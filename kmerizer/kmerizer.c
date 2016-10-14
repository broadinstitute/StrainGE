#include <Python.h>
#include <python2.7/numpy/ndarrayobject.h>

#define A 0
#define C 1
#define G 2
#define T 3
#define RC(b) ((b) ^ 3)

typedef uint64_t kmer_t;
typedef int64_t count_t;

int kmerizer_kmerize_internal(int k, char *seq, int length, kmer_t *kmers) {
  int kcount = 0;

  kmer_t mask = ((kmer_t) 1 << (2 * k)) - 1;
  int shift = 2 * (k - 1);

  int n = 0;
  kmer_t fw = 0;
  kmer_t rc = 0;

  for (int i = 0; i < length; ++i) {
    char b = toupper(seq[i]);
    kmer_t value;
    switch (b) {
    case 'A':
      value = A;
      break;
    case 'C':
      value = C;
      break;
    case 'G':
      value = G;
      break;
    case 'T':
      value = T;
      break;
    default:
      fw = rc = n = 0;
      continue;
    }
    fw = ((fw << 2) & mask) | value;
    rc = ((rc >> 2) & mask) | (RC(value) << shift);
    if (++n >= k) {
      kmer_t kmer = (fw > rc) ? fw : rc;
      *kmers++ = kmer;
      kcount++;
    }
  }
  return kcount;
}

static PyObject *
kmerizer_kmerize(PyObject *self, PyObject *args)
{
  int k;
  char *seq;
  int length;

  if (!PyArg_ParseTuple(args, "is#", &k, &seq, &length))
    return NULL;

  int max_kmers = length + 1 - k;

  kmer_t *kmers = malloc(max_kmers * sizeof(kmer_t));
  int kcount = kmerizer_kmerize_internal(k, seq, length, kmers);

  npy_intp dim = kcount;
  PyObject *result = PyArray_SimpleNew(1, &dim, NPY_ULONG);
  kmer_t *kmerp = PyArray_GETPTR1(result, 0);
  bcopy(kmers, kmerp, kcount * sizeof(kmer_t));
  free(kmers);

  return result;
}

static PyObject *
kmerizer_kmerize_into_array(PyObject *self, PyObject *args)
{
  int k;
  char *seq;
  int length;
  PyObject *array;
  int offset;

  if (!PyArg_ParseTuple(args, "is#Oi", &k, &seq, &length, &array, &offset))
    return NULL;

  int max_kmers = length + 1 - k;

  kmer_t *kmers = malloc(max_kmers * sizeof(kmer_t));
  int kcount = kmerizer_kmerize_internal(k, seq, length, kmers);

  kmer_t *kmerp = PyArray_GETPTR1(array, offset);
  bcopy(kmers, kmerp, kcount * sizeof(kmer_t));
  free(kmers);
  return Py_BuildValue("i", kcount);
}


/* for two sorted arrays, count the number of common elements */
int
kmerizer_count_common_internal(kmer_t *kmers1, count_t size1, kmer_t *kmers2, count_t size2) {
  /* walk through sorted arrays in parallel */
  int kcount, i1, i2;
  for (kcount = 0, i1 = 0, i2 = 0; i1 < size1 && i2 < size2; ) {
    kmer_t kmer1 = kmers1[i1];
    kmer_t kmer2 = kmers2[i2];
    /*printf("count=%d i1=%d i2=%d k1=%lx k2=%lx c1=%d c2=%d\n", kcount, i1, i2, kmer1, kmer2, c1[i1], c2[i2]);*/
    if (kmer1 == kmer2) {
      /* in both */
      ++kcount;
      ++i1;
      ++i2;
    } else if (kmer1 < kmer2) ++i1;
    else ++i2;
  }
  return kcount;
}


static PyObject *
kmerizer_count_common(PyObject *self, PyObject *args)
{
  PyObject *kmers1, *kmers2;
  if (!PyArg_ParseTuple(args, "OO", &kmers1, &kmers2))
    return NULL;
  npy_intp k1size = *PyArray_DIMS(kmers1);
  npy_intp k2size = *PyArray_DIMS(kmers2);
  kmer_t *k1 = PyArray_GETPTR1(kmers1, 0);
  kmer_t *k2 = PyArray_GETPTR1(kmers2, 0);
  int kcount = kmerizer_count_common_internal(k1, k1size, k2, k2size);
  return Py_BuildValue("i", kcount);
}

static PyObject *
kmerizer_merge_counts(PyObject *self, PyObject *args)
{
  PyObject *kmers1, *count1, *kmers2, *count2;

  if (!PyArg_ParseTuple(args, "OOOO", &kmers1, &count1, &kmers2, &count2))
    return NULL;
  
  npy_intp k1size = *PyArray_DIMS(kmers1);
  npy_intp k2size = *PyArray_DIMS(kmers2);

  kmer_t *k1 = PyArray_GETPTR1(kmers1, 0);
  kmer_t *k2 = PyArray_GETPTR1(kmers2, 0);
  count_t *c1 = PyArray_GETPTR1(count1, 0);
  count_t *c2 = PyArray_GETPTR1(count2, 0);

  int kcount, i1, i2;;

  kcount = k1size + k2size - kmerizer_count_common_internal(k1, k1size, k2, k2size);

  /* allocate output arrays */
  npy_intp dim = kcount;
  PyObject *kmersOut = PyArray_SimpleNew(1, &dim, NPY_ULONG);
  kmer_t *kout = PyArray_GETPTR1(kmersOut, 0);

  PyObject *countsOut = PyArray_SimpleNew(1, &dim, NPY_LONG);
  count_t *cout = PyArray_GETPTR1(countsOut, 0);

  /* walk through sorted arrays in parallel */
  for (kcount = 0, i1 = 0, i2 = 0; i1 < k1size && i2 < k2size; ) {
    kmer_t kmer1 = k1[i1];
    kmer_t kmer2 = k2[i2];
    /*printf("count=%d i1=%d i2=%d k1=%lx k2=%lx c1=%d c2=%d\n", kcount, i1, i2, kmer1, kmer2, c1[i1], c2[i2]);*/
    if (kmer1 == kmer2) {
      /* in both */
      kout[kcount] = kmer1;
      cout[kcount] = c1[i1++] + c2[i2++];
    } else if (kmer1 < kmer2) {
      /* only in kmers1 */
      kout[kcount] = kmer1;
      cout[kcount] = c1[i1++];
    } else {
      /* only in kmers2 */
      kout[kcount] = kmer2;
      cout[kcount] = c2[i2++];
    }
    ++kcount;
  }
  /* leftovers in kmers1 */
  while (i1 < k1size) {
    kout[kcount] = k1[i1];
    cout[kcount++] = c1[i1++];
  }
  /* leftovers in kmers2 */
  while (i2 < k2size) {
    kout[kcount] = k2[i2];
    cout[kcount++] = c2[i2++];
  }
  
  return Py_BuildValue("(NN)", kmersOut, countsOut);
}

static PyObject *KmerizerError;

static PyMethodDef KmerizerMethods[] = {
      {"kmerize",  kmerizer_kmerize, METH_VARARGS,
       "Kmerize sequence returning new numpy array."},
      {"kmerize_into_array",  kmerizer_kmerize_into_array, METH_VARARGS,
       "Kmerize sequence into existing numpy array."},
      {"merge_counts",  kmerizer_merge_counts, METH_VARARGS,
       "Merge and sum  two sets of kmer and count arrays."},
      {"count_common",  kmerizer_count_common, METH_VARARGS,
       "Count common elements in two sorted kmer arrays."},
      {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initkmerizer(void)
{
  PyObject *m;

  m = Py_InitModule("kmerizer", KmerizerMethods);
  if (m == NULL)
    return;

  import_array();

  KmerizerError = PyErr_NewException("kmerizer.error", NULL, NULL);
  Py_INCREF(KmerizerError);
  PyModule_AddObject(m, "error", KmerizerError);
}
