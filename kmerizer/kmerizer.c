#include <Python.h>
#include <python2.7/numpy/ndarrayobject.h>

#define A 0
#define C 1
#define G 2
#define T 3
#define RC(b) ((b) ^ 3)

typedef uint64_t kmer_t;


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

static PyObject *KmerizerError;

static PyMethodDef KmerizerMethods[] = {
      {"kmerize",  kmerizer_kmerize, METH_VARARGS,
       "Kmerize sequence returning new numpy array."},
      {"kmerize_into_array",  kmerizer_kmerize_into_array, METH_VARARGS,
       "Kmerize sequence into existing numpy array."},
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
