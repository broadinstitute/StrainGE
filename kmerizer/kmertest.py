import kmerizer

seq = 'agcttttcattctgactgcaacgggcaatatgtctctgtgtggattaaaaaaagagtgtc'

k15 = kmerizer.kmerize(15, seq)
k15.sort()
print k15

k23 = kmerizer.kmerize(23, seq)
k23.sort()
print k23

k31 = kmerizer.kmerize(31, seq)
k31.sort()
print k31
