from distutils.core import setup, Extension

setup (name = 'kmertools',
       version = '1.0',
       author = 'Bruce Walker',
       author_email = 'bruce@w1bw.us',
       description = 'Kmer generation and manipulation modules',
       py_modules = ['kmertools'],
       ext_modules = [Extension('kmerizer', sources = ['kmerizer.c'])])





