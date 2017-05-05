from distutils.core import setup, Extension
import numpy

setup (name = 'kmertools',
       version = '1.0',
       author = 'Bruce Walker',
       author_email = 'bruce@w1bw.us',
       description = 'Kmer generation and manipulation modules',
       py_modules = ['kmertools'],
       ext_modules = [Extension('kmerizer', sources = ['kmerizer.c'])],
       include_dirs=[numpy.get_include()])





