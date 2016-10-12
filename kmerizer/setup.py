from distutils.core import setup, Extension

module1 = Extension('kmerizer',
                    sources = ['kmerizer.c'])

setup (name = 'kmerizer',
              version = '1.0',
              description = 'This is a kmer generation package',
              ext_modules = [module1])



