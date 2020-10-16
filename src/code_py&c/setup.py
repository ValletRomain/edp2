from distutils.core import setup, Extension

module1 = Extension('godunov',
                    sources = ['godunov.c'])

setup (name = 'Godunov',
       version = '1.0',
       description = 'This is a demo godunov package',
       ext_modules = [module1])
