#!/usr/bin/env python
try:
    from setuptools import setup
    args = {}
except ImportError:
    from distutils.core import setup
    args = dict(scripts=['jediepcserver.py'])
    print("""\
*** WARNING: setuptools is not found.  Using distutils...
""")

setup(name='Fortpy',
      version='1.0',
      description='Fortran Parsing, Unit Testing and Intellisense',
      author='Conrad W Rosenbrock',
      author_email='rosenbrockc@gmail.com',
      url='',
      # install_requires=[
      #     "argparse",
      #     "pyparsing",
      # ],
      packages=['fortpy', 'fortpy.parsers', 'fortpy.isense', 'fortpy.testing',
            ],
     )
