#!/usr/bin/env python
try:
    from setuptools import setup
    args = {}
except ImportError:
    from distutils.core import setup
    print("""\
*** WARNING: setuptools is not found.  Using distutils...
""")

setup(name='Fortpy',
      version='1.0.3',
      description='Fortran Parsing, Unit Testing and Intellisense',
      author='Conrad W Rosenbrock',
      author_email='rosenbrockc@gmail.com',
      url='',
      install_requires=[
          "argparse",
          "pyparsing",
          "paramiko",
      ],
      packages=['fortpy', 'fortpy.parsers', 'fortpy.isense', 'fortpy.testing',
            ],
     )
