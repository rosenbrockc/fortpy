#!/usr/bin/env python
try:
    from setuptools import setup
    args = {}
except ImportError:
    from distutils.core import setup
    print("""\
*** WARNING: setuptools is not found.  Using distutils...
""")

from setuptools import setup
try:
    from pypandoc import convert
    read_md = lambda f: convert(f, 'rst')
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    read_md = lambda f: open(f, 'r').read()

setup(name='Fortpy',
      version='1.0.4',
      description='Fortran Parsing, Unit Testing and Intellisense',
      long_description=read_md('README.md'),
      author='Conrad W Rosenbrock',
      author_email='rosenbrockc@gmail.com',
      url='https://github.com/rosenbrockc/fortpy',
      license='MIT',
      install_requires=[
          "argparse",
          "pyparsing",
          "python-dateutil",
          "paramiko",
      ],
      packages=['fortpy', 'fortpy.parsers', 'fortpy.isense', 'fortpy.testing',
                'fortpy.templates', 'fortpy.interop', 'fortpy.scripts',
                'fortpy.printing' ],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'Natural Language :: English',
          'License :: OSI Approved :: MIT License',          
          'Operating System :: MacOS',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
     )
