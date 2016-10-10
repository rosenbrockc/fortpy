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
      version='1.7.7',
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
          "termcolor",
          "numpy",
          "matplotlib",
          "scipy",
          "tqdm"
      ],
      packages=['fortpy', 'fortpy.parsers', 'fortpy.isense', 'fortpy.testing',
                'fortpy.templates', 'fortpy.interop',
                'fortpy.printing', 'fortpy.stats' ],
      scripts=['fortpy/scripts/compare.py', 'fortpy/scripts/convert.py', 'fortpy/scripts/runtests.py',
               'fortpy/scripts/analyze.py', 'fortpy/scripts/parse.py', 'fortpy/scripts/ftypes.py',
               'fortpy/scripts/bestprac.py'],
      package_data={'fortpy': ['isense/builtin.xml']},
      include_package_data=True,
      classifiers=[
          'Development Status :: 5 - Production/Stable',
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
