[![PyPI](https://img.shields.io/pypi/v/fortpy.svg)](https://pypi.python.org/pypi/fortpy/)

# FORTPY

## Python Emacs Intellisense and Unit Testing Support for Fortran

Fortpy is a python based parsing, unit testing and auto-complete framework for supporting Fortran 2003 including object oriented constructs. Auto-completion integration currently only available for emacs. Here are some of the features:

### Python Package Installation

Fortpy requires the following python packages:
- [python-epc](https://github.com/tkf/python-epc)
- argparse
- [paramiko](https://github.com/paramiko/paramiko) if you will be editing over SSH using tramp in emacs and still want auto-complete support.
- pyparsing
- dateutil

If you install Fortpy using `pip install` the dependencies will automatically get installed.

    pip install fortpy

### Fortpy Configuration

In most real-usage scenarios, the out-of-the-box support for multiple libraries is incomplete. A [configuration file](https://github.com/rosenbrockc/fortpy/wiki/Fortpy-Global-Configuration) can be created that gives additional information about source code folders to parse, server information for SSH editing over tramp and auto-completion configuration settings. After you have created a `config.xml` file on your system, you need to configure an environment variable to tell Fortpy where to find it:

`export FORTPY_CONFIG="~/path/to/config.xml"`

NOTE: Environment variables for emacs are set when it first starts; adding a config.xml file will not affect the emacs isense support until you restart emacs.

## Screenshots

Here are some of the things you can do once Fortpy is integrated with Emacs using fortpy.el:

![Automatic Signature Suggestions](../master/docs/screenshots/signature.png "Help with call signatures of functions and subroutines.")

Help with call signatures of functions and subroutines.

![Embedded Member Suggestions](../master/docs/screenshots/completion.png "Completion suggestions for both methods and variables embedded in user-defined types.")

Completion suggestions for both methods and variables embedded in user-defined types.

![Bracket Complete Embedded Methods](../master/docs/screenshots/bracket_complete.png "Documentation strings for methods embedded in user-defined types")

Documentation strings for methods embedded in user-defined types.
 
## Features

### Auto-complete Support
- Integration with emacs using virtualenv, epc and auto-complete.
- Context-specific suggestions for variable names and available executables.

See the [isense package documentation](https://github.com/rosenbrockc/fortpy/wiki/Intellisense-Package)

### XML Documentation Standard

Fortpy uses an [XML documentation standard](https://github.com/rosenbrockc/fortpy/wiki/XML-Documentation-Standard) that is derived from the Microsoft XML code documentation standard. As such, many of the same tags and attributes are available for decorating the methods, types and variables in a Fortran module. Because the parsers were designed to work with Fortran 2003, they understand embedded methods in user types that point to other methods in the module. See the wiki page on documenting your code to work with Fortpy.

![XML Documentation Example](../master/docs/screenshots/xml_docs.png "XML documentation standard allows for complex documentation strings and structures.")

### Automated Unit Testing

Because Fortpy is aware of the structure and dependencies between modules, it can generate automated unit tests for individual functions. You can specify the tests to run in some XML documentation decorating the method to test. When the unit testing scripts are run, a fortran program gets generated and all the necessary dependencies are automatically included. After running the program, the framework tests the output and generates a similarity score. [See the full documentation](https://github.com/rosenbrockc/fortpy/wiki/Unit-Testing-Package).

### Code Version Support

Fortpy's unit testing framework makes use of XML templates that specify the structure of input/output files to run and test. Using these same templates, the interoperability package allows input and ouput files to be converted between versions of the Fortran program. Since most Fortran programs run on simple text files for input, the interop package also gives support for XML input files to be used; scripts convert the XML input files to the correct plain-text format at run-time.
