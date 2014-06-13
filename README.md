FORTPY
======

Python Emacs Intellisense and Unit Testing Support for Fortran
------

Fortpy is a python based parsing, unit testing and auto-complete framework for supporting Fortran 2003 including object oriented constructs. Auto-completion integration currently only available for emacs. Here are some of the features:

### Auto-complete Support
- Integration with emacs using virtualenv, epc and auto-complete.
- Completion suggestions for both methods and variables embedded in user-defined types.
- Signature prompts when calling functions or subroutines.
- Context-specific suggestions for variable names and available executables.

### XML Documentation Standard

Fortpy uses an XML documentation standard that is derived from the Microsoft XML code documentation standard. As such, many of the same tags and attributes are available for decorating the methods, types and variables in a Fortran module. Because the parsers were designed to work with Fortran 2003, they understand embedded methods in user types that point to other methods in the module. See the wiki page on documenting your code to work with Fortpy.

### Automated Unit Testing

