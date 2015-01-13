FORTPY: Revision History
======

Revision 1.3.5
------

- Fixed a bug with function fitting introduced with Revision 1.3.5.
- Added fit-related dependent variable removal from `self.curargs["dependents"]` when `rmfit` is used.
- Fixed issue #28

Revision 1.3.4
------

- Postfix functions weren't being applied anymore to the data because of the filter inclusion in each dependent variable. Fixed.

Revision 1.3.3
------

- Added property editing for all variables at once using `*` as the variable name.
- Fixed a bug when setting plot line widths.
- Added additional support for postfix functions that return `None` values.
- Fixed bug where line plots weren't including lines by default if any marker settings were given.

Revision 1.3.2
------

- Added option to set x and y axis plot limits.
- Added option for dependent variables to be plotted either with scatter *or* with plot. In the future it will be easier to add other options if we need to.
- Fixed some minor bugs from previous revisions uncovered during the enhancements.

Revision 1.3.1
------

- Allowed adding dependent variables from multiple filters for plotting and tabulating.
- Added plot options for `markers`, `lines` and `ticks`.
- Added font options for labels, axes, title and legend.
- Added the option to save and load `font` and `ticks` options to file; this allows publication quality settings to be saved and loaded into any analysis group easily.
- This revision breaks the saved sessions from previous versions of the shell. Increment minor version.
- Fixed a problem in `fortpy.f90` where tab-delimited files were not being read in correctly.

Revision 1.2.17
------

- Added the option `[None]` as a possible label to remove a series from the plot legend in the test analysis shell.
- Added ranges for specifying the cases in a unit test so that developers don't have to type out series of tests with sequential integers.

Revision 1.2.16
------

- Debugged the handling of multi-dimensional data for arrays restricted in one of its dimensions combined with an aggregation function.
- Added handling for fitting parameters that have values below `1e-2` so that their equations are displayed in scientific notation instead of plain decimal.
- Added commands `rmcolor`, `rmfit`, `rmlabel` and `rmpostfix` for removing these variable to value mappings.

Revision 1.2.15
------

- Added shell restarts on unhandled exceptions; the default is to restart 10 times before saving the shell and exiting.
- Added handling so that multiple dimensional data will still be attempted, but with a warning; this allows data that is 1D in certain dimensions to still be plotted/tabulated with an appropriate aggregation function across test cases.
- Added automatic fitting of linear and exponential functions to the plots.
- Added colors and plot legend/title settings/commands to the shell.
- Debugged the shell history getting extremely large from shell restarts.
- Added the option to view/clear the last unhandled exception caught by the shell.

Revision 1.2.14
------

- Fixed a bug in the test case filter that was matching nothing, so that no results ever made it to plots/tables.

Revision 1.2.13
------

- Added output of exception message in generalized unhandled exception trapper.
- Added `numpy` and `matplotlib` as pre-reqs for the fortpy installation.

Revision 1.2.12
------

- Fixed a bug introduced in revision 1.2.11 with path completions involving `~`.
- Added validation for commands that require variables to be set; previously it caused an unhandled exception.
- Added global exception handling for commands. Instead of just dying, the session gets saved automatically to `#fortpy.shell#` before exiting when an unhandled exception occurs.

Revision 1.2.11
------

- Added support for directory navigation withing fortpy shell using `cd`, `pwd`, `ls`.
- Fixed support for relative directories such as `..` when selecting files or folders.

Revision 1.2.10
------

- Fixed a bug in `scripts/analyze.py` for loading the shell history. For first time users with no history file, the script broke with a file not found error.

Revision 1.2.9
------

- Fixed a bug with imports in `/scripts/analyze.py` when the `-pypath` parameter isn't used.
- Added full documentation to the fortpy unit test analysis shell in `scripts/analyze.py`.

Revision 1.2.8
------

- Fixed bug in `scripts/compare.py` that broke when the file comparer was updated during a previous upgrade of the unit testing framework. Ironically, the lack of unit tests for fortpy didn't catch the bug introduced in the framework for performing fortran unit tests.
- Added the module `testing/parser.py` to parse test case results and performing generic functions such as plotting, creating tables, etc.
- Added the fortpy shell `scripts/analyze.py` that interfaces with test results to tabulate results, create plots, investigate failures, etc.
- Fixed a bug that raised an exception when the model output files are empty.
- Fixed a bug where input files with more than 500 characters per line don't get all the data imported. Upped the value to 5000.
- Fixed a bug in the test input/output file generator where the generic value writer didn't handle the case of a `long` in python.
- Added overflow error checking for the make files when compiling in DEBUG mode. Changed the default behavior of the `scripts/runtests.py` script to always run the fortpy unit tests in DEBUG mode unless the `-nodebug` flag is specified.

Revision 1.2.7
------

- Fixed a bug that added the `pointer` modifier to class variables even if allocatable or pointer was already present as a modifier. The problem was an `or` that should have been an `and`.
- Decided that the previous revision 1.1.6 should have included a minor increment because of the 1000+ lines of new/changed code. Switching this revision to 1.2.7 instead of 1.1.7.

Revision 1.1.6
------

- Fixed bug with uninitialized list of parsed programs as part of trying to support fortran PROGRAM files.
- Fixed bug for goto definition lookups not working in cases where they should.
- Fixed bug where global definitions for variables that are never in a parameter list can't be auto-generated.
- Fixed bug for variable value assignments from embedded methods in derived types.
- Fixed bug where allocation for local variables in the unit testing wasn't been generated for embedded methods.
- Fixed bug in variable dimension and default value parsing. The regex was stopping at the first `)` and then stripping all other parentheses in the string. The end result was that the property `D` was wrong because it only counted `,` instances in the string and when dimensions were specified with functions, there were lots of extra commas.
- Added support for ragged array input files. If the file is 2D but a 1D array of objects needs to be initialized from each line in the file (with lines of differing length), then `ragged="true"` can be added to the `<value>` tag and the ragged initialization is handled automatically.
- Added support for `<part>` tags so that specific dimensions of arrays and arrays of objects can be initialized separately. This includes wildcard support for input file copying and initialization of higher-dimensional arrays from a series of 2D arrays.
- Added support for allocate-only assignments of variables. An `<assignment>` tag can now have just a `name` and `allocate` attribute. Unless the variable is a derived-type pointer, the `allocate` attribute must explicitly specify the dimensions of the array to allocate.
- Added support for mixed file, function and embedded value assignments. This is especially useful in the context of ragged arrays and part assignments. For example, an array of objects where each elements needs to be initialized from a 2D input file via an embedded function is possible using a combination of `<part>` and `<value>` tags where the `<value>` tag has both `file` and `embedded` attributes specified.
- Fixed a bug where the `allocatable` or `pointer` modifiers were not being automatically added to `class` and `type` variables with `regular="true"`.

Revision 1.1.5
------

- Fixed issue #20
- Fixed a "No summary for element." error in the intellisense that was returning no summary when one actually existed.

Revision 1.1.4
------

- Fixed an error in `fortpy.f90` where a failing condition on the first part of an if statement that used `.and.` wasn't exiting the condition without evaluating the second part. This seems to be a quirk of the gfortran compiler.
- Fixed issue #19
- Fixed issue #18
- Fixed issue #17

Revision 1.1.3
------

- Reversed the renaming of the variables kinds (e.g. `dp`) in `fortpy.f90` because they were breaking unit tests that previously functioned correctly. Opened an issue to address the problem.

Revision 1.1.2
------

- Fixed a bug where parameters with `regular="true"` weren't being found because the dictionary had lower-case keys and the variable name was upper-case.
- Fixed parser bug for displaying and saving modifiers on executables.
- Fixed parser bug for determining the dependency modules.
- Added support for identifying modules that require pre-processing before compilation.
- Added configuration section `<libraries>` to fortpy `config.xml`. Allows libraries to be explicitly specified for inclusion when linking the unit testing executable.
- Added rule to the auto-generated makefile for creating a shared library of the executable being unit tested and all of its dependencies.

Revision 1.1.1
------

- Fixed parsing issues when identifying dependency calls inside of executables. This also fixed some isense issues caused by broken dependencies definitions.

Revision 1.1.0
------

- Added multiple cache support so that unit testing doesn't affect the "production" cache.
- Fixed versioning issue where the cache was not being invalidated when incompatible changes are made to the structure of classes that get serialized.
- Minor version release because of additional unit testing support and new support for interfaces.

Revision 1.0.13
------

- Added feature for auto-generating the global variable for custom types whose embedded methods are being tested. Usually, a `regular="true"` attribute is required to auto-generate a global variable. However, when embedded types are being called, they seldomly get commented with `<parameter>` tags because it is well understood that the first argument is the type/class instance.
- Fixed issue #9
- Added limited support for interfaces. Generic interfaces are supported with isense now. Although the parsers will handle the operator and assignment interfaces, support is not yet provided for actually using the symbols defined in those interfaces.
- Fixed issue #11
- Added more intelligent completions for embedded methods in derived types. Previously, all members and module procedures in an embedded type were returned with %-completion; now if the symbol is preceded by `call`, only embedded subroutine definitions are returned.

Revision 1.0.12
------

- Debugged error where the first parameter to an embedded procedure on a derived type was still being explicitly passed in to the embedded procedure. Since the compiler does that automatically, it was choking.

Revision 1.0.11
------

- Debugged the first real unit test case of an embedded method unit test.
- Debugged some additional allocation problems with auto-generated parameters in the unit testing.

Revision 1.0.10
------

- Fixed a bug that was preventing deferred shape arrays from being allocated automatically by fortpy in unit testing.

Revision 1.0.9
------

- Fixed a bug when comparing two simple representations where the second file has more entries than the first one.
- Fixed a misrepresentation of match success in the testing. The framework was using `common_match` in the `percent_match` routine which made wrong results seem a lot better.
- Added coloring to the terminal output for all packages.

Revision 1.0.8
------

- Fixed bug in `fortpy.pysave` that made 2d saving not produce tabular input.
- Fixed bug that crashed the program when `<isense>` was not present in the `config.xml` file.
- Added support for automated profiling of methods being tested using the -profile argument to `scripts/runtests.py`.

Revision 1.0.7
------

- Added `timed="true"` attribute support to the `<test>` tag.
- Fixed bug where an empty folder being cleaned would raise an exception on `remove()`
- Added a check for executables being tested to make sure they have a `public` modifier or are in the list of public symbols in the parent module.
- Added debug compile option to unit testing.
- Fixed automatic allocation errors for pointers.
- Fixed fortpy.f90 line and value counting routines.
- Added support for arrays of derived types; still get suggestions for embedded members if operating on an element of an array of derived types.
- Fixed a bug in `parsers.types` where "endtype" was not recognized.
- Added `<isense>` section to `config.xml` with an option to uppercase the builtin Fortran names
  by default with the suggestions.

Revision 1.0.6
------

- Fixed issue with fortpy/templates directory not being installed.

Revision 1.0.5
------

- Fixed issue with the builtin.xml file not being included in the distribution.

Revision 1.0.4
------

- Added support for gfortran compilation in unit testing
- Altered testing XML template to use `<test>` instead of `<outcome>`
  - Allow multiple inputs using `<input>` tags instead of comma-separated list.
  - Allow multiple output comparisons with `<output>` tags instead of comma-separated list.
  - Added _description_ attribute to `<test>` tag.
- Constant-input mode for testing. When a `<line>` template is specified for an `<input>` tag, the call to the main method being unit tested can be repeated for each line in that input file.
- Altered testing XML template to use `<assignment>` instead of `<instance>`.
- Variable value assignments directly from files.
  - `<assignment>` tag now accepts additional attributes _folder_, _file_ and _rename_.
- Major overhaul of underlying classes for the unit testing to make future feature enhancements much easier. Added the `elements` module in the `isense` package.
- Added support for built-in functions using information documented from a textbook that is now in `builtin.xml`. Added the `builtin` module to parse and cache that XML documentation.