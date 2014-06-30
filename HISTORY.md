FORTPY: Revision History
======

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