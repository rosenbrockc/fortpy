"""This module generates the contents of the fortpy.f90 file, which has
lots of duplication by nature of the interfaces that it implements.
"""
kmap = {
    "integer": {"dp": "fli", "sp": "fsi"},
    "complex": {"dp": "fdp", "sp": "fsp"},
    "real": {"dp": "fdp", "sp": "fsp"},
    "character": {"string": "len=*"}
}
linemax = 250000
"""Specifies the maximum length that any line in a data file can have and still
be supported by fortpy.
"""
prec = {
    "integer": {None: "'('// adjustl(FMT) // 'i25)'",
                "dp": "'('// adjustl(FMT) // 'i50)'",
                "sp": "'('// adjustl(FMT) // 'i5)'"},
    "complex": {"dp": "'('// adjustl(FMT) // 'e22.12)'", "sp": "'('// adjustl(FMT) // 'e22.12)'"},
    "real": {"dp": "*", "sp": "*"},
    "logical": {None: "*"},
    "character": {None: "'(A)'"}
}
"""The fortran format strings for writing the values for each data type and kind
to file.
"""
xsuffix = {
    None: "",
    "_p": "p",
    "_f": "f"
}
"""Suffixes for the executable names indexed by the pointer/fixed dimension suffixes.
"""
defaults = {
    "integer": "0",
    "real": "0",
    "complex": 0,
    "logical": ".false.",
    "character": "''"
}
def spacer(n):
    """Returns a string with 'n' spaces."""
    if n > 0:
        return ''.join([' ']*n)
    else:
        return ""

def fpy_read(D, dtype, kind, suffix=None):
    """Generates the fortran code for the fpy_read interface module procedure
    of the specified type and dimensionality.
    
    :arg D: the integer dimensionality of the variable that is having its
      value set from the file.
    :arg dtype: the fortran data type that the subroutine will be reading data
      for. Can be ["real", "integer", "character", "logical", "complex"].
    :arg kind: the precision/kind information for the type. Can be one of
      ["dp", "sp"] for double/single precision real, integer or complex types.
    :arg suffix: whether to generate the module procedure for allocatable (None),
      pointer ("_p"), or fixed ("_f) variables.
    """
    modifiers = {
        None: ", allocatable",
        "_p": ", pointer",
        "_f": ""
    }
    common = {
        "dtype": dtype,
        "D": D,
        "modify": modifiers[suffix] if D > 0 else "",
        "Dx": "({})".format(','.join([":"]*D)) if D >0 else "",
        "maxlen": linemax,
        "suffix": xsuffix[suffix],
        "skind": kind if kind is not None else ""
    }

    if dtype in defaults:
        common["default"] = defaults[dtype]
    else:
        common["default"] = "0"
        
    if kind is not None and dtype in kmap:
        common["kind"] = "({})".format(kmap[dtype][kind])
    else:
        common["kind"] = ""
    
    #Add the specific sections for each of the dimensionalities that we read.
    if D==0:
        common["vars"] = "integer :: nlines, nvalues"
    elif D in [1, 2]:
        common["vars"] = "integer :: nlines, nvalues, i"
    else:
        fmtstr = "integer :: dims({0:d}), {1}, indices({0:d})"
        common["vars"] = fmtstr.format(D, ', '.join(["i{}".format(i+1) for i in range(D)]))

    if suffix == "_p":
        common["dealloc"] = "    if (associated(variable)) variable => null()\n"
    elif D > 0 and suffix is None:
        common["dealloc"] = "    if (allocated(variable)) deallocate(variable)\n"
    else:
        common["dealloc"] = ""
        
    if D in [0,1,2]:
        if dtype == "character":
            common["analyze"] = "    call fpy_linevalue_count(filename, commentchar, nlines, nvalues, .true.)\n"
        else:
            common["analyze"] = "    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)\n"
    else:
        common["analyze"] = ""

    if D==0:
        lines = ["if ((nvalues .gt. 1) .or. (nlines /= 1)) then",
                 "  write(*,*) \"Cannot read a single value from \", filename",
                 "  write(*,*) \"Found \", nlines, \" lines and \", nvalues, \" values\"",
                 "  if (present(success_)) success_ = .false.",
                 "  if (strict) stop",
                 "end if\n"]
        common["warning"] = '    '+'\n    '.join(lines)
    elif D==1:
        lines = ["if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then",
                 '  write(*,*) "Cannot read a vector value from ", filename',
                 '  write(*,*) "Found ", nlines, " lines and ", nvalues, " values"',
                 "  if (present(success_)) success_ = .false.",
                 "  if (strict) stop",
                 "end if\n"]
        if suffix == "_f":
            lines.extend(["",
                          "if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then",
                          "  write(*,*) \"File data dimensions don't match fixed variable shape \", shape(variable)",
                          "  write(*,*) \"Fortpy sees data dimensions in '\", filename, \"' as \", nlines, nvalues",
                          "  if (present(success_)) success_ = .false.",
                          "  if (strict) stop",
                          "end if\n"])
        common["warning"] = '    '+'\n    '.join(lines)
    elif D==2 and suffix == "_f":
        lines = ["if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then",
                 "  write(*,*) \"File data dimensions don't match fixed variable shape \", shape(variable)",
                 "  write(*,*) \"Fortpy sees data dimensions in '\", filename, \"' as \", nlines, nvalues",
                 "  if (present(success_)) success_ = .false.",
                 "  if (strict) stop",
                 "end if\n"]
        common["warning"] = '    '+'\n    '.join(lines)
    else:
        common["warning"] = ""

    if D==1 and suffix in [None, "_p"]:
        lines = ["if (nlines .gt. 1) then",
                 "  allocate(variable(nlines))",
                 "else",
                 "  allocate(variable(nvalues))",
                 "end if",
                 "variable = {default}\n".format(**common)]
        common["allocate"] = '    '+'\n    '.join(lines)
    elif D==2 and suffix in [None, "_p"]:
        common["allocate"] = ("    allocate(variable(nlines, nvalues))\n"
                              "    variable = {default}\n".format(**common))
    elif D==0:
        common["allocate"] = "    variable = {}\n".format(common["default"])
    else:
        common["allocate"] = ""

    if D==0:
        common["initial"] = ""    
    elif D in [1,2]:
        common["initial"] = "\n      i=1"
    else:
        lines = ["\n      do",
                 '  read(funit, "(A)", iostat=ioerr) line',
                 "  if (ioerr == 0) then",
                 "    cleaned = trim(adjustl(line))",
                 "    if (len(cleaned) .gt. 1) then",
                 "      if (cleaned(1:2) .eq. '##') then",
                 "        read(cleaned(3:) ,*) dims",
                 "        exit",
                 "      end if",
                 "    end if",
                 "  else",
                 "    write(*,*) \"No shape information found for {0}D data file '\", filename, \"'\"",
                 "    stop",
                 "  end if",
                 "end do",
                 ""]
        if suffix in [None, "_p"]:
            dims = ', '.join(["dims({})".format(i+1) for i in range(D)])
            lines.append("allocate(variable({}))".format(dims))
        lines.append("variable = {default}".format(**common))
        lines.append("")

        common["initial"] = '\n      '.join(lines).format(D)

    if D==0:
        lines = ["if (cleaned(1:1) /= commentchar) then",
                 "  read(line, *) variable",
                 "end if"]
        common["readvar"] = '\n            '.join(lines)
    elif D==1:
        lines = ["if (cleaned(1:1) /= commentchar) then",
                 "  if (nlines .gt. 1) then",
                 "    read(line, *) variable(i)",
                 "    i = i+1",
                 "  else",
                 "    read(line, *) variable",
                 "  end if",
                 "end if"]
        common["readvar"] = '\n            '.join(lines)
    elif D==2:
        lines = ["if (cleaned(1:1) /= commentchar) then",
                 "  read(line, *) variable(i,:)",
                 "  i = i+1",
                 "end if"]
        common["readvar"] = '\n            '.join(lines)
    else:
        lines = ["if (cleaned(1:2) .eq. '##') then",
                 "  read(cleaned(3:), *) indices"]
        ivars = ["i{}".format(i+1) for i in range(D-1)]
        for i in range(D-2):
            lines.append("  i{0} = indices({0})".format(i+1))
        lines.append("  i{} = 1".format(D-1))
        
        lines.extend(["elseif (cleaned(1:1) /= commentchar) then",
                      "  read(line, *) variable({0},:)",
                      "  i{1} = i{1}+1",
                      "end if"])
        common["readvar"] = '\n            '.join(lines).format(','.join(ivars), D-1)

    if dtype=="character":
        common["xname"] = "fpy_read_{dtype}_{D}d{suffix}".format(**common)
    else:
        common["xname"] = "fpy_read_{dtype}{skind}_{D}d{suffix}".format(**common)

    if suffix is None:
        common["noexist"] = "return"
    elif suffix == "_p":
        #Set pointers to null if they weren't set on construction in their definition.
        common["noexist"] = "variable => null()\n      return"
    else:
        common["noexist"] = "variable = {}\n      return".format(common["default"])
        
    template = """  subroutine {xname}(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    {dtype}{kind}{modify}, intent(inout) :: variable{Dx}

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    {vars}
    character({maxlen:d}) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      {noexist}
    end if

{dealloc}{analyze}{warning}{allocate}

    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then{initial}
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            {readvar}
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine {xname}
    """
    return (common["xname"], template.format(**common))

def fpy_interface_sub(fpy, dtype, kind, suffix=None):
    """Generates the interface for reading/writing values for variables of the specified suffix
    for dimensions 0 through 7.

    :arg fpy: the generating function that accepts the dimension and type information
      and returns the (xnames, subtext).
    :arg suffix: one of [None, "_f", "_p"] corresponding to "allocatable", fixed dimension,
      and "pointer" type variables.
    """
    allsubs = []
    xnames = []
    drange = list(range(8)) if suffix is None else list(range(1,8))
    for D in drange:
        xname, sub = fpy(D, dtype, kind, suffix)
        xnames.append(xname)
        allsubs.append(sub)

    return (xnames, '\n'.join(allsubs))

def fpy_interface(fpy, static, interface, typedict):
    """Splices the full list of subroutines and the module procedure list
    into the static.f90 file.

    :arg static: the string contents of the static.f90 file.
    :arg interface: the name of the interface *field* being replaced.
    :arg typedict: the dictionary of dtypes and their kind and suffix combos. 
    """
    modprocs = []
    subtext = []
    for dtype, combos in list(typedict.items()):
        for tcombo in combos:
            kind, suffix = tcombo
            xnames, sub = fpy_interface_sub(fpy, dtype, kind, suffix)
            modprocs.extend(xnames)
            subtext.append(sub)
        subtext.append("\n")

    #Next, chunk the names of the module procedures into blocks of five
    #so that they display nicely for human readability.
    from fortpy.printing.formatting import present_params
    splice = static.replace(interface, present_params(modprocs, 21))
    return splice.replace(interface.replace("py", "xpy"), ''.join(subtext))

def pysave(D, dtype, kind, suffix=None):
    """Generates the fortran code for the pysave interface module procedure
    of the specified type and dimensionality.
    
    :arg D: the integer dimensionality of the variable that is having its
      value set from the file.
    :arg dtype: the fortran data type that the subroutine will be reading data
      for. Can be ["real", "integer", "logical", "complex", "character"].
    :arg kind: the precision/kind information for the type. Can be one of
      ["dp", "sp"] for double/single precision real, integer or complex types.
    :arg suffix: not actually used by the routine; in the signature for compatibility
      with the automatic generation of the code.
    """
    common = {
        "dtype": dtype,
        "D": D,
        "Dx": "({})".format(','.join([":"]*D)) if D >0 else "",
        "suffix": xsuffix[suffix],
        "skind": kind if kind is not None else ""
    }

    #For character, we use the auto-length arrays.
    if dtype in prec and kind in prec[dtype]:
        common["fmt"] = prec[dtype][kind]
    else:
        common["fmt"] = "*"

    if kind is not None and dtype in kmap:
        common["kind"] = "({})".format(kmap[dtype][kind])
    else:
        common["kind"] = ""

    if D > 0:
        inds = ', '.join(["i{}".format(i+1) for i in range(D-1)])
        if D > 1:
            inds = ', ' + inds
        common["vars"] = "integer :: dims({}){}\n".format(D, inds)
    else:
        common["vars"] = ""

    if D > 0:
        lines = ["dims = shape(variable)",
                 "write(FMT, *) dims({})".format(D),
                 "if (dims({}) .eq. 0) return".format(D)]
        common["initial"] = '\n    '.join(lines)
    else:
        common["initial"] = "write(FMT, *) 1"
        
    if D==0:
        if dtype=="character":
            common["write"] = "write(fileunit, {fmt}) trim(adjustl(variable))".format(**common)
        else:
            common["write"] = "write(fileunit, {fmt}) variable".format(**common)
    else:
        lines = []
        if D > 2:
            #Add the shape declaration at the top of the file.
            lines.append("write(fileunit, \"(A, {}i15)\") \"## \", dims".format(D))
        #Next, we make a do loop for each of the dimensions except the last and then
        #write the lines of the last 2D array for each slice.
        for i in range(D-1):
            lines.append("{0}do i{1}=1, dims({1})".format(spacer(2*i), i+1))
            if i == D-3 and D > 2:
                inds = ', '.join(["i{}".format(j+1) for j in range(D-2)])
                fmtstr = "{}write(fileunit, \"(A, {}i15)\") \"## \", {}, 0, 0"
                lines.append(fmtstr.format(spacer(2*(i+1)), D, inds))

        sinds = ', '.join(["i{}".format(j+1) for j in range(D-1)])
        if D > 1:
            sinds = sinds + ', '
        if dtype=="character":
            fmtstr = "{0}call char_write_trimmed(variable({2}:))"
        else:
            fmtstr = "{}write(fileunit, {}) variable({}:)"
        lines.append(fmtstr.format(spacer(2*(D-1)), common["fmt"], sinds))

        #Finally, add the 'end do' statements to close the loops. 
        for i in range(D-2,-1,-1):
            lines.append("{0}end do".format(spacer(2*i), i+1))

        common["write"] = '\n    '.join(lines)

    if dtype=="character":
        common["xname"] = "pysave_{dtype}_{D}d{suffix}".format(**common)
    else:
        common["xname"] = "pysave_{dtype}{skind}_{D}d{suffix}".format(**common)
    
    template = """  subroutine {xname}(variable, filename)
    character(len=*), intent(in) :: filename
    {dtype}{kind}, intent(in) :: variable{Dx}
    character(20) :: FMT
    {vars}
    {initial}

    call file_open(filename, len(filename), '{dtype}')
    {write}
    call file_close()
  end subroutine {xname}
    """
    return (common["xname"], template.format(**common))

interfaces = {
    "__fpy_read__": {
        "integer": [(None, None), ("sp", None), ("dp", None)],
        "real": [("sp", None), ("dp", None)],
        "complex": [("sp", None), ("dp", None)],
        "logical": [(None, None)],
        "character": [("string", None)]
    },
    "__fpy_read_f__": {
        "integer": [(None, "_f"), ("sp", "_f"), ("dp", "_f")],
        "real": [("sp", "_f"), ("dp", "_f")],
        "complex": [("sp", "_f"), ("dp", "_f")],
        "logical": [(None, "_f")],
        "character": [("string", "_f")]
    },
    "__fpy_read_p__": {
        "integer": [(None, "_p"), ("sp", "_p"), ("dp", "_p")],
        "real": [("sp", "_p"), ("dp", "_p")],
        "complex": [("sp", "_p"), ("dp", "_p")],
        "logical": [(None, "_p")],
        "character": [("string", "_p")]
    },
    "__pysave__": {
        "integer": [(None, None), ("sp", None), ("dp", None)],
        "real": [("sp", None), ("dp", None)],
        "complex": [("sp", None), ("dp", None)],
        "logical": [(None, None)],
        "character": [("string", None)]
    }
}
"""Dict tracking which subroutines to generate for specific data types,
kinds and allocatable, pointer or fixed dimension types.
"""
generators = {
    "__fpy_read__": fpy_read,
    "__fpy_read_f__": fpy_read,
    "__fpy_read_p__": fpy_read,
    "__pysave__": pysave
}
"""Dict indicating which python function generates the fortran subroutines
for the given interface field.
"""

from os import path
import sys
sys.path.insert(0, path.expanduser("~/codes/fortpy-dist/"))
import fortpy
from fortpy.utility import get_fortpy_templates_dir

templates = get_fortpy_templates_dir()
statpath = path.join(templates, "static.f90")
with open(statpath) as f:
    static = f.read()

for iface, idict in list(interfaces.items()):
    if iface in generators:
        static = fpy_interface(generators[iface], static, iface, idict)
    else:
        raise ValueError("No generator for field {}".format(iface))

static = static.replace("__version__", fortpy.__version__)
fortpath = path.join(templates, "fortpy.f90")
with open(fortpath, 'w') as f:
    f.write(static)
