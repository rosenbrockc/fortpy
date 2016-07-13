"""Generates a f90 module file for a code directory with interfaces to save
any user-derived type variables (or arrays of those variables) using the
auto-class functionality of the unit tests.
"""
from fortpy import msg

warnings = {}
"""Keeps track of the user types that have already produced warnings."""

def _append_member(member, result):
    """Appends the element to the result list *only* if it can actually be
    handled by the auxiliary module.
    """
    global warnings
        
    if "private contents" in member.customtype.modifiers:
        if member.kind.lower() not in warnings:
            msg.warn("User type {} skipped ".format(member.kind) +
                     "because members are private.", 0)
            warnings[member.kind.lower()] = 1
        else:
            warnings[member.kind.lower()] += 1
        return

    if member.customtype.name.lower() not in member.customtype.module.publics:
        if member.kind.lower() not in warnings:
            msg.warn("Skipping user type {} ".format(member.kind) +
                     "because it is not marked as public in its module.")
            warnings[member.kind.lower()] = 1
        else:
            warnings[member.kind.lower()] += 1
        return

    result.append(member)

def _el_members(element, result):
    """In connection with _find_members() searches a single code element for
    user-derived type variable declarations.
    """
    for memkey, member in element.members.items():
        if member.customtype is not None and not any([r.matched(member) for r in result]):
            _append_member(member, result)
    if hasattr(element, "parameters"):
        for memkey, member in element.parameters.items():
            if member.customtype is not None and not any([r.matched(member) for r in result]):
                _append_member(member, result)

    #Handle modules that have executables and the nested executables
    #inside of subroutines.
    if hasattr(element, "executables") and type(element).__name__ != "CustomType":
        for nested in element.executables.values():
            _el_members(nested, result)
    
def _find_members(parser, coderoot):
    """Searches the list of modules in the CodeParser instance to find all
    user-derived variable declarations. Filters the list to find the unique
    ones for which saving subroutines should be generated.
    """
    result = []
    for modname, module in parser.modules.items():
        if (modname not in ["fortpy", "fpy_auxiliary"] and
            coderoot.lower() in module.filepath.lower()):
            _el_members(module, result)
            for utype in module.types.values():
                _el_members(utype, result)

    return sorted(result, key=(lambda m: m.kind))

def _print_member(member):
    """Prints the essential details of a member that is having its saving subroutine 
    generated in the auxiliary f90 file.
    """
    if member.D > 0:
        dimension = "({})".format(member.dimension)
    else:
        dimension = "[scalar]"
    if "allocatable" in member.modifiers:
        modifiers = ", allocatable"
    elif "pointer" in member.modifiers:
        modifiers = ", pointer"
    else:
        modifiers = ""
    return "{}{} :: {}".format(member.strtype, modifiers, dimension)

def _get_suffix(member):
    """Returns the fpy_read suffix for the specified member.
    """
    if "pointer" in member.modifiers:
        suffix = "_p"
    elif "allocatable" not in member.modifiers and member.D > 0:
        suffix = "_f"
    else:
        suffix = ""

    return suffix
        
def _get_rsubroutine_recursive(classer, common):
    """Gets the recursive code subroutine for saving a variable with members
    of the same type as itself.
    """
    variable = classer.variable
    lines = []
    varnames = []
    for memname, member in variable.customtype.members.items():
        if member.is_custom:
            if member.kind.lower() == member.customtype.name.lower() and member.D==0:
                xr = "call fpy_read(prefix//'-{0}', '#', variable_{0}_, fpy_success)"
                lines.append(xr.format(member.name))
                lines.append("if (fpy_success .and. variable_{}_ .gt. 0) then".format(memname))
                lines.append("  variable%{0} => pointer_stack(variable_{0}_)".format(memname))
                lines.append("end if")
                varnames.append("variable_{}_".format(memname))
            else:
                xr = "call auxread_{varname}{D}d{suffix}(variable%{memname}, prefix//'-{memname}', pointer_stack, stack)"
                lines.append(xr.format(**{"varname": member.customtype.name,
                                          "D": member.D,
                                          "memname": member.name,
                                          "suffix": "_p" if "pointer" in member.modifiers else ""}))
        else:
            lines.append("call fpy_read{1}(prefix//'-{0}', '#', variable%{0})".format(memname, _get_suffix(member)))
    common["read"] = '\n    '.join(lines)
    common["vars"] = "integer :: {}".format(', '.join(varnames)) if len(varnames) > 0 else ""
    common["init"] = '\n    '.join(["{} = 0".format(v) for v in varnames])

    template = """  recursive subroutine {xname}_(var_iloc, stack, pointer_stack, prefix)
    integer, intent(in) :: var_iloc
    type(fpy_address), allocatable, intent(inout) :: stack(:)
    {dtype}{kind}, pointer, intent(inout) :: pointer_stack(:)
    character(len=*), intent(in) :: prefix

    logical :: fpy_success
    {dtype}{kind}, pointer :: variable
    {vars}
    if (fpy_verbose > 0) write (*,*) "Start recursive read for {xname}_ with '"//prefix//"'."
    {init}

    variable => pointer_stack(var_iloc)
    {read}
    if (fpy_verbose > 0) write (*,*) "End recursive read for {xname}_ with '"//prefix//"'."
  end subroutine {xname}_

  subroutine {xname}(variable, folder, multi_stack, rstack, nested_)
    {dtype}{kind}, pointer, intent(inout) :: variable
    character(len=*), intent(in) :: folder
    type(fpy_address), optional, allocatable, intent(in) :: rstack(:)
    {dtype}{kind}, optional, pointer, intent(in) :: multi_stack(:)
    logical, optional, intent(in) :: nested_

    type(fpy_address), allocatable :: stack(:)
    {dtype}{kind}, pointer :: pointer_stack(:)
    character(len=:), allocatable :: lfolder
    integer :: iloc
    logical :: nested

    if (fpy_verbose > 0) write (*,*) "Start recursive read for {xname} in '"//folder//"'."
    if (present(nested_)) then
      nested = nested_
    else
      nested = .false.
    end if

    if ((present(multi_stack) .and. present(rstack)) .or. nested) then
       lfolder = folder
       if (present(multi_stack) .and. present(rstack)) then
         stack = rstack
         iloc = fpy_address_index(stack, lfolder)
         pointer_stack => multi_stack
       end if
    else
       lfolder = folder//'_'
       call fpy_read_address(stack, lfolder)
       iloc = 1
       allocate(pointer_stack(size(stack)))
    end if

    call {xname}_(iloc, stack, pointer_stack, lfolder)
    variable = pointer_stack(iloc)
    if (fpy_verbose > 0) write (*,*) "End recursive read for {xname} in '"//lfolder//"'."
  end subroutine {xname}

"""
    return ((common["xname"], ), template.format(**common))

def _get_rsubroutine(classer):
    """Gets the code to generate a subroutine for the single member specified.
    """
    member = classer.variable
    common = {
        "varname": member.customtype.name,
        "D": member.D,
        "Dx": "({})".format(','.join([":"]*member.D)) if member.D >0 else "",
        "dtype": member.dtype if member.dtype.lower() != "class" else "type",
        "kind": "({})".format(member.kind),
        "hierarchy": member.customtype.full_name
    }
    common["xname"] = "auxread_{varname}{D}d".format(**common)
    return _get_rsubroutine_nested(classer, common)
    
def _get_rsubroutine_nested(classer, common):
    """Gets the code for by calling the 0d cases in a nested loop."""
    member = classer.variable
    if member.D == 0:
        if member.customtype.recursive:
            return _get_rsubroutine_recursive(classer, common)
        else:
            return _get_rsubroutine_flat(classer, common)
    else:
        vsplice = ', '.join(["fpy{}".format(i) for i in range(member.D)])
        vlines = ["integer :: {}".format(vsplice)]
        vlines.append("character(100) :: pslist")
        vlines.append("type(fpy_address), allocatable :: stack(:)")
        vlines.append("integer, allocatable :: drange(:)")
        vlines.append("{dtype}{kind}, pointer :: tempvar".format(**common))
        common["vars"] = '\n    '.join(vlines)
        
        lines = []
        spacing = "    "

        maxvlen = max([len(v) for v in member.customtype.fixedvar])
        fstr = "'{{0: <{}}}'".format(maxvlen)
        fixedvars = [fstr.format(v) for v in member.customtype.fixedvar]
        sfixedvar = "(/ {} /)".format(', '.join(fixedvars))
        lines.append("{}call fpy_get_drange(drange, folder//'-', {})".format(spacing, sfixedvar))
        asplice = ', '.join(["drange({})".format(i+1) for i in range(member.D)])
        lines.append("{}if (.not. allocated(drange)) return".format(spacing))
        lines.append("{}allocate(variable({}))".format(spacing, asplice))
        
        for i in range(member.D):
            lines.append("{}do fpy{}=1, size(variable, {})".format(spacing, i, i+1))
            spacing += "  "

        nname = "auxread_{varname}0d".format(**common)
        lines.append("{}call fpy_period_join_indices(pslist".format(spacing) +
                     ", (/ {} /), {})".format(vsplice, member.D))
        lines.append("{}tempvar => variable({})".format(spacing, vsplice))
        xr = "{}call {}(tempvar, folder//'-'//trim(adjustl(pslist)), pointer_stack, stack, .true.)"
        lines.append(xr.format(spacing, nname))

        for i in range(member.D, 0, -1):
            spacing = spacing[:-2]
            lines.append("{}end do".format(spacing))
        common["read"] = '\n'.join(lines)
        
        template = """  subroutine {xname}_p(variable, folder, multi_stack, rstack, nested_)
    character(len=*), intent(in) :: folder
    {dtype}{kind}, pointer, intent(inout) :: variable{Dx}
    {dtype}{kind}, pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), optional, intent(in) :: rstack(:)
    logical, optional, intent(in) :: nested_

    {dtype}{kind}, pointer :: pointer_stack(:)
    {vars}
    logical :: nested

    if (fpy_verbose > 0) write (*,*) "Start nested read for {xname}_p in '"//folder//"'."
    if (present(nested_)) then
      nested = nested_
    else
      nested = .false.
    end if

    if (present(multi_stack) .and. present(rstack)) then
      stack = rstack
      pointer_stack => multi_stack
    else
      call fpy_read_address(stack, folder)
      allocate(pointer_stack(size(stack)))
    end if

{read}
    if (fpy_verbose > 0) write (*,*) "End nested read for {xname}_p in '"//folder//"'."
  end subroutine {xname}_p

  subroutine {xname}(variable, folder, multi_stack, rstack, nested_)
    character(len=*), intent(in) :: folder
    {dtype}{kind}, allocatable, target, intent(inout) :: variable{Dx}
    {dtype}{kind}, pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), optional, intent(in) :: rstack(:)
    logical, optional, intent(in) :: nested_

    {dtype}{kind}, pointer :: lvar{Dx}
    logical :: nested

    if (fpy_verbose > 0) write (*,*) "Start nested read for {xname} in '"//folder//"'."
    if (present(nested_)) then
      nested = nested_
    else
      nested = .false.
    end if

    if (allocated(variable)) then
      lvar => variable
    else
      lvar => null()
    end if
    call {xname}_p(lvar, folder, multi_stack, rstack, nested)
    variable = lvar
    if (fpy_verbose > 0) write (*,*) "End nested read for {xname} in '"//folder//"'."
  end subroutine {xname}
"""
        return ((common["xname"], common["xname"]+"_p"), template.format(**common))
    
def _get_rsubroutine_flat(classer, common):
    """Gets the *non* recursive code subroutine for saving a variable with members
    of the same type as itself.
    """
    variable = classer.variable
    lines = []
    for memname, member in variable.customtype.members.items():
        if member.is_custom:
            if member.customtype is None:
                msg.err("Unable to find Type instance for : {}".format(member.definition()))
                continue
            if "allocatable" in member.modifiers:
                d0 = ', '.join(['0' for i in range(member.D)])
                lines.append("allocate(variable%{}({}))".format(memname, d0))
            xr = "call auxread_{varname}{D}d{suffix}(variable%{memname}, lfolder//'-{memname}', nested_=.true.)"
            lines.append(xr.format(**{"varname": member.customtype.name,
                                      "D": member.D,
                                      "memname": member.name,
                                      "suffix": _get_suffix(member)}))
        else:
            lines.append("call fpy_read{1}(lfolder//'-{0}', '#', variable%{0})".format(member.name, _get_suffix(member)))
    common["read"] = '\n    '.join(lines)
    
    template = """  subroutine {xname}(variable, folder, multi_stack, rstack, nested_)
    character(len=*), intent(in) :: folder
    {dtype}{kind}, intent(inout) :: variable{Dx}
    {dtype}{kind}, pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), allocatable, intent(in), optional :: rstack(:)
    logical, optional, intent(in) :: nested_

    character(len=:), allocatable :: lfolder
    logical :: nested

    if (fpy_verbose > 0) write (*,*) "Start flat read for {xname} in '"//folder//"'."

    if (present(nested_)) then
      nested = nested_
    else
      nested = .false.
    end if

    if ((present(multi_stack) .and. present(rstack)) .or. nested) then
       lfolder = folder
    else
       lfolder = folder//'_'
    end if

    {read}
    if (fpy_verbose > 0) write (*,*) "End flat read for {xname} in '"//folder//"'."
  end subroutine {xname}
"""
    return ((common["xname"],), template.format(**common))

def _get_wsubroutine_recursive(classer, common):
    """Gets the recursive code subroutine for saving a variable with members
    of the same type as itself.
    """
    lines = []
    classer.code(lines, "vars", "  ")
    classer.code(lines, "init", "  ")
    common["vars"] = '\n  '.join(lines)

    lines = []
    classer.code(lines, "assign", "  ")
    common["write"] = '\n  '.join(lines)

    template = """  recursive subroutine {xname}_(variable, prefix, stack, wrote)
    {dtype}{kind}, pointer, intent(in) :: variable{Dx}
    type(fpy_address), allocatable, intent(inout) :: stack(:)
    character(len=*), intent(in) :: prefix
    logical, intent(out) :: wrote
    integer :: ploc
    type(fpy_address) :: tadd
    character(len=:), allocatable :: nprefix

  {vars}

    tadd = fpy_get_address(c_loc(variable), prefix)
    ploc = address_loc(stack, tadd)
    if (ploc .eq. -1) then
       call append_address(stack, tadd)
    else
       call pysave(ploc, prefix)
       return
    end if

  {write}
  end subroutine {xname}_

  subroutine {xname}(variable, pfolder, nested_, wrote_)
    character(len=*), intent(in) :: pfolder
    {dtype}{kind}, pointer, intent(in) :: variable{Dx}
    logical, optional, intent(in) :: nested_
    logical, optional, intent(out) :: wrote_
    type(fpy_address), allocatable :: stack_(:)
    
    character(len=:), allocatable :: folder
    logical :: wrote

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('{hierarchy}', folder//'.fpy.type')
    end if

    allocate(stack_(0))
    call {xname}_(variable, folder, stack_, wrote)
    call fpy_save_addresses(stack_, folder)

    if (present(wrote_)) then
      wrote_ = wrote
    end if
  end subroutine {xname}
"""
    return ((common["xname"],), template.format(**common))

def _get_wsubroutine(classer):
    """Gets the code to generate a subroutine for the single member specified.
    """
    member = classer.variable
    common = {
        "varname": member.customtype.name,
        "D": member.D,
        "Dx": "({})".format(','.join([":"]*member.D)) if member.D >0 else "",
        "dtype": member.dtype if member.dtype.lower() != "class" else "type",
        "kind": "({})".format(member.kind),
        "hierarchy": member.customtype.full_name
    }
    common["xname"] = "auxsave_{varname}{D}d".format(**common)
    return _get_wsubroutine_nested(classer, common)
    
def _get_wsubroutine_nested(classer, common):
    """Gets the code for D>0 by calling the 0d cases in a nested loop."""
    member = classer.variable
    if member.D == 0:
        if member.customtype.recursive:
            return _get_wsubroutine_recursive(classer, common)
        else:
            return _get_wsubroutine_flat(classer, common)
    else:
        vsplice = ', '.join(["fpy{}".format(i) for i in range(member.D)])
        vlines = ["integer :: {}".format(vsplice)]
        vlines.append("character(100) :: pslist")
        common["vars"] = '\n    '.join(vlines)
        
        lines = []
        spacing = "    "
        for i in range(member.D):
            lines.append("{}do fpy{}=1, size(variable, {})".format(spacing, i, i+1))
            spacing += "  "

        nname = "auxsave_{varname}0d".format(**common)
        lines.append("{}call fpy_period_join_indices(pslist".format(spacing) +
                     ", (/ {} /), {})".format(vsplice, member.D))
        lines.append("{}call {}(variable({}), folder//'-'//trim(adjustl(pslist)), .true., wrote)".format(spacing, nname, vsplice))

        for i in range(member.D, 0, -1):
            spacing = spacing[:-2]
            lines.append("{}end do".format(spacing))
        common["write"] = '\n'.join(lines)
        
        template = """  subroutine {xname}_p(variable, pfolder, nested_, wrote_)
    character(len=*), intent(in) :: pfolder
    {dtype}{kind}, pointer, intent(in) :: variable{Dx}
    logical, optional, intent(in) :: nested_
    logical, optional, intent(out) :: wrote_
    logical :: wrote
    character(len=:), allocatable :: folder
    {vars}

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('{hierarchy}', folder//'.fpy.type')
    end if

{write}
    if (present(wrote_)) then
      wrote_ = wrote
    end if
  end subroutine {xname}_p

  subroutine {xname}(variable, folder, wrote_)
    character(len=*), intent(in) :: folder
    {dtype}{kind}, allocatable, target, intent(in) :: variable{Dx}
    logical, optional, intent(out) :: wrote_
    {dtype}{kind}, pointer :: lvar{Dx}
    logical :: wrote
    lvar => variable
    call {xname}_p(lvar, folder, wrote)

    if (present(wrote_)) then
      wrote_ = wrote
    end if
  end subroutine {xname}
"""
        return ((common["xname"], common["xname"]+"_p"), template.format(**common))
    
def _get_wsubroutine_flat(classer, common):
    """Gets the *non* recursive code subroutine for saving a variable with members
    of the same type as itself.
    """
    #Just use the AutoClasser to initialize the variables it needs for
    #the save to be successful.
    lines = []
    classer.code(lines, "vars", "  ")
    classer.code(lines, "init", "  ")
    common["vars"] = '\n  '.join(lines)

    lines = []
    classer.code(lines, "assign", "  ")
    common["write"] = '\n  '.join(lines)

    template = """  subroutine {xname}(variable, pfolder, nested_, wrote_)
    character(len=*), intent(in) :: pfolder
    {dtype}{kind}, intent(in) :: variable{Dx}
    logical, optional, intent(out) :: wrote_
    logical, optional, intent(in) :: nested_
    logical :: wrote
    character(len=:), allocatable :: folder
  {vars}

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('{hierarchy}', folder//'.fpy.type')
    end if
    wrote = .false.
  {write}
    if (present(wrote_)) then
      wrote_ = wrote
    end if
  end subroutine {xname}
"""
    return ((common["xname"],), template.format(**common))

def _generate_single(classers, write=True):
    """Generates the subroutines for reading *or* writing using the specified set
    of Auto-classers.
    """
    routines = []
    xnames = []
    modules = []
    modcode = []
    scalars = {}
    
    for c in classers:        
        if write:
            txname, code = _get_wsubroutine(c)
        else:
            txname, code = _get_rsubroutine(c)
            
        addcode = False
        for xname in txname:
            if xname.lower() not in xnames:
                xnames.append(xname.lower())
                addcode = True
        if addcode:
            routines.append(code)
            
        if write and c.variable.customtype.module.name not in modules:
            modules.append(c.variable.customtype.module.name)
            modcode.append("  use {}".format(modules[-1]))
        skey = c.variable.customtype.name.lower()
        #Sometimes we have multi-D variables of a certain type, but no scalar for that type.
        #Since we need the scalar to perform the multi-D saves, we keep track of them.
        if c.variable.D == 0 and skey not in scalars:
            scalars[skey] = c

    for c in classers:
        skey = c.variable.customtype.name.lower()
        if c.variable.D > 0 and skey not in scalars:
            #Make a copy of the variable, set its dimensionality to scalar and then process it.
            from copy import copy
            from fortpy.testing.elements import AutoClasser
            cvar = copy(c.variable)
            ccopy = AutoClasser(cvar, c.folder, c.varfile, c.cases, c.coderoot, True)
            ccopy.variable.D = 0
            ccopy.variable.dimension = None

            if write:
                txname, code = _get_wsubroutine(ccopy)
            else:
                txname, code = _get_rsubroutine(ccopy)

            addcode = False
            for xname in txname:
                if xname.lower() not in xnames:
                    xnames.append(xname.lower())
                    addcode = True
            if addcode:
                routines.append(code)

            scalars[c.variable.customtype.name.lower()] = ccopy

    return (xnames, modcode, routines, modules)

def _prepare_dir(parser, modnames, auxdir):
    """Prepares the directory with the necessary dependencies to compile the auxiliary
    module.

    :arg parser: a fortpy.code.CodeParser for accessing the dependency modules.
    :arg modnames: a list of module names that need to be included in the directory.
    """
    #We need to see whether to include the pre-compiler directive or not.
    precompile = False
    for needed in modnames:
        if parser.modules[needed].precompile:
            precompile = True
            break
        
    modnames.append("fpy_auxiliary")

    from os import path
    from fortpy.interop.make import makefile
    lines = []
    makepath = path.join(auxdir, "Makefile.fpy_aux")
    makefile("fpy_aux", modnames, makepath, "fpy_auxiliary.all", precompile,
             parser=parser, executable="so", makefpyaux=True)

def _compile(auxdir, compiler=None, debug=False, profile=False):
    """Compiles the auxiliary module to generate a .mod and a .o file.

    :arg auxdir: the auxiliary directory to compile in.
    :arg tversion: the template version of fpy_auxiliary.f90 from fortpy. Needed
      to determine if the code-specific module is out of date.
    :arg compiler: the key of the compiler from the compilers.xml or None for default.
    """
    from fortpy.testing.compilers import compile_general
    return compile_general(auxdir, compiler, "fpy_aux", debug, profile, vupdates=["fpy_auxiliary"], quiet=True)

def _should_recompile(auxdir, parser, modules, compiler):
    """Determines whether the fpy_auxiliary module should be rewritten and recompiled.
    """
    from os import path
    from shutil import copy
    
    recompile = False
    for modulename in modules:
        module = parser.modules[modulename]
        auxpath = path.join(auxdir, path.split(module.filepath)[1])
        if not path.isfile(auxpath):
            copy(module.filepath, auxdir)
            recompile = True
        else:
            fmtime = path.getmtime(module.filepath)
            xmtime = path.getmtime(auxpath)
            if xmtime < fmtime:
                recompile = True
                copy(module.filepath, auxdir)

    #Also check the version numbers of the template fpy_auxiliary.f90 and the one present
    #in our directory (if it exists).
    fpyaux = path.join(auxdir, "fpy_auxiliary.f90")
    if path.isfile(fpyaux):
        from fortpy.testing.compilers import template_version, get_fortpy_version
        tversion = template_version(compiler, "fpy_auxiliary.f90")
        xversion = get_fortpy_version(compiler, fpyaux)
        recompile = recompile or (xversion != tversion)
    else:
        recompile = True

    return recompile

def generate(parser, coderoot, stagedir, compiler=None, debug=False, profile=False,
             strict=False, docompile=True):
    """Generates the f90 module file for all members referenced in the parser's
    modules list.
    """
    import fortpy
    from fortpy.utility import get_fortpy_templates_dir
    from fortpy.testing.elements import AutoClasser
    
    members = _find_members(parser, coderoot)
    folder = "filename"
    classers = [AutoClasser(m, folder, m.name, [], coderoot, True) for m in members]
    
    from os import path
    templates = get_fortpy_templates_dir()
    statpath = path.join(templates, "fpy_auxiliary.f90")
    with open(statpath) as f:
        static = f.read()
        
    from fortpy.printing.formatting import present_params
    xnames, modcode, routines, modules = _generate_single(classers)

    from fortpy.code import order_module_dependencies
    modules = order_module_dependencies(modules, parser)

    from os import mkdir, path
    if docompile:
        auxdir = path.join(stagedir, "fortpy.aux")
        if not path.isdir(auxdir):
            mkdir(auxdir)
    else:
        auxdir = coderoot

    if docompile and not _should_recompile(auxdir, parser, modules, compiler):
        msg.okay("The current version of fpy_auxiliary.f90 is up-to-date.")
        from fortpy.testing.compilers import replace
        target = replace(auxdir + ".[c]", compiler)
        return (0, True, target)
    
    static = static.replace("__auxsave__", present_params(xnames, 21))
    static = static.replace("__aux_uses__", '\n'.join(modcode))        
    static = static.replace("__version__", fortpy.__version__)
    static = static.replace("__fxauxsave__", '\n'.join(routines))

    xnames, modcode, routines, dmods = _generate_single(classers, False)
    static = static.replace("__auxread__", present_params(xnames, 21))
    static = static.replace("__fxauxread__", '\n'.join(routines))

    fortpath = path.join(auxdir, "fpy_auxiliary.f90")
    with open(fortpath, 'w') as f:
        f.write(static)

    if docompile:
        _prepare_dir(parser, modules, auxdir)
        return _compile(auxdir, compiler, debug, profile)
