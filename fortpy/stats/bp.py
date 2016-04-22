"""Methods for testing a code library against Fortran best practices to
help uncover subtle bugs that took a while for us to track down. See
especially http://www.cs.rpi.edu/~szymansk/OOF90/bugs.html"""

def _exec_check_pointers(executable):
    """Checks the specified executable for the pointer condition that not
    all members of the derived type have had their values set.

    Returns (list of offending members, parameter name).
    """
    oparams = []
    pmembers = {}
    xassigns = map(lambda x: x.lower().strip(), executable.external_assignments())

    def add_offense(pname, member):
        """Adds the specified member as an offender under the specified parameter."""
        if pname not in oparams:
            oparams.append(pname)
        if pname not in pmembers:
            pmembers[pname] = [member]
        else:
            pmembers[pname].append(member)

    def check_buried(executable, pname, member):
        """Checks whether the member has its value changed by one of the dependency
        subroutines in the executable.
        """
        for d in executable.dependencies:
            if pname in d.argnames:
                pindex = d.argnames.index(pname)
                dtarget = d.target
                if dtarget is not None:
                    mparam = dtarget.ordered_parameters[pindex]
            
    for pname, param in executable.parameters.items():
        if param.direction == "(out)" and param.is_custom:
            utype = param.customtype
            if utype is None:
                continue
            for mname, member in utype.members.items():
                key = "{}%{}".format(pname, mname).lower().strip()
                if key not in xassigns:
                    #We also need to check the dependency calls to other, buried subroutines.
                    compname = "{}%{}".format(pname, mname).lower()
                    if executable.changed(compname) is None:
                        add_offense(pname, member)
                    
    return (oparams, pmembers)

def _type_check_pointers(utype):
    """Checks the user-derived type for non-nullified pointer array declarations
    in its base definition.

    Returns (list of offending members).
    """
    result = []
    for mname, member in utype.members.items():
        if ("pointer" in member.modifiers and member.D > 0 and
            (member.default is None or "null" not in member.default)):
            result.append(member)

    return result
    
def check_pointers(parser, codedir=None, mfilter=None, recursive=False):
    """Checks the modules in the specified code parser to see if they
    have common, but subtle, pointer bugs in:

    1. subroutines with a parameter of intent(out) and user-derived type
      must* set *all* members of that parameter or they will have an
      *undefined* status.
    2. pointer-type arrays that are not nullified are set to a valid target
      will return 'T' when passed to `associated`. Best practice is to nullify
      pointer arrays in user-derived types as the default value on those types.
    
    :arg parser: [fortpy.code.CodeParser] with the modules to search *already loaded*.
    :arg codedir: specify the full path to the library whose modules should be searched,
      just another way to filter which modules are generating the warnings.
    :arg mfilter: filter to apply to module names; can use the wildcard standard
      from bash.
    """
    from fnmatch import fnmatch
    from fortpy.msg import std, set_verbosity, info

    set_verbosity(0)
    W1 = "   {} '{}' does not set the value of members '{}' in parameter '{}'."
    W2 = "   Type '{}' does not nullify members '{}' on creation."

    offenders = {}
    for (modname, module) in parser.modules.items():
        if not recursive and codedir is not None and not codedir.lower() in module.filepath.lower():
            continue
        if mfilter is not None and not fnmatch(module.name.lower(), mfilter.lower()):
            continue

        #Test the first condition above for all subroutines in the module; also handle
        #the recursively defined subroutines.
        hprinted = False
        for xname, xvalue in module.executables.items():
            oparams, pmembers = _exec_check_pointers(xvalue)
            if len(oparams) > 0:
                if not hprinted:
                    info("Best practice suggestions: {}".format(module.filepath))
                    hprinted = True
                    
                for oparam in oparams:
                    plist = ', '.join([p.name for p in pmembers[oparam]])                
                    std(W1.format(type(xvalue).__name__, xname, plist, oparam), 0)
                offenders[xvalue.full_name] = (oparams, pmembers)

        for tname, tvalue in module.types.items():
            result = _type_check_pointers(tvalue)
            if len(result) > 0:
                if not hprinted:
                    info("Best practice suggestions: {}".format(module.filepath))
                    hprinted = True

                plist = ', '.join([p.name for p in result])
                std(W2.format(tname, plist), 0)            
                offenders[xvalue.full_name] = result

    return offenders
