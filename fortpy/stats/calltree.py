"""Analyzes executables to find out which other executables anywhere
in the loaded libraries call that one. This is repeated recursively
to generate possible trees that end up in calling the one at the end.
"""
#NOTE: this was originally intended to augment the setting of breakpoints
#in the interactive debugging; I realized that only the explicitly
#specified pre-reqs ever get called first, so that we only need list those.
def _exec_callers(xinst, result):
    """Adds the dependency calls from the specified executable instance
    to the results dictionary.
    """
    for depkey, depval in xinst.dependencies.items():
        if depval.target is not None:
            if depval.target.name in result:
                if xinst not in result[depval.target.name]:
                    result[depval.target.name].append(xinst)
            else:
                result[depval.target.name] = [xinst]

    for xname, xvalue in xinst.executables:
        _exec_callers(xvalue, result)

def _module_callers(parser, modname, result):
    """Adds any calls to executables contained in the specified module.
    """
    if modname in result:
        #We have already processed this module.
        return
    
    module = parser.get(modname)
    mresult = {}
    if module is not None:
        for xname, xinst in module.executables():
            _exec_callers(xinst, mresult)
        result[modname] = mresult
        
        for depkey in module.dependencies:
            depmod = depkey.split('.')[0].lower()
            _module_callers(parser, depmod, result)

def tree(parser, startmod):
    """Returns the call tree for all modules in the library that are
    linked to the specified starting module.
    """
    result = {}
    _module_callers(parser, startmod, result)
    return result

def _call_fan(branch, calls, executable):
    """Appends a list of callees to the branch for each parent
    in the call list that calls this executable.
    """
    #Since we don't keep track of the specific logic in the executables
    #it is possible that we could get a infinite recursion of executables
    #that keep calling each other.
    if executable in branch:
        return
    branch.append(executable)
    
    if executable.name in calls:
        for caller in calls[executable.name]:
            twig = []
            _call_fan(twig, calls, caller)
            branch

def callers(parser, executable):
    """Returns a list of 'module.executable' that call the specified
    executable (i.e. is part of their dependency list).
    """
    calls = tree(parser, executable.module)
    return calls
    # if executable.name in calls:
    #     stack = calls[executable.name]
    #     result = []
    #     while len(stack) > 0:
    #         branch = [executable]
    #         caller = stack.pop()
    #         branch.append(caller)
    #         if caller.name in calls:                
