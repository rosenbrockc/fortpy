"""Methods for gathering statistics on the unit test coverage of a library.
"""
def testability(parser, module):
    """Returns a dictionary of testability scores, indexed by executable
    name with value {"score", "dscore", "pscore"}.
    """
    execs = _list_executables(parser, module)
    result = {}
    for xname, xinst in execs.items():
        pscore, dscore, total = _testability(xinst)
        if pscore is None:
            continue
        result[xinst.full_name] = {"score": total, "pscore": pscore, "dscore": dscore}
    return result   

def _testability(xinst, xall=None):
    """Returns a testability score for the specified executable instance.
    This takes into account the depth of the dependency chains, number of
    parameters and the test coverage of the dependency calls.

    :arg xall: when True, scores are even returned for executables that have
      unit tests defined; otherwise, only non-tested ones are scored.
    """
    analysis = _analyze(xinst)
    if (xall is None and ((analysis["ntests"] > 0 and analysis["ncases"] > 0)
        or analysis["skip"])):
        return (None, None, None)
    
    #Only the input parameters are relevant for the testability score.
    inpars = [n for n, p in analysis["params"].items() if "in" in p["direction"]]
    #For the dependencies, we are interested in how many there are and
    #which of those have good coverage with unit tests already.
    covered = []
    uncovered = []
    
    for kdep, vdepl in xinst.dependencies.items():
        for vdep in vdepl:
            if isinstance(vdep, str):
                continue
            if vdep.target is not None:
                danal = _analyze(vdep.target)
                if (danal["ntests"] > 0 and danal["ncases"] > 0):
                    #This dependency has unit tests defined for it already, so
                    #the hierarchy is preserved if the calling method is tested.
                    covered.append(vdep)
                else:
                    iparams = [a for (a, p) in zip(vdep.argnames, vdep.target.ordered_parameters)
                               if "in" in p.direction]
                    uncovered.extend([p for p in iparams if p not in uncovered and p not in inpars])

    #Now we assign scores to the number of input parameters and dependencies.
    from math import exp
    #For dependencies, the covered ones offer no *penalty* to running the tests
    #whereas the uncovered ones do. Also, if the covered ones use any of the input
    #parameters that the calling method uses, there is a high likelihood that the
    #data for those parameters already exists. For uncovered dependencies, we only
    #care about them if they have input parameters.
    dscore = -(1. - exp(-len(uncovered)/5.))
    overlap = []
    for cdep in covered:
        overlap.extend([a for a in cdep.argnames if a in inpars and a not in overlap])

    #If there are no input parameters, then this score will be 1.0; as more
    #input parameters are added, it becomes more difficult to test.
    pscore = exp(-(len(inpars) - len(overlap))/5.)

    return pscore, dscore, pscore + dscore

def _list_executables(parser, module):
    """Returns the dictionary of executables for the specified module *name*
    in the code parser.
    """
    mkey = module.lower()
    if mkey not in parser.modules:
        parser.load_dependency(mkey, True, True)

    if mkey in parser.modules:
        return parser.modules[mkey].executables
    else:
        raise ValueError("The code parser cannot locate module '{}'".format(mkey) +
                         " for analyzing unit testing statistics.")

def summary(parser, module):
    """Returns a dictionary and text summarizing the coverage of the
    specified module.

    :arg parser: the fortpy.code.CodeParser instance.
    :arg module: the name of the module to return a summary for.
    """
    execs = _list_executables(parser, module)
    result = {}
    for xname, xinst in execs.items():
        analysis = _analyze(xinst)
        result[xinst.full_name] = (analysis, _describe(xinst, analysis))
    return result

def _analyze(anexec):
    """Analyzes the specified executable instance to get unit testing
    coverage statistics.
    """
    #Examine the test_group instance to see if we have any testable
    #instances. If we do, take a look at the cases to determine coverage.
    #Setup a defaults dictionary with the stats fields we may find
    #information about.
    result = {
        "ntests": 0,
        "tests": {},
        "ncases": 0,
        "summary": False,
        "params": {},
        "nparams": 0,
        "skip": False
    }
    """Dictionary of results for analysis of the unit test coverage for
    this executable. Keys:
      ntests: number of unit tests defined for the executable.
      tests: dictionary keyed by unit test identifier; describes the details
        of coverage for each unit test (cases, etc.).
      ncases: the total number of test cases across all tests.
      summary: true if the docstring summary exists.
      params: dictionary of summary, regularity and constraint data for each
        parameter in the executable.
      nparams: the number of parameters with summary tags.
    """
    if anexec.test_group is not None:
        tg = anexec.test_group
        for tname, tinst in tg.tests.items():
            result["tests"][tname] = {
                #The 1 for no cases listed is for the default case that exists whenever
                #a test is set to run checks on model outputs
                "cases": 1 if tinst.cases is None else len(tinst.cases),
                "inputs": len(tinst.inputs),
                "outputs": len(tinst.outputs),
                "assignments": len(tinst.methods),
                "globals": len(tinst.variables),
                "targets": len(tinst.targets),
                "enabled": tinst.runchecks and tinst.execute,
                "executes": tinst.execute,
                "timed": tinst.timed
            }
            result["ncases"] += result["tests"][tname]["cases"] if tinst.runchecks else 0
            result["ntests"] += 1 if result["tests"][tname]["enabled"] else 0

    #<skip enabled="true">
    result["skip"] = any([d.xml.tag == "skip" and "enabled" in d.xml.attrib
                          and d.xml.attrib["enabled"].lower() == "true"
                          for d in anexec.docstring])
    result["summary"] = anexec.summary != "No summary for element."
    result["nparams"] = 0
    result["nparams.in"] = 0
    for pinst in anexec.ordered_parameters:
        hassum = pinst.summary != "No summary for element.",
        result["params"][pinst.name] = {
            "summary": hassum,
            "regular": any(["regular" in d.xml.attrib and d.xml.attrib["regular"].lower() == "true"
                            for d in pinst.docstring]),
            "bounded": any(["range" in d.xml.attrib for d in pinst.docstring]),
            "direction": pinst.direction
        }
        if hassum:
            result["nparams"] += 1
        if "in" in pinst.direction:
            result["nparams.in"] += 1

    pgood = result["nparams"] == len(anexec.parameters)
    if result["summary"]:
        if pgood:
            docs = "GOOD"
        else:
            docs = "OK"
    else:
        docs = "NO"
    result["docsum"] = docs
            
    return result 

fmtstr = "{0:<40s} | {1:^5d} | {2:^9.2f} | {3:^4s} |"
headstr = "{0:<40s} | {1:^5s} | {2:^9s} | {3:^4s} |"
dheader = headstr.format(*("Executable Identifier", "Tests", "Case/Test", "Docs"))
def _describe(anexec, analysis):
    """Returns a text description of the analysis for the unit tests of the
    specified executable instance.
    
    :arg analysis: the result from calling _analyze(anexec).
    """
    #The text description just has the mass totals for everything.
    #Since it will be tabulated, we just return a tuple of the values
    #that need to be passed to the format string.
    pcase = float(analysis["ncases"]) / (analysis["ntests"] if analysis["ntests"] > 0 else 1)
    tup = (anexec.full_name, analysis["ntests"], pcase, analysis["docsum"])
    return fmtstr.format(*tup)
