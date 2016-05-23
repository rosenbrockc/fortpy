#!/usr/bin/env python
from termcolor import cprint
import argparse

def print_result(testkey, percent, time, common):
    """Prints the specified result to the terminal with coloring based
    on how successful it was.
    """
    printkey = testkey.replace("|", " | ")
    if percent > .99:
        color = "green"
    elif percent > .50:
        color = "yellow"
    else:
        color = "red"

    text = "RESULT: {0} \n\t{1:.2%} success ({3:.2%} common) in {2} ms\n"
    cprint(text.format(printkey, percent, time, common), color)

def _get_compilers():
    """Returns a list of compilers from the command-line arguments."""
    complist = []
    if args["compiler"] and args["compiler"] != "*":
        complist.append(args["compiler"])
    elif args["compiler"] == "*":
        from fortpy.testing.compilers import compilers
        for c in compilers:
            complist.append(c)
    else:
        complist.append("gfortran")
    return complist
    
def do_testing(args):
    """Runs the unit tests for all the modules in the code directory."""
    from fortpy.msg import set_quiet
    from fortpy.testing.tester import UnitTester
    set_quiet(args["quiet"])
    
    t = UnitTester(args["stagedir"], args["verbose"], args["templates"], args["fortpy"],
                   args["rerun"], debug=(not args["nodebug"]), profile=args["profile"],
                   strict=args["strict"], quiet=args["quiet"])

    #We only have to write the testing folder once; it gets copied for all
    #the remaining tests that need to be run for different compilers.
    t.writeall(args["codedir"])

    complist = _get_compilers()
    totalperc = 0
    totaltest = 0
    for c in complist:
        result = t.runall(c)
        print("")

        for idk in result:
            totaltime = result[idk].totaltime
            if totaltime is not None:
                timestr = "{0:.4f}".format(totaltime*1000)
            else:
                timestr = "<untimed>"
            print_result(idk, result[idk].percent, timestr, result[idk].common)
            totalperc += result[idk].percent
            totaltest += 1

    #This section for exit codes helps the continuous integration server to know what's
    #going on with all the test results.

    if totaltest == 0 and totalperc == 0:
        _exit_code(0, "No Tests")
    else:
        score = totalperc/totaltest
        
    if score == 1.:
        _exit_code(0, "Success")
    elif score < 1.:
        _exit_code(2, "Failure")
    else:
        _exit_code(3, "Didn't Run")

def _exit_code(i, prefix):
    """Informs the user of the exit code and then exits."""
    print("{1}: exiting with code {0:d}".format(i, prefix))
    exit(i)
        
def _get_parser(codedir, modules=None):
    """Gets a CodeParser instance with specified modules loaded from file or cache."""
    from fortpy.code import CodeParser
    parser = CodeParser()
    parser.verbose = args["verbose"]
    
    if modules is None:
        files = {}
        parser.scan_path(codedir, files)
        for f in files:
            filepath = files[f]
            parser.parse(filepath, True, True)
    else:
        for modname in modules:
            parser.load_dependency(modname, True, True)
        
    return parser
            
def do_auxiliary():
    """Generates an auxiliary f90 file with an interface to save all the user-derived
    type variables referenced anywhere in the code.
    """
    from os import path
    from fortpy.testing.auxiliary import generate
    codedir = path.abspath(path.expanduser(args["codedir"]))
    parser = _get_parser(codedir)

    complist = _get_compilers()
    for c in complist:
        generate(parser, codedir, args["stagedir"], c, debug=(not args["nodebug"]),
                 profile=args["profile"], strict=args["strict"], docompile=args["compileaux"])
    
#Create a parser so that the script can receive arguments
parser = argparse.ArgumentParser(description="Fortpy Automated Unit Testing Tool")

#Add arguments to decide which of the systems and penalties to process.
parser.add_argument("codedir", help="Specify the path to the directory of code files to run tests for.")
parser.add_argument("-stagedir", help="Sets the directory in which to stage the unit tests.")
parser.add_argument("-templates", help="Specify the path to the folder that houses the XML templates.")
parser.add_argument("-outfile", help="Specify a path to save the comparison reports to.")
parser.add_argument("-verbose", help="Sets whether the comparison output is verbose.", type=int, default=1)
parser.add_argument("-mode", help="Sets the strictness of the comparison.", default="default")
parser.add_argument("-fortpy", help="The path to the fortpy templates directory.")
parser.add_argument("-rerun",
                    help=("Specifies a filter for *module* names whose unit tests will be rerun. "
                          "When a test is rerun, it is recompiled and tested, even if the code base "
                          "has not changed since the last test. Value '*' reruns the unit "
                          "tests of *all* modules in the code directory."))
parser.add_argument("-compiler", help="Specify the compiler to use for the unit testing")
parser.add_argument("-pypath", help="Specify a path to add to sys.path before running the tests.")
parser.add_argument("-nodebug", 
                    help=("Compile the executables with DEBUG=false; the default behavior is to "
                          "always compile with DEBUG on to check bounds and overflow errors etc."), 
                    action="store_true")
parser.add_argument("-profile", help="Compile and link with profiling enabled. Analyze profile.",
                    action="store_true")
parser.add_argument("-quiet", help="Run in quiet mode; only essential output appears on stdout.",
                    action="store_true")
parser.add_argument("-auxiliary", action="store_true",
                    help=("Generate an auxiliary f90 module with interfaces to save user-derived "
                          "type variables."))
parser.add_argument("-strict", action="store_true",
                    help="Enable all warnings for the compilers.")
parser.add_argument("-compileaux", action="store_true",
                    help=("Also compile the fpy_auxiliary.f90 into .o, .mod and .so library. "
                          "Requires -stagedir to be specified."))

if __name__ == "__main__":
    #Parse the args from the commandline that ran the script, call initialize
    args = vars(parser.parse_args())

    #We added this argument for debugging installations. That way we can make changes
    #without doing a pip install each time; just put the path to the repo root in 'pypath'
    if args["pypath"]:
        import sys
        sys.path.append(args["pypath"])

    import fortpy
    if args["pypath"]:
        #Change to unit-testing mode so that we don't compete with the live cache.
        from fortpy import settings
        settings.use_test_cache = True

    from fortpy import msg
    msg.set_verbosity(args["verbose"])
        
    testing = not (args["auxiliary"])
    if testing:
        do_testing(args)
    elif args["auxiliary"]:
        do_auxiliary()
