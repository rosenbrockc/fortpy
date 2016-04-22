#!/usr/bin/env python
import argparse

def parse():
    """Parses all the modules in the library specified by the script args.
    """
    from fortpy.code import CodeParser
    c = CodeParser()
    
    if args["verbose"]:
        c.verbose = True

    f90files = {}
    c.scan_path(args["source"], f90files)
    for fname, fpath in f90files.items():
        if fname not in c._modulefiles:
            c._modulefiles[fname] = []
        c._parse_from_file(fpath, fname, args["recursive"], args["recursive"], False)

    return c

def _check_pointers(parser):
    """Checks the pointer best-practice conditions."""
    from fortpy.stats.bp import check_pointers
    check_pointers(parser, args["source"], args["filter"], args["recursive"])

checks = {
    "pointer": _check_pointers
}

#Create a parser so that the script can receive arguments
parser = argparse.ArgumentParser(description="Fortpy Best Practices Analyzer")

#Add arguments to decide which of the systems and penalties to process.
parser.add_argument("source", help="Specify the path to the source directory.")
parser.add_argument("-filter", help="Only analyze modules that match this filter.")
parser.add_argument("-check", nargs="+", help="Specify the best-practice checks to run.",
                    choices=["pointers", "all"], required=True)
parser.add_argument("-verbose", help="Sets the level of verbosity for comparison output.", type=int,
                    default=0)
parser.add_argument("-pypath", help="Specify a path to add to sys.path before running the tests.")
parser.add_argument("-recursive", help=("Also analyze the libraries that the specified 'source' "
                                        "directory depends on."), action="store_true")

#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_args())

#We added this argument for debugging installations. That way we can make changes
#without doing a pip install each time; just put the path to the repo root in 'pypath'
if args["pypath"]:
    import sys
    sys.path.insert(0, args["pypath"])

import fortpy
parser = parse()
from fortpy.msg import okay

if "all" in args["check"]:
    for ckey in checks:
        okay("Running {} checks".format(ckey))
        checks[ckey](parser)
else:
    for ckey in args["check"]:
        if ckey in checks:
            info("Running {} checks".format(ckey))
            checks[ckey](parser)
