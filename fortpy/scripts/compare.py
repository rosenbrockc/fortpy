#!/usr/bin/env python
import argparse
def initialize():
    from fortpy.msg import set_verbosity
    set_verbosity(args["verbose"])
    
    comparer = FileComparer(template_folder=args["templates"])
    result = comparer.compare(args["source"], args["target"], args["xmltemplate"], args["mode"])

    if args["save"]:
        from os import path
        fullpath = path.abspath(path.expanduser(args["save"]))
        with open(fullpath, "w") as f:
            f.write(print_compare_result(result, args["verbose"]))
    else:
        print((print_compare_result(result, args["verbose"])))

    return result

#Create a parser so that the script can receive arguments
parser = argparse.ArgumentParser(description="Fortpy File Comparison Tool")

#Add arguments to decide which of the systems and penalties to process.
parser.add_argument("-source", help="Specify the path to the source (templated) file.", required=True)
parser.add_argument("-target", help="Specify the path to the file to compare to source.", required=True)
parser.add_argument("-templates", help="Specify the path to the folder that houses the XML templates.")
parser.add_argument("-save", help="Specify a path to save the comparison report to.")
parser.add_argument("-verbose", help="Sets the level of verbosity for comparison output.", type=int,
                    default=0)
parser.add_argument("-mode", help="Sets the strictness of the comparison.", default="default")
parser.add_argument("-xmltemplate", help="Specify the name of the XML template to use.")
parser.add_argument("-pypath", help="Specify a path to add to sys.path before running the tests.")

#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_args())

#We added this argument for debugging installations. That way we can make changes
#without doing a pip install each time; just put the path to the repo root in 'pypath'
if args["pypath"]:
    import sys
    sys.path.insert(0, args["pypath"])

from fortpy.testing.comparer import FileComparer
from fortpy.testing.results import print_compare_result
result = initialize()
