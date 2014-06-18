#!/apps/python/2.7.3/bin/python
from fortpy.testing.tester import UnitTester
import argparse

def initialize():
    t = UnitTester(args["stagedir"], args["verbose"], args["templates"], args["fortpy"],
                   args["rerun"])
    t.writeall(args["codedir"])
    result = t.runall()

    for idk in result:
        print("RESULT: {0} = {1:.2%} successful".format(idk, result[idk].percent))

#Create a parser so that the script can receive arguments
parser = argparse.ArgumentParser(description="Fortpy File Comparison Tool")

#Add arguments to decide which of the systems and penalties to process.
parser.add_argument("codedir", help="Specify the path to the directory of code files to run tests for.")
parser.add_argument("-stagedir", help="Sets the directory in which to stage the unit tests.", required=True)
parser.add_argument("-templates", help="Specify the path to the folder that houses the XML templates.")
parser.add_argument("-outfile", help="Specify a path to save the comparison reports to.")
parser.add_argument("-verbose", help="Sets whether the comparison output is verbose.", action="store_true")
parser.add_argument("-mode", help="Sets the strictness of the comparison.", default="default")
parser.add_argument("-fortpy", help="The path to the fortpy templates directory.", 
                    default="~/pythonpkg/fortpy/templates")
parser.add_argument("-rerun", help="If specified, the tests are re-run with minimal compilation.",
                    action="store_true")

#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_args())
initialize()
