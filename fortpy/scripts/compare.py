#!/usr/bin/env python
import argparse
from fortpy.testing.comparer import FileComparer
from fortpy.testing.results import print_compare_result

def initialize():
    if args["templates"]:
        comparer = FileComparer(args["templates"])
    else:
        comparer = FileComparer()

    if not args["xmltemplate"]:
        args["xmltemplate"] = "{}.xml".format(source.split("/")[-1])
    result = comparer.compare(args["source"], args["target"], args["xmltemplate"], args["mode"])

    if args["save"]:
        with open(args["save"], "w") as f:
            f.write(print_compare_result(result, args["verbose"]))
    else:
        print(print_compare_result(result, args["verbose"]))

#Create a parser so that the script can receive arguments
parser = argparse.ArgumentParser(description="Fortpy File Comparison Tool")

#Add arguments to decide which of the systems and penalties to process.
parser.add_argument("-source", help="Specify the path to the source (templated) file.", required=True)
parser.add_argument("-target", help="Specify the path to the file to compare to source.", required=True)
parser.add_argument("-templates", help="Specify the path to the folder that houses the XML templates.")
parser.add_argument("-save", help="Specify a path to save the comparison report to.")
parser.add_argument("-verbose", help="Sets whether the comparison output is verbose.", action="store_true")
parser.add_argument("-mode", help="Sets the strictness of the comparison.", default="default")
parser.add_argument("-xmltemplate", help="Specify the name of the XML template to use.")

#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_args())
initialize()
