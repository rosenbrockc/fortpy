#!/usr/bin/env python
import argparse
from fortpy.interop.converter import InputConverter

def initialize():
    converter = InputConverter(args["templates"])
    if args["target"]:
        #We can only convert a single file, take the first one in the list
        converter.convert(args["source"][0], args["version"], args["target"])
    else:
        for path in args["source"]:
            converter.convert(path, args["version"])

#Create a parser so that the script can receive arguments
parser = argparse.ArgumentParser(description="Fortpy File Conversion Tool")

#Add arguments to decide which of the systems and penalties to process.
parser.add_argument("source", help="Specify the path(s) to the source XML input file.", nargs="+")
parser.add_argument("-target", help="Specify the path to the file to convert.")
parser.add_argument("-version", help="Specify the target version of the converted input files.", default="1",
                    type=int)
parser.add_argument("-templates", help="Specify the path to the folder that houses the XML templates.", required=True)
parser.add_argument("-verbose", help="Sets whether the comparison output is verbose.", action="store_true")

#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_args())
initialize()

