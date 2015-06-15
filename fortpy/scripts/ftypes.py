#!/usr/bin/env python
import argparse
import sys
from os import path

def _parse():
    """Parses the specified Fortran source file from which the wrappers will
    be constructed for ctypes.
    """
    if not args["reparse"]:
        settings.use_filesystem_cache = False

    c = CodeParser()    
    if args["verbose"]:
        c.verbose = True

    if args["reparse"]:
        c.reparse(args["source"])
    else:
        c.parse(args["source"])

    return c

def writef90():
    """Writes the F90 module file from the specified code parser for *all* the modules
    inside of it.
    """
    cparser = _parse()
    dependencies = ["{}_c".format(m) for m in cparser.modules]
    originals = list(cparser.modules.keys())
    writers = {}
    for modname, module in list(cparser.modules.items()):
        msg.info("Writing wrapper module for {}".format(modname))
        fwriter = f90.WrapperModule(module, args["library"], args["f90"], args["link"])
        writers[modname] = fwriter
        fwriter.write_f90()
        msg.okay("Finished writing {}_c.f90".format(modname))
        
    #We can compile all these modules together into a single shared library.
    writer = list(writers.values())[0]
    code = writer.make(remake=args["remake"], dependencies=dependencies, compiler=args["compiler"],
                       debug=args["debug"], profile=args["profile"])
    if code == 0:
        for modname in originals:
            msg.info("Writing python wrapper module {}.py".format(modname))
            writers[modname].write_py()
            msg.okay("Finished writing {}.py".format(modname))
        
parser = argparse.ArgumentParser(description="Fortpy ctypes Interoperability Tool")

parser.add_argument("source", help="Specify the path to the source file to parse.")
parser.add_argument("-verbose", help="Sets whether the comparison output is verbose.", action="store_true")
parser.add_argument("-reparse", help="Overwrite the cached version of the module.", action="store_true")
parser.add_argument("-pypath", help="Specify a path to add to sys.path before running the tests.")
parser.add_argument("-f90", help="Specify the directory to save the Fortran modules to.")
parser.add_argument("-python", help="Specify the directory to save the Python modules to.")
parser.add_argument("-library", help="Specify the name of the library to create with ftypes.", required=True)
parser.add_argument("-link", help=("Specify the name of a library with the compiled, original code that "
                                   "can be used to create the wrapper shared library."))
parser.add_argument("-remake", help=("When writing and compiling the wrapper library, force a re-compile "
                                     "even if the wrapper module already exists."), action="store_true")
parser.add_argument("-compiler", default="gfortran", help="Specify the compiler to use on the wrapper module.")
parser.add_argument("-debug", help="Compile the wrapper shared library with debug flags.",
                    action="store_true")
parser.add_argument("-profile", help="Compiler the wrapper shared library with gprof profiling enabled.",
                    action="store_true")
#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_args())

#We added this argument for debugging installations. That way we can make changes
#without doing a pip install each time; just put the path to the repo root in 'pypath'
if args["pypath"]:
    import sys
    sys.path.append(args["pypath"])

import fortpy
from fortpy.code import CodeParser
from fortpy import settings
from fortpy import msg
from fortpy.interop import ftypes as f90

writef90()
