#!/usr/bin/env python
import argparse
import sys
from os import path

def parse(args):
    if not args["reparse"]:
        settings.use_filesystem_cache = False
    settings.unit_testing_mode = True

    c = CodeParser()
    
    if args["verbose"]:
        c.verbose = True

    if args["reparse"]:
        fname = args["source"].split("/")[-1].lower()
        c._modulefiles[fname] = []
        c._programfiles[fname] = []
        c._parse_from_file(args["source"], fname, False, False, False)
    elif args["scan"]:
        files = {}
        c.scan_path(args["source"], files)
        
    else:
        c.parse(args["source"])

    if not args["scan"]:
        #Since this is for unit testing, we will access the "private" variables.
        for fname in c._modulefiles:
            for moduledat in c._modulefiles[fname]:
                if args["verbose"] > 2:
                    print(c.modules[moduledat])
                else:
                    print(moduledat)
                
            # for progdat in c._programfiles[fname]:
            #     if args["verbose"]:
            #         print c.programs[progdat]
            #     else:
            #         print progdat

        return c
    else:
        from fortpy import msg
        msg.info("Valid code files surviving .fpyignore:")
        for fpath in files:
            msg.okay(fpath)
        return files

if __name__ == "__main__":    
    #Create a parser so that the script can receive arguments
    parser = argparse.ArgumentParser(description="Fortpy File Parsing Unit Testing Tool")

    #Add arguments to decide which of the systems and penalties to process.
    parser.add_argument("source", help="Specify the path to the source file/folder to parse.")
    parser.add_argument("-verbose", help="Sets whether the comparison output is verbose.", type=int, default=0)
    parser.add_argument("-reparse", help="Overwrite the cached version of the module.", action="store_true")
    parser.add_argument("-pypath", help="Specify a path to add to sys.path before running the tests.")
    parser.add_argument("-scan", action="store_true",
                        help="Scan the specified directory and show which code files would be run.")

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

    if args["verbose"]:
        msg.set_verbosity(args["verbose"])
        
    cparser = parse(args)
