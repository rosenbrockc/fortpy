import os
import sys
import argparse

sys.path.append("/Users/trunks/codes/fortpy-dist")
import fortpy
import fortpy.isense

def get_path(partial):
    """Returns the full path to a file to use during the testing.

    :arg partial: the partial path relative to the working directory.
    """
    cwd = os.getcwd()
    return os.path.join(cwd, partial)

def get_script(edit):
    """Gets an isense.Script object for the specified edit.

    :arg edit: a list of [partial path, line, char] for the script.
    """
    fullpath = get_path(edit[0])
    with open(fullpath) as f:
        edited = f.read()
    return fortpy.isense.Script(edited, edit[1], edit[2], fullpath)

def get_completion_tests(key= "monte"):
    """Returns a list of functions to test and cursor locations to test
    at for the specified script object.
    """
    scripts = {
        "monte": [["in_function_call", 189, 30],["completions", 121, 12],
                  ["bracket_complete", 133, 25], ["in_function_call", 141, 9],
                  ["bracket_complete", 122, 12], ["bracket_complete", 191, 9],
                  ["bracket_complete", 311, 16], ["in_function_call", 312, 29],
                  ["completions", 313, 11], ["in_function_call", 314, 19],
                  ["in_function_call", 315, 13], ["completions", 316, 14],
                  ["goto_definitions", 103, 24]],
        "deriv": [["completions", 207, 4], ["bracket_complete", 208, 4],
                  ["completions", 209, 8], ["in_function_call", 330, 22],
                  ["completions", 345, 8], ["in_function_call", 52, 27],
                  ["goto_definitions", 52, 27]],
        "cefit": [["in_function_call", 3087, 30], ["in_function_call", 264, 17]]
        }
    if key in scripts:
        return scripts[key]
    else:
        return []

source = ["tests/isense/monte_orderparam.f90", 212, 9]

def test_original(key="monte"):
    """Tests the standard completions *excluding* the real-time update
    functionality for the source file. Returns a 
    """
    paths = {
        "monte": "tests/isense/original/monte_orderparam.f90",
        "deriv": "tests/isense/original/derivative_structure_generator.f90",
        "cefit": "tests/isense/original/ce_fit.f90"
    }
    fullpath = get_path(paths[key])
    
    with open(fullpath) as f:
        original = f.read()
    settings.real_time_update = False
    settings.unit_testing_mode = True

    results = []
    for test in get_completion_tests(key):
        results.append(test)
        try:
            script = fortpy.isense.Script(original, test[1], test[2], fullpath)
            result = getattr(script, test[0])()
            results.append(result)        
        except:
            print([test, fullpath, key])
            raise

        #Also examine the description for in function calls and completions.
        if test[0] in ["in_function_call", "completions"] and len(result) > 0:
            if type(result) == list:
                results.append(result[0].description)
                results.append(result[0].fulldoc)

        if test[0] == "goto_definitions":
            results.append(result.fulldoc)

    return results

import fortpy.settings as settings
settings.real_time_update = False

#Create a parser so that the script can receive arguments
parser = argparse.ArgumentParser(description="Fortpy Automated Unit Testing Tool")

#Add arguments to decide which of the systems and penalties to process.
parser.add_argument("-isense", action="store_true",
                    help=("Run the tests for the isense auto-completion; otherwise run tests"
                          " on the unit testing framework."))

#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_known_args()[0])

#if args["isense"]:
settings.real_time_update = False
keys = ["monte", "deriv", "cefit"]
for key in keys:
    with open(get_path("tests/isense/compl_{}.tmp".format(key)), 'w') as f:
        f.write("\n\n".join([str(i) for i in test_original(key)]))

def test_real_time():
    """Tests all the real-time update machinery for the isense package."""
    #Get the unmodified original file to make changes to.
    original = get_script(source)

    #Technically, we should be able to run each of these edits in
    #turn, which each make a single adjustment.
    results = []
    for edit in edits:
        single = []
        try:
            single.append(edit[0:3])
            if edit[3] is not None and edit[4] in edit[3]:
                single.append(getattr(edit[3][edit[4]], edit[5]))
            elif edit[3] is not None:
                single.append("'{0}' not found in {1}".format(edit[4], 
                                                            list(edit[3].keys())))
            script = get_script(edit)
            result = script.in_function_call()
            if edit[3] is not None:
                single.append(getattr(edit[3][edit[4]], edit[5]))
            single.append(result)

        except:
            print(single)
            raise
        results.append("\n".join([str(i) for i in single]))

    return results

import fortpy.isense.cache as cache
cache.parser().reparse(get_path(source[0]))
module = cache.parser().modules["monte_orderparam"]
edits = [["tests/isense/edits/signature.f90", 208, 30, 
          module.executables, "single_corr", "paramorder"], 
         ["tests/isense/edits/docstring.f90", 208, 30, None],
         ["tests/isense/edits/variable.f90", 206, 9, 
          module.executables, "flip_after", "members"], 
         ["tests/isense/edits/new_element.f90", 142, 4, 
          module.executables, "new_method", "paramorder" ]]

settings.real_time_update = True
with open(get_path("tests/isense/rtupdate.tmp"), 'w') as f:
    f.write('\n\n')
    f.write("\n\n".join([str(i) for i in test_real_time()]))
