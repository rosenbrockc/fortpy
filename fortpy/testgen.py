"""This module provides a generic method for saving input and output files to use
with the unit testing framework.
"""
def _get_filename(absdir, name, infile=True, case=None):
    """Gets the file name for the specified variable name and case.

    :arg absdir: the expanded folder path that will contain the file.
    :arg name: the name of the variable to use.
    :arg infile: if true, '.in.case' will be used; otherwise '.out.case'.
    :arg case: the identifier for the variable set that the variable belongs to.
    """
    from os import path
    if infile == True:
        suffix = "{}.in".format(name.lower())
    elif infile==False:
        suffix = "{}.out".format(name.lower())
    else:
        suffix = name.lower()

    if case is None:
        return path.join(absdir, suffix)
    else:
        return "{}.{}".format(suffix, case)

def print_value(value):
    """Determines the type of the specified value and formats it with high precision for
    reading in by fortran (or for later comparison by python).
    """
    if isinstance(value, int):
        return "%15d" % value
    elif isinstance(value, float):
        return "%.09e" % value
    elif isinstance(value, complex):
        raise ValueError("Complex data types not supported yet for auto-test generation "
                         "in fortpy.")
    elif isinstance(value, list) or isinstance(value, tuple):
        return '   '.join(map(print_value, value))

def write_generic(value, filepath=None):
    """Writes the specified variable value to the file, formatting it for use with Fortran
    and taking its dimensionality into account.

    :arg value: the variable value to save into 'filepath'.
    :arg filepath: the full, expanded path to save the variable to.
    """
    #We basically just test for the dimensionality of each variable first and then
    #call the relevant routine for it.
    if isinstance(value, list) or isinstance(value, tuple):
        #Check if it is 1D or 2D.
        if len(value) > 0 and (isinstance(value[0], list) or isinstance(value[0], tuple)):
            contents = '\n'.join(map(print_value, value))
        else:
            contents = print_value(value)
    else:
        contents = print_value(value)

    if filepath is not None:
        with open(filepath, 'w') as f:
            f.write(contents)
    else:
        return contents

def load(names, folder, infile=True, case=None):
    """Loads the value of a previously saved variable from file using the fortpy test
    generation naming standards.

    :arg names: a list of variable names whose values should be loaded from file.
    :arg folder: the path to the folder to load files from.
    :arg infile: when true, the suffix '.in.case' is assumed; otherwise '.out.case'.
    :arg case: the identifier for the set of variables to load.
    """
    from os import path
    absdir = path.expanduser(folder)    
    result = []

    for name in names:
        absfile = path.join(absdir, _get_filename(absdir, name, infile, case))
        varval = []
        if not path.isfile:
            result.append(None)
        else:
            with open(absfile) as f:
                for line in f:
                    values = list(map(eval, line.strip().split()))
                    if len(values) == 1:
                        varval.append(values[0])
                    else:
                        varval.append(values)

            if len(varval) == 1:
                result.append(varval[0])
            else:
                result.append(varval)

    return result

def save(values, names, folder, infile=True, case=None, overwrite=False):
    """Saves the list of variable values to separate files so that they can be used
    by the automated unit testing framework in fortpy.

    :arg values: A list of variables whose values should be saved.
    :arg names: The name prefix to give each file, one for each variable in 'values'.
    :arg folder: the path to the folder to save the files in.
    :arg infile: Specifies whether the file should have ".in.{case}" or ".out.{case}" appended
      to its file name.
    :arg case: the identifier to use for this *set* of variable values.
    :arg overwrite: unless True, if a file already exists with that name, the variable's value
      will not be overwritten.
    """
    from os import path
    absdir = path.expanduser(folder)
    if not path.isdir(absdir):
        from os import makedirs
        makedirs(absdir)

    for val, name in zip(values, names):
        absfile = path.join(absdir, _get_filename(absdir, name, infile, case))
        if overwrite or not path.isfile(absfile):
            write_generic(val, absfile)
