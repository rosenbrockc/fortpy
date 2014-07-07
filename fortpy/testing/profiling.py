from .. import msg
from os import path, system
#What I learned. gprof is a tool that does exactly what we want and would
#work with both ifort and gfortran (both using the -pg option for compiling
#and for linking. This is supposed to produce a gmon.out file in the directory
#where the executable ran. gprof analyzes the output file and prints nice
#tables with function calls and can even limit the output to just the function
#we are interested in.

#The problem is that Apple's proprietary junk makes it not work. It works on
#linux, but not Mac OSX. The first problem seems to be that the executables
#don't even produce a 'gmon.out' file to analyze. The second is that we don't
#have gprof. The internet says that Apple wants you to use their instruments
#application; but that requires you to attach to a long-running process, which
#then gets analyzed. We just want what gprof does.
gprof = "gprof {} -F {} -b > profiling.out"

def profile(testpath, identifier, exepath, compiler):
    """Performs the profiling steps for the specified unit test specification
    based on the compiler.

    :arg testpath: the path to the folder in which the test was run.
    :arg identifier: a module.method string identifying the method being tested.
    :arg testid: the identifier for the specific test being run.
    :arg exepath: the full path to the executable that was run.
    """
    if compiler == "ifort":
        function = _after_ifort(identifier)
    else:
        function = _after_gfortran(identifier)

    command = gprof.format(exepath, function)
    fullcmd = "cd {}; {}".format(testpath, command)
    system(fullcmd)

    profpath = path.join(testpath, "profiling.out")
    if path.isfile(profpath):
        msg.okay("PROFILING: successfully executed for {}".format(identifier))

def _after_ifort(identifier):
    """Performs post-execution cleanup of the ifort profiling data in the test
    directory and extracts information related exclusively to the method being
    unit tested.

    :arg identifier: a module.method string identifying the method being tested.
    """
    return identifier.lower().replace(".", "_mp_") + "_"

def _after_gfortran(identifier):
    """Same post-execution cleanup as _after_ifort() but for gfortran using
    the gprof executable.
    
    :arg identifier: a module.method string identifying the method being tested.
    """
    #Pattern __group_theory_MOD_grouper
    return "__" + identifier.lower().replace(".", "_MOD_")
