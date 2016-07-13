"""Module for interacting with the compiler settings for running multiple
unit tests and handling input/output file name modifications.
"""
from fortpy import msg

active = None
"""The name of the active compiler being used by fortpy at the moment.
"""
_default = None
"""The name of the default compiler to use if no other is specified.
"""
compilers = None
"""A dict of the compilers that were specified in the 'compilers.xml'
file (if it exists).
"""
def detect_compiler(exepath):
    """Determines the underlying family of the compiler.

    :arg exepath: the full path to the compiler executable.
    """
    #We just run the --version command and see if the first line
    #includes a "gfortran" label. Otherwise, we assume it is ifort.
    from os import waitpid, path
    from subprocess import Popen, PIPE

    #Append any trailing / that we need to get rsync to work correctly.
    xpath = path.abspath(path.expanduser(exepath))
    pgfort = Popen("{} --version".format(xpath), close_fds=True,
                    shell=True, executable="/bin/bash", stdout=PIPE, stderr=PIPE)
    waitpid(pgfort.pid, 0)
    #Redirect the output and errors so that we don't pollute stdout.
    output = pgfort.stdout.readlines()
    error = pgfort.stderr.readlines()
    if len(output) == 0:
        raise ValueError("Unable to determine the family of compiler at {}.".format(exepath))
    pgfort.stdout.close()
    pgfort.stderr.close()
    
    if "gnu" in output[0].decode("UTF-8").lower():
        return "gfortran"
    else:
        return "ifort"

class Compiler(object):
    """Represents a compiler installed on the machine to use when 
    compiling the unit tests to be run.
    """
    def __init__(self, xml):
        """:arg xml: the <compiler> child tag to create the instance from.
        """
        self.path = None
        """The full path to the compiler executable."""
        self.key = None
        """The string identifier to use in file names and execution 
        directory names.
        """
        self.name = None
        """The short name of the compiler that the developer refers to.
        """
        self.default = False
        """Specifies whether this compiler is set to be the default 
        compiler when none is specified explicitly.
        """
        self.family = None
        """The type of the compiler; one of ["ifort", "gfortran"]."""
        self._parse_xml(xml)

    def _parse_xml(self, xml):
        from fortpy.utility import get_attrib
        self.path = get_attrib(xml, "path", "compiler")
        self.key = get_attrib(xml, "key", "compiler")
        self.name = get_attrib(xml, "name", "compiler")
        self.default = get_attrib(xml, "default", cast=bool, default=False)
        #Auto-detect the underlying type of the compiler.
        self.family = detect_compiler(self.path)

def replace(text, compiler, family=False):
    """Replaces all instances of [c] in the specified text with the key of
    the compiler if it exists.

    :arg text: the text to replace [c] or [f] in.
    :arg compiler: the lowered text key of the compiler to use.
    :arg family: when True, the string "[f]" is replaced instead of "[c]"
    """
    srep = "[f]" if family else "[c]"
    sval = ""
    if compiler in compilers:
        if family:
            sval = compilers[compiler].family[0]
        else:
            sval = compilers[compiler].key
    else:
        if _default is not None and _default in compilers:
            if family:
                sval = compilers[_default].family[0]
            else:
                sval = compilers[_default].key
        elif compiler == "gfortran":
            #This is the default that we specified for fortpy, just replace
            #the text with the key 'g'. This is true for both family and
            #compiler replacements.
            sval = 'g'

    return text.replace(srep, sval)

def family(compiler):
    """Returns the compiler family (ifort or gfortran) for the compiler
    with the specified name.
    """
    if compiler in compilers:
        return compilers[compiler].family
    else:
        #Check if a default is specified, then return that.
        if _default is not None and _default in compilers:
            return compilers[_default].family
        else:
            return "gfortran"

def executor(compiler):
    """Returns the full path to the compiler to use when executing actual
    compile commands.
    """
    if compiler in compilers:
        return compilers[compiler].path
    else:
        #Check if a default is specified, then return that.
        if _default is not None and _default in compilers:
            return compilers[_default].path
        else:
            return "gfortran"
    
def _load_xml():
    global compilers
    if compilers is not None:
        return
    
    import fortpy.config
    import sys
    config = sys.modules["config"]
    compilers = {}
    global _default
    
    if config.compilers is not None:
        from os import path
        import xml.etree.ElementTree as ET
        #Make sure the file exists and then import it as XML and read the values out.
        uxpath = path.abspath(path.expanduser(config.compilers))
        if path.isfile(uxpath):
            tree = ET.parse(uxpath)
            root = tree.getroot()
            for child in root:
                if child.tag == "compiler":
                    compiler = Compiler(child)
                    compilers[compiler.name] = compiler
                    if compiler.default:
                        _default = compiler.name

            if _default is None and len(compilers) > 0:
                _default = list(compilers.keys())[0]

_templatev = {}
"""Holds the version number of various files needed for compilation."""
        
def template_version(compiler, filename):
    """Returns the version number of the specified file."""
    global _templatev
    if filename not in _templatev:
        from os import path
        from fortpy.utility import get_fortpy_templates_dir
        tempath = path.join(get_fortpy_templates_dir(), filename)
        _templatev[filename] = get_fortpy_version(compiler, tempath, attribute="codeversion")

    return _templatev[filename]
                
def get_fortpy_version(compiler, fortpath, recursed=False, attribute="version"):
    """Gets the fortpy version number from the first line of the specified file."""
    result = []
    #If the file doesn't exist yet, we don't try to find the version information.
    from os import path
    ext = path.splitext(fortpath)[1]
    cext = "" if compiler is None else replace(".[c]", compiler)
    if not path.isfile(fortpath) or ext in [".o", ".mod", cext]:
        if path.isfile(fortpath + '.v'):
            return get_fortpy_version(compiler, fortpath + '.v', True)
        else:
            return result

    with open(fortpath) as f:
        for line in f:
            try:
                lt = line.index("<")
                vxml = "<doc>{}</doc>".format(line[lt::])
            except ValueError:
                vxml = ""
            break
        
    if "<fortpy" in vxml:
        import xml.etree.ElementTree as ET
        x = list(ET.XML(vxml))
        if len(x) > 0:
            try:
                result = list(map(int, x[0].attrib[attribute].split(".")))
            except ValueError:
                fstr = "can't extract version information from '{}' in '{}'."
                msg.err(fstr.format(x[0].attrib[attribute], fortpath))
            except KeyError:
                fstr = "fortpy tag does not contain attribute '{}': {} ({})"
                msg.err(fstr.format(attribute, x[0].attrib, fortpath))

    if len(result) == 0 and not recursed:
        return get_fortpy_version(compiler, fortpath + '.v', True)
    else:
        return result

def _compile_simple(compiler, modnames, folder):
    """Compiles the specified list of modules (in specified order) to produce a .mod and a .o
    for each of them.
    """
    for modname in modnames:
        msg.info("Compiling {0}.mod and {0}.o for {1}".format(modname, compiler))
    
    from os import waitpid, path
    from subprocess import Popen, PIPE
    codefiles = ' '.join(["{}.f90".format(m) for m in modnames])
    command = "cd {0}; {1} {2}; {1} -c {2}".format(folder, executor(compiler), codefiles)
    pcompile = Popen(command, shell=True, executable="/bin/bash",
                     stdout=PIPE, stderr=PIPE, close_fds=True)
    pcompile.stdout.close()
    pcompile.stderr.close()
    waitpid(pcompile.pid, 0)

def _vupdate_compiled_module(compiler, modname, folder, tversion=None, rename=True):
    """Compiles a .mod and .o for the specified compiler and specified.
    """
    from os import path
    opath = path.join(folder, "{}.o".format(modname))
    mpath = path.join(folder, "{}.mod".format(modname))
    if path.isfile(opath) and path.isfile(mpath):
        if rename:
            from shutil import move
            nopath = path.join(folder, replace("{}.o.[c]".format(modname), compiler))
            nmpath = path.join(folder, replace("{}.mod.[c]".format(modname), compiler))
            move(opath, nopath)
            move(mpath, nmpath)
        else:
            nopath = opath
            nmpath = mpath
            
        #Create the version files so we can keep track of the compiled versions.
        vpaths = [nopath, nmpath]
        for np in vpaths:
            if tversion is None:
                tversion = get_fortpy_version(compiler, np)
            if len(tversion) == 0:
                #Try to get the template version for the module.
                tversion = template_version(compiler, modname + ".f90")

            vp = np + ".v"
            with open(vp, 'w') as f:
                f.write('#<fortpy version="{}" />'.format('.'.join(map(str, tversion))))
    else:
        msg.err("Unable to find {0}.o and {0}.mod.".format(modname))

def _ensure_fileversion(compiler, modname, folder, target, trycompile=True):
    """Makes sure that the module's f90, mod and o files are up to date with the template
    version. If they aren't compile and copy them.

    :arg compiler: the compiler key from compilers.xml
    :arg modname: the name of the module to check file versions for and move (e.g. "fortpy").
    :arg folder: the folder that contains the up-to-date, "template" version of the module.
    :arg target: the folder to copy the compiled files to.
    :arg trycompile: if the codefile has not been compiled yet, or if the version is out of 
      date, should the code try a simple compile?
    """
    from os import path
    codefile = "{}.f90".format(modname)
    compfiles = ["{}.{}".format(modname, ext) for ext in ["o", "mod"]]
    
    tversion = template_version(compiler, codefile)
    for sdfile in compfiles:
        fdfile = replace(sdfile + ".[c]", compiler)
        ftarget = path.join(target, sdfile)
        dversion = get_fortpy_version(compiler, ftarget)

        if not path.isfile(ftarget) or dversion != tversion:
            source = path.join(folder, fdfile)
            sversion = get_fortpy_version(compiler, source)
            if trycompile and (not path.isfile(source) or sversion != tversion):
                _compile_simple(compiler, [modname], folder)
                _vupdate_compiled_module(compiler, modname, folder, tversion)
            elif not path.isfile(source):
                msg.warn("{} does not exist.".format(source))
                continue
            elif sversion != tversion:
                msg.warn("{} has an old version number.".format(source))

            from fortpy.utility import symlink
            symlink(source, ftarget)
            #If the file is a binary, we need to save a .v with version
            #information as well for the next time we want to copy it.
            pre, ext = path.splitext(ftarget)
            if ext in [".o", ".so", ".mod"]:
                with open(ftarget + '.v', 'w') as f:
                    f.write("# <fortpy version=\"{}\" />".format('.'.join(map(str, tversion))))

def compile_general(folder, compiler=None, identifier=None, debug=False, profile=False,
                    quiet=False, moptions=None, inclfortpy=True, vupdates=None, strict=False):
    """Same as for compile() but the folder is assumed to be generic and a copy of it
    is made for the specific compiler that we are dealing with.
    """
    #Because we often run the tests for multiple compiler versions, we need
    #a copy of the execution directory that was setup for the testing.
    from fortpy.utility import copytree
    from os import path
    target = replace(folder + ".[c]", compiler)
    copytree(folder, target)
    
    code, success = compile(target, compiler, identifier, debug, profile, quiet, moptions,
                            inclfortpy, vupdates, strict)
    return (code, success, target)
                    
def compile(folder, compiler=None, identifier=None, debug=False, profile=False,
            quiet=False, moptions=None, inclfortpy=True, vupdates=None,
            strict=False, inclfpyaux=False):
    """Runs the makefile in the specified folder to compile the 'all' rule.

    :arg vupdates: a list of module names for which the output .mod and .o files
      should have version information attached.
    """
    if inclfortpy:
        #Before we can compile the library, we need to make sure that we have a fortpy
        #.mod and .o compiled with the *same* compiler version specified.
        from fortpy.utility import get_fortpy_templates_dir
        _ensure_fileversion(compiler, "fortpy", get_fortpy_templates_dir(), folder)
        
    options = ""
    if debug:
        options += " DEBUG=true"
    if profile:
        options += " GPROF=true"
    if strict:
        options += " STRICT=true"    

    if moptions is not None:
        for opt in moptions:
            options += " {}".format(opt)
        
    if identifier is not None:
        codestr = "cd {}; make -f 'Makefile.{}' F90='{}' FAM='{}'" + options
        command = codestr.format(folder, identifier, executor(compiler), family(compiler))
    else:
        codestr = "cd {}; make F90='{}' FAM='{}'" + options
        command = codestr.format(folder, executor(compiler), family(compiler))
        
    #If we are running in quiet mode, we don't want the compile information
    #to post to stdout; only errors should be redirected. This means we need
    #to wrap the execution in a subprocess and redirect the std* to PIPE
    from os import waitpid, path
    from subprocess import Popen, PIPE
    pcompile = Popen(command, shell=True, executable="/bin/bash", stdout=PIPE, stderr=PIPE,
                     close_fds=True)
    waitpid(pcompile.pid, 0)
    
    if not quiet:
        output = [x.decode('utf8') for x in pcompile.stdout.readlines()]
        msg.std(''.join(output))
    #else: #We don't need to get these lines since we are purposefully redirecting them.
    error = [x.decode('utf8') for x in pcompile.stderr.readlines()]
    code = len(error)
    if code != 0:
        msg.err(''.join(error))

    pcompile.stdout.close()
    pcompile.stderr.close()
    #It turns out that the compiler still returns a code of zero, even if the compile
    #failed because the actual compiler didn't fail; it did its job properly. We need to
    #check for the existence of errors in the 'compile.log' file.
    lcount = 0
    errors = []
    log = path.join(folder, "compile.{}.log".format(identifier if identifier is not None else "default"))
    with open(log) as f:
        for line in f:
            lcount += 1
            if lcount > 21 and lcount < 32:
                errors.append(line)
            elif lcount > 21:
                break

    if len(errors) > 0:
        #There are 21 lines in the compile.log file when everything runs correctly
        #Overwrite code with a bad exit value since we have some other problems.
        code = 1
        #We want to write the first couple of errors to the console and give them the
        #option to still execute if the compile only generated warnings.
        msg.warn("compile generated some errors or warnings:")
        msg.blank()
        msg.info(''.join(errors))

    if vupdates is not None:
        for modname in vupdates:
            _vupdate_compiled_module(compiler, modname, folder, rename=False)
        
    return (code, len(errors)==0)
                
_load_xml()
