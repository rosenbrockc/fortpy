"""Module for interacting with the compiler settings for running multiple
unit tests and handling input/output file name modifications.
"""
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
    pgfort = Popen("{} --version".format(xpath),
                    shell=True, executable="/bin/bash", stdout=PIPE, stderr=PIPE)
    waitpid(pgfort.pid, 0)
    #Redirect the output and errors so that we don't pollute stdout.
    output = pgfort.stdout.readlines()
    error = pgfort.stderr.readlines()
    if len(output) == 0:
        raise ValueError("Unable to determine the family of compiler at {}.".format(exepath))

    if "gnu" in output[0].lower():
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

_load_xml()
