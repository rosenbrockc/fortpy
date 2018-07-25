"""Exposes objects and methods for library-level configuration settings.
"""
libs = {}
"""dict: keys are paths to code directories that have configuration settings
specified; values are :class:`Library` instances.
"""
def option(library, name, default=None):
    """Returns the value of the option for a given library.

    Args:
        library (str): full path to the library code root.
        name (str): name of the option/attribute to get.
        default: value to return if the library doesn't define that option.
    """
    result = default
    if library not in libs:
        L = Library(library)
    else:
        L = libs[library]
        
    if hasattr(L, name):
        result = getattr(L, name)
    return result

class Library(object):
    """Represents a set of library-level configuration settings for unit testing
    etc.

    Args:
        libpath (str): path to the code root of the library. The configuration
          file should be called `.fortpy.ini` in that directory.

    Attributes:
        libpath (str): path to the code root of the library. The configuration
          file should be called `.fortpy.ini` in that directory.
    """
    def __init__(self, libpath):
        from os import path
        self.libpath = path.abspath(path.expanduser(libpath))
        self.options = {}

        self._parse()
        
    def __getattr__(self, name):
        if name in self.options:
            return self.options[name]

    def _parse_flags(self, parser, ftype):
        """Parses a set of directory/file flags.
        
        Args:
            ftype (str): one of ['I', 'L'] for `-I` and `-L` flags in compilers and
              linkers. 
        """
        from os import path
        secname = "flags.{}".format(ftype)
        if parser.has_section(secname):
            self.options[ftype] = []
            for ipath, incl in parser.items(secname):
                if incl != "1":
                    continue
                
                apath = path.abspath(path.expanduser(ipath))
                if path.isfile(apath) or path.isdir(apath):
                    self.options[ftype].append(apath)
                
    def _parse(self):
        """Parses the configuration file; sets the library globally.
        """
        #Add the reference to this library configuration to the global
        #settings.
        global libs
        libs[self.libpath] = self

        from six.moves.configparser import ConfigParser
        from os import path
        filepath = path.join(self.libpath, ".fortpy.ini")
        if path.isfile(filepath):
            parser = ConfigParser()
            parser.readfp(open(filepath))

        if parser.has_section("modules.external"):
            self.options["externals"] = parser.options("modules.external")

        self._parse_flags(parser, "I")
        self._parse_flags(parser, "L")

        if parser.has_section("flags.libs"):
            self.options["libs"] = []
            for lname, incl in parser.items("flags.libs"):
                if incl == "1":
                    self.options["libs"].append(lname)
