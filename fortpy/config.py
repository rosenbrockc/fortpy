#This config module handles retrieval of global variable values for fortpy from environment
#variables or a config XML file.
import types
from sys import modules

class _config(types.ModuleType):
    def __init__(self):
        """A configuration class to store global variables under a configuration
        module name. Exposes values as properties to allow arbitrary variable
        storage via XML as well as static oft-used variables."""
        self._vardict = {}
        self._initialized = False
        self.getenvar("FORTPY_CONFIG")

        if self.implicit_XML is not None:
#            print "Loading FORTPY config variables from " + self.implicit_XML
            self.load_xml(self.implicit_XML)
            self._initialized = True      

    @property
    def implicit_XML(self):
        """Returns the path to the XML file that stores the config info."""
        return self.property_get("FORTPY_CONFIG")

    @property    
    def codes(self):
        """Returns a list of additional code folders to look in for files."""
        return self.property_get("codes")

    @property
    def mappings(self):
        """Returns a list of module mappings for modules that aren't named the
        same as the files they are in."""
        return self.property_get("mappings")

    def property_get(self, key):
        if key in self._vardict:
            return self._vardict[key]
        else:
            return None

    def getenvar(self, envar):
        from os import getenv
        """Retrieves the value of an environment variable if it exists."""
        if getenv(envar) is not None:
            self._vardict[envar] = getenv(envar)

    def load_xml(self, filepath):
        """Loads the values of the configuration variables from an XML path."""
        from os import path
        import xml.etree.ElementTree as ET
        #Make sure the file exists and then import it as XML and read the values out.
        uxpath = path.expanduser(filepath)
        if path.isfile(uxpath):
            tree = ET.parse(uxpath)
            root = tree.getroot()

            for child in root:
                if child.tag == "codes":
                    self._load_codes(child)
                elif child.tag == "mappings":
                    self._load_mapping(child)
      
    def _load_codes(self, tag):
        """Extracts all the paths to additional code directories to 
        be considered.

        :arg tag: the ET tag for the <codes> element."""
        codes = []
        for code in tag:
            if code.tag == "trunk":
                codes.append(code.attrib["value"])

        self._vardict["codes"] = codes

    def _load_mapping(self, tag):
        """Extracts all the alternate module name mappings to 
        be considered.

        :arg tag: the ET tag for the <mappings> element."""
        mappings = {}
        for mapping in tag:
            if mapping.tag == "map":
                mappings[mapping.attrib["module"]] = mapping.attrib["file"]

        self._vardict["mappings"] = mappings

#modules["config"] = _config()
