import re
from fortpy.elements import Interface

class InterfaceParser(object):
    """Parses generic, assignment and operator interfaces in Fortran modules and programs."""
    def __init__(self, docparser):
        self.docparser = docparser
        self.setup_regex()

    def setup_regex(self):
        """Regex definitions for parsing the code elements."""
        self._RX_INTERFACE = (r"\n\s*interface\s+(?P<name>[a-z0-9_]+)(\s\((?P<symbol>[.\w+*=/-]+)\))?"
                              r"(?P<contents>.+?)"
                              r"end\s*interface\s+(?P=name)?")
        self.RE_INTERFACE = re.compile(self._RX_INTERFACE, re.I | re.DOTALL)
        
    def parse(self, parent):
        """Parses any interfaces out of the module preamble which is the pre-contains section
        of the program or module.
        
        :arg parent: the ancle.elements.Module instance that will own the interfaces that get
          parsed from the module.
        """
        for iface in self.RE_INTERFACE.finditer(parent.preamble):
            name = iface.group("name")
            symbol = iface.group("symbol")
            contents = iface.group("contents")
            result = Interface(name, [], parent, symbol)
            self._parse_procedures(contents, result)

            #Set the regex start and end char indices
            result.start, result.end = parent.absolute_charindex(iface.string, iface.start(),
                                                                 iface.end())
            self.update_docs(result, parent)
            parent.interfaces[name.lower()] = result

    def update_docs(self, iface, module):
        """Updates the documentation for the specified interface using the module predocs."""        
        #We need to look in the parent module docstrings for this types decorating tags.
        key = "{}.{}".format(module.name, iface.name)
        if key in module.predocs:
            iface.docstring = self.docparser.to_doc(module.predocs[key][0], iface.name)
            iface.docstart, iface.docend = (module.predocs[key][1], module.predocs[key][2])

    def _parse_procedures(self, contents, iface):
        """"Parses the list of procedures or signatures defined in the generic, operator
        or assignment interface.
        
        :arg contents: the text between the interface ... end interface keywords.
        :arg iface: the fortpy.elements.Interface instance that will own the parsed
          procedures/signatures.
        """
        procs = contents.split("module procedure")
        stripped = [p.strip() for p in procs if p.strip() != ""]
        for embed in stripped:
            #We want to first extract the name of the procedure, then find its element
            #instance to append to the Interface instance's procedures list.
            methods = re.split(",\s*", embed.replace("&\n", ""))
            keys = ["{}.{}".format(iface.module.name, m) for m in methods]
            iface.procedures.extend(keys)
