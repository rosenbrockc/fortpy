import re
from ..elements import ValueElement

class VariableParser(object):
    """Extracts members, locals and parameter definitions from fortran code."""
    def __init__(self):
        self.setup_regex()

    def setup_regex(self):
        "Sets up all the patterns and compiled regexes for extracting variables."""        
        #Regex for finding members of a type
        self._RX_MEMBERS = r"^\s*(?P<type>character|real|type|logical|integer|class)" + \
                           r"(?P<kind>\([A-Za-z0-9_]+\))?" + \
                           r",?(?P<modifiers>[^:]+)?\s*::\s*(?P<names>[^\n!]+)" #Removed $ from end.
        self.RE_MEMBERS = re.compile(self._RX_MEMBERS, re.M | re.I)

        self._RX_MULTIDEF = r"(?P<name>[^(,=]+)(?P<dimension>\([^)]+\))?(\s?=\s*(?P<default>.+))?"
        self.RE_MULTIDEF = re.compile(self._RX_MULTIDEF, re.I)

    def parse(self, string, parent):
        """Parses all the value code elements from the specified string."""
        result = {}
        for member in self.RE_MEMBERS.finditer(string):
            mems = self._process_member(member, parent)
            #The regex match could contain multiple members that were defined
            #on the same line in the code file.
            for onemem in mems:
                result[onemem.name.lower()] = onemem
        return result

    def _process_member(self, member, parent):
        """Extracts all the member info from the regex match; returns a ValueElements."""
        #The modifiers regex is very greedy so we have some cleaning up to do
        #to extract the mods.
        modifiers = member.group("modifiers")
        if modifiers is not None:
            modifiers = re.split("[,\s]+", modifiers.strip())
            if "" in modifiers:
                modifiers.remove("")

        dtype = member.group("type")
        kind = member.group("kind")
        names = member.group("names")
        
        #If there are multiple vars defined on this line we need to return
        #a list of all of them.
        result = []

        #They might have defined multiple vars on the same line
        for vmatch in self.RE_MULTIDEF.finditer(names.strip()):
            name = vmatch.group("name").strip()
            default = vmatch.group("default")
            dimension = vmatch.group("dimension")
            #Now construct the element and set all the values, then add it in the results list.
            result.append(ValueElement(name, modifiers, dtype, kind, default, dimension, parent))

        return result
