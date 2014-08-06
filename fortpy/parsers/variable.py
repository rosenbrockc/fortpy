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
                           r",?(?P<modifiers>[ \w\t:,()]+)?::\s*(?P<names>[^\n!]+)" #Removed $ from end.
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
            #Unfortunately, the dimension can also be specified as a modifier and
            #the dimensions can include variable names and functions. This introduces
            #the possibility of nested lists.
            if "dimension" in modifiers:
                start, end = self._get_dim_modifier(modifiers)
                dimension = modifiers[start+1:end]
                dimtext = modifiers[modifiers.index("dimension"):end+1]
                modifiers = re.split(",\s*", modifiers.replace(dimtext, "").strip())
                modifiers.append("dimension")
            else:
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

    def _get_dim_modifier(self, modifiers):
        """Extracts the dimension information from the string of modifiers extracted by
        the regex.
        
        :arg modifiers: the list of modifiers identified by the regex.
        """
        suffix = modifiers.split("dimension")[1]
        start = modifiers.index("dimension") + len("dimension")
        #We use a stack to monitor how many parenthesis we have traversed in the string.
        #Once we reach the closing parenthesis, we know that we have the dimension info.
        stack = []
        args = []
        for i in range(len(suffix)):
            if suffix[i] == '(':
                stack.append(i + start)
            elif suffix[i] == ')':
                args.append((stack.pop(), i + start))

        #The last entry in args should be the indices of the entire dimension expression
        return args[-1]
