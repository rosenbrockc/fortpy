import re
from ..elements import CustomType, TypeExecutable

class TypeParser(object):
    """Extracts custom type definitions from fortran code files."""
    def __init__(self, vparser, docparser):
        self.vparser = vparser
        self.docparser = docparser
        self.setup_regex()

    def setup_regex(self):
        """Sets up the patterns and compiled regex objects for parsing types."""
        #Regex for matching the entire body of the type and getting top-level modifiers.
        self._RX_TYPE = r"\s+type(?P<modifiers>,\s+(public|private))?(\s+::)?\s+(?P<name>[A-Za-z0-9_]+)" + \
                        r"(?P<contents>.+?)end\s+type(\s+(?P=name))?"
        self.RE_TYPE = re.compile(self._RX_TYPE, re.DOTALL | re.I)
        #Regex for finding if the type is private
        self._RX_PRIV = "private.+?contains"
        self.RE_PRIV = re.compile(self._RX_PRIV, re.DOTALL | re.I)
        #Regex for finding methods buried in a type declaration.
        self._RX_EXEC = r"^\s*(?P<modifiers>[^:]+)\s+::\s+(?P<name>[A-Za-z0-9_]+)" + \
                        r"(\s+=>\s+(?P<points>[A-Za-z0-9_]+))?$"
        self.RE_EXEC = re.compile(self._RX_EXEC, re.M | re.I)
        #Regex for getting text after contains statement
        self._RX_CONTAINS = "\n\s*contains(?P<remainder>.+)"
        self.RE_CONTAINS = re.compile(self._RX_CONTAINS, re.DOTALL | re.I)
        
    def parse(self, module):
        """Extracts all the types from the specified module body."""
        matches = self.RE_TYPE.finditer(module.contents)
        result = {}
        for match in matches:
            name = match.group("name")
            modifiers = match.group("modifiers")
            if modifiers is not None:
                cleanmods = re.split("[\s,]+", modifiers.strip())
            else:
                cleanmods = []

            contents = match.group("contents")
            result[name] = self._process_type(name, cleanmods, contents, module, match)
            if "public" in result[name].modifiers:
                module.publics[name] = 1
        module.types = result

    def _process_type(self, name, modifiers, contents, module, match):
        """Processes a regex match of a type's contents."""
        #First, we need to see if the types children are private.
        if self.RE_PRIV.search(contents):
            modifiers.append("private contents")
        
        #Next, we need to parse out all the members of the type and their docstrings
        members = self.vparser.parse(contents, None)
        #Last of all, parse out all the executables including the finalizer
        execs = self._process_execs(contents, module.name)
        
        #Now we can create the type code element and handle the member docstrings
        t = CustomType(name, modifiers, members, execs, module)

        #Set the regex start and end char indices
        t.start, t.end = module.absolute_charindex(match.string, match.start(),
                                                             match.end())

        #Update the parent for embedded members and executables
        for key in t.members.keys():
            t.members[key].parent = t
        for key in t.executables.keys():
            t.executables[key].parent = t

        #We need to look in the parent module docstrings for this types decorating tags.
        key = "{}.{}".format(module.name, name)
        if key in module.predocs:
            t.docstring = self.docparser.to_doc(module.predocs[key], name)

        #Extract the docstrings from the type body and associate them with their members
        memdocs = self.docparser.parse_docs(contents, name)
        if name in memdocs:
            docs = self.docparser.to_doc(memdocs[name], name)
            self.docparser.process_memberdocs(docs, t)

        return t
        
    def _process_execs(self, contents, modulename):
        """Extracts all the executable methods that belong to the type."""
        #We only want to look at text after the contains statement
        match = self.RE_CONTAINS.search(contents)
        result = {}
        #It is possible for the type to not have any executables
        if match is not None:
            exectext = match.group("remainder")
            if exectext:
                for mexec in self.RE_EXEC.finditer(exectext):
                    name = mexec.group("name")
                    modifiers = mexec.group("modifiers")
                    pointsto = mexec.group("points")
                    if pointsto is not None:
                        pointsto = "{}.{}".format(modulename, pointsto)
                    e = TypeExecutable(name, modifiers, None, pointsto)
                    result[name] = e
                    #Executables at the type level don't have any documentation because each
                    #one points to a module-level executable that will be documented.
        return result
