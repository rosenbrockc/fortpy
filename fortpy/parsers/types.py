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
        self._RX_TYPE = r"\n\s*type(?P<modifiers>,\s+(public|private))?(\s*::)?\s+(?P<name>[A-Za-z0-9_]+)" + \
                        r"(?P<contents>.+?)end\s*type(\s+(?P=name))?"
        self.RE_TYPE = re.compile(self._RX_TYPE, re.DOTALL | re.I)
        #This regex is the same as RE_TYPE, only the contents are removed from the definition.
        self._RX_SIG =  r"type(?P<modifiers>,\s+(public|private))?(\s+::)?\s+(?P<name>[A-Za-z0-9_]+)"
        self.RE_SIG = re.compile(self._RX_SIG, re.I)

        #Regex for finding if the type is private
        self._RX_PRIV = "private.+?(contains)?"
        self.RE_PRIV = re.compile(self._RX_PRIV, re.DOTALL | re.I)
        #Regex for finding methods buried in a type declaration.
        self._RX_EXEC = r"^\s*(?P<modifiers>[^:]+)\s+::\s+(?P<name>[A-Za-z0-9_]+)" + \
                        r"(\s+=>\s+(?P<points>[A-Za-z0-9_]+))?$"
        self.RE_EXEC = re.compile(self._RX_EXEC, re.M | re.I)
        #Regex for getting text after contains statement
        self._RX_CONTAINS = "\n\s*contains(?P<remainder>.+)"
        self.RE_CONTAINS = re.compile(self._RX_CONTAINS, re.DOTALL | re.I)
        
    def parse_signature(self, statement, element, module=None):
        """Parses the specified line as a new version of the signature for 'element'.

        :arg statement: the string that has the new signature.
        :arg element: the code element whose signature will be changed.
        :arg module: for real-time module update, the module who the new element
          will be child to.
        """
        #Types don't have their signatures updated in any way that affects
        #the intellisense other than a name change.
        tmatch = self.RE_SIG.match(statement)
        result = (None, None, None)

        if tmatch is not None:
            name = tmatch.group("name")

            #We also need to overwrite the modifiers list with the new signature
            modifiers = tmatch.group("modifiers")
            if modifiers is not None:
                cleanmods = re.split("[\s,]+", modifiers.strip())
            else:
                cleanmods = []
            if module is None:
                element.update_name(name)
                element.modifiers.extend(cleanmods)
            else:
                result = (CustomType(name, modifiers, {}, module), 
                          tmatch.start(), tmatch.end())
                                  
        return result

    def parse_line(self, statement, element, mode):
        """As part of real-time update, parses the statement and adjusts the attributes
        of the specified CustomType instance to reflect the changes.

        :arg statement: the lines of code that was added/removed/changed on the 
          element after it had alread been parsed. The lines together form a single
          continuous code statement.
        :arg element: the CustomType instance to update.
        :arg mode: 'insert', or 'delete'.
        """
        if element.incomplete:
            #We need to check for the end_token so we can close up the incomplete
            #status for the instance.
            if element.end_token in statement:
                element.incomplete = False
                return

        #This method deals with updating the *body* of the type declaration. The only
        #possible entries in the body are member variable declarations and type
        #executable definitions.
        self._process_execs_contents(statement, element.module.name, element, mode)      
        self._rt_parse_members(statement, element, mode)
        
    def _rt_parse_members(self, statement, element, mode):
        """As part of parse_line() check for type members in the statement."""
        if mode == "delete":
            self._rt_members_delete(element, statement)
        elif mode == "insert":
            self._rt_members_add(element, statement)

    def _rt_members_add(self, element, statement):
        """Finds all the member declarations in 'statement' and adds the
        corresponding instances to element.members."""
        members = self.vparser.parse(statement, None)
        for member in members:
            single = members[member]
            single.parent = element
            element.members[member] = single

    def _rt_members_delete(self, element, statement):
        """Finds all the member declarations in 'statement' and removes the
        corresponding instances from element.members."""
        removals = self.vparser.parse(statement, None)
        for member in removals:
            if member in element.members:
                del element.members[member]

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
            result[name.lower()] = self._process_type(name, cleanmods, contents, module, match)
            if "public" in result[name.lower()].modifiers:
                module.publics[name.lower()] = 1

        #Set the types we found in the module and then move the embedded
        #ones into their correct parent executables.
        module.types = result
        module.update_embedded("types")

    def _process_type(self, name, modifiers, contents, module, match):
        """Processes a regex match of a type's contents."""
        #First, we need to see if the types children are private.
        if self.RE_PRIV.search(contents):
            modifiers.append("private contents")

        #Next, we need to parse out all the members of the type and their docstrings
        members = self.vparser.parse(contents, None)
        
        #Now we can create the type code element and handle the member docstrings
        t = CustomType(name, modifiers, members, module)
        #parse out all the executables including the finalizer
        execs = self._process_execs(contents, module.name, t)

        #Set the regex start and end char indices
        t.start, t.end = module.absolute_charindex(match.string, match.start(),
                                                             match.end())

        #Update the parent for embedded members and executables
        for key in list(t.members.keys()):
            t.members[key].parent = t
        for key in list(t.executables.keys()):
            t.executables[key].parent = t

        #Extract the docstrings from the type body and associate them with their members
        memdocs = self.docparser.parse_docs(contents, t)
        if name in memdocs:
            docs = self.docparser.to_doc(memdocs[name][0], name)
            self.docparser.process_memberdocs(docs, t)

        return t

    def update_docs(self, t, module):
        """Updates the documentation for the specified type using the module predocs."""        
        #We need to look in the parent module docstrings for this types decorating tags.
        key = "{}.{}".format(module.name, t.name)
        if key in module.predocs:
            t.docstring = self.docparser.to_doc(module.predocs[key][0], t.name)
            t.docstart, t.docend = (module.predocs[key][1], module.predocs[key][2])

    def _process_execs(self, contents, modulename, atype, mode="insert"):
        """Extracts all the executable methods that belong to the type."""
        #We only want to look at text after the contains statement
        match = self.RE_CONTAINS.search(contents)

        #It is possible for the type to not have any executables
        if match is not None:
            exectext = match.group("remainder")
            self._process_execs_contents(exectext, modulename, atype, mode)

    def _process_execs_contents(self, exectext, modulename, atype, mode="insert"):
        if exectext:
            for mexec in self.RE_EXEC.finditer(exectext):
                name = mexec.group("name")
                if mode == "insert":
                    modifiers = mexec.group("modifiers")
                    pointsto = mexec.group("points")
                    if pointsto is not None:
                        pointsto = "{}.{}".format(modulename, pointsto)
                    e = TypeExecutable(name, modifiers, None, pointsto)
                    atype.executables[name] = e

                elif mode == "delete":
                    if name in atype.executables:
                        del atype.executables[name]

        #Executables at the type level don't have any documentation because each
        #one points to a module-level executable that will be documented.
