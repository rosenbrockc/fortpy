import re
from ..elements import Module
from .docstring import DocStringParser
from .variable import VariableParser
from .types import TypeParser
from .executable import ExecutableParser

class ModuleParser(object):
    """Extracts modules from fortran code files."""
    def __init__(self):
        self.setup_regex()
        self.vparser = VariableParser()
        self.docparser = DocStringParser()
        self.tparser = TypeParser(self.vparser, self.docparser)
        self.xparser = ExecutableParser(self.vparser, self.docparser)
    
    def setup_regex(self):
        """Sets up compiled regex objects for parsing code elements."""
        #Regex for extracting modules from the code
        self._RX_MODULE = r"(\n|^)\s*module\s+(?P<name>[a-z0-9_]+)(?P<contents>.+?)end\s*module"
        self.RE_MODULE = re.compile(self._RX_MODULE, re.I | re.DOTALL)
        #Regex for use statements in a module
        self._RX_USE = r"^\s*use\s+(?P<name>[^,]+?)(,\s+only:(?P<only>[A-Za-z0-9_\s,]+?))?$"
        self.RE_USE = re.compile(self._RX_USE, re.I | re.M)
        #Regex for finding if the module is private
        self._RX_PRIV = "private.+?(type|contains)"
        self.RE_PRIV = re.compile(self._RX_PRIV, re.DOTALL | re.I)
        #Regex for finding publcily labeled members declared using public keyword.
        self._RX_PUBLIC = r"\n\s*public\s+(?P<methods>[A-Za-z0-9_,\s&\n]+)"
        self.RE_PUBLIC = re.compile(self._RX_PUBLIC, re.I)
        #Regex for finding text before type or contains declarations that may contian members.
        self._RX_MEMBERS = "(?P<preamble>.+?)(type|contains)"
        self.RE_MEMBERS = re.compile(self._RX_MEMBERS, re.DOTALL | re.I)

    def parse(self, string, parent):
        """Extracts any modules from the specified fortran code file."""
        #First, get hold of the docstrings  for all the modules so that we can
        #attach them as we parse them.
        moddocs = self.docparser.parse_docs(string)
        #Now look for modules in the file and then match them to their decorators.
        matches = self.RE_MODULE.finditer(string)
        result = []
        for rmodule in matches:
            name = rmodule.group("name")
            contents = rmodule.group("contents")
            module = self._process_module(name, contents, parent, rmodule)
            #Check whether the docparser found docstrings for the module.
            if name in moddocs:                
                module.docstring = self.docparser.to_doc(moddocs[name], name)
            result.append(module)
        return result

    def _process_publics(self, contents):
        """Extracts a list of public members, types and executables that were declared using
        the public keyword instead of a decoration."""
        matches = self.RE_PUBLIC.finditer(contents)
        result = {}
        for public in matches:
            methods = public.group("methods")
            for item in re.split(r"[\s&\n,]+", methods.strip()):
                self._dict_increment(result, item)
        return result

    def _process_module(self, name, contents, parent, match):
        """Processes a regex match for a module to create a CodeElement."""
        #First, get hold of the name and contents of the module so that we can process the other
        #parts of the module.
        modifiers = []

        #We need to check for the private keyword before any type or contains declarations
        if self.RE_PRIV.search(contents):
            modifiers.append("private")
        #The only other modifier for modules ought to be implicit none
        if re.search("implicit\s+none", contents):
            modifiers.append("implicit none")

        #Next, parse out the dependencies of the module on other modules
        dependencies = self._parse_use(contents)
        publics = self._process_publics(contents)

        #We can now create the CodeElement
        result = Module(name, modifiers, dependencies, publics, contents, parent)
        result.start = match.start()
        result.end = match.end()
        result.refstring = match.string

        #It is possible for the module to have members, parse those
        self._parse_members(contents, result)
        #Now we just need to handle the custom types and executables
        self.tparser.parse(result)
        self.xparser.parse(result)
        
        return result            

    def _parse_use(self, string):
        """Extracts use dependencies from the innertext of a module."""
        result = {}
        for ruse in self.RE_USE.finditer(string):
            name = ruse.group("name").strip()
            if ruse.group("only"):
                only = ruse.group("only").split(",")
                for method in only:
                    key = "{}.{}".format(name, method.strip())
                    self._dict_increment(result, key)
            else:
                self._dict_increment(result, name)
        return result

    def _dict_increment(self, dictionary, key):
        """Increments the value of the dictionary at the specified key."""
        if key in dictionary:
            dictionary[key] += 1
        else:
            dictionary[key] = 1

    def _parse_members(self, contents, module):
        """Extracts any module-level members from the code. They must appear before
        any type declalations."""
        match = self.RE_MEMBERS.match(contents)
        preamble = match.group("preamble")
        #Get a dictionary of all the members in this module body
        members = self.vparser.parse(preamble, module)

        #The docstrings for these members will appear as member tags in the same
        #preamble text. Extract them.
        module.predocs = self.docparser.parse_docs(preamble, module.name)
        if module.name in module.predocs:
            #We can only do member docstrings if the module had internal docstrings
            #that may refer to members.
            memdocs = self.docparser.to_doc(module.predocs[module.name], module.name)
            remainingdocs = self.docparser.process_memberdocs(memdocs, module)
            module.predocs[module.name] = remainingdocs
