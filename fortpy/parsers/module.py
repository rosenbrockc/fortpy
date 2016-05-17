import re
from ..elements import Module
from .docstring import DocStringParser
from .variable import VariableParser
from .types import TypeParser
from .executable import ExecutableParser
from .interface import InterfaceParser

class ModuleParser(object):
    """Extracts modules from fortran code files."""
    def __init__(self):
        self.setup_regex()
        self.vparser = VariableParser()
        self.docparser = DocStringParser()
        self.iparser = InterfaceParser(self.docparser)
        self.tparser = TypeParser(self.vparser, self.docparser)
        self.xparser = ExecutableParser(self.vparser, self.docparser)
    
    def setup_regex(self):
        """Sets up compiled regex objects for parsing code elements."""
        #Regex for extracting modules from the code
        self._RX_MODULE = r"(\n|^)\s*module\s+(?P<name>[a-z0-9_]+)(?P<contents>.+?)end\s*module"
        self.RE_MODULE = re.compile(self._RX_MODULE, re.I | re.DOTALL)
        self._RX_PROGRAM = r"(\n|^)\s*program\s+(?P<name>[a-z0-9_]+)(?P<contents>.+?)end\s*program"
        self.RE_PROGRAM = re.compile(self._RX_PROGRAM, re.I | re.DOTALL)

        #Regex for use statements in a module
        self._RX_USE = r"^\s*use\s+(?P<name>[^,]+?)(\s*,\s+only\s*:(?P<only>[A-Za-z0-9_\s,]+?))?$"
        self.RE_USE = re.compile(self._RX_USE, re.I | re.M)
        #Regex for finding if the module is private
        self._RX_PRIV = "private.+?(type|contains)"
        self.RE_PRIV = re.compile(self._RX_PRIV, re.DOTALL | re.I)
        #Regex for finding publcily labeled members declared using public keyword.
        self._RX_PUBLIC = r"\n\s*public\s+(?P<methods>[A-Za-z0-9_,\s&\n]+)"
        self.RE_PUBLIC = re.compile(self._RX_PUBLIC, re.I)
        #Regex for finding text before type or contains declarations that may contian members.
        self._RX_MEMBERS = "(?P<preamble>.+?)(\s+type[,\s]|contains)"
        self.RE_MEMBERS = re.compile(self._RX_MEMBERS, re.DOTALL | re.I)

        self._RX_PRECOMP = r"#endif"
        self.RE_PRECOMP = re.compile(self._RX_PRECOMP, re.I)

    def rt_update(self, statement, element, mode, linenum, lineparser):
        """Performs a real-time update of the specified statement that is in the body of the
        module.
        
        :arg statement: the lines of code that was added/removed/changed on the 
          element after it had alread been parsed. The lines together form a single
          continuous code statement.
        :arg element: the Module instance to update.
        :arg mode: 'insert', 'replace', or 'delete'.
        """
        #First find out if we are passed the CONTAINS section separating executables from
        #the type and member definitions. In order for this code to be reached, the lines
        #that are being changed are *between* definitions that already exist. The likelihood
        #is that they are *new* definitions of members, types or executables.
        if linenum <= element.contains_index:
            #we only have to look for type and member definitions.
            self._rt_parse_members(statement, element, mode)
            self._rt_parse_types(statement, element, mode, lineparser)
        else:
            #we only have to deal with executables.
            self._rt_parse_execs(statement, element, mode, lineparser)

    def _rt_parse_execs(self, statement, element, mode, lineparser):
        """As part of parse_line(), checks for new executable declarations in the statement."""
        #This is the same deal as with _rt_parse_types() below.
        if mode == "insert":
            enew, start, end = self.xparser.parse_signature(statement, element, element)
            if enew is not None:
                enew.start, enew.end = lineparser.absolute_charindex(statement, start, end)
                enew.incomplete = True
                element.executables[enew.name.lower()] = enew
                lineparser.additions.append((enew, element))

    def _rt_parse_types(self, statement, element, mode, lineparser):
        """As part of parse_line(), checks for new type declarations in the statement."""
        if mode == "insert":
            #Since we got to this point, there is *no* code element that owns the current
            #line which is being replaced; we are merely checking to see if the new line
            #being entered matches a type definition, do the same thing as "insert"
            tnew, start, end = self.tparser.parse_signature(statement, element, element)
            #We need to set sensible boundaries for 'start' and 'end' so that if the lines
            #immediately after this one are member definitions (for example) they get
            #associated correctly with this type.
            if tnew is not None:
                tnew.start, tnew.end = lineparser.absolute_charindex(statement, start, end)
                tnew.incomplete = True
                element.types[tnew.name.lower()] = tnew
                lineparser.additions.append((tnew, element))
#        elif mode == "delete":
            #This line is not part of any existing code element. Pointless to do anything
            #with delete since it probably is just a member that got removed.

    def _rt_parse_members(self, statement, element, mode):
        """As part of parse_line(), checks for member declarations in the statement."""
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

    def parse(self, string, parent, module=True, filepath=None):
        """Extracts modules *and* programs from a fortran code file.

        :arg string: the contents of the fortran code file.
        :arg parent: the instance of CodeParser that will own the return Module.
        :arg module: when true, the code file will be searched for modules; otherwise
          it will be searched for programs.
        """
        if module:
            result = self._parse_modules(string, parent, filepath)
        else:
            result = self._parse_programs(string, parent, filepath)
        return result
    
    def _parse_programs(self, string, parent, filepath=None):
        """Extracts a PROGRAM from the specified fortran code file."""
        #First, get hold of the docstrings  for all the modules so that we can
        #attach them as we parse them.
        moddocs = self.docparser.parse_docs(string)
        #Now look for modules in the file and then match them to their decorators.
        matches = self.RE_PROGRAM.finditer(string)
        result = []
        for rmodule in matches:
            name = rmodule.group("name").lower()
            contents = re.sub("&[ ]*\n", "", rmodule.group("contents"))
            module = self._process_module(name, contents, parent, rmodule, filepath)
            #Check whether the docparser found docstrings for the module.
            if name in moddocs:                
                module.docstring = self.docparser.to_doc(moddocs[name][0], name)
                module.docstart, module.docend = module.absolute_charindex(string, moddocs[name][1],
                                                                           moddocs[name][2])
            result.append(module)
        return result

    def _parse_modules(self, string, parent, filepath=None):
        """Extracts any modules from the specified fortran code file."""
        #First, get hold of the docstrings  for all the modules so that we can
        #attach them as we parse them.
        moddocs = self.docparser.parse_docs(string)
        #Now look for modules in the file and then match them to their decorators.
        matches = self.RE_MODULE.finditer(string)
        result = []
        for rmodule in matches:
            name = rmodule.group("name").lower()
            contents = re.sub("&[ ]*\n", "", rmodule.group("contents"))
            module = self._process_module(name, contents, parent, rmodule, filepath)
            #Check whether the docparser found docstrings for the module.
            if name in moddocs:                
                module.docstring = self.docparser.to_doc(moddocs[name][0], name)
                module.docstart, module.docend = module.absolute_charindex(string, moddocs[name][1],
                                                                           moddocs[name][2])

            #Before we append the module to the list, we need to update its list
            #of publics if it hasn't explicitly been declared as private.
            module.all_to_public()                
            result.append(module)
        return result

    def _process_publics(self, contents):
        """Extracts a list of public members, types and executables that were declared using
        the public keyword instead of a decoration."""
        matches = self.RE_PUBLIC.finditer(contents)
        result = {}
        start = 0
        for public in matches:
            methods = public.group("methods")
            #We need to keep track of where the public declarations start so that the unit
            #testing framework can insert public statements for those procedures being tested
            #who are not marked as public
            if start == 0:
                start = public.start("methods")
            for item in re.split(r"[\s&\n,]+", methods.strip()):
                if item.lower() in ["interface", "type", "use"]:
                    #We have obviously reached the end of the actual public
                    #declaration in this regex match.
                    break
                self._dict_increment(result, item.lower())
        return (result, start)

    def _process_module(self, name, contents, parent, match, filepath=None):
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
        publics, pubstart = self._process_publics(match.string)

        #We can now create the CodeElement
        result = Module(name, modifiers, dependencies, publics, contents, parent)
        if filepath is not None:
            result.filepath = filepath.lower()
        result.start = match.start()
        result.end = match.end()
        result.refstring = match.string
        result.set_public_start(pubstart)
        if self.RE_PRECOMP.search(contents):
            result.precompile = True

        self.xparser.parse(result)
        self.tparser.parse(result)

        #It is possible for the module to have members, parse those
        self._parse_members(contents, result)
        self.iparser.parse(result)

        #Now we can update the docstrings for the types. They rely on data
        #extracted during parse_members() which is why they have to run
        #separately over here.
        for t in result.types:
            self.tparser.update_docs(result.types[t], result)

        return result            

    def _parse_use(self, string):
        """Extracts use dependencies from the innertext of a module."""
        result = {}
        for ruse in self.RE_USE.finditer(string):
            #We also handle comments for individual use cases, the "only" section
            #won't pick up any comments.
            name = ruse.group("name").split("!")[0].strip()
            if name.lower() == "mpi":
                continue
            
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
        #We need to get hold of the text before the module's main CONTAINS keyword
        #so that we don't find variables from executables and claim them as
        #belonging to the module.
        icontains = module.contains_index
        ichar = module.charindex(icontains, 0)
        module.preamble = module.refstring[:ichar]

        #Get a dictionary of all the members in this module body
        #We only want to look at variable definitions before the first type
        lowest = ichar
        remove = [] #Will use later below, see next comment
        for t in module.types:
            remove.append((module.types[t].start, module.types[t].end))
            if module.types[t].start < lowest:
                lowest = module.types[t].start

        module.members.update(self.vparser.parse(contents[:lowest-(module.start + 10 + len(module.name))], module))

        #The docstrings for these members will appear as member tags in the same
        #preamble text. We can't use the entire preamble for this because member
        #docs inside of a type declaration will show up as belonging to the
        #module, when in fact, they don't.
        remove.sort(key=lambda tup: tup[0])
        retain = []
        cur_end = 0
        for rem in remove:
            signature = module.refstring[rem[0]+1:rem[1]].index("\n") + 2
            keep = module.refstring[cur_end:rem[0] + signature]
            cur_end = rem[1]
            retain.append(keep)

        #If there weren't any types in the module, we still want to get at the docs in
        #the preamble.
        if len(remove) == 0:
            retain = module.preamble

        docsearch = "".join(retain)
        module.predocs = self.docparser.parse_docs(docsearch, module)
        if module.name in module.predocs:
            #We can only do member docstrings if the module had internal docstrings
            #that may  to members.
            memdocs = self.docparser.to_doc(module.predocs[module.name][0], module.name)
            remainingdocs = self.docparser.process_memberdocs(memdocs, module)
            module.predocs[module.name] = remainingdocs
