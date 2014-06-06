from difflib import Differ

class CodeBlock(object):
    """Represents a contiguous block of code in the buffer text.

    :arg start: the starting line number of this block in the file.
    :arg end: the ending line number.
    :arg contents: a list of strings representing the continuous block
      of lines of code.
    :arg iscached: specifies whether the lines belong to the cached
      version of the code or the buffer version.
    :arg action: a +, - or ? signifying how the lines should be treated.
    """
    def __init__(self, start, end, contents, action):
        self.start = start
        self.end = end
        self.contents = contents
        self.action = ""
        self.decoratables = {}

    @property
    def isline(self):
        """Specifies whether this block is actually just a single line."""
        return self.end - self.start == 0

class LiveDiffer(object):
    """Parses out the differences between the cached version of a
    code file and the one sent in from emacs to see what has been
    added or removed since the last parsing.

    :arg codeparser: an instance of the code parser that has cached
      versions of the other modules already parsed.
    :arg source: the string representing the current state of the
      buffer in emacs that should be diffed with the cached module.
    """
    def __init__(self, codeparser):
        self.parser = codeparser
        self.d = Differ()

    def diff(self, source, cached):
        """Returns a list of CodeBlock objects that describe which changes to
        parse in the new source code in order to update the cached module
        representation.

        :arg source: the string representing the current state of the
          buffer in emacs that should be diffed with the cached module.
        :arg cached: the string representation of the cached module.
        """
        lines = self._get_diff(source, cached)
        return self._split_diff(lines)        
        
    def _split_diff(self, lines):
        """Breaks the diff results up by who was affected, associates
        contiguous block of code together.

        :arg lines: the lines returned from the Differ.
        """
        result = []
        #Keep track of which line we are on in each file.
        linec = 0 #Cached line count
        liner = 0 #Source line count

        #These variables keep track of continuous blocks of code that are
        #either additions or removals
        contig = []
        #When true, the current block of contig lines is in the emacs (source)
        #file contents.
        src_contig = False

        for i in range(len(lines)):
            l = lines[i]
            code = l[:2]

            #Safely get the codes of the next and previous lines so we can handle
            #the way that ? is presented.
            if i + 1 < len(lines):
                coden = lines[i + 1][:2] #The code for the next line
            else:
                coden = ""
            if i - 1 > 0:
                codep = lines[i - 1][:2] #The code for the previous line
            else:
                codep = ""

            if coden == "? ":
                #The code is in both files; the previous entry was the line in
                #cached file, the next entry is the line it should be replaced
                #with. We will need to re-parse the line using the new code to
                #update possible var references etc.
                linec += 1
            elif codep == "? ":
                liner += 1
                result.append(CodeBlock(linec, linec, l, "?"))
            elif code == "+ ":
                liner += 1
                self._get_block_result(linec, contig, "-", src_contig, result)
                src_contig = True
            elif code == "- ":
                linec += 1
                self._get_block_result(liner, contig, "+", not src_contig, result)
                src_contig = False
            elif code == "  ":
                linec += 1
                liner += 1

        return result

    def _get_block_result(self, linenum, contig, mode, src_contig, result):
        """Appends a code block for the specified mode and contig list if required."""
        #See what the last contig block was made of (cached/source files)
        if src_contig:
            contig.append(l)
        else:
            #We had a continuous list of lines that were opposite type from
            #the opposite file. We need to create block object to store them.
            b = CodeBlock(linec - len(contig), linec, contig, "-")
            result.append(b)
            #Reset the contig block tracker to be in the opposite mode.
            contig = [ l ]

    def _get_diff(self, source, cached):
        """Gets a list of lines from the Differ that represents what
        has changed between the source and the cached module.

        :arg source: the string representing the current state of the
          buffer in emacs that should be diffed with the cached module.
        :arg path: the full path to the file being edited in the buffer.
        """
        if cached != "":
            return list(self.d.compare(cached, source))
        else:
            return []
       
class LineParser(object):
    """Parses individual lines or blocks of code to retrieve code element
    representations from them. Exposes low level regex methods for type
    testing of strings."""
    def __init__(self, codeparser):
        self.parser = codeparser
        self.modulep = codeparser.modulep
        self.docparser = self.modulep.docparser
        self.tparser = self.modulep.tparser
        self.xparser = self.modulep.xparser
        self.modlines = {}

    def _get_element(self, linenum, module):
        """Gets the code element who owns the line number specified.

        :arg linenum: the number of the line to test ownership for.
        :arg module: the module whose children will be tested.
        """
        if not module.name in self.modlines:
            self._load_modlines(module)

        #Just cycle through the line numbers in order from the dictionary
        #and see which is the first one.
        element = None
        end = 0
        for i in sorted(self.modlines[module.name].keys()):
            if i <= linenum:
                end, element = self.modlines[module.name][i]
                if linenum <= end:
                    break
                else:
                    #Something is up, the line number doesn't fall within
                    #the full body of this code element.
                    element = None
            
        if element is not None:
            return element
        else:
            #It must be between code elements, which makes it part of
            #the parent modules contents, return the parent.
            return module

    def _load_modlines(self, module):
        """Loads the line numbers spanned by each code element in the
        module into a dictionary in modlines under the module name."""
        result = {}
        #We will save the *start* line numbers as the keys and then
        #the end line numbers and object as a tuple in the values.
        #We need to look at types and executables.
        for xkey in module.executables:
            x = module.executables[xkey]
            result[x.start] = (x.end, x)

        for tkey in module.types:
            t = module.types[tkey]
            result[t.start] = (t.end, t)

        self.modlines[module.name] = result

    def is_decoratable(self, string):
        """Tests whether the specified string is a decoratable code element. 
        Uses the docstring decorator regex. Can be a module, type, subroutine
        function."""
        return self.docparser.RE_DECOR.match(string)

    def is_terminator(self, string, match):
        """Determines whether the specifed string is a terminator for
        the decorator match."""
        ftype = match.group("functype")
        return "end {}".format(ftype) in string

    def get_terminator(self, match):
        """Gets an appropriate terminator string based on the decorator
        match for a code element."""
        ftype = match.group("functype")
        name = match.group("name")
        return "  end {} {}".format(ftype, name)

class ModuleUpdater(object):
    """Updates the representations of the fortran modules in memory using
    new source code supplied from the emacs buffer."""
    def __init__(self, codeparser):
        self.parser = codeparser
        self.differ = LiveDiffer(codeparser)
        self.linep = LineParser(codeparser)

        self._actions = {
            "+": self._handle_new,
            "-": self._handle_kill,
            "?": self._handle_edit
        }

    def update(self, source, path):
        """Updates all references in the cached representation of the module
        with their latest code from the specified source code string that
        was extracted from the emacs buffer.

        :arg source: the string representing the current state of the
          buffer in emacs that should be diffed with the cached module.
        :arg path: the full path to the file being edited in the buffer.
        """
        #First, use the differ to get a list of the changes that need to
        #be made. Then cycle through the changes and handle each one.
        module = self._get_module(path)
        changes = self.differ.diff(source, module.refstring)
        for change in changes:
            #Before we can process the change, we need to make sure that
            #the code is clean. For example, if a subroutine keyword has
            #been specified to mark the start of the routine, we need to
            #make sure that it has a matching end subroutine, even if the
            #user hasn't typed one in yet.
            self._clean_block(change)
            
            #Now we can process the action based on how it changes the
            #cached representations of the code elements.
            self._actions[change.action](change, module)
            
    def _handle_new(self, change, module):
        """Handles new code that was added to the code file in the buffer
        but wasn't in the cached representation of the module."""
        #The block changes are either a single line to be added or a whole
        #new code element to insert. If there are decoratables in the block
        #we will do those separately.
        elements = []
        if len(change.decoratables.keys()) > 0:
            #The executable parser will ignore anything in the text that doesn't
            #match a subroutine or function. We can just join all the contents
            #of this block and pass it in as if it were the contents of the module
            #after the contains keyword.
            contents = "".join(change.contents)
            parsed = []

            #The regex match on the decoratable tells us which parser needs to
            #be called to do the parsing.
            for dec in change.decoratables:
                mtype = change.decoratables[dec].group("element")
                if mtype is not None and mtype not in parsed:
                    parser = self._get_parser(mtype)
                    parser.parse_block(contents, module)
                    parsed.append(mtype)
        else:
            #These are just ad-hoc lines to be added. We need to find out
            #which code element the line numbers belong to. If there is
            #more than one line in the change's contents, then the change
            #forms a contiguous block, so we only need to check the first
            #line of the set for where it belongs.
            element = self.linep._get_element(change.start, module)
            for line in change.contents:
                parser = self._get_parser(element)
                if parser is not None:
                    parser.parse_line(line, element, "+")

    def _handle_kill(self, change, module):
        """Handles the removal of some code from the cached representation
        of the module."""
        
    def _handle_edit(self, change, module):
        """Handles a line of code in the cached representation that is
        different in the buffer being edited."""
        
    def _get_parser_s(self, element):
        """Returns the appropriate parser based on the string type of the element."""
        if element == "subroutine" or element == "function":
            return self.linep.xparser
        elif element == "type":
            return self.linep.tparser
        elif element == "module":
            return self.linep.modulep
        else:
            return None            

    def _get_parser(self, element):
        """Returns the appropriate parser based on the type of the element."""
        #The most common type of addition will be to an executable's contents
        #then to the parent module and finally to a derived type definition.
        if type(element) == type(""):
            return self._get_parser_s(element)
        elif isinstance(element, Executable):
            return self.linep.xparser
        elif isinstance(element, Module):
            return self.linep.modulep
        elif isinstance(element, CustomType):
            return self.linep.tparser
        else:
            return None

    def _clean_block(self, block):
        """Ensures that the specified block is ready to be parsed by
        inserting missing terminating "end" statements."""
        hasdecor = []
        terminates = []

        for i in range(len(block.contents)):
            l = block.contents[i]
            m = self.linep.is_decoratable(l)
            if m is not None:
                hasdecor.append(m)
                #We keep track of the decoratable status of the line so we don't
                #have to test whether it is decoratable again later.
                block.decoratables[i] = m
            if len(hasdecor) > len(terminates) and \
               self.linep.is_terminator(hasdecor[len(terminates)]):
                terminates.append(True)

        #Now see if we need to append any terminators
        if len(terminates) < len(hasdecor):
            for i in range(len(terminates), len(hasdecor)):
                block.contents.append(self.linep.get_terminator(hasdecor[i]))
    
    def _get_module(self, path):
        """Gets the cached version of the module.

        :arg path: the full path to the file being edited in the buffer.
        """ 
        #We need to determine the module from the file name of the path
        #There could be mappings for the file name to module name.
        name = path.split("/")[-1]
        module = None

        for mapkey in self.parser.mappings:
            if self.parser.mappings[mapkey] == name.lower():
                module = self.parser.mappings[mapkey]
                break

        if module is None:
            module = name.replace(".f90", "")
            
        #The parser would already have the module loaded before calling
        #a differ. We can assume that it is there.
        if module in self.parser.modules:
            return self.parser.modules[module]
        else:
            print "ERROR: there is no cached module to diff with."
            return None

