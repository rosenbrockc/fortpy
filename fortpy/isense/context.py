from .. import msg
from . import cache
from fortpy.elements import Function, Subroutine, ValueElement, Executable, Module
import re

#This module looks at the position of the cursor in the source code
#and determines what kind of symbol string needs to be passed back
#to the Script class for the intellisense.
class UserContext(object):
    """Represents the context of the user in the Fortran code file
    so that intellisense can be taylored to the logical position 
    in the file.

    :attr modulename: the name of the module the user is editing.
    :attr section: either 'vars', 'types', 'execution', or 'contains'
      depending on whether the cursor is still in the section where
      variables are declared, special user-defined types are declared
      or (if in a PROGRAM) the execution section. 'contains' refers to
      subroutine/function declarations after the CONTAINS keyword.
    :attr el_type: the type of the element the user is in at the moment
      can be any of the types in fortpy.elements.
    :attr el_name: the name of the element that their cursor is in. Is
      set to be the module if UserContext.section is not 'types' or 'contains'
    :attr el_section: the section within the element where the cursor
      is placed. Can be 'params', 'vars' or 'body' if the cursor is in
      the argument list for a subroutine/function, in the variable
      declaration section or in the body of the function.
    :attr el_call: set to either 'sub', 'fun', 'arith', 'assign' or 'name'
      if the cursor is on the symbol or argument list for a subroutine or 
      function, or if the cursor is in some arithmetic.
    """
    def __init__(self, source, pos, fullpath, parser_key="default"):
        """Creates an instance of the UserContext."""
        self._strsource = source
        self._source = source.splitlines()
        self._orig_path = fullpath
        self.pos = pos
        self.parser = cache.parser(parser_key)

        #Initialize the lazy variables
        self._symbol = None
        self._full_symbol = None
        self._short_symbol = None
        self._short_full_symbol = None
        self.element = None
        self._call_index = None
        self._cachedstr = None
        self._exact_match = None

        #Setup all the regular expressions by grabbing cached versions so
        #we don't have to recompile them the whole time.
        self._setup_regex()
        self._contextualize()

    def __str__(self):
        output = []
        output.append("\n\nSYMBOLS: {}; {}".format(self.symbol, self.full_symbol))
        output.append("SHORT: {}; {}".format(self.short_symbol, self.short_full_symbol))
        output.append("MODULE: " + self.modulename)
        output.append("SECTION: " + self.section)
        tel = self.el_type or self.element.name
        output.append("TYPE: {}".format(str(tel)))
        output.append("NAME: " + self.el_name)
        output.append("EL SECTION: " + self.el_section)
        output.append("CALL :" + self.el_call)
        return "\n".join(output)

    @property
    def refstring(self):
        """Returns the *unsplit* source code string for the current code context."""
        return self._strsource

    @property
    def bufferstr(self):
        """Returns the full string of the file's current state in the buffer."""
        return self._source

    @property
    def cachedstr(self):
        """Returns the full string of the file contents from the cache for
        the file that we are currently providing intellisense for."""
        if self._cachedstr is None:
            if self.module is not None:
                refstring = self.module.refstring
                self._cachedstr = refstring.splitlines()
            else:
                self._cachedstr = []

        return self._cachedstr

    @property
    def current_line(self):
        """Gets the text of the line under the cursor."""
        return self._source[self.pos[0]]

    @property
    def exact_match(self):
        """Returns the symbol under the cursor looking both directions as part
        of a definition lookup for an exact match.
        """
        #We don't have to worry about grouping or anything else fancy. Just
        #loop through forward and back until we hit a character that can't be
        #part of a variable or function name.
        if self._exact_match is None:
            i = self.pos[1] - 1
            start = None
            end = None
            line = self.current_line
            terminators = ['(', ')', '\n', ' ', '=', '%', ',']
            while i >= 0 and start is None:
                if line[i] in terminators:
                    start = i + 1
                i -= 1

            i = self.pos[1]
            while i < len(line) and end is None:
                if line[i] in terminators:
                    end = i
                i += 1

            self._exact_match = line[start:end].lower()

        return self._exact_match

    @property
    def short_symbol(self):
        """Gets the symbol for the current cursor *excluding* the
        last character under the cursor."""
        if self._short_symbol is None:
            self._short_symbol = self._symbol_extract(cache.RE_CURSOR, False)

        return self._short_symbol

    @property
    def short_full_symbol(self):
        """Gets the full symbol excluding the character under the cursor."""
        if self._short_full_symbol is None:
            self._short_full_symbol = self._symbol_extract(cache.RE_FULL_CURSOR,
                                                           False, True)

        return self._short_full_symbol

    @property
    def symbol(self):
        """Gets the symbol under the current cursor."""
        if self._symbol is None:
            self._symbol = self._symbol_extract(cache.RE_CURSOR)
            
        return self._symbol

    @property
    def full_symbol(self):
        """Returns the symbol under the cursor AND additional contextual
        symbols in the case of %-separated lists of type members."""
        if self._full_symbol is None:
            self._full_symbol = self._symbol_extract(cache.RE_FULL_CURSOR, brackets=True)

        return self._full_symbol
    
    def _symbol_extract(self, regex, plus = True, brackets=False):
        """Extracts a symbol or full symbol from the current line, 
        optionally including the character under the cursor.

        :arg regex: the compiled regular expression to use for extraction.
        :arg plus: when true, the character under the cursor *is* included.
        :arg brackets: when true, matching pairs of brackets are first removed
          before the regex is run.
        """
        charplus = self.pos[1] + (1 if plus else -1)
        consider = self.current_line[:charplus][::-1]

        #We want to remove matching pairs of brackets so that derived types
        #that have arrays still get intellisense.
        if brackets==True:
            #The string has already been reversed, just run through it.
            rightb = []
            lastchar = None
            for i in range(len(consider)):
                if consider[i] == ")":
                    rightb.append(i)
                elif consider[i] == "(" and len(rightb) > 0:
                    lastchar = i
                    rightb.pop()

            if lastchar is not None:
                consider = '%' + consider[lastchar+1:]

        rematch = regex.match(consider)
        if rematch is not None:
            return rematch.group("symbol")[::-1]
        else:
            return ""

    @property
    def call_arg_index(self):
        """Determines the index of the parameter in a call list using
        string manipulation and context information."""
        #The function name we are calling should be in el_name by now
        if self._call_index is None:
            if (self.el_section == "body" and
                self.el_call in [ "sub", "fun" ]):
                #Get hold of the element instance of the function being called so
                #we can examine its parameters.
                fncall = self.el_name
                if fncall in self.current_line[:self.pos[1]]:
                    args = self.current_line[:self.pos[1]].split(fncall)[1]
                else:
                    args = ""

                #This handles the case of the bracket-complete.
                if args == "":
                    return 0
                #The nester requires each start bracket to have an end bracket
                if args[-1] != ")":
                    args += ")"

                #Pyparsing handles calls where functions are being called as
                #the values for parameters like function(a, fun(c,d), r).
                try:
                    nl = cache.nester.parseString(args).asList()[0]
                    clean = [n for n in nl if not isinstance(n, list) and n != ","]
                    #We need to specially handle the case where they have started typing the
                    #name of the first parameter but haven't put a ',' in yet. In that case, 
                    #we are still on the first argument.
                    if clean[-1][-1] != ',':
                        self._call_index = len(clean) - 1
                    else:
                        self._call_index = len(clean)
                except:
                    msg.warn("ARG INDEX: lookup failed on bracket parsing.")
            
        return self._call_index

    def _setup_regex(self):
        """Sets up the constant regex strings etc. that can be used to
        parse the strings for determining context."""
        self.RE_COMMENTS = cache.RE_COMMENTS
        self.RE_MODULE = cache.RE_MODULE
        self.RE_TYPE = cache.RE_TYPE
        self.RE_EXEC = cache.RE_EXEC
        self.RE_MEMBERS = cache.RE_MEMBERS
        self.RE_DEPEND = cache.RE_DEPEND

    def _contextualize(self):
        """Finds values for all the important attributes that determine the 
        user's context."""
        line, column = self.pos
        #Get the module top-level information
        self._get_module(line)
        if self.module is None:
            return
        
        #Use the position of the cursor in the file to decide which
        #element we are working on.
        self.element = self.module.get_element(line, column)

        #Now all that's left is to contextualize the line that the
        #cursor is on.
        self._deep_match(line, column)
            
    def _get_module(self, line):
        """Finds the name of the module and retrieves it from the parser cache."""
        #Finding the module name is trivial; start at the beginning of
        #the module and iterate lines until we find the module.
        for sline in self._source:
            if len(sline) > 0 and sline.strip()[0] != "!":
                rmatch = self.RE_MODULE.match(sline)
                if rmatch is not None:
                    self.modulename = rmatch.group("name")
                    break
        else:
            #We don't even have the start of a module in this code file
            return
            
        #Before we carry on with the rest of the context, find the separating
        #CONTAINS keyword so we know whether to look for types or subs/funcs.
        #If the code parser hasn't ever parsed this module, parse it now.
        self.parser.isense_parse(self._orig_path, self.modulename)
        self.module = self.parser.modules[self.modulename]

        if line > self.module.contains_index:
            self.section = "contains"
        #else:
            #We are in the pre-contains section of types and module members
            #We have to do some more complex searching.

    def _deep_match(self, line, column):
        """Checks the contents of executables, types and modules for member
        definitions and updates the context."""
        #Now we just try each of the possibilities for the current line
        if self._match_member(line, column):
            self.el_section = "vars"
            self.el_type = ValueElement
            self.el_name = self.col_match.group("names")
            return
             
        if isinstance(self.element, Executable) or isinstance(self.element, Module):
            #We are inside of a subroutine or function definition
            #It is either params, vars or body. We already tested for variable
            #declarations in the match_member test. Check now to see if we are
            #on the line that defines the function
            self._match_exec(line)
            if self.col_match:
                self.el_section = "params"
                self.el_call = "assign"
                return

            #This regex incorrectly grabbing things like 'if' is functions because
            #they do actually like like functions... We need to filter the list 
            #with special keywords before we claim victory. TODO
            self.col_match = self.RE_DEPEND.match(self._source[line])
            if self.col_match:
                self.el_section = "body"
                self.el_name = self.col_match.group("exec")
                if self.col_match.group("sub") and "call" in self.col_match.group("sub"):
                    self.el_call = "sub"
                    self.el_type = Subroutine
                else:
                    self.el_call = "fun"
                    self.el_type = Function

                return
        #If we are inside a type, we either got the member variable declaration
        #already, or it is a pointer to a method inside of the module. Add the
        #test for that later. TODO
#        if isinstance(self.element, CustomType):
            
        self.el_section = "body"
        self.el_type = None
        #We just need to figure out what kind of a call is being made
        #at this column position, the only thing left is
        if " = " in self._source[line]:
            eqi = self._source[line].index("=")
            if column < eqi:
                self.el_call = "name"
            else:
                self.el_call = "assign"
            self.el_name = self._source[line].split("=")[0].strip()
        elif re.match("\s*call\s+", self._source[line]):
            self.el_call = "sub"
            self.el_name = self._source[line].split("call ")[1]
        else:
            #The only thing left (and the hardest to nail down) is
            #the arithmetic catch all type.
            self.el_call = "arith"
            self.el_name = "[arith]"

    def _match_exec(self, i):
        """Looks at line 'i' for a subroutine or function definition."""
        self.col_match = self.RE_EXEC.match(self._source[i])
        if self.col_match is not None:
            if self.col_match.group("codetype") == "function":
                self.el_type = Function
            else:
                self.el_type = Subroutine
            self.el_name = self.col_match.group("name")
        
            return True
        else:
            return False

    def _match_member(self, i, column):
        """Looks at line 'i' to see if the line matches a module member def."""
        self.col_match = self.RE_MEMBERS.match(self._source[i])
        if self.col_match is not None:
            if column < self._source[i].index(":"):
                self.el_call = "name"
            else:
                self.el_call = "assign"
            
            return True
        else:
            return False

    def _match_type(self, i):
        """Looks at line 'i' to see if the line matches a module user type def."""
        self.col_match = self.RE_TYPE.match(self._source[i])
        if self.col_match is not None:
            self.section = "types"
            self.el_type = CustomType
            self.el_name = self.col_match.group("name")
           
            return True
        else:
            return False
