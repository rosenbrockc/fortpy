from .. import msg
import re
from ..elements import Subroutine, Function, Dependency, Executable, Module, TypeExecutable
import pyparsing

class ExecutableParser(object):
    """Parses subroutine and function definitions from a module body."""
    
    def __init__(self, vparser, docparser):
        self.setup_regex()
        self.vparser = vparser
        self.docparser = docparser
        self.nester = pyparsing.nestedExpr('(', ')')

        #A list of all the intrinsic functions in fortran 2003
        self._intrinsic = [ k.strip().lower() for k in self._intrinsic_functions()]

    def setup_regex(self):
        """Sets up compiled regex objects for parsing the executables from a module."""
        self._RX_CONTAINS = r"^\s*contains[^\n]*?$"
        self.RE_CONTAINS = re.compile(self._RX_CONTAINS, re.M | re.I)

        #Setup a regex that can extract information about both functions and subroutines
        self._RX_EXEC = r"\n[ \t]*((?P<type>character|real|type|logical|integer)?" + \
                        r"(?P<kind>\([a-z0-9_]+\))?)?((?P<modifiers>[\w, \t]+?))?[ \t]*" + \
                        r"(?P<codetype>subroutine|function)\s+(?P<name>[^(]+)" + \
                        r"\s*\((?P<parameters>[^)]*)\)(?P<result>\sresult\([a-z0-9_]+\))?" + \
                        r"(?P<contents>.+?)end\s*(?P=codetype)\s+(?P=name)"
        self.RE_EXEC = re.compile(self._RX_EXEC, re.DOTALL | re.I)
        #Regex for the signature is almost identical to the full executable, but it doesn't
        #look for any contents after the parameter list.
        self._RX_SIG =  r"((?P<type>character|real|type|logical|integer)?" + \
                        r"(?P<kind>\([a-z0-9_]+\))?)?(,?(?P<modifiers>[^\n]+?))?\s*" + \
                        r"(?P<codetype>subroutine|function)\s+(?P<name>[^(]+)" + \
                        r"\s*\((?P<parameters>[^)]*)\)"
        self.RE_SIG = re.compile(self._RX_SIG, re.I)

        #The contents of the executable already have any & line continuations
        #removed, so we can use a simple multiline regex.
        self._RX_ASSIGN = r"^(?P<assignee>[^!<=\n/]+?)=[^\n]+?$"
        self.RE_ASSIGN = re.compile(self._RX_ASSIGN, re.M)

        self._RX_DEPEND = r"^\s*(?P<sub>call\s+)?(?P<exec>[a-z0-9_%]+\s*\([^\n]+)$"
        self.RE_DEPEND = re.compile(self._RX_DEPEND, re.M | re. I)

        self._RX_DEPCLEAN = r"(?P<key>[a-z0-9_%]+)\("
        self.RE_DEPCLEAN = re.compile(self._RX_DEPCLEAN, re.I)

        self._RX_CONST = '[^"\']+(?P<const>["\'][^\'"]+["\'])'
        self.RE_CONST = re.compile(self._RX_CONST)

        self._RX_COMMENTS = r'\s*![^\n"]+?\n'
        self.RE_COMMENTS = re.compile(self._RX_COMMENTS)

    def parse_signature(self, statement, element, module=None):
        """Parses the specified line as a new version of the signature for 'element'.

        :arg statement: the string that has the new signature.
        :arg element: the code element whose signature will be changed.
        """
        #If the signature changes, the user might not have had a chance to add the
        #detailed member information for it yet. Here
        #we will just update the modifiers and attributes. Also, since all the mods
        #etc. will be overwritten, we don't need to handle replace separately.
        smatch = self.RE_SIG.match(statement)
        result = (None, None, None)
        eresult = None

        if smatch is not None:
            name = smatch.group("name").strip()
            modifiers = smatch.group("modifiers") or []
            codetype = smatch.group("codetype")

            #If the exec is a function, we also may have a type and kind specified.
            if codetype.lower() == "function":
                dtype = smatch.group("type")
                kind = smatch.group("kind")
                if module is None:
                    element.update(name, modifiers, dtype, kind)
                else:
                    eresult = Function(name, modifiers, dtype, kind, module)
            else:
                if module is None:
                    element.update(name, modifiers)
                else:
                    eresult = Subroutine(name, modifiers, module)

            #The parameter sets are actually driven by the body of the executable 
            #rather than the call signature. However, the declarations will be
            #interpreted as members if we don't add the parameters to the ordered
            #list of parameter names. Overwrite that list with the new names.
            params = re.split("[\s,]+", smatch.group("parameters").lower())
            if eresult is None:
                element.paramorder = params
            else:
                eresult.paramorder = params

            result = (eresult, smatch.start(), smatch.end())

        return result
                        
    def parse_line(self, statement, element, mode):
        """Parses the contents of the specified line and adds its representation
        to the specified element (if applicable).

        :arg statement: the lines of code that was added/removed/changed on the 
          element after it had alread been parsed. The lines together form a single
          continuous code statement.
        :arg element: the Subroutine or Function instance to update.
        :arg mode: 'insert', or 'delete'.
        """
        if element.incomplete:
            #We need to check for the end_token so we can close up the incomplete
            #status for the instance.
            if element.end_token in statement:
                element.incomplete = False
                return

        #The line can either be related to an assignment or a dependency since
        #those are the only relevant contents that aren't local variable definitions
        #However if it is a local var definition, we need to process it with any
        #doctags etc that it used. NOTE: doctags are parsed separately by the line
        #parser at a higher level because of the XML.
        self._process_assignments(element, statement, mode)
        self._process_dependencies(element, statement, mode)
        self._parse_members(statement, element, element.paramorder, mode)

    def parse(self, module):
        """Extracts all the subroutine and function definitions from the specified module."""
        #Because of embedded types, we have to examine the entire module for
        #executable definitions.
        self.parse_block(module.refstring, module, module, 0)

        #Now we can set the value of module.contains as the text after the start of
        #the *first* non-embedded executable.
        min_start = len(module.refstring)
        for x in module.executables:
            if module.executables[x].start < min_start:
                min_start = module.executables[x].start

        module.contains = module.refstring[min_start::]

    def parse_block(self, contents, parent, module, depth):
        """Extracts all executable definitions from the specified string and adds
        them to the specified parent."""
        for anexec in self.RE_EXEC.finditer(contents):
            x = self._process_execs(anexec, parent, module)
            parent.executables[x.name.lower()] = x
            if  isinstance(parent, Module) and "public" in x.modifiers:
                parent.publics[x.name.lower()] = 1
            
            #To handle the embedded executables, run this method recursively
            self.parse_block(x.contents, x, module, depth + 1)

        #Now that we have the executables, we can use them to compile a string
        #that includes only documentation *external* to the executable definitions
        #Because we enforce adding the name to 'end subroutine' statements etc.
        #all the embedded executables haven't been parsed yet.
        if len(parent.executables) > 0:
            remove = []
            for x in parent.executables:
                remove.append((parent.executables[x].start, parent.executables[x].end))

            remove.sort(key=lambda tup: tup[0])
            retain = []
            cur_end = 0
            for rem in remove:
                if "\n" in contents[rem[0]+1:rem[1]]:
                    signature = contents[rem[0]+1:rem[1]].index("\n") + 2
                    keep = contents[cur_end:rem[0] + signature]
                    cur_end = rem[1]
                    retain.append(keep)

            #Now we have a string of documentation segments and the signatures they
            #decorate that only applies to the non-embedded subroutines
            docsearch = "".join(retain)
            docblocks = self.docparser.parse_docs(docsearch, parent)

            #Process the decorating documentation for the executables including the
            #parameter definitions.
            for x in parent.executables:
                self._process_docs(parent.executables[x], docblocks, 
                                   parent, module, docsearch)

    def _process_execs(self, execmatch, parent, module):
        """Processes the regex match of an executable from the match object."""
        #Get the matches that must be present for every executable.
        name = execmatch.group("name").strip()
        modifiers = execmatch.group("modifiers")
        if modifiers is None:
            modifiers = []
        else:
            modifiers = re.split(",[ \t]*", modifiers)
        codetype = execmatch.group("codetype")
        params = re.split("[\s,]+", execmatch.group("parameters"))
        contents = execmatch.group("contents")

        #If the exec is a function, we also may have a type and kind specified.
        if codetype.lower() == "function":
            dtype = execmatch.group("type")
            kind = execmatch.group("kind")
            result = Function(name, modifiers, dtype, kind, parent)
        else:
            result = Subroutine(name, modifiers, parent)

        #Set the regex start and end char indices
        result.start, result.end = module.absolute_charindex(execmatch.string, execmatch.start(),
                                                             execmatch.end())
        result.contents = contents

        #Now we can handle the rest which is common to both types of executable
        #Extract a list of local variables
        self._parse_members(contents, result, params)
        if isinstance(result, Function):
            #We need to handle the syntax for defining a function with result(variable)
            #to specify the return type of the function.
            if execmatch.group("result") is not None:
                resvar = execmatch.group("result").lower().split("result(")[1].replace(")", "")
            else:
                resvar = None
            result.update_dtype(resvar)
        
        #Fortran allows lines to be continued using &. The easiest way
        #to deal with this is to remove all of those before processing
        #any of the regular expressions
        decommented = "\n".join([ self._depend_exec_clean(l) for l in contents.split("\n") ])
        cleaned = re.sub("&\s*", "", decommented)

        #Finally, process the dependencies. These are calls to functions and
        #subroutines from within the executable body
        self._process_dependencies(result, cleaned)
        self._process_assignments(result, cleaned)
        
        return result

    def _process_assignments(self, anexec, contents, mode="insert"):
        """Extracts all variable assignments from the body of the executable.

        :arg mode: for real-time update; either 'insert', 'delete' or 'replace'.
        """
        for assign in self.RE_ASSIGN.finditer(contents):
            assignee = assign.group("assignee").strip()
            target = re.split(r"[(%\s]", assignee)[0].lower()

            #We only want to include variables that we know are in the scope of the
            #current executable. This excludes function calls etc. with optional params.
            if target in self._intrinsic:
                continue

            if target in anexec.members or \
               target in anexec.parameters or \
               (isinstance(anexec, Function) and target.lower() == anexec.name.lower()):
                if mode == "insert":
                    anexec.add_assignment(re.split(r"[(\s]", assignee)[0])
                elif mode == "delete":
                    #Remove the first instance of this assignment from the list
                    try:
                        index = element.assignments.index(assign)
                        del element.assignments[index]
                    except ValueError:
                        #We didn't have anything to remove, but python
                        pass                    

    def _process_dependencies(self, anexec, contents, mode="insert"):
        """Extracts a list of subroutines and functions that are called from
        within this executable.

        :arg mode: specifies whether the matches should be added, removed
          or merged into the specified executable.
        """
        #At this point we don't necessarily know which module the executables are
        #in, so we just extract the names. Once all the modules in the library
        #have been parsed, we can do the associations at that level for linking.
        for dmatch in self.RE_DEPEND.finditer(contents):
            isSubroutine = dmatch.group("sub") is not None
            if "!" in dmatch.group("exec"):
                execline = self._depend_exec_clean(dmatch.group("exec"))
            else:
                execline = "(" + dmatch.group("exec").split("!")[0].replace(",", ", ") + ")"

            if not "::" in execline:
                try:
                    dependent = self.nester.parseString(execline).asList()[0]
                except:
                    msg.err("parsing executable dependency call {}".format(anexec.name))
                    msg.gen("\t" + execline)
            
                #Sometimes the parameter passed to a subroutine or function is 
                #itself a function call. These are always the first elements in
                #their nested lists.
                self._process_dependlist(dependent, anexec, isSubroutine, mode)

    def _depend_exec_clean(self, text):
        """Cleans any string constants in the specified dependency text to remove
        embedded ! etc. that break the parsing.
        """
        #First remove the escaped quotes, we will add them back at the end.
        unquoted = text.replace('""', "_FORTPYDQ_").replace("''", "_FORTPYSQ_")
        for cmatch in self.RE_CONST.finditer(unquoted):
            string = cmatch.string[cmatch.start():cmatch.end()]
            newstr = string.replace("!", "_FORTPYEX_")
            unquoted = unquoted.replace(string, newstr)

        requote = unquoted.split("!")[0].replace("_FORTPYDQ_", '""').replace("_FORTPYSQ_", "''")
        result = "(" + requote.replace("_FORTPYEX_", "!").replace(",", ", ") + ")"
        return result

    def _process_dependlist(self, dependlist, anexec, isSubroutine, mode="insert"):
        """Processes a list of nested dependencies recursively."""
        for i in range(len(dependlist)):
            #Since we are looping over all the elements and some will
            #be lists of parameters, we need to skip any items that are lists.
            if isinstance(dependlist[i], list):
                continue

            key = dependlist[i].lower()
            if len(dependlist) > i + 1:
                has_params = isinstance(dependlist[i + 1], list)
            else:
                has_params = False

            #Clean the dependency key to make sure we are only considering valid
            #executable symbols.
            cleanmatch = self.RE_DEPCLEAN.match(key + "(")
            if not cleanmatch:
                continue
            else:
                key = cleanmatch.group("key")

            #We need to handle if and do constructs etc. separately
            if key not in ["then", ",", "", "elseif"] and has_params \
               and not "=" in key and not ">" in key:
                if key in ["if", "do"]:
                    self._process_dependlist(dependlist[i + 1], anexec, False, mode)
                else:
                    #This must be valid call to an executable, add it to the list
                    #with its parameters and then process its parameters list
                    #to see if there are more nested executables
                    if mode == "insert":
                        self._add_dependency(key, dependlist, i, isSubroutine, anexec)
                    elif mode == "delete":
                        #Try and find a dependency already in the executable that 
                        #has the same call signature; then remove it.
                        self._remove_dependency(dependlist, i, isSubroutine, anexec)

                    self._process_dependlist(dependlist[i + 1], anexec, False, mode)

    def _remove_dependency(self, dependlist, i, isSubroutine, anexec):
        """Removes the specified dependency from the executable if it exists
        and matches the call signature."""
        if dependlist[i] in anexec.dependencies:
            all_depends = anexec.dependencies[dependlist[i]]
            if len(all_depends) > 0:
                clean_args = all_depends[0].clean(dependlist[i + 1])
                for idepend in range(len(all_depends)):
                    #Make sure we match across all relevant parameters
                    if (all_depends[idepend].argslist == clean_args
                        and all_depends[idepend].isSubroutine == isSubroutine):
                        del anexec.dependencies[dependlist[i]][idepend]
                        #We only need to delete one, even if there are multiple
                        #identical calls from elsewhere in the body.
                        break

    def _add_dependency(self, key, dependlist, i, isSubroutine, anexec):
        """Determines whether the item in the dependency list is a valid function
        call by excluding local variables and members."""
        #First determine if the reference is to a derived type variable
        lkey = key.lower()
        if "%" in key:
            #Find the type of the base variable and then perform a tree
            #search at the module level to determine if the final reference
            #is a valid executable
            base = key.split("%")[0]
            ftype = None
            if base in anexec.members and anexec.members[base].is_custom:
                ftype = anexec.members[base].kind
            elif base in anexec.parameters and anexec.parameters[base].is_custom:
                ftype = anexec.parameters[base].kind

            if ftype is not None:
                end = anexec.module.type_search(ftype, key)
                if end is not None and isinstance(end, TypeExecutable):
                    #We have to overwrite the key to include the actual name of the type
                    #that is being referenced instead of the local name of its variable.
                    tname = "{}%{}".format(ftype, '%'.join(key.split('%')[1:]))
                    d = Dependency(tname, dependlist[i + 1], isSubroutine, anexec)
                    anexec.add_dependency(d)
                    
        elif lkey not in ["for", "forall", "do"]:
            #This is a straight forward function/subroutine call, make sure that
            #the symbol is not a local variable or parameter, then add it
            if not lkey in anexec.members and not lkey in anexec.parameters \
               and not lkey in self._intrinsic:                
                #One issue with the logic until now is that some one-line statements
                #like "if (cond) call subroutine" don't trigger the subroutine flag
                #of the dependency.
                if dependlist[i-1] == "call":
                    isSubroutine = True
                d = Dependency(dependlist[i], dependlist[i + 1], isSubroutine, anexec)
                anexec.add_dependency(d)
                
    def _process_docs(self, anexec, docblocks, parent, module, docsearch):
        """Associates the docstrings from the docblocks with their parameters."""
        #The documentation for the parameters is stored outside of the executable
        #We need to get hold of them from docblocks from the parent text
        key = "{}.{}".format(parent.name, anexec.name)
        if key in docblocks:
            docs = self.docparser.to_doc(docblocks[key][0], anexec.name)
            anexec.docstart, anexec.docend = (docblocks[key][1], docblocks[key][2])
            self.docparser.process_execdocs(docs, anexec, key)
        #else: the module didn't have any docstrings for this executable...

    def _parse_members(self, contents, anexec, params, mode="insert"):
        """Parses the local variables for the contents of the specified executable."""
        #First get the variables declared in the body of the executable, these can
        #be either locals or parameter declarations.
        members = self.vparser.parse(contents, anexec)

        #If the name matches one in the parameter list, we can connect them
        for param in list(params):
            lparam = param.lower()
            if lparam in members:
                if mode == "insert" and not lparam in anexec.parameters:
                    anexec.add_parameter(members[lparam])
                elif mode == "delete":
                    anexec.remove_parameter(members[lparam])
            
        #The remaining members that aren't in parameters are the local variables
        for key in members:
            if mode == "insert":
                if not key.lower() in anexec.parameters:
                    anexec.members[key] = members[key]
            elif mode == "delete" and key in anexec.members:
                del anexec.members[key]

        #Next we need to get hold of the docstrings for these members
        if mode == "insert":
            memdocs = self.docparser.parse_docs(contents, anexec)
            if anexec.name in memdocs:
                docs = self.docparser.to_doc(memdocs[anexec.name][0], anexec.name)
                self.docparser.process_memberdocs(docs, anexec)

            #Also process the embedded types and executables who may have
            #docstrings just like regular executables/types do.
            self.docparser.process_embedded(memdocs, anexec)

    def _intrinsic_functions(self):
        """Returns a list of fortran intrinsic functions."""
        base = set(['PRESENT', 'ABS ', 'AIMAG ', 'AINT ', 'ANINT ', 'CEILING ', 'CMPLX ', 'CONJG ', 'DBLE ',
                'DIM ', 'DPROD ', 'FLOOR ', 'INT ', 'MAX ', 'MIN ', 'MOD ', 'MODULO ', 'NINT ', 'REAL ',
                'SIGN ', 'ACOS ', 'ASIN ', 'ATAN ', 'ATAN2 ', 'COS ', 'COSH ', 'EXP ', 'LOG ', 'LOG10 ',
                'SIN ', 'SINH ', 'SQRT ', 'TAN ', 'TANH ', 'ACHAR ', 'ADJUSTL ', 'ADJUSTR ', 'CHAR ',
                'IACHAR ', 'ICHAR ', 'INDEX ', 'LEN_TRIM ', 'LGE ', 'LGT ', 'LLE ', 'LLT ', 'REPEAT ',
                'SCAN ', 'TRIM ', 'VERIFY ', 'LEN ', 'KIND ', 'SELECTED_INT_KIND ', 'SELECTED_REAL_KIND ',
                'LOGICAL ', 'DIGITS ', 'EPSILON ', 'HUGE ', 'MAXEXPONENT ', 'MINEXPONENT ', 'PRECISION ',
                'RADIX ', 'RANGE ', 'TINY ', 'BIT_SIZE ', 'BTEST ', 'IAND ', 'IBCLR ', 'IBITS ', 'IBSET ',
                'IEOR ', 'IOR ', 'ISHFT ', 'ISHFTC ', 'NOT ', 'TRANSFER ', 'EXPONENT ', 'FRACTION ',
                'NEAREST ', 'RRSPACING ', 'SCALE ', 'SET_EXPONENT ', 'SPACING ', 'DOT_PRODUCT ', 'MATMUL ',
                'ALL ', 'ANY ', 'COUNT ', 'MAXVAL ', 'MAXVAL ', 'MINVAL ', 'MINVAL ', 'PRODUCT ', 'PRODUCT ',
                'SUM ', 'SUM ', 'ALLOCATED ', 'LBOUND ', 'SHAPE ', 'SIZE ', 'UBOUND ', 'MERGE ', 'PACK ',
                'SPREAD ', 'UNPACK ', 'RESHAPE ', 'CSHIFT ', 'EOSHIFT ', 'TRANSPOSE ', 'MAXLOC ', 'MAXLOC ',
                'MINLOC ', 'MINLOC ', 'ASSOCIATED ', 'NULL ', 'COMMAND_ARGUMENT_COUNT ', 'GET_COMMAND ',
                'CPU_TIME ', 'DATE_AND_TIME ', 'ZONE, VALUES])', 'MVBITS ', 'LEN, TO, TOPOS)', 'RANDOM_NUMBER ',
                'RANDOM_SEED ', 'SYSTEM_CLOCK ', 'AINT ', 'ALOG ', 'ALOG10 ', 'AMAX0 ', 'REAL ', 'AMAX1 ',
                'AMIN0 ', 'REAL ', 'AMIN1 ', 'AMOD ', 'ANINT ', 'CABS ', 'CCOS ', 'CEXP ', 'CHAR ', 'CLOG ',
                'CSIN ', 'CSQRT ', 'DABS ', 'DACOS ', 'DASIN ', 'DATAN ', 'DATAN2 ', 'DCOS ', 'DCOSH ',
                'DDIM ', 'DEXP ', 'DINT ', 'DLOG ', 'DLOG10 ', 'DMAX1 ', 'DMIN1 ', 'DMOD ', 'DNINT ', 'DSIGN ',
                'DSIN ', 'DSINH ', 'DSQRT ', 'DTAN ', 'DTANH ', 'FLOAT ', 'REAL ', 'IABS ', 'IDIM ', 'IDINT ',
                'INT ', 'IDNINT ', 'NINT ', 'IFIX ', 'INDEX ', 'ISIGN ', 'MAX0 ', 'MAX1 ', 'INT ', 'MIN0 ', 
                'MIN1 ', 'INT ', 'SNGL ', 'ISO_C_BINDING', 'C_LOC',
                'C_ASSOCIATED', 'C_PTR_1', 'C_PTR_2', 'C_F_POINTER', 'CAXPY', 'DAXPY', 'SAXPY', 'ZAXPY',
                'CCOPY', 'DCOPY', 'SCOPY', 'ZCOPY', 'CDOTC', 'CDOTU', 'DDOT', 'SDOT', 'ZDOTC', 'ZDOTU', 'CSCAL',
                'DSCAL', 'SSCAL', 'ZSCAL', '-xia', 'DINTERVAL', 'DIVIX', 'INF', 'INTERVAL', 'ISEMPTY', 'MAG',
                'MID', 'MIG', 'NDIGITS', 'QINTERVAL', 'SINTERVAL', 'SUP', 'VDABS', 'VDACOS', 'VDASIN', 'VDATAN',
                'VDATAN2', 'VDCEILING', 'VDCOS', 'VDCOSH', 'VDEXP', 'VDFLOOR', 'VDINF', 'VDINT', 'VDISEMPTY',
                'VDLOG', 'VDLOG10', 'VDMAG', 'VDMID', 'VDMIG', 'VDMOD', 'VDNINT', 'VDSIGN', 'VDSIN', 'VDSINH',
                'VDSQRT', 'VDSUP', 'VDTAN', 'VDTANH', 'VDWID', 'VQABS', 'VQCEILING', 'VQFLOOR', 'VQINF',
                'VQINT', 'VQISEMPTY', 'VQMAG', 'VQMID', 'VQMIG', 'VQNINT', 'VQSUP', 'VQWID', 'VSABS', 'VSACOS',
                'VSASIN', 'VSATAN', 'VSATAN2', 'VSCEILING', 'VSCOS', 'VSCOSH', 'VSEXP', 'VSFLOOR', 'VSINF',
                'VSINT', 'VSISEMPTY', 'VSLOG', 'VSLOG10', 'VSMAG', 'VSMID', 'VSMIG', 'VSMOD', 'VSNINT',
                'VSSIGN', 'VSSIN', 'VSSINH', 'VSSQRT', 'VSSUP', 'VSTAN', 'VSTANH', 'VSWID', 'WID', 'CLOC',
                'COMPL', 'COT', 'CSMG', 'DSHIFTL', 'DSHIFTR', 'EQV', 'FCD', 'IBCHNG', 'ISHA', 'ISHC', 'ISHL',
                'LEADZ', 'LENGTH', 'LOC', 'NEQV', 'POPCNT', 'POPPAR', 'SHIFT', 'SHIFTA', 'SHIFTL',
                'SHIFTR', 'TIMEF', 'UNIT', 'XOR', 'MPI_SIZEOF', 'malloc', 'realloc', 'open', 'close',
                    'allocate', 'deallocate', 'write', 'flush', '.not.present', 'equal', 'do', 'if'])
        return [p.strip() for p in base]
