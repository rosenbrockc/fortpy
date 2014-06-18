from fortpy.elements import Executable, Subroutine, DocElement
import re

class GlobalDeclaration(object):
    """Represents a declaration to have a global variable instance available
    for execution of a unit test.

    :attr element: the docstring element representing the <global> tag."""
    def __init__(self, element):
        self.element = element

    @property
    def dependency(self):
        """Returns the name of a variable that this variable needs in order
        to be initialized or assigned a value or ""."""
        if "default" in self.attributes:
            return self.attributes["default"].lower()
        else:
            return ""

    def compare(self, element):
        """Determines whether the specified element defines the same variable
        type and name as this one.

        The two are considered equal if they have the same:
         - name
         - type
         - dimensionality
         - allocatable/pointer modifiers
        ."""
        if self.attributes["name"].lower() != element.attributes["name"].lower():
            return False
        elif not self._kind_check(element):
            return False
        elif not self._modifier_check(element):
            return False
        else:
            return self._dimensions_check(element)
            
    def _kind_check(self, element):
        """Checks whether the kind declaration of the element matches this global."""
        if ("kind" not in self.attributes and "kind" not in element.attributes) or \
           (self.attributes["kind"].lower() == element.attributes["kind"].lower()):
            return True
        else:
            return False

    def _dimensions_check(self, element):
        """Checks whether the dimensions of this declaration match the specified
        element."""
        if not ("dimensions" in self.attributes and "dimensions" in element.attributes):
            return True
        elif "dimensions" in self.attributes and "dimensions" in element.attributes:
            #The dimension text has to match perfectly. If variables names are specified
            #for bounds, we have no way of knowing whether the sizes are the same before
            #runtime. However, we can do some cleanup befor comparing.
            match = True
            selfdim = self.attributes["dimensions"].lower().split(",")
            eldim = element.attributes["dimensions"].lower().split(",")

            i = 0
            #We only need to compare dimensions until one fails
            while match and i < len(selfdim):
                if selfdim[i].strip() != eldim[i].strip():
                    match = False
                i += 1
                
            return match
        else:
            return False

    def _modifier_check(self, element):
        """Checks whether this global var and the specified element match in crucial
        modifier definitions (i.e. pointer, allocatable)."""
        #We need both to have modifiers or not have modifiers
        if not ("modifiers" in self.attributes and "modifiers" in element.attributes):
            return True
        elif "modifiers" in self.attributes and "modifiers" in element.attributes:
            selfmods = self.attributes["modifiers"]
            elmods = element.attributes["modifiers"]
            #At first it seemed like this would be an important check, but a pre-req
            #could theoretically allocate pointer and then it could be passed into
            #a routine that just needed a value.
            #if ("pointer" in selfmods and "pointer" not in elmods) or \
            #   ("pointer" not in selfmods and "pointer" in elmods):
            #    return False
            if ("allocatable" in selfmods and "allocatable" not in elmods) or \
               ("allocatable" not in selfmods and "allocatable" in elmods):
                return False
            else:
                return True
        else:
            return False
        
    @property
    def attributes(self):
        """Returns a dictionary of attributes for this global variable declaration."""
        return self.element.attributes

    def definition(self):
        """Returns the fortran declaration that defines a global variable."""
        result = []
        if "type" not in self.attributes or "name" not in self.attributes:
            print("FATAL: required variable for execution missing some attributes. {}".format(
                self.attributes))
            exit(1)

        result.append(self.attributes["type"])
        if "kind" in self.attributes and self.attributes["kind"] is not None:
            result.append("({})".format(self.attributes["kind"]))

        if "modifiers" in self.attributes:
            mods = self.attributes["modifiers"].split(",")
            smods = self._clean_mods(mods)
        else:
            smods = ""

        if smods != "":
            result.append(", " + smods)
        result.append(" :: ")
        result.append(self.attributes["name"])
        
        if "dimensions" in self.attributes:
            result.append("({})".format(self.attributes["dimensions"]))

        if "default" in self.attributes:
            result.append(" ={}".format(self.attributes["default"]))
            
        return "  " + "".join(result)

    def _clean_mods(self, mods):
        """Returns only those modifiers that belong in a definition statement."""
        result = []
        if "pointer" in mods:
            result.append("pointer")
        if "allocatable" in mods:
            result.append("allocatable")

        return ", ".join(result)

    def initialization(self):
        """Returns the fortran declaration to initialize the global variable's value."""
        #See which of the options for initializing the variable were specified
        #in the docstring.
        if "default" in self.attributes:
            result = None #We handle default values at the definition level.
        elif "mode" in self.attributes:
            result = "call " + self.attributes["name"] + \
                     "%set_mode({})".format(self.attributes["mode"])
        else:
            result = None
        #else: result = "! {} had no initialization value specified".format(self.attributes["name"])

        if result is not None:
            return "  " + result  
        else:
            return result

    def range_check(self):
        """Determines if this variable has a range check specified and
        writes the fortran code to compute the check."""
        if "range" in self.attributes:
            #The fortpy module has an interface to handle any sorts of types we throw at
            #it, so we just need to make the call.
            rangecall = "call fpyin_range({},{},{})"
            #TODO, we need to create the range variables in the program and assign there
            #values if range was specified. 2D arrays need a 2D vector of ranges for
            #each dimension. Need to adjust the fortpy module to accomodate. Call the vars
            #by paramname_min etc.
        else:
            return None

class Assignment(object):
    """Represents an instance value assignment as part of a pre-req chain."""
    def __init__(self, element):
        self.element = element
        self.methods = []

    @property
    def name(self):
        """Returns the name of the variable that needs its value changed."""
        if "name" in self.element.attributes:
            return self.element.attributes["name"]
        else:
            return ""

    @property
    def attributes(self):
        """Provides one-level-up access to the XML elements attributes collection."""
        return self.element.attributes
        
    def code(self):
        """Generates the code to accomplish the variable assignment."""
        if "conditionals" in self.attributes:
            return self._code_conditional()
        elif "value" in self.attributes:
            return "  {} = {}".format(self.name, self.attributes["value"])
        else:
            return ""

    def _code_conditional(self):
        """Returns code to complete a conditional assignment."""
        #Generate some if statements to handle all the conditions and values.
        ifs = self.attributes["conditionals"].split("|")
        template = "  {}if ({}) then"
        lhs = "{} = ".format(self.attributes["name"])

        first = ifs[0].strip().split(":")
        allifs = [ template.format("", first[0]) ]
        allifs.append("    " + lhs + first[1])
        

        for anif in ifs[1::]:
            clean = anif.strip().split(":")
            if not "else:" in clean:
                #This must be an elseif clause. The template already has the
                #if part, so just add in the else.
                allifs.append(template.format("else", clean[0]))
            else:
                allifs.append("else")        
            allifs.append("    " + lhs + clean[1])

        allifs.append("  endif")
        return "\n".join(allifs)

class MethodFinder(object):
    """Class for recursively finding methods for executing pre-req methods
    that takes recursive dependencies into account.

    :arg identifier: the module.executable that this method finder represents
    :arg parser: the code parser instance for inter-module access.
    :arg element: the DocElement that specified this method as a pre-req.
    :arg fatal_if_missing: specified whether the framework chokes if it can't
      find a method that is referenced as a pre-req.
    """
    def __init__(self, identifier, parser, element, fatal_if_missing = True):
        self.identifiers = identifier.split(".")
        self.name = identifier
        self.methods = []
        self._parser = parser
        self.element = element

        self._module = None
        self._fatal = fatal_if_missing
        self._tests = []

        self.executable = self._find_executable()

        #Recursively get all the pre-requisite for this method
        if element is None or \
           not ("terminate" in self.element.attributes and self.element.attributes["terminate"] == "true"):
            self._get_prereqs()
        
    @property
    def attributes(self):
        """Returns the dictionary of attributes from the underlying DocElement."""
        if self.element is not None:
            return self.element.attributes
        else:
            return {}

    def _get_prereqs(self):
        """Compiles a list of subroutines that need to run. Handles recursive calls."""
        for test in self._tests:
            if test.doctype == "prereq" and "method" in test.attributes:
                identify = test.attributes["method"]
                self.methods.append(MethodFinder(identify, self._parser, test, self._fatal))

            elif test.doctype == "instance":
                if "method" in test.attributes:
                    #This is slightly more complicated because the name of the instance
                    #is the name of a parameter, NOT a type. We need to look up the type
                    #from the list of parameters in this executable
                    param = test.attributes["name"]
                    if param in self.executable.parameters:
                        ptype = self.executable.parameters[param]
                        identify = "{}.{}".format(ptype, test.attributes["method"])
                        self.methods.append(MethodFinder(identify, self._parser, test, self._fatal))
                    elif self._fatal:
                        print("FATAL: the parameter {} for instance method {} execution missing.".format(
                            param, test.attributes["method"]))
                elif "value" in test.attributes or "conditionals" in test.attributes:
                    #This is a variable assignment. Sometimes variable values are changed
                    #between calls to functions. The order in which prereq/instance elements
                    #appear defines when these assignments take place. We will treat a var
                    #value assignment as a method to simplify the implementation.
                    self.methods.append(Assignment(test))
                    
                
    def _find_executable(self):
        """Finds the executable code element from the code parser instance."""
        #First, we need to determine which module the identifier points to.
        #See if the module is in the code parser
        result = None

        if self.module is not None:
            #We want to use the code parsers tree search to get the executable
            #If the identifier refers to a type method, we need to use the
            #type search
            if "%" in self.identifiers[1]:
                tbase = self.identifiers[1].split("%")[0].strip()
                element = self._parser.type_search(tbase, self.identifiers[1], self._module)
                if isinstance(element, Executable):
                    result = element                   
            else:
                #This is just a standard exectable inside the module
                if self.identifiers[1] in self._module.executables:
                    result = self._module.executables[self.identifiers[1]]
                elif self._fatal:
                    print("FATAL: a specified executable was not found in the module: {}".format(
                        self.identifiers[1]))
                    
        #If we have a result, get a pointer to the tests so we don't have to go one level
        #deep all the time.
        if result is not None:
            self._tests = result.tests

        return result

    @property
    def module(self):
        """Returns the module name for this method."""
        if self._module is None:
            if not self.identifiers[0] in self._parser.modules:
                #Try a dependency search for the specific module
                self._parser.load_dependency(self.identifiers[0], True, True, False)
                #If the module was found, we can set the value of the local variable
                #otherwise we will end the program by default
                if self.identifiers[0] in self._parser.modules:
                    self._module = self._parser.modules[self.identifiers[0]]
                elif self._fatal:
                    print("FATAL: the module for a pre-requisite method could not be located: {}".format(
                        self.identifiers[0]))
                    exit -1
                    
            if self.identifiers[0] in self._parser.modules:
                self._module = self._parser.modules[self.identifiers[0]]
        
        return self._module        

class MethodWriter(object):
    """Class for writing the code to execute a set of dependent pre-requisite
    methods.

    :arg method: a string identifying the method to test module.executable.
    :arg parser: an instance of the code parser for inter-module access.
    :arg fatal_if_missing: determines whether the framework chokes when it can't
       find a code element in the code parser that is referenced in the docstrings.
    """
    def __init__(self, method, parser, fatal_if_missing = True):
        self.method = MethodFinder(method, parser, None, fatal_if_missing)
        self._ordered = []
        self._method_dict = {}
        self._globals = {}
        self._ordered_globals = []
        self._uses = {}

        self._tracker = 0
        self._order_dependencies(self.method)

    def lines(self):
        """Writes all the lines needed to declare and initialize variables and
        call the methods etc. to test the executable."""
        #First we write the global definitions, and then the initializations.
        result = []
        result.append("")

        for glob in self.globaldefs:
            result.append(glob.definition())
            sinit = glob.initialization()
            if sinit is not None:
                result.append(sinit)

        #Add some whitespace between these declarations and the method calls
        result.append("")

        #Add calls to the methods and functions and any variable assignments
        #in between.
        for methodk in self._ordered:
            method = self._method_dict[methodk]
            result.append(self._write_method(method))

        result.append("")

        return "\n".join(result)

    def uses(self):
        """Gets a list of module.executable dependencies that need to appear in the
        uses clauses of the fortran program."""
        #This dictionary has module names as keys and a list of module executables
        #as the value.
        if len(list(self._uses.keys())) == 0:     
            for methodk in self._method_dict:
                if not isinstance(self._method_dict[methodk], Assignment):
                    executable = self._method_dict[methodk].executable
                    self._process_uses_dependency(executable, executable.parent, self._uses)

        return self._uses

    def _process_uses_dependency(self, method, module, result):
        """Adds the specified dependency to the result dictionary, ignoring duplicates."""
        if module.name in result:
            if method.name not in result[module.name]:
                result[module.name].append(method.name)
        else:
            result[module.name] = [ method.name ]                
        
    def _write_method(self, method):
        """Returns the call to a subroutine or function or assignment as a string."""
        if isinstance(method, Assignment):
            #This is just a simple variable value assignment.
            return method.code()
        else:
            #This is either a call to a subroutine or a function. We need to make
            #sure that we handle any mapping tags or call tags in the docstrings
            prefix = "call " if isinstance(method.executable, Subroutine) else ""
            spacing = len(list(prefix)) + len(list(method.executable.name)) + 1

            if "paramlist" in method.attributes:
                return "  {}{}({})".format(prefix, method.executable.name,
                                           self._present_params(re.split(",\s*", method.attributes["paramlist"]),
                                                                spacing))
            else:
                #We need to construct it manually from the parameters in the method
                #Create a dictionary of name mappings from the testing docstrings
                mappings = {}
                for test in method.executable.tests:
                    if test.doctype == "mapping" and "name" in test.attributes and \
                       "target" in test.attributes:
                        mappings[test.attributes["target"]] = test.attributes["name"]

                #Now we can construct the actual list of parameters to use in the call
                #and substitute mappings where appropriate.
                calllist = []
                for param in method.executable.ordered_parameters:
                    if param in mappings:
                        calllist.append(mappings[param])
                    else:
                        calllist.append(param.name)
                return "  {}{}({})".format(prefix, method.executable.name,
                                           self._present_params(calllist, spacing))
        
    def _present_params(self, paramlist, spacing = 0):
        """Creates the (paramlist) for a method call formatted nicely for calls
        with lots of parameters."""
        #The +2 is spacing is for the tab indent at the start of the line.
        #The +3 is for indent and the extra parenthesis at the start of the call.
        line = []
        length = 0
        result = []
        for param in paramlist:
            extra = len(list(param))
            if length + extra + 2 + spacing > 90:
                result.append(", ".join(line) + ", &")
                line = [ param ]
                length = extra + 2
            else:
                line.append(param)
                length += extra + 2

        #Add on the remaining bits of the line
        result.append(", ".join(line))

        return "\n{}".format(" ".join([ "" for i in range(spacing + 3)])).join(result)

    @property
    def globaldefs(self):
        """Gets a list of all the global variables required by all the methods
        in the entire executable."""
        if len(self._ordered_globals) == 0:
            #We haven't extracted and ordered the globals yet. Do it.
            self._get_globals()
            
            #Now we want to examine the list of globals and see if any of them depend
            #on each other via default value assignments or conditional assignments.
            for name in self._globals:
                glob = self._globals[name]
                if glob.dependency != "" and glob.dependency in self._ordered_globals:
                    gindex = self._ordered_globals.index(glob.dependency)
                    self._ordered_globals.insert(gindex + 1, name)
                else:
                    self._ordered_globals.append(name)

        return [ self._globals[n] for n in self._ordered_globals ]

    def _get_globals(self):
        """Compiles a list of global variables declarations to add to the fortran program."""
        #Parameters marked as "regular" get defined using the assumed requirements
        #for their parameter definitions in the methods. All others have to have a
        #global keyword to define them inside the testing docblock.

        #If two globals are defined with the same name but different types, the 
        #framework should choke.

        #We want to sweep over all methods and pull the global tags out of their modules
        #and executable testing groups. Each time we need to make sure we haven't added
        #it already. The order in which they are declared and initialized could matter 
        #(since the initialization could depend on an existing variable).
        checked_modules = []

        for kmethod in self._ordered:
            #First check its modules global statements
            method = self._method_dict[kmethod]
            if not isinstance(method, Assignment):
                module = method.module
                if module.name not in checked_modules:
                    self._get_test_globals(module.tests)
                    checked_modules.append(module.name)

                #Now, check the executables global statements and regular parameter declarations
                self._get_test_globals(method.executable.tests)
                self._get_param_globals(method.executable)

    def _get_param_globals(self, executable):
        """Extracts global declarations from the parameter regular="true" attribute
        of paramater decorations."""
        for name in executable.parameters:
            param = executable.parameters[name]
            for doc in param.docstring:
                if doc.doctype == "parameter" and \
                   "regular" in doc.attributes and doc.attributes["regular"] == "true":
                    glob = self._global_from_param(executable, name)
                    if glob is not None:
                        self._global_add(name, glob)

    def _global_from_param(self, executable, name):
        """Creates a global DocElement from the specified executable's parameter "name"."""
        if name in executable.parameters:
            param = executable.parameters[name]
            result = DocElement(None, None)
            result.doctype = "AUTOPARAM"
            result.attributes["name"] = name
            result.attributes["type"] = param.dtype
            self._global_clean_param(result, "kind", param.kind)
            self._global_clean_param(result, "modifiers", ", ".join(param.modifiers))
            self._global_clean_param(result, "dimensions", param.dimension)
            self._global_clean_param(result, "default", param.default)
            return result
        else:
            return None

    def _global_clean_param(self, result, name, value):
        """Adds the specified parameter name and value to the result if the value is not None."""
        if value is not None:
            result.attributes[name] = value

    def _get_test_globals(self, tests):
        """Extracts a list of the global declarations from the specified tests
        and adds them to the writer's globals dictionary if they are unique."""
        for doc in tests:
            if doc.doctype == "global" and "name" in doc.attributes:
                name = doc.attributes["name"]
                self._global_add(name, doc)

    def _global_add(self, name, doc):
        """Adds a global declaration for the specified DocElement if it isn't already
        represented in the list."""
        if name not in self._globals:
            self._globals[name] = GlobalDeclaration(doc)
        else:
            #We need to make sure that it is unique compared to the existing
            #one. If it isn't, stop the execution. If it is, we don't need
            #to add it again.
            existing = self._globals[name]
            if not existing.compare(doc):
                print("FATAL: variables in the calling namespace have the same name," + \
                    " but different types: \n{}{}".format(existing.element, doc))
                exit(1)
    
    def _order_dependencies(self, method):
        """Determines the optimal order to execute the methods in so that we don't
        have any runtime errors."""
        #We need to use a list and a dictionary so that we have the actual objects in
        #a list at the end. If we only keep track of names, we won't know
        #where to find the executable code element instances without going
        #through the list another time.

        #Now we need to recursively take care of all its children
        for anexec in method.methods:
            if anexec.name in self._ordered:
                #This executable is already in the list! We need to see if it
                #is above the current entry or below it. Get the first occurrenc
                #in the list and see if it is less than where we want to insert it.
                existx = self._ordered.index(anexec.name)
                if existx > index:
                    #We want to insert this method in at the index AND remove the other
                    #entry so it doesn't get called twice.
                    self._ordered.remove(anexec.name)
                    self._order_dependencies(anexec)
                #elif existx <= index: we don't need to do anything since it will already
                #be executed before it is needed.
            else:
                #We can just add it in at the right place since it isn't in the list yet.
                self._order_dependencies(anexec)

        #The last thing to do is add the executable to the list before its parent
        #since it is a pre-req of its parent
        self._ordered.insert(self._tracker, method.name)
        if method.name not in self._method_dict:
            self._method_dict[method.name] = method
        #We need to make sure that the methods are inserted in the right place across multiple
        #recursive calls to this function. This is the index handling for that.
        self._tracker = self._ordered.index(method.name) + 1
