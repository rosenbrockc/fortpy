from fortpy.elements import Executable, Subroutine
from fortpy.docelements import DocElement
from fortpy.testing.elements import MethodFinder, Assignment
import re

class MethodWriter(object):
    """Class for writing the code to execute a set of dependent pre-requisite
    methods.

    :arg method: a string identifying the method to test format: "module.executable".
    :arg parser: an instance of the code parser for inter-module access.
    :arg fatal_if_missing: determines whether the framework chokes when it can't
       find a code element in the code parser that is referenced in the docstrings.
    """
    def __init__(self, method, parser, fatal_if_missing = True):
        self.method = MethodFinder(method, parser, None, fatal_if_missing)

        #Dictionary of dependencies for the methods that need be executed in the
        #pre-req chain. Keys are the module names and values are the method names.
        self._uses = {}
        #The pre-req methods required for the main method being unit tested are all
        #in the _method_dict; the order they need to be run in is in _ordered.
        self._ordered = []
        self._method_dict = {}
        #The variables that need to be declared at the top of the program are in
        #_globals; the order that they need to appear in is in _ordered_globals.
        self._ordered_globals = []
        self._globals = {}
        #_tracker gets used by _order_dependencies() to arrange the methods in the
        #correct order so that we don't gen any runtime errors.
        self._tracker = 0
        self._order_dependencies(self.method)

    @property
    def tests(self):
        """Returns a dictionary of the tests available to the method writer
        from the underlying XML documentation.
        """
        if self.method.group is not None:
            return self.method.group.tests
        else:
            return {}

    def copy(self, coderoot, testroot, case):
        """Performs the copy operation on any assignments that derive their values
        from files.

        :arg coderoot: the full path to the folder that houses all the code files.
        :arg testroot: the full path to the folder that the parent test is being
          performed in.
        :arg case: if a specific case of the same test is being performed, the
          case identifier to use for string formatting."""
        for methodk in self._ordered:
            method = self._method_dict[methodk]
            if isinstance(method, Assignment):
                method.copy(coderoot, testroot, case)
        
    def lines(self, testid):
        """Writes all the lines needed to declare and initialize variables and
        call the methods etc. to test the executable.

        :arg testid: the identifier of the test to run.
        """
        result = []
        result.append("")
        self._code_vars(result, "  ", testid)

        #If we have to repeat the main method calls multiple times, we need a variable
        #to loop over to simplify things.
        if self.tests[testid].constant:
            result.append("  integer :: fpy_repeat")

        self._code_init(result, "  ", testid)

        #Before we can run the tests, we need to validate the input selection for
        #constant-input testing. These make sure that if the input values come from
        #multiple files, that we have the same number of input values from each file.
        self.method.group.code(testid, "validate", result, "  ")
        self._code_once_off(result, "  ")

        #If the test has any constant-input specifcations, we need to run a loop
        #over each set of values in the input files.
        if self.tests[testid].constant:
            repeats = self.tests[testid].repeats
            result.append("  do fpy_repeat = 1, {}".format(repeats))

        #The group 'before' code must run once before each call to the main methods
        #It is not sensitive to being inside or outside the loop. Placing it here
        #makes sure it gets appended correctly either may.

        #After the method has run, we need to do some target saving of any varibles
        #that the developer finds interesting.
        if self.tests[testid].constant:
            self.method.group.code(testid, "before", result, "    ")
            self._code_repeats(result, "    ")
            self.method.group.code(testid, "after", result, "    ")
            result.append("    !---------CLEANUP SEPARATOR---------------------------")
        else:
            self.method.group.code(testid, "before", result, "  ")
            self._code_repeats(result, "  ")
            self.method.group.code(testid, "after", result, "  ")
            result.append("  !---------CLEANUP SEPARATOR---------------------------")

        #The group code for 'after' would have saved the values of any variables that
        #we needed, so we can deallocate now.
        self._code_assignments(result, "after", "  ")

        if self.tests[testid].constant:
            result.append("  end do")

        #All that's left now is to do the cleanup of open files etc.
        result.append("")
        self.method.group.code(testid, "final", result, "  ")

        return "\n".join(result)

    def _code_repeats(self, lines, spacer):
        """Appends lines for calls to methods and variable assignments that
        need to happen *every time* the main method is called."""
        #See the comment in _code_once_off about the arrangement of methods
        #and assignments via ordered dependencies.
        for methodk in self._ordered:
            method = self._method_dict[methodk]
            if isinstance(method, MethodFinder) and method.repeats:
                method.code(lines, "call", spacer)
            elif isinstance(method, Assignment):
                method.code(lines, "before", spacer)

    def _code_once_off(self, lines, spacer):
        """Appends lines for calls to methods and variable assignments that 
        only need to happen once in a specific order before the main method
        being unit tested is called.
        """
        #Add calls to the methods and functions and any variable assignments
        #Since assignments can rely on the values of variables altered by
        #method calls, we have to reassign the values of variables every time
        #the loop executes.
        for methodk in self._ordered:
            method = self._method_dict[methodk]
            if isinstance(method, MethodFinder) and not method.repeats:
                method.code(lines, "call", spacer)
            elif isinstance(method, Assignment):
                method.code(lines, "assign", spacer)

    def _code_init(self, lines, spacer, testid):
        """Appends the code for initializing variables required by the outcomes
        testing and the variable assignments."""
        #Add some whitespace between these declarations and the method calls
        #Add any initialization code for the constant-input test cases.
        lines.append("")
        self.method.group.code(testid, "init", lines, spacer)
        self._code_assignments(lines, "init", spacer)

    def _code_vars(self, lines, spacer, testid):
        """Appends code to declare variables for all the methods and pre-reqs
        needed to run the unit test."""
        #First we write the global definitions from the unit testing specification
        #and then the variable declarations required by the constant-input files.
        for glob in self.globaldefs:
            lines.append(glob.definition())
            sinit = glob.initialization()
            if sinit is not None:
                lines.append(sinit)

        lines.append("")
        self.method.group.code(testid, "vars", lines, spacer)
        
        #Add any variable declarations from the variable assignments that use
        #file contents.
        self._code_assignments(lines, "vars", spacer)

        #We also have declarations for functions that need to have their output values
        #tested by the framework.
        for methodk in self._ordered:
            method = self._method_dict[methodk]
            if isinstance(method, MethodFinder):
                method.code(lines, "vars", spacer)

    def _code_assignments(self, lines, position, spacer):
        """Appends the relevant lines from all assignments in this unit test based
        on the current position in the code file.

        :arg position: one of ['vars', 'init', 'assign', 'before', 'after'].
        """
        for methodk in self._ordered:
            method = self._method_dict[methodk]
            if isinstance(method, Assignment):
                method.code(lines, position, spacer)

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
                    mod_group = module.test_group
                    if mod_group is not None:
                        self._globals.update(module.test_group.variables)
                    checked_modules.append(module.name)

                #Now, check the executables global statements and regular parameter declarations
                self._globals.update(method.variables)

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
