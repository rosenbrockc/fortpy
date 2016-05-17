from fortpy.elements import Executable, Subroutine
from fortpy.docelements import DocElement
from fortpy.testing.elements import MethodFinder, Assignment
import re

class MethodWriter(object):
    """Class for writing the code to execute a set of dependent pre-requisite
    methods.

    :arg mainid: a string identifying the method to test format: "module.executable".
    :arg parser: an instance of the code parser for inter-module access.
    :arg testgen: the TestGenerator instance that owns this method writer.
    :arg testids: a list of all the testids that we are creating executables for.
    :arg fatal_if_missing: determines whether the framework chokes when it can't
       find a code element in the code parser that is referenced in the docstrings.
    """
    def __init__(self, method, parser, testgen, fatal_if_missing = True, stagedir=None):
        self.mainid = method
        self.parser = parser
        self.method = None
        """The MethodFinder instance for the *current* testid whose driver file is
        being constructed."""
        self.testgen = testgen
        self.finders = {}
        """Dictionary of MethodFinder instances indexed by the testid that they are
        constructed for."""
        self.method_dicts = {}
        """Dictionary of the methods for setting up the unit tests, indexed by the
        testid they are sorted for."""
        self.ordered = {}
        """Dictionary of the method orders indexed by testid."""
        self.globals = {}
        """Dictionary of dictionaries for all globals needed by each of the test
        specifications. Indexed by testid."""
        self.ordered_globals = {}
        """The variables that need to be declared at the top of the program are in
        globals[testid]; the order that they need to appear in is in ordered_globals[testid].
        """
        self.stagedir = stagedir
        """The script-argument staging directory (if any) to override the default staging
        directory in the XML."""
        
        self._uses = {}
        """Dictionary of dependencies for the methods that need be executed in the
        pre-req chain. Keys are the module names and values are the method names."""
        self._tracker = 0
        """_tracker gets used by _order_dependencies() to arrange the methods in the
        correct order so that we don't gen any runtime errors."""

        #Initialize the dictionary of all possible test specifications using a basic
        #version of the method finder.
        finder = MethodFinder(self.mainid, parser, None, None, basic=True)
        if finder.group is not None:
            self.group = finder.group
            """The TestingGroup instance for the executable being unit tested."""
            self.tests = finder.group.tests
            """Returns a dictionary of the tests available to the method writer
            from the underlying XML documentation.
            """
        else:
            self.group = None
            self.tests = {}
        
        for testid in self.tests:
            self.finders[testid] = MethodFinder(self.mainid, parser, None, testid,
                                                fatal_if_missing, True)
            self.method_dicts[testid] = {}
            self.ordered[testid] = []
            self._tracker = 0
            self._order_dependencies(self.finders[testid], testid)

            self.globals[testid] = {}
            self.ordered_globals[testid] = []

    def reset(self, testid, coderoot):
        """Resets the executable writer to work with the new testid under the same
        unit test identifier.
        """
        self.method = self.finders[testid]
        self.method.group.set_finder(self.method)
        self.method.group.coderoot = coderoot

    def copy(self, coderoot, testroot, case, testid, compiler):
        """Performs the copy operation on any assignments that derive their values
        from files.

        :arg coderoot: the full path to the folder that houses all the code files.
        :arg testroot: the full path to the folder that the parent test is being
          performed in.
        :arg case: if a specific case of the same test is being performed, the
          case identifier to use for string formatting.
        :arg compiler: the name of the compiler being used for the unit tests.
        """
        for methodk in self.ordered[testid]:
            method = self.method_dicts[testid][methodk]
            if isinstance(method, Assignment):
                method.copy(coderoot, testroot, case, compiler,
                            self.testgen.xgenerator.libraryroot)

    def setup(self, testid, testroot):
        """Creates any folders for auto-class variables that are needed to save
        output to.
        """
        #With the auto-class system, we also need to make sure that the test targets
        #get to create their directories for complex variables.
        for target in self.finders[testid].test.targets:
            target.init(testroot)
        
    def lines(self, testid, coderoot):
        """Writes all the lines needed to declare and initialize variables and
        call the methods etc. to test the executable.

        :arg testid: the identifier of the test to run.
        """
        self.reset(testid, coderoot)
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
        self._code_once_off(result, "  ", testid)

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
            self._code_repeats(result, "    ", testid)
            self.method.group.code(testid, "after", result, "    ")
            result.append("    !---------CLEANUP SEPARATOR---------------------------")
        else:
            self.method.group.code(testid, "before", result, "  ")
            self._code_repeats(result, "  ", testid)
            self.method.group.code(testid, "after", result, "  ")
        
        #The group code for 'after' would have saved the values of any variables that
        #we needed, so we can deallocate now.
        self._code_assignments(testid, result, "after", "  ")

        if self.tests[testid].constant:
            result.append("  end do")

        #All that's left now is to do the cleanup of open files etc.
        result.append("")
        self.method.group.code(testid, "final", result, "  ")

        #Also, if we were timing this method, save its timing results. We only do
        #timing of the main method being tested.
        self.method.code(result, "final", "  ", testid)

        return "\n".join(result)

    def _check_exec_public(self, anexec):
        """Makes sure that the specified executable is marked as public in the module
        file that houses it.

        :arg anexec: the instance of fortpy.elements.Executable that needs to be public.
        """
        module = anexec.module
        if module is not None:
            if not ("public" in anexec.modifiers or anexec.name.lower() in module.publics):
                target = self.testgen.get_module_target(module.name)
                with open(target) as f:
                    contents = f.readlines()

                #When tests get re-run a bunch of times, we might repeatedly insert the
                #public declaration. Also, once pre-reqs are chained up, there is a 
                #possibility that multiple executables will insert their public declaration
                #just above the first one in unmodified module.

                #One thing we do know is that the public declarations have to appear before
                #members and type definitions in the module.
                insert = True
                keywords = ["type", "interface", "contains"]
                for line in contents[module.public_linenum[0]::]:
                    lline = line.lower().strip()
                    if anexec.name.lower() in lline:
                        insert = False
                        break

                    if (any([re.match(r"\b{}\b".format(b), lline, re.I) for b in keywords])):
                        break
                    if (re.match(r"[\w\s\d]::[\w\s]", lline)):
                        break

                if insert:
                    contents.insert(module.public_linenum[0], "public {}\n".format(anexec.name))
                    with open(target, 'w') as f:
                        f.writelines(contents)
        else:
            raise ValueError("Can't find the module for {};".format(anexec.name) + 
                             "unable to check public declaration for unit testing.")            

    def _code_repeats(self, lines, spacer, testid):
        """Appends lines for calls to methods and variable assignments that
        need to happen *every time* the main method is called."""
        #See the comment in _code_once_off about the arrangement of methods
        #and assignments via ordered dependencies.
        for methodk in self.ordered[testid]:
            method = self.method_dicts[testid][methodk]
            if isinstance(method, MethodFinder) and method.repeats:
                #There is a one-to-one correspondence betweer a MethodFinder and the
                #executable it represents. The only executables that make their way into
                #MethodFinder instances are either pre-reqs or the actual method. We need
                #to make sure that the relevant module files have the executable marked
                #as 'public' so that the unit test executable compiles.
                
                #By now, the files should have been copied over already for the module.
                self._check_exec_public(method.executable)
                method.code(lines, "call", spacer, testid)
            elif isinstance(method, Assignment):
                method.code(lines, "before", spacer)

    def _code_once_off(self, lines, spacer, testid):
        """Appends lines for calls to methods and variable assignments that 
        only need to happen once in a specific order before the main method
        being unit tested is called.
        """
        #Add calls to the methods and functions and any variable assignments
        #Since assignments can rely on the values of variables altered by
        #method calls, we have to reassign the values of variables every time
        #the loop executes.
        for methodk in self.ordered[testid]:
            method = self.method_dicts[testid][methodk]
            if isinstance(method, MethodFinder) and not method.repeats:
                #See the comment in the repeats coding about marking executables as public
                #automatically in the unit test folder.
                self._check_exec_public(method.executable)
                method.code(lines, "call", spacer, testid)
            elif isinstance(method, Assignment):
                method.code(lines, "assign", spacer)

    def _code_init(self, lines, spacer, testid):
        """Appends the code for initializing variables required by the outcomes
        testing and the variable assignments."""
        #Add some whitespace between these declarations and the method calls
        #Add any initialization code for the constant-input test cases.
        lines.append("")
        self.method.group.code(testid, "init", lines, spacer)
        self._code_assignments(testid, lines, "init", spacer)

    def _code_vars(self, lines, spacer, testid):
        """Appends code to declare variables for all the methods and pre-reqs
        needed to run the unit test."""
        #First we write the global definitions from the unit testing specification
        #and then the variable declarations required by the constant-input files.
        for glob in self.globaldefs(testid):
            lines.append(glob.definition())
            sinit = glob.initialization()
            if sinit is not None:
                lines.append(sinit)

        lines.append("")
        self.method.group.code(testid, "vars", lines, spacer)
        
        #Add any variable declarations from the variable assignments that use
        #file contents.
        self._code_assignments(testid, lines, "vars", spacer)

        #We also have declarations for functions that need to have their output values
        #tested by the framework; this also adds the variable for recording the execution
        #time of the unit-tested method.
        for methodk in self.ordered[testid]:
            method = self.method_dicts[testid][methodk]
            if isinstance(method, MethodFinder):
                method.code(lines, "vars", spacer, testid)

    def _code_assignments(self, testid, lines, position, spacer):
        """Appends the relevant lines from all assignments in this unit test based
        on the current position in the code file.

        :arg position: one of ['vars', 'init', 'assign', 'before', 'after'].
        """
        for methodk in self.ordered[testid]:
            method = self.method_dicts[testid][methodk]
            if isinstance(method, Assignment):
                method.code(lines, position, spacer)

    @property
    def autoclass(self):
        """Returns True if any of the test specifications requires auto-class support.
        """
        return any([t.autoclass for t in self.tests.values()])
            
    def uses(self):
        """Gets a list of module.executable dependencies that need to appear in the
        uses clauses of the fortran program."""
        #This dictionary has module names as keys and a list of module executables
        #as the value.
        if len(list(self._uses.keys())) == 0:
            #Check if any of the tests requires the fpy auxiliary module. If so, add
            #it in.
            if self.autoclass:
                self._uses["fpy_auxiliary"] = []
            
            for testid in self.finders:
                for methodk in self.method_dicts[testid]:
                    if not isinstance(self.method_dicts[testid][methodk], Assignment):
                        executable = self.method_dicts[testid][methodk].executable
                        self._process_uses_dependency(executable, executable.parent, self._uses)

                #We also need to add dependencies for variables whose types have custom kinds
                #so that we include the relevant modules. This is because an executable that
                #consumes a derived type doesn't necessarily have that derived type in declared
                #in its *same* module.
                for glob in self.globaldefs(testid):
                    if glob.kind is not None:
                        self._find_global_dependency(glob, self._uses, testid)

                #Next, we check that the output type of the function (if the unit tested executable
                #is a function).
                if type(self.finders[testid].executable).__name__ == "Function":
                    if self.finders[testid].executable.kind is not None:
                        self._find_global_dependency(self.finders[testid].executable, self._uses,
                                                     testid)

        return self._uses

    def _find_global_dependency(self, glob, result, testid):
        """Makes sure that the name of the module and derived type for the kind of
        a global variable is in the dependencies list.

        :arg glob: an instance of GlobalDeclaration to extract kind information from.
        :arg result: the dictionary of module dependencies (keys) and executables, types
          or parameters to use (list of values).
        """
        #We need to make sure that the module who has the kind declaration is
        #included in self._uses. However, some kinds need to be ignored like len=*, digits
        #etc. Variable names cannot have a digit as their first character. We are only
        #interested in variable names.
        if "len" not in glob.kind.lower() and not re.match("\d", glob.kind[0]):
            #Find the module that declares this kind as:
            # 1) derived type
            # 2) parameter/member
            #The good news is that in order for the module being unit tested to compile,
            #it must have a 'use' statement for the derived type. We can just search from
            #that module with a tree find. Also, at the end of the day, it will have to
            #be declared as public, in order to be used in other modules.
            found, foundmod = self.parser.tree_find(glob.kind.lower(), self.finders[testid].module,
                                                    "publics")
            if found is not None:
                if foundmod.name in result:
                    if glob.kind not in result[foundmod.name]:
                        result[foundmod.name].append(glob.kind)
                else:
                    result[foundmod.name] = [ glob.kind ]
            else:
                print((glob.definition()))
                raise ValueError("Cannot find dependency for variable kind '{}'".format(glob.kind))
    
    def _process_uses_dependency(self, method, module, result):
        """Adds the specified dependency to the result dictionary, ignoring duplicates."""
        if module.name in result:
            if method.name not in result[module.name]:
                result[module.name].append(method.name)
        else:
            result[module.name] = [ method.name ]                
    
    def globaldefs(self, testid):
        """Gets a list of all the global variables required by all the methods
        in the entire executable for this testid."""
        if len(self.ordered_globals[testid]) == 0:
            #We haven't extracted and ordered the globals yet. Do it.
            self._get_globals(testid)
            
            #Now we want to examine the list of globals and see if any of them depend
            #on each other via default value assignments or conditional assignments.
            for name in self.globals[testid]:
                glob = self.globals[testid][name]
                if glob.dependency != "" and glob.dependency in self.ordered_globals[testid]:
                    gindex = self.ordered_globals[testid].index(glob.dependency)
                    self.ordered_globals[testid].insert(gindex + 1, name)
                else:
                    self.ordered_globals[testid].append(name)

        return [ self.globals[testid][n] for n in self.ordered_globals[testid] ]

    def _get_globals(self, testid):
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

        for kmethod in self.ordered[testid]:
            #First check its modules global statements
            method = self.method_dicts[testid][kmethod]
            if not isinstance(method, Assignment):
                module = method.module
                if module.name not in checked_modules:
                    mod_group = module.test_group
                    if mod_group is not None:
                        self.globals[testid].update(mod_group.variables)
                    checked_modules.append(module.name)

                #Now, check the executables global statements and regular parameter declarations
                self.globals[testid].update(method.variables)

    def _order_dependencies(self, method, testid):
        """Determines the optimal order to execute the methods in so that we don't
        have any runtime errors."""
        #We need to use a list and a dictionary so that we have the actual objects in
        #a list at the end. If we only keep track of names, we won't know
        #where to find the executable code element instances without going
        #through the list another time.

        #Now we need to recursively take care of all its children
        for anexec in method.methods:
            if anexec.writekey in self.ordered[testid]:
                #This executable is already in the list! We need to see if it
                #is above the current entry or below it. Get the first occurrenc
                #in the list and see if it is less than where we want to insert it.
                existx = self.ordered[testid].index(anexec.writekey)
                if existx > self._tracker:
                    #We want to insert this method in at the index AND remove the other
                    #entry so it doesn't get called twice.
                    self.ordered[testid].remove(anexec.writekey)
                    self._order_dependencies(anexec, testid)
                #elif existx <= index: we don't need to do anything since it will already
                #be executed before it is needed.
            else:
                #We can just add it in at the right place since it isn't in the list yet.
                self._order_dependencies(anexec, testid)

        #The last thing to do is add the executable to the list before its parent
        #since it is a pre-req of its parent
        self.ordered[testid].insert(self._tracker, method.writekey)
        if method.writekey not in self.method_dicts[testid]:
            self.method_dicts[testid][method.writekey] = method
        #We need to make sure that the methods are inserted in the right place across multiple
        #recursive calls to this function. This is the index handling for that.
        self._tracker = self.ordered[testid].index(method.writekey) + 1
