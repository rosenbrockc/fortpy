from .generator import TestGenerator
from ..code import CodeParser, secondsToStr
from os import system, path, mkdir
from time import clock
from shutil import copy, copyfile
from fortpy.code import secondsToStr
from .comparer import FileComparer
from fortpy.testing.results import print_compare_result

class ExecutionResult(object):
    """The result of running the executable on the system, NOT the
    result of the unit test file comparisons.

    :arg execfolder: the folder in which the executable ran.
    :arg exitcode: the system exit code after the process terminated.
    :arg runtime: a python datetime for the execution run time.
    :arg tester: an OutcomeTester that can test the outcome of this
      execution.
    :arg case: the case identifier if this was part of a series of
      cases run for the same outcome.
    """
    def __init__(self, folder, exitcode, runtime, tester, case = None):
        self.folder = folder
        self.code = exitcode
        self.runtime = runtime
        self.tester = tester
        self.case = case

    def _compare_time(self):
        """Determines if the runtime for this execution exceeded the
        maximum in the tester."""
        if self.tester.runtime is not None and \
           self.tester.unit is not None:
            #Get the runtime in a nice-to-compare format 0:00:00.000
            stime = secondsToStr(self.runtime)
            #we have limits specified as minutes or hours in a range
            #No unit test should/can take more than 24 hours.
            utime = [ int(u) for u in stime.split(":") ]
            if self.tester.unit == "h":
                return self._compare_list_time(self.tester.runtime, utime[0])
            elif self.tester.unit == "m":
                return self._compare_list_time(self.tester.runtime, utime[1])  
        else:
            return True

    def _compare_list_time(self, bounds, actual):
        """Sees if the actual time integer value lies between the
        integer list bounds [min,max]."""
        return actual >= bounds[0] and actual <= bounds[1]

    def test(self, caseid, uresult):
        """Tests the outcome of this executable using its tester attribute."""
        #The first thing to check is the contents of the files
        self.tester.test(caseid, self, uresult)
        #All that's left is to check the runtime against its max.
        if not self._compare_time():
            uresult.overtimes[caseid] = (self.runtime, self.tester.runtime)        

class ValueCompareResult(object):
    """The result of comparing the contents of a file to an explicit
    value hard-coded in the testing XML tags."""
    def __init__(self, path, value):
        self.path = path
        self.value = value
        self.equal = self._compare_values()

    def _compare_values(self):
        """Tests the equality of the file contents in path to the
        explicit value specified."""
        result = False

        #We only deal with the first line in the file, anything more
        #complicated should be handled with file templates.
        with open(self.path) as f:
            line = f.readline()

        try:
            if eval(line) == eval(self.value):
                result = True
        except ValueError:
            result = False
        
        return result

class OutcomeTester(object):
    """Performs outcomes tests for a single 'outcome' DocString.

    :arg doc: the DocString for the 'outcome' tag.
    :arg codefolder: the path to the folder that has all the modules from which
      the unit tests were built.
    :arg comparer: an instance of FileComparer to compare output files
      to the model ones.
    """
    def __init__(self, doc, codefolder, comparer, verbose):
        self.doc = doc
        self.comparer = comparer
        self.verbose = verbose

        #Extract the paths to the source input and output folders
        self.source = doc.attributes["sourcepath"].replace(".", codefolder)
        if "comparepath" in doc.attributes:
            self.compare = doc.attributes["comparepath"].replace(".", codefolder)
        else:
            self.compare = self.source
        
        #Extract the targets that need to be compared and the files to compare
        #them to. Since multiple targets can be specified, we need to look at 
        #each one to see if it is a file or a variable target

        #The corresponding entries in the compare attribute must be either files
        #or values too.
        if "target" not in doc.attributes or "compare" not in doc.attributes:
            print("WARNING: either targets or model output files were missing for {}".format(doc))
            self.targets = []
            self.models = []
        else:
            self.targets = doc.attributes["target"].split(",")
            self.models = doc.attributes["compare"].split("?")

        #Get all the other optional attributes
        self._parse_xml()

    def _parse_xml(self):
        """Extracts XML attributes related to comparison modes and tolerances
        that are used to analyze the file comparison."""
        #Extract the doc attributes related to the comparison performance
        if "mode" in self.doc.attributes:
            self.mode = self.doc.attributes["mode"]
        else:
            self.mode = "default"

        if "runtime" in self.doc.attributes:
            runtime = list(self.doc.attributes["runtime"])
            self.unit = runtime.pop()
            self.runtime = [ int(t) for t in "".join(runtime).split("-") ]
        else:
            self.runtime = None
            self.unit = None

        if "tolerance" in self.doc.attributes:
            self.tolerance = float(self.doc.attributes["tolerance"])
        else:
            self.tolerance = 1
        
    def test(self, caseid, xresult, uresult):
        """Checks the output of the execution for the specified outcome doctag.

        :arg xresult: an ExecutionResult instance that has information about the
          execution whose outcome is being tested.
        :arg uresult: the unittest's TestResult that holds information from all
          executions that were run for different cases.
        """        
        #Now compare the files with their model outputs
        compresults = []
        for i in range(len(self.models)):
            paths = self._run_get_paths(self.targets[i], self.models[i], xresult.folder,
                                        self.compare, xresult.case)
            #The third entry in the tuple identifies it as an explicit value
            #comparison as opposed to a file comparison.
            if paths[2]:
                result = self._run_compare_file(paths[0], paths[1])
                if result is None or result.common_match < self.tolerance:
                    self._add_failure(caseid, result, uresult)
            else:
                result = self._run_compare_var(paths[0], eval(paths[1]))
                if result is None or not result.equal:
                    self._add_failure(caseid, result, uresult)
            
            compresults.append(result)

        #Add the outcome to the results dictionary
        uresult.outcomes[caseid] = compresults

    def _add_failure(self, caseid, result, uresult):
        """Adds the specified result as a failure for the caseid."""
        if caseid in uresult.failures:
            uresult.failures[caseid].append(result)
        else:
            uresult.failures[caseid] = [ result ]

    def _run_compare_var(self, exe, value, doc):
        """Compares the value of a variable written to an output file with an
        explicit value."""
        return ValueCompareResult(exe, value)

    def _run_compare_file(self, exe, model):
        """Compares an output file from an executable with its model output using
        settings in the doctag.

        :arg exe: the full path to the output file from the unit test.
        :arg model: the full path to the model output file to compare.
        """
        #First get the comparison results, then analyze them using the tolerances
        #etc. from the doctag.
        result = self.comparer.compare(exe, model, self.mode)
        #Write the results out to file as a record. If the result is none create
        #a file that says so. The default file name is the output file name with
        #an extra extension of .compare
        resultpath = exe + ".compare"
        with open(resultpath, "w") as f:  
            if result is not None:          
                f.write(print_compare_result(result, self.verbose))
            else:
                f.write("The result comparison failed. Check the unit test console output.")
            
        return result

    def _run_get_paths(self, target, compare, exefolder, source, case):
        """Gets the file paths to the executable output and model output
        files that need to be compared. If the model output is not a file
        the third entry in the tuple will be False.

        :arg target: the target file or variable that needs to be compared.
        :arg compare: the name of the file in the source folder to which
          the target will be compared.
        :arg execfolder: the folder in which the executable ran, also where
          the output files from the executable reside.
        :arg source: the path to the source folder where the model outputs are.
        """
        #Keep track of whether we are comparing a variable value file to a hard-coded
        #value or to another output file.
        isvar = False

        if target[0] == ".":
            exe = path.join(exefolder, target[2::]).format(case)
        else:
            #This is a variable, get the auto-generated filename
            filename = target.replace("%", ".")
            exe = path.join(exefolder, filename)

        if compare[0] == ".":
            mod = path.join(source, compare[2::].format(case))
        else:
            isvar = True
            mod = compare

        return (exe, mod, not isvar)     

class TestResult(object):
    """Represents a set of unit test results.

    :arg identifier: the module.method identifier for this unit test.
    :attr cases: a dictionary of ExecutionResult objects with detail
      on how the system execution went for each case. Key is caseId
    :attr outcomes: a dictionary of CompareResult objects with detail
      on how similar the new output files are to the model ones. Keys
      are the same caseIds used in self.cases.
    :attr compiled: specifies whether the executable associated with
      this identifier compiled successfully.
    """
    def __init__(self, identifier):
        self.identifier = identifier
        self.cases = {}
        self.outcomes = {}
        self.compiled = None
        self.warnings = []
        self.failures = {}
        self.overtimes = {}

    @property
    def percent(self):
        """Returns the percent success of the entire unit test."""
        total = 0
        common = 0

        for caseid in self.outcomes:
            case = self.outcomes[caseid]
            for result in case:
                if isinstance(result, ValueCompareResult):
                    if result.equal:
                        total += 1
                else:
                    if result is not None:
                        total += result.common_match

        return float(total) / (len(list(self.cases.keys())) + self.failure_count)
            
    @property
    def failure_count(self):
        """Gets the number of outcomes that failed because of tolerance
        or misformatted or no output files."""
        total = 0
        for caseid in self.failures:
            for result in self.failures[caseid]:
                total += 1

        return total            

class UnitTester(object):
    """Performs automatic unit testing of individual subroutines and
    functions based on XML doctags found in the code files.

    :arg libraryroot: the path to folder in which to stage the tests.
    """
    def __init__(self, libraryroot, verbose = False, compare_templates = None, 
                 fortpy_templates = "~/pythonpkg/fortpy/templates", rerun = False):
        self.libraryroot = path.abspath(libraryroot)
        self.parser = CodeParser()
        self.parser.verbose = verbose
        self.fortpy_templates = path.abspath(fortpy_templates)
        self.tgenerator = TestGenerator(self.parser, self.libraryroot, self.fortpy_templates,
                                        rerun)
        
        #A flag to track whether the generator has already written
        #the executables.
        self._written = False
        #The user's raw value for the compare_templates directory
        self._compare_templates = compare_templates

    def writeall(self, codefolder):
        """Writes all the unit test executables that are new or modified
        as well as the makefiles for all subroutines in all modules."""
        #The test generator already loops over all modules in the code
        #parser and does all the heavy-lifting. We need to pre-load any
        #modules that we are interested in testing. Only those loaded
        #when write() is first called will have their tests processed.
        self._codefolder = path.abspath(codefolder)
        if self._compare_templates is not None:
            self.compare_templates = path.abspath(self._compare_templates)
        else:
            self.compare_templates = path.join(self._codefolder, "templates/")

        #We will load all the modules in the code folder specified and
        #then run the test generator.
        files = {}
        self.parser.scan_path(self._codefolder, files)
        for f in files:
            filepath = files[f]
            self.parser.parse(filepath, True, True)

        #Now that we have loaded all the codefiles in the path, we can
        #generate the unittest executables
        self.tgenerator.write()
        self._written = True

    def runall(self):
        """Compiles and runs each new or modified unit test executable.
        After the run is complete, the outcomes are checked for consistency."""
        if self._written:
            #We will likely need a file comparer for the outcomes testing
            self.comparer = FileComparer(self.compare_templates)
            #Run them each individually and return a dictionary of all the
            #test results
            result = {}
            for identifier in self.tgenerator.tests_to_run:
                result[identifier] = self._run_single(identifier)

            return result
        else:
            print("WARNING: you can't run tests until the executables have been written. Exiting.")
            return None
       
    def _run_single(self, identifier):
        """Runs all unit test cases for the specified identifier."""
        #Just initialize a result object and populate its properties
        #using the various _run_* methods.
        result = TestResult(identifier)
        result.compiled = self._run_compile(identifier)
        if result.compiled:
            self._run_exec(identifier, result)

        return result

    def _run_compile(self, identifier):
        """Compiles the executable that was created for the specified identifier,
        returns True if the compile was successful."""
        #Find the target folder that has the executables etc then run
        #make and check the exit code.
        target = path.join(self.libraryroot, identifier)
        print("\n\n")
        code = system("cd {}; make F90=ifort".format(target))
        print("\n")
        return code == 0

    def _run_exec(self, identifier, result):
        """Runs the executable for unit test for the specified identifier
        for each of the outcomes specified in the doctags."""
        #Get the home path of the executable. A sub-folder for tests
        #needs to be created. For tests that have input and output files
        #a home/tests/case folder gets created and the source files
        #get copied.

        #Create the folder for staging the tests.
        tests = path.join(self.libraryroot, identifier, "tests")
        if not path.exists(tests):
            mkdir(tests)
        
        #Now we need determine which tests to run from the outcomes and folder tags.
        kmodule, kmethod = identifier.lower().split(".")
        module = self.parser.modules[kmodule]
        method = module.executables[kmethod]

        #Get the absolute path to the executable that we created
        exepath = path.join(self.libraryroot, identifier, "run.x")

        for doc in method.tests:
            if doc.doctype == "outcome"  and "target" in doc.attributes and \
               "sourcepath" in doc.attributes:
                self._run_folder(doc, tests, result, exepath)

        #Now that we have run all of the executables, we can analyze their
        #output to see if it matches.
        for case in result.cases:
            xres = result.cases[case]
            xres.test(case, result)

    def _run_folder(self, doc, testsfolder, result, exepath):
        """Runs the executable for the sources in the folder doctag.

        :arg doc: the docstring 'outcome' that motivates the execution.
        :arg testsfolder: the path to the unit tests unique staging folder.
        :arg result: a TestResult instance that execution results can be
          added to as they complete.
        """
        #The target path from the tag is applicable for all test types.
        source = doc.attributes["sourcepath"].replace(".", self._codefolder)
        #Renames is a linked list of target names for the files in sources.
        if "rename" in doc.attributes:
            renames = doc.attributes["rename"].split(",")
        else:
            renames = []

        system("cd {}".format(source))

        #We can use the same outcome tester afterwards for all of the cases.
        tester = OutcomeTester(doc, self._codefolder, self.comparer, self.parser.verbose)

        #The execution can either be case-based or once-off.
        if "cases" in doc.attributes and "sources" in doc.attributes:
            #We need to run the executable multiple times, once for each case
            #Each case has input files specified relative to the code folder.
            cases = doc.attributes["cases"].split(",")
            for case in cases:
                caseid = "/".join([doc.attributes["sourcepath"].replace(".",""), case])
                if not caseid in result.cases:
                    inputs = doc.attributes["sources"].split(",")
                    #Get the paths to the input files needed to run each case
                    copies = [ path.join(source, i.format(case)) for i in inputs ]

                    #Make a separate directory for the case and copy all its inputs.
                    casepath = path.join(testsfolder, case)
                    if not path.exists(casepath):
                        mkdir(casepath)

                    #Copy all the input files we need to run this case. If a rename
                    #attribute was specified, change the name of the target files.
                    for i in range(len(copies)):
                        f = copies[i]
                        if len(renames) > 0 and i < len(renames):
                            new = path.join(casepath, renames[i])
                            copyfile(f, new)
                        else:
                            copy(f, casepath)

                    print("Executing run.x for case {}".format(caseid))

                    #We need to time the execution as one of the possible outcome tests.
                    start_time = clock()                    
                    code = system("cd {}; {}".format(casepath, exepath))
                    result.cases[caseid] = ExecutionResult(casepath, code, 
                                                           clock() - start_time, tester, case)
                    #Add some whitespace for readability between tests
                    print("")
                else:
                    result.warnings.append("Duplicate CASES specified for unit testing: {}".format(caseid))
        else:
            print("Executing run.x in explicit folder {}".format(target))
            start_time = clock()                    
            code = system("cd {}; {}".format(testsfolder, exepath))
            result.cases[doc.attributes["sourcepath"].replace(".","")] = ExecutionResult(
                testsfolder, code, clock() - start_time, tester)
            #Add some whitespace for readability between tests
            print("")        
