import fortpy
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

    :arg folder: the folder in which the executable ran.
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
        if self.tester.testspec.runtime is not None and \
           self.tester.testspec.unit is not None:
            #Get the runtime in a nice-to-compare format 0:00:00.000
            stime = secondsToStr(self.runtime)
            #we have limits specified as minutes or hours in a range
            #No unit test should/can take more than 24 hours.
            utime = [ int(u) for u in stime.split(":") ]
            if self.tester.testspec.unit == "h":
                return self._compare_list_time(self.tester.testspec.runtime, utime[0])
            elif self.tester.testspec.unit == "m":
                return self._compare_list_time(self.tester.testspec.runtime, utime[1])  
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
            uresult.overtimes[caseid] = (self.runtime, self.tester.testspec.runtime)

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
    """Performs outcomes tests for a single 'test' DocString.

    :arg testspec: the TestSpecification for the <test> tag.
    :arg codefolder: the path to the folder that has all the modules from which
      the unit tests were built.
    :arg comparer: an instance of FileComparer to compare output files
      to the model ones.
    """
    def __init__(self, testspec, codefolder, comparer, verbose):
        self.testspec = testspec
        self.codefolder = codefolder
        self.comparer = comparer
        self.verbose = verbose
        
    def test(self, caseid, xresult, uresult):
        """Checks the output of the execution for the specified test specification.

        :arg xresult: an ExecutionResult instance that has information about the
          execution whose outcome is being tested.
        :arg uresult: the unittest's TestResult that holds information from all
          executions that were run for different cases.
        """        
        if not self.testspec.testable:
            return

        #Compare the files with their model outputs
        compresults = []
        for i in range(len(self.testspec.targets)):
            target = self.testspec.targets[i]
            exepath, outvar, isfile = self._run_get_paths(target, xresult.folder)

            #The third entry in the tuple identifies it as an explicit value
            #comparison as opposed to a file comparison.
            if isfile:
                result = self._run_compare_file(exepath, outvar, self.codefolder, caseid)
                if result is None or result.common_match < outvar.tolerance:
                    self._add_failure(caseid, result, uresult)
            else:
                result = self._run_compare_var(exepath, outvar)
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

    def _run_compare_var(self, exepath, outvar):
        """Compares the value of a variable written to an output file with an
        explicit value.

        :arg exepath: the full path to the file created by the unit test executable.
        :arg outvar: the TestOutput instance with testing specifications.
        """
        return ValueCompareResult(exe, outvar.value)

    def _run_compare_file(self, exepath, outvar, coderoot, caseid):
        """Compares an output file from an executable with its model output using
        settings in the doctag.

        :arg exe: the full path to the output file from the unit test.
        :arg outvar: the TestOutput instance with testing specifications.
        :arg coderoot: the full path to the folder that has all the code files.
        :arg caseid: the identifier for the specific test case being tested.
        """
        #First get the comparison results, then analyze them using the tolerances
        #etc. from the doctag.
        targetpath = outvar.abspath(coderoot, caseid)
        result = self.comparer.compare(exepath, targetpath, outvar.template, outvar.mode)

        #Write the results out to file as a record. If the result is none create
        #a file that says so. The default file name is the output file name with
        #an extra extension of .compare
        resultpath = exepath + ".compare"
        with open(resultpath, "w") as f:  
            if result is not None:          
                f.write(print_compare_result(result, self.verbose))
            else:
                f.write("The result comparison failed. Check the unit test console output.")
            
        return result

    def _run_get_paths(self, target, exefolder):
        """Gets the file paths to the executable output and model output
        files that need to be compared. If the model output is not a file
        the third entry in the tuple will be False.

        :arg target: the target file or variable that needs to be compared.
        :arg execfolder: the folder in which the executable ran, also where
          the output files from the executable reside.
        :arg caseid: the identifier for the specific test case that ran.
        """
        if target.name[0] == ".":
            exe = path.join(exefolder, target.name[2::])
        else:
            #This is a variable, get the auto-generated filename
            exe = path.join(exefolder, target.varfile)

        #In order to run the comparison, we must have an <output> tag to compare
        #the target to.
        if target.compareto in self.testspec.outputs:
            outvar = self.testspec.outputs[target.compareto]
        else:
            raise ValueError("Target's compareto='{}' ".format(target.compareto) + 
                             "does not match any <output> tag.")

        return (exe, outvar, outvar.filemode)

class TestResult(object):
    """Represents a set of unit test results.

    :arg identifier: the module.method identifier for this unit test.
    :arg testid: the identifier of the <test> tag for the specific test that this
      result represents.
    :attr cases: a dictionary of ExecutionResult objects with detail
      on how the system execution went for each case. Key is caseId
    :attr outcomes: a dictionary of CompareResult objects with detail
      on how similar the new output files are to the model ones. Keys
      are the same caseIds used in self.cases.
    :attr compiled: specifies whether the executable associated with
      this identifier compiled successfully.
    :attr failures: CompareResult or ValueResult instances that failed the
      tests outlined by the relevant test specification. Keys are the
      caseIds used by all the other dicts.
    :attr overtimes: a tuple of (actual runtime, desired runtime range) for
      tests that were outside of the range for runtime.
    """
    def __init__(self, identifier, testid):
        self.identifier = identifier
        self.testid = testid
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
                 fortpy_templates=None, rerun = False, compiler="gfortran"):
        self.libraryroot = path.abspath(libraryroot)
        self.parser = CodeParser()
        self.parser.verbose = verbose
        self._set_fortpy_dir(fortpy_templates)
        self.tgenerator = TestGenerator(self.parser, self.libraryroot, 
                                        self.fortpy_templates, rerun)
        self.compiler = compiler
        self._compiler_exists()
        
        #A flag to track whether the generator has already written
        #the executables.
        self._written = False
        #The user's raw value for the compare_templates directory
        self._compare_templates = compare_templates

    @property
    def tests(self):
        """Returns a dictionary of all the tests that need to be run for the
        current configuration of the test generator."""
        return self.tgenerator.xgenerator.writer.tests

    @property
    def writer(self):
        """Returns the underlying executable writer that has a list of ordered
        method/assignment dependencies."""
        return self.tgenerator.xgenerator.writer
        
    def _compiler_exists(self):
        """Tests whether the specified compiler is available on the machine. If
        it isn't give an error and exit."""
        if self.which(self.compiler) is None:
            print("ERROR: compiler {} not found. Exiting.".format(self.compiler))
            exit(1)
        
    def which(self, program):
        """Tests whether the specified program is anywhere in the environment
        PATH so that it probably exists."""
        import os
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None

    def get_fortpy_templates_dir(self):
        """Gets the templates directory from the fortpy package."""
        fortdir = path.dirname(fortpy.__file__)
        return path.join(fortdir, "templates")

    def _set_fortpy_dir(self, fortpy_templates=None):
        """Sets the directory path for the fortpy templates. If no directory
        is specified, use the default one that shipped with the package.
        """
        #If they didn't specify a custom templates directory, use the default
        #one that shipped with the package.
        if fortpy_templates is not None:
            self.fortpy_templates = path.abspath(fortpy_templates)
        else:
            self.fortpy_templates = self.get_fortpy_templates_dir()

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
        self.tgenerator.write(self._codefolder)
        self._written = True

    def runall(self):
        """Compiles and runs each new or modified unit test executable.
        After the run is complete, the outcomes are checked for consistency."""
        if self._written:
            #We will likely need a file comparer for the outcomes testing
            self.comparer = FileComparer(self.fortpy_templates, self.compare_templates)
            #Run them each individually and return a dictionary of all the
            #test results
            result = {}
            for composite_id in self.tgenerator.tests_to_run:
                identifier, testid = composite_id.split("|")
                result[composite_id] = self._run_single(identifier, testid)

            return result
        else:
            print("WARNING: you can't run tests until the executables have been written. Exiting.")
            return None
       
    def _run_single(self, identifier, testid):
        """Runs all unit test cases for the specified identifier."""
        #Just initialize a result object and populate its properties
        #using the various _run_* methods.
        result = TestResult(identifier, testid)
        result.compiled = self._run_compile(identifier, testid)
        if result.compiled:
            self._run_exec(identifier, testid, result)

        return result

    def _run_compile(self, identifier, testid):
        """Compiles the executable that was created for the specified identifier,
        returns True if the compile was successful."""
        #Find the target folder that has the executables etc then run
        #make and check the exit code.
        target = path.join(self.libraryroot, identifier)
        print("\n\n")
        code = system("cd {}; make -f 'Makefile.{}' F90={}".format(target, testid, self.compiler))
        print("\n")
        return code == 0

    def _run_exec(self, identifier, testid, result):
        """Runs the executable for unit test for the specified identifier
        for each of the outcomes specified in the doctags."""
        #Get the home path of the executable. A sub-folder for tests
        #needs to be created. For tests that have input and output files
        #a home/tests/testid.case folder gets created and the source files
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
        exepath = path.join(self.libraryroot, identifier, "{}.x".format(testid))
        #Since we have already split out all the tests that need to be run and 
        #we have a 'testid' for the current test to run, just run that test.
        self._run_folder(self.tests[testid], tests, result, exepath)

        #Now that we have run all of the executables, we can analyze their
        #output to see if it matches.
        for case in result.cases:
            xres = result.cases[case]
            xres.test(case, result)

    def _run_folder(self, testspec, testsfolder, result, exepath):
        """Runs the executable for the sources in the folder doctag.

        :arg testspec: a TestSpecification instance for this unit test.
        :arg testsfolder: the path to the unit tests unique staging folder.
        :arg result: a TestResult instance that execution results can be
          added to as they complete.
        :arg expath: the full path to the executable to run in the folder.
        """
        #We can use the same outcome tester afterwards for all of the cases.
        tester = OutcomeTester(testspec, self._codefolder, self.comparer, self.parser.verbose)

        #The execution can either be case-based or once-off.
        if testspec.cases is not None and len(testspec.inputs) > 0:
            #We need to run the executable multiple times, once for each case
            #Each case has input files specified relative to the code folder.
            for case in testspec.cases:
                caseid = "{}.{}".format(testspec.identifier, case)
                if not caseid in result.cases:
                    #Make a separate directory for the case and copy all its inputs.
                    casepath = path.join(testsfolder, caseid)
                    if not path.exists(casepath):
                        mkdir(casepath)

                    for i in testspec.inputs:
                        i.copy(self._codefolder, casepath, case)
                    #Also copy any assignment file dependencies.
                    self.writer.copy(self._codefolder, casepath, case)
                    #Clean the testing folder to remove any target variable output files
                    #from any earlier test runs.
                    testspec.clean(casepath)

                    print("Executing {}.x for case {}".format(testspec.identifier, case))

                    #We need to time the execution as one of the possible outcome tests.
                    start_time = clock()                    
                    code = system("cd {}; {}".format(casepath, exepath))
                    result.cases[caseid] = ExecutionResult(casepath, code, clock() - start_time, 
                                                           tester, case)

                    #Add some whitespace for readability between tests
                    print("")
                else:
                    result.warnings.append("Duplicate CASES specified for unit testing: {}".format(caseid))
        else:
            #Create a folder for this test specification to run in.
            testpath = path.join(testsfolder, testspec.identifier)
            if not path.exists(testpath):
                mkdir(testpath)
            print("Executing {}.x in explicit folder {}".format(testspec.identifier, testpath))

            #Copy across all the input files we need to run.
            for i in testspec.inputs:
                i.copy(self._codefolder, testpath)
            #Also copy any assignment file dependencies.
            self.writer.copy(self._codefolder, testpath, "")
            #Clean the testing folder to remove any target variable output files
            #from any earlier test runs.
            testspec.clean(testpath)


            start_time = clock()                              
            code = system("cd {}; {}".format(testpath, exepath))
            result.cases[testspec.identifier] = ExecutionResult(
                testpath, code, clock() - start_time, tester)

            #Add some whitespace for readability between tests
            print("")        
