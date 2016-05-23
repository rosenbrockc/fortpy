from .. import msg
import fortpy
from .generator import TestGenerator
from ..code import CodeParser, secondsToStr
from os import system, path, mkdir, remove
from time import clock
import time
from shutil import copy, copyfile
from fortpy.code import secondsToStr
from .comparer import FileComparer
from fortpy.testing.results import print_compare_result
from fortpy.testing import profiling

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
        """Tests the outcome of this executable using its tester attribute.

        :arg uresult: the overall TestResult for the test.
        """
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
        import re
        with open(self.path) as f:
            for line in f:
                if not re.match("^\s*#", line):
                    break

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
    def __init__(self, testspec, codefolder, comparer, verbose, compiler):
        self.testspec = testspec
        self.codefolder = codefolder
        self.comparer = comparer
        self.verbose = verbose
        self.compiler = compiler
        """Active compiler for the current test whose outcome is being run."""
        
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
            exepath, outvar, isfile, exists = self._run_get_paths(target, xresult.folder, caseid)
            #Exists tells us if we have the necessary output to compare with. If we don't
            #we can't complete the comparison.
            if not exists:
                self._add_failure(caseid, None, uresult)
                continue

            #The third entry in the tuple identifies it as an explicit value
            #comparison as opposed to a file comparison.
            if isfile and not outvar.autoclass:
                result = self._run_compare_file(exepath, outvar, self.codefolder, caseid)
                if result is None or result.common_match < outvar.tolerance:
                    self._add_failure(caseid, result, uresult)
            elif outvar.autoclass:
                result = self._run_compare_autoclass(exepath, outvar, self.codefolder, caseid)
                if result.common_match < outvar.tolerance:
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
        return ValueCompareResult(exepath, outvar.value)

    def _run_compare_autoclass(self, exepath, outvar, coderoot, caseid):
        """Compares an output *folder* from an executable (for a variable saved using
        the autoclass feature) with its model output using settings in the doctag.

        :arg exe: the full path to the output folder from the unit test.
        :arg outvar: the TestOutput instance with testing specifications.
        :arg coderoot: the full path to the folder that has all the code files.
        :arg caseid: the identifier for the specific test case being tested.
        """
        #We need to get a list of all the files in the model and test folders
        #for the variable and then compare each of them with _run_compare_file.
        #Luckily, there are no nested folders, the files are saved linearly and
        #the filenames contain the recursive complexity of the variable.
        from os import walk, path, mkdir
        modelpath = outvar.abspath(coderoot, caseid, self.compiler)
        mfiles = []
        xfiles = []
        for (dirpath, dirnames, filenames) in walk(modelpath):
            mfiles.extend(filenames)
            break

        for (dirpath, dirnames, filenames) in walk(exepath):
            xfiles.extend(filenames)
            break

        #Create a directory for all the .compare files of the member variables
        #in the auto-classed variable.
        compath = exepath + ".compare"
        if not path.isdir(compath):
            mkdir(compath)

        #Hopefully we get a one-to-one match, otherwise we have to record the
        #set difference as failures.
        onlym = []
        mx = []
        summary = {}
        for m in mfiles:
            if m[0] != "_" or ".fpy" in m:
                #Ignore the files that don't follow the convention; otherwise the
                #statistics will be messed up.
                continue
            if m in xfiles:
                xpath = path.join(exepath, m)
                mpath = path.join(modelpath, m)
                mxres = self.comparer.compare(xpath, mpath, outvar.template, outvar.mode)
                mx.append(mxres)
                summary[m] = (mxres.percent_match, mxres.common_match)

                #Write a comparison report for this particular variable.
                cpath = path.join(compath, m)
                with open(cpath, "w") as f:
                    if mxres is not None:
                        f.write(print_compare_result(mxres, self.verbose))
                    else:
                        f.write("The result comparison failed. Check the unit test console output.")
            else:
                onlym.append(m)
        onlyx = [f for f in xfiles if (f[0] == "_" and ".fpy" not in f and f not in mfiles)]

        #It turns out to be useful to have a summary file that lists the percentage
        #for each file comparison. Otherwise, it is hard to track down where the errors
        #have occurred.
        from fortpy.testing.results import ACResult
        reltot = len(summary) + len(onlyx) + len(onlym)
        result = ACResult(mx, onlym, onlyx, outvar.actolerance)
        slines = ["{0:.2%} success ({1:.2%} common) TOTAL MATCH".format(result.percent_match, result.common_match),
                  "",
                  "OUTPUT : {}".format(exepath),
                  "MODEL  : {}".format(modelpath),
                  "",
                  "Files in both MODEL and TEST directories ({}/{}):".format(len(summary), reltot)]
        for mfile, match in summary.items():
            slines.append("{0} => {1:.2%} success ({2:.2%} common)".format(mfile, match[0], match[1]))

        if len(onlym) > 0:
            slines.append("")
            slines.append("Files only present in the MODEL output ({}/{}):".format(len(onlym), reltot))
            slines.extend(onlym)
        if len(onlyx) > 0:
            slines.append("")
            slines.append("Files only present in the TEST output ({}/{}):".format(len(onlyx), reltot))
            slines.extend(onlyx)

        sfile = path.join(compath, "summary")
        with open(sfile, 'w') as f:
            f.write('\n'.join(slines))

        return result

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
        targetpath = outvar.abspath(coderoot, caseid, self.compiler)
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

    def _run_get_paths(self, target, exefolder, caseid):
        """Gets the file paths to the executable output and model output
        files that need to be compared. If the model output is not a file
        the third entry in the tuple will be False.

        :arg target: the target file or variable that needs to be compared.
        :arg execfolder: the folder in which the executable ran, also where
          the output files from the executable reside.
        :arg caseid: the identifier for the specific test case that ran.
        """
        #In order to run the comparison, we must have an <output> tag to compare
        #the target to.
        if target.compareto in self.testspec.outputs:
            outvar = self.testspec.outputs[target.compareto]
        else:
            raise ValueError("Target's compareto='{}' ".format(target.compareto) + 
                             "does not match any <output> tag.")

        if target.name[0] == ".":
            exe = path.join(exefolder, target.name[2::])
        elif outvar.autoclass:
            exe = path.join(exefolder, target.varfile)
        else:
            #This is a variable, get the auto-generated filename
            exe = path.join(exefolder, target.varfile)

        if outvar.filemode:
            abspath = outvar.abspath(self.codefolder, caseid, self.compiler)
            exists = ((not outvar.autoclass and path.isfile(abspath)) or
                      (outvar.autoclass and path.isdir(abspath)))
            if not exists:
                if outvar.autoclass:
                    msg.warn("Auto-class model output folder {} ".format(abspath) +
                             "does not exist for case {}".format(caseid))
                else:
                    msg.warn("Model output file {} does not exist for case {}".format(abspath, caseid))
            return (exe, outvar, True, exists)
        else:
            return (exe, outvar, False, True)

class TestResult(object):
    """Represents a set of unit test results.

    :arg identifier: the module.method identifier for this unit test.
    :arg testid: the identifier of the <test> tag for the specific test that this
      result represents.
    :attr cases: a dictionary of ExecutionResult objects with detail
      on how the system execution went for each case. Key is caseId
    :attr paths: a dictionary of paths to the folders that were executed for each case.
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
        self.paths = {}
        self.outcomes = {}
        self.compiled = None
        self.target = None
        """The full path to the directory containing the compiled executable
        to run."""

        self.warnings = []
        self.failures = {}
        self.overtimes = {}

    def clean(self, case):
        """Dereferences the test results' memory to the representations of
        the output and only retains the percent and common matches etc. for
        aggregate statistics. Should only be called after the *.compare files
        have been written out for the test case.
        """
        if case in self.outcomes:
            if isinstance(self.outcomes[case], list):
                for result in self.outcomes[case]:
                    if result is not None:
                        result.clean()
            else:
                self.outcomes[case].clean()

        #Technically, the self.failures dict also has pointers to test result
        #instances. But these will also be in the outcomes dictionary, so they
        #will have been cleaned already.
        
    @property
    def totaltime(self):
        """Returns the total amount of time spent running *only* the unit tested
        executable for *all* cases in this result."""
        time = 0
        for caseid in self.paths:
            casetime = self.runtime(caseid)
            if casetime is not None:
                time += casetime

        if time == 0:
            return None
        else:
            return time

    def runtime(self, caseid):
        """Gets the absolute amount of time spent running *only* the executable
        being unit tested for this result."""
        time = None

        if caseid in self.paths:
            timepath = path.join(self.paths[caseid], "fpytiming.out")
            if path.isfile(timepath):
                with open(timepath) as f:
                    time = float(f.readlines()[1].strip())

        return time

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
                        total += 1/len(case)
                else:
                    if result is not None:
                        total += result.percent_match/len(case)

        if len(self.cases) == 0:
            return 0
        else:
            return float(total) / (len(self.cases) + self.failure_count)
            
    @property
    def common(self):
        """Returns the common match percent success."""
        total = 0
        common = 0

        for caseid in self.outcomes:
            case = self.outcomes[caseid]
            for result in case:
                if isinstance(result, ValueCompareResult):
                    if result.equal:
                        total += 1/len(case)
                else:
                    if result is not None:
                        total += result.common_match/len(case)
        if len(self.cases) == 0:
            return 0
        else:
            return float(total) / (len(self.cases) + self.failure_count)

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
    def __init__(self, libraryroot=None, verbose=False, compare_templates=None,
                 fortpy_templates=None, rerun=None , compiler=None,
                 debug=False, profile=False, quiet=False, strict=False):
        self.parser = CodeParser()
        self.parser.verbose = verbose
        from fortpy.utility import set_fortpy_templates
        set_fortpy_templates(self, fortpy_templates)
        self.tgenerator = TestGenerator(self.parser, libraryroot, self.fortpy_templates, self, rerun)
        self.compiler = None
        """The full path to the compiler to use for the unit tests.
        """
        self.quiet = quiet
        """Specifies whether the tester should run in quiet mode, which only prints
        essential information for the unit tests to stdout.
        """
        self.strict = strict
        """Specifies whether to compile files with *all* warnings enabled."""
        self.debug = debug == True
        self.set_compiler(compiler)
        self.profile = self._profiler_exists(profile)
        
        #A flag to track whether the generator has already written
        #the executables.
        self._written = False
        #The user's raw value for the compare_templates directory
        self._compare_templates = compare_templates

    def get_fortpy_version(self, fortpath):
        """Gets the fortpy version number from the first line of the specified file."""
        from fortpy.testing.compilers import get_fortpy_version
        return get_fortpy_version(self.compiler, fortpath) 

    def tests(self, identifier):
        """Returns a dictionary of all the tests that need to be run for the
        specified module.executable identifier."""
        if identifier in self.tgenerator.xtests:
            return self.tgenerator.xtests[identifier]
        else:
            return {}

    def writer(self, identifier):
        """Returns the underlying executable writer that has a list of ordered
        method/assignment dependencies."""
        if identifier in self.tgenerator.xwriters:
            return self.tgenerator.xwriters[identifier]
        
    def libraryroot(self, identifier):
        """Returns the absolute path to the staging directory for the unit tests with
        the specified testid.
        """
        if identifier in self.tgenerator.xgenerator.folders:
            return self.tgenerator.xgenerator.folders[identifier]

    def _profiler_exists(self, profile):
        """Tests whether we have gprof available to do the profiling of the methods
        to unit test.
        
        :arg profile: the value specified in the UnitTester constructor.
        """
        if profile==True:
            from fortpy.utility import which
            gprof = which("gprof")
            if gprof is None:
                msg.err("gprof is required to run profiling with fortpy.")
                exit(1)
            else:
                return True
        else:
            return False

    def set_compiler(self, compiler):
        """Sets the compiler to use for the unit testing of this code parser.
        """
        if compiler is not None:
            self.compiler = compiler
            self._compiler_exists()

    def _compiler_exists(self):
        """Tests whether the specified compiler is available on the machine. If
        it isn't give an error and exit."""
        from fortpy.testing.compilers import compilers
        if self.compiler in compilers:
            #Overwrite the *name* of the compiler with its full path; since
            #fortpy assumes that self.compiler is the name of a valid executable
            #this will still work correctly.
            compiler = compilers[self.compiler].path
        else:
            compiler = self.compiler

        from fortpy.utility import which
        if which(compiler) is None:
            msg.err("compiler {} not found. Exiting.".format(self.compiler))
            exit(1)
        
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

    def runall(self, compiler=None):
        """Compiles and runs each new or modified unit test executable.
        After the run is complete, the outcomes are checked for consistency.

        :arg compiler: the name of the compiler in 'compilers.xml' to use.
        """
        #We will likely need a file comparer for the outcomes testing
        self.comparer = FileComparer(self.fortpy_templates, self.compare_templates)

        if self._written:
            self.set_compiler(compiler)

            from fortpy.testing.compilers import replace
            from fortpy.testing.auxiliary import generate
            from fortpy.utility import copy

            #Run them each individually and return a dictionary of all the
            #test results
            result = {}
            fpyauxs= []            
            for composite_id in self.tgenerator.tests_to_run:
                identifier, testid = composite_id.split("|")
                #Compile and copy fpy_auxiliary if it isn't in the identifiers directory yet.
                source = path.join(self.libraryroot(identifier), identifier)
                target = replace(source + ".[c]", self.compiler)
                if self.writer(identifier).autoclass and identifier not in fpyauxs:
                    code, success, fpytarget = generate(self.parser, self._codefolder,
                                                        self.libraryroot(identifier), self.compiler,
                                                        self.debug, self.profile, self.strict)
                    opath = path.join(fpytarget, "fpy_auxiliary.o")
                    mpath = path.join(fpytarget, "fpy_auxiliary.mod")
                    copy(opath, target)
                    copy(mpath, target)
                    fpyauxs.append(identifier)
                
                oneresult = self._run_single(identifier, testid, source)
                if oneresult is not None:
                    result[composite_id] = oneresult

            return result
        else:
            msg.warn("you can't run tests until the executables have been written. Exiting.")
            return None
       
    def _run_single(self, identifier, testid, staging):
        """Runs all unit test cases for the specified identifier."""
        #Just initialize a result object and populate its properties
        #using the various _run_* methods.
        result = TestResult(identifier, testid)
        result.compiled, result.target = self._run_compile(identifier, testid, staging)
        if result.compiled:
            self._run_exec(identifier, testid, result)

        if self.tests(identifier)[testid].runchecks:
            #Only return a test result if the checks were actually run.
            return result
        
    def _run_compile(self, identifier, testid, staging):
        """Compiles the executable that was created for the specified identifier,
        returns True if the compile was successful."""
        #We need to make sure that the fpy_auxiliary .mod and .o files are available
        #if the test uses auto-class functionality.
        from os import path
        from fortpy.testing.compilers import compile_general          
        msg.blank()
        code, success, target = compile_general(staging, self.compiler, testid, self.debug,
                                                self.profile, self.quiet, strict=self.strict)

        if not success:
            #There are 21 lines in the compile.log file when everything runs correctly
            #Overwrite code with a bad exit value since we have some other problems.
            code = 1

            #If the executable exists, we could still prompt them to run it (in case the
            #additional lines were just warnings).
            exe = path.join(target, "{}.x".format(testid))
            if path.isfile(exe):
                choice = raw_input("\nWould you still like to run the executable? ").lower()
                code = 0 if "y" in choice else code
                if "n" in choice:
                    msg.err("Unit testing terminated by user.")
                    exit(0)
            else:
                msg.err("Could not compile executable {}.x".format(testid))
                exit(-1)

        return code == 0, target

    def _run_exec(self, identifier, testid, result):
        """Runs the executable for unit test for the specified identifier
        for each of the outcomes specified in the doctags."""
        if not self.tests(identifier)[testid].execute:
            #We don't want to carry on with this execution at all. User-specified
            #override.
            return

        #Get the home path of the executable. A sub-folder for tests
        #needs to be created. For tests that have input and output files
        #a home/tests/testid.case folder gets created and the source files
        #get copied.

        #Create the folder for staging the tests.
        tests = path.join(result.target, "tests")
        if not path.exists(tests):
            mkdir(tests)
        
        #Now we need determine which tests to run from the outcomes and folder tags.
        kmodule, kmethod = identifier.lower().split(".")
        module = self.parser.modules[kmodule]
        method = module.executables[kmethod]

        #Get the absolute path to the executable that we created
        exepath = path.join(result.target, "{}.x".format(testid))

        #Since we have already split out all the tests that need to be run and 
        #we have a 'testid' for the current test to run, just run that test.
        self._run_folder(self.tests(identifier)[testid], tests, result, exepath,
                         self.writer(identifier))

        if not self.tests(identifier)[testid].runchecks:
            return
        
        #Now that we have run all of the executables, we can analyze their
        #output to see if it matches.
        from tqdm import tqdm
        pbar = tqdm(result.cases)
        for case in pbar:
            if "." in case:
                pbar.set_description("Checking {}".format(case.split(".")[1]))
            else:
                #The default case doesn't have any extra id.
                pbar.set_description("Checking {}".format(case))
                
            xres = result.cases[case]
            #This next step generates large representations in memory of all the output
            #files that need to be compared. If we don't clean it up here, then a test
            #with thousands of cases could use lots of memory. Or if you have 100 executables
            #each with 1000 test cases, then you have 100,000 representations of model and
            #actual output files, with their differences, etc.
            xres.test(case, result)

            #Now, we should dereference all the test results for this case since we have
            #the percent matches and the actual comparisions should be written out to
            #*.compare files in the individual test case directories.
            result.clean(case)

    def _run_folder(self, testspec, testsfolder, result, exepath, testwriter):
        """Runs the executable for the sources in the folder doctag.

        :arg testspec: a TestSpecification instance for this unit test.
        :arg testsfolder: the path to the unit tests unique staging folder.
        :arg result: a TestResult instance that execution results can be
          added to as they complete.
        :arg expath: the full path to the executable to run in the folder.
        :arg testwriter: the MethodWriter instance that generated the test
          that will be run.
        """
        #We can use the same outcome tester afterwards for all of the cases.
        tester = OutcomeTester(testspec, self._codefolder, self.comparer, 
                               self.parser.verbose, self.compiler)

        #The execution can either be case-based or once-off.
        msg.okay("Executing {}.x in {}".format(testspec.identifier, testsfolder))
        if testspec.cases is not None:
            #We need to run the executable multiple times, once for each case
            #Each case has input files specified relative to the code folder.
            from tqdm import tqdm
            pbar = tqdm(testspec.cases) if not self.quiet else testspec.cases
            for case in pbar:
                caseid = "{}.{}".format(testspec.identifier, case)
                if not self.quiet:
                    pbar.set_description("Running {}".format(case))
                if not caseid in result.cases:
                    #Make a separate directory for the case and copy all its inputs.
                    casepath = path.join(testsfolder, caseid)
                    self._execute_testpath(testspec, testwriter, casepath, exepath, 
                                           result, tester, caseid, case)
                else:
                    result.warnings.append("Duplicate CASES specified for unit testing:" + 
                                           " {}".format(caseid))
        else:
            #Create a folder for this test specification to run in.
            testpath = path.join(testsfolder, testspec.identifier)
            self._execute_testpath(testspec, testwriter, testpath, exepath, result, 
                                   tester, testspec.identifier)

    def _execute_testpath(self, testspec, testwriter, testpath, exepath, result, 
                          tester, caseid, case=""):
        """Executes the unit test in the specified testing folder for 'case'."""
        if not path.exists(testpath):
            mkdir(testpath)

        #Copy across all the input files we need to run.
        for i in testspec.inputs:
            i.copy(self._codefolder, testpath, case, self.compiler,
                   self.tgenerator.xgenerator.libraryroot)
        #Also copy any assignment file dependencies.
        testwriter.copy(self._codefolder, testpath, case, testspec.identifier, self.compiler)
        #Clean the testing folder to remove any target variable output files
        #from any earlier test runs.
        testspec.clean(testpath)
        testwriter.setup(testspec.identifier, testpath)

        #If the testspec needs auto-class support, write the case to file.
        if testspec.autoclass:
            with open(path.join(testpath, "fpy_case"), 'w') as f:
                f.write('"{}"'.format(case))
        
        #Save the path to the folder for execution in the result.
        result.paths[caseid] = testpath
        start_time = clock()                              
        from os import waitpid
        from subprocess import Popen, PIPE
        command = "cd {}; {} > .fpy.x.out".format(testpath, exepath)
        prun = Popen(command, shell=True, executable="/bin/bash", stdout=PIPE, stderr=PIPE)
        waitpid(prun.pid, 0)        
        if not self.quiet:
            output = prun.stdout.readlines()
            if len(output) > 0:
                msg.std(''.join(output))
        #else: #We don't need to get these lines since we are purposefully redirecting them.
        error = prun.stderr.readlines()
        if len(error) > 0:
            if self.quiet:
                msg.info("With Executable at {}".format(exepath), 1)
            msg.err('\n  '+'  '.join(error))
        code = len(error)
        
        if case == "":
            result.cases[caseid] = ExecutionResult(testpath, code, clock() - start_time, tester)
        else:
            result.cases[caseid] = ExecutionResult(testpath, code, clock() - start_time, tester, case)
        self._write_success(testpath, code)

        #See if we need to run the post-execution profiling
        if self.profile:
            profiling.profile(testpath, testspec.testgroup.method_fullname, 
                              exepath, self.compiler)

    def _write_success(self, testpath, code):
        """Creates a SUCCESS file in the specified testpath if code==0 that has
        the time of the last execution. If code != 0, any existing SUCCESS file
        is deleted.
        """
        sucpath = path.join(testpath, "SUCCESS")
        if code == 0:
            with open(sucpath, 'w') as f:
                f.write(time.strftime("%c"))
        else:
            if path.isfile(sucpath):
                remove(sucpath)
