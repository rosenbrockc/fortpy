from .. import msg
from .executable import ExecutableGenerator
from shutil import copy
import xml.etree.ElementTree as ET
import os
import datetime
import dateutil.parser
from fortpy.code import config

class TestGenerator(object):
    """Generates automatic unit tests based on docstrings in the fortran
    doc elements defined in the code files.
    :arg parser: an instance of the code parser that has code elements
      with all the docstrings.
    :arg tester: an instance of UnitTester handling the overall testing
      automation.
    :arg libraryroot: the path to the folder that will contain all the
      unit test code files and executables.
    :arg fortpy_templates: the path to the fortpy templates folder to
      copy dependencies from.
    :arg rerun: specifies whether to re-run the tests (i.e. don't recopy
      and re-compile everything, just re-run the tests using existing
      code files. May fail if something has changed.
    :attr xgenerator: an instance of ExecutableGenerator to generate the
      .f90 program files for performing the unit tests.
    :attr dependfiles: a list of additional files to be copied from the
      fortpy templates directory that are needed for the unit testing.
    :attr xtests: a dictionary of the test collections for each of the
      executables visited by this generator.
    :attr xwriters: a dictionary of the writers from each executable that
      was visited by the generator.
    """
    def __init__(self, parser, libraryroot, fortpy_templates, tester, rerun=None):
        self.parser = parser
        self.tester = tester
        self.xgenerator = ExecutableGenerator(parser, libraryroot, self)
        if rerun is not None:
            self.rerun = rerun.lower()
        else:
            self.rerun = rerun

        self.dependfiles = []
        self.xtests = {}
        """Dictionary of the dictionary of test specifications indexed by test identifier.
        """
        self.xwriters = {}
        """Dictionary of the MethodWriter instances for each test run, indexed by
        test identifier.
        """
        self.archive = None
        """A dictionary of the unit tests' code files and when they were last compiled
        for unit testing. Used to determine if the underlying code has changed since the
        test was last run.
        """
        self._fortpy = fortpy_templates

        #Stores the identifiers of unit tests whose files changed so they
        #need to be recompiled and executed.
        self._changed = []

    @property
    def tests_to_run(self):
        """Returns the identifiers for those unit tests that are new or
        needed to be recreated because of code changes."""
        return self._changed

    def write(self, codefolder):
        """Creates a fortran program for each subroutine in the code parsers
        modules lists that tests the subroutine/function.
        :arg codefolder: the full path to the folder in which the code files
          reside that tests will be run for.
        """
        #We need to enumerate over a *copy* of the keys list since the list of 
        #modules is likely to change during the execution as dependencies
        #are found and loaded.
        currentlist = list(self.parser.modules.keys())
        lcoderoot = codefolder.lower()

        for mkey in currentlist:
            #We only want to perform unit tests for executables that are in the
            #same code folder as the one being executed by the unit tester.
            module = self.parser.modules[mkey]
            if lcoderoot in module.filepath.lower():
                self._write_module(module, codefolder)
        
    def _write_module(self, module, coderoot):
        """Generates the fortran programs for all executables in the module
        code element specified."""      
        for execkey in module.executables:
            anexec = module.executables[execkey]
            #We need to check whether this executable has any outcomes defined
            #for a unit test. If it doesn't, we just skip it.
            found = False
            i = 0

            #The executable must have a testing group *and* be marked as public
            #in its parent module (either via a public modifier or via the public
            #keyword in the module).
            if anexec.test_group is not None:
                self._write_executable(module, anexec, coderoot)
                
    def get_module_target(self, module):
        """Gets the full path to the file that houses the specified module."""
        if module not in self.parser.mappings:
            return os.path.join(self.xgenerator.folder, module + ".f90")
        else:
            return os.path.join(self.xgenerator.folder, self.parser.mappings[module])        

    def _write_executable(self, module, executable, coderoot):
        """Generates the fortran program for the specified executable code
        element."""
        #Each test written needs separated from the others
        msg.std("{}".format(type(executable).__name__.lower()) + 
                " {}".format(executable.name))

        #The way this works is that the latest copy of each module that
        #this executable needs to run is copied to a staging folder
        #where it can be compiled.
        identifier = "{}.{}".format(module.name, executable.name)
        from os import path
        newstage = self.xgenerator.reset(identifier, path.dirname(module.filepath))
        self._xml_get(identifier, newstage)
        needs = self.xgenerator.needs()

        #Now we need to check whether the files it depends on have changed
        #since the last time we wrote and compiled. If re-run is specified,
        #force the difference check to be True so that it re-writes the driver.
        if self.rerun is not None and (self.rerun == "*" or self.rerun in identifier.lower()):
            different = True
        else:
            different = False
            
        if identifier in self.archive:
            previous = self.archive[identifier]
        else:
            previous = {}

        #As we copy files, we need to keep track of the last date they were
        #modified so that we only retest things that change.
        files = {}
        for needk in needs:
            #Ignore references to fortpy since we include it automatically.
            #Also ignore MPI since it isn't something we need to worry about.
            if needk == "fortpy" or needk=="fpy_auxiliary":
                continue

            needed = self.parser.modules[needk]
            moddate = modification_date(needed.filepath)
            #We also consider a module altered if its XML file has changed since
            #we last parsed it.
            xmlpath = needed.xmlpath
            if os.path.isfile(xmlpath):
                xmoddate = modification_date(xmlpath)
                if xmoddate > moddate:
                    moddate = xmoddate

            #Get the path to the code file in the executable directory so that
            #we can copy it over if it doesn't exist.
            target = self.get_module_target(needk)

            if needed.filepath not in previous or \
               (needed.filepath in previous and previous[needed.filepath] < moddate) or \
               not os.path.exists(target):
                msg.info("   COPY {}".format(needed.compile_path))
                copy(needed.compile_path, self.xgenerator.folder)
                different = True
            files[needed.filepath] = moddate

        #We also need to copy across any dependency files that don't exist
        #These don't ever change, so we only need to check for existence
        for dfile in self.dependfiles:
            target = os.path.join(self.xgenerator.folder, dfile)
            dversion = self.tester.get_fortpy_version(target)
            tversion = self.tester.template_version(dfile)
            if config.symlink and (os.path.isfile(target) or os.path.islink(target)):
                os.remove(target)
                
            if not os.path.exists(target) or dversion != tversion:
                source = os.path.join(self._fortpy, dfile)
                if config.symlink:
                    os.symlink(source, target)
                else:
                    copy(source, self.xgenerator.folder)
                different = True

        #We also need to rewrite the files if the user deleted the testid.f90
        #file from the directories to force a re-write.
        if not different:
            for testid in self.xgenerator.writer.tests:
                codefile = os.path.join(self.xgenerator.folder, "{}.f90".format(testid))
                if not os.path.exists(codefile):
                    msg.std("Missing codefile in {}; re-run.".format(self.xgenerator.folder))
                    different = True
                    break

        #All the code files needed for compilation are now in the directory.
        #Create the executable file and the makefile for compilation
        if different:
            if len(self.xgenerator.writer.tests) > 1:
                msg.okay("   CASES:")
            #It is possible that multiple tests were defined for the executable
            #being unit tested. We need to write a *.f90 PROGRAM file for each
            #test scenario *and* a separate makefile for the executable.
            for testid in self.xgenerator.writer.tests:
                self.xgenerator.write(testid, coderoot)
                self.xgenerator.makefile(testid)
                self._changed.append("{}|{}".format(identifier, testid))
                self.xtests[identifier] = self.xgenerator.writer.tests
                self.xwriters[identifier] = self.xgenerator.writer
                spacer = "  " if len(self.xgenerator.writer.tests) > 1 else ""
                msg.info("   {}Wrote Test: {}".format(spacer, testid))

            #Overwrite the file date values for this executable in the archive
            #Also, save the archive in case something goes wrong in the next
            #Executable writing.
            self.archive[identifier] = files
            self._xml_save(identifier)
        else:            
            msg.gen("   IGNORED: code hasn't changed.".format(executable.name))

    def _xml_get(self, identifier, reset=False):
        """Returns an XML tree for the document that tracks dates for code
        files and unit tests."""
        if self.archive is not None and not reset:
            #We only load this once per reset of the test generator's executable generator.
            return
        
        target = os.path.join(self.xgenerator.folders[identifier], "archive.xml")
        self.archive = {}

        if os.path.exists(target):
            el = ET.parse(target).getroot()

            for test in el:
                files = {}
                for f in test:
                    files[f.attrib["path"]] = dateutil.parser.parse(f.attrib["modified"])
                self.archive[test.attrib["name"]] = files

    def _xml_save(self, identifier):
        """Saves the archive dictionary to XML."""
        root = ET.Element("archive")
        for testk in self.archive:
            subel = ET.SubElement(root, "unittest", attrib={ "name": testk })
            for f in self.archive[testk]:
                single = self.archive[testk][f]
                fileel = ET.SubElement(subel, "file", attrib={ "path": f, "modified": single.isoformat() })

        tree = ET.ElementTree(root)
        xmlpath = os.path.expanduser(os.path.join(self.xgenerator.folders[identifier], "archive.xml"))
        tree.write(xmlpath)

def modification_date(filename):
    t = os.path.getmtime(filename)
    return datetime.datetime.fromtimestamp(t)
