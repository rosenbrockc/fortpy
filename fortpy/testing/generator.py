from .. import msg
from .executable import ExecutableGenerator
from shutil import copy
import xml.etree.ElementTree as ET
import os
import datetime
import dateutil.parser

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
        self.libraryroot = libraryroot
        self.xgenerator = ExecutableGenerator(parser, libraryroot, self)
        if rerun is not None:
            self.rerun = rerun.lower()
        else:
            self.rerun = rerun

        self.dependfiles = [ "fortpy.f90", "Makefile.ifort", "Makefile.gfortran" ]
        self.xtests = {}
        self.xwriters = {}
        self._fortpy = fortpy_templates

        #Stores the identifiers of unit tests whose files changed so they
        #need to be recompiled and executed.
        self._changed = []
        #Load the file dates for previous tests into self.archive
        self._xml_get()

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
                self._write_module(module)
        
    def _write_module(self, module):
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
                self._write_executable(module, anexec)
                
    def get_module_target(self, module):
        """Gets the full path to the file that houses the specified module."""
        if module not in self.parser.mappings:
            return os.path.join(self.xgenerator.folder, module + ".f90")
        else:
            return os.path.join(self.xgenerator.folder, self.parser.mappings[module])        

    def _write_executable(self, module, executable):
        """Generates the fortran program for the specified executable code
        element."""
        #Each test written needs separated from the others
        msg.blank()

        #The way this works is that the latest copy of each module that
        #this executable needs to run is copied to a staging folder
        #where it can be compiled.
        identifier = "{}.{}".format(module.name, executable.name)
        self.xgenerator.reset(identifier, self.libraryroot, self.rerun)
        needs = self.xgenerator.needs()

        #Now we need to check whether the files it depends on have changed
        #since the last time we wrote and compiled.
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
            if needk == "fortpy":
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
                msg.info("COPY {}".format(needed.compile_path))
                copy(needed.compile_path, self.xgenerator.folder)
                different = True
            files[needed.filepath] = moddate

        #We also need to copy across any dependency files that don't exist
        #These don't ever change, so we only need to check for existence
        for dfile in self.dependfiles:
            target = os.path.join(self.xgenerator.folder, dfile)
            dversion = self.tester.get_fortpy_version(target)
            tversion = self.tester.template_version(dfile)
            if not os.path.exists(target) or dversion != tversion:
                source = os.path.join(self._fortpy, dfile)
                msg.info("COPY: {}".format(source))
                copy(source, self.xgenerator.folder)
                different = True

        #We also need to rewrite the files if the user deleted the testid.f90
        #or testid.x files from the directories to force a re-write.
        if not different:
            for testid in self.xgenerator.writer.tests:
                codefile = os.path.join(self.libraryroot, identifier, "{}.f90".format(testid))
                xfile = os.path.join(self.libraryroot, identifier, "{}.x".format(testid))
                if not os.path.exists(codefile) or not os.path.exists(xfile):
                    different = True
                    break

        #All the code files needed for compilation are now in the directory.
        #Create the executable file and the makefile for compilation
        if different:
            msg.okay("UNITTEST: writing executable(s) for {}".format(type(executable).__name__) + 
                     " {}".format(executable.name))
            #It is possible that multiple tests were defined for the executable
            #being unit tested. We need to write a *.f90 PROGRAM file for each
            #test scenario *and* a separate makefile for the executable.
            for testid in self.xgenerator.writer.tests:
                self.xgenerator.write(testid)
                self.xgenerator.makefile(testid)
                self._changed.append("{}|{}".format(identifier, testid))
                self.xtests[identifier] = self.xgenerator.writer.tests
                self.xwriters[identifier] = self.xgenerator.writer
                msg.info("\tWROTE TEST: {}\n".format(testid))

            #Overwrite the file date values for this executable in the archive
            #Also, save the archive in case something goes wrong in the next
            #Executable writing.
            self.archive[identifier] = files
            self._xml_save()
        else:            
            msg.gen("UNITTEST: ignored '{}' because code hasn't changed.".format(executable.name))

    def _xml_get(self):
        """Returns an XML tree for the documont that tracks dates for code
        files and unit tests."""
        target = os.path.join(self.libraryroot, "archive.xml")
        self.archive = {}

        if os.path.exists(target):
            el = ET.parse(target).getroot()

            for test in el:
                files = {}
                for f in test:
                    files[f.attrib["path"]] = dateutil.parser.parse(f.attrib["modified"])
                self.archive[test.attrib["name"]] = files

    def _xml_save(self):
        """Saves the archive dictionary to XML."""
        root = ET.Element("archive")
        for testk in self.archive:
            subel = ET.SubElement(root, "unittest", attrib={ "name": testk })
            for f in self.archive[testk]:
                single = self.archive[testk][f]
                fileel = ET.SubElement(subel, "file", attrib={ "path": f, "modified": single.isoformat() })

        tree = ET.ElementTree(root)
        xmlpath = os.path.expanduser(os.path.join(self.libraryroot, "archive.xml"))
        tree.write(xmlpath)

def modification_date(filename):
    t = os.path.getmtime(filename)
    return datetime.datetime.fromtimestamp(t)
