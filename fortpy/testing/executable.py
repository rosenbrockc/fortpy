from .. import msg
from fortpy.testing.method import MethodWriter
from os import path, mkdir, remove
from datetime import datetime

class ExecutableGenerator(object):
    """Generates a fortran executable to perform unit tests for a given
    subroutine or function.

    :arg parser: an instance of the code parser for inter-module access
    :arg folder: the folder in which to generate the code files and execs.
    :arg testgen: the TestGenerator instance that owns this executable generator.
    """
    def __init__(self, parser, folder, testgen):
        #We need the code parser because the pre-reqs specified may be in other
        #modules and may have pre-reqs of their own. This way, we can find
        #them all easily.
        self.parser = parser
        self.folder = folder
        self.testgen = testgen

        self._needs = None

    def needs(self):
        """Returns a list of all the modules that this executable needs to
        run correctly."""
        if self._needs is None:
            self._needs = self._calc_needs()
        return self._needs

    def _calc_needs(self):
        """Calculates a list of all the modules that this executable needs to
        run correctly."""
        result = []        
        for modk in self.writer.uses():
            if modk not in result:
                result.append(modk)

        #We also need to look up the dependencies of each of these modules
        recursed = list(result)
        for i in range(len(result)):
            module = result[i]
            self._process_module_needs(module, i, recursed)

        return recursed

    def _process_module_needs(self, module, i, result):
        """Adds the module and its dependencies to the result list."""
        #Some code might decide to use the fortpy module methods for general
        #development, ignore it since we know it will be present in the end.
        if module == "fortpy":
            return

        #See if the parser has alread loaded this module.
        if module not in self.parser.modules:
            self.parser.load_dependency(module, True, True, False)

        #It is possible that the parser couldn't find it, if so
        #we can't create the executable!
        if module in self.parser.modules:
            modneeds = self.parser.modules[module].needs
            for modn in modneeds:
                if modn not in result:
                    #Since this module depends on the other, insert the other
                    #above it in the list.
                    result.insert(i, modn)
                else:
                    x = result.index(modn)
                    if x > i:
                        #We need to move this module higher up in the food chain
                        #because it is needed sooner.
                        result.remove(modn)
                        result.insert(i, modn)

                newi = result.index(modn)
                self._process_module_needs(modn, newi, result)
        else:
            raise ValueError("unable to find module {}.".format(module))

    def reset(self, identifier, libraryroot, rerun=None):
        """Resets the writer to work with a new executable."""        
        #A method writer has all the guts needed to generate the executable
        #We are just wrapping it and making sure the variable values get
        #recorded to temporary files for comparison etc.
        self.identifier = identifier
        self.folder = path.expanduser(path.join(libraryroot, identifier))
        self._needs = None

        #Create the directory for the executable files to be copied and written to.
        if not path.exists(self.folder):
            msg.okay("EXEC DIR: create {}".format(self.folder))
            mkdir(self.folder)

        #If re-run is specified, delete the fortpy.f90 file to force
        #a quick re-compile and run of the tests.
        if rerun is not None and (rerun == "*" or rerun in identifier.lower()):
            fortpath = path.join(self.folder, "fortpy.f90")
            if path.exists(fortpath):
                remove(fortpath)
            
        self.writer = MethodWriter(identifier, self.parser, self.testgen)

    def write(self, testid):
        """Writes the fortran program file for the executable specified.

        :arg testid: the identifier of the test to construct the executable for.
        """
        lines = []
        identifier = self.writer.tests[testid].identifier

        #First off, we need to start the program and set the module dependencies.
        lines.append("!!<summary>Auto-generated unit test for {}\n".format(
            self.identifier))
        lines.append("!!using FORTPY. Generated on {}.\n".format(datetime.now()))
        lines.append("!!{}</summary>\n".format(self.writer.tests[testid].description))
        lines.append("PROGRAM UNITTEST_{}\n".format(self.writer.finders[testid].executable.name))
        lines.append(self._get_uses(testid))

        #Next add the variable declarations and initializations and the calls
        #to execute the pre-req methods and the one we are trying to test.
        lines.append(self.writer.lines(testid))

        lines.append("\nEND PROGRAM UNITTEST_{}".format(self.writer.finders[testid].executable.name))

        with open(path.join(self.folder, "{}.f90".format(identifier)), 'w') as f:
            f.writelines(lines)
        
    def _get_mapping(self, mapped):
        """Gets the original file name for a module that was mapped
        when the module name does not coincide with the file name
        that the module was defined in."""
        if mapped in self.parser.mappings:
            return self.parser.mappings[mapped]
        else:
            return mapped + ".f90"

    def makefile(self, identifier):
        """Generates a makefile to create the unit testing executable
        for the specified test identifier.

        :arg identifier: the id of the test that this executable should be made for.
        """
        lines = []

        #Append the general variables
        lines.append("EXENAME\t\t= {}.x".format(identifier))
        lines.append("SHELL\t\t= /bin/bash")
        lines.append("UNAME\t\t= $(shell uname)")
        lines.append("HOSTNAME\t= $(shell hostname)")
        lines.append("LOG\t\t= compile.log")
        lines.append("")

        #Now the standard entries for ifort. We will just have the ifort include
        #file so that the MPI and other options can be tested to.
        allneeds = self.needs()
        lines.append(self._make_compiler_include(allneeds))
        lines.append(".SILENT:")
        lines.append("")

        #Append all the dependent modules to the makefile
        lines.append("LIBMODULESF90\t= \\")
        #Copy over the fortpy module in case we need it.
        lines.append("\t\tfortpy.f90 \\")
        for modk in allneeds[0:len(allneeds)-1]:
            if modk != "fortpy":
                lines.append("\t\t{} \\".format(self._get_mapping(modk)))
        lines.append("\t\t{}".format(self._get_mapping(allneeds[-1])))

        lines.append("MAINF90\t\t= {}.f90".format(identifier))
        lines.append("SRCF90\t\t= $(LIBMODULESF90) $(MAINF90)")
        lines.append("OBJSF90\t\t= $(SRCF90:.f90=.o)")
        lines.append("SLIBF90\t\t= $(LIBMODULESF90:.f90=.o)")
        lines.append("")

        #Add explicitly defined libraries that should be included when linking
        #the unit testing executable.
        linklibs = self._add_explicit_includes(lines)
        lines.append("")

        #We need to add the error handling commands to make debugging compiling easier.
        lines.append(self._make_error())
        lines.append("")

        lines.append("all:	info $(EXENAME)")
        lines.append(self._make_info())
        lines.append(self._make_exe(linklibs, identifier))
        lines[-1] += "	make -f Makefile.{}".format(identifier)

        makepath = path.join(self.folder, "Makefile.{}".format(identifier))
        with open(makepath, 'w') as f:
            f.writelines("\n".join(lines))

    def _add_explicit_includes(self, lines):
        """Adds any relevant libraries that need to be explicitly included according
        to the fortpy configuration file. Libraries are appended to the specified
        collection of lines. Returns true if relevant libraries were added.
        """
        from fortpy import config
        import sys
        includes = sys.modules["config"].includes
        linklibs = False

        if len(includes) > 0:
            lines.append("LIBS\t\t= \\")
            for library in includes:
                addlib = False
                if "modules" in library:
                    #We need to loop over the modules specified for the library and see
                    #if any of them are in our list of modules.
                    for libmod in library["modules"]:
                        if libmod.lower() in allneeds:
                            addlib = True
                            break
                else:
                    addlib = True

                if addlib:
                    linklibs = True
                    lines.append("\t\t{} \\".format(library["path"]))

        return linklibs

    def _make_compiler_include(self, allneeds):
        """Returns the include statement for the compiler to use."""
        #We need to see whether to include the pre-compiler directive or not.
        precompile = False
        for needed in allneeds:
            if self.parser.modules[needed].precompile:
                precompile = True
                break

        base = """ifeq ($(F90),ifort)
  include Makefile.ifort{}
else
ifeq ($(F90),gfortran)
  include Makefile.gfortran{}
else
  include Makefile.error
endif
endif"""
        if precompile:
            insert = ["\n  FFLAGS += -fpp -save-temps -heap-arrays", "\n  FFLAGS += -cpp"]
            return base.format(*insert)
        else:
            return base.format("", "")

    def _make_error(self):
        """Generates script for compile-time error handling."""
        return """
# Error handling
NEWFILE		= \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
ERR		= ******************************* ERROR *******************************
SHOW_LOG	= ( perl -pi -e 's/ [Ee]rror \#/\\n\\n\\n$(ERR)\\n*** error \#/' $(LOG); perl -pi -e 's/^\# 1 \"/\\n\\n$(NEWFILE)\\n\\n\\n/' $(LOG); grep -n -A3 -E "$(ERR)|$(NEWFILE)" $(LOG) )
"""

    def _make_exe(self, linklibs, identifier):
        """Generates the script to run the compiling."""
        linktxt = "$(LIBS) " if linklibs else ""
        base = """
$(EXENAME): $(OBJSF90)
	-rm $(EXENAME) 2> /dev/null
	echo -n "Linking... "
	-$(F90) $(LDFLAGS) -o $(EXENAME) $(OBJSF90) {0}>> $(LOG) 2>> $(LOG)
	echo "done."
	if test -e $(EXENAME); then echo "Produced executable: $(EXENAME)"; else $(SHOW_LOG); echo "Error."; fi

$(OBJSF90): %.o: %.f90
	echo -n "Compiling: $^... "
	-$(F90) -c $(FFLAGS) $^ >> $(LOG) 2>> $(LOG)
	echo "done."

{1}.so: $(SLIBF90)
	-rm $(EXENAME).so 2> /dev/null
	echo -n "Creating shared library..."
	-$(F90) -shared -fPIC $(FFLAGS) -o {1}.so $(SLIBF90) {0}>> $(LOG) 2>> $(LOG)

clean:
	-rm *.o *.mod *.i90 $(EXENAME)
remake:
	-rm *.o *.mod *.i90 $(EXENAME)
"""
        return base.format(linktxt, identifier)

    def _make_info(self):
        """Generates the script for displaying compile-time info."""
        module, method = self.identifier.split(".")
        return """
info: 
	echo -e "\\nCompile time:" > $(LOG)
	date >> $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo "                     FORTPY"                           | tee -a $(LOG)
	echo "               >>> version 1.1 <<<                    "| tee -a $(LOG)         
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo -e "Compiling on system  : $(UNAME)"                    | tee -a $(LOG)
	echo -e "             machine : $(HOSTNAME)"                 | tee -a $(LOG)
	echo "Compiling for module : {0}"                            | tee -a $(LOG)         
	echo "              method : {1}"                            | tee -a $(LOG)         
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo -e "DEBUG mode\\t:\\t$(DEBUG)"                          | tee -a $(LOG)
	echo -e "GPROF mode\\t:\\t$(GPROF)"                          | tee -a $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo "F90    : $(F90)"                                       | tee -a $(LOG)
	echo "FFLAGS : $(FFLAGS)"                                    | tee -a $(LOG)
	echo "LDFLAGS: $(LDFLAGS)"                                   | tee -a $(LOG)
	echo "MKLpath:$(MKL)"                                        | tee -a $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo ""                                                      | tee -a $(LOG)

""".format(module, method)

    def _get_uses(self, testid):
        """Gets a list of use statements to add to the program code."""
        #The writer can extract a dictionary of module dependencies for us.
        alluses = self.writer.uses()
        #Now we just need to generate the code statements.
        uselist = []
        for module in alluses:
            if module == self.writer.finders[testid].module.name:
                uselist.append("use {}".format(module))
            else:
                uselist.append("use {}, only: {}".format(module, ", ".join(alluses[module])))

        #Last of all, we need to append the module that handles interaction with
        #the fortran results and the python testing framework.
        if "fortpy" not in alluses:
            uselist.append("use fortpy\n")

        return self._tabjoin(uselist)

    def _tabjoin(self, values):
        """Joins a list of values with \t\n."""
        return "  {}".format("\n  ".join(values))
