from fortpy.testing.method import MethodWriter
from os import path, mkdir, remove
from datetime import datetime

class ExecutableGenerator(object):
    """Generates a fortran executable to perform unit tests for a given
    subroutine or function.

    :arg parser: an instance of the code parser for inter-module access
    :arg folder: the folder in which to generate the code files and execs.
    """
    def __init__(self, parser, folder):
        #We need the code parser because the pre-reqs specified may be in other
        #modules and may have pre-reqs of their own. This way, we can find
        #them all easily.
        self.parser = parser
        self.folder = folder

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
        for i in range(len(result)):
            module = result[i]
            self._process_module_needs(module, i, result)

        return result

    def _process_module_needs(self, module, i, result):
        """Adds the module and its dependencies to the result list."""
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
            print("FATAL: unable to find module {}.".format(module))
            exit(1)

    def reset(self, identifier, libraryroot, rerun = False):
        """Resets the writer to work with a new executable."""        
        #A method writer has all the guts needed to generate the executable
        #We are just wrapping it and making sure the variable values get
        #recorded to temporary files for comparison etc.
        self.identifier = identifier
        self.folder = path.expanduser(path.join(libraryroot, identifier))

        #Create the directory for the executable files to be copied and written to.
        if not path.exists(self.folder):
            print("EXEC DIR: create {}".format(self.folder))
            mkdir(self.folder)

        #If re-run is specified, delete the fortpy.f90 file to force
        #a quick re-compile and run of the tests.
        if rerun:
            fortpath = path.join(self.folder, "fortpy.f90")
            if path.exists(fortpath):
                remove(fortpath)
            
        self.writer = MethodWriter(identifier, self.parser)

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
        lines.append("PROGRAM UNITTEST_{}\n".format(self.writer.method.executable.name))
        lines.append(self._get_uses())

        #Next add the variable declarations and initializations and the calls
        #to execute the pre-req methods and the one we are trying to test.
        lines.append(self.writer.lines(testid))

        lines.append("\nEND PROGRAM UNITTEST_{}".format(self.writer.method.executable.name))

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

        #We are only going to use ifort for the testing framework since we have it
        #First add stuff for the C compilers
        lines.append("CC\t\t= gcc")
        lines.append("CFLAGS\t\t= -g -O3 -fPIC")
        lines.append("")

        #Now the standard entries for ifort. We will just have the ifort include
        #file so that the MPI and other options can be tested to.
        lines.append(self._make_compiler_include())
        lines.append(".SILENT:")
        lines.append("")

        #Append all the dependent modules to the makefile
        lines.append("LIBMODULESF90\t= \\")
        #Copy over the fortpy module in case we need it.
        lines.append("\t\tfortpy.f90 \\")
        allneeds = self.needs()
        for modk in allneeds[0:len(allneeds)-1]:
            lines.append("\t\t{} \\".format(self._get_mapping(modk)))
        lines.append("\t\t{}".format(self._get_mapping(allneeds[-1])))

        lines.append("MAINF90\t\t= {}.f90".format(identifier))
        lines.append("SRCF90\t\t= $(LIBMODULESF90) $(MAINF90)")
        lines.append("OBJSF90\t\t= $(SRCF90:.f90=.o)")
        lines.append("")
        lines.append("SRCC\t\t= timing.c")
        lines.append("HEADC\t\t= timing.h")
        lines.append("OBJSC\t\t= timing.o")

        #We need to add the error handling commands to make debugging compiling easier.
        lines.append(self._make_error())
        lines.append("")

        lines.append("all:	info $(EXENAME)")
        lines.append(self._make_info())
        lines.append(self._make_exe())

        makepath = path.join(self.folder, "Makefile.{}".format(identifier))
        with open(makepath, 'w') as f:
            f.writelines("\n".join(lines))

    def _make_compiler_include(self):
        """Returns the include statement for the compiler to use."""
        return """ifeq ($(F90),ifort)
  include Makefile.ifort
else
ifeq ($(F90),gfortran)
  include Makefile.gfortran
else
  include Makefile.error
endif
endif"""

    def _make_error(self):
        """Generates script for compile-time error handling."""
        return """
# Error handling
NEWFILE		= \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
ERR		= ******************************* ERROR *******************************
SHOW_LOG	= ( perl -pi -e 's/ [Ee]rror \#/\\n\\n\\n$(ERR)\\n*** error \#/' $(LOG); perl -pi -e 's/^\# 1 \"/\\n\\n$(NEWFILE)\\n\\n\\n/' $(LOG); grep -n -A3 -E "$(ERR)|$(NEWFILE)" $(LOG) )
"""

    def _make_exe(self):
        """Generates the script to run the compiling."""
        return """
$(EXENAME): $(OBJSF90) $(OBJSC)
	-rm $(EXENAME) 2> /dev/null
	echo -n "Linking... "
	-$(F90) $(LDFLAGS) -o $(EXENAME) $(OBJSC) $(OBJSF90) >> $(LOG) 2>> $(LOG)
	echo "done."
	if test -e $(EXENAME); then echo "Produced executable: $(EXENAME)"; else $(SHOW_LOG); echo "Error."; fi

$(OBJSF90): %.o: %.f90
	echo -n "Compiling: $^... "
	-$(F90) -c $(FFLAGS) $^ >> $(LOG) 2>> $(LOG)
	echo "done."

$(OBJSC): %.o: %.c $(HEADC)
	echo -n "Compiling: $^... "
	$(CC) $(CFLAGS) -c $^
	echo "done."
"""

    def _make_info(self):
        """Generates the script for displaying compile-time info."""
        return """
info: 
	echo -e "\\nCompile time:" > $(LOG)
	date >> $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo "                     FORTPY"                           | tee -a $(LOG)
	echo "               >>> version 1.0 <<<                    "| tee -a $(LOG)         
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo -e "Compiling on system  : $(UNAME)"                    | tee -a $(LOG)
	echo -e "             machine : $(HOSTNAME)"                 | tee -a $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo -e "DEBUG mode\\t:\\t$(DEBUG)"                            | tee -a $(LOG)
	echo -e "MPI mode\\t:\\t$(UNCLE_MPI)"                          | tee -a $(LOG)
	echo -e "openMP mode\\t:\\t$(UNCLE_OMP)"                       | tee -a $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo "F90    : $(F90)"                                       | tee -a $(LOG)
	echo "FFLAGS : $(FFLAGS)"                                    | tee -a $(LOG)
	echo "CC     : $(CC)"                                        | tee -a $(LOG)
	echo "CFLAGS : $(CFLAGS)"                                    | tee -a $(LOG)
	echo "LDFLAGS: $(LDFLAGS)"                                   | tee -a $(LOG)
	echo "MKLpath:$(MKL)"                                        | tee -a $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo ""                                                      | tee -a $(LOG)

"""

    def _get_uses(self):
        """Gets a list of use statements to add to the program code."""
        #The writer can extract a dictionary of module dependencies for us.
        alluses = self.writer.uses()
        #Now we just need to generate the code statements.
        uselist = []
        for module in alluses:
            if module == self.writer.method.module.name:
                uselist.append("use {}".format(module))
            else:
                uselist.append("use {}, only: {}".format(module, ", ".join(alluses[module])))

        #Last of all, we need to append the module that handles interaction with
        #the fortran results and the python testing framework.
        uselist.append("use fortpy\n")

        return self._tabjoin(uselist)

    def _tabjoin(self, values):
        """Joins a list of values with \t\n."""
        return "  {}".format("\n  ".join(values))
