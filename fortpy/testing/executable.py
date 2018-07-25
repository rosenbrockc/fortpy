from .. import msg
from fortpy.testing.method import MethodWriter
from os import path, mkdir, remove
from datetime import datetime
from fortpy.interop.make import makefile

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
        self.libraryroot = folder
        """Specifies a staging directory to override the one specified in the <group>
        tag (if one was specified there)."""
        self.testgen = testgen
        self.folders = {}
        """Dictionary of full paths to the staging directory for each unit test
        being generated by the framework.
        """
        self.folder = None
        """The full path to the *current* staging directory for the executable whose
        tests are being generated.
        """
        self._needs = None

    def needs(self):
        """Returns a list of all the modules that this executable needs to
        run correctly."""
        if self._needs is None:
            from fortpy.code import order_module_dependencies
            self._needs = order_module_dependencies(self.writer.uses(), self.parser)
        return self._needs

    def reset(self, identifier, coderoot):
        """Resets the writer to work with a new executable."""        
        #A method writer has all the guts needed to generate the executable
        #We are just wrapping it and making sure the variable values get
        #recorded to temporary files for comparison etc.
        self.writer = MethodWriter(identifier, self.parser, self.testgen,
                                   stagedir=self.libraryroot)
        self.identifier = identifier
        if self.libraryroot is not None:
            self.folders[identifier] = path.abspath(self.libraryroot)
        elif self.writer.group is not None and self.writer.group.staging is not None:
            #Determine the absolute path using the relative path specified in the group
            #staging directory attribute.
            from fortpy.tramp import coderelpath
            self.folders[identifier] = coderelpath(coderoot, self.writer.group.staging)
        else:
            raise ValueError("Can't determine the staging directory for running the tests.")
        
        self._needs = None
        #The calling methods need to know whether the context (i.e. the staging folder)
        #of this executable writer has changed. This variable tracks that.
        newstage = self.folder != self.folders[identifier]
        self.folder = path.join(self.folders[identifier], identifier)

        #Create the directory for the executable files to be copied and written to.
        if not path.exists(self.folder):
            msg.okay("EXEC DIR: create {}".format(self.folder))
            mkdir(self.folder)

        return newstage
                
    def write(self, testid, coderoot):
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
        lines.append("  implicit none\n")

        #Next add the variable declarations and initializations and the calls
        #to execute the pre-req methods and the one we are trying to test.
        lines.append(self.writer.lines(testid, coderoot))

        lines.append("\nEND PROGRAM UNITTEST_{}".format(self.writer.finders[testid].executable.name))

        with open(path.join(self.folder, "{}.f90".format(identifier)), 'w') as f:
            f.writelines(lines)
        
    def makefile(self, identifier):
        """Generates a makefile to create the unit testing executable
        for the specified test identifier.

        :arg identifier: the id of the test that this executable should be made for.
        """
        allneeds = self.needs()

        #We need to see whether to include the pre-compiler directive or not.
        precompile = False
        for needed in allneeds:
            if needed=="fortpy" or needed=="fpy_auxiliary":
                continue
            #Skip the external modules in pre-compiled libraries.
            if needed.lower() in self.parser.externals:
                continue
            
            if self.parser.modules[needed].precompile:
                precompile = True
                break

        lines = []
        #Determine whether the compilation should produce headers on the output screen.
        from fortpy.msg import will_print
        verbose = will_print(2)
        makepath = path.join(self.folder, "Makefile.{}".format(identifier))
        makefile(identifier, allneeds, makepath, self.identifier, precompile,
                 parser=self.parser, inclfpyaux=self.writer.autoclass, verbose=verbose)
        
    def _get_uses(self, testid):
        """Gets a list of use statements to add to the program code."""
        #The writer can extract a dictionary of module dependencies for us.
        alluses = self.writer.uses()
        #Now we just need to generate the code statements.
        uselist = []
        for module in alluses:
            if module == self.writer.finders[testid].module.name:
                uselist.append("use {}".format(module))
            elif module == "fpy_auxiliary":
                uselist.append("use fpy_auxiliary")
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
