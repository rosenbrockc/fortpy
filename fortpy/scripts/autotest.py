"""Shell interface for automating the creation of unit tests with
minimal user input."""
import cmd
import argparse

class WizardShell(cmd.Cmd):
    """Shell for interacting with fortpy unit test results."""
    def __init__(self, start=None):
        """
        :arg start: the starting code folder/file to auto-test."""
        cmd.Cmd.__init__(self)
        self.intro = 'Welcome to the fortpy unit test automation shell. Type help or ? to list commands.\n'
        self.prompt = '(fortpy) '
        self.stagedir = None
        """The directory to stage the tests in."""
        self.folders = {}
        """Dict of the default folder paths (relative to code root) for each of the modules
        that are automated.
        """
        self.wizards = {}
        """The automation.Wizard instances that stores the result of the user interaction.
        Dict keys are the code library paths.
        """
        self.active = None
        """The name of the active module being automated."""
        self.parser = None
        """An instance of fortpy.code.CodeParser for interacting with code elements."""
        self._maxerr = 5
        """The maximum number of times that the command loop can error out before the script
        quits."""
        self._errcount = 0
        """The number of times the shell has caught an unhandled exception so far."""
        self.lasterr = None
        """The last unhandled exception caught by the shell."""
        self.xauto = None
        """The fortpy.elements.Executable instance of the executable that we are automating
        tests for at the moment.
        """
        self.tauto = None
        """The fortpy.testing.elements.TestSpecification instance that is being edited
        for the currently active self.xauto.
        """
        self.pauto = None
        """The fortpy.elements.ValueElement instance of the parameter that we are automating
        data for at the moment.
        """
        self.caseauto = None
        """The identifier of the case that is currently being edited in the active self.tauto.
        """
        global input
        try: input = raw_input
        except NameError: pass
        self.input = input
        """Access to the input function (py3+) or raw_input function (py2.x) for getting
        values from the user.
        """
        
        if start is not None:
            self._load_library(start)
        
    def _fixed_width_info(self, lines, pfun=None):
        """Prints the specified string as information with fixed width of 80 chars."""
        if pfun is None:
            pfun = msg.info
        for string in lines:
            for line in [string[i:i+80] for i in range(0, len(string), 80)]:
                pfun(line)
            msg.blank()

    def _check_parser(self):
        """Checks whether we have a functional code parser."""
        if self.parser is None:
            msg.err("Cannot set a module for unit-testing. First load one with 'auto'.")
            return False
        else:
            return True

    def _redirect_split(self, args):
        """Determines whether the specified shell args have a redirect specification for
        the output. Returns a tuple of (usable args, filename, append).
        """
        if ">>" in args:
            append = True
            usable, filename = args.split(">>")
        elif ">" in args:
            append = False
            usable, filename = args.split(">")
        else:
            append = filename = None
            usable = args

        if filename is not None:
            return (usable, filename.strip(), append)
        else:
            return (usable, filename, append)

    def _redirect_output(self, value, filename=None, append=None, printfun=None):
        """Outputs the specified value to the console or a file depending on the redirect
        behavior specified.

        :arg value: the string text to print or save.
        :arg filename: the name of the file to save the text to.
        :arg append: when true, the text is appended to the file if it exists.
        """
        if filename is None:
            if printfun is None:
                print(value)
            else:
                printfun(value)
        else:
            if append:
                mode = 'a'
            else:
                mode = 'w'

            from os import path
            fullpath = path.abspath(filename)
            with open(filename, mode) as f:
                f.write(value + '\n')

    def _xclear(self):
        """Clears the terminal and re-prints the current status of the executable
        being automated.
        """
        self.do_clear(None)
        msg.info(self.xauto.signature.split('(')[0])
        self._fixed_width_info([self.xauto.summary], msg.std)
        msg.blank(1,0)

    def _pclear(self):
        """Clears the terminal and re-prints the current status of the parameter
        being automated.
        """
        self.do_clear(None)
        msg.info(str(self.pauto))
        self._fixed_width_info([self.pauto.summary], msg.std)
        msg.blank(1,0)
                
    def do_help(self, arg):
        """Sets up the header for the help command that explains the background on how to use
        the script generally. Help for each command then stands alone in the context of this
        documentation. Although we could have documented this on the wiki, it is better served
        when shipped with the shell.
        """
        if arg == "":
            lines = [("This shell works with the fortpy unit testing framework to automate the "
                      "*creation* of unit tests (in addition to automating the execution). Since "
                      "the fortpy unit tests rely on the existence of XML files with details on "
                      "how to perform the tests, this amounts to automating the creation of those "
                      "files."),
                     ("In most real-life codes, there is a small set of input parameters for the "
                      "code as a whole. The multitude of other executables in the various modules "
                      "consume the output of previously run executables as inputs. This shell "
                      "automates the connecting of the various executables' parameters together "
                      "so that duplication is avoided, and so that adding one extra test-case to "
                      "the first executables in the chain automatically addes test cases to others "
                      "that depend on that one.")]
            
            self._fixed_width_info(lines)
        cmd.Cmd.do_help(self, arg)
        
    def _set_def_prompt(self):
        """Sets the default prompt to match the currently active unit test."""
        if len(self.active) > 15:
            self.prompt = "(fpy:{})".format(self.active[0:16])
        else:
            self.prompt = "(fpy:{})".format(self.active)
    def _complete_modules(self, text, line, istart, iend):
        """Returns a completion list of possible, parsed unit tests that could be loaded or
        otherwise interacted with.
        """
        if not self._check_parser():
            return []
        if text == "":
            return list(self.parser.modules.keys())
        else:
            return [m for m in self.parser.modules if m.startswith(text)]

    def _complete_executables(self, text, line, istart, iend):
        """Returns a completion list of executables for a 'module.executable'.
        """
        if "." in text:
            modname, rest = text.strip().split('.')
            module = self.parser.get(modname)
            if module is not None:
                if rest != "":
                    kset = [x for x in module.executables.keys()
                                   if x.startswith(rest)]
                else:
                    kset = module.executables.keys()
                return ["{}.{}".format(modname, x) for x in kset]
            else:
                msg.warn("The module '{}' does not exist or can't be found.".format(modname))
        else:
            suggest = self._complete_modules(text, line, istart, iend)
            if len(suggest) > 1:
                return suggest
            else:
                return [suggest[0] + "."]

    def _complete_params(self, text, line, istart, iend):
        """Returns the completion list for the parameters of the active executable.
        """
        from fortpy.elements import Executable
        if self.xauto is None or not isinstance(self.xauto, Executable):
            return []
        
        if text == "":
            return list(self.xauto.parameters.keys())
        else:
            return [p for p in self.xauto.parameters if p.startswith(text)]
            
    def _get_ordered_xsuggestions(self):
        """Returns a list of executables in the active module, ordered by ease
        of testability.
        """
        from fortpy.stats.testing import testability
        scores = testability(self.parser, self.active)
        from operator import itemgetter
        skeys = sorted([(k, v["score"]) for (k,v) in scores.items()], key=itemgetter(1), reverse=True)
        return [(k, scores[k]) for (k, v) in skeys]
        
    def do_testability(self, arg):
        """Prints a summary of testability scores for the active module.
        """
        if self.active is not None:
            if arg.strip() == "":
                fmtstr = "{0:<40s} | {1:^10.2f} | {2:^12.2f} | {3:^5.2f} |"
                headstr = "{0:<40s} | {1:^10s} | {2:^12s} | {3:^5s} |"
                dheader = headstr.format(*("Executable Identifier", "Dep. Score", "Param. Score", "Total"))
            else:
                fmtstr = "{0:^3d} | {1:<40s} | {2:^10.2f} | {3:^12.2f} | {4:^5.2f} |"
                headstr = "{0:^3s} | {1:<40s} | {2:^10s} | {3:^12s} | {4:^5s} |"
                dheader = headstr.format(*("#", "Executable Identifier", "Dep. Score",
                                           "Param. Score", "Total"))
                
            msg.info(dheader)
            msg.std(''.join(['-' for i in range(len(dheader))]))
            ce, cs, cw, co = msg.cenum["cerrs"], msg.cenum["cstds"], msg.cenum["cwarn"], msg.cenum["cokay"]

            i = 0
            for (sname, score) in self._get_ordered_xsuggestions():
                stotal = score["score"]
                pscore = score["pscore"]
                if pscore > 0.95:
                    pcol = co
                elif pscore > 0.5:
                    pcol = cw
                else:
                    pcol = ce

                dscore = score["dscore"]
                if abs(dscore) < 0.05:
                    dcol = co
                elif abs(dscore) < 0.5:
                    dcol = cw
                else:
                    dcol = ce

                if stotal > 0.95:
                    tcol = co
                elif stotal > 0.5:
                    tcol = cw
                else:
                    tcol = ce

                if arg.strip() == "":
                    cols = (cs, dcol, pcol, tcol)
                    text = fmtstr.format(sname, dscore, pscore, stotal)
                else:
                    cols = (cs, cs, dcol, pcol, tcol)
                    text = fmtstr.format(i, sname, dscore, pscore, stotal)
                msg.arb(text, cols, '|')
                i += 1
        else:
            msg.warn("There is no active module to print a testability summary for.")
    def help_testability(self):
        lines = [("Prints a summary of the automation testability for all the executables "
                  "in the currently active module."),
                 ("See Also: 'set' for setting the currently active module.")]
        self._fixed_width_info(lines)        
            
    def do_summary(self, arg):
        """Prints a summary of the unit test coverage for the active module.
        """
        from fortpy.stats.testing import summary, dheader
        if self.active is not None:
            msg.info(dheader)
            msg.info(''.join(['-' for i in range(len(dheader))]))
            ce, cs, cw, co = msg.cenum["cerrs"], msg.cenum["cstds"], msg.cenum["cwarn"], msg.cenum["cokay"]
            execsum = summary(self.parser, self.active)
            for xname in sorted(execsum.keys()):
                analysis, text = execsum[xname]
                if analysis["skip"]:
                    cols = (cs, cs, cs, cs)
                else:
                    if analysis["ncases"] >= 2*analysis["ntests"]:
                        casecol = co
                    elif analysis["ncases"] >= analysis["ntests"]:
                        casecol = cw
                    else:
                        casecol = ce

                    if analysis["ntests"] > 0:
                        testcol = co
                    else:
                        testcol = ce

                    if analysis["docsum"] == "GOOD":
                        doccol = co
                    elif analysis["docsum"] == "OK":
                        doccol = cw
                    else:
                        doccol = ce
                    cols = (testcol, testcol, casecol, doccol)

                msg.arb(text, cols, '|')
        else:
            msg.warn("There is no active module to print a coverage summary for.")
    def help_summary(self):
        lines = [("Prints a summary of the unit testing coverage for all the executables "
                  "in the currently active module."),
                 ("See Also: 'set' for setting the currently active module.")]
        self._fixed_width_info(lines)        
        
    def do_set(self, arg):
        """Sets the specified f90 file to be the active one for automation.
        """
        if not self._check_parser():
            return
        if arg not in self.parser.modules:
            self.parser.load_dependency(arg, True, True)
            
        if arg in self.parser.modules:
            self.active = arg               
            #Change the prompt so that they know which module is being automated.
            self._set_def_prompt()
        else:
            msg.err("The module '{}' is not valid or cannot be found.".format(arg))
    def complete_set(self, text, line, istart, iend):
        return self._complete_modules(text, line, istart, iend)        
    def help_set(self):
        lines = [("Sets the currently active module whose executables are having their "
                  "unit tests automated."),
                 ("EXAMPLE \"set utilities\" sets the automation target to "
                  "the 'utilities' module."),
                 ("NOTE: you can see which module is currently active by looking at the prompt.")]
        self._fixed_width_info(lines)

    def _load_library(self, arg):
        """Loads the specified code library and f90 file (if specified) as the targets of
        the automated test generation.
        """
        from os import path
        fullpath = path.abspath(path.expanduser(arg))
        folder = None
        if path.isfile(fullpath):
            folder, f90 = path.split(fullpath)
        if path.isdir(fullpath):
            folder, f90 = fullpath, None
            
        if folder.lower() not in self.wizards:
            self.wizards[folder.lower()] = Wizard(fullpath)

        if self.parser is None:
            from fortpy.code import CodeParser
            self.parser =CodeParser()

        allfiles = {}
        self.parser.scan_path(folder, allfiles)
        for fname, fpath in allfiles.items():
            self.parser.parse(fpath, True, True)

        if f90 is not None:
                if f90 in self.parser._modulefiles:
                    modname = self.parser._modulefiles[f90]
                    self.do_set(modname)
                else:
                    msg.warn("The parser couldn't find a valid fortran module in {}".format(arg))
                
        if folder is None:
            msg.err("The file/folder {} does not exist.".format(fullpath))

    def do_parameter(self, value):
        """Prompts for the hookup of the parameter with the specified name.
        """
        #We want to print the summary information for the parameter (if any
        #exists) so that they know what they are working with.
        if value in self.xauto.parameters:
            self.pauto = self.xauto.parameters[value]
            if "in" in self.pauto.direction:
                from fortpy.testing.automation import prompt_input
                prompt_input(self, self.pauto)
            else:
                
            self._pclear()
        else:
            msg.warn("can't locate parameter '{}' in '{}'.".format(value, self.xauto.name))
        
    def complete_parameter(self, text, line, istart, iend):
        return self._complete_params(text, line, istart, iend)
    def help_parameter(self):
        lines = [("Starts a prompt for setting the data targets for input and "
                  "output parameters for the active executable."),
                 ("EXAMPLE \"parameter inmat\" starts an automation prompt for the "
                  "parameter 'inmat'."),
                 ("NOTE: you can see which executable is currently active by looking at the prompt.")]
        self._fixed_width_info(lines)
            
    def _prompt_executable(self, value):
        """Prompts the user to set the unit tests up for the executable with the specified
        numerical id or name.
        """
        if isinstance(value, int):
            suggest = self._get_ordered_xsuggestions()
            if value < len(suggest):
                xkey = suggest[value][0].lower()
                self.xauto = self.parser.get_executable(xkey)
                if self.xauto is None:
                    msg.warn("the executable '{}' could not be found.".format(xkey))

            else:
                msg.warn("the choice {} is not valid; must be less than {}.".format(value, len(suggest)))
        else:
            self.xauto = self.parser.get_executable(value)
            if self.xauto is None:
                msg.warn("the executable '{}' could not be found.".format(value))

        #Finally, set the prompt so the user knows what they are editing.
        if self.xauto is not None:
            self._xclear()
            self.prompt = "({}:{})".format(self.active, self.xauto.name[0:min(15, len(self.xauto.name))])
            
    def do_auto(self, arg):
        """Starts an interactive automation session for the code directory specified in args."""
        if arg != "":
            self._prompt_executable(arg)
        else:
            self.do_clear(None)
            msg.std("Please select an executable to automate:")
            msg.std(" - A low parameter score means more input parameters to define data for.")
            msg.std(" - A low dependency score means that the executable calls other ")
            msg.std("   executables in the module for which unit tests have not yet ")
            msg.std("   been defined.")
            msg.blank(1, -1)
            self.do_testability("A")
            msg.blank(1,-1)

            value = None
            i = 0
            while value is None and i < 3:
                choice = self.input("Your choice: ")
                import re
                if re.match("[\d]+", choice.strip()):
                    value = int(choice)
                else:
                    msg.warn("Invalid choice; please try again", prefix=False)
                i += 1

            if value is None:
                return
            else:
                self._prompt_executable(value)
    def complete_auto(self, text, line, istart, iend):
        return self._complete_executables(text, line, istart, iend)
    def help_auto(self):
        lines = [("Starts an interactive automation session for the active module. If no "
                  "'module.executable' key is provided, a list of executables without any "
                  "tests is displayed in order of decreasing ease of testability."),
                 ("EXAMPLE: \"auto bcs.do_wrapped\" prompts for setting up the tests for "
                  "the 'do_wrapped' executable in module 'bcs'.")]
        self._fixed_width_info(lines)
            
    def do_load(self, arg):
        """Loads the code files in the specified library and (optionally) sets the active module.
        """
        self._load_library(arg)
            
    def complete_load(self, text, line, istart, iend):
        import glob
        from os import path

        result = []
        command = line.split()[0]
        rest = line[len(command):len(line)].strip()
        tilde = path.expanduser("~")

        for p in glob.glob(path.expanduser(rest)+'*'):
            if "~" in rest:
                suggest = p.replace(tilde, "~")
            else:
                suggest = p
            if path.isdir(p):
                result.append(suggest + "/")
            else:
                result.append(suggest)

        if rest == text:
            return result
        else:
            return [text + r[len(rest):len(r)] for r in result]
    def help_load(self):
        lines = [("Loads the modules in the specified code library so that their unit "
                  "tests can be automated. Sets the default module for automation to "
                  "the one specified."),
                 ("EXAMPLE \"auto symlib/src/utilities.f90\" parses all the modules in "
                  "symlib 'src/' directory and sets the 'utilities.f90' file as the active "
                  "module whose tests are being automated.")]
        self._fixed_width_info(lines)

    def _level_up(self):
        """Exits the current automation context if we are any layers deep.
        """
        if self.pauto is not None:
            self.pauto = None
            self._xclear()
            return

        if self.xauto is not None:
            self.xauto = None
            print()
            self.do_clear(None)
            self._set_def_prompt()
            return

        print()
        return True
        
    def do_quit(self, arg):
        """Exit the fortpy test automation shell. Any unsaved results will be lost."""
        return self._level_up()

    def do_EOF(self, arg):
        """Exit the fortpy test automation shell. Any unsaved results will be lost."""
        return self._level_up()

    def do_history(self, arg):
        import readline
        usable, filename, append = self._redirect_split(arg)
        if usable == "":
            for i in range(1, readline.get_current_history_length()+1):
                print(("{}: {}".format(i, readline.get_history_item(i))))
        elif usable == "clear":
            readline.clear_history()
        elif "limit" in usable:
            try:
                length = int(usable.split()[1])
                readline.set_history_length(length)
            except ValueError:
                msg.err("The maximum history length you entered is not a valid integer.")
    def complete_history(self, text, line, istart, iend):
        possible = ["clear", "limit"]
        if text == "":
            return possible
        else:
            return [p for p in possible if p.startswith(text)]
    def help_history(self):
        lines = [("Commands for interacting with the fortpy shell history. The history "
                  "functions similarly to the bash shell and commands that work there "
                  "(such as Ctrl-R for reverse i-search) will work in the fortpy shell "
                  "as well. At the end of a fortpy session, the history is saved to a "
                  "file called 'history' in the fortpy cache directory specified in the "
                  "fortpy settings."),
                 ("EXAMPLE: \"history\" lists all the items currently in the history."),
                 ("EXAMPLE: \"history clear\" deletes all the items from the history."),
                 ("EXAMPLE: \"history limit 10000\" limits the number of items stored in "
                  "the history to 10000. By default, sequential duplicates are not stored  "
                  "in the fortpy shell history.")]
        self._fixed_width_info(lines)

    @property
    def histpath(self):
        """Returns the full path to the console history file."""
        from os import path
        from fortpy import settings
        return path.join(settings.cache_directory, "autotest-history")

    def preloop(self):
        cmd.Cmd.preloop(self)
        #We need to restore the console history if it exists.
        import readline
        from os import path
        if path.isfile(self.histpath) and self.lasterr is None:
            readline.read_history_file(self.histpath)

    def postloop(self):
        cmd.Cmd.postloop(self)
        #Save the readline console history to the fortpy cache for the next session.
        import readline
        readline.write_history_file(self.histpath)

    def _store_lasterr(self):
        """Stores the information about the last unhandled exception."""
        from sys import exc_info
        from traceback import format_exception
        e = exc_info()
        self.lasterr = '\n'.join(format_exception(e[0], e[1], e[2]))

    # def cmdloop(self):
    #     try:
    #         cmd.Cmd.cmdloop(self)
    #     except Exception as exsimple:
    #         msg.err(exsimple.message)
    #         self._store_lasterr()
    #         if self._errcount < self._maxerr:
    #             self._errcount += 1
    #             msg.err("The shell has caught {} unhandled exceptions so far.\n".format(self._errcount) + 
    #                     "When that value reaches {}, the shell will save a ".format(self._maxerr) + 
    #                     "recovery file and exit.")
    #             self.postloop()
    #             self.cmdloop()
    #         else:
    #             self.do_save("#fortpy.shell#")
    #             msg.err("Something unexpected happened. The shell has died. Your session "
    #                     "has been saved as '#fortpy.shell#' in the current directory.")
        
    def precmd(self, line):
        """Makes sure that the command specified in the line is valid given the current
        status of loaded unit tests and analysis group.
        """
        if line == "":
            return ""
        command = line.split()[0]
        if "!" in command:
            value = command.split("!")[1]
            try:
                ihist = int(value)
                import readline
                if ihist <= readline.get_current_history_length():
                    return readline.get_history_item(ihist)
                else:
                    msg.warn("The specified history item {} ".format(ihist) + 
                             "does not exist in the history.")
                    return ""
            except ValueError:
                #This clearly isn't a bash style recursion of a past history item.
                #Just run the command as it was originally entered.
                return line
        else:
            return line

    def emptyline(self):
        """Prevents the last non-empty command from being run when the line is empty."""
        msg.info("No command entered. Type help or ? for command listing.")

    def do_clear(self, arg):
        """Clears the screen."""
        from os import system
        system("clear")

    def do_shell(self, arg):
        """Executes a bash command from within this shell."""
        from os import system
        system(arg)

    def do_cd(self, arg):
        """Imitates the bash shell 'cd' command."""
        from os import chdir, path
        fullpath = path.abspath(path.expanduser(arg))
        if path.isdir(fullpath):
            chdir(fullpath)
        else:
            msg.err("'{}' is not a valid directory.".format(arg))
    def complete_cd(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)

    def do_ls(self, arg):
        """Imitates the bash 'ls' command."""
        self.do_shell("ls")

    def do_pwd(self, arg):
        """Imitates the bash 'pwd' command."""
        self.do_shell("pwd")

    def do_error(self, arg):
        usable, filename, append = self._redirect_split(arg)
        if usable == "":
            if self.lasterr is not None:
                self._redirect_output(self.lasterr, filename, append, msg.warn)
            else:
                msg.okay("No uncaught exceptions on record.")
        elif usable == "clear":
            self.lasterr = None
            msg.okay("Cleared the last uncaught exception.")            
    def complete_error(self, text, line, istart, iend):
        possible = ["clear"]
        return [p for p in possible if p.startswith(text)]
    def help_error(self):
        lines = [("Prints the last uncaught exception generated in the shell. If 'clear' is "
                  "passed as a single argument, the last exception is cleared from the shell.")]
        self._fixed_width_info(lines)
        
parser = argparse.ArgumentParser(description="Fortpy Automated Test Result Analyzer")
parser.add_argument("-pypath", help="Specify a path to add to sys.path before running the tests.")
parser.add_argument("source", help="Specify the code library (and/or file) to auto-test.")
#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_args())

#We added this argument for debugging installations. That way we can make changes
#without doing a pip install each time; just put the path to the repo root in 'pypath'
if args["pypath"]:
    import sys
    sys.path.append(args["pypath"])
from fortpy import msg
from fortpy.testing.automation import Wizard

if __name__ == '__main__':
    msg.set_verbosity(0)
    shell = WizardShell(args["source"])
    shell.cmdloop()
