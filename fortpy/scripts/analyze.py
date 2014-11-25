"""This module/script provides an interactive console interface to the test results
from multiple cases of an individual test."""
import cmd
import argparse

class FortpyShell(cmd.Cmd):
    """Shell for interacting with fortpy unit test results."""
    def __init__(self):
        cmd.Cmd.__init__(self)
        self.intro = 'Welcome to the fortpy unit test analysis shell. Type help or ? to list commands.\n'
        self.prompt = '(fortpy) '
        self.tests = {}
        """A dictionary of all the unit test directories that were parsed so far and their analyses."""
        self.active = None
        """The name of the 'module.executable' that is current active in the shell."""
        self.group = None
        """The name of the analysis group in the current unit test being worked on."""
        self._template_args = {
            "xscale": None,
            "yscale": None,
            "tfilter": None,
            "threshold": 1.,
            "xlabel": None,
            "ylabel": None,
            "functions": {},
            "independent": None,
            "dependents": [],
            "headings": []
        }
        """A dictionary of template arguments that can be set to affect the plotting/tabulating 
        behavior of the shell for a unit test analysis group. This gets duplicated for each
        analysis group and unit test worked on."""
        self.args = {}
        """A list of all the commands that are implemented in the fortpy shell."""
        self._test_cmds = ["inputs", "outputs", "timing", "group", "examine"]
        """A hard-coded list of the shell commands that require a valid test to be loaded into
        the console session. Used for validation of commands to prevent unhandled exceptions
        crashing the console. Anything that uses properties 'live', 'allvars' or 'allprops'
        relies on a valid unit test being loaded.
        """
        self._group_cmds = ["xlabel", "ylabel", "filter", "rmfilter", "threshold", "dep",
                            "indep", "rmdep", "vars", "postfix", "plot", "logplot", "loglogplot",
                            "table", "failures"]
        """As for self._test_cmds but for a valid analysis group. Anything needing property
        'curargs' relies on a valid analysis group.
        """

    @property
    def live(self):
        """Returns the currently active unit test Analysis instance's details"""
        return self.tests[self.active].details

    @property
    def allvars(self):
        """Returns all file names that can be part of the variables in the analysis."""
        return self.tests[self.active].allvars

    @property
    def allprops(self):
        """Returns the dictionary of possible properties by file name."""
        return self.tests[self.active].props

    @property
    def curargs(self):
        """Returns the dictionary of arguments for plotting/tabulating of the active
        unit test and analysis group.
        """
        return self.args[self.active][self.group]

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

    def _print_dict_files(self, case, target):
        """Prints the file names and available properties in the specified test case's 'target'
        dictionary.
        
        :arg case: the identifier of the test case to print contents for.
        :arg target: one of ['inputs', 'outputs', 'percent']
        """
        if case in self.live and target in self.live[case]:
            lines = []
            data = self.live[case][target]
            skeys = list(sorted(data.keys()))
            for key in skeys:
                props = ", ".join(sorted(data[key].keys()))
                lines.append("'{}': {}".format(key, props))
            return '\n'.join(lines)
        else:
            msg.err("Can't find the result set {} in {}".format(target, case))            

    def _complete_cases(self, text, line, istart, iend):
        """Returns the completion list of possible test cases for the active unit test."""
        if text == "":
            return self.live.keys()
        else:
            return [c for c in self.live if c.startswith(text)]

    def do_inputs(self, arg):
        """Prints a list of the available input files to analyze in the specified test case."""
        usable, filename, append = self._redirect_split(arg)
        result = self._print_dict_files(usable, "inputs")
        if result is not None:
            self._redirect_output(result, filename, append, msg.info)
    def complete_inputs(self, text, line, istart, iend):
        return self._complete_cases(text, line, istart, iend)

    def do_outputs(self, arg):
        """Prints a list of the available output files to analyze in the specified test case."""
        usable, filename, append = self._redirect_split(arg)
        result = self._print_dict_files(usable, "outputs")
        if result is not None:
            self._redirect_output(result, filename, append, msg.info)
    def complete_outputs(self, text, line, istart, iend):
        return self._complete_cases(text, line, istart, iend)

    def do_timing(self, arg):
        """Prints the number of milliseconds spent executing the specified test case."""
        usable, filename, append = self._redirect_split(arg)
        if usable in self.live and "timing" in self.live[usable]:
            if self.live[usable]["timing"] is not None:
                result = "{0:.4f} ms".format(self.live[arg]["timing"]*1000)
                self._redirect_output(result, filename, append, msg.okay)
            else:
                msg.warn("No timing data available for {}.".format(arg))
        else:
            msg.err("The test case {} does not exist.".format(arg))
    def complete_timing(self, text, line, istart, iend):
        return self._complete_cases(text, line, istart, iend)

    def _complete_tests(self, text, line, istart, iend):
        """Returns a completion list of possible, parsed unit tests that could be loaded or
        otherwise interacted with.
        """
        if text == "":
            return self.tests.keys()
        else:
            return [t for t in self.tests if t.startswith(text)]
        
    def _set_arg_generic(self, argid, arg, cast=str):
        """Sets the value of the argument with the specified id using the argument passed
        in from the shell session.
        """
        usable, filename, append = self._redirect_split(arg)
        if usable != "":
            self.curargs[argid] = cast(usable)
        if argid in self.curargs:
            result = "{}: '{}'".format(argid.upper(), self.curargs[argid])
            self._redirect_output(result, filename, append, msg.info)

    def do_xlabel(self, arg):
        """Sets the label of the x-axis for plotting."""
        self._set_arg_generic("xlabel", arg)
    def do_ylabel(self, arg):
        """Sets the label of the y-axis for plotting."""
        self._set_arg_generic("ylabel", arg)
    def do_filter(self, arg):
        """Sets the filter for the test cases to include in the plot/table by name. Only those
        test cases that include this text are included in plots, tables etc."""
        self._set_arg_generic("tfilter", arg)
    def complete_filter(self, text, line, istart, iend):
        return self._complete_cases(text, line, istart, iend)

    def do_rmfilter(self, arg):
        """Removes the test case filter that limits which results are included in plots/tables.
        """
        self.curargs["tfilter"] = None
        self._set_arg_generic("tfiletr", "")

    def do_threshold(self, arg):
        """Specify a success threshold (percent as float) that output files must attain before 
        they can be included in the plots/tables."""
        try:
            self._set_arg_generic("threshold", arg, float)
        except ValueError:
            err("The specified threshold value is not a valid float percentage. Try 1. for 100%")

    def _complete_vars(self, text):
        #Variables can come from inputs or outputs.
        if text == "":
            return self.allvars+["timing"]
        else:
            vlist = [v + "|" for v in self.allvars if v.startswith(text)]
            if text in "timing":
                vlist.append("timing")
            return vlist

    def _complete_props(self, var, text):
        if var not in self.allprops:
            return []

        if text == "":
            #We arbitrarily check the first test case to determine what properties are
            #available.
            return self.allprops[var]
        else:
            return [p for p in self.allprops[var] if p.startswith(text)]

    def _complete_fullvar(self, text, line, istart, iend):
        #Determine if we have a bar in the text, if we do then we are completing attributes;
        #otherwise we are completing variables.
        value = line.split()
        if len(value) == 1:
            value = ""
        else:
            value = value[-1]

        if "|" in value:
            var, prop = value.split("|")
            return self._complete_props(var, prop)
        else:
            return self._complete_vars(value)

    def _validate_var(self, var):
        """Validates the form of the specified variable."""
        if var == "timing":
            return True
        else:
            if "|" in var:
                varname, prop = var.split("|")
                return prop != ""
            else:
                return False
            
    def do_indep(self, arg):
        """Sets the name and attribute of the independent variable for plotting/tabulating functions."""
        if not self._validate_var(arg):
            msg.err("Variable {} is not a valid file name and property combination.".format(arg))
        else:
            self.curargs["independent"] = arg
    def complete_indep(self, text, line, istart, iend):
        return self._complete_fullvar(text, line, istart, iend)

    def do_dep(self, args):
        """Adds the name and attribute of a dependent variable to the list 
        for plotting/tabulating functions."""
        for arg in args.split():
            if not self._validate_var(arg):
                msg.err("Variable {} is not a valid file name and property combination.".format(arg))
                continue
            if arg not in self.curargs["dependents"]:
                self.curargs["dependents"].append(arg)
    def complete_dep(self, text, line, istart, iend):
        return self._complete_fullvar(text, line, istart, iend)

    def do_rmdep(self, args):
        """Removes dependent variables currently set for plotting/tabulating etc."""
        for arg in args.split():
            if arg in self.curargs["dependents"]:
                self.curargs["dependents"].remove(arg)
    def complete_rmdep(self, text, line, istart, iend):
        if text == "":
            return self.curargs["dependents"]
        else:
            return [v for v in self.curargs["dependents"] if v.startswith(text)]

    def do_vars(self, arg):
        """Prints the current settings for dependent and independent variables in the shell."""
        usable, filename, append = self._redirect_split(arg)
        result = []
        result.append("INDEPENDENT: '{}'".format(self.curargs["independent"]))
        result.append("DEPENDENTS")
        for i in range(len(self.curargs["dependents"])):
            result.append("  {}. '{}'".format(i, self.curargs["dependents"][i]))

        self._redirect_output('\n'.join(result), filename, append, msg.info)

    def do_postfix(self, arg):
        """Sets the function to apply to the values of a specific variable before plotting
        or tabulating values.
        """
        usable, filename, append = self._redirect_split(arg)
        sargs = usable.split()
        if len(sargs) == 1 and sargs[0] == "list":
            result = []
            skeys = list(sorted(self.curargs["functions"].keys()))
            for key in skeys:
                result.append("'{}' => {}".format(key, self.curargs["functions"][key]))
            self._redirect_output('\n'.join(result), filename, append, msg.info)
        elif len(sargs) == 2:
            var, fxn = sargs
            if not self._validate_var(var):
                msg.err("Variable '{}' is not a valid variable|property combination.")
        
            from importlib import import_module
            lib = fxn.split(".")
            fun = lib.pop()
            numpy = import_module('.'.join(lib))
            if not hasattr(numpy, fun):
                msg.err("Function '{}' is not a valid numpy function.")
            else:
                self.curargs["functions"][var] = fxn
                #Give the user some feedback so that they know it was successful.
                self.do_postfix("list")
    def complete_postfix(self, text, line, istart, iend):
        #We just need to look at the lib string and then return the valid functions in it
        els = line.split()
        if len(els) == 1 or (len(els) == 2 and line[-1] != " "):
            varlist = self._complete_fullvar(text, line, istart, iend)
            if text in "list" and "|" not in line:
                varlist.append("list")
            return varlist
        elif line[-1] == " ":
            return ["numpy."]
        else:
            if "." not in els[-1]:
                return ["numpy."]
            libs = els[-1].split(".")
            from importlib import import_module
            partial = libs.pop()
            prefix = '.'.join(libs)
            module = import_module(prefix)
            alldir = [a for a in dir(module) if a[0] != "_"]

            if partial == "":
                result = alldir
            else:
                result = [a for a in alldir if a.startswith(partial)]                

            return [prefix + '.' + a for a in result]

    def _plot_generic(self, filename=None):
        """Plots the current state of the shell, saving the value to the specified file
        if specified.
        """
        #Since the filename is being passed directly from the argument, check its validity.
        if filename == "":
            filename = None

        if self.curargs["xlabel"] is None:
            #Set a default x-label since we know what variable is being plotted.
            self.curargs["xlabel"] = "Value of '{}' (unknown units)".format(self.curargs["independent"])
        args = self.curargs
        a = self.tests[self.active]
        a.plot(args["independent"], args["dependents"], args["threshold"], args["xlabel"], 
               args["ylabel"], args["tfilter"], filename, args["functions"], 
               args["xscale"], args["yscale"])

    def do_logplot(self, arg):
        """Plots the current state of the shell's independent vs. dependent variables on the
        same set of axes with the y-scale set to logarithmic. Give filename to save to as
        argument or leave blank to show.
        """        
        usable, filename, append = self._redirect_split(arg)
        self.curargs["xscale"] = None
        self.curargs["yscale"] = "log"
        self._plot_generic(filename)
    def complete_logplot(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)

    def do_loglogplot(self, arg):
        """Plots the current state of the shell's independent vs. dependent variables on the
        same set of axes with the x- and y-scale set to logarithmic. Give filename to save to as
        argument or leave blank to show.
        """        
        usable, filename, append = self._redirect_split(arg)
        self.curargs["xscale"] = "log"
        self.curargs["yscale"] = "log"
        self._plot_generic(filename)
    def complete_loglogplot(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)

    def do_plot(self, arg):
        """Plots the current state of the shell's independent vs. dependent variables on the
        same set of axes. Give filename to save to as argument or leave blank to show.
        """
        usable, filename, append = self._redirect_split(arg)
        self.curargs["xscale"] = None
        self.curargs["yscale"] = None
        self._plot_generic(filename)
    def complete_plot(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)

    def _set_def_prompt(self):
        """Sets the default prompt to match the currently active unit test."""
        if len(self.active) > 15:
            module, executable = self.active.split(".")
            self.prompt = "({}*.{}*:{})".format(module[0:6], executable[0:6], self.group)
        else:
            self.prompt = "({}:{})".format(self.active, self.group)

    def do_set(self, arg):
        """Sets the specified 'module.executable' to be the active test result to interact with.
        """
        if arg in self.tests:
            self.active = arg
            #Create a default argument set and analysis group for the current plotting
            if arg not in self.args:
                self.args[arg] = {"default": dict(self._template_args)}
                self.group = "default"
            else:
                self.group = list(self.curargs.keys())[0]
                
            #Change the prompt so that they know which unit test is being edited.
            self._set_def_prompt()
        else:
            msg.err("The test case '{}' is not valid.".format(arg))
    def complete_set(self, text, line, istart, iend):
        return self._complete_tests(text, line, istart, iend)        

    def do_save(self, arg):
        """Saves the session variables to a file so that the same analysis can be continued later."""
        #We need to save the full paths to the staging directories for the tests that are loaded
        #so far; then when the session is restored, we can reparse the results.
        from os import path
        import json
        fullpath = path.expanduser(arg)

        data = {
            "tests": [path.abspath(a.stagedir) for a in self.tests.values()],
            "args": self.args
        }

        with open(fullpath, 'w') as f:
            json.dump(data, f)
        msg.okay("Saved current session to {}".format(fullpath))
    def complete_save(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)

    def do_load(self, arg):
        """Loads a saved session variables, settings and test results to the shell."""
        from os import path
        import json
        fullpath = path.expanduser(arg)
        if path.isfile(fullpath):
            with open(fullpath) as f:
                data = json.load(f)

            #Now, reparse the staging directories that were present in the saved session.
            for stagepath in data["tests"]:
                self.do_parse(stagepath)
            self.args = data["args"]
    def complete_load(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)
        
    def do_parse(self, arg, fullparse=False):
        """Parse the test results from the specified directory and load them under the name
        of 'module.executable ' that they were created with. E.g. parse classes.polya/
        """
        from os import path
        fullpath = path.expanduser(arg)
        if path.isdir(fullpath):
            if arg[-1] == "/":
                end = -2
            else:
                end = -1
            case = fullpath.split("/")[end]
            self.tests[case] = Analysis(fullpath, fullparse)
            self.do_set(case)
        else:
            msg.err("The folder {} does not exist.".format(fullpath))
    def complete_parse(self, text, line, istart, iend):
        import glob
        from os import path
        result = []
        for p in glob.glob(path.expanduser(text)+'*'):
            if path.isdir(p):
                result.append(p + "/")
            else:
                result.append(p)
        return result

    def do_fullparse(self, arg):
        """Parses a set of unit tests *and* the contents of their input and output files
        so that test results can be plotted/tabulated against the data in the input/output
        files. Takes much longer than the regular 'parse' command if the data files are large
        since all the data has to be loaded to memory.
        """
        self.do_parse(arg, True)
    def complete_fullparse(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)

    def do_reparse(self, arg):
        """Reparses the currently active unit test to get the latest test results loaded
        to the console.
        """
        #We just get the full path of the currently active test and hit reparse.
        full = arg == "full"
        from os import path
        fullpath = path.abspath(self.tests[self.active].stagedir)
        self.tests[self.active] = Analysis(fullpath, full)        
    def complete_reparse(self, text, line, istart, iend):
        possible = ["full"]
        if text == "":
            return possible
        else:
            return [p for p in possible if p.startswith(text)]

    def do_group(self, arg):
        """Creates a new analysis group with unique settings for plotting/tabulating etc.
        or switches the active group to the specified name.
        """
        if arg == "":
            msg.err("Can't switch to empty group {}.")
            return
        elif arg not in self.args[self.active]:
            self.args[self.active][arg] = dict(self._template_args)
            msg.okay("Created analysis group '{}'.".format(arg))

        self.group = arg
        self._set_def_prompt()
    def complete_group(self, text, line, istart, iend):
        if text == "":
            return self.args[self.active].keys()
        else:
            return [k for k in self.args[self.active] if k.startswith(text)]

    def do_table(self, arg):
        """Prints the set of values for the independent vs. dependent variables in the
        active unit test and analysis group as a table.
        """
        usable, filename, append = self._redirect_split(arg)
        a = self.tests[self.active]
        args = self.curargs
        result = a.table(args["independent"], args["dependents"], args["threshold"],
                         args["headings"], args["tfilter"], args["functions"])
        if result is not None:
            self._redirect_output(result, filename, append, msg.info)
    def complete_table(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)

    def do_failures(self, arg):
        """Prints a list of test cases that failed for the current unit test and analysis
        group settings. To only check failure on specific output files, set the list of
        files to check as arguments.
        """
        usable, filename, append = self._redirect_split(arg)
        a = self.tests[self.active]
        args = self.curargs
        outfiles = usable.split()
        if len(outfiles) == 0:
            outfiles = None
        result = a.failures(outfiles, args["threshold"], args["tfilter"])
        self._redirect_output(result, filename, append, msg.info)
    def complete_failures(self, text, line, istart, iend):
        #We only want to return those output files in the current unit test who *also*
        #have .compare files to check status on.
        if text == "":
            return [f for f in self.tests[self.active].comparable]
        else:
            return [f for f in self.tests[self.active].comparable if f.startswith(text)]

    def do_examine(self, arg):
        """Opens a unit test case's .out.compare file to examine the verbose comparison
        report across values.
        """
        #We use their default editor (if it has been set); otherwise we can't do much of
        #anything and issue a warning.
        from os import getenv, path, system
        testcase, output = arg.split()
        target = path.join(self.tests[self.active].stagedir, "tests", testcase, 
                           "{}.compare".format(output))
        if getenv("EDITOR") is not None:
            system("`$EDITOR {}`".format(target))
        else:
            msg.warn("$EDITOR not set in environment. Can't open {}".format(target))

    def complete_examine(self, text, line, istart, iend):
        #First we need to complete on just the test case identifier; then we complete
        #on the name of the .out file.
        args = line.split()
        if len(args) == 1 or (len(args) == 2 and line[-1] != ' '):
            return self._complete_cases(text, line, istart, iend)
        elif len(args) >= 2:
            if len(args) == 2:
                return self.complete_failures("", line, istart, iend)
            else:
                return self.complete_failures(args[2], line, istart, iend)

    def do_quit(self, arg):
        """Exit the fortpy test analysis shell. Any unsaved results will be lost."""
        print()
        return True

    def do_EOF(self, arg):
        """Exit the fortpy test analysis shell. Any unsaved results will be lost."""
        print()
        return True

    def do_history(self, arg):
        import readline
        usable, filename, append = self._redirect_split(arg)
        if usable == "":
            for i in range(1, readline.get_current_history_length()+1):
                print("{}: {}".format(i, readline.get_history_item(i)))
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
        result = []
        result.append("Fortpy Shell Console History Command Args")
        result.append("  <default>: lists the contents of the history.")
        result.append("  clear: removes all items from the history.")
        result.append("  limit: sets the maximum number of items to retain in the history.")
        return '\n'.join(result)

    @property
    def histpath(self):
        """Returns the full path to the console history file."""
        from os import path
        from fortpy import settings
        return path.join(settings.cache_directory, "history")

    def preloop(self):
        cmd.Cmd.preloop(self)
        #We need to restore the console history if it exists.
        import readline
        readline.read_history_file(self.histpath)

    def postloop(self):
        cmd.Cmd.postloop(self)
        #Save the readline console history to the fortpy cache for the next session.
        import readline
        readline.write_history_file(self.histpath)

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
            if command in self._test_cmds or command in self._group_cmds:
                #We have to test the active unit test for both the test commands and the
                #group commands, since the group commands rely on the active unit test.
                if self.active is None or self.active not in self.tests:
                    msg.err("The specified command '{}' requires an ".format(command) + 
                            "active unit test to be loaded. Use 'parse' or 'load'.")
                    return ""
                elif (command in self._group_cmds and (self.group is None or 
                                                       self.group not in self.args[self.active])):
                    msg.err("No valid analysis group is active. Use 'group' to create " 
                            "one or mark an existing one as active.")
                    return ""
                else:
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
        from os import system
        system(arg)

parser = argparse.ArgumentParser(description="Fortpy Automated Test Result Analyzer")
parser.add_argument("-pypath", help="Specify a path to add to sys.path before running the tests.")
           
#Parse the args from the commandline that ran the script, call initialize
args = vars(parser.parse_args())

#We added this argument for debugging installations. That way we can make changes
#without doing a pip install each time; just put the path to the repo root in 'pypath'
if args["pypath"]:
    import sys
    sys.path.append(args["pypath"])
    from fortpy import msg
    from fortpy.testing.parser import Analysis

if __name__ == '__main__':
    FortpyShell().cmdloop()
