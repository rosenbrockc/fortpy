#!/usr/bin/env python
"""This module/script provides an interactive console interface to the test results
from multiple cases of an individual test."""
import cmd
import argparse
import re
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
            "version": 7,
            "xscale": None,
            "yscale": None,
            "tfilter": ["*"],
            "threshold": 1.,
            "functions": {},
            "independent": None,
            "dependents": [],
            "headings": [],
            "fits": {},
            "labels": {},
            "colors": {},
            "fonts": {},
            "markers": {},
            "ticks": {},
            "lines": {},
            "plottypes": {},
            "limits": {},
            "twinplots": {},
            "legend": {}
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
        self._group_cmds = ["label", "filter", "rmfilter", "threshold", "dep",
                            "indep", "rmdep", "vars", "postfix", "failures", "headings"]
        """As for self._test_cmds but for a valid analysis group. Anything needing property
        'curargs' relies on a valid analysis group.
        """
        self._var_cmds = ["plot", "logplot", "loglogplot", "table", "fit", "color", "label",
                          "font", "fontsave", "fontload", "ticks", "lines", "markers", "limit"]
        """List of commands that require an independent variable to be set and at least one
        dependent variable to be set."""
        self._maxerr = 10
        """The maximum number of times that the command loop can error out before the script
        quits."""
        self._errcount = 0
        """The number of times the shell has caught an unhandled exception so far."""
        self.lasterr = None
        """The last unhandled exception caught by the shell."""
        self._version_check = {}
        """Keeps track of which unit test and analysis group argument dictionaries have been checked
        for version upgrade when the shell loads each one."""
        self._possible_cols = {
            "blue": "b",
            "green": "g",
            "red": "r",
            "cyan": "c",
            "magenta": "m",
            "yellow": "y",
            "black": "k",
            "white": "w"
        }
        """A dict of possible matplotlib colors and their corresponding codes."""
        self.legend_locs = {'upper_right'  : 1,
                            'upper_left'   : 2,
                            'lower_left'   : 3,
                            'lower_right'  : 4,
                            'right'        : 5,
                            'center_left'  : 6,
                            'center_right' : 7,
                            'lower_center' : 8,
                            'upper_center' : 9,
                            'center'       : 10}
        """Dict of legend locations and their corresponding numerical values."""
        self.legend_opts = {'loc': self.legend_locs, 'numpoints': int, 'markerscale': float, 'markerfirst': bool,
                            'scatterpoints': int, 'scatteryoffsets': (list, float),
                            'fontsize': float, 'borderpad': float, 'labelspacing': float, 'handlelength': float,
                            'handleheight': float, 'handletextpad': float, 'borderaxespad': float,
                            'columnspacing': float, 'ncol': int, 'fancybox': bool, 'shadow': bool,
                            'title': str, 'framealpha': float, 'frameon': bool}
        """Dict of legend options and their corresponding types expected by the plotter."""
        from matplotlib.markers import MarkerStyle
        self._possible_markers = MarkerStyle.markers 
        """A dict of possible marker styles to use on plots."""

        from matplotlib.lines import Line2D
        self._possible_linestyles = Line2D.lineStyles
        """A dict of possible line styles for plots with lines."""

    def do_help(self, arg):
        """Sets up the header for the help command that explains the background on how to use
        the script generally. Help for each command then stands alone in the context of this
        documentation. Although we could have documented this on the wiki, it is better served
        when shipped with the shell.
        """
        if arg == "":
            lines = [("The fortpy unit testing analysis shell makes it easy to analyze the results "
                      "of multiple test cases, make plots of trends and tabulate values for use in "
                      "other applications. This documentation will provide an overview of the basics. "
                      "Use 'help <command>' to get specific command help."),
                     ("Each fortpy shell session can hold the results of multiple unit tests. You can "
                      "load a unit test's results into the session using one of the 'parse' commands. "
                      "Once the test is loaded you can tabulate and plot results by setting test case "
                      "filters ('filter'), and independent and dependent variables ('indep', 'dep')."
                      "Switch between different unit tests loaded into the session using 'set'."),
                     ("To make multiple plots/tables for the same unit test, create new analysis "
                      "groups ('group'). "
                      "Each group has its own set of properties that can be set (e.g. variables, plot "
                      "labels for axes, filters for test cases, etc.) The possible properties that affect "
                      "each command are listed in the specific help for that command."),
                     ("You can save the state of a shell session using 'save' and then recover it at "
                      "a later time using 'load'. When a session is re-loaded, all the variables and "
                      "properties/settings for plots/tables are maintained and the latest state of the "
                      "unit test's results are used. A console history is also maintained with bash-like "
                      "commands (e.g. Ctrl-R for reverse history search, etc.) across sessions. You can "
                      "manipulate its behavior with 'history'.")]
            
            self._fixed_width_info(lines)
        cmd.Cmd.do_help(self, arg)

    def _fixed_width_info(self, lines):
        """Prints the specified string as information with fixed width of 80 chars."""
        for string in lines:
            for line in [string[i:i+80] for i in range(0, len(string), 80)]:
                msg.info(line)
            msg.blank()

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
        #This is where we update the arguments to new dictionary versions when the code
        #is updated with new features etc. so that the older, serialized sessions don't
        #break.
        if (self.active not in self._version_check or 
            self.group not in self._version_check[self.active]):
            if ("version" not in self.args[self.active][self.group] or
                self.args[self.active][self.group]["version"] < self._template_args["version"]):
                for key in self._template_args:
                    if key not in self.args[self.active][self.group]:
                        self.args[self.active][self.group][key] = self._template_args[key]
                self.args[self.active][self.group]["version"] = self._template_args["version"]

            if self.active not in self._version_check:
                self._version_check[self.active] = [self.group]
            else:
                self._version_check[self.active].append(self.group)

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
            return list(self.live.keys())
        else:
            return [c for c in self.live if c.startswith(text)]

    def do_inputs(self, arg):
        usable, filename, append = self._redirect_split(arg)
        result = self._print_dict_files(usable, "inputs")
        if result is not None:
            self._redirect_output(result, filename, append, msg.info)
    def complete_inputs(self, text, line, istart, iend):
        return self._complete_cases(text, line, istart, iend)
    def help_inputs(self):
        lines = [("Prints a list of the available input files to analyze in the specified "
                  "test case. In general, these are files with '.in' in the file name. If "
                  "your unit test XML groups don't rename input files with the '.in' suffix, "
                  "you may have difficulty with the analysis."),
                 ("EXAMPLE: \"inputs standard.1\" lists all the input files available for analysis "
                  "in the 'standard.1' test case. See also 'outputs'.")]
        self._fixed_width_info(lines)

    def do_outputs(self, arg):
        usable, filename, append = self._redirect_split(arg)
        result = self._print_dict_files(usable, "outputs")
        if result is not None:
            self._redirect_output(result, filename, append, msg.info)
    def complete_outputs(self, text, line, istart, iend):
        return self._complete_cases(text, line, istart, iend)
    def help_outputs(self):
        lines = [("Prints a list of the available output files to analyze in the specified "
                  "test case. In general, these are files with '.out' in the file name. If "
                  "your unit test XML groups don't use output files with the '.out' suffix, "
                  "you may have difficulty with the analysis."),
                 ("EXAMPLE: \"outputs standard.1\" lists all the output files available for analysis "
                  "in the 'standard.1' test case. See also: 'inputs'.")]
        self._fixed_width_info(lines)

    def do_timing(self, arg):
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
    def help_timing(self):
        lines = [("Prints the number of milliseconds (ms) that a specific test case took to run. "
                  "This value *only* corresponds to the total time that the actual method being "
                  "unit tested took to run; all other setup/take down methods such as input/output "
                  "reading/writing etc. are *not* included in this value."),
                 ("EXAMPLE: \"timing standard.1\" prints the total run time for the method using "
                  "data from the 'standard.1' test case.")]
        self._fixed_width_info(lines)

    def _complete_tests(self, text, line, istart, iend):
        """Returns a completion list of possible, parsed unit tests that could be loaded or
        otherwise interacted with.
        """
        if text == "":
            return list(self.tests.keys())
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

    def do_headings(self, arg):
        usable, filename, append = self._redirect_split(arg)
        if usable != "":
            headings = usable.split("|")
            self.curargs["headings"] = [h.strip() for h in headings]
        result = "HEADINGS: " + " | ".join(self.curargs["headings"])
        self._redirect_output(result, filename, append, msg.info)        
    def help_headings(self):
        lines = [("Sets the headings to use when printing the table of values for the "
                  "current variable set (i.e. independent and dependent variables). It "
                  "is your responsibility to make sure that the number of headings matches "
                  "the number of variables and that they appear in the same order as they "
                  "are in the session. You can see the order they will appear in the table "
                  "by using the 'vars' command. Add the variables as a '|'-separated list "
                  "of strings; as such bars are not allowed in the heading names."),
                 ("EXAMPLE: \"headings Group Size|Time (ms)|Result\" sets the headings for "
                  "the output table to 'Group Size', 'Time (ms)' and 'Result' respectively."),
                 ("See also: 'table'.")]
        self._fixed_width_info(lines)        

    def do_filter(self, arg):
        """Sets the filter for the test cases to include in the plot/table by name. Only those
        test cases that include this text are included in plots, tables etc."""
        if arg == "list":
            msg.info("TEST CASE FILTERS")
            for f in self.curargs["tfilter"]:
                if f == "*":
                    msg.info("  * (default, matches all)")
                else:
                    msg.info("  " + f)
        elif arg not in self.curargs["tfilter"]:
            self.curargs["tfilter"].append(arg)
            self.do_filter("list")        
    def complete_filter(self, text, line, istart, iend):
        return self._complete_cases(text, line, istart, iend)
    def help_filter(self):
        lines = [("Adds a test case filter to exclude specific test case results from plots "
                  "and tables or failure reports. For example, if I wanted to plot the timing "
                  "for a set of test cases that all have 'standard.fg*' as their name (with * "
                  "a wildcard character), then I could use the following code: "),
                 ("EXAMPLE: \"filter standard.fg*\". If the wildcard character is not present, "
                  "the filter must match one of the test cases *exactly*; otherwise there will be "
                  "no results included in any of the plots/tables. You can remove a filter with "
                  "the 'rmfilter' command.")]
        self._fixed_width_info(lines)

    def do_rmfilter(self, arg):
        """Removes the test case filter that limits which results are included in plots/tables.
        """
        if arg in self.curargs["tfilter"]:
            if arg == "*":
                msg.warn("The default filter cannot be removed.")
            else:
                self.curargs["tfilter"].remove(arg)
                self.do_filter("list")
    def help_rmfilter(self):
        lines = [("Removes the test case filter to exclude specific test case results from plots "
                  "and tables or failure reports. See 'filter'.")]
        self._fixed_width_info(lines)
    def complete_rmfilter(self, text, line, istart, iend):
        return [f for f in self.curargs["tfilter"] if f.startswith(text)]

    def do_threshold(self, arg):
        """Specify a success threshold (percent as float) that output files must attain before 
        they can be included in the plots/tables."""
        try:
            self._set_arg_generic("threshold", arg, float)
        except ValueError:
            err("The specified threshold value is not a valid float percentage. Try 1. for 100%")
    def help_threshold(self):
        lines = [("Except in the case of the 'failures' report, all plots/tables only use those "
                  "test case results that are considered successful. Most unit tests specify an "
                  "output file that needs to be compared to a model output file. When the comparison "
                  "is performed, a percentage value is calculated for the overall match of the output "
                  "with the model output. This comparison is saved in a '.out.compare' file."),
                 ("Success is defined by a comparison percentage that exceeds the 'threshold' value. "
                  "Specify a decimal percentage (e.g. 0.78 for 78%) to change the threshold from its "
                  "default value of 1. (i.e. 100% success)."),
                 ("EXAMPLE: \"threshold 0.66\" sets the minimum success required to 66%. If no "
                  "files were specified for comparison, the test case is considered successful if "
                  "the program exited without any unhandled exceptions (this creates a 'SUCCESS' file "
                  "with the last successful run time in it.")]
        self._fixed_width_info(lines)

    def _complete_filters(self, text):
        return [f + "/" for f in self.curargs["tfilter"] if f.startswith(text)]

    def _complete_vars(self, text):
        #Variables can come from inputs or outputs.
        if text == "":
            return self.allvars+["timing"]
        else:
            vlist = [v + "|" for v in self.allvars if v.startswith(text)]
            if text in "timing":
                if any("timing" in f for f in self.curargs["fits"]):
                    vlist.append("timing|")
                else:
                    vlist.append("timing")
            return vlist

    def _complete_props(self, var, text):
        if var not in self.allprops and not any("timing" in f for f in self.curargs["fits"]):
            return []

        extra = []
        for fitvar in self.curargs["fits"]:
            if var in fitvar and var != "timing":
                extra.append(fitvar.split("|")[1] + ".fit")
            elif var == "timing":
                extra.append("fit")

        if var not in self.allprops:
            return [e for e in extra if e.startswith(text)]

        if text == "":
            return self.allprops[var] + extra
        else:
            return [p for p in (self.allprops[var]+extra) if p.startswith(text)]

    def _complete_fullvar(self, text, line, istart, iend, wfilter=True):
        value = line.split()
        if len(value) == 1:
            value = ""
        else:
            value = value[-1]

        if "/" in value or not wfilter:
            if wfilter:
                filt, rest = value.split("/")
            else:
                #This handles the case where we don't want to complete on filters.
                rest = value
                
            #Determine if we have a bar in the text, if we do then we are completing attributes;
            #otherwise we are completing variables.
            if "|" in rest:
                var, prop = rest.split("|")
                if wfilter:
                    return ["{}/{}|{}".format(filt, var, p) for p in self._complete_props(var, prop)]
                else:
                    return ["{}|{}".format(var, p) for p in self._complete_props(var, prop)]
            else:
                if wfilter:
                    return ["{}/{}".format(filt, v) for v in self._complete_vars(rest)]
                else:
                    return self._complete_vars(rest)
        else:
            return self._complete_filters(value)

    def _validate_var(self, var):
        """Validates the *form* only of the specified variable."""
        if "/" in var:
            filt, rest = var.split("/")
            if "|" in rest:
                varname, prop = rest.split("|")
                return prop != ""
            elif rest == "timing":
                return True
            else:
                return False
        else:
            return True
            
    def do_indep(self, arg):
        """Sets the name and attribute of the independent variable for plotting/tabulating functions."""
        if not self._validate_var(arg):
            msg.err("Variable {} is not a valid file name and property combination.".format(arg))
        else:
            self.curargs["independent"] = arg
    def complete_indep(self, text, line, istart, iend):
        return self._complete_fullvar(text, line, istart, iend, False)
    def help_indep(self):
        lines = [("Sets the variable and property combination for the *independent* variable in plots "
                  "and tables, etc. There can only be a single independent variable for any plot/table. "
                  "The format of the variable specification is filename|property and can be auto-"
                  "completed with <tab>. The properties are extracted from the file when it is parsed. "
                  "Depending on the variable property you choose, you either plot values extracted from "
                  "*all* test cases across the whole unit test, OR only array-type values from multiple "
                  "files within the *same* test case (see the note about exceptions below)."),
                 ("Possible properties are:\n"
                  "  - 'width': the number of values per row for 2D square arrays.\n"
                  "  - 'depth': the number of rows of data in the file.\n"
                  "  - 'shape': the (depth, width) of the data as a tuple.\n"
                  "  - 'value': the scalar value for a file that has only a single value in it.\n"
                  "  - 'rowvals': a list of the values for each row. Requires aggregation function.\n"
                  "  - 'colvals': a list of the values for each column. Requires aggregation function"),
                 ("For the 'rowvals' and 'colvals' properties, the unit test needs to be parsed using "
                  "the 'fullparse' command. These properties allow the contents of an input/output file "
                  "for a single unit test to be plotted against some other data from a different file."
                  "Choose an aggregation function (using the 'postfix' command) like numpy.mean or "
                  "numpy.sum to turn each list of values into a single number. EXCEPTION: As long as the "
                  "dimensionality of the data allows the property to be reduced to a single number "
                  "when the aggregation function is applied, you could also plot 'rowvals' over all "
                  "the unit tests. For example, if all the 'rowvals' for a given filter are just 1D "
                  "arrays and you take a sum, then the values are still aggregatable."),
                 ("EXAMPLE: \"indep group.in|depth\" chooses the number of rows in each 'group.in' "
                  "file across all the unit test's cases to be the x-value for any plots or tables."),
                 ("EXAMPLE: \"indep generators.in|rowvals\" selects the set of row values from the "
                  "'generators.in' file to be plotted. To work, you would also need to *first* "
                  "set the 'postfix' function to a valid aggregation function."),
                 ("NOTE: use tab completion to make sure that the variable names and properties you "
                  "select are valid. The auto-completion is specific to the file you choose, so it "
                  "only reflects what is possible."),
                 ("See Also: 'dep'.")]
        self._fixed_width_info(lines)

    def do_dep(self, args):
        """Adds the name and attribute of a dependent variable to the list 
        for plotting/tabulating functions."""
        vals = args.split()
        twin = None
        if len(vals) > 1:
            if len(vals) == 2:
                var, plot = vals
            elif len(vals) == 3:
                var, plot, twin = vals
        else:
            var = vals[0]
            plot = None

        if not self._validate_var(var):
            msg.err("Variable {} is not a valid file name and property combination.".format(var))
        else:
            if var in self.curargs["dependents"]:
                self.curargs["dependents"].remove(var)
            self.curargs["dependents"].append(var)
            if plot is not None:
                self.curargs["plottypes"][var] = plot
            if twin is not None:
                self.curargs["twinplots"][var] = twin
    def complete_dep(self, text, line, istart, iend):
        els = line.split()
        if len(els) == 1 or (len(els) == 2 and line[-1] != " "):
            if len(els) == 1:
                part = ""
            else:
                part = els[-1]
            return self._complete_fullvar(part, line, istart, iend)
        elif (len(els) == 2 and line[-1] == " ") or (len(els) == 3 and line[-1] != " "):
            keys = ["line", "scatter"]
            if line[-1] == " ":
                return keys
            else:
                return [c for c in keys if c.startswith(text)]
        else:
            twins = ["twinx", "twiny"]
            if line[-1] == " ":
                return twins
            else:
                return [c for c in twins if c.startswith(text)]
    def help_dep(self):
        lines = [("Adds one variable as a *dependent* variable for plotting or tabulating. "
                  "A description of the format of the variables and properties can be found by typing "
                  "\"help indep\" in the shell. The dependent variables follow similar "
                  "conventions and the same limitations apply regarding the use of 'rowvals' and "
                  "'colvals' properties. One addition for the dependent variables, however, is that "
                  "each variable can have its own filter specification. The filter must have been "
                  "added using the 'filter' command. This allows dependent variables from different "
                  "sets of test cases to be plotted on the same plot."),
                 ("EXAMPLE: \"dep */concs.in|width line\" adds 'concs.in|depth' as a "
                  "dependent variable whose values will be extracted across all test cases"
                  " that match the filter specified. The '*' filter matches *all* test cases. The "
                  "line specification sets the plot-type to line instead of scatter."),
                 ("NOTE: use tab completion to make sure that the variable names and properties you "
                  "select are valid. The auto-completion is specific to the file you choose, so it "
                  "only reflects what is possible.")]
        self._fixed_width_info(lines)

    def do_rmdep(self, args):
        """Removes dependent variables currently set for plotting/tabulating etc."""
        for arg in args.split():
            if arg in self.curargs["dependents"]:
                self.curargs["dependents"].remove(arg)
            if arg in self.curargs["plottypes"]:
                del self.curargs["plottypes"][arg]
            if arg in self.curargs["twinplots"]:
                del self.curargs["twinplots"][arg]
            if arg in self.curargs["colors"]:
                del self.curargs["colors"][arg]
            if arg in self.curargs["labels"]:
                del self.curargs["labels"][arg]
            if arg in self.curargs["markers"]:
                del self.curargs["markers"][arg]
            if arg in self.curargs["lines"]:
                del self.curargs["lines"][arg]
                
    def complete_rmdep(self, text, line, istart, iend):
        els = line.split()
        if len(els) == 1:
            return self.curargs["dependents"]
        else:
            if text != els[1]:
                return [text + v.replace(els[1], "") for v in self.curargs["dependents"] 
                        if v.startswith(els[1])]
            else:
                return [v for v in self.curargs["dependents"] if v.startswith(els[1])]
    def help_rmdep(self):
        lines = [("Removes the specified variable from the list of dependent variables for "
                  "the tables or plots generated using the shell.")]
        self._fixed_width_info(lines)

    def do_vars(self, arg):
        usable, filename, append = self._redirect_split(arg)
        result = []
        result.append("INDEPENDENT: '{}'".format(self.curargs["independent"]))
        result.append("DEPENDENTS")
        for i in range(len(self.curargs["dependents"])):
            var = self.curargs["dependents"][i]
            plottype = "scatter"
            twinplot = "" if var not in self.curargs["twinplots"] else "-{}".format(self.curargs["twinplots"][var])
            if var in self.curargs["plottypes"]:
                plottype = self.curargs["plottypes"][var]
            result.append("  {}. '{}' ({}{})".format(i, var, plottype, twinplot))
        self._redirect_output('\n'.join(result), filename, append, msg.info)
    def help_vars(self):
        lines = [("Prints the current status of the independent and dependent variables for the "
                  "current analysis group.")]
        self._fixed_width_info(lines)

    def _print_map_dict(self, argkey, filename, append):
        """Prints a dictionary that has variable => value mappings."""
        result = []
        skeys = list(sorted(self.curargs[argkey].keys()))
        for key in skeys:
            result.append("'{}' => {}".format(key, self.curargs[argkey][key]))
        self._redirect_output('\n'.join(result), filename, append, msg.info)

    def do_postfix(self, arg):
        """Sets the function to apply to the values of a specific variable before plotting
        or tabulating values.
        """
        usable, filename, append = self._redirect_split(arg)
        sargs = usable.split()
        if len(sargs) == 1 and sargs[0] == "list":
            self._print_map_dict("functions", filename, append)
        elif len(sargs) >= 2:
            defvars = self._postfix_varlist("postfix " + arg)
            for var in defvars.values():
                if not self._validate_var(var):
                    msg.err("Variable '{}' is not a valid variable|property combination.")
                    return

            fxn = arg.split()[-1]
            if ":" not in fxn:
                msg.err("{} is not a valid postfix function expression.")
                self.help_postfix()
                return

            modvar = fxn.split(":")[0]
            if modvar not in defvars:
                msg.err("Invalid postfix function: variable '{}' not defined.".format(modvar))
                return
            
            defvars["lambda"] = fxn
            self.curargs["functions"][defvars[modvar]] = defvars
            #Give the user some feedback so that they know it was successful.
            self.do_postfix("list")
    def _postfix_varlist(self, line):
        """Returns a dictionary of the global variable names (keys) and their local 
        name for the lambda function.
        """
        els = line.split()
        result = {}
        if len(els) >= 2:
            defvars = [v for v in els if "=" in v]
            varlist = []
            for dvar in defvars:
                gvar, lvar = dvar.split("=")
                if lvar != "":
                    result[lvar] = gvar
        return result
    
    def complete_postfix(self, text, line, istart, iend):
        els = line.split()
        if len(els) >= 2:
            if ":" in els[-1]:
                #We are completing the lambda function. Find a list of the variables that
                #they have defined so far.
                defvars = self._postfix_varlist(line)
                if len(defvars) == 0:
                    return []
                else:
                    #We want to complete for the numpy functions and the variable names
                    #that we know they have defined.
                    comptext = re.split(r"[()\]\-\^\\/[+*:]+", els[-1])[-1]
                    left = els[-1][:els[-1].rfind(comptext)]
                    if "." in comptext:
                        libs = comptext.split(".")
                        from importlib import import_module
                        partial = libs.pop()
                        prefix = '.'.join(libs)
                        module = import_module(prefix)
                        alldir = [a for a in dir(module) if a[0] != "_"]

                        if partial == "":
                            result = alldir
                        else:
                            result = [a for a in alldir if a.startswith(partial)]
                        return [left + prefix + '.' + a for a in result]
                    else:
                        varlist = list(defvars.keys()) + ["numpy.", "math."]
                        return [left + a for a in varlist if a.startswith(comptext)]
            else:
                #We are compiling the list of variables to be used in the lambda.
                defvars = self._postfix_varlist(line)
                if "=" in els[-1] and line[-1] != " ":
                    #They are choosing a variable name, return the generic variable name.
                    gvar, lvar = els[-1].split("=")
                    if len(lvar) == 0:
                        varlist = [els[-1] + chr(97+len(defvars))]
                    else:
                        varlist = [els[-1] + " "]
                elif "=" in els[-1] and line[-1] == " ":
                    #Auto-complete on the variable name.
                    defvars = self._postfix_varlist(line)
                    varlist = self._complete_fullvar("", "postfix ", 0, 0)
                    varlist.extend([v + ':' for v in defvars])
                else:
                    varlist = self._complete_fullvar(els[-1], "postfix {}".format(els[-1]).format(els[-1]), 0, 0)
                    if len(varlist) == 1 and varlist[0] == els[-1]:
                        varlist = [els[-1] + '=']
                    varlist.extend([v + ':' for v in defvars if v.startswith(els[-1])])
                    if "list".startswith(els[-1]):
                        varlist.append("list")
                return varlist
        else:
            varlist = self._complete_fullvar(text, line, istart, iend)
            if text in "list" and "|" not in line:
                varlist.append("list")
            return varlist
                
    def help_postfix(self):
        lines = [("Sets the postfix function for a specific variable (either dependent or "
                  "independent). The postfix function is applied to the variable's value "
                  "before plotting. Any numpy function can be selected for the postfix."),
                 ("EXAMPLE \"postfix group.in|rowvals numpy.mean\" sets the postfix function "
                  "to take the mean value of each row in the 'group.in' file for the plotting.")]
        self._fixed_width_info(lines)

    def do_rmpostfix(self, arg):
        """Removes a postfix function from a variable. See 'postfix'."""
        altered = False
        if arg in self.curargs["functions"]:
            del self.curargs["functions"][arg]
            altered = True
        elif arg == "*":
            for varname in list(self.curargs["functions"].keys()):
                del self.curargs["functions"][varname]
            altered = True
        if altered:
            self.do_postfix("list")
    def complete_rmpostfix(self, text, lines, istart, iend):
        return [p for p in list(self.curargs["functions"].keys()) if p.startswith(text)]

    def do_fit(self, arg):
        usable, filename, append = self._redirect_split(arg)
        sargs = usable.split()
        if len(sargs) == 1 and sargs[0] == "list":
            self._print_map_dict("fits", filename, append)
        elif len(sargs) == 2:
            var, fxn = sargs
            if not self._validate_var(var):
                msg.err("Variable '{}' is not a valid variable|property combination.")
            else:
                self.curargs["fits"][var] = fxn
                #Give the user some feedback so that they know it was successful.
                self.do_fit("list")
                #Also add the fit to the dependent variables automatically, since
                #it is extremely likely that they will want to plot it.
                if "timing" in var:
                    self.do_dep(var + "|fit")
                else:
                    self.do_dep(var + ".fit")
    def complete_fit(self, text, line, istart, iend):
        els = line.split()
        if len(els) == 1 or (len(els) == 2 and line[-1] != " "):
            varlist = self._complete_fullvar(text, line, istart, iend)
            if text in "list" and "|" not in line:
                varlist.append("list")
            return varlist
        elif line[-1] == " ":
            return ["exp", "lin"]
        else:
            return [f for f in ["exp", "lin"] if f.startswith(els[-1])]
    def help_fit(self):
        lines = [("Sets a function to fit a dependent variable's data to relative to its "
                  "independent variable. The values are fit using one of the common fitting "
                  "function types and can then be plotted on the same curve as the original "
                  "data or tabulated. To add the fit as a variable, use the 'fit' property "
                  "of the variable."),
                 ("Possible values are 'exp', 'lin'."),
                 ("EXAMPLE \"fit group.in|rowvals exp\" allows the aggregated 'rowvals' for "
                  "the 'group.in' file to be fitted by an exponential curve. This adds the "
                  "property 'fit' to 'group.in' as a property.")]
        self._fixed_width_info(lines)

    def do_rmfit(self, arg):
        """Removes a fit function from a variable. See 'fit'."""
        if arg in self.curargs["fits"]:
            del self.curargs["fits"][arg]        
        #We also need to remove the variable entry if it exists.
        if "timing" in arg:
            fitvar = "{}|fit".format(arg)
        else:
            fitvar = "{}.fit".format(arg)
        if fitvar in self.curargs["dependents"]:
            self.curargs["dependents"].remove(fitvar)
    def complete_rmfit(self, text, lines, istart, iend):
        return [p for p in list(self.curargs["fits"].keys()) if p.startswith(text)]

    def _complete_deps(self, els, line, addlist=True):
        if len(els) == 1:
            part = ""
        else:
            part = els[1]

        varlist = [v for v in self.curargs["dependents"] if v.startswith(part)]
        if part in "list" and "|" not in line and addlist:
            varlist.append("list")
        if part == "":
            varlist.append("*")
        return varlist

    def do_label(self, arg):
        usable, filename, append = self._redirect_split(arg)
        sargs = usable.split()
        if len(sargs) == 1 and sargs[0] == "list":
            self._print_map_dict("labels", filename, append)
        else:
            var, label = sargs[0], ' '.join(sargs[1:len(sargs)])
            options = ["plot", "x", "y", "z", "x-twin", "y-twin"]
            if var not in options and not self._validate_var(var):
                msg.err("Variable '{}' is not a valid variable|property combination.")
            else:
                if var == "*":
                    for depvar in self.curargs["dependents"]:
                        self.curargs["labels"][depvar] = label
                else:
                    self.curargs["labels"][var] = label
                #Give the user some feedback so that they know it was successful.
                self.do_label("list")
    def complete_label(self, text, line, istart, iend):
        els = line.split()
        if len(els) == 1 or (len(els) == 2 and line[-1] != " "):
            varlist = self._complete_deps(els, line)
            varlist.extend(["plot", "x", "y", "z", "x-twin", "y-twin"])
            if len(els) == 1:
                return varlist
            else:
                return [p for p in varlist if p.startswith(els[-1])]
        else:
            return []
    def help_label(self):
        lines = [("Sets the text for the legend of a specific variable in the plot. "
                  "When setting the label, it will show up exactly as typed, you don't "
                  "need to use quotes around strings etc. To prevent a label from being "
                  "applied to a specific variable, set its value to [None]."),
                 ("EXAMPLE: \"label group.in|width Group Size\" "
                  "sets the legend label to 'Group Size'.")]
        self._fixed_width_info(lines)

    def do_rmlabel(self, arg):
        """Removes a label mapping from a variable. See 'label'."""
        if arg in self.curargs["labels"]:
            del self.curargs["labels"][arg]        
    def complete_rmlabel(self, text, lines, istart, iend):
        return [p for p in list(self.curargs["labels"].keys()) if p.startswith(text)]

    def do_color(self, arg):
        usable, filename, append = self._redirect_split(arg)
        sargs = usable.split()
        if len(sargs) == 1 and sargs[0] == "list":
            self._print_map_dict("colors", filename, append)
        elif len(sargs) == 2:
            var, col = sargs
            if not self._validate_var(var):
                msg.err("Variable '{}' is not a valid variable|property combination.")
            else:
                if var == "*":
                    for depvar in self.curargs["dependents"]:
                        self.curargs["colors"][depvar] = col
                else:
                    self.curargs["colors"][var] = col
                #Give the user some feedback so that they know it was successful.
                self.do_color("list")
    def complete_color(self, text, line, istart, iend):
        els = line.split()
        if len(els) == 1 or (len(els) == 2 and line[-1] != " "):
            return self._complete_deps(els, line)
        else:
            keys = list(self._possible_cols.keys()) + ["0.", "#"]
            if line[-1] == " ":
                return keys
            else:
                return [c for c in keys if c.startswith(text)]
    def help_color(self):
        lines = [("Sets the color of a specific variable in a plot."),
                 ("EXAMPLE: \"color group.in|depth blue\" sets the plot color for the variable "
                  "to be blue for any plots in the current analysis group.")]
        self._fixed_width_info(lines)

    def do_rmcolor(self, arg):
        """Removes a color mapping from a variable. See 'color'."""
        if arg in self.curargs["colors"]:
            del self.curargs["colors"][arg]        
    def complete_rmcolor(self, text, lines, istart, iend):
        return [p for p in list(self.curargs["colors"].keys()) if p.startswith(text)]

    def _get_matplot_dict(self, option, prop, defdict):
        """Returns a copy of the settings dictionary for the specified option in 
        curargs with update values where the value is replaced by the key from 
        the relevant default dictionary.

        :arg option: the key in self.curargs to update.
        :arg defdict: the default dictionary whose keys should be used when values match.
        """
        cargs = self.curargs[option]
        result = cargs.copy()
        for varname in cargs:
            if prop in cargs[varname]:
                name = cargs[varname][prop]
                for key, val in list(defdict.items()):
                    if val == name:
                        cargs[varname][prop] = key
                        break

        return result

    def _plot_generic(self, filename=None):
        """Plots the current state of the shell, saving the value to the specified file
        if specified.
        """
        #Since the filename is being passed directly from the argument, check its validity.
        if filename == "":
            filename = None

        if "x" not in self.curargs["labels"]:
            #Set a default x-label since we know what variable is being plotted.
            self.curargs["labels"]["x"] = "Value of '{}' (unknown units)".format(self.curargs["independent"])
        args = self.curargs
        a = self.tests[self.active]
        self._make_fits()
        
        #Before we can pass the markers in, we need to translate from keys to values so
        #that matplotlib understands.
        markdict = self._get_matplot_dict("markers", "marker", self._possible_markers)
        linedict = self._get_matplot_dict("lines", "style", self._possible_linestyles)

        #Set the remaining arguments to have the right keyword name.
        args["savefile"] = filename
        args["markers"] = markdict
        args["lines"] = linedict
        
        a.plot(**args)

    def do_logplot(self, arg):
        """Plots the current state of the shell's independent vs. dependent variables on the
        same set of axes with the y-scale set to logarithmic. Give filename to save to as
        argument or leave blank to show.
        """        
        usable, filename, append = self._redirect_split(arg)
        self.curargs["xscale"] = None
        self.curargs["yscale"] = "log"
        self._plot_generic(filename)
    def help_logplot(self):
        lines = [("Plots the behavior of the dependent variables as functions of the independent "
                  "variable for the current analysis group. The y-scale is set to logarithmic. To "
                  "save the plot to a file, redirect the output."),
                 ("You can control the appearance of the plot using the following commands: "
                  "'xlabel', 'ylabel'. Consider using separate analysis groups for each plot "
                  "you need to create so that you don't waste time overwriting property values. "
                  "You can save a sessions settings using 'save' if you want to return and "
                  "redo the plots later or make adjustments (highly recommended)."),
                 ("EXAMPLE \"logplot > plot.pdf\" creates a log plot of the data and saves it "
                  "to the specified PDF file."),
                 ("See also: 'loglogplot', 'plot'.")]
        self._fixed_width_info(lines)
        
    def do_semilogxplot(self, arg):
        """Plots the current state of the shell's independent vs. dependent variables on the
        same set of axes with the y-scale set to logarithmic. Give filename to save to as
        argument or leave blank to show.
        """        
        usable, filename, append = self._redirect_split(arg)
        self.curargs["xscale"] = "log"
        self.curargs["yscale"] = None
        self._plot_generic(filename)
    def help_semilogxplot(self):
        lines = [("Plots the behavior of the dependent variables as functions of the independent "
                  "variable for the current analysis group. The x-scale is set to logarithmic. To "
                  "save the plot to a file, redirect the output."),
                 ("You can control the appearance of the plot using the following commands: "
                  "'xlabel', 'ylabel'. Consider using separate analysis groups for each plot "
                  "you need to create so that you don't waste time overwriting property values. "
                  "You can save a sessions settings using 'save' if you want to return and "
                  "redo the plots later or make adjustments (highly recommended)."),
                 ("EXAMPLE \"semilogxplot > plot.pdf\" creates a log plot of the data and saves it "
                  "to the specified PDF file."),
                 ("See also: 'loglogplot', 'plot', 'logplot'.")]
        self._fixed_width_info(lines)

    def do_loglogplot(self, arg):
        """Plots the current state of the shell's independent vs. dependent variables on the
        same set of axes with the x- and y-scale set to logarithmic. Give filename to save to as
        argument or leave blank to show.
        """        
        usable, filename, append = self._redirect_split(arg)
        self.curargs["xscale"] = "log"
        self.curargs["yscale"] = "log"
        self._plot_generic(filename)
    def help_loglogplot(self):
        lines = [("Plots the behavior of the dependent variables as functions of the independent "
                  "variable for the current analysis group. The y-scale is set to logarithmic and "
                  "the x-scale is also set to logarithmic. To save the plot to a file, redirect "
                  "the output."),
                 ("You can control the appearance of the plot using the following commands: "
                  "'xlabel', 'ylabel'. Consider using separate analysis groups for each plot "
                  "you need to create so that you don't waste time overwriting property values. "
                  "You can save a sessions settings using 'save' if you want to return and "
                  "redo the plots later or make adjustments (highly recommended)."),
                 ("EXAMPLE \"loglogplot > plot.pdf\" creates a log-log plot of the data and saves it "
                  "to the specified PDF file."),
                 ("See also: 'logplot', 'plot', 'semilogx'.")]
        self._fixed_width_info(lines)

    def do_plot(self, arg):
        """Plots the current state of the shell's independent vs. dependent variables on the
        same set of axes. Give filename to save to as argument or leave blank to show.
        """
        usable, filename, append = self._redirect_split(arg)
        self.curargs["xscale"] = None
        self.curargs["yscale"] = None
        self._plot_generic(filename)
    def help_plot(self):
        lines = [("Plots the behavior of the dependent variables as functions of the independent "
                  "variable for the current analysis group. To save the plot to a file, redirect "
                  "the output."),
                 ("You can control the appearance of the plot using the following commands: "
                  "'xlabel', 'ylabel'. Consider using separate analysis groups for each plot "
                  "you need to create so that you don't waste time overwriting property values. "
                  "You can save a sessions settings using 'save' if you want to return and "
                  "redo the plots later or make adjustments (highly recommended)."),
                 ("EXAMPLE \"plot > plot.pdf\" creates a plot of the data and saves it "
                  "to the specified PDF file."),
                 ("See also: 'logplot', 'loglogplot', 'semilogx'.")]
        self._fixed_width_info(lines)

    def _set_def_prompt(self):
        """Sets the default prompt to match the currently active unit test."""
        if len(self.active) > 15:
            ids = self.active.split(".")
            if len(ids) > 2:
                module, executable, compiler = ids
            else:
                module, executable = ids
                compiler = "g"
            self.prompt = "({}*.{}*.{}:{})".format(module[0:6], executable[0:6], compiler, self.group)
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
    def help_set(self):
        lines = [("Sets the currently active unit test whose test cases are being analyzed or "
                  "plotted. You can add multiple unit tests to the session by using the 'parse' "
                  "or 'fullparse' commands. When any of those commands is used, this command "
                  "is also called automatically to set the session to use the newly parsed "
                  "unit test as the active one. This command toggles between them."),
                 ("EXAMPLE \"set classes.polya\" sets the unit tests for the executable 'polya' in "
                  "the 'classes' module to be the active unit test."),
                 ("NOTE: you can see which unit test is currently active by looking at the prompt.")]
        self._fixed_width_info(lines)

    def do_save(self, arg):
        """Saves the session variables to a file so that the same analysis can be continued later."""
        #We need to save the full paths to the staging directories for the tests that are loaded
        #so far; then when the session is restored, we can reparse the results.
        from os import path
        import json
        fullpath = path.expanduser(arg)

        data = {
            "tests": [path.abspath(a.stagedir) for a in list(self.tests.values())],
            "args": self.args
        }

        with open(fullpath, 'w') as f:
            json.dump(data, f)
        msg.okay("Saved current session to {}".format(fullpath))
    def complete_save(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)
    def help_save(self):
        lines = [("Saves the current test analysis session to disk so it can be resumed later. "
                  "The current state of all variables, settings and a list of the unit tests that "
                  "were parsed into the session are also saved. When the session is re-loaded "
                  "later using the 'load' command, the *latest* state of all the unit tests is "
                  "re-parsed from the disk; i.e. no results are cached, just pointers to the results, "
                  "which are re-parsed on load."),
                 ("EXAMPLE: \"save analysis.json\" saves the session to the specified file. The "
                  "session is always serialized in JSON."),
                 ("See also: 'load'")]
        self._fixed_width_info(lines)

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
    def help_load(self):
        lines = [("Loads a previously saved shell session from disk. For details, see 'save'."),
                 ("EXAMPLE \"load analysis.json\" loads the saved session.")]
        self._fixed_width_info(lines)
        
    def do_parse(self, arg, fullparse=False):
        """Parse the test results from the specified directory and load them under the name
        of 'module.executable ' that they were created with. E.g. parse classes.polya/
        """
        from os import path
        fullpath = path.abspath(path.expanduser(arg))
        if path.isdir(fullpath):
            if fullpath[-1] == "/":
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

    def help_parse(self):
        lines = [("Parses the test results in the sub-folders of the specified unit test "
                  "staging directory. A staging directory is the one that contains all the "
                  "*.f90 files, the Makefile.* files, etc. as well as a folder called tests/ "
                  "that holds the results of the individual test case runs. You can call this "
                  "command multiple times to reload a specific unit test. However, if the "
                  "unit test you are trying to re-parse is the active one, it is probably "
                  "quicker to use the 'reparse' command, which reparses the active unit test."),
                 ("EXAMPLE \"parse classes.polya/\" parses all the test case runs in the "
                  "unit test for the 'polya' executable that belongs to the 'classes' module.")]
        self._fixed_width_info(lines)

    def do_fullparse(self, arg):
        """Parses a set of unit tests *and* the contents of their input and output files
        so that test results can be plotted/tabulated against the data in the input/output
        files. Takes much longer than the regular 'parse' command if the data files are large
        since all the data has to be loaded to memory.
        """
        self.do_parse(arg, True)
    def complete_fullparse(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)
    def help_fullparse(self):
        lines = [("As for the 'parse' command, but also parses the *contents* of the input "
                  "and output files in each test case folder. This allows the shell to plot "
                  "data from a specific test case easily, which can be useful for generating "
                  "plots for papers if the unit testing framework was used to debug and run "
                  "the code. See also 'parse'."),
                 ("EXAMPLE \"fullparse classes.polya/\" parses all the test case runs in the "
                  "unit test for the 'polya' executable that belongs to the 'classes' module. "
                  "Also parses the *contents* of the input and output files.")]
        self._fixed_width_info(lines)

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
    def help_reparse(self):
        lines = [("Performs either 'parse' or 'fullparse' for the *active* unit test "
                  "in the shell."),
                 ("EXAMPLE \"(classes.polya:default) reparse\" reparses the test cases for "
                  "the 'polya' executable in the 'classes' module. That is the active unit "
                  "test as can be seen by the prompt beside the command."),
                 ("EXAMPLE \"(classes.polya:default) reparse full\" same as the above example, "
                  "but uses the 'fullparse' command instead."),
                 ("See also: 'parse', 'fullparse'")]
        self._fixed_width_info(lines)

    def do_group(self, arg):
        """Creates a new analysis group with unique settings for plotting/tabulating etc.
        or switches the active group to the specified name.
        """
        from copy import deepcopy
        vals = arg.split()
        if len(vals) == 0 or vals[0] not in ["list", "duplicate", "add", "switch", "remove"]:
            self.help_group()
            return

        if vals[0] == "add":
            if vals[1] not in self.args[self.active]:
                self.args[self.active][vals[1]] = self._template_args.copy()
                msg.okay("Created analysis group '{}'.".format(vals[1]))
                self.do_group("switch {}".format(vals[1]))
            else:
                msg.info("Group '{}' already exists. Switching to it.".format(vals[1]))
                self.do_group("switch {}".format(vals[1]))
        elif vals[0] == "switch":
            if vals[1] in self.args[self.active]:
                self.group = vals[1]
                self._set_def_prompt()
            else:
                msg.warn("The group '{}' does not exist.".format(vals[1]))
        elif vals[0] == "duplicate" and len(vals) == 3 and vals[1] in self.args[self.active]:
            self.args[self.active][vals[2]] = deepcopy(self.args[self.active][vals[1]])
            msg.okay("Duplicated analysis group '{}' into '{}'.".format(vals[1], vals[2]))
            self.do_group("switch {}".format(vals[2]))
        elif vals[0] == "list":
            for key in self.args[self.active]:
                msg.info(key)
        elif vals[0] == "remove" and vals[1] in self.args[self.active]:
            del self.args[self.active][vals[1]]
            self.do_group("list")
            
    def complete_group(self, text, line, istart, iend):
        els = line.split()
        options = ["list", "duplicate", "add", "switch", "remove"]
        if len(els) == 1:
            return options
        elif len(els) == 2 and line[-1] != " ":
            return [o for o in options if o.startswith(text)]
        elif len(els) == 3 or (len(els) == 2 and line[-1] == " "):
            if text == "":
                return list(self.args[self.active].keys())
            else:
                return [k for k in self.args[self.active] if k.startswith(text)]
    def help_group(self):
        lines = [("Adds/duplicates/switches the active analysis group. When a new group is "
                  "created, it gets a new set of default variables and properties that "
                  "can be set independently of any other analysis groups in the session. "
                  "Groups are most useful for defining multiple plots that can be made from "
                  "the test cases of a single unit test."),
                 ("EXAMPLE: \"group add newplot\" creates a new analysis group called 'newplot' that "
                  "has its own set of plot properties. All the variables (independent and "
                  "dependent) have to be set since a new group is completely blank."),
                 ("EXAMPLE: \"group duplicate default copy\" creates an identical copy of the "
                  "variables in the 'default' group under the new group name 'copy'.")]
        self._fixed_width_info(lines)

    def _make_fits(self):
        """Generates the data fits for any variables set for fitting in the shell."""
        a = self.tests[self.active]
        args = self.curargs
        #We need to generate a fit for the data if there are any fits specified.
        if len(args["fits"]) > 0:
            for fit in list(args["fits"].keys()):
                a.fit(args["independent"], fit, args["fits"][fit], args["threshold"], args["functions"])

    def do_table(self, arg):
        """Prints the set of values for the independent vs. dependent variables in the
        active unit test and analysis group as a table.
        """
        usable, filename, append = self._redirect_split(arg)
        a = self.tests[self.active]
        args = self.curargs
        self._make_fits()
        result = a.table(args["independent"], args["dependents"], args["threshold"],
                         args["headings"], args["functions"])
        if result is not None:
            self._redirect_output(result, filename, append, msg.info)
    def complete_table(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)
    def help_table(self):
        lines = [("Generates a table of data points that would be plotted were any of the "
                  "plot commands to be called. The first column is the independent variable, "
                  "all other columns belong to the dependent variables in the order in which "
                  "they were added to the analysis group. The output can be redirected to a file "
                  "for use by other applications."),
                 ("EXAMPLE: \"table > vartable.dat\" saves a table of the data to the specified "
                  "file. Headings for the table can be set using the 'headings' command.")]
        self._fixed_width_info(lines)

    def do_failures(self, arg):
        """Prints a list of test cases that failed for the current unit test and analysis
        group settings. To only check failure on specific output files, set the list of
        files to check as arguments.
        """
        usable, filename, append = self._redirect_split(arg)
        a = self.tests[self.active]
        args = self.curargs
        splitargs = usable.split()
        if len(splitargs) > 0:
            tfilter = splitargs[0]
        else:
            tfilter = "*"

        outfiles = None
        if len(splitargs) > 1:
            outfiles = splitargs[1:len(splitargs)]

        result = a.failures(outfiles, args["threshold"], tfilter)
        self._redirect_output(result, filename, append, msg.info)
    def complete_failures(self, text, line, istart, iend):
        els = line.split()
        if len(els) == 1 or (len(els) == 2 and line[-1] != " "):
            if len(els) == 1:
                part = ""
            else:
                part = els[1]
            return [f for f in self.curargs["tfilter"] if f.startswith(part)]
        else:
            #We only want to return those output files in the current unit test who *also*
            #have .compare files to check status on.
            if text == "":
                return [f for f in self.tests[self.active].comparable]
            else:
                return [f for f in self.tests[self.active].comparable if f.startswith(text)]
    def help_failures(self):
        lines = [("Searches the current unit test for test cases that failed, either because "
                  "they don't meet the threshold criteria, or because the executable generated "
                  "an uncaught exception. The failures are listed as a table and may be saved "
                  "to an external file using a redirect. To examine a failure because of low "
                  "percent accuracy, use the 'examine' command."),
                 ("You can filter the failures that are searched by specifying output files "
                  "that have accompanying '.out.compare' comparison reports. In that case, only "
                  "the comparisons with sub-threshold values from the specified files will be "
                  "reported."),
                 ("You can also filter by specifying a test case filter that only returns test "
                  "cases who match the string set using the 'filter' command. Any time the "
                  "failures are filtered, the 2nd argument is *always* the test case filter, all "
                  "subsequent arguments should be output file names to filter failures on."),
                 ("EXAMPLE: \"failures * polya.out\" searches for test cases that failed because "
                  "of low accuracy in the 'polya.out' file. This will include tests that never "
                  "finished because of exceptions.")]
        self._fixed_width_info(lines)

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
    def help_examine(self):
        lines = [("Opens a '.out.compare' file from a specific test case inside a "
                  "unit test to see which parts of the test comparison caused a failure. "
                  "The file is opened with the text program specified by the $EDITOR "
                  "environment variable. If the variable is not set, the command will "
                  "display a warning and exit."),
                 ("EXAMPLE: \"examine standard.1 polya.out\" opens the 'polya.out.compare' "
                  "file for examination inside the $EDITOR program. Use the tab-completion "
                  "functionality when selecting the test case and file to make sure that "
                  "it exists.")]
        self._fixed_width_info(lines)

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
        return path.join(settings.cache_directory, "history")

    def preloop(self):
        cmd.Cmd.preloop(self)
        #We need to restore the console history if it exists.
        import readline
        readline.set_completer_delims(' \t\n')
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

    def cmdloop(self):
#        try:
        cmd.Cmd.cmdloop(self)
        # except Exception as exsimple:
        #     msg.err(exsimple.message)
        #     self._store_lasterr()
        #     if self._errcount < self._maxerr:
        #         self._errcount += 1
        #         msg.err("The shell has caught {} unhandled exceptions so far.\n".format(self._errcount) + 
        #                 "When that value reaches {}, the shell will save a ".format(self._maxerr) + 
        #                 "recovery file and exit.")
        #         self.postloop()
        #         self.cmdloop()
        #     else:
        #         self.do_save("#fortpy.shell#")
        #         msg.err("Something unexpected happened. The shell has died. Your session "
        #                 "has been saved as '#fortpy.shell#' in the current directory.")

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

            if command in self._var_cmds:
                #We need to make sure that we have variables set.
                if self.curargs["independent"] is None or len(self.curargs["dependents"]) == 0:
                    msg.err("This command requires an independent variable to be set and "
                            "at least one dependent variable.\n See 'dep' and 'indep' commands.")
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

    def do_font(self, arg):
        opts = arg.split()
        if len(opts) == 0:
            self.help_font()
            return
        if opts[0] == "family":
            value = ' '.join(opts[1:len(opts)])
            if value == "rm":
                del self.curargs["fonts"]["family"]
            else:
                self.curargs["fonts"]["family"] = value
            self.do_font("list")
        elif opts[0] in ["xticks", "xlabel", "yticks", "ylabel", "legend", "title",
                         "x-twin-ticks", "x-twin-label", "y-twin-ticks", "y-twin-label"]:
            if opts[1] == "rm":
                del self.curargs["fonts"][opts[0]]
            else:
                props = {}
                for p in opts[1:len(opts)]:
                    prop, val = p.split("=")
                    props[prop] = val
                self.curargs["fonts"][opts[0]] = props
            self.do_font("list")
        elif opts[0] == "list":
            for key in self.curargs["fonts"]:
                if key == "family":
                    msg.info("FONT FAMILY: {}".format(self.curargs["fonts"]["family"]))
                else:
                    msg.info("{} FONTS:".format(key.upper()))
                    for prop in self.curargs["fonts"][key]:
                        msg.info("  {} => {}".format(prop, self.curargs["fonts"][key][prop]))
    def help_font(self):
        lines = [("Sets the font settings on plots for the current analysis group. Possible "
                  "font options are:\n "
                  "- family: font family for all fonts in the plot.\n"
                  "- xticks, yticks: font settings for the digits marking ticks on the axes.\n"
                  "- xlabel, ylabel: font settings for the axes labels.\n"
                  "- legend: font settings for the legend box.\n"
                  "- list: show the font settings for the entire plot.\n"
                  "- x-twin-ticks, y-twin-ticks: ticks fonts for twin axes.\n"
                  "- x-twin-label, y-twin-label: font settings for twin axes labels."),
                 ("For 'axes', 'labels' and 'legend' font settings, the size, weight, variant "
                  "and style of the font can be set. Use the tab completion to select valid "
                  "values. To type without tab completion, specify the property value set as "
                  "\"property=value\", for example \"size=large\"."),
                 ("To delete a font specification for one of the options above, specify 'rm'. "
                  "EXAMPLE: \"font family rm\" or \"font yticks rm\"."),
                 ("Since it is common to use the same font settings for plots across multiple "
                  "papers, you can save and load a set of font settings for an analysis group "
                  "using the 'fontsave' and 'fontload' commands.")]
        self._fixed_width_info(lines)
    def complete_font(self, text, line, istart, iend):
        els = line.split()
        if len(els) == 1 or (len(els)==2 and line[-1] != " "):
            #We use the font command for multiple options, list the possibilities
            #in the second position.
            if len(els) == 2:
                part = els[1]
            else:
                part = ""
            options = ["family", "xlabel", "ylabel", "xticks", "yticks", "list", "legend", "title",
                       "x-twin-ticks", "x-twin-label", "y-twin-ticks", "y-twin-label"]
            return [o for o in options if o.startswith(part)]
        else:
            #We need to look at the second element to decide how to complete.
            if line[-1] == " ":
                part = ""
            else:
                part = els[-1]

            #We want to allow repitition of font options on the same line of the shell.
            #If the last element has an equal sign, then we are completing the font
            #property.
            if "list" not in line and "family" not in line and len(els) > 2:
                cfonts = {
                    "size": ['xx-small', 'x-small', 'small', 'medium', 'large',
                             'x-large', 'xx-large'],
                    "variant": ['normal', 'small-caps'],
                    "style": ['normal', 'italic', 'oblique'],
                    "weight": ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'],
                }
                if "=" in part:
                    prop, val = part.split("=")
                    return [p for p in cfonts[prop] if p.startswith(val)]
                else:
                    if part in cfonts:
                        return [part + "="]
                    else:
                        return [p for p in cfonts if p.startswith(part)]            
            elif part != "list":
                propkeys = ["size", "variant", "style", "weight", "rm"]
                compl = {
                    "family": ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace', 'rm'],
                    "xticks": propkeys,
                    "yticks": propkeys,
                    "xlabel": propkeys,
                    "ylabel": propkeys,
                    "legend": propkeys,
                    "title": propkeys
                }
                #See if we have already specified which font option we are editing.
                return [p for p in compl[els[1]] if p.startswith(part)]

    def do_fontsave(self, arg):
        """Saves the session variables to a file so that the same analysis can be continued later."""
        #We need to save the full paths to the staging directories for the tests that are loaded
        #so far; then when the session is restored, we can reparse the results.
        from os import path
        import json
        fullpath = path.expanduser(arg)
        data = {
            "fonts": self.curargs["fonts"],
            "ticks": self.curargs["ticks"]
        }
        with open(fullpath, 'w') as f:
            json.dump(data, f)
        msg.okay("Saved current font settings to {}".format(fullpath))
    def complete_fontsave(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)
    def help_save(self):
        lines = [("Saves the current font settings for plots to the specified file."),
                 ("EXAMPLE: \"fontsave pubfonts.json\" saves the font settings to the specified "
                  "file. The session is always serialized in JSON."),
                 ("See also: 'fontload'")]
        self._fixed_width_info(lines)

    def do_fontload(self, arg):
        from os import path
        import json
        fullpath = path.expanduser(arg)
        if path.isfile(fullpath):
            with open(fullpath) as f:
                data = json.load(f)
            self.curargs["fonts"] = data["fonts"]
            self.curargs["ticks"] = data["ticks"]
            msg.okay("Loaded font settings from {}".format(fullpath))
    def complete_fontload(self, text, line, istart, iend):
        return self.complete_parse(text, line, istart, iend)
    def help_load(self):
        lines = [("Loads previously saved font settings from disk. For details, see 'fontsave'. "
                  "NOTE: the font settings in the active analysis group will be overwritten."),
                 ("EXAMPLE \"fontload pubfonts.json\" loads the saved session.")]
        self._fixed_width_info(lines)

    def _complete_dep_keyval(self, line, propkeys, opts):
        els = line.split()
        if len(els) == 1 or (len(els) == 2 and line[-1] != " "):
            return self._complete_deps(els, line)
        else:
            #We need to look at the second element to decide how to complete.
            if line[-1] == " ":
                part = ""
            else:
                part = els[-1]

            if len(els) > 2:
                if "=" in part:
                    prop, val = part.split("=")
                    if prop in opts:
                        return [p for p in opts[prop] if p.startswith(val)]
                else:
                    if part in propkeys and part != "rm":
                        return [part + "="]
                    else:
                        return [p for p in propkeys if p.startswith(part)]
            elif part != "list":
                return [p for p in propkeys if p.startswith(part)]        

    def _do_dep_keyval_single(self, option, varname, vals):
        if vals[1] == "rm":
            if varname in self.curargs[option]:
                del self.curargs[option][varname]
        else:
            props = [v.split("=") for v in vals[1:len(vals)]]
            if varname in self.curargs[option]:
                #Update the existing values or insert new ones where necessary.
                target = self.curargs[option][varname]
            else:
                target = {}

            for pname, pval in props:
                target[pname] = pval
            if varname not in self.curargs[option]:
                self.curargs[option][varname] = target

    def _do_dep_keyval(self, arg, option, listfun, heading):
        vals = arg.split()
        if len(vals) >= 2:
            varname = vals[0]
            if varname == "*":
                for depvar in self.curargs["dependents"]:
                    self._do_dep_keyval_single(option, depvar, vals)
            else:
                self._do_dep_keyval_single(option, varname, vals)
            listfun("list")
        elif arg == "list":
            for var in self.curargs[option]:
                msg.info("'{}' {} SETTINGS".format(var, heading.upper()))
                for pname, pval in list(self.curargs[option][var].items()):
                    msg.info("  {} => {}".format(pname, pval))
        else:
            msg.warn("You haven't specified a valid variable and property combination.")

    def do_markers(self, arg):
        self._do_dep_keyval(arg, "markers", self.do_markers, "marker")
    def help_markers(self):
        lines = [("Sets the style of the markers for one of the dependent variables. "
                  "Use <tab> completion to select the variable (using filter, file name, "
                  "and property) and then select marker options from the lists. You can "
                  "remove the marker settings for a variable by specifying 'rm' as the "
                  "property name."),
                 ("EXAMPLE \"markers */concs.in|depth size=2 marker=tri_down fill=full\" "
                  "sets the size to 2 points^2, the marker to be an upside-down triangle "
                  "and the filling to be full.")]
        self._fixed_width_info(lines)
    def complete_markers(self, text, line, istart, iend):
        propkeys = ["marker", "size", "fill", "rm"]        
        opts = {
            "marker": list(self._possible_markers.values()),
            "fill": ('full', 'left', 'right', 'bottom', 'top', 'none')
        }
        return self._complete_dep_keyval(line, propkeys, opts)

    def do_lines(self, arg):
        self._do_dep_keyval(arg, "lines", self.do_lines, "line")
    def help_lines(self):
        lines = [("Sets the width and style of lines drawn on plots for a specific variable."
                  "Use <tab> completion to select the variable (using filter, file name, "
                  "and property) and then select line options from the lists. You can "
                  "remove the line settings for a variable by specifying 'rm' as the "
                  "property name."),
                 ("EXAMPLE \"lines */concs.in|depth.fit width=2 style=_draw_dashed\" "
                  "sets the line width to 2 points and the style to dashed.")]
        self._fixed_width_info(lines)
    def complete_lines(self, text, line, istart, iend):
        propkeys = ["style", "width", "rm"]        
        opts = {
            "style": list(self._possible_linestyles.values())
        }
        return self._complete_dep_keyval(line, propkeys, opts)

    def do_ticks(self, arg):
        els = arg.split()
        if len(els) > 0 and arg != "list":
            cast = ["length", "width", "pad"]
            dprops = {}
            delete = False
            for entry in els:
                if "=" in entry:
                    prop, val = entry.split("=")
                    if prop in cast:
                        dprops[prop] = float(val)
                    else:
                        if prop == "axis":
                            #Handle the twin-x and twin-y. A little tricky because when
                            #we think about handling the twin-x ticks, we really mean
                            #that we want to adjust the y-axis.
                            if "-" in val:
                                dprops["axis"] = "y" if val == "x-twin" else "x"
                            else:
                                dprops["axis"] = val
                        else:
                            dprops[prop] = val
                elif entry == "rm":
                    delete = True

            if "axis" in dprops:
                axis = dprops["axis"]
            else:
                axis = "both"
            if "which" in dprops:
                which = dprops["which"]
            else:
                which = "major"

            key = "{}|{}".format(axis, which)
            if delete and key in self.curargs["ticks"]:
                del self.curargs["ticks"][key]
            else:
                self.curargs["ticks"][key] = dprops
            self.do_ticks("list")
        elif arg == "list":
            for key in self.curargs["ticks"]:
                axis, which = key.split("|")
                msg.info("AXIS: {} and WHICH: {}".format(axis.upper(), which.upper()))
                for prop, val in list(self.curargs["ticks"][key].items()):
                    msg.info("  {} => {}".format(prop, val))                    
    def help_ticks(self):
        lines = [("Sets the properties for drawing the tick marks on plots for the active "
                  "analysis group. Use the tab completion to select properties and valid "
                  "values for the ticks. To style the x and y axes differently, use the "
                  "'axis' property with different settings each time. The settings are stored "
                  "according to the 'axis' and 'which' properties and are applied in an"
                  " arbitrary order. You can "
                  "use \"axis=both\" to set common settings and then specialize with the other "
                  "possibilities."),
                 ("To remove settings for an axis-which combination, specify the axis and "
                  "which keywords and then use 'rm' as a keyword. E.g. \"ticks rm\" would "
                  "remove the settings for the *major* ticks on *both* axes (the default "
                  "axis and which values). You could also use \"ticks rm axis=x which=minor\" "
                  "to remove settings for that set.")]
        self._fixed_width_info(lines)
    def complete_ticks(self, text, line, istart, iend):
        bopts = ("on", "off")
        opts = {
            "axis": ("x", "y", "both", "x-twin", "y-twin"),
            "which": ("major", "minor", "both"),
            "direction": ("in", "out", "inout"),
            "color": self._possible_cols,
            "bottom": bopts,
            "top": bopts,
            "left": bopts,
            "right": bopts,
            "labelbottom": bopts, 
            "labeltop": bopts,
            "labelleft": bopts,
            "labelright": bopts,
            "reset": ("true", "false")
        }
        propkeys = list(opts.keys()) + ["length", "width", "pad", "list", "rm"]
        
        els = line.split()
        if line[-1] == " ":
            part = ""
        else:
            part = els[-1]

        if '=' in part:
            prop, val = part.split('=')
            if prop in opts:
                return [p for p in opts[prop] if p.startswith(val)]
        else:
            if part in propkeys:
                return [part + '=']
            else:
                return [p for p in propkeys if p.startswith(part)]

    def do_limit(self, arg):
        import re
        vals = arg.split()
        if len(vals) >= 2:
            if vals[0] in ["x", "y", "z", "x-twin", "y-twin"]:
                if re.match("[\d.]+", vals[1]):
                    self.curargs["limits"][vals[0]] = tuple(map(float, vals[1:3]))
                elif vals[1] == "rm":
                    del self.curargs["limits"][vals[0]]
                elif vals[1] == "auto":
                    self.curargs["limits"][vals[0]] = "auto"
                self.do_limit("list")
            else:
                msg.err("Only 'x', 'y', 'z', 'x-twin' and 'y-twin' are valid axis designations.")                
        elif arg == "list":
            for dim in self.curargs["limits"]:
                if isinstance(self.curargs["limits"][dim], tuple):
                    msg.info("{0} LIMITS: {1[0]}-{1[1]}".format(dim.upper(), self.curargs["limits"][dim]))
                else:
                    msg.info("{0} LIMITS: {1}".format(dim.upper(), self.curargs["limits"][dim]))
        else:
            msg.warn("Enter the axis type ('x', 'y', 'z', 'x-twin', 'y-twin') and the"
                     " start and end values. E.g. \"limit x 0 20\"")
    def help_limit(self):
        lines = [("Sets the limiting values for the axes on the plot."),
                 ("EXAMPLE \"limit x 0 20\" sets the x-axis to only plot from zero to twenty.")]
        self._fixed_width_info(lines)
    def complete_limit(self, text, line, istart, iend):
        els = line.split()
        if len(els) == 1 or (len(els) == 2 and line[-1] != ' '):
            axes = ["x", "y", "z", "list", "x-twin", "y-twin"]
            if len(els) == 2:
                return [a for a in axes if a.startswith(els[1])]
            else:
                return axes
        elif (len(els) == 2 and line[-1] == " ") or len(els) == 3:
            opts = ["<float>", "rm", "auto"]
            if len(els) == 3:
                return [o for o in opts if o.startswith(text)]
            else:
                return opts

    def do_figsize(self, arg):
        vals = arg.split()
        if len(vals) >= 2:
            self.curargs["figsize"] = tuple(map(float, vals))
        elif arg == "list":
            if "figsize" in self.curargs:
                msg.info("FIGSIZE: {0[0]}x{0[1]}".format(self.curargs["figsize"]))
            else:
                msg.info("FIGSIZE: {0[0]}x{0[1]}".format((7,5)))
        else:
            self.help_figsize()
    def help_figsize(self):
        lines = [("Sets the size of the figure for the plot."),
                 ("EXAMPLE \"figsize 10 20\" sets the figure size to 10x20 (WxH) *inches*.\n"
                  "EXAMPLE \"figsize 10 20 0.8\" also sets the figure padding to 0.8 (fraction of text width).")]
        self._fixed_width_info(lines)
    def complete_figsize(self, text, line, istart, iend):
        els = line.split()
        if len(els) == 1 or (len(els) == 2 and line[-1] != ' '):
            options = ["list"]
            if len(els) == 2:
                return [a for a in options if a.startswith(els[1])]
            else:
                return options

    def do_legend(self, arg):
        usable, filename, append = self._redirect_split(arg)
        vals = usable.split()
        if len(vals) == 1 and vals[0] == "list":
            self._print_map_dict("legend", filename, append)
        elif len(vals) >= 1:
            if len(vals) >= 2 and vals[0] == "rm":
                if vals[1] == "*":
                    for pkey in list(self.curargs["legend"].keys()):
                        del self.curargs["legend"][pkey]
                else:
                    for pkey in vals[1:]:
                        if pkey in self.curargs["legend"]:
                            del self.curargs["legend"][pkey]
            else:
                for pkey in vals:
                    if "=" in pkey:
                        prop, val = pkey.split("=")
                        if val == "none":
                            value = None
                        elif val in ["true", "false"]:
                            value = val == "true"
                        else:
                            #Check what the plotter expects and cast the value.
                            dst = self.legend_opts[prop]
                            if isinstance(dst, tuple):
                                #We are expecting a list, cast each value and return a list.
                                value = list(map(dst[1], val[1:-1].split(",")))
                            elif isinstance(dst, dict):
                                value = val.replace("_", " ")
                            else:
                                value = dst(val)
                        self.curargs["legend"][prop] = value
                    else:
                        msg.err("'{}' is an invalid property-value combination.".format(pkey))
            self.do_legend("list")
    def help_legend(self):
        lines = [("Sets the properties for displaying the legend."),
                 ("EXAMPLE \"legend pos=right shadow=true\" sets the legend position to the right "
                  "of the figure and includes a shadow around the box.")]
    def complete_legend(self, text, line, istart, iend):
        suggest = {
            str: ["<string>"],
            bool: ["true", "false"],
            int: ["<int>"],
            float: ["<float>"]
        }
        defaults = {
            "markerfirst": True,
            "ncol": 1
        }
        els = line.split()
        if line[-1] == " ":
            if len(els) <= 2:
                return list(self.legend_opts.keys()) + ["list", "rm"]
            else:
                return list(self.legend_opts.keys())
        elif len(els) >= 2:
            if els[1] == "rm":
                #Return the list of properties that are set so far.
                return [k for k in self.curargs["legend"] if k.startswith(els[-1])] + ["*"]
            else:
                #We specify keyword value combinations.
                if "=" in els[-1] and line[-1] != " ":
                    prop, val = els[-1].split("=")
                    if prop not in self.legend_opts:
                        msg.warn("'{}' is not a valid legend property.".format(prop))
                    else:
                        skey = self.legend_opts[prop]
                        if isinstance(skey, dict):
                            return ["{}={}".format(prop, k) for k in skey if k.startswith(val)]
                        elif isinstance(skey, tuple):
                            return ["{0}=[{1},{1},...]".format(prop, str(skey[1]).replace("type ", "").replace("'", "")),
                                    "{}=none".format(prop)]
                        else:
                            if skey in defaults:
                                oplist = suggest[skey] + [defaults[skey]]
                            else:
                                oplist = suggest[skey] + ["none"]
                            return ["{}={}".format(prop, o) for o in oplist if o.startswith(val)]
                else:
                    oplist = [p + '=' for p in self.legend_opts if p.startswith(els[-1])]
                    if "list".startswith(els[-1]) and len(els) == 2:
                        oplist.append("list")
                    if "rm".startswith(els[-1]) and len(els) == 2:
                        oplist.append("rm")
                    return oplist
            
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
