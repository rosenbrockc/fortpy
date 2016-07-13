"""Module for analyzing the results of multiple test cases against each other
for a single unit test.
"""
from fortpy import msg
def _fit_exp(x, a, b, c):
    from numpy import exp
    return a * exp(b * x) + c

def _fit_lin(x, a, b):
    return a*x + b

class Analysis(object):
    """Represents the results of analyzing all the test case directories."""
    def __init__(self, stagedir, fullparse=False):
        """Analyzes the results from the test cases in the specified staging directory.

        :arg stagedir: the full path to the unit test directory.
        :arg fullparse: when true the entire contents of input/output files is also parsed
          and interpreted so that it is available for plotting along with other results. This
          increases the parsing time significantly.
        """
        from os import path
        self.stagedir = path.expanduser(stagedir)
        self.infiles = []
        """A list of input file names common to all the directories, and which can
        act as (in)dependent variables for further plotting or analysis.
        """
        self.outfiles = []
        """A list of output file names common to all the directories, and which can
        act as (in)dependent variables for further plotting or analysis.
        """
        self.details = {}
        """The dictionary with all the detailed analysis of each test case directory.
        Created by self.parse()."""
        self.props = {}
        """A dictionary of properties available for each file name. Assumes that the
        first file of a given name found has the same properties available as all
        subsequent ones."""
        self.comparable = []
        """A list of file names for which we also have a 'filename.compare' file with
        the comparison results from the unit test runs.
        """
        self.fits = {}
        """A dict of the fitting results for a given set of indpendent and dependent
        variables."""

        self.parse(fullparse)

        self._vars = None

    @property
    def allvars(self):
        """Return a list of possible file names that can be used as variables during analysis."""
        if self._vars is None:
            self._vars = []
            self._vars.extend(self.outfiles)
            self._vars.extend(self.infiles)
        return self._vars

    def _test_percent(self, testcase, outfile):
        """Retrieves the percent success for the specified output file and test case.

        :arg testcase: the dictionary of parsed test result data to use in analysis.
        :arg outfile: a specific output file present in self.outfiles that should have
          its percent accuracy checked for success.        
        """
        result = None
        comparefile = "{}.compare".format(outfile)
        if comparefile in testcase["percents"]:
            result = testcase["percents"][comparefile]
        else:
            #Since the user isn't explicitly comparing the result, it must be true.
            result = 1.0
        return result        

    def _is_successful(self, testcase, outfile=None, threshold=1.):
        """Determines whether the dictionary representing the specified test case was
        successful or not and whether it can be used.

        :arg testcase: the dictionary of parsed test result data to use in analysis.
        :arg outfile: a specific output file present in self.outfiles that should have
          its percent accuracy checked for success.
        :arg threshold: float specifying the minimum percentage required for the test to
          considered successful with respect to the specified outfile.
        """
        result = False
        if "lastrun" in testcase:
            #We know that at least it completed, so there may be a 'outfile.compare'
            #report for the output that describes how well it did. It is also possible
            #that the output file wasn't set for comparison, but was being created by
            #fortpy for another purpose.
            if outfile is not None and ".out" in outfile:
                percent = self._test_percent(testcase, outfile)
                result = (percent is not None and percent >= threshold)
            else:
                result = True
            
        return result

    def _failures_data(self, variables=None, threshold=1., tfilter=None):
        """Returns a list of test cases that failed and their completion percentage for
        each variable, optionally filtered.

        :arg variable: (list) filter the output file names whose values are considered.
        :arg threshold: the percentage value that defines success. Values below the threshold
          are considered failures.
        :arg tfilter: only test cases whose names contain this value are included in the search.
        """
        #Failures only include the outputs in the case folders.
        result = {}
        for case in self.details:
            if not self._case_filter(case, tfilter):
                continue

            rset = {}
            for outfile in self.details[case]["outputs"]:
                if variables is not None and outfile not in variables:
                    continue
                if not self._is_successful(self.details[case], outfile, threshold):
                    rset[outfile] = self._test_percent(self.details[case], outfile)

            if len(rset) > 0:
                result[case] = rset

        return result

    def failures(self, variables=None, threshold=1., tfilter=None):
        """Formats a list of test cases who failed.

        :arg variable: (list) filter the output file names whose values are considered.
        :arg threshold: the percentage value that defines success. Values below the threshold
          are considered failures.
        :arg tfilter: only test cases whose names contain this value are included in the search.
        """
        #We just need to format the failures data into a nice interactive list format.
        if tfilter is None:
            title = "Failures for *all* Test Cases"
        else:
            title = "Failures for Test Cases Matching '{}'".format(tfilter)
        printed = [title, ''.join(['-']*50)]
        data = self._failures_data(variables, threshold, tfilter)
        keys = list(sorted(data.keys()))
        #In order to facilitate console interaction with the results, we need to number the
        #output values.
        for i in range(len(keys)):
            results = ["{}. {}:".format(i+1, keys[i].upper())]
            fkeys = list(sorted(data[keys[i]].keys()))
            for j in range(len(fkeys)):
                if data[keys[i]][fkeys[j]] is not None:
                    results.append("   {0}: '{1}' => {2:.2f}%".format(j+1, fkeys[j], 
                                                                      data[keys[i]][fkeys[j]]*100))
                else:
                    results.append("   {0}: '{1}' => None".format(j+1, fkeys[j]))
            printed.append('\n'.join(results) + '\n')

        return '\n'.join(printed)

    def _case_filter(self, caseid, tfilter):
        """Determines whether the specified filter matches the case.

        :arg caseid: the identifier for the specific test case.
        :arg tfilter: the case filter; if '*' is present in the filter, then the case
          matches true as long as it contains the filter text; otherwise the filter needs
          to match exactly.
        """
        if tfilter is None or tfilter == "*":
            return True
        else:
            from fnmatch import fnmatch
            return (("*" in tfilter and fnmatch(caseid, tfilter)) or
                    (not "*" in tfilter and tfilter == caseid))

    def _get_data(self, variable, fullvar, order=None, threshold=1., tfilter=None,
                  independent=None, x=None):
        """Returns a raw list of the valid data points for the specified variable.

        :arg variable: a string indentifying the variable's filename and the
          appropriate property of the data in the file. For example: "concs.in|depth" would
          use the depth of the data in the file 'concs.in' as its value for each test case.
        :arg fullvar: the full variable name including the filter and property spec.
        :arg threshold: a float value specifying the minimum level of success that the output
          file must attain before it can be used in the plot.
        :arg independent: the name of the independent variable including property.
        :arg x: a list of values for the *independent* variable to use when evaluating
          a *fitted* curve to a variable.
        """
        from numpy import ndarray
        target = None
        if variable != "timing" or "timing|" in variable:
            varfile, attribute = variable.split("|")
        else:
            target = varfile = "timing"
            attribute = None

        if target is None:
            if varfile in self.infiles:
                target = "inputs"
            elif varfile in self.outfiles:
                target = "outputs"
            elif varfile == "timing" and attribute is not None:
                target = "timing"
            else:
                raise ValueError("Cannot locate the values for variable {}".format(variable))

        values = []
        cases = []
        #If we have already ruled out some of the test cases while doing the independent
        #variable, we want to only use those cases.
        if order is not None:
            details = order
        else:
            details = self.details

        for testcase in details:
            if not self._case_filter(testcase, tfilter):
                #We don't consider test cases that are blocked by the filter.
                continue

            value = 0
            if target == "timing" and self._is_successful(self.details[testcase]) and attribute is None:
                value = self.details[testcase][target]
            else:
                if (self._is_successful(self.details[testcase], varfile, threshold)):
                    if target != "timing" and varfile in self.details[testcase][target]:
                        data = self.details[testcase][target][varfile]
                    else:
                        data = None

                    if data is not None and attribute in data:
                        value = data[attribute]
                    elif attribute == "fit" and independent is not None and x is not None:
                        #Return the value of the *fitted* function for the corresponding independent
                        #variable value; this is one-to-one for tabulating. For plotting we don't use
                        #these values.
                        key = "{}${}".format(independent, varfile)
                        if len(values) < len(x) and key in self.fits:
                            ival = x[len(values)]
                            value = self.fits[key]["function"](ival)

            if isinstance(value, ndarray) or value != 0 or order is not None:
                values.append(value if value is not None else 0)
                cases.append(testcase)

        return (values, cases)

    def _get_data_series(self, independent, dependents, threshold, functions=None):
        """Used by plot() and table() to get the x and y data series that will be plotted
        or tabulated against each other. See the arguments on those methods.
        """
        if (("rowvals" in independent or "colvals" in independent) or
            any([("rowvals" in d or "colvals" in d) for d in dependents])):
            #We need to make sure that the test case filter they specified returns only a
            #single result so that we get a reasonable plot. This will happen if the data
            #is 1D in the other dimension or if the filter only has a single case. We
            #issue a warning if the filter doesn't match up, so they can check their data.
            for dep in dependents:
                tfilter, depvar = dep.split("/")
                if tfilter is None or sum([(1 if self._case_filter(d, tfilter) else 0)
                                           for d in self.details]) > 1:
                    msg.warn("Plotting aggregated data for more than one test case. Check results \n"
                             "for consistency and completeness.")
                    break

        rawvars = {}
        oseries = []
        icases = {}
        def _load_raw(variable, rawvars):
            if variable in rawvars:
                return
            
            tfilter, depvar = variable.split("/")
            fullindvar = "{}/{}".format(tfilter, independent)
            cases = None if fullindvar not in icases else icases[fullindvar]
            if fullindvar not in rawvars:
                x, cases = self._get_data(independent, fullindvar, None, threshold, tfilter)
                if (len(x) == 1 and isinstance(x[0], list) and 
                    ("rowvals" in independent or "colvals" in independent)):
                    x = x[0]
                rawvars[fullindvar] = x

            if cases is not None:
                x = rawvars[fullindvar]
                ypts, names = self._get_data(depvar, variable, cases, threshold, tfilter, independent, x)
                if (len(ypts) == 1 and isinstance(ypts[0], list) and
                    ("rowvals" in variable or "colvals" in variable)):
                    ypts = ypts[0]
                rawvars[variable] = ypts
            return (fullindvar, cases)

        for variable in dependents:
            fullindvar, cases = _load_raw(variable, rawvars)
            icases[fullindvar] = cases
            oseries.append((fullindvar, variable))                    
        #Next, we need to apply any postfix functions and then check that the data sizes
        #are commensurate (i.e. same number of points for dependent and independent series).
        ys = []
        xs = []        
        if functions is not None:
            for postfix, fixdef in functions.items():
                fxn = None
                vorder = [v for v in fixdef.keys() if v != "lambda"]
                values = {}
                if postfix not in rawvars:
                    msg.err("Variable '{}' to be postfixed has no data series.".format(postfix))
                    break
                
                N = len(rawvars[postfix])
                for lvar, gvar in fixdef.items():
                    if lvar == "lambda":
                        #This is the function definition for each item, create a lambda for it.
                        import numpy
                        import math
                        try:
                            evals = gvar.split(":")[1]
                            lambstr = "lambda {}: {}".format(", ".join(vorder), evals)
                            fxn = eval(lambstr)
                        except:
                            msg.err("Could not evaluate function '{}'.".format(gvar))
                    else:
                        if gvar not in rawvars:
                            _load_raw(gvar, rawvars)
                        if gvar in rawvars:
                            values[lvar] = rawvars[gvar]
                        else:
                            emsg = "Variable '{}' in the lambda function postfix for '{}' is missing data."
                            msg.err(emsg.format(gvar, postfix))

                #Now if we have a valid function defined and data to operate on, just
                #evaluate the postfix one value at a time.
                if fxn is not None and N > 0:
                    for i in range(N):
                        args = [values[v][i] for v in vorder]
                        pval = fxn(*args)
                        rawvars[postfix][i] = pval

        #Now we can use the adjusted/postfixed raw values to construct the data series
        for indepvar, depvar in oseries:
            xs.append((rawvars[indepvar], indepvar))
            #We need the arrays to be the same length for plotting; if we are following
            #an existing independent variable, we still need to append zero to the list
            #for the values we don't have.
            if len(rawvars[depvar]) != len(rawvars[indepvar]):
                msg.err("Can't coerce an array to a single value without an aggregation "
                        "function such as numpy.mean or numpy.sum")
            ys.append((rawvars[depvar], depvar))
            
        return (xs, ys)

    def _get_font(self, dfont, fonts, option):
        """Gets a copy of the default font specified and updates its properties based on the
        values in the fonts dictionary.

        :arg option: the name of the plot item whose font is being changed. Possible values are
          ['axes', 'labels', 'legend'].
        """
        if fonts is not None and option in fonts:
            fdict = fonts[option]
            xfont = dfont.copy()
            if "size" in fdict:
                xfont.set_size(fdict["size"])
            if "style" in fdict:
                xfont.set_style(fdict["style"])
            if "variant" in fdict:
                xfont.set_variant(fdict["variant"])
            if "weight" in fdict:
                xfont.set_weight(fdict["weight"])

            return xfont
        else:
            return dfont

    def _get_markers(self, dargs, markers):
        m = markers.get("marker")
        f = markers.get("fill", "full")
        #This is a bug in matplotlib that I submitted an issue on at github.
        #For the meanwhile, only full fill marker styles. So we have this hack:
        if f == "none":
            dargs["facecolors"] = f
        dargs["edgecolors"] = dargs["color"]
        # When the bug is fixed, we can use this instead MarkerStyle(marker=m, fillstyle=f)
        dargs["marker"] = m

    def _reset_lineplot(self, dargs):
        if "marker" in dargs and "linestyle" not in dargs:
            #Set the default line style to solid line.
            dargs["linestyle"] = "-"
        if "facecolors" in dargs:
            dargs["markerfacecolor"] = dargs["facecolors"]
            del dargs["facecolors"]
        if "s" in dargs:
            dargs["markersize"] = dargs["s"]
            del dargs["s"]
        #Set the marker edge color if markers are specified; otherwise they show up gray.
        if "marker" in dargs:
            dargs["markeredgecolor"] = dargs["color"]
            if "edgecolors" in dargs:
                del dargs["edgecolors"]

    def plot(self, independent=None, dependents=None, threshold=1.,
             savefile=None, functions=None, xscale=None, yscale=None,
             colors=None, labels=None, fonts=None, markers=None, lines=None, ticks=None,
             plottypes=None, limits=None, twinplots=None, legend=None, **kwargs):
        """Plots the specified dependent variables as functions of the independent one.

        :arg independent: a string indentifying the independent variable's filename and the
          appropriate property of the data in the file. For example: "concs.in|depth" would
          use the depth of the data in the file 'concs.in' as its value for each test case.
        :arg dependent: a list of strings that use the same format as the string for the
          'independent' variable.
        :arg threshold: a float value specifying the minimum level of success that the output
          file must attain before it can be used in the plot.
        :arg xlabel: the label for the x-axis of the plot.
        :arg ylabel: the label for the y-axis of the plot.
        :arg functions: a dictionary of functions to apply to the values of each variable. E.g.
          {"group.in|depth": {lambda dict}}. See the help text help_postfix() for details.
        """
        if independent is None or dependents is None:
            raise ValueError("Must specify at least the variables to plot.")
        
        import matplotlib.pyplot as plt
        from matplotlib import cm, rc
        from numpy import linspace, array
        from itertools import cycle
        from matplotlib.font_manager import FontProperties

        xs, ys = self._get_data_series(independent, dependents, threshold, functions)
        
        #Set the default font for the plots, or the custom one.
        if fonts is not None and "family" in fonts:
            rc('font',family=fonts["family"])
        else:
            rc('font',family='Times New Roman')
            
        layout = None
        if "figsize" in kwargs:
            if len(kwargs["figsize"]) > 2:
                layout = {"pad": kwargs["figsize"][2]}
                
        fig = plt.figure(figsize=(None if "figsize" not in kwargs else kwargs["figsize"][0:2]), tight_layout=layout)
        dfont = FontProperties()

        ax = fig.add_subplot(111)
        rbcolors = cm.rainbow(linspace(0, 1, len(ys)))
        cycols = cycle(rbcolors)

        axx = None
        axy = None
        if twinplots is not None:
            if "twinx" in twinplots.values():
                axx = ax.twinx()
            if "twiny" in twinplots.values():
                axy = ax.twiny()

        if labels is not None:
            if "x" in labels:
                ax.set_xlabel(labels["x"], fontproperties=self._get_font(dfont, fonts, "xlabel"))
            if "y" in labels:
                ax.set_ylabel(labels["y"], fontproperties=self._get_font(dfont, fonts, "ylabel"))
            if "plot" in labels:
                plt.title(labels["plot"], fontproperties=self._get_font(dfont, fonts, "title"))
            if "x-twin" in labels and axx is not None:
                axx.set_ylabel(labels["x-twin"], fontproperties=self._get_font(dfont, fonts, "x-twin-label"))
            if "y-twin" in labels and axy is not None:
                axy.set_xlabel(labels["y-twin"], fontproperties=self._get_font(dfont, fonts, "y-twin-label"))
                
        if xscale is not None:
            ax.set_xscale(xscale)
        if yscale is not None:
            ax.set_yscale(yscale)
                
        for xset, yset in zip(xs, ys):
            x, xlabel = xset
            y, ylabel = yset

            #Set up a dictionary with all the kwargs for the plotting options
            #that may be common to both kinds of plots.
            dargs = {}
            if colors is not None and ylabel in colors:
                dargs["color"] = colors[ylabel]
            else:
                dargs["color"] = next(cycols)

            if markers is not None and ylabel in markers:
                self._get_markers(dargs, markers[ylabel])
                if "size" in markers[ylabel]:
                    dargs["s"] = float(markers[ylabel]["size"])
            if "s" not in dargs:
                dargs["s"] = 15
                    
            if lines is not None and ylabel in lines:
                if "style" in lines[ylabel]:
                    dargs["linestyle"] = lines[ylabel]["style"]
                if "width" in lines[ylabel]:
                    dargs["linewidth"] = float(lines[ylabel]["width"])

            if labels is not None and ylabel in labels:
                dargs["label"] = None if labels[ylabel].lower() == "[none]" else labels[ylabel] 
            else:
                dargs["label"] = ylabel

            if twinplots is not None and ylabel in twinplots:
                if twinplots[ylabel] == "twinx":
                    liveax = axx
                elif twinplots[ylabel] == "twiny":
                    liveax = axy
            else:
                liveax = ax
            try:
                if plottypes is not None and ylabel in plottypes:
                    lineplot = plottypes[ylabel] == "line"
                else:
                    lineplot = None

                if ylabel[len(ylabel)-4:len(ylabel)] in ["|fit", ".fit"]:
                    varfile = ylabel[0:len(ylabel)-4]
                    key = "{}${}".format(independent, varfile)
                    allx = linspace(min(x), max(x), 50)
                    if labels is not None and ylabel in labels:
                        dargs["label"] = labels[ylabel].format(self._format_fit(key))
                        if dargs["label"].lower() == "[none]":
                            dargs["label"] = None
                    else:
                        dargs["label"] = self._format_fit(key)

                    self._reset_lineplot(dargs)
                    liveax.plot(allx, self.fits[key]["function"](allx), **dargs)
                elif lineplot:
                    self._reset_lineplot(dargs)
                    from operator import itemgetter
                    sdata = array(sorted(zip(x,y),key=itemgetter(0)))
                    liveax.plot(sdata[:,0], sdata[:,1], **dargs)
                else:
                    size = dargs["s"]
                    del dargs["s"]
                    liveax.scatter(x, y, s=size, **dargs)
            except ValueError:
                msg.err("The values for {} can't be log-plotted.".format(ylabel))

        #Set the tick labels font sizes.
        def tickfonts(fonts, key, labels):
            if fonts is not None and key in fonts:
                for label in labels:
                    tfont = self._get_font(dfont, fonts, key)
                    label.set_fontproperties(tfont)
        tickfonts(fonts, "xticks", ax.get_xticklabels())
        tickfonts(fonts, "yticks", ax.get_yticklabels())
        if axx is not None:
            tickfonts(fonts, "x-twin-ticks", axx.get_yticklabels())
        if axy is not None:
            tickfonts(fonts, "y-twin-ticks", axy.get_xticklabels())

        if ticks is not None and len(ticks) > 0:
            for key in ticks:
                if "reset" in ticks[key]:
                    ticks[key]["reset"] = ticks[key]["reset"] == "true"
                if "x-twin" in key and axx is not None:
                    axx.tick_params(**ticks[key])
                elif "y-twin" in key and axy is not None:
                    axy.tick_params(**ticks[key])
                else:
                    ax.tick_params(**ticks[key])

        for dim in ["x", "y", "z", "x-twin", "y-twin"]:
            if limits is not None and dim in limits:
                if dim == "x-twin" and axx is not None:
                    if isinstance(limits[dim], tuple):
                        axx.set_ylim(limits[dim])
                    elif limits[dim] == "auto":
                        axx.set_ylim(auto=True)
                elif dim == "y-twin" and axy is not None:
                    if isinstance(limits[dim], tuple):
                        axy.set_ylim(limits[dim])
                    elif limits[dim] == "auto":
                        axy.set_ylim(auto=True)                        
                elif "-" not in dim:
                    if isinstance(limits[dim], tuple):
                        getattr(ax, "set_{}lim".format(dim))(limits[dim])
                    elif limits[dim] == "auto":
                        getattr(ax, "set_{}lim".format(dim))(auto=True)

        if len(dependents) > 1:
            if legend is not None:
                plt.legend(prop=self._get_font(dfont, fonts, "legend"), **legend)
            else:
                plt.legend(loc='upper right', prop=self._get_font(dfont, fonts, "legend"))
        if savefile is None:
            plt.show(block=False)
        else:
            plt.savefig(savefile)

    def _format_heading(self, heading):
        """Formats the specified heading text to conform to the table layout."""
        if "/" in heading:
            filt, rest = heading.split("/")
            cases = filt.split(".")
            fstr = '.'.join(cases[1:len(cases)])
            heading = "{}/{}".format(fstr, rest)
        if len(heading) > 15:
            return "{0:^15}".format(heading[0:15])
        else:
            return "{0:^15}".format(heading)

    def table(self, independent, dependents, threshold=1., headings=None, functions=None):
        """Outputs the values for the specified dependent and independent variables in table
        format with the independent variable in the first column.

        :arg independent: a string indentifying the independent variable's filename and the
          appropriate property of the data in the file. For example: "concs.in|depth" would
          use the depth of the data in the file 'concs.in' as its value for each test case.
        :arg dependent: a list of strings that use the same format as the string for the
          'independent' variable.
        :arg threshold: a float value specifying the minimum level of success that the output
          file must attain before it can be used in the plot.
        :arg headings: a list of strings to identify each of the variables; there should be
          one for each of the dependent variables plus another for the independent variable.
          The order of the headings should be [independent, dependents...] where the dependents
          are ordered the same way as they appear in 'dependents'.
        :arg functions: a dictionary of functions to apply to the values of each variable. E.g.
          {"group.in|depth", "log"} would apply the numpy.log function to each value extracted
          from the depth property of the data in 'group.in' before tabulating it. If the value
          is not a string, it must be callable with some argument.
        """
        from fortpy.testgen import print_value
        xs, ys = self._get_data_series(independent, dependents, threshold, functions)

        data = []
        header = []
        for xset, yset in zip(xs, ys):
            x, xlabel = xset
            y, ylabel = yset
            data.append(x)
            header.append(self._format_heading(xlabel))
            data.append(y)
            header.append(self._format_heading(ylabel))

        printed = []
        if headings is not None and len(headings) > 0:
            header = [self._format_heading(h) for h in headings]
        printed.append(' | '.join(header))
        printed.append('-'*(17*len(header)))

        i = 0
        done = False
        while not done:
            done = True
            row = []
            for col in range(len(data)):
                if i < len(data[col]):
                    row.append(str(print_value(data[col][i])))
                    done = False
                else:
                    row.append(' '*17)
            i += 1
            printed.append('   '.join(row))
                
        return '\n'.join(printed)

    def _format_fit(self, key):
        """Formats the function fit at the specified key for printing/plotting."""
        if key in self.fits:
            params = self.fits[key]["params"]
            model = self.fits[key]["model"]
            if any(params < 1e-2):
                mdict = {
                    "exp": "{0:.2e}exp({1:.2e}*x){2:+.2e}",
                    "lin": "{0:.2e}*x{1:+.2e}"
                }
            else:
                mdict = {
                    "exp": "{0:.2f}exp({1:.2f}*x){2:+.2f}",
                    "lin": "{0:.2f}*x{1:+.2f}"
                }
            return mdict[model].format(*params)
        else:
            msg.warn("Couldn't format the fit {}; fit not found.".format(key))
            return key

    def fit(self, independent, dependent, model, threshold=1., functions=None):
        """Tries to fit the data selection for the specified dependent variable against the
        indepedent one using the specified model. Arguments are similar to self.table(). Returns
        a tuple of (fitting parameters, 1 standard deviation errors).

        :arg model: a description of the function to fit to. Possible values are ['exp', 'lin'].
        """
        key = "{}${}".format(independent, dependent)
        xs, ys = self._get_data_series(independent, [dependent], threshold, functions)
        if xs[0][0] is None:
            return

        if len(xs[0][0]) != len(ys[0][0]):
            msg.err("Can't fit data unless we have the same number of values for the independent "
                    "and dependent variables. " 
                    "len({}) == {}; len({}) == {}".format(independent, len(xs[0][0]), 
                                                          dependent, len(ys[0][0])))
            return

        from scipy.optimize import curve_fit
        from numpy import sqrt, diag
        mdict = {
            "exp": _fit_exp,
            "lin": _fit_lin
        }
        if not model in mdict:
            msg.err("Cannot fit using model '{}'; unknown function.".format(model))
            return

        func = mdict[model]
        popt, pcov = curve_fit(func, xs[0][0], ys[0][0])
        perr = sqrt(diag(pcov))
        
        self.fits[key] = {
            "params": popt,
            "covarmat": pcov,
            "function": (lambda x: func(x, *popt)),
            "model": model,
            "error": perr
        }

    def _analyze_data(self, fullpath, name, fullparse):
        """Analyzes the data in the specified file to determine its shape and some other
        properties. Returns a dictionary with the properties.
        """
        from os import path
        import re
        result = {}
        rxcomment = re.compile("^\s*#")
        width = 0
        depth = 0
        values = []
        ragged = False

        with open(path.join(fullpath, name)) as f:
            for line in f:
                if not rxcomment.match(line):
                    if fullparse or len(values) == 0:
                        values.append(list(map(eval, line.split())))
                    if width == 0:
                        width = len(values[-1])
                    elif width != len(values[-1]):
                        ragged = True
                    depth += 1

        if len(values) > 1 or len(values[0]) > 1:
            from numpy import array
            values = array(values)
            result["rowvals"] = values
            if not ragged:
                result["colvals"] = array([values[:,i] for i in range(len(values[0]))])
        elif len(values) == 1 and len(values[0]) == 1:
            result["value"] = values[0][0]

        result["depth"] = depth
        result["width"] = width
        result["shape"] = (depth, width)
        return result

    def _process_test(self, fullpath, files, fullparse):
        """Processes the contents of the specified test directory to get statistics
        about input/output and result files.
        """
        from os import path
        result = {}
        #The fortpy files we are interested in are 'fpytiming.out' and 'SUCCESS'.
        if "fpytiming.out" in files:
            with open(path.join(fullpath, "fpytiming.out")) as f:
                abstime = float(f.readlines()[1])
            result["timing"] = abstime

        if "SUCCESS" in files:
            from dateutil import parser
            with open(path.join(fullpath, "SUCCESS")) as f:
                lastrun = parser.parse(f.readlines()[0])
            result["lastrun"] = lastrun

        #The rest of the files should have either "in" or "out" in them. We just
        #analyze the contents to determine the shape of the files. Later we could
        #add some statistics about mean and variance of the date etc.
        inputs = {}
        outputs = {}
        percents = {}
        for name in files:
            if name == "SUCCESS" or name == "fpytiming.out":
                continue
            if ".in" in name:
                analysis = self._analyze_data(fullpath, name, fullparse)
                if name not in self.infiles:
                    self.infiles.append(name)
                if name not in self.props:
                    self.props[name] = list(analysis.keys())
                inputs[name] = analysis
            elif ".out.compare" in name:
                percent = None
                with open(path.join(fullpath, name)) as f:
                    i = 0
                    for line in f:
                        if i == 3:
                            percent = float(line.split()[0].split("%")[0])/100.
                            break
                        i += 1
                percents[name] = percent
                short = name.replace(".compare", "")
                if short not in self.comparable:
                    self.comparable.append(short)                    
            elif ".out" in name:
                analysis = self._analyze_data(fullpath, name, fullparse)
                if name not in self.outfiles:
                    self.outfiles.append(name)
                if name not in self.props:
                    self.props[name] = list(analysis.keys())
                outputs[name] = analysis

        result["inputs"] = inputs
        result["outputs"] = outputs
        result["percents"] = percents
        return result

    def parse(self, fullparse):
        """Parses all the test result in the specified staging directory.

        :arg fullparse: when true the entire contents of input/output files is also parsed
          and interpreted so that it is available for plotting along with other results. This
          increases the parsing time significantly.
        """
        #We need to look at each of the sub-directories in 'stagedir/tests'.
        from os import path, walk
        testpath = path.join(self.stagedir, "tests")
        result = {}
        for root, dirs, files in walk(testpath):
            if root.lower() != testpath.lower():
                result[root.split("/")[-1]] = self._process_test(root, files, fullparse)
        self.details = result
