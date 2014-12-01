"""Module for analyzing the results of multiple test cases against each other
for a single unit test.
"""
from fortpy import msg
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
                results.append("   {0}: '{1}' => {2:.2f}%".format(j+1, fkeys[j], 
                                                                  data[keys[i]][fkeys[j]]*100))
            printed.append('\n'.join(results) + '\n')

        return '\n'.join(printed)

    def _case_filter(self, caseid, tfilter):
        """Determines whether the specified filter matches the case.

        :arg caseid: the identifier for the specific test case.
        :arg tfilter: the case filter; if '*' is present in the filter, then the case
          matches true as long as it contains the filter text; otherwise the filter needs
          to match exactly.
        """
        if tfilter is None:
            return True
        else:
            from fnmatch import fnmatch
            return (("*" in tfilter and fnmatch(caseid, tfilter)) or
                    (not "*" in tfilter and tfilter == caseid))

    def _get_data(self, variable, order=None, threshold=1., tfilter=None, functions=None):
        """Returns a list of the valid data points for the specified variable.

        :arg variable: a string indentifying the variable's filename and the
          appropriate property of the data in the file. For example: "concs.in|depth" would
          use the depth of the data in the file 'concs.in' as its value for each test case.
        :arg threshold: a float value specifying the minimum level of success that the output
          file must attain before it can be used in the plot.
        """
        from numpy import ndarray
        target = None
        if variable != "timing":
            varfile, attribute = variable.split("|")
        else:
            target = varfile = "timing"

        if target is None:
            if varfile in self.infiles:
                target = "inputs"
            elif varfile in self.outfiles:
                target = "outputs"
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

        #Handle cases where the values need to have functions applied to them before tabulating
        #or printing the values.
        if functions is not None and variable in functions:
            if isinstance(functions[variable], str) or isinstance(functions[variable], unicode):
                from importlib import import_module
                module = functions[variable].split(".")
                fname = module.pop()
                numpy = import_module('.'.join(module))
                fx = getattr(numpy, fname)
            else:
                fx = functions[variable]
        else:
            fx = None

        for testcase in details:
            if not self._case_filter(testcase, tfilter):
                #We don't consider test cases that are blocked by the filter.
                continue

            value = 0
            if target == "timing" and self._is_successful(self.details[testcase]):
                value = self.details[testcase][target]
            else:
                if (self._is_successful(self.details[testcase], varfile, threshold) and
                    varfile in self.details[testcase][target]):
                    data = self.details[testcase][target][varfile]
                    if attribute in data:
                        value = data[attribute]

            #We need the arrays to be the same length for plotting; if we are following
            #an existing independent variable, we still need to append zero to the list
            #for the values we don't have.
            if isinstance(value, ndarray) or value != 0 or order is not None:
                if fx is not None:
                    values.append(fx(value))
                else:
                    #Handle the case where we are plotting an entire row as a single point
                    #and an aggregrate function should have been specified
                    if isinstance(value, list):
                        msg.err("Can't coerce an array to a single value without an aggregation "
                                "function such as numpy.mean or numpy.sum")
                        return ([0], "Array Aggregation Error")
                    else:
                        values.append(value)
                cases.append(testcase)

        return (values, cases)

    def _get_data_series(self, independent, dependents, threshold, tfilter, functions=None):
        """Used by plot() and table() to get the x and y data series that will be plotted
        or tabulated against each other. See the arguments on those methods.
        """
        if (("rowvals" in independent or "colvals" in independent) or
            any([("rowvals" in d or "colvals" in d) for d in dependents])):
            #We need to make sure that the test case filter they specified returns only a
            #single result so that we get a reasonable plot.
            if tfilter is None or sum([(1 if self._case_filter(d, tfilter) else 0)
                                       for d in self.details]) > 1:
                msg.err("Can't plot/tabulate array-valued data for more than one test case at a time.")
                return (None, None)

        x, cases = self._get_data(independent, None, threshold, tfilter)
        if (len(x) == 1 and isinstance(x[0], list) and 
            ("rowvals" in independent or "colvals" in independent)):
            x = x[0]
        ys = []
        for variable in dependents:
            ypts, names = self._get_data(variable, cases, threshold, tfilter, functions)
            if (len(ypts) == 1 and isinstance(ypts[0], list) and
                ("rowvals" in variable or "colvals" in variable)):
                ypts = ypts[0]
            ys.append((ypts, variable))

        return (x, ys)

    def plot(self, independent, dependents, threshold=1., xlabel=None, ylabel=None,
             tfilter=None, savefile=None, functions=None, xscale=None, yscale=None):
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
        :arg tfilter: specifies a string that must be present in each of the test case names in
          order for the test case to be considered.
        :arg functions: a dictionary of functions to apply to the values of each variable. E.g.
          {"group.in|depth": "log"} would apply the numpy.log function to each value extracted
          from the depth property of the data in 'group.in' before plotting it. If the value
          is not a string, it must be callable with some argument.
        """
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from numpy import linspace
        from itertools import cycle
        x, ys = self._get_data_series(independent, dependents, threshold, tfilter)
        if x is None:
            return        

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if xscale is not None:
            ax.set_xscale(xscale)
        if yscale is not None:
            ax.set_yscale(yscale)

        colors = cm.rainbow(linspace(0, 1, len(ys)))
        cycols = cycle(colors)

        for y, label in ys:
            try:
                ax.scatter(x, y, s=10, color=next(cycols), label=label)
            except ValueError:
                msg.err("The values for {} can't be log-plotted.".format(label))
        plt.legend(loc='upper left');
        if savefile is None:
            plt.show()
        else:
            plt.savefig(savefile)

    def _format_heading(self, heading):
        """Formats the specified heading text to conform to the table layout."""
        if len(heading) > 15:
            return "%15s" % heading[0:15]
        else:
            return "%15s" % heading

    def table(self, independent, dependents, threshold=1., headings=None,
              tfilter=None, functions=None):
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
        :arg tfilter: specifies a string that must be present in each of the test case names in
          order for the test case to be considered.
        :arg functions: a dictionary of functions to apply to the values of each variable. E.g.
          {"group.in|depth", "log"} would apply the numpy.log function to each value extracted
          from the depth property of the data in 'group.in' before tabulating it. If the value
          is not a string, it must be callable with some argument.
        """
        from fortpy.testgen import write_generic
        x, ys = self._get_data_series(independent, dependents, threshold, tfilter, functions)
        if x is None:
            return

        data = [x]
        header = [self._format_heading(independent)]
        for y, label in ys:
            data.append(y)
            header.append(self._format_heading(label))

        from numpy import array
        printed = [write_generic(zip(*data))]
        if headings is not None and len(headings) > 0:
            header = [self._format_heading(h) for h in headings]
        printed.insert(0, ' | '.join(header))
        printed.insert(1, ''.join(['-']*(17*len(header))))
        return '\n'.join(printed)

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
                        values.append(map(eval, line.split()))
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
