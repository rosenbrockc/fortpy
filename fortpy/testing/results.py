class ACResult(object):
    """The result of an auto-class folder comparison across all nested
    variables recursively.

    :arg compared: a list of FileComparer results for files that were
      actually compared because both files existed.
    :arg onlym: files that were only in the model folder.
    :arg onlyx: files that were only in the test execution folder.
    :arg tolerance: the maximum error allowed for the any individual 
      file comparison result.
    """
    def __init__(self, compared, onlym, onlyx, tolerance):
        self.compared = compared
        self.onlym = onlym
        self.onlyx = onlyx
        self.tolerance = tolerance
        
        self.failures = []
        """A list of all the file comparisons that failed or had similarities
        below the threshold.
        """

        #These are lazy versions of the properties so that the variables
        #they depend on don't need to be stored indefinitely.
        self._has_data = None
        self._percent_match = None
        self._common_match = None

    def clean(self):
        """Dereferences pointers to test results and keeps only the aggregate
        properties (percent_match, common_match, has_data).
        """
        #Make sure that all the properties are calculated.
        self.percent_match + self.common_match
        self.compared = None
        self.onlym = None
        self.onlyx = None
        
    @property
    def has_data(self):
        """Returns True if this result had actual data to compare (compared
        to two empty sets of data).
        """
        if self._has_data is None:
            self._has_data = len(self.compared) > 0 or len(self.onlym) > 0 or len(self.onlyx) > 0
        return self._has_data
    
    @property
    def percent_match(self):
        """Returns a value indicating how similar the two variables are."""
        if self._percent_match is None:
            total = sum([c.percent_match for c in self.compared if c.percent_match is not None])
            if self.has_data:
                self._percent_match = total/(len(self.compared) + len(self.onlym) + len(self.onlyx))
            else:
                self._percent_match = 0
        return self._percent_match

    @property
    def common_match(self):
        """Returns a value indicating how similar the overlapping entries
        of the two dictionaries are."""
        if self._common_match is None:
            total = 0.
            for cres in self.compared:
                if cres is None or cres.common_match < self.tolerance:
                    self.failures.append(cres)
                elif cres is not None:
                    total += cres.common_match

            if self.has_data:
                self._common_match = total/(len(self.compared) + len(self.onlym) + len(self.onlyx))
            else:
                self._common_match = 0
        return self._common_match
            
class DictResult(object):
    """The result of a dictionary comparison between two dicts.

    :attr common: the number of keys and values that matched.
    :attr different: a list of the values where the keys matched
      but the values didn't.
    :attr only1: the keys and values that were only in the first.
    :attr only2: the keys and values that were only in the second.
    """
    def __init__(self, dict1, dict2, label, outcomes):
        self.common = 0
        self.different = []
        self.only1 = {}
        self.only2 = {}
        self.label = label
        self._matched = []

        #Keep a pointer to the original dictionaries for easier reference
        self.dict1 = dict1
        self.dict2 = dict2

        #Also keep a pointer to the outcomes so we know how to calculate
        #the match percent.
        self.outcomes = outcomes

        #These are lazy versions of the properties so that the variables
        #they depend on don't need to be stored indefinitely.
        self._has_data = None
        self._percent_match = None
        self._common_match = None

    def __str__(self):
        if self.dict1 is None:
            return "Dereferenced fortpy.testing.results.DictResult<{}>".format(id(self))
        else:
            return print_dict_result(self)

    def clean(self):
        """Dereferences pointers to test results and keeps only the aggregate
        properties (percent_match, common_match, has_data).
        """
        #Make sure that all the properties are calculated.
        self.percent_match + self.common_match
        self.different = None
        self.only1 = None
        self.only2 = None
        self._matched = None
        self.dict1 = None
        self.dict2 = None
    
    @property
    def has_data(self):
        """Returns True if this result had actual data to compare (compared
        to two empty sets of data).
        """
        if self._has_data is None:
            self._has_data = not (len(self.dict1) == 0 and len(self.dict2) == 0)
        return self._has_data        
    
    def add_common(self, key):
        """Adds the specified key to the list of keys AND values that
        matched between the dictionaries."""
        self._matched.append(key)
        self.common += 1

    def _ignore_total(self, aniter):
        """Counts the number of non-ignorable items in the iterable."""
        total = 0
        for k in aniter:
            if not self.outcomes.can_ignore(self.label, k):
                total += 1
        return total

    def _diff_total(self):
        """Determines how many of the differences between common elements
        matter in light of the outcomes."""
        total = 0
        for t in self.different:
            if not self.outcomes.can_ignore(self.label, t[0]):
                total += 1
        return total

    @property
    def percent_match(self):
        """Returns a value indicating how similar the two dictionaries are."""
        if self._percent_match is None:
            if self.outcomes.can_ignore(self.label):
                self._percent_match = 1
            else:
                #This is the total number of dictionary entries that can *not* be ignored.
                intersect = set(self.dict1.keys()).intersection(set(self.dict2.keys()))
                total = self._ignore_total(intersect)
                #This is how many of the keys that actually had the same values cannot be ignored.
                common = self._ignore_total(self._matched)
                if total == 0 and common == 0:
                    #Everything was told to be ignored, so this is a perfect match!
                    self._percent_match = 1
                else:
                    self._percent_match = float(common) / total
        return self._percent_match
                
    @property
    def common_match(self):
        """Returns a value indicating how similar the overlapping entries
        of the two dictionaries are."""
        if self._common_match is None:
            if self.outcomes.can_ignore(self.label):
                self._common_match = 1
            else:
                common = self._ignore_total(self._matched)
                diff = self._ignore_total([ k[0] for k in self.different ])
                if common + diff != 0:
                    self._common_match = float(common) / (common + diff)
                else:
                    #There was nothing to compare, return 1
                    self._common_match = 1
        return self._common_match
    
class ListResult(object):
    """The result of comparing two lists."""
    def __init__(self, list1, list2, label=None, outcomes=None):
        self.common = 0
        self.different = []
        self.label = label
        self.outcomes = outcomes
        self.sharedkey = False

        #Keep a pointer to orig references.
        self.list1 = list1
        self.list2 = list2

        #These are lazy versions of the properties so that the variables
        #they depend on don't need to be stored indefinitely.
        self._has_data = None
        self._percent_match = None
        self._common_match = None

    def __str__(self):
        if self.list1 is None:
            return "Dereferenced fortpy.testing.results.ListResult<{}>".format(id(self))
        else:
            return print_list_result(self)

    def clean(self):
        """Dereferences pointers to test results and keeps only the aggregate
        properties (percent_match, common_match, has_data).
        """
        #Make sure that all the properties are calculated.
        self.percent_match + self.common_match
        self.different = None
        self.list1 = None
        self.list2 = None
    
    @property
    def has_data(self):
        """Returns True if this result had actual data to compare (compared
        to two empty sets of data).
        """
        if self._has_data is None:
            self._has_data = not (len(self.list1) == 0 and len(self.list2) == 0)
        return self._has_data

    @property
    def percent_match(self):
        """Returns a value indicating how similar the two lists are."""
        if self._percent_match is None:
            if self.outcomes is not None and self.outcomes.can_ignore(self.label):
                self._percent_match = 1
            else:
                total = len(self.list1) + len(self.list2)
                if total > 0:
                    self._percent_match = float(2 * self.common) / total
                else:
                    #If both the lists had no elements, it is a perfect match.
                    self._percent_match = 1
        return self._percent_match

    @property
    def common_match(self):
        """Returns the same value as the percent match since lists
        don't have a definition of "common"."""
        if self._common_match is None:
            if not self.sharedkey:
                self._common_match = self.percent_match
            else:
                self._common_match = 1
        return self._common_match

class CompatibilityResult(object):
    """The result of a compatibility comparison between two blocks."""
    def __init__(self, label, outcomes):
        self.common = 0
        self.key_errors = []
        self.total = 1
        self.label = label

        #Also keep a pointer to the outcomes so we know how to calculate
        #the match percent.
        self.outcomes = outcomes

        #These are lazy versions of the properties so that the variables
        #they depend on don't need to be stored indefinitely.
        self._common_match = None

    def clean(self):
        """Dereferences pointers to test results and keeps only the aggregate
        properties (percent_match, common_match, has_data).
        """
        #Make sure that all the properties are calculated.
        self.percent_match + self.common_match
        self.key_errors = None

    def __str__(self):
        if self.key_errors is None:
            return "Dereferenced fortpy.testing.results.CompatibilityResult<{}>".format(id(self))
        else:
            return print_compat_result(self)

    #TODO, the compatibility results don't have very good granularity because
    #it is a pain to implement. If it is ever needed, we can come back.

    @property
    def percent_match(self):
        return float(2 * self.common) / self.total

    @property
    def common_match(self):
        """Returns a value indicating how similar the overlapping entries
        of the two blocks are."""
        if self._common_match is None:
            self._common_match = float(2 * self.common) / (self.total - len(self.key_errors))
        return self._common_match

def accumulate_matches(results, outcomes = None, ispercent = True):
    """Calculates the average of all the entries matches in the dict."""
    total = 0
    count = 0
    for rkey in results:
        if outcomes is not None and type(rkey) == type(""):
            #We need to count the results only if they shouldn't be ignored.
            if "." in rkey:
                nkey = rkey.split(".")[0]
            else:
                nkey = rkey
            
            if not outcomes.can_ignore(nkey):
                r = results[rkey]
            else:
                r = None
        else:
            r = results[rkey]

        if r is not None:
            if ispercent:
                total += r.percent_match
            else:
                total += r.common_match
            count += 1

    if count != 0:
        return total / count
    elif total == 0:
        #If both are zero, the body is empty for both cases and it is a perfect
        #match (two empty sets).
        return 1
    else:
        return 0

class BlockResult(object):
    """Represents the result of comparing two body blocks."""
    def __init__(self, outcomes, key = None, index = None, template = None):
        self.key = key
        self.index = index
        self.results = {}

        if template is not None and template.key is not None:
            self.key = "{} == {}".format(template.key, key)

        #Also keep a pointer to the outcomes so we know how to calculate
        #the match percent.
        self.outcomes = outcomes

        #These are lazy versions of the properties so that the variables
        #they depend on don't need to be stored indefinitely.
        self._has_results = None
        self._percent_match = None
        self._common_match = None

    def clean(self):
        """Dereferences pointers to test results and keeps only the aggregate
        properties (percent_match, common_match, has_data).
        """
        #Make sure that all the properties are calculated.
        self.percent_match + self.common_match
        bool(self.has_results)
        self.results = None

    def __str__(self):
        if self.results is None:
            return "Dereferenced fortpy.testing.results.BlockResult<{}>".format(id(self))
        else:
            return print_block_result(self, False)

    @property
    def has_data(self):
        """Returns True if this result had actual data to compare (compared
        to two empty sets of data).
        """
        return self.has_results
    
    @property 
    def has_results(self):
        """Returns True if this block result has any child results."""
        if self._has_results is None:
            self._has_results = len(list(self.results.keys())) > 0
        return self._has_results

    @property
    def percent_match(self):
        """Calculates the combined percent match for the whole
        block from its constituent results."""
        if self._percent_match is None:
            self._percent_match = accumulate_matches(self.results, self.outcomes)
        return self._percent_match

    @property
    def common_match(self):
        """Calculates the combined match for the block by taking
        only those overlapping keys into account from named and compatibility
        comparisons."""
        if self._common_match is None:
            self._common_match = accumulate_matches(self.results, self.outcomes, False)
        return self._common_match

class BodyResult(object):
    """The result of comparing all the body entries using the template."""
    def __init__(self, rep1, rep2):
        self.blocks = {}
        self.only1 = {}
        self.only2 = {}
        self.rep1 = rep1
        self.rep2 = rep2
        """Pointer to the parent file representation that has template information
        and the contents of the files in python variables.
        """

        #These are lazy versions of the properties so that the variables
        #they depend on don't need to be stored indefinitely.
        self._has_data = None
        self._percent_match = None
        self._common_match = None

    def clean(self):
        """Dereferences pointers to test results and keeps only the aggregate
        properties (percent_match, common_match, has_data).
        """
        #Make sure that all the properties are calculated.
        for block in self.blocks.values():
            block.clean()
            
        self.percent_match + self.common_match
        self.only1 = None
        self.only2 = None
        self.rep1 = None
        self.rep2 = None

    def __str__(self):
        if self.rep1 is None:
            return "Dereferenced fortpy.testing.results.BodyResult<{}>".format(id(self))
        else:
            return print_body_result(self)

    @property
    def has_data(self):
        """Returns True if this result had actual data to compare (compared
        to two empty sets of data).
        """
        if self._has_data is None:
            self._has_data = not (len(self.blocks) == 0 and len(self.only1) == 0 and len(self.only2) == 0)
        return self._has_data
        
    @property
    def percent_match(self):
        """Calculates the combined percent match for the whole
        body from its constituent results."""
        if self._percent_match is None:
            #We need to reduce the match based on how many blocks actually
            #overlapped in the body.
            overlap = len(list(self.blocks.keys()))
            total = overlap + len(list(self.only1.keys())) + len(list(self.only2.keys()))
            if total == 0: #Since both have zero entries, we have a perfect match
                self._percent_match = 1
            else:
                self._percent_match = accumulate_matches(self.blocks) * float(overlap) / total
        return self._percent_match
    
    @property
    def common_match(self):
        """Calculates the combined match for the body by taking
        only those overlapping blocks into account from named and compatibility
        comparisons."""
        if self._common_match is None:
            self._common_match = accumulate_matches(self.blocks, False)
        return self._common_match

class CompareResult(object):
    """The results from comparing two files using a file template."""
    def __init__(self, preamble, body, file1, file2):
        self.preamble = preamble
        self.body = body
        self.file1 = file1
        self.file2 = file2

    def clean(self):
        """Dereferences pointers to test results and keeps only the aggregate
        properties (percent_match, common_match, has_data).
        """
        self.preamble.clean()
        self.body.clean()

    def __str__(self):
        if self.body.rep1 is None:
            return "Dereferenced fortpy.testing.results.CompareResult<{}>".format(id(self))
        else:
            return print_compare_result(self)

    @property
    def has_data(self):
        """Returns True if this result had actual data to compare (compared
        to two empty sets of data).
        """
        return self.preamble.has_results or self.body.has_data

    @property
    def percent_match(self):
        """Calculates the combined percent match for the whole
        file from its constituent results."""
        if self.preamble.has_results and self.body.has_data:
            return (self.preamble.percent_match + self.body.percent_match) / 2
        elif self.preamble.has_results:
            return self.preamble.percent_match
        else:
            return self.body.percent_match

    @property
    def common_match(self):
        """Calculates the combined match for the file by taking
        only those overlapping keys into account from named and compatibility
        comparisons."""
        if self.preamble.has_results and self.body.has_data:
            return (self.preamble.common_match + self.body.common_match) / 2
        elif self.preamble.has_results:
            return self.preamble.common_match
        else:
            return self.body.common_match


def print_compare_result(result, verbose = False):
    """Presents the results of a file comparison operation for
    printing to console or saving to file."""
    if isinstance(result, CompareResult):
        return _print_r_header(result) + _print_body(result, verbose)
    else:
        return print_list_result(result)

def _print_r_header(result):
    """Returns global information from the comparison."""
    lines = []
    lines.append("\n----------------------------------------------------------------------")
    lines.append("FORTPY File Comparison Report")
    lines.append(print_matches(result))
    lines.append("----------------------------------------------------------------------")
    lines.append("")

    lines.append("SOURCE: {}".format(result.file1))
    lines.append("TARGET: {}".format(result.file2))
    lines.append("\n")

    return "\n".join(lines)

def _print_body(result, verbose = False):
    lines = []
    lines.append("----------------------------------------------------------------------")
    lines.append("FILE PREAMBLE COMPARISON DETAILS")
    if result.preamble.has_results:
        lines.append(print_matches(result.preamble))
    else:
        lines.append(" No Preamble Defined in Template")
    lines.append("----------------------------------------------------------------------\n")

    if result.preamble.has_results:
        lines.append(print_block_result(result.preamble, False, verbose))
    lines.append(print_body_result(result.body, verbose))

    return "\n".join(lines)
        

def print_body_result(result, verbose = False):
    """Returns a string representation of the block result."""
    lines = []
    lines.append("\n----------------------------------------------------------------------")
    lines.append("FILE BODY COMPARISON DETAILS")
    lines.append(print_matches(result))
    lines.append("----------------------------------------------------------------------\n")

    for block in sorted(result.blocks.keys()):
        if result.blocks[block].common_match < 1:
            lines.append(print_block_result(result.blocks[block], True, verbose))

    lines.extend(print_body_onlys(result.only1, result.rep1, 1, verbose))
    lines.append('\n')
    lines.extend(print_body_onlys(result.only2, result.rep2, 2, verbose))
            
    return "\n".join(lines)

def print_body_onlys(only, rep, i, verbose=False):
    """Returns a string representation of those entries in the body that were in only one
    or the other, but not both.
    """
    lines = []
    if len(only) > 0:
        lines.append("ENTRIES ONLY IN FILE {}\n".format(i))
        for key, block in only.items():
            for lineid, linevals in block.items():
                for lineval in linevals:
                    for vkey in sorted(lineval.named.keys()):
                        compkey = "{}.{}".format(lineid, vkey)
                        if isinstance(rep.template.key, list):
                            add = compkey in rep.template.key
                        else:
                            add = compkey == rep.template.key
                        if add:
                            lines.append("{}: {}".format(vkey, str(lineval.named[vkey])))
                    lines.append("")
    return lines
                        
def print_block_result(result, separator = True, verbose=False):    
    """Returns a string representation of the block result."""
    lines = []
    if result.key is not None:
        lines.append("BLOCK: key '{}'".format(result.key))
    elif result.index is not None:
        lines.append("BLOCK: index {}".format(result.index))
    else:
        lines.append("BLOCK")        
        
    lines.append(print_matches(result))
    lines.append("")

    for d in sorted(result.results.keys()):
        c = result.results[d]
        ignore = result.outcomes.can_ignore(d)
        if isinstance(c, DictResult):
            lines.append(print_dict_result(c, d, verbose, ignore))
        elif isinstance(c, ListResult):
            lines.append(print_list_result(c, d, verbose, ignore))
        elif isinstance(c, CompatibilityResult):
            lines.append(print_compat_result(c, d, verbose, ignore))
        
    if separator:
        lines.append("\n--------------------------------------------------")

    return "\n".join(lines)

def print_list_result(result, label = "", verbose = False, ignored = False):
    """Returns a string representation of the list result."""
    lines = []
    signore = "(I) " if ignored else ""

    if len(result.list1) == 1 and len(result.list2) == 1:
        if isinstance(result.list1[0], float):
            #Had to add this so that the compare report showed enough significant
            #figures for the human to understand the difference.
            fmtstr = "MISMATCH ({0:.12f} vs. {1:.12f})"
        else:
            fmtstr = "MISMATCH ({} vs. {})"
        match = fmtstr.format(result.list1[0], result.list2[0]) if result.common == 0 else "MATCH"
        lines.append("{}ITEM\t{} - {}".format(signore, label, match))            
    else:
        lines.append("{}LIST\t{} : {}".format(signore, label, print_matches(result)))
        if verbose:
            lines.append("\t{}/{} matches, {} items vs. {} items".format(
                result.common, result.common + len(result.different),
                len(result.list1), len(result.list2)))
            lines.append("\tDIFF: " + result.different.__str__())

    return "\n".join(lines)

def print_dict_result(result, label = "", verbose = False, ignored = False):
    """Returns a string representation of the dict result."""
    lines = []
    signore = "(I) " if ignored else ""

    lines.append("{}DICT\t{} - {}".format(signore, label, print_matches(result)))
    if verbose:
        lines.append("\t{}/{} matches, {} keys vs. {} keys".format(
            result.common, result.common + len(result.different),
            len(list(result.dict1.keys())), len(list(result.dict2.keys()))))
        lines.append("\tDIFF: " + result.different.__str__())
    return "\n".join(lines)
    
def print_compat_result(result, label = "", verbose = False, ignored = False):    
    """Returns a string representation of the compatibility result."""
    lines = []
    signore = "(I) " if ignored else ""

    lines.append("{}COMPAT\t{} - {}".format(signore, label, print_matches(result)))
    if verbose:
        lines.append("\t{}/{} matches/key errors".format(result.common,
                                                         len(result.key_errors)))
        
    return "\n".join(lines)

def print_matches(result):
    """Prints the specified result in decimal percent format as a string."""
    if result.has_data:
        return "{0:.2%} ({1:.2%} RAW)".format(result.common_match, result.percent_match)
    else:
        return "Excluded from Comparison by Template"
