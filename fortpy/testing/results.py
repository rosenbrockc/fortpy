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

    def __str__(self):
        return print_dict_result(self)

    def add_common(self, key):
        """Adds the specified key to the list of keys AND values that
        matched betweend the dictionaries."""
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

    @property
    def percent_match(self):
        """Returns a value indicating how similar the two dictionaries are."""
        if self.outcomes.can_ignore(self.label):
            return 1
        else:
            total = self._ignore_total(self.dict1) + self._ignore_total(self.dict2)
            common = self._ignore_total(self._matched)
            if total == 0 and common == 0:
                return 1
            else:
                return float(2 * common) / total

    @property
    def common_match(self):
        """Returns a value indicating how similar the overlapping entries
        of the two dictionaries are."""
        if self.outcomes.can_ignore(self.label):
            return 1
        else:
            common = self._ignore_total(self._matched)
            diff = self._ignore_total([ k[0] for k in self.different ])
            if common + diff != 0:
                return float(common) / (common + diff)
            else:
                #There was nothing to compare, return 1
                return 1
    
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

    def __str__(self):
        return print_list_result(self)

    @property
    def percent_match(self):
        """Returns a value indicating how similar the two lists are."""
        if self.outcomes is not None and self.outcomes.can_ignore(self.label):
            return 1
        else:
            total = len(self.list1) + len(self.list2)
            if total > 0:
                return float(2 * self.common) / total
            else:
                #If both the lists had no elements, it is a perfect match.
                return 1

    @property
    def common_match(self):
        """Returns the same value as the percent match since lists
        don't have a definition of "common"."""
        if not self.sharedkey:
            return self.percent_match
        else:
            return 1

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

    def __str__(self):
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
        return float(2 * self.common) / (self.total - len(self.key_errors))

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

    def __str__(self):
        return print_block_result(self, False)

    @property 
    def has_results(self):
        """Returns True if this block result has any child results."""
        return len(list(self.results.keys())) > 0

    @property
    def percent_match(self):
        """Calculates the combined percent match for the whole
        block from its constituent results."""
        return accumulate_matches(self.results, self.outcomes)

    @property
    def common_match(self):
        """Calculates the combined match for the block by taking
        only those overlapping keys into account from named and compatibility
        comparisons."""
        return accumulate_matches(self.results, self.outcomes, False)

class BodyResult(object):
    """The result of comparing all the body entries using the template."""
    def __init__(self, body1, body2):
        self.blocks = {}
        self.only1 = {}
        self.only2 = {}

    def __str__(self):
        return print_body_result(self)

    @property
    def percent_match(self):
        """Calculates the combined percent match for the whole
        body from its constituent results."""
        #We need to reduce the match based on how many blocks actually
        #overlapped in the body.
        overlap = len(list(self.blocks.keys()))
        total = overlap + len(list(self.only1.keys())) + len(list(self.only2.keys()))
        return accumulate_matches(self.blocks) * float(overlap) / total

    @property
    def common_match(self):
        """Calculates the combined match for the body by taking
        only those overlapping blocks into account from named and compatibility
        comparisons."""
        return accumulate_matches(self.blocks, False)

class CompareResult(object):
    """The results from comparing two files using a file template."""
    def __init__(self, preamble, body, file1, file2):
        self.preamble = preamble
        self.body = body
        self.file1 = file1
        self.file2 = file2

    def __str__(self):
        return print_compare_result(self)

    @property
    def percent_match(self):
        """Calculates the combined percent match for the whole
        file from its constituent results."""
        if self.preamble.has_results:
            return (self.preamble.percent_match + self.body.percent_match) / 2
        else:
            return self.body.percent_match

    @property
    def common_match(self):
        """Calculates the combined match for the file by taking
        only those overlapping keys into account from named and compatibility
        comparisons."""
        if self.preamble.has_results:
            return (self.preamble.common_match + self.body.common_match) / 2
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

    return "\n".join(lines)

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
        match = "MISMATCH ({} vs. {})".format(result.list1[0], result.list2[0]) if result.common == 0 else "MATCH"
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
    return "{0:.2%} ({1:.2%} RAW)".format(result.common_match, result.percent_match)
