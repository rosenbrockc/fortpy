from .. import msg
from .results import DictResult, ListResult
import xml.etree.ElementTree as ET
import re

class TemplateContents(object):
    """Represents the contents of a file template for a single version.

    :attr comments: the character that designates the rest of the line
      as a comment in the file being templated.
    :attr stop: possible values are EOF|finite. Specifies how the body
      template will be iterated over.
    :attr comparisons: a dictionary with mode types as keys and compare
      rules as values (type templates.FieldComparisons).
    :attr key: as the body is iterated, values are added to a list in
      the order that they appear. If a key is specified, the value of
      the field key will be used to determine which body blocks represent
      the same data.
    """
    def __init__(self):
        self.comments = ""
        self.stop = "EOF"
        self.comparisons = {}
        self.outcomes = {}
        self.key = None

        #The actual specification in the template for how many times to
        #repeat the body block read. Can be a number or a variable name.
        self._count = None
        #The parsed count value. If it was a number, this will be a valid
        #python number type; otherwise a variable name.
        self._pcount = None

        #The order in which the line identifiers will be encountered in the
        #actual output file for each body block. List of line identifiers.
        self._bodyorder = []
        #The dictionary of template lines (type FileLine) that define the
        #lines that make up a single body block.
        self._body = {}
        #As for bodyorder, but for the lines in the template preamble
        self._preambleorder = []
        #As for body, but defines a single preamble block whose stored variable
        #values are available to all the body blocks.
        self._preamble = {}

    def add_body(self, fline):
        """Adds the specified file line to the body of the template.

        :arg fline: a FileLine object representing the next line that
          will be encountered when the file body is parsed.
        """
        self._bodyorder.append(fline.identifier)
        self._body[fline.identifier] = fline

    def add_preamble(self, fline):
        """Adds the specified file line to the preamble of the template.

        :arg fline: a FileLine object representing the next line that
          will be encountered when the file body is parsed.
        """
        self._preambleorder.append(fline.identifier)
        self._preamble[fline.identifier] = fline
        
    @property
    def body(self):
        """Returns the body FileLines in order of appearance."""
        return [ self._body[n] for n in self._bodyorder ]

    @property
    def preamble(self):
        """Returns the preamble FileLines in order of appearance."""
        return [ self._preamble[n] for n in self._preambleorder ]

    @property
    def count(self):
        """Returns the number of lines that this template should be used for."""
        if self._pcount is None:
            if self._count is not None:
                if self._count.isdigit():
                    self._pcount = int(self._count)
                else:
                    self._pcount = self._count
            else:
                self._pcount = 1

        return self._pcount

class LineValues(object):
    """Represents the results of parsing a line using a FileLine.parse().

    :attr named: a dictionary of values extracted from a line that were
      named according to the list in the 'names' attribute. Key is variable
      name, value is the value.
    :attr stored: a dictionary of values that were specified for storing so
      that other lines in the same block/file can access their values.
    :attr values: the raw values extracted from the line and cast to the
      data types specified by the 'type' attribute.
    """
    def __init__(self):
        self.named = {}
        self.stored = {}
        self.values = []

class LineGroup(object):
    """Represents a logical grouping of <line> entries in an output template
    that should be repeated a *variable* number of times.
    """
    def __init__(self, xml):
        self.identifier = "default" if "name" not in xml.attrib else xml.attrib["name"]
        if "repeat" not in xml.attrib:
            raise ValueError("'repeat' is a required attribute of <group> in an output template.")
        else:
            #We have to do some fancy things here with stored/named values etc.
            self.count = xml.attrib["repeat"]
        self.comment = None if "comment" not in xml.attrib else xml.attrib["comment"]

        self.lines = []
        for child in xml:
            if (child.tag == "line" or child.tag == "lines") and "id" in child.attrib:
                fline = FileLine(child)
                self.lines.append(fline)

        self.line_counts = None
        """List of integer line counts, one for each FileLine in self.lines.
        """
        self._line_cum = None
        """Cumulative number of times to use each FileLine. Makes it easier to
        pick which one we should be using.
        """

    def update_counts(self, lcounts, gcount):
        """Updates the list of counts describing how often each FileLine should be repeated.

        :arg counts: a list of integer values, one for each FileLine in self.lines.
        """
        self.line_counts = counts
        self._line_cum = [sum(a[0:i]) for i in range(1, len(self.lines))]

    def parse(self, line, i):
        """Parses the specified line using the relevant FileLine object, based on the global line
        counter 'i'.
        """
        #i is zero-based. However, once we reach sum(self.line_counts), we need to repeat
        #the line templates again. igroup is the index within the group (instead of global).
        igroup = i % sum(self.line_counts)
        iline = [0 if c < i else 1 for c in self._line_cum].index(1)
        return self.lines[iline].parse(line)
                
class FileLine(object):
    """Represents a template line specification in a file.

    :arg element: the XML 'line' tag element.
    :attr compatibility: a dictionary that maps variable names in one version
      to different names in another version of the template.
    """
    def __init__(self, element, parent):
        self.xml = element
        self.compatibility = {}
        self.defaults = {}
        self.parent = parent
        """The FileTemplate that this FileLine belongs to."""

        #Overwrite makes the default values get used *even if* a value was
        #specified in the dictionary for write mode.
        self.overwrite = False

        #The names that should be used for the parsed values of
        #a particular line in the file being compared.
        self._names = None
        self._raw_names = None

        #The names that should be used for 'global' variables whose values
        #will be available to the entire block/file.
        self._stores = None

        #The actual specification in the template for how many times to
        #repeat the line read. Can be a number or a variable name.
        self._count = None
        #The parsed count value. If it was a number, this will be a valid
        #python number type; otherwise a variable name.
        self._pcount = None

        self._load_xml()

    @property
    def count(self):
        """Returns the number of lines that this template should be used for."""
        if self._pcount is None:
            if self._count is not None:
                if self._count.isdigit():
                    self._pcount = int(self._count)
                else:
                    self._pcount = self._count
            else:
                self._pcount = 1

        return self._pcount

    @property
    def unique_names(self):
        """Returns a list of all the named variables where each variable only
        appears once, even if it is multi-valued.
        """
        return [n.split("*")[0] for n in self._raw_names]

    def write(self, valuedict, version, stored):
        """Creates a string representation for this line template using the
        specified values as part of output file conversion.

        :arg valuedict: the dictionary of values from the version being 
          converted.
        :arg version: the version number of the values from the version
          being converted.
        """
        result = []
        count = self.count
        if type(count) == type("") and count in stored:
            try:
                count = int(stored[count])
            except ValueError:
                msg.warn("Can't understand how to use {} for count".format(count))
                return
        else:
            count = 1

        if self.identifier in valuedict and not self.overwrite:
            values = valuedict[self.identifier]

            for i in range(count):
                outvals = []
                if self._raw_names is None:
                    #There aren't any named variables, so we just write the
                    #values directly to the line.
                    outvals.append(self._write_values(values[i].values))
                else:
                    outvals.extend(self._write_values_generic(values[i].named, version))

                result.append(" ".join(outvals))

        elif self.identifier in self.defaults:
            #We are going to use defaults. If there need to be multiple entries
            #use the same default value for all of them
            if type(self.defaults[self.identifier]) == type({}):
                value = " ".join(self._write_values_generic(self.defaults[self.identifier],
                                                           version))                    
            else:
                value = self.defaults[self.identifier]

            for i in range(count):
                result.append(value)
    
        return "\n".join(result)

    def _write_values_generic(self, values, version):
        """Creates a list of elements to write for this line using a generic
        dict of values to work from."""
        result = []
        for name in self._raw_names:
            sname = name.split("*")[0]
            value = self._write_find_values(sname, values)
            if value is None and version in self.compatibility:
                value = self._write_compat_values(sname, version, values)
            if value is not None:
                result.append(value)

        return result

    def _write_compat_values(self, name, version, values):
        """Returns a string representing the values obtained when compatibility
        is taken into account between versions.
        """
        usename = None
        for oldname in self.compatibility[version]:
            if self.compatibility[version][oldname] == name:
                usename = oldname
                break

        if usename is not None:
            return self._write_find_values(usename, values)

    def _write_find_values(self, name, values):
        """Searches for the value to use for the specified variable; first looks
        in 'values', then in defaults for this line.
        """
        if name in values:
            if hasattr(values[name], "values"):
                return self._write_values(values[name].values)
            else:
                return self._write_values(values[name])
        elif name in self.defaults:
            return self._write_values(self.defaults[name])
        elif (self.identifier in self.defaults and 
              name in self.defaults[self.identifier]):
            return self._write_values(self.defaults[self.identifier][name])
        else:
            return None
        
    def _write_values(self, values):
        """Returns a string representing the specified values."""
        if type(values) == type([]):
            return " ".join([str(v) for v in values])
        else:
            return str(values)

    def _load_xml(self):
        """Examines XML element to extract file line info."""        
        #We can handle multiple lines with the same class and template info.
        self.multiple = self.xml.tag == "lines"
        if self.multiple and "count" in self.xml.attrib:
            self._count = self.xml.attrib["count"]

        #Get the mandatory attributes first.
        self.identifier = self.xml.attrib["id"]
        self.dtypes = re.split(",\s*", self.xml.attrib["type"])
        self.values = re.split(",\s*", self.xml.attrib["values"])

        #Handle default value specifiers for the output conversion capability
        if "default" in self.xml.attrib:
            defaults = re.split(",\s*", self.xml.attrib["default"])
            innerdict = {}
            for d in defaults:
                if "=" in d:
                    name, value = d.split("=")
                    innerdict[name] = value

            if len(list(innerdict.keys())) == 0:
                self.defaults[self.identifier] = d
            else:
                self.defaults[self.identifier] = innerdict

        #See which of the optional attribs are in the element
        if "overwrite" in self.xml.attrib:
            self.overwrite = self.xml.attrib["overwrite"] == "true"
        if "store" in self.xml.attrib:
             self._stores = re.split(";\s*", self.xml.attrib["store"])
        if "names" in self.xml.attrib:
            #The template allows them to repeat names using a *[int] notation
            #If the same name appears multiple times, the values are grouped
            #into a single liste under that name when values are extracted.
            self._names = []
            self._raw_names = re.split(",\s*", self.xml.attrib["names"])

            for n in self._raw_names:
                if "*" in n:
                    name, times = n.split("*")
                    for t in range(int(times)):
                        self._names.append(name)
                else:
                    self._names.append(n)

        #The line(s) element may have some children for compatibility
        kids = list(self.xml)
        if len(kids) > 0:
            for kid in kids:
                if kid.tag == "compatibility":
                    self._load_compat_xml(kid)

    def _load_compat_xml(self, element):
        """Extracts XML data from a compatibility tag in the line element."""
        for vtag in element:
            #Each child of compatibility is a version element that describes
            #mappings between version names of values.
            versions = xml_get_versions(vtag)
            for v in versions:
                if not v in self.compatibility:
                    self.compatibility[v] = {}
                #Add all the mappings from this version tag to the list.
                mappings = {}
                for amap in re.split(",\s*", vtag.attrib["mappings"]):
                    source, target = amap.split("=")
                    mappings[source.strip()] = target.strip()
                self.compatibility[v][vtag.attrib["id"]] = mappings

    def parse(self, line):
        """Parses a line from an actual file using the rules in
        this line's template definition."""
        #Initialize the result of this parsing operation.
        result = LineValues()
        #First, we split on whitespace to get all the elements in the line
        raw = line.strip().split()

        #Loop over the total number of known values and cast them to the right type
        k = 0
        for i in range(len(self.values)):
            #If there are a variable number of entries for this value
            #just use up all the rest
            if self.values[i] == "*":
                loop = list(range(k, len(raw)))
                namek = k
            else:
                loop = list(range(int(self.values[i])))
            #If there are a variable number of entries, they would be stored together
            #as a list under a certain name. Use this as a clearing house. After we have
            #stored the values, we can extend the results list.
            current = []

            for j in loop:
                if k >= len(raw):
                    if self.parent is not None:
                        emsg = "Specified known value index '{}/{}' exceeds line value count ({}). Using template '{}'."
                        msg.err(emsg.format(k, len(loop)-1, raw, self.parent.filepath))
                    else:
                        msg.err("Specified known value index '{}/{}' exceeds line value count ({}).".format(k, len(loop)-1, raw))
                val = raw[k]
                dtype = self.dtypes[i]
                try:
                    if dtype == "int":
                        current.append(int(val))
                    elif dtype == "float":
                        current.append(float(val))
                    else:
                        current.append(val)
                except ValueError:
                    msg.err("[{}] could not parse value '{}' of type '{}'.\n".format(
                        self.identifier, val, dtype))

                #If names were specified for the values, we need to populate the dict
                #now
                if self._names is not None:
                    if self.values[i] != "*":
                        if self._names[k] not in result.named:
                            result.named[self._names[k]] = current[j]
                        else:
                            if type(result.named[self._names[k]]) == type([]):
                                result.named[self._names[k]].append(current[j])
                            else:
                                result.named[self._names[k]] = [ result.named[self._names[k]], 
                                                                 current[j] ] 
                k += 1
            
            #Now that the look is over, if we were naming variables, we want 
            #to save the rest of the current
            #values list under the name.
            result.values = current

            if self.values[i] == "*":
                if self._names is not None:
                    result.named[self._names[namek]] = current            
                #We used up all the values, save the user from themselves
                break

        #Now that we have determined all the values, we can store the ones
        #that need to be used later. We import operator so that the eval()
        #can work properly.
        import operator
        if self._stores is not None:
            for s in self._stores:
                name, target = s.split("=")
                if "$" in target:
                    store = eval(target.replace("$", "result.values"))
                elif re.match("\d+", target) is not None:
                    store = eval(target)
                result.stored[name] = store

        return result

class LineComparer(object):
    """Compares values for specific names between dictionaries
    using specified operators and tolerances.

    :attr name: the name of the value in the line that will be compared.
    :attr element: the xml element that defined the comparison.
    """
    def __init__(self, name, element):
        self.name = name
        self.numeric = False

        if "operator" in element.attrib:
            self.operator = element.attrib["operator"]
        else:
            self.operator = "equals"
        if "tolerance" in element.attrib:
            self.tolerance = element.attrib["tolerance"]
            if self.tolerance[0].isdigit():
                #We are working with a number, just eval() it and use it in
                #a finite difference comparison.
                try:
                    self.tolerance = eval(self.tolerance)
                    self.numeric = True
                except ValueError:
                    msg.warn("tolerance for comparison {} ".format(element.attrib["id"]) + 
                             "should be a number but can't be evaluated.")
                    self.tolerance = None
        else:
            self.tolerance = None

        self._comparer = {
            "equals": self._compare_equals,
            "finite": self._compare_finite
        }

    @property
    def isdict(self):
        """Determines whether this comparer compares only dictionaries."""
        return self.name is not None

    def compare(self, value1, value2):
        """Compares a value in the two dicts/values according to the settings
        in this line comparer. Returns True if they match within the
        specified tolerance."""
        #If anything doesn't match up, we just say they don't match.
        result = False
        if self.isdict:
            if self.name in value1 and self.name in value2:
                result = self._comparer[self.operator](value1, value2, self.isdict)
        #We can't process regular values with finite differences unless
        #we have numeric tolerances.
        elif self.operator == "equals" or self.numeric:
            result = self._comparer[self.operator](value1, value2, self.isdict)
        
        return result

    def _compare_equals(self, value1, value2, isdict = True):
        """Determines if the two values are equal."""
        if isdict:
            return value1[self.name] == value2[self.name]
        else:
            return value1 == value2

    def _compare_finite(self, value1, value2, isdict = True):
        """Determines if the two values are equal within the tolerance."""
        #First we need to check if the tolerance is a number or a reference
        #to a variable in the dictionary.
        if self.tolerance is not None:
            if self.numeric:
                if isdict:
                    return value1[self.name] - value2[self.name] <= self.tolerance
                else:
                    return value1 - value2 <= self.tolerance
            else:
                #Look for the values in the dictionaries and build a dynamic
                #tolerance value. This won't be reached unless isdict==True
                try:
                    s1 = eval(self.tolerance.replace("$", "value1"))
                    return value1[self.name] - value2[self.name] <= s1
                except ValueError:
                    msg.warn("could not generate dynamic tolerance for comparison" +
                             "{} and tolerance {}""".format(self.name, self.tolerance))
                    return False
        #We can't perform a finite difference calculation unless a tolerance
        #was specified.
        else:
            return self._compare_equals(value1, value2, isdict)

class FieldComparisons(object):
    """Represents instructions on how to compare fields that have
    random variance."""
    def __init__(self, element = None):
        self.compares = {}
        if element is not None:
            self._load_xml(element)

    def _load_xml(self, element):
        """Extracts all of the child field comparisons from a
        comparisons tag in the template."""
        for child in element:
            self._load_compare_xml(child)
        
    def _load_compare_xml(self, element):
        """Extracts comparison information for a single compare entry 
        in a comparison."""
        if "id" in element.attrib:
            identifier = element.attrib["id"]
            if "." in identifier:
                line, name = identifier.split(".")
                if not line in self.compares:
                    self.compares[line] = {}
                self.compares[line][name] = LineComparer(name, element)
            else:
                self.compares[identifier] = LineComparer(None, element)
                
    def compare_d(self, dict1, dict2, key, outcomes):
        """Compares all values in the two dictionaries using any
        comparison rules for the specific field specified in the
        template.
        
        :arg key: the identifier of the line that these lists are
          associated with.
        :arg outcoms: a TemplateOutcomes with information on how to
          intepret comparison results.
        """
        #Initialize a list result. The basic initialization is common to
        #both logic branches.
        result = DictResult(dict1, dict2, key, outcomes)

        if key in self.compares:
            self._compare_dict(dict1, dict2, self.compares[key], result)
        else:
            self._compare_dict(dict1, dict2, None, result)

        return result

    def _compare_dict(self, dict1, dict2, line_comparer, result):
        """Compares the values of two dictionaries by key."""
        #First compare the keys in the first with the second.
        for key in dict1:
            if key not in dict2:
                result.only1[key] = dict1[key]
            else:
                #Use a flag because of the complicated logic tree
                compared = False

                if line_comparer is not None:
                    if type(line_comparer) == type({}) and key in line_comparer:
                        ctrue = line_comparer[key].compare(dict1, dict2)
                        compared = True
                    elif isinstance(line_comparer, LineComparer):
                        ctrue = line_comparer.compare(dict1[key], dict2[key])
                        compared = True

                if not compared:
                    ctrue = dict1[key] == dict2[key]

                if ctrue:
                    result.add_common(key)
                else:
                    result.different.append((key, dict1[key], dict2[key]))

        #Now, see if the second has anything not in the first
        for key in dict2:
            if key not in dict1:
                result.only2[key] = dict2[key]

    def compare_l(self, list1, list2, key, outcomes):
        """Compares the values at corresponding indices in the lists
        using any comparison rules defined.

        :arg key: the identifier of the line that these lists are
          associated with.
        """
        #Initialize a list result. The basic initialization is common to
        #both logic branches.
        result = ListResult(list1, list2, key, outcomes)
        elcount = min([len(list1), len(list2)])

        if key in self.compares:
            #The key cannot have any variable names, otherwise it would be
            #a dictionary
            self._compare_list(list1, list2, self.compares[key], result)
        else:
            #We only do equality comparison on each element.
            self._compare_list(list1, list2, None, result)

        return result

    def _compare_list(self, list1, list2, line_comparer, result):
        """Performs the element-wise list comparison using the specified
        comparer and appending the outcomes to the result."""        
        elcount = min([len(list1), len(list2)])
        for i in range(elcount):
            if line_comparer is not None:
                ctrue = line_comparer.compare(list1[i], list2[i])
            else:
                if isinstance(list1[i], float):
                    #Override the default behavior of finite precision comparisions for
                    #float values to be the precision default in fortpy.
                    ctrue = (list1[i]-list2[i] < 1e-13)
                else:
                    ctrue = list1[i] == list2[i]
                    
            if not ctrue:
                result.different.append((i, list1[i], list2[i]))
            else:
                result.common += 1

class TemplateOutcomes(object):
    """Represents a set of rules to use in determining outcomes of 
    comparisons between files.

    :arg element: the xml 'outcomes' element that this class handles.
    :attr ignore: a dictionary of data elements that should be ignored
      in comparison interpretation. Keys are line identifiers in the
      preamble or body block. Values are lists of variable names in
      the line that should be ignored.
    """
    def __init__(self, element = None):
        self.ignore = {}

        #The original list of line.name ignore directives.
        self._ignore = []

        if element is not None:
            self._load_xml(element)

    def can_ignore(self, line, name = None):
        """Determines whether the line and/or name should be ignored."""
        if name is None:
            return line in self.ignore and self.ignore[line] is None
        else:
            if line in self.ignore:
                return self.ignore[line] is None or name in self.ignore[line]
            else:
                return False

    def _load_xml(self, element):
        """Extracts all relevant tags from the parent 'outcome' element."""
        for child in element:
            if child.tag == "ignore" and "id" in child.attrib:
                self._ignore.append(child.attrib["id"])

        self._parse_ignores()

    def _parse_ignores(self):
        """Extracts line id and line names as lists in dictionaries
        to make ignore checking faster."""
        for i in self._ignore:
            if "." in i:
                line, name = i.split(".")
                if line in self.ignore:
                    if self.ignore[line] is None:
                        self.ignore[line] = [ name ]
                    else:
                        self.ignore[line].append(name)
                else:
                    self.ignore[line] = [ name ]
            else:
                self.ignore[i] = None

class FileTemplate(object):
    """Represents an XML template defining multiple versions
    of the same file for comparison.

    :arg filepath: the path to the XML template file to load.
    :attr contents: a dictionary with key=version# and value 
      a TemplateContents() for that specific version.
    """
    def __init__(self, filepath):
        self.filepath = filepath
        self.contents = {}
        self._xml_load()

    def _xml_load(self):
        """Loads the XML file and splits the template entries
        based on version numbers."""
        with open(self.filepath) as f:
            lines = f.read()
        from fortpy.utility import XML_fromstring
        root = XML_fromstring(lines, self.filepath)

        #The first element must be a fortpy with attribute template
        #otherwise give a message about it being an invalid template file
        if "mode" in root.attrib and root.attrib["mode"] == "template":
            #See if we have multiple versions to work with or if it is a
            #straight laced template
            versions = xml_get_versions(root)
            #Create a contents object for each version specified in the file.
            for v in versions:
                if v not in self.contents:
                    self.contents[v] = TemplateContents()                    
            self._xml_load_versions(root)
        else:
            msg.warn("The XML template {} is not".format(self.filepath) + 
                     " a valid fortpy template file.")

    def _xml_load_versions(self, root):
        """Loads the template from XML tracking important version information."""
        #Creating a dictionary like this is a slick way to handle multiple cases.
        methods = {
            "preamble": self._xml_v_lineparent,
            "body": self._xml_v_body,
            "comments": self._xml_v_comments,
            "comparisons": self._xml_v_comparisons,
            "outcomes": self._xml_v_outcomes
            }
        #Run the extraction method for the relevant tag
        for child in root:
            methods[child.tag](child)

        #There has to be a FieldComparer to process file comparisons. If one wasn't
        #explicitly specified, create a default one that does regular comparisons.
        defaultfc = FieldComparisons()
        for vkey in self.contents:
            if "default" not in self.contents[vkey].comparisons:
                self.contents[vkey].comparisons["default"] = defaultfc

        #THere also has to be a TemplateOutcomes.
        defaultto = TemplateOutcomes()
        for vkey in self.contents:
            if "default" not in self.contents[vkey].outcomes:
                self.contents[vkey].outcomes["default"] = defaultto

    def _xml_v_outcomes(self, element):
        """Extracts outcome information from the template."""
        versions = xml_get_versions(element)
        modes = xml_get_modes(element)

        to = TemplateOutcomes(element)
        for v in versions:
            for m in modes:
                self.contents[v].outcomes[m] = to

    def _xml_v_comparisons(self, element):
        """Extracts comparison information for specific fields from a comparison element."""
        versions = xml_get_versions(element)
        modes = xml_get_modes(element)

        fc = FieldComparisons(element)
        for v in versions:
            for m in modes:
                self.contents[v].comparisons[m] = fc

    def _xml_v_body(self, element):
        """Extracts the body attributes and child lines."""
        #We are looking for a specification on how to handle reading the body in
        versions = xml_get_versions(element)
        for v in versions:
            if "stop" in element.attrib:
                self.contents[v].stop = element.attrib["stop"]
            if "count" in element.attrib:
                self.contents[v]._count = element.attrib["count"]
            if "key" in element.attrib:
                self.contents[v].key = element.attrib["key"]
                if "," in self.contents[v].key:
                    self.contents[v].key = re.split(",\s*", self.contents[v].key)
            
        self._xml_v_lineparent(element)
            
    def _xml_v_comments(self, element):
        """Extracts the comments information from the specified comments element."""
        versions = xml_get_versions(element)
        for v in versions:
            self.contents[v].comments = element.text

    def _xml_v_lineparent(self, element):
        """Extracts the line-type elements from the specified parent element."""
        for child in element:
            versions = xml_get_versions(child)
            if (child.tag == "line" or child.tag == "lines") and "id" in child.attrib:
                fline = FileLine(child, self)
            elif child.tag == "group":
                fline = LineGroup(child)
            else:
                msg.warn("non line-type tag in <{0}>\n{1}\n</{0}>".format(element.tag, element.text))

            for v in versions:
                if element.tag == "preamble":
                    self.contents[v].add_preamble(fline)
                elif element.tag == "body":
                    self.contents[v].add_body(fline)
        
def xml_get_modes(element):
    """Returns a list of comparison modes declared in the XML element."""
    if "mode" in element.attrib:
        return re.split(",\s*", element.attrib["mode"])
    else:
        return [ "default" ]

def xml_get_versions(element):
    """Returns a list of versions referenced in the XML element."""
    if "versions" in element.attrib:
        return [ int(n.strip()) for n in  element.attrib["versions"].split(",") ]
    else:
        #There are not multiple versions, so this is the first!
        return [ 1 ]                
