import os
from .. import msg
import re
import xml.etree.ElementTree as ET
from fortpy.testing.comparer import FileComparer

class TemplateLine(object):
    """Represents a single line in the template file and how to format it.

    :arg element: the XML element that defines this line in the file.
    :arg group: the [group] that this line belongs to.
    """
    def __init__(self, element, group, commentchar):
        self.identifier = element.attrib["id"]

        #First look at the "mandatory" attributes and assign defaults if missing
        if "type" in element.attrib:
            self.dtype = re.split(",\s*", element.attrib["type"])
        else:
            msg.warn("no type specified for {}. Assuming string.".format(self.identifier))
            self.dtype = [ "string" ]

        #Values specifies how many variable values are present in the file
        if "values" in element.attrib:
            self.values = re.split(",\s*", element.attrib["values"])
            i = 0
            for i in range(len(self.values)):
                if self.values[i].isdigit():
                    self.values[i] = int(self.values[i])
        elif "from" not in element.attrib:
            msg.warn("no value count specified for {}. Assuming *.".format(self.identifier))
            self.values = [ "*" ]
        else:
            self.values = []

        #Handle all the optional attributes
        if "choices" in element.attrib:
            self.choices = re.split(",\s*", element.attrib["choices"])
        else:
            self.choices = []
        if "comment" in element.attrib:
            self.comment = "{} {}".format(commentchar, element.attrib["comment"])
        else:
            self.comment = ""
        if "default" in element.attrib:
            self.default = eval(element.attrib["default"])
        else:
            self.default = None
        #How from works: if an element has a from attribute, it is included in
        #the plaintext file after conversion but does *not* appear in the XML
        #file that is being converted. It grabs its value from another group
        #or line whose id is the from attribute's value.
        if "from" in element.attrib:
            self.fromtag = element.attrib["from"]
        else:
            self.fromtag = None
        #Related to from, this operator specifies how the value should be generated
        #using the line/group whose id is the from attribute's value.
        if "operator" in element.attrib:
            self.operator = element.attrib["operator"]
        else:
            self.operator = "count"

        self.group = group
        self._nvalues = None
        self._caster = {
            "int": self._cast_int,
            "float": self._cast_float,
            #We want to use the same syntax even though we do nothing with strings
            "string": lambda s: s 
        }

    @property
    def nvalues(self):
        """Returns the number of values recorded on this single line. If the
        number is variable, it returns -1."""
        if self._nvalues is None:
            self._nvalues = 0
            for val in self.values:
                if type(val) == type(int):
                    self._nvalues += val
                else:
                    self._nvalues = -1
                    break

        return self._nvalues

    def write(self, valuedict):
        """Returns the lines that this template line should add to the input file."""
        if self.identifier in valuedict:
            value = valuedict[self.identifier]
        elif self.default is not None:
            value = self.default
        elif self.fromtag is not None and self.fromtag in valuedict:
            if self.operator == "count":
                value = len(valuedict[self.fromtag])
            else:
                msg.err("referenced 'from' attribute/operator {} not in xml dictionary.".format(self.fromtag))
                exit(1)
        else:
            msg.err("a required line {} had no value or default specified.".format(self.identifier))
            exit(1)            

        #Before we generate the result, validate the choices if they exist
        if len(self.choices) > 0:
            for single in value:
                if str(single) not in self.choices:
                    msg.warn("failed choices validation for {} in {} (line {})".format(
                        single, self.choices, self.identifier))

        result = []
        #Get the string representation of the value
        if isinstance(value, list):
            sval = " ".join([ str(val) for val in value])
        else:
            sval = str(value)

        if self.comment != "" and (self.nvalues < 0 or self.nvalues > 5):
            #We will put the comments on a separate line from the actual values.
            result.append(self.comment)
            result.append(sval)
        else:
            result.append("{} {}".format(sval, self.comment))
                    
        return result

    def parse(self, element):
        """Parses the contents of the specified XML element using template info.

        :arg element: the XML element from the input file being converted.
        """
        result = []
        if element.text is not None and element.tag == self.identifier:
            l, k = (0, 0)
            raw = element.text.split()
            while k < len(self.values):
                dtype = self.dtype[k]
                if isinstance(self.values[k], int):
                    for i in range(self.values[k]):
                        result.append(self._caster[dtype](raw[i + l]))
                    l += self.values[k]
                    k += 1
                else:
                    #This is a variable argument line, just use up the rest
                    #of them as the type of the current line
                    rest = [ self._caster[dtype](val) for val in raw[l::] ]
                    result.extend(rest)
                    break
        else:
            msg.warn("no results for parsing {} using line {}".format(element.tag, self.identifier))

        return result

    def _cast_int(self, value):
        """Returns the specified value as int if possible."""
        try:
            return int(value)
        except ValueError:
            msg.err("Cannot convert {} to int for line {}.".format(value, self.identifier))
            exit(1)

    def _cast_float(self, value):
        """Returns the specified value as float if possible."""
        try:
            return float(value)
        except ValueError:
            msg.err("Cannot convert {} to float for line {}.".format(value, self.identifier))
            exit(1)

class TemplateGroup(object):
    """Represents a logical grouping of line templates.

    :arg element: the XML group element to parse.
    :arg commentchar: the character(s) that specify comment lines. Used when
      inserting comments beside lines in the plaintext file.
    """
    def __init__(self, element, commentchar):
        self.identifier = element.attrib["name"]
        self.order = []
        self.lines = {}

        if "comment" in element.attrib:
            self.comment = "{} {}".format(commentchar, element.attrib["comment"])
        else:
            self.comment = ""
        if "repeat" in element.attrib:
            self.repeat = element.attrib["repeat"]
        else:
            self.repeat = None

        self._load(element, commentchar)
        
    def _load(self, element, commentchar):
        """Loads all the child line elements from the XML group element."""
        for child in element:            
            if "id" in child.attrib:
                tline = TemplateLine(child, self, commentchar)
                self.order.append(tline.identifier)
                self.lines[tline.identifier] = tline
            else:
                msg.warn("no id element in {}. Ignored. (group._load)".format(child))

    def parse(self, element):
        """Extracts the values from the specified XML element that is being converted."""
        #All the children of this element are what we are trying to parse.
        result = []
        for child in element:
            if child.tag in self.lines:
                values = { child.tag: self.lines[child.tag].parse(child) }
                result.append(values)

        return result

    def write(self, valuedict):
        """Generates the lines for the converted input file using the specified
        value dictionary."""
        result = []
        if self.identifier in valuedict:
            values = valuedict[self.identifier]
        else:
            return result

        if self.comment != "":
            result.append(self.comment)

        if self.repeat is not None and type(values) == type([]):
            if self.repeat.isdigit():
                for i in range(int(self.repeat)):  
                    result.extend(self._write_iterate(values[i]))
            else:
                #We are repeating for as many values as we have in the value
                #entry for the group in the dictionary.
                for value in values:
                    result.extend(self._write_iterate(value))
        elif type(values) == type({}):
            #This group doesn't get repeated, so the values variable must
            #be a dictionary, just run it once.
            result = self._write_iterate(values)

        return result

    def _write_iterate(self, values):
        """Generates the lines for a single pass through the group."""
        result = []
        for key in self.order:
            result.append(self.lines[key].write(values))

        if len(result) > 1:
            return result
        else:
            return result[0]
        
class TemplateContents(object):
    """The contents of an XML input template.

    :attr order: a list of id attributes from the lines in the template file
      that preserves the order in which the lines showed up in the file.
    :attr entries: a dictionary of TemplateLine and TemplateGroup instances
      for the corresponding lines and groups in the template. Dict keys are
      the identifiers in the order list.
    :attr comment: the character(s) at the start of a line that specify it as
      a comment line."""
    def __init__(self):
        self.order = []
        self.entries = {}
        self.comment = "#"

class FileTemplate(object):
    """Represents an XML template that specifies how to format an input/output
    file using a dictionary of keyed values.

    :arg path: the full path to the XML template file to load.
    """
    def __init__(self, path, name, direction="input"):
        self.name = name
        self.path = os.path.expanduser(path)
        self.versions = {}
        self.direction = direction
        self._load()
        
    def _load(self):
        """Extracts the XML template data from the file."""
        if os.path.exists(self.path):
            root = ET.parse(self.path).getroot()
            if (root.tag == "fortpy" and "mode" in root.attrib and
                root.attrib["mode"] == "template" and "direction" in root.attrib and
                root.attrib["direction"] == self.direction):
                #First, we need instances of the template contents for each of the
                #versions listed in the fortpy tag.
                for v in _get_xml_version(root):
                    self.versions[v] = TemplateContents()
                #Now we can update the contents objects using the XML data.
                self._load_entries(root)

                #See if a custom name was specified for the auto-converted
                #files.
                if "autoname" in root.attrib:
                    self.name = root.attrib["autoname"]
            else:
                msg.err("the specified template {} ".format(self.path) + 
                        "is missing the mode and direction attributes.")
                exit(1)
        else:
            msg.err("could not find the template {}.".format(self.path))
            exit(1)

    def parse(self, root):
        """Returns a dictionary of values extracted from the root of the
        specified XML file. It is assumed that the file is an input/output
        file to be converted into plaintext. As such the file should only
        specify a single version number."""
        #Use the first element in the versions list since there should only be one.
        v = _get_xml_version(root)[0]
        result = {}

        for child in root:         
            if child.tag in self.versions[v].entries:
                entry = self.versions[v].entries[child.tag]
                #Entry can be either a line or a group. Both objects have a parse
                #method that returns a list of values. In the line's case, the
                #list is the values from that line. For the group, it is a list
                #of dictionaries, a dictionary for each tag name.
                result[child.tag] = entry.parse(child)

        return result

    def write(self, valuedict, version):
        """Generates the lines for the converted input file from the valuedict.

        :arg valuedict: a dictionary of values where the keys are ids in the
          template and the values obey their template rules.
        :arg version: the target version of the output file.
        """
        result = []

        if version in self.versions:
            for tag in self.versions[version].order:
                entry = self.versions[version].entries[tag]
                result.extend(entry.write(valuedict))

        return result

    def _load_entries(self, root):
        """Loads all the child entries of the input template from the 
        specified root element."""
        mdict = {
            "comments": self._comment,
            "line": self._line,
            "group": self._group
        }
        for entry in root:
            mdict[entry.tag](entry)

    def _comment(self, element):
        """Extracts the character to use for comments in the input file."""
        for v in _get_xml_version(element):
            self.versions[v].comment = element.text

    def _line(self, element):
        """Parses the XML element as a single line entry in the input file."""
        for v in _get_xml_version(element):
            if "id" in element.attrib:
                tline = TemplateLine(element, None, self.versions[v].comment)
                self.versions[v].entries[tline.identifier] = tline
                self.versions[v].order.append(tline.identifier)
            else:
                msg.warn("no id element in {}. Ignored. (_line)".format(element))
    
    def _group(self, element):
        """Parses the XML element as a group of [unknown] number of lines."""
        for v in _get_xml_version(element):
            if "name" in element.attrib:
                g = TemplateGroup(element, self.versions[v].comment)
                self.versions[v].entries[g.identifier] = g
                self.versions[v].order.append(g.identifier)
            else:
                msg.warn("no name element in {}. Ignored. (_group)".format(element))            
        
def _get_xml_version(element):
    """Extracts a list of versions that an xml element references. Returns
    a [ 1 ] list if there isn't a versions attribute."""
    if "versions" in element.attrib:
        result = [ int(v) for v in re.split(",\s*", element.attrib["versions"]) ]
    else:
        result = [ 1 ]

    return result

class FileConverter(object):
    """Converts XML-based input/output files into non-keyword based ones.
    
    :arg template_dir: the path to the directory containing input file templates.
    """
    def __init__(self, template_dir):
        self.template_dir = os.path.expanduser(template_dir)
        self.templates = {}

    def convert(self, path, version, target = None):
        """Converts the specified file using the relevant template.

        :arg path: the full path to the file to convert.
        :arg version: the new version of the file.
        :arg target: the optional path to save the file under. If not
          specified, the file is saved based on the template file name.
        """
        #Get the template and values out of the XML input file and
        #write them in the format of the keywordless file.
        values, template = self.parse(path)
        lines = template.write(values, version)

        #Finally, write the lines to the correct path.
        if target is None:
            target = os.path.join(os.path.dirname(path), template.name)

        with open(os.path.expanduser(target), 'w') as f:
            f.write("\n".join(lines))

    def parse(self, path):
        """Extracts a dictionary of values from the XML file at the specified path."""
        #Load the template that will be used for parsing the values.
        expath, template, root = self._load_template(path)
        if expath is not None:
            values = template.parse(root)

        return (values, template)
        
class OutputConverter(object):
    """Converts plain-text output files between versions."""
    def __init__(self, template_dir):
        self.comparer = FileComparer(os.path.expanduser(template_dir))

    def convert(self, path, version, target):
        """Converts the specified source file to a new version number."""
        source = self.comparer.get_representation(path)
        lines = [ '# <fortpy version="{}"></fortpy>\n'.format(version) ]

        for line in self.comparer.template.contents[version].preamble:
            lines.append(line.write(source.preamble, source.version, source.stored) + "\n")

        for line in self.comparer.template.contents[version].body:
            for valueset in source.body:
                lines.append(line.write(valueset, source.version, source.stored) + "\n")

        with open(os.path.expanduser(target), 'w') as f:
            f.writelines(lines)

class InputConverter(FileConverter):
    """Converts XML-based input files into non-keyword based ones.
    
    :arg template_dir: the path to the directory containing input file templates.
    """
    def __init__(self, template_dir):
        super(InputConverter, self).__init__(template_dir)

    def _load_template(self, path):
        #First we extract the file name for the template or look for it
        #in the root element. The naming convention is to use .xin.xml
        #as the extension. If we replace the .xin.xml by .in.xml it 
        #should cover most cases.
        expath = os.path.expanduser(path)
        root = ET.parse(expath).getroot()
        if root.tag == "fortpy" and "mode" in root.attrib and \
           root.attrib["mode"] == "input":
            #This is a valid input file.
            if "template" in root.attrib:
                template = root.attrib["template"]
            else:
                template = os.path.split(expath)[1].replace(".xin.xml", ".in.xml")

            tpath = os.path.join(self.template_dir, template)
            name = template.replace(".xml","")
            self.templates[template] = FileTemplate(tpath, name)
            return (expath, self.templates[template], root)
        else:
            msg.warn("the input file {} is missing the mode attribute.".format(path))
            return None
