import re

class DocGroup(object):
    """Represents a list of DocElements that have been logically grouped together."""
    
    def __init__(self, XMLElement, decorates = None):
        self.xml = XMLElement
        self.decorates = decorates
        self.attributes = self.xml.attrib
        self.doctype = "group"

    def __str__(self):
        return "GROUP: {}\nAttributes: {}\nDecorates: {}\n".format(self.name,
                                                                self.attributes, self.decorates)
    @property
    def name(self):
        """Gets the name of this group if it exists."""
        if "name" in self.attributes:
            return self.attributes["name"]
        else:
            return None

class DocElement(object):
    """Represents a docstring enabled element in a code file."""
    
    def __init__(self, XMLElement, parser, decorates = None, parent = None):
        """Initializes the docstring element by parsing out the contents and references.

         - XMLElement: the element from the XML tree for the docstring element.
         - parser: an instance of the DocStringParser() with compiled re objects.
        """
        if XMLElement is not None:
            if XMLElement.text is not None:
                self.contents = re.sub("\s+", " ", XMLElement.text.replace("\n", ""))
            else:
                self.contents = ""
            self.doctype = XMLElement.tag
            self.xml = XMLElement
        else:
            self.contents = ""
            self.doctype = ""
            self.xml = None

        self.references = []        
        #Group is the parent of the docstring (i.e. group), NOT the code element it decorates
        self.group = parent
        
        #Decorates is the code element that the docstring refers to. This is common to all code
        #elements but is only set at a higher level by the code element.
        self.decorates = decorates
        self.attributes = {}
        if XMLElement is not None:
            self.parse(parser, XMLElement)

    def __getstate__(self):
        """Cleans up the object so that it can be pickled without any pointer
        issues."""
        odict = self.__dict__.copy() # copy the dict since we change it
        if self.group is not None and not isinstance(self.group, str):
            odict['parent_name'] = self.group.name
        elif isinstance(self.group, str):
            odict['parent_name'] = self.group
        else:
            odict['parent_name'] = None

        del odict['group']
        return odict

    def __setstate__(self, dict):
        self.group = None
        self.__dict__.update(dict)

    @property
    def pointsto(self):
        """Returns the name of the variable, parameter, etc. that this
        docstring points to by checking whether a "name" attribute is
        present on the docstring."""
        if "name" in self.attributes:
            return self.attributes["name"].lower()
        else:
            return None

    def children(self, doctype):
        """Returns the child XML tag that matches the specified doctype.
        """
        if self.xml is not None:
            return self.xml.findall(doctype)
        else:
            return []

    @staticmethod
    def format_dimension(dimension):
        """Formats the specified <dimension> XML tag for string output."""
        result = ""
        if "type" in dimension.attrib:
            result += "[R" if dimension.attrib["type"] == "row" else "[C"
        else:
            result += "[C"
        if "index" in dimension.attrib:
            result += ":" + dimension.attrib["index"]

        result += "] " + re.sub("\s+", " ", dimension.text.replace("\n", " "))
        return result

    def parse(self, parser, xml):
        """Parses the rawtext to extract contents and references."""
        #We can only process references if the XML tag has inner-XML
        if xml.text is not None:
            matches = parser.RE_REFS.finditer(xml.text)
            if matches:
                for match in matches:
                    #Handle "special" references to this.name and param.name here.
                    self.references.append(match.group("reference"))

        #We also need to get all the XML attributes into the element
        for key in list(xml.keys()):
            self.attributes[key] = xml.get(key)
            
    def __str__(self):
        lines = ["  {}: {}".format(self.doctype.upper(), self.contents)]
        if len(self.attributes) > 0:
            strattr = []
            for key, value in list(self.attributes.items()):
                strattr.append("{} -> {}".format(key, value))
            lines.append("    Attributes:" + "; ".join(strattr))

        if self.decorates is not None:
            lines.append("    Decorates: {}".format(self.decorates))
            
        return "\n".join(lines)
