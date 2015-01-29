# Using the magic encoding
# -*- coding: utf-8 -*-
from .. import msg
import xml.etree.ElementTree as ET
import re
from ..docelements import DocElement, DocGroup
import os

class DocStringParser(object):
    """Parses the XML tags from the custom docstrings in our fortran code."""

    def __init__(self):
        self.setup_regex()

    def setup_regex(self):
        """Sets up the patterns and regex objects for parsing the docstrings."""
        #Regex for grabbing out valid XML tags that represent known docstrings that we can work with.
        self.keywords = [ "summary", "usage", "errors", "member", "group", "local", 
                          "comments", "parameter" ]
        #Regex for extracting the contents of docstrings minus the !! and any leading spaces.
        self._RX_DOCS = "^\s*!!(?P<docstring>.+?)$"
        self.RE_DOCS = re.compile(self._RX_DOCS, re.M)
        #Regex for handling cross references in the documentation
        self._RX_REFS = r"@CREF\[(?P<reference>[^\]]+)\]"
        self.RE_REFS = re.compile(self._RX_REFS)
        #Regex to match first lines of declarations for code elements that can be
        #decorated by docstrings.
        self._RX_DECOR = (r"((?P<type>character|real|type|logical|integer)?"
                          r"(?P<kind>\([a-z0-9_]+\))?)?(,?(?P<modifiers>[^\n]+?))?"
                          r"\s*(?P<codetype>subroutine|function|type|module|interface)\s+(?P<name>[^(]+)")
        self.RE_DECOR = re.compile(self._RX_DECOR, re.I)
        #Regex for getting the docstrings decorating one or more modules in a code file,
        #Since they aren't contained inside any other code element, we can't just use
        #the normal docblocks routines.
        self._RX_MODDOCS = (r"^(?P<docstring>\s*!!.+?)\n\s*module\s+(?P<name>[A-Za-z0-9_]+)"
                            ".+?end\s+module(\s+(?P=name))?")
        self.RE_MODDOCS = re.compile(self._RX_MODDOCS, re.DOTALL | re.I)

    def parse_docs(self, string, container = None):
        """Parses the docstrings from the specified string that is the contents of container.

        Returns a dictionary with keys as parent.code_element_name and the values
        a list of XML elements for corresponding docstrings.

        :arg container: the instance of the element who owns the string.
        """
        result = {}
        if container is None:
            #We are working with the code file at the module level. Extract the module
            #docstrings and XML and return the dictionary with module names as keys.
            for module in self.RE_MODDOCS.finditer(string):
                docstring = re.sub("\s*!!", "", module.group("docstring"))
                doctext = "<doc>{}</doc>".format(re.sub("\n", "\s", docstring))
                try:
                    docs = ET.XML(doctext)
                    #Get the name of the module to use as the key and then add the list
                    #of XML docstrings to the result.
                    key = module.group("name")
                    if not key in result:
                        result[key] = [list(docs), module.start(), module.end()]
                    else:
                        result[key][0].extend(list(docs))
                except ET.ParseError:
                    msg.err(doctext)
        else:
            #This is the text content of a code element that was buried inside of the module.
            #Get all the docblocks and the items they decorate from this parent.
            result = self._parse_docblocks(string, container)
            
        return result
        
    def _process_docgroup(self, group, code_el, add=True):
        """Explodes the group members into a list; adds the group to the
        specified code element and updates the group value for each
        of the docstring elements in the group.

        :arg add: when true, docgroups must be unique in the code element;
          otherwise, existing groups are overwritten."""
        if group.name in code_el.groups and add:
            msg.warn("duplicate group names in code element {}".format(code_el.name))
        else:
            code_el.groups[group.name] = group

        kids = self.to_doc(list(group.xml), group.decorates)
        for child in kids:
            child.group = group.name

        return kids

    def process_execdocs(self, docs, anexec, key, add=True):
        """Associates parameter documentation with parameters for the executable
        and any remaining docs with the executable itself.

         - key: the module.executable identifier for the function or subroutine.
        """
        #Paramdocs has a list of docstrings for summary, usage, parameters, etc.
        #check which belong to parameters and associate them, otherwise append
        #them to the executable.
        for doc in docs:
            if doc.doctype == "parameter":
                if doc.pointsto is not None and doc.pointsto in anexec.parameters:
                    if add:
                        anexec.parameters[doc.pointsto].docstring.append(doc)
                    else:
                        anexec.parameters[doc.pointsto].overwrite_docs(doc)
                else: 
                    #the parameter docstring is orphaned, give a warning.
                    wmsg = "the docstring for parameter '{}' had no corresponding " + \
                          "parameter in the executable definition for '{}'."
                    msg.warn(wmsg.format(doc.pointsto, anexec))
            elif doc.doctype == "group":
                kids = self._process_docgroup(doc, anexec)
                if add:
                    anexec.docstring.extend(kids)
                else:
                    for kid in kids:
                        anexec.overwrite_docs(kid)
            else:
                #The docstring must be for the executable
                if add:
                    anexec.docstring.append(doc)
                else:
                    anexec.overwrite_docs(doc)
                    
    def process_embedded(self, xlist, anexec, add=True):
        """Processes the specified xml list and executable to link *embedded*
        types and executables to their docstrings.

        :arg xlist: a list of XML elements returned by parse_docs().
        :arg add: when true, docstrings are only appended, never overwritten.
        """
        #Keep track of the changes that took place in the lengths of the 
        #docstrings that got added/updated on the elements children.
        delta = 0
        for t in anexec.types:
            key = "{}.{}".format(anexec.name, t)
            if key in xlist:
                docs = self.to_doc(xlist[key][0], t)
                self.process_memberdocs(docs, anexec.types[t], add)
                anexec.types[t].docstart = xlist[key][1]
                delta += xlist[key][2] - anexec.types[t].docend
                anexec.types[t].docend = xlist[key][2]

        for iexec in anexec.executables:
            key = "{}.{}".format(anexec.name, iexec)
            if key in xlist:
                docs = self.to_doc(xlist[key][0], t)
                self.process_memberdocs(docs, anexec.executables[iexec], add)
                anexec.executables[iexec].docstart = xlist[key][1]
                delta += xlist[key][2] - anexec.executables[iexec].docend
                anexec.executables[iexec].docend = xlist[key][2]

        if not add:
            return delta

    def process_memberdocs(self, docs, codeEl, add=True):
        """Associates member type DocElements with their corresponding members
        in the specified code element. The element must have a dictionary of
        members already."""        
        #Now we need to associate the members with their docstrings
        #Some of the members may be buried inside a group tag and
        #need to be handled separately.
        remainingdocs = []
        expandeddocs = []

        #Process any groups that are in the doc list.
        for doc in docs:
            if isinstance(doc, DocGroup):
                kids = self._process_docgroup(doc, codeEl, add)
                expandeddocs.extend(kids)
            else:
                expandeddocs.append(doc)

        for doc in expandeddocs:
            #Process the docstring, if it doesn't belong to a member
            #we will add it to the list of unassigned docstrings,
            #these most likely point to type declarations.           
            if not self._process_docstrings(doc, codeEl.members, add):
                remainingdocs.append(doc)
        return remainingdocs

    def _process_docstrings(self, doc, members, add=True):
        """Adds the docstrings from the list of DocElements to their
        respective members.

        Returns true if the doc element belonged to a member."""
        if ((doc.doctype == "member" or doc.doctype == "local") and 
            doc.pointsto is not None and 
            doc.pointsto in members):
            if add:
                members[doc.pointsto].docstring.append(doc)
            else:
                members[doc.pointsto].overwrite_docs(doc)
            return True
        else:
            return False

    def to_doc(self, xmllist, decorates):
        """Converts the specified xml list to a list of docstring elements."""
        result = []
        for xitem in xmllist:
            if xitem.tag != "group":
                #The docstring allows a single string to point to multiple
                #names in a comma-separated list in the names attribute.
                if "name" in list(xitem.keys()):
                    names = re.split("[\s,]+", xitem.get("name"))
                    for name in names:                        
                        #Once we have created the DocElement, we need to override
                        #its name attribute (which will have the comma-separated
                        #list) with the single name
                        docel = DocElement(xitem, self, decorates)
                        docel.attributes["name"] = name
                        result.append(docel)
                else:
                    #This docstring doesn't have a name attribute, just add it
                    result.append(DocElement(xitem, self, decorates))
            else:
                docel = DocGroup(xitem, decorates)
                result.append(docel)
        return result
        
    def _parse_docblocks(self, string, container):
        """Parses all the docstrings out of the specified string.

        Returns a dictionary of docstrings with the key as parent.code_element_name
        and the value a list of XML elements that contain docstrings.
        """
        #The easiest way to do this is to look at one line at a time and see if it is a docstring
        #When we find a group of docstrings that suddenly ends, the next item is the code element
        #that they were decorating (which may or may not be pertinent).
        current = []
        docblocks = {}
        docstart = 0

        for line in string.split("\n"):
            match = self.RE_DOCS.match(line)
            if match is not None:
                current.append(match.group("docstring"))
                if len(current) == 1:
                    #This was the first docstring of a new documentation block.
                    docend = docstart + len(line) + 1  # +1 for \n removed by split()
                else:
                    #We already have some docstrings in the block, update start/end
                    docend += len(line) + 1
            else:
                #See if we were previously working on a docstring block or not.
                if len(current) > 0:
                    #Save all the current docstrings in the blocks dictionary
                    #under the name of the code element in this line.
                    key = self._parse_docline(line, container)
                    #If the docblock has multiple XML tags at the same depth, the XML
                    #parser will scream. Wrap everything in a doc tag.
                    doctext = "<doc>{}</doc>".format(" ".join(current))
                    try:
                        docs = ET.XML(doctext)
                        if not key in docblocks:
                            #Let the docstart and docend *always* be absolute 
                            #character references.
                            absstart, absend = container.module.absolute_charindex(string, 
                                                                                   docstart,
                                                                                   docend-len(line))
                            docblocks[key] = [list(docs), absstart, absend]
                        else:
                            docblocks[key][0].extend(list(docs))
                    except ET.ParseError:
                        msg.err(doctext)

                    #Reset the list of current docstrings
                    current = []
                    docstart = docend + len(line) + 1
                else:
                    #We need to keep track of the line lengths for docstart/end.
                    docstart += len(line) + 1

        return docblocks

    def _parse_docline(self, line, container):
        """Parses a single line of code following a docblock to see if
        it as a valid code element that can be decorated. If so, return
        the name of the code element."""
        match = self.RE_DECOR.match(line)
        if match is not None:
            return "{}.{}".format(container.name, match.group("name"))
        else:
            return container.name

    def parsexml(self, xmlstring, modules):
        """Parses the docstrings out of the specified xml file."""
        result = {}

        xmlroot = ET.fromstring(xmlstring)
        if xmlroot.tag == "fortpy" and "mode" in xmlroot.attrib and \
           xmlroot.attrib["mode"] == "docstring":
            #We fill the dictionary with decorates names as keys and lists
            #of the xml docstring elements as values.
            for child in xmlroot:
                xmltags = []
                if child.tag == "decorates" and "name" in child.attrib:
                    decorates = child.attrib["name"]
                    xmltags.extend(list(child))
                elif "decorates" in child.attrib:
                    decorates = child.attrib["decorates"]
                    xmltags.append(child)

                if decorates in result:
                    result[decorates].extend(xmltags)
                else:
                    result[decorates] = xmltags

            #Loop through all the docstrings we found and team them up with
            #their respective module members.
            self._xml_update_modules(result, modules)

    def _xml_update_modules(self, xmldict, modules):
        """Updates the docstrings in the specified modules by looking for
        docstrings in the xmldict."""
        for kdecor in xmldict:
            modname, memname = kdecor.split(".")
            if modname in modules:
                module = modules[modname]
                #We only need to check the members, types and executables
                memname = memname.lower()
                if memname in module.members:
                    docs = self.to_doc(xmldict[kdecor], modname)
                    self.process_memberdocs(docs, module)
                elif memname in module.types:
                    member = module.types[memname]
                    docs = self.to_doc(xmldict[kdecor], memname)
                    member.docstring.extend(docs)
                elif memname in module.executables:
                    member = module.executables[memname]
                    docs = self.to_doc(xmldict[kdecor], memname)
                    self.process_execdocs(docs, member, kdecor)
                else:
                    msg.warn("orphaned docstring. No member {} in module {}.".format(
                        memname, modname))
            else:
                msg.warn("orphaned docstring from XML docfile for {}".format(kdecor))

    def rt_update_module(self, xmldict, module):
        """Updates the members, executables and types in the specified module
        to have the latest docstring information from the xmldict.
        """
        #This keeps track of how many character were added/removed by
        #updating the docstrings in xmldict.
        delta = 0
        for kdecor in xmldict:
            if "." in kdecor:
                modname, memname = kdecor.split(".")
            else:
                modname, memname = module.name, None

            if module.name == modname:
                #This tag is relevant to the specified module. Continue
                xlist, docstart, docend = xmldict[kdecor]

                #We only need to check the members, types and executables
                #For executables and types, we need to update the docstart and
                #docend attributes since their docstrings must come as a single
                #block immediately preceding the signature, so that our values
                #from the updater will be correct.
                if memname in module.types:
                    member = module.types[memname]
                    docs = self.to_doc(xlist, memname)
                    member.docstring = docs
                    delta += self._rt_update_docindices(member, docstart, docend)
                elif memname in module.executables:
                    member = module.executables[memname]
                    docs = self.to_doc(xlist, memname)
                    self.process_execdocs(docs, member, kdecor, False)
                    delta += self._rt_update_docindices(member, docstart, docend)
                else:
                    #Since it didn't point to anything else, it must be for the
                    #members of the module.
                    docs = self.to_doc(xlist, modname)
                    self.process_memberdocs(docs, module, False)               
                
        return delta

    def _rt_update_docindices(self, element, docstart, docend):
        """Updates the docstart, docend, start and end attributes for the 
        specified element using the new limits for the docstring."""
        #see how many characters have to be added/removed from the end
        #of the current doc limits.
        delta = element.docend - docend
        element.docstart = docstart
        element.docend = docend
        element.start += delta
        element.end += delta

        return delta
