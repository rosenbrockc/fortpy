import fortpy
import os
import re
import xml.etree.ElementTree as ET
from fortpy.elements import ValueElement, Subroutine, Function
from fortpy.docelements import DocElement
import fortpy.config
import sys
config = sys.modules["config"]

def load(parser, serializer):
    """Returns a dictionary of builtin functions for Fortran. Checks the
    cache first to see if we have a serialized version. If we don't, it
    loads it from the XML file.

    :arg parser: the DocParser instance for parsing the XML tags.
    :arg serializer: a Serializer instance from the CodeParser to cache
      the loaded XML file.
    """
    fortdir = os.path.dirname(fortpy.__file__)
    xmlpath = os.path.join(fortdir, "isense", "builtin.xml")
    if not os.path.isfile(xmlpath):
        return {}

    changed_time = os.path.getmtime(xmlpath)
    cached = serializer.load_module("builtin.xml", changed_time)
    if cached is None:
        result = _load_builtin_xml(xmlpath, parser)
        serializer.save_module("builtin.xml", result, changed_time)
    else:
        result = cached

    return result

def _isense_builtin_symbol(symbol):
    """Returns the symbol with modifications according to settings in the
    <isense> configuration section.
    """
    if "builtin" in config.isense and "uppercase" in config.isense["builtin"]:
        if config.isense["builtin"]["uppercase"] == "true":
            return symbol.upper()
        else:
            return symbol.lower()
    else:
        #Make it lowercase by default
        return symbol.lower()

def _load_builtin_xml(xmlpath, parser):
    """Loads the builtin function specifications from the builtin.xml file.

    :arg parser: the DocParser instance for parsing the XML tags.
    """
    #First we need to get hold of the fortpy directory so we can locate
    #the isense/builtin.xml file.
    result = {}

    el = ET.parse(xmlpath).getroot()
    if el.tag == "builtin":
        for child in el:
            anexec = _parse_xml(child, parser)
            result[anexec.name.lower()] = anexec

    return result
                
def _parse_xml(child, parser):
    """Parses the specified child XML tag and creates a Subroutine or
    Function object out of it."""
    name, modifiers, dtype, kind = _parse_common(child)

    #Handle the symbol modification according to the isense settings.
    name = _isense_builtin_symbol(name)

    if child.tag == "subroutine":
        parent = Subroutine(name, modifiers, None)
    elif child.tag == "function":
        parent = Function(name, modifiers, dtype, kind, None)

    if parent is not None:
        for kid in child:
            if kid.tag == "parameter":
                _parse_parameter(kid, parser, parent)
            elif kid.tag == "summary":
                _parse_summary(kid, parser, parent)
            elif kid.tag == "usage":
                _parse_usage(kid, parser, parent)

    return parent

def _parse_summary(tag, parser, parent):
    """Parses a <summary> tag and adds it the Executable parent instance.
    
    :arg parser: an instance of DocParser to create the DocElement with.
    """
    summary = DocElement(tag, parser, parent)
    parent.docstring.append(summary)

def _parse_usage(tag, parser, parent):
    """Parses a <usage> tag and adds it the Executable parent instance.
    
    :arg parser: an instance of DocParser to create the DocElement with.
    """
    usage = DocElement(tag, parser, parent)
    parent.docstring.append(usage)

def _parse_parameter(tag, parser, parent):
    """Creates a ValueElement instance from the specified <parameter> XML tag."""
    name, modifiers, dtype, kind = _parse_common(tag)
    if "default" in tag.attrib:
        default = tag.attrib["default"]
    else:
        default = None
    if "dimension" in tag.attrib:
        dimension = tag.attrib["dimension"]
    else:
        dimension = None

    result = ValueElement(name, modifiers, dtype, kind, default, dimension, parent)
    doc = DocElement(tag, parser, result)
    result.docstring.append(doc)
    parent.add_parameter(result)

def _parse_common(tag):
    """Returns a tuple of (name, modifiers, dtype, kind)
    for the specified tag. Any missing attributes will have values of None.
    """
    if "modifiers" in tag.attrib:
        modifiers = re.split(",\s*", tag.attrib["modifiers"].strip())
        if "" in modifiers:
            modifiers.remove("")
    else:
        modifiers = None

    if "name" in tag.attrib:
        name = tag.attrib["name"]
    if "type" in tag.attrib:
        dtype = tag.attrib["type"]
    else:
        dtype = None
    if "kind" in tag.attrib:
        kind = tag.attrib["kind"]
    else:
        kind = None

    return (name, modifiers, dtype, kind)
