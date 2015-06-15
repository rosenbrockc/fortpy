#This module handles the printing of XML docstrings into paragraph format
#for output during isense documentation lookups or while generating HTML
#documentation from the code.
from ..elements import ValueElement

def format(element):
    """Formats all of the docstrings in the specified element and its
    children into a user-friendly paragraph format for printing.

    :arg element: an instance of fortpy.element.CodeElement.
    """
    result = []
    if type(element).__name__ in ["Subroutine", "Function"]:
        _format_executable(result, element)
    elif type(element).__name__ == "CustomType":
        _format_type(result, element)
    elif isinstance(element, ValueElement):
        _format_value_element(result, element)

    return '\n'.join(result)

def _format_executable(lines, element, spacer=""):
    """Performs formatting specific to a Subroutine or Function code
    element for relevant docstrings.
    """
    rlines = []
    rlines.append(element.signature)
    _format_summary(rlines, element)

    rlines.append("")
    rlines.append("PARAMETERS")
    for p in element.ordered_parameters:
        _format_value_element(rlines, p)

    rlines.append("")
    _format_generic(rlines, element, ["summary"])

    #Subroutines can have embedded types and functions which need to be handled.
    if len(element.types) > 0:
        rlines.append("\nEMBEDDED TYPES")
        for key, value in list(element.types.items()):
            _format_type(rlines, value, "  ")

    if len(element.executables) > 0:
        rlines.append("\nEMBEDDED EXECUTABLES")
        for key, value in list(element.executables.items()):
            _format_executable(rlines, value, "  ")

    lines.extend([spacer + l for l in rlines])

def _format_type(lines, element, spacer=""):
    """Formats a derived type for full documentation output."""
    rlines = []
    rlines.append(element.signature)
    _format_summary(rlines, element)

    rlines.append("")
    _format_generic(rlines, element, ["summary"])

    if len(element.executables) > 0:
        rlines.append("\nEMBEDDED PROCEDURES")
        for key, value in list(element.executables.items()):
            rlines.append("  {}".format(value.__str__()))
            target = value.target
            if target is not None:
                _format_executable(rlines, target, "  ")

    if len(element.members) > 0:
        rlines.append("\nEMBEDDED MEMBERS")
        for key, value in list(element.members.items()):
            _format_value_element(rlines, value, "  ")

    lines.extend([spacer + l for l in rlines])    

def _format_value_element(lines, element, spacer=""):
    """Formats a member or parameter for full documentation output."""
    lines.append(spacer + element.definition())
    _format_summary(lines, element)
    _format_generic(lines, element, ["summary"])

def _format_summary(lines, element, spacer="  "):
    """Adds the element's summary tag to the lines if it exists."""
    summary = spacer + element.summary
    if element.summary == "":
        summary = spacer + "No description given in XML documentation for this element."
    lines.append(summary)

def _format_generic(lines, element, printed, spacer=""):
    """Generically formats all remaining docstrings and custom XML
    tags that don't appear in the list of already printed documentation.

    :arg printed: a list of XML tags for the element that have already
      been handled by a higher method.
    """
    for doc in element.docstring:
        if doc.doctype.lower() not in printed:
            lines.append(spacer + doc.__str__())
