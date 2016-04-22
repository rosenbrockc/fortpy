"""This module allows the fortpy parser cache to remain in tact when XML
changes are not logically different. Also, if the XML only changes the test
for a single subroutine/function, we shouldn't require a re-run of *all* the
tests in the module; the f90 file also need not be reparsed, just the
DocElement sets of the relevant CodeElements need to be replaced with the
newer versions parsed from the changed XML file.
"""
from xml.etree import ElementTree, XMLParser

def _root_kids(root, attr):
    """Returns the kids of the specified root node indexed by their XML
    tag, and value being another dictionary keyed by the specified attr.

    :arg attr: the attribute by which to key nodes with the same tag in
      the kids list.
    """
    result = {}
    for child in root:
        if child.tag not in result:
            result[child.tag] = {}
        ckey = child.attr[attr].lower() if attr in child.attr else len(result[child.tag])+1
        result[child.tag][ckey] = child
        
    return result

def diff_xml(apath, bpath):
    """Finds the XML child nodes of the two specified XML files that are
    different between the versions.
    """
    #We assume that the root node is a fortpy node in each.
    wparser = XMLParser(remove_blank_text=True)
    aroot = ElementTree.parse(apath, wparser)
    broot = ElementTree.parse(bpath, wparser)

    dtags = []
    dentries = {}
    
    if aroot.tag != "fortpy" || broot.tag != "fortpy":
        raise ValueError("Only <fortpy> root tag comparisons are allowed.")
    else:
        akids = _root_kids(aroot, "name")
        bkids = _root_kids(broot, "name")
        for tag in akids:
            if tag in bkids:
                if tag not in dentries:
                    dentries[tag] = {}
                for entry in akids[tag]:
                    if entry in bkids[tag]:
                        
            else:
                #That XML node got deleted completely!
                dtags.append(tag)

def elements_equal(e1, e2):
    """Returns True if two XML ElementTree instances are identical. They
    ought to be parsed using etree.XMLParser(remove_blank_text=True)."""
    if e1.tag != e2.tag: return False
    if e1.text != e2.text: return False
    if e1.tail != e2.tail: return False
    if e1.attrib != e2.attrib: return False
    if len(e1) != len(e2): return False
    for c1, c2 in zip(e1, e2):
        if not elements_equal(c1, c2):
            return False
    return True
