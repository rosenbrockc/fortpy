"""Utility functions used by multiple modules in the package."""
def get_attrib(xml, name, tag=None, cast=str, default=None):
    """Returns the specified attribute from the XML element, raising a ValueError
    if it is not availaible.
    
    :arg xml: the XMLElement instance to get the attribute from.
    :arg name: the name of the attribute in xml.attrib dictionary.
    :arg tag: the name of the tag to display with the ValueError if the attribute is missing.
    """
    if name in xml.attrib:
        return cast(xml.attrib[name])
    elif default is not None:
        return default
    elif tag is not None:
        raise ValueError("'{}' is a required attribute of <{}> tag.".format(name, tag))

def copyfile(src, dst, verbose=False):
    """Copies the specified source file to destination *file* if it is newer
    or does not yet exist.
    """
    from os import waitpid
    from subprocess import Popen, PIPE
    prsync = Popen("rsync {}-t -u {} {}".format("-v " if verbose else "", src, dst),
                   shell=True, executable="/bin/bash", stdout=PIPE, stderr=PIPE)
    waitpid(prsync.pid, 0)
    
    #Redirect the output and errors so that we don't pollute stdout.
    error = prsync.stderr.readlines()

    if len(error) > 0:
        from fortpy.msg import warn
        warn("Error while copying {} using rsync.\n\n{}".format(src, '\n'.join(error)))
        
def copy(src, dst, verbose=False):
    """Copies the specified source file to destination directory if it is newer
    or does not yet exist.
    """
    from os import waitpid, path
    from subprocess import Popen, PIPE
    desti = path.join(dst, "")
    if not path.isdir(desti):
        from os import mkdir
        mkdir(desti)

    prsync = Popen("rsync -t -u {} {}".format(src, desti),
                   shell=True, executable="/bin/bash", stdout=PIPE, stderr=PIPE)
    waitpid(prsync.pid, 0)
    
    #Redirect the output and errors so that we don't pollute stdout.
    #output = prsync.stdout.readlines()
    error = prsync.stderr.readlines()

    if len(error) > 0:
        from fortpy.msg import warn
        warn("Error while copying {} using rsync.\n\n{}".format(src, '\n'.join(error)))
    
def copytree(src, dst):
    """Recursively copies the source directory to the destination
    only if the files are newer or modified by using rsync.
    """
    from os import path, waitpid
    from subprocess import Popen, PIPE

    #Append any trailing / that we need to get rsync to work correctly.
    source = path.join(src, "")
    desti = path.join(dst, "")
    if not path.isdir(desti):
        from os import mkdir
        mkdir(desti)
    
    prsync = Popen("rsync -t -u -r {} {}".format(source, desti),
                    shell=True, executable="/bin/bash", stdout=PIPE, stderr=PIPE)
    waitpid(prsync.pid, 0)
    #Redirect the output and errors so that we don't pollute stdout.
    #output = prsync.stdout.readlines()
    error = prsync.stderr.readlines()

    if len(error) > 0:
        from fortpy.msg import warn
        warn("Error while copying {} using rsync.\n\n{}".format(source, '\n'.join(error)))

def get_fortpy_templates_dir():
    """Gets the templates directory from the fortpy package."""
    import fortpy
    from os import path
    fortdir = path.dirname(fortpy.__file__)
    return path.join(fortdir, "templates")

def set_fortpy_templates(obj, fortpy_templates=None):
    """Sets the directory path for the fortpy templates. If no directory
    is specified, use the default one that shipped with the package.
    """
    #If they didn't specify a custom templates directory, use the default
    #one that shipped with the package.
    from os import path
    if fortpy_templates is not None:
        obj.fortpy_templates = path.abspath(path.expanduser(fortpy_templates))
    else:
        from fortpy.utility import get_fortpy_templates_dir
        obj.fortpy_templates = get_fortpy_templates_dir()

def get_dir_relpath(base, relpath):
    """Returns the absolute path to the 'relpath' taken relative to the base
    directory.

    :arg base: the base directory to take the path relative to.
    :arg relpath: the path relative to 'base' in terms of '.' and '..'.
    """
    from os import path
    xbase = path.abspath(path.expanduser(base))
    if not path.isdir(xbase):
        return path.join(xbase, relpath)
    
    if relpath[0:2] == "./":
        return path.join(xbase, relpath[2::])
    else:
        from os import chdir, getcwd
        cd = getcwd()
        chdir(xbase)
        result = path.abspath(relpath)
        chdir(cd)
        return result

import io
import itertools as IT
import xml.etree.ElementTree as ET

def wrap_line(line, limit=None, chars=80):
    """Wraps the specified line of text on whitespace to make sure that none of the
    lines' lengths exceeds 'chars' characters.
    """
    result = []
    builder = []
    length = 0
    if limit is not None:
        sline = line[0:limit]
    else:
        sline = line
        
    for word in sline.split():
        if length <= chars:
            builder.append(word)
            length += len(word) + 1
        else:
            result.append(' '.join(builder))
            builder = [word]
            length = 0

    result.append(' '.join(builder))
    return result

def x_parse_error(err, content, source):
    """Explains the specified ParseError instance to show the user where the error
    happened in their XML.
    """
    lineno, column = err.position
    if "<doc>" in content:
        #Adjust the position since we are taking the <doc> part out of the tags
        #since we may have put that in ourselves.
        column -= 5
    start = lineno - (1 if lineno == 1 else 2)
    lines = []
    tcontent = content.replace("<doc>", "").replace("</doc>", "").split('\n')
    for context in IT.islice(tcontent, start, lineno):
        lines.append(context.strip())
    last = wrap_line(lines.pop(), column)
    lines.extend(last)
    
    caret = '{:=>{}}'.format('^', len(last[-1]))
    if source is not None:
        err.msg = '\nIn: {}\n{}\n{}\n{}'.format(source, err, '\n'.join(lines), caret)
    else:
        err.msg = '{}\n{}\n{}'.format(err, '\n'.join(lines), caret)

    raise err

def XML(content, source=None):
    """Parses the XML text using the ET.XML function, but handling the ParseError in
    a user-friendly way.
    """
    try:
        tree = ET.XML(content)
    except ET.ParseError as err:
        x_parse_error(err, content, source)
    return tree

def XML_fromstring(content, source=None):
    """Parses the XML string into a node tree. If an ParseError exception is raised,
    the error message is formatted nicely to show the badly formed XML to the user.
    """
    try:
        tree = ET.fromstring(content)
    except ET.ParseError as err:
        x_parse_error(err, content, source)
    return tree

def which(program):
    """Tests whether the specified program is anywhere in the environment
    PATH so that it probably exists."""
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

