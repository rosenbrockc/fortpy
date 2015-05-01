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

def copytree(src, dst):
    """Recursively copies the source directory to the destination
    only if the files are newer or modified by using rsync.
    """
    from os import path
    from os import waitpid
    from subprocess import Popen, PIPE

    #Append any trailing / that we need to get rsync to work correctly.
    source = path.join(src, "")
    
    prsync = Popen("rsync -t -u -r {} {}".format(source, dst),
                    shell=True, executable="/bin/bash", stdout=PIPE, stderr=PIPE)
    waitpid(prsync.pid, 0)
    #Redirect the output and errors so that we don't pollute stdout.
    output = prsync.stdout.readlines()
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
    if relpath[0:2] == "./":
        return path.join(xbase, relpath[2::])
    else:
        from os import chdir, getcwd
        cd = getcwd()
        chdir(xbase)
        result = path.abspath(relpath)
        chdir(cd)
        return result
