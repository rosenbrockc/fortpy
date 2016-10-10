from .. import msg
from fortpy.testing.templates import FileLine
from fortpy.docelements import DocElement, DocGroup
from fortpy.printing.formatting import present_params
from os import path, remove
from fortpy.utility import copyfile
import re
def _expand_cases(casestr):
    """Returns a list of case identifiers from the shorthand string.
    """
    # For specifying the cases in a unit test, we should allow ranges like
    # "standard.cr[1-12]" so that the developer doesn't need to enter each
    # of the cases separately. We should still allow a comma-separated list
    # of cases, but each must allow the shorthand notation.
    rawcases = re.split(",\s*", casestr)
    if "[" in casestr:
        cases = []
        rxcase = r"(?P<prefix>[^[]*)\[(?P<range>\d+-\d+)](?P<suffix>.*)"
        recase = re.compile(rxcase)
        for craw in rawcases:
            m = recase.match(craw)
            if m:
                prefix = m.group("prefix")
                vals = list(map(int, m.group("range").split("-")))
                suffix = m.group("suffix")
                if prefix is None:
                    prefix = ""
                if suffix is None:
                    suffix = ""
                for v in range(vals[0], vals[1]+1):
                    cases.append("{}{}{}".format(prefix, v, suffix))
            else:
                cases.append(craw)
        return cases
    else:
        return rawcases

class TestSource(object):
    """Represents a link between the output of one unit test that is being used
    as input for another unit test.
    """
    def __init__(self, compkey, testspec):
        """Construct the test source.
        
        :arg compkey: module.executable.parameter:testid:cases.
        :arg testspec: the instance of TestSpecification for the <value> or <input>
          tag using the 'testsource" attribute.
        """
        subkeys = compkey.split(":")
        self.key = subkeys[0].strip().split('.')
        """The 'module.executable.parameter' directive."""
        self.testid = None if len(subkeys) < 2 else subkeys[1]
        """The identifier of the test whose output will be used."""
        self.cases = None
        """The list of test cases to use inputs for. One test case will be created
        for each distinct test source case using either 'one-to-one' or 'cartesian'.
        """
        if len(subkeys) == 3:
            if subkeys[2] in ["match", "product"]:
                self.cases = subkeys[2]
            else:
                self.cases = _expand_cases(subkeys[2])

        self.testspec = testspec
        """The TestSpecification instance of the target whose input is being defined.
        """        
        self.source = self._source_parameter()
        """The CodeElement instance of the parameter whose output is being used.
        """
        self.target = self._source_testoutput(*self.source)
        """The TestTarget instance for the variable whose output is being used as
        the test source.
        """

    @property
    def varname(self):
        """Returns the name of the local file *variable* that will contain
        the full path to the data to use.
        """
        return "fpy_{}_ts".format(self.source[1].name)
    
    @property
    def filevar(self):
        """Returns the code to reference a trimmed, clean version of the contents
        of the variable containing the full path to the data.
        """
        return "trim(adjustl({}))".format(self.varname)

    @property
    def filename(self):
        """Returns the name of the local file that will contain
        the full path to the data to use.
        """
        return "fpy_{}_ts".format(self.source[1].name)

    def _is_relevant(self, testcase, sourcespec):
        """Returns True if the specified target case identifier matches
        the specification in self.cases.
        """
        #First check that the test names match.
        if self.testid is not None:
            #We are filtering on a specific test id *and* the existence of that
            #test in the source specification.
            testok = self.testid == sourcespec.identifier
        else:
            #We just need to make sure that the test exists
            testok = True

        if not testok:
            return False
        
        if self.cases == "every":
            #Check whether the specified test case matches any of those for
            #which we have sample output.
            if (testcase in sourcespec.cases or (testcase == "" and sourcespec.cases is None)):
                return True
        if isinstance(self.cases, list) and len(self.cases) > 0:
            #Check that the testcase is in our filtering list *and* that it exists
            #in the set of output (source) test cases.
            return testcase in self.cases and testcase in sourcespec.cases
        return False

    def path(self, coderoot, case, compiler=None, stagedir=None):
        """Returns the full path to the *source* file that this test source
        is going to use (either directly or with copy).

        :arg stagedir: the path to the staging directory to use (override).
        """
        from os import path
        from fortpy.testing.compilers import replace
        identifier = replace("{}.{}.[c]".format(xinst.module.name, xinst.name), compiler)

        if not self._is_relevant(case, self.target.testspec):
            return None
        
        if stagedir is not None:
            folder = path.abspath(stagedir)
        elif xinst.group is not None and xinst.group.staging is not None:
            #Determine the absolute path using the relative path specified in the group
            #staging directory attribute.
            from fortpy.tramp import coderelpath
            folder = coderelpath(coderoot, xinst.group.staging)

        #This is the path to the source output folder that we need to copy from
        if case != "":
            folder = path.join(folder, identifier, "tests", self.target.testspec.identifier)
        else:
            caseid = "{}.{}".format(self.target.testspec.identifier, case)
            folder = path.join(folder, identifier, "tests", caseid)

        if path.isdir(folder):
            source = path.join(folder, self.target.varfile)
            if path.isfile(source):
                return source
            else:
                raise ValueError("File '{}' for 'testsource' does not exist.".format(source))
        else:
            raise ValueError("Can't find source folder for 'testsource' at '{}'".format(folder))
    
    def copy(self, coderoot, testroot, case, compiler=None, stagedir=None):
        """Copies the test-source file to the specified testing directory.
        
        :arg stagedir: the path to the staging directory to use (override).
        """
        source = self.path(coderoot, case, compiler, stagedir)
        if source is not None:
            copyfile(source, testroot)
        
    def _source_testoutput(self, executable, parameter):
        """Returns the test output source file information for the specified
        parameter from a unit test.
        """
        if executable.test_group is None:
            return None
        result = None
        for tname, tinst in executable.test_group.tests:
            for target in tinst.targets:
                if (parameter == target.name and
                    (target.testspec.identifier == self.testid or self.testid is None)):
                    result = target
                    break

        return result
        
    def _source_parameter(self):
        """Returns the parameter CodeElement instance defined by the testsource
        attribute if it was specified. Result is a tuple of (executable, parameter)
        instances.
        """
        if len(self.key) == 3:
            emsg = "Couldn't find the test-source {} '{}' for test '{}' in '{}'."
            tid, xid = self.testspec.identifier, self.testspec.testgroup.executable.name
            mname, xname, pname = self.key
            parser = self.testspec.testgroup.executable.module.parent
            
            if mname not in parser.modules:
                parser.load_dependency(mname, True, True)
            if mname in parser.modules:
                if xname in parser.modules[mnam].executables:
                    xinst = parser.modules.executables[xname]
                    if pname in xinst.parameters:
                        return (xinst, xinst.parameters[pname])
                    else:
                        raise ValueError(emsg.format("parameter", xname, tid, xid))
                else:
                    raise ValueError(emsg.format("executable", xname, tid, xid))
            else:
                raise ValueError(emsg.format("module", mname, tid, xid))
        else:
            raise ValueError("The 'testsource' attribute requires a list of 3 "
                             "'.'-separated identifiers 'module.executable.parameter'. "
                             "Specified value was '{}'.".format(testsource))               
    
class GlobalDeclaration(object):
    """Represents a declaration to have a global variable instance available
    for execution of a unit test.

    :attr element: the docstring element representing the <global> tag."""
    def __init__(self, element):
        self.element = element

    @property
    def dependency(self):
        """Returns the name of a variable that this variable needs in order
        to be initialized or assigned a value or "".
        """
        if "default" in self.attributes:
            return self.attributes["default"].lower()
        else:
            return ""

    @property
    def position(self):
        """Returns the position of a <global> declaration in a testing group relative
        to the test-local <global> definitions.
        """
        if "position" in self.attributes:
            return self.attributes["position"]
        else:
            return "before"
        
    @property
    def dimension(self):
        """Returns the dimension specification from the <global> tag."""
        if "dimensions" not in self.attributes:
            return None
        else:
            return self.attributes["dimensions"]

    @property
    def kind(self):
        """Returns the kind of the variable from the <global> tag declaration, if 
        it exists.
        """
        if "kind" in self.attributes:
            return self.attributes["kind"]
        else:
            return None

    @property
    def D(self):
        """Returns the integer number of dimensions that this variable
        is declared as having."""
        if self.dimension is None:
            return 0
        else:
            return self.dimension.count(",") + 1

    @property
    def ignore(self):
        """Specifies whether this variable should be ignored when passed in to argument lists
        and for definitions etc.
        """
        return "ignore" in self.attributes and self.attributes["ignore"] == "true"

    def value_elem(self, module):
        """Returns a fortpy.element.ValueElement to represent this GlobalDeclaration.
        """
        from fortpy.elements import ValueElement
        a = self.attributes
        default = None if "default" not in a else a["default"]
        modifiers = [] if "modifiers" not in a else a["modifiers"]
        return ValueElement(a["name"], modifiers, a["type"], self.kind,
                            default, self.dimension, module, self.D)
    
    def compare(self, element):
        """Determines whether the specified element defines the same variable
        type and name as this one.

        The two are considered equal if they have the same:
         - name
         - type
         - dimensionality
         - allocatable/pointer modifiers
        ."""
        if self.ignore:
            return False
        elif self.attributes["name"].lower() != element.attributes["name"].lower():
            return False
        elif not self._kind_check(element):
            return False
        elif not self._modifier_check(element):
            return False
        else:
            return self._dimensions_check(element)
            
    def _kind_check(self, element):
        """Checks whether the kind declaration of the element matches this global."""
        if ("kind" not in self.attributes and "kind" not in element.attributes) or \
           ("kind" in self.attributes and "kind" in element.attributes
            and self.attributes["kind"].lower() == element.attributes["kind"].lower()):
            return True
        else:
            return False

    def _dimensions_check(self, element):
        """Checks whether the dimensions of this declaration match the specified
        element."""
        if not ("dimensions" in self.attributes and "dimensions" in element.attributes):
            return True
        elif "dimensions" in self.attributes and "dimensions" in element.attributes:
            #The dimension text has to match perfectly. If variables names are specified
            #for bounds, we have no way of knowing whether the sizes are the same before
            #runtime. However, we can do some cleanup befor comparing.
            match = True
            selfdim = self.attributes["dimensions"].lower().split(",")
            eldim = element.attributes["dimensions"].lower().split(",")

            i = 0
            #We only need to compare dimensions until one fails
            while match and i < len(selfdim):
                if selfdim[i].strip() != eldim[i].strip():
                    match = False
                i += 1
                
            return match
        else:
            return False

    def _modifier_check(self, element):
        """Checks whether this global var and the specified element match in crucial
        modifier definitions (i.e. pointer, allocatable)."""
        #We need both to have modifiers or not have modifiers
        if not ("modifiers" in self.attributes and "modifiers" in element.attributes):
            return True
        elif "modifiers" in self.attributes and "modifiers" in element.attributes:
            selfmods = self.attributes["modifiers"]
            elmods = element.attributes["modifiers"]
            #At first it seemed like this would be an important check, but a pre-req
            #could theoretically allocate pointer and then it could be passed into
            #a routine that just needed a value.
            #if ("pointer" in selfmods and "pointer" not in elmods) or \
            #   ("pointer" not in selfmods and "pointer" in elmods):
            #    return False
            if ("allocatable" in selfmods and "allocatable" not in elmods) or \
               ("allocatable" not in selfmods and "allocatable" in elmods):
                return False
            else:
                return True
        else:
            return False
        
    @property
    def attributes(self):
        """Returns a dictionary of attributes for this global variable declaration."""
        if isinstance(self.element, DocElement):
            return self.element.attributes
        else:
            #This handles the case that we are constructing from an XML element rather
            #than a DocElement
            return self.element.attrib

    def definition(self):
        """Returns the fortran declaration that defines a global variable.

        :arg altname: specifies an alternative name to use for this variable; all
          other variable properties remain the same.
        """
        if self.ignore:
            return "  ! Variable '{}' was set to be ignored.".format(self.attributes["name"])
        
        result = []
        if "type" not in self.attributes or "name" not in self.attributes:
            msg.err("required variable for execution missing some attributes." + 
                    " {}".format(self.attributes))
            exit(1)

        result.append(self.attributes["type"])
        if "kind" in self.attributes and self.attributes["kind"] is not None:
            if self.attributes["type"] == "character" and "len=*" in self.attributes["kind"]:
                msg.warn("Assumed length character array defaulting to kind of 'len=100'. "
                         "Use a <global> tag to override this default.")
                result.append("(100)")
            else:
                result.append("({})".format(self.attributes["kind"]))

        dimalloc = False
        if "dimensions" in self.attributes:
            normal = True
            if isinstance(self.element, DocElement) and self.element.doctype == "AUTOPARAM":
                #We used to force the user to explicitly declare a <global> tag if the
                #dimensions were fixed; however, we should be able to predict whether it
                #is necessary: if the dimensions only include digits etc. then it can be
                #declared as is; otherwise, we replace it with ':'.
                complexdim=[re.match("^\d+", d) is None for d in re.split(',\s*', self.attributes["dimensions"])]
                if any(complexdim):
                    normal = False
            if normal:
                dimstr = "({})".format(self.attributes["dimensions"])
            else:
                dimstr = "({})".format(','.join([':' for i in range(self.attributes["D"])]))
                dimalloc = True
        else:
            dimstr = ""
        if "modifiers" in self.attributes:
            mods = re.split(",\s*", self.attributes["modifiers"])
            smods = self._clean_mods(mods)
            if dimalloc and "allocatable" not in smods and "pointer" not in smods:
                smods += (", " if len(smods) > 0 else "") + "allocatable"
        else:
            smods = ""

        if smods != "":
            result.append(", " + smods)
        result.append(" :: ")
        result.append(self.attributes["name"])

        if dimstr != "":
            result.append(dimstr)
        if "default" in self.attributes:
            result.append(" ={}".format(self.attributes["default"]))
            
        return "  " + "".join(result)

    def _clean_mods(self, mods):
        """Returns only those modifiers that belong in a definition statement."""
        result = []
        if "pointer" in mods:
            result.append("pointer")
        if "allocatable" in mods:
            result.append("allocatable")

        return ", ".join(result)

    def initialization(self):
        """Returns the fortran declaration to initialize the global variable's value."""
        #See which of the options for initializing the variable were specified
        #in the docstring.
        if "default" in self.attributes:
            result = None #We handle default values at the definition level.
        elif "mode" in self.attributes:
            result = "call " + self.attributes["name"] + \
                     "%set_mode({})".format(self.attributes["mode"])
        else:
            result = None
        #else: result = "! {} had no initialization value specified".format(self.attributes["name"])

        if result is not None:
            return "  " + result  
        else:
            return result

    def range_check(self):
        """Determines if this variable has a range check specified and
        writes the fortran code to compute the check."""
        if "range" in self.attributes:
            #The fortpy module has an interface to handle any sorts of types we throw at
            #it, so we just need to make the call.
            rangecall = "call fpyin_range({},{},{})"
            #TODO, we need to create the range variables in the program and assign there
            #values if range was specified. 2D arrays need a 2D vector of ranges for
            #each dimension. Need to adjust the fortpy module to accomodate. Call the vars
            #by paramname_min etc.
            #We don't need to use fpyin_range b/c it just calls any/all built-in values.
        else:
            return None

class AutoClasser(object):
    """Class for generating read/save calls in Fortran for complex
    user-derived types. Ragged arrays are *not* supported.

    :arg variable: an instance of ValueElement to auto-class.
    :arg folder: the path to the *parent* folder that houses the 'varfile' folder with
      the contents of the variable.
    :arg varfile: the name of the folder to create (relative to 'folder') to contain
      the files for the variable.
    :arg cases: the TestSpecification instance's cases that will run.
    :arg coderoot: the full path to the code folder that houses that executable being
      unit tested.
    """
    def __init__(self, variable, folder, varfile, cases, coderoot=None, generic=False):
        self.variable = variable
        self.folder = folder
        self.varfile = varfile
        self.cases = cases
        self.coderoot = coderoot
        self._ofolder = folder
        """The original folder specification before it got overwritten with the
        relative one in the fortran program.
        """
        self.generic = generic
        """When true, instead of outputing the specific variable name for save and read
        operations, the generic 'variable' will be used as the name.
        """
        
        self._wcode = None
        """The fortran code to *save* the variable to its folder."""
        self._rcode = None
        """The fortran code to *read* a variable from folder so long as the
        files are named using the correct convention.
        """
        self._vars = []
        """The variables needed to perform the read/write operations; essentially
        the same list whether it is input/output.
        """
        self._is_vcoded = None
        """Specifies whether this auto-classer has already output code to define
        the variables (for cases where the same variable has its value read,
        modified and then saved using auto-class)."""
        self.tree_order = None
        """The order in which the tree's dimensionalities will be saved in the
        .fortpy.analysis file.
        """
        self._read_dims = None
        """Has the actual file dimensionality for each variable in self.tree_order as
        extracted from the first auto-class folder that gets scanned. There should be
        consistency between the folders to match the strongly-typed variables."""

    def set_folder(self, folder):
        """Sets the folder for the auto-class relative to the *code root*.
        """
        if folder is None:
            raise ValueError("Can't set folder to 'None'.")

        self._ofolder = folder
        from fortpy.tramp import coderelpath
        relpath = coderelpath(self.coderoot, folder).format(r"'//trim(adjustl(fpy_case))//'")
        if self.coderoot in relpath:
            self.folder = relpath.replace(self.coderoot, "fpy_coderoot//'")
        else:
            self.folder = relpath
        
    @property
    def acroot(self):
        """The folder path to the auto-class target for read/write, *relative*
        to 1) code directory for reads; 2) test directory for writes.
        """
        if self.varfile != "":
            return "{}/{}/".format(self.folder, self.varfile)
        else:
            return self.folder + '/'

    def abspath(self, root):
        """Returns the absolute path to the auto-class folder relative to the root
        directory specified: 1) code root for reads; 2) test execution directory
        for writes.
        """
        from fortpy.tramp import coderelpath
        return coderelpath(root, self.acroot)

    def xcode(self, lines, position, spacer, write=True, first=True, testsource=None):
        """Generates the code to save/read a variable from its auto-class folder.
        This version uses the `auxsave` and `auxread` interfaces from `fpy_auxiliary`
        instead of generating all the save and read statements directly.
        """
        #Everything is handled internally by the auxiliary module. We don't have
        #any additional variables/initialization outside of the ones handled
        #globally by the <assignment> and <target> objects.
        acroot = testsource if testsource is not None else self.acroot
        if position == "vars" and "fpy_coderoot" in acroot:
            rootstr = "{}character({}), parameter :: fpy_coderoot = '{}'"
            lines.append(rootstr.format(spacer, len(self.coderoot), self.coderoot))
        if position == "save" and write:
            if type(self.variable).__name__ == "ValueElement":
                fstr = "{}call auxsave({}, '{}')"
            elif type(self.variable).__name__ == "Function":
                fstr = "{}call auxsave({}_fpy, '{}')"
            lines.append(fstr.format(spacer, self.variable.name, acroot))
        elif position == "assign" and not write:
            fstr = "{}call auxread({}, '{}')"
            lines.append(fstr.format(spacer, self.variable.name, acroot))
    
    def _scan_folder(self, spath):
        """Creates a tree representation of the specified folder under the assumption
        that it contains file data for an auto-class read.

        :arg spath: the full path to the source folder to scan.
        """
        from os import walk, path
        import re
        sources = []
        for (dirpath, dirnames, filenames) in walk(spath):
            sources.extend(filenames)
            break
            
        tree = {}
        for source in sources:
            if source[0] != "_":
                #Ignore files that don't match the convention.
                continue

            parts = source.split("-")
            branch = None

            for i, p in enumerate(parts):
                if not re.match("[\d.]+", p):
                    #This is a variable name, *not* a dimension specification
                    if branch is None:
                        if p not in tree:
                            tree[p] = {"0": []}
                        branch = tree[p]
                    else:
                        if p not in branch:
                            branch[p] = {"0": []}
                        branch = branch[p]
                else:
                    branch["0"].append(p)

        return tree

    def _prep_single(self, relpath):
        """Performs an analysis and saves the necessary files for the folder
        at the specified path.
        """
        tree = self._scan_folder(relpath)
        _read_dims = self._analyze_tree(tree)
        if self._read_dims is None:
            self._read_dims = _read_dims
        if self.tree_order is None:
            #We only need to save the order for the first folder that we encounter.
            self.tree_order = []
            saveorder = True
        else:
            saveorder = False

        def _write_val(f, value):
            if isinstance(value, list):
                f.write("{}\n".format(' '.join(map(str, value))))
            else:
                f.write("{}\n".format(value))
            
        anpath = path.join(relpath, ".fortpy.analysis")
        with open(anpath, 'w') as f:
            if saveorder:
                for key, value in _read_dims.items():
                    self.tree_order.append(key)
                    _write_val(f, value)
                with open(path.join(relpath, ".fortpy.tree"), 'w') as g:
                    g.write('\n'.join(self.tree_order))    
            else:
                for key in self.tree_order:
                    if key in _read_dims:
                        _write_val(f, _read_dims[key])
                    else:
                        msg.warn("Folder '{}' is missing member variable data '{}'".format(relpath, key))
                    
    def prep_read(self):
        """Sets this auto-classer up to handle reading-in from a folder
        by analyzing the contents of the folder.
        """
        def _load_tree_order(relpath):
            """Attempts to load the tree order from file for the specified auto-class
            directory.
            """
            if self.tree_order is None:
                trpath = path.join(relpath, ".fortpy.tree")
                if path.isfile(trpath):
                    with open(trpath) as f:
                        self.tree_order = f.read().split('\n')

        def _load_read_dims(relpath):
            """Attempts to load the variable dimensionality from file for the specified
            auto-class directory.
            """
            if self.tree_order is None:
                return
            
            if self._read_dims is None:
                anpath = path.join(relpath, ".fortpy.analysis")
                if path.isfile(anpath):
                    self._read_dims = {}
                    with open(anpath) as f:
                        lines = f.readlines()
                    for itree, v in enumerate(self.tree_order):
                        self._read_dims[v] = list(map(int, lines[itree].split()))
                        if len(self._read_dims[v]) == 1 and self._read_dims[v][0] == 0:
                            self._read_dims[v] = self._read_dims[v][0]
                        
        #For handling multiple cases, it gets a little more complex; we need to
        #scan every folder that matches a case. All the folders need to use the
        #same tree order for consistency with the fortran driver.
        from fortpy.tramp import coderelpath
        from os import path
        if self.cases is not None and "{}" in self._ofolder:
            proclist = []
            for case in self.cases:
                relpath = coderelpath(self.coderoot, self._ofolder.format(case))
                if path.isdir(relpath):
                    anpath = path.join(relpath, ".fortpy.analysis")
                    if not path.isfile(anpath):
                        proclist.append(relpath)

                    #One of the folders will already have a tree order if we have ever
                    #processed this variable's folders before.
                    _load_tree_order(relpath)
                    _load_read_dims(relpath)

            for relpath in proclist:
                self._prep_single(relpath)
        else:
            relpath = coderelpath(self.coderoot, self._ofolder)
            _load_tree_order(relpath)
            _load_read_dims(relpath)
            self._prep_single(relpath)
                    
    def _analyze_tree(self, tree=None, prefix=""):
        result = {}            
        for key, value in tree.items():
            if key != "0":
                if prefix != "":
                    xkey = "{}.{}".format(prefix, key)
                else:
                    xkey = key
                result[xkey] = self._analyze_branch(key, value)
                result.update(self._analyze_tree(value, xkey))

        return result
                
    def _analyze_branch(self, memname, branch):
        if len(branch["0"]) == 0:
            return 0
        else:
            #Make sure that the dimensionality is consistent across the files.
            bdim = None
            D0 = None
            for leaf in branch["0"]:
                D = tuple(map(int, leaf.split('.')))
                if bdim is None:
                    bdim = [0]*len(D)
                    D0 = D

                if len(D0) != len(D):
                    raise ValueError("Inconsistent dimensionality for member "
                                     "'{}': {}".format(memname, branch["0"]))

                #Now we just need to extract the maximum value in each dimension.
                for i in range(len(D0)):
                    if D[i] > bdim[i]:
                        bdim[i] = D[i]

            return bdim
        
    def code(self, lines, position, spacer, write=True, first=True):
        """Appends the code required to save/read the value of this auto-class to
        file-folder if the specified position is appropriate.
        """
        #We assume that whatever is calling the auto-classer has already figured
        #out whether the position is correct for the read/write code.
        #Since the procesing generates code for both the read/write *and* the vars, we
        #need to run it irrespective.
        if self._wcode is None and write:
            self._process(spacer, True)
        elif self._rcode is None and not write:
            self._process(spacer, False)

        if position not in ["vars", "init"]:
            if write:
                lines.extend(self._wcode)
            else:
                lines.extend(self._rcode)        

        #For declaring the variables, we only do it once per application, even if we
        #are performing both a read and a write.
        if position == "vars":
            if self._is_vcoded is None:
                lines.append("{}!V: {} auto-class support variables".format(spacer, self.variable.name))
                for v in self._vars:
                    if "_acps" in v:
                        lines.append("{}character(100) :: {}".format(spacer, v))
                    elif "_wrote" in v:
                        lines.append("{}logical :: {}".format(spacer, v))
                    else:
                        lines.append("{}integer :: {}".format(spacer, v))
                self._is_vcoded = []

            if not write and "read" not in self._is_vcoded:
                #Also append on the ragged-array for handling the dimensionality of
                #the separate members in the auto-class folder.
                lines.append("{}class(fpy_vararray), allocatable :: {}_fpy_acdims(:)".format(spacer, self.variable.name))
                if first:
                    lines.append("{}integer :: ifpy_ac".format(spacer))
                    rootstr = "{}character({}), parameter :: fpy_coderoot = '{}'"
                    lines.append(rootstr.format(spacer, len(self.coderoot), self.coderoot))
                self._is_vcoded.append("read")
            if write:
                self._is_vcoded.append("write")

        if position == "init" and not write:
            #Read in the contents of the fortpy analysis into the ragged array
            quote = "'" if "fpy_coderoot" not in self.acroot else ""
            fstr = "{}call autoclass_analyze({}{}.fortpy.analysis', {}_fpy_acdims)"
            lines.append(fstr.format(spacer, quote, self.acroot, self.variable.name))
        
    def _process_autovar(self, variable, context, filecontext, treecontext, spacer,
                         write, depth=0):
        """Creates a system to recursively save the contents of a single variable
        that is a user-derived type with multiple members.

        :arg variable: an instance of ValueElement representing the type of the variable
          to save.
        :arg context: for nested derived types like 'type1%type2', variable would be 'type2'
          and the context would be 'type1%'.
        :arg vslice: if the derived type is part of an array, the slice of the array to
          generate the pysave for. This will be a list of loop variables.
        :arg filecontext: the context of the *filename* that the variable values are being
          saved to.
        :arg spacer: whitespace to prepend to each line of code generated.
        """
        result = []
        spacing = spacer
        #Find the effective dimension; since fortpy can save multi-D arrays for simple
        #data types, those don't have any effective dimensions. 
        effD = variable.D if variable.is_custom else 0
        prefix = context.replace("%", "_") + variable.name
        loopvars = ["{}_ac{}".format(prefix, i+1) for i in range(effD)]
        lslice = ', '.join(loopvars)
        pslist = "{}_acps{}".format(prefix, depth+1)

        if treecontext == "" and variable.name.lower() == self.variable.name.lower():
            ntreecontext = "_"
        else:
            ntreecontext = "{}.{}".format(treecontext, variable.name.lower())

        if not write:
            #Next, we see which line in the analysis file has the dimensionality
            #information for us to set the limits on the loops.
            if ntreecontext in self.tree_order:
                itree = self.tree_order.index(ntreecontext) + 1
                limit = "{}_fpy_acdims({})%items({{}})".format(self.variable.name, itree)
            else:
                limit = "0"
        else:
            limit = "size({}{}, {})"

        for i, loopvar in enumerate(loopvars):
            if write:
                limit = limit.format(context, variable.name, i+1)
            else:
                limit = limit.format(i+1)
            result.append("{}do {}=1, {}".format(spacing, loopvar, limit))
            spacing = spacing + '  '

        if effD > 0:
            indices = "(/ {} /)".format(lslice)
            result.append("{}call fpy_period_join_indices({}, ".format(spacing, pslist) +
                          "{}, {})".format(indices, effD))
            loopvars.append(pslist)

        #Since the folder that all the files gets saved in is already named with the
        #highest-level variable's name, we don't need to add that to the path.
        if self.generic and variable.name.lower() == self.variable.name.lower() and depth==0:
            if (variable.kind is not None and self.variable.kind is not None and
                variable.kind.lower() == self.variable.kind.lower()):
                filevarname = ""
            else:
                filevarname = "'"
        elif variable.name.lower() == self.variable.name.lower() and depth==0:
            filevarname = "_"
        else:
            filevarname = variable.name

        if (variable.is_custom and depth==0):
            custype = variable.customtype
            for member in custype.members.values():
                if effD > 0:
                    if self.generic:
                        ncontext = "{}variable({})%".format(context, lslice)
                    else:
                        ncontext = "{}{}({})%".format(context, variable.name, lslice)
                    nfilecontext = "{}{}-'//adjustl(trim({}))//'-".format(filecontext, filevarname, pslist)
                else:
                    if self.generic:
                        ncontext = "{}variable%".format(context)
                    else:
                        ncontext = "{}{}%".format(context, variable.name)
                    nfilecontext = "{}{}-".format(filecontext, filevarname)

                icode, ivars = self._process_autovar(member, ncontext, nfilecontext, ntreecontext,
                                                     spacing, write, depth+1)
                result.extend(icode)
                loopvars.extend(ivars)
        else:
            #Yay, we are finally down to the standard type level! Just call it with pysave
            if effD > 0:
                varname = "{}{}({})".format(context, variable.name, lslice)
                filename = "{}{}-'//trim(adjustl({}))".format(filecontext, filevarname, pslist)
            else:
                varname = "{}{}".format(context, variable.name)
                filename = "{}{}'".format(filecontext, filevarname)

            if write:
                if "allocatable" in variable.modifiers:
                    result.append("{}if (allocated({}{})) then".format(spacing, context, variable.name))
                    spacing += "  "
                if "pointer" in variable.modifiers:
                    result.append("{}if (associated({}{})) then".format(spacing, context, variable.name))
                    spacing += "  "
                if (variable.is_custom and 
                    self.variable.customtype.name.lower() == variable.kind.lower()):
                    rl = []
                    rl.append("nprefix = {}".format(filename))
                    rl.append("call auxsave_{}0d_({}, nprefix, stack, wrote)".format(self.variable.customtype.name,
                                                                                 varname))
                    result.extend([spacing + l for l in rl])
                elif variable.is_custom:
                    if "private contents" not in variable.customtype.modifiers:
                        loopvars.append("{}_wrote".format(variable.name))
                        result.append("{}call auxsave({}, {}, .true., {}_wrote)".format(spacing, varname, filename, variable.name))
                        result.append("{}wrote = wrote .or. {}_wrote".format(spacing, variable.name))
                        result.append("{}if (.not. {}_wrote) then".format(spacing, variable.name))
                        result.append("{}  call pysave(.false., {}//'-.fpy.blank')".format(spacing, filename))
                        result.append("{}end if".format(spacing))
                else:
                    if filename[-2]=="'":
                        print(filename, filecontext, filevarname)
                    result.append("{}call pysave({}, {})".format(spacing, varname, filename))
                    result.append("{}wrote = .true.".format(spacing))

                if "allocatable" in variable.modifiers or "pointer" in variable.modifiers:
                    spacing = spacing[0:-2]
                    result.append("{}end if".format(spacing))
            else:
                #We need to see whether we have allocatable, pointer or fixed dimension variable
                #and call the relevant interface.
                fstr = "{}call fpy_read{}({}, '#', {})"
                if "pointer" in variable.modifiers:
                    result.append(fstr.format(spacing, "_p", filename, varname))
                elif ("allocatable" not in variable.modifiers and variable.D > 0 and
                      not re.match("[\w]+", variable.dimension)):
                    result.append(fstr.format(spacing, "_f", filename, varname))
                else:
                    result.append(fstr.format(spacing, "", filename, varname))                                

        for i in range(effD):
            spacing = spacing[0:-2]
            result.append("{}end do".format(spacing))
        return (result, loopvars)
                
    def _process(self, spacer, write):
        """Processes a single variable as an auto-class for reading/writing with a single
        folder full of files for its members.
        """
        if not write and self.tree_order is None:
            self.prep_read()

        quote = "'" if "fpy_coderoot" not in self.acroot and not self.generic else ""
        root = (self.acroot if not self.generic else
                ("prefix//'" if self.variable.customtype.recursive else "folder//'"))
        icode, ivars = self._process_autovar(self.variable, "", quote + root, "", spacer, write)
        if write:
            self._wcode = icode
        else:
            self._rcode = []
            #We need to handle all the allocations for the members based on the actual size of the data
            #unless they are standard types that will be allocated by the fortpy interfaces.
            flines = []
            acdims = "{}_fpy_acdims".format(self.variable.name)
            flines.append("do ifpy_ac=1, size({}, 1)".format(acdims))
            for itree, v in enumerate(self.tree_order):
                if isinstance(self._read_dims[v], list):
                    varname = v.replace("_", self.variable.name).replace('.', '%')
                    dslice = ','.join(["{}(ifpy_ac)%items({})".format(acdims, i+1)
                                       for i in range(len(self._read_dims[v]))])
                    alloc = "{}({})".format(varname, dslice)
                    flines.append("  if ({}(ifpy_ac)%items(1) .ne. 0) allocate({})".format(acdims, alloc))
            flines.append("end do")
            self._rcode.extend([spacer + l for l in flines])
            self._rcode.extend(icode)

        if len(self._vars) == 0 and len(ivars) > 0:
            self._vars = ivars
        
class AssignmentValue(object):
    """Represents a value that can be assigned to the variable represented
    by an Assignment instance.

    :attr xml: the xml <value> tag that this object was created from.
    :attr parent: the instance of Assignment that owns this value.
    :attr identifier: the unique identifier of the this value.
    :attr folder: the folder path relative to the code folder where the
      file for setting a variable value is found.
    :attr filename: the name of the file in the source folder.
    :attr rename: the new name the file should have in the execution folder.
    :attr constant: if specified, the constant value that the variable will
      be set to during the unit test.
    :attr embedded: the name of an embedded procedure in a derived type to
      run. An error is thrown if the variable is not a derived type.
    :attr function: exact fortran code for a function call to set the vars value.
    :attr repeats: if the variable will have its value set as part of a 
      constant-input mode program, this specifies whether it should keep being
      set during each iteration of the main method being tested.
    :attr prereq: if not None/False, the embedded method on a derived type will be treated
      along with all its prereqs as part of the execution chain. This value should be the
      test identifier of the embedded subroutine's test spec to use.
    :attr D: the dimensionality of the data in the file that is setting the variable value.
    :attr ragged: when true, the lines in a 2D array file are treated individually with
      some other (embedded or function) value assignment.
    :attr dtype: for ragged array files, the data type of the values in each row.
    :attr kind: the kind of each value in the data rows of the ragged array.
    """
    def __init__(self, xmltag, parent):
        """Initializes the assignment value using the <value> tag."""
        self.xml = xmltag
        self.parent = parent
        self.identifier = None
        self.folder = None
        self.filename = None
        self.rename = None
        self.constant = None
        self.embedded = None
        self.function = None
        self.repeats = None
        self.prereqs = None
        self.paramlist = None
        self.ragged = False
        self.dtype = None
        self.kind = None
        self.commentchar = "#"
        self.member = None
        """If this value should be assigned to the member of a derived type,
        then this is the name of that member.
        """
        self.suffix = None
        """Overrides the default period joined list of loop variables for <part>
        or ragged assignment from files.
        """
        self.autoclass = False
        """Specifies whether the value comes from an auto-class structured folder.
        """
        self.testsource = None
        """Specifies the module.executable.parameter name whose output from a previous
        unit test run will be used as input for this value.
        """
        
        self._derived_type = None
        self._codes = {
            "vars": self._code_vars,
            "init": self._code_init,
            "assign": self._code_assign,
            "after": self._code_after,
            "before": self._code_before
        }

        if xmltag is not None:
            self._parse_xml()

    @property
    def iid(self):
        """Returns the globally unique identifier for this value."""
        return "{}_{}".format(self.parent.name, self.identifier)
       
    @property
    def xname(self):
        """Returns the name of the file to use taking renaming into account."""
        if self.rename is not None:
            return self.rename
        else:
            return self.filename

    def livefile(self, testroot, case):
        """Returns the full path to the file that should be *created* in the
        specified test directory as part of this assignment. This is only
        useful for file-based assignments.
        """
        if self.rename is not None:
            target = path.join(testroot, self.rename)
        else:
            target = path.join(testroot, self.filename.format(case))
        return target
        
    def copy(self, coderoot, testroot, case, compiler, stagedir=None):
        """Copies the input files needed for this value to set a variable.

        :arg coderoot: the full path to folder with the code files.
        :arg testroot: the full path the folder where the test is running.
        :arg case: the case id for multi-case testing.
        :arg compiler: the name of the compiler being used for the unit tests.
        :arg stagedir: the path to the staging directory to use (override).
        """
        #We only want to do the copy if we are assigning a value from a file.
        if (self.folder is not None and self.filename is not None and not self.autoclass
            and self.testsource is None):
            from fortpy.tramp import coderelpath
            from fortpy.testing.compilers import replace
            from fortpy.utility import symlink
            relpath = coderelpath(coderoot, self.folder)
            source = replace(path.join(relpath, self.filename.format(case)), compiler)
            source = replace(source, compiler, True)
            
            #For the cases where multiple files specify the values for different
            #parts of an array, the <value> file specification will have a wildcard
            #at the position where the array index id will go. We just copy *all* the
            #input files over that match that definition.
            if "*" in source:
                import glob
                for filename in glob.glob(source):
                    suffix = ".{}".format(case)
                    if filename[-len(suffix)::] == suffix:
                        target = path.join(testroot, filename[0:len(filename)-len(suffix)].split("/")[-1])
                        symlink(filename, target)
            else:
                target = self.livefile(testroot, case)
                symlink(source, target)

        if (self.testsource):
            #We don't want to duplicate lots of files all over the system. If they already
            #exist in an output folder, we should just tell Fortran to get the files from
            #there. For the copy operation, we set the *path* to the actual data in a file
            #with the variable name and 'ts': 'fpy_var_ts'.
            target = self.livefile(testroot, case)
            datpath = self.testsource.path(coderoot, case, compiler, stagedir)
            if path.isfile(target):
                with open(target) as f:
                    rewrite = datpath == f.readlines()[1]
            else:
                rewrite = True

            if rewrite:
                with open(target, 'w') as f:
                    f.write('#<fortpy mode="testsource" copy="false" />\n')
                    f.write(datpath)

    def check_prereqs(self, finder):
        """Checks whether this value element requires an embedded method to be
        added to the prereq chain of the method writer."""
        if self.embedded is not None:
            #The embedded method could use a pointer so that the method we are
            #calling is *not* the name of the actual method. We need to find the
            #instance of the TypeExecutable to locate its actual target.
            var = self.parent.variable
            target = var.kind
            module = finder.module
            self._derived_type, typemod = self.parent.parser.tree_find(target.lower(), module, "types")

            if self._derived_type is None:
                raise ValueError("The type for embedded method {} cannot be found.".format(self.embedded))

            if self.prereqs:
                typex = self._derived_type.executables[self.embedded.lower()]
                if self.typex.pointsto is not None:
                    key = "{}.{}".format(self._derived_type.module, typex.pointsto)
                else:
                    key = "{}.{}".format(self._derived_type.module, typex.name)

                if self.prereqs == True:
                    #Let the MethodFinder use the first testing group test.
                    finder.add_prereq(key, self.parent.element, None)
                else:
                    finder.add_prereq(key, self.parent.element, self.prereqs)
 
    def code(self, lines, position, spacer, slices=None, first=False):
        """Appends the code lines to initialize the parent variable.

        :arg lines: the list of strings to append to.
        :arg position: one of ['vars', 'init', 'assign', 'before', 'after'] indicating 
          where in the fortran program the code will be appended.
        :arg slices: for array value assignments, the specific indices to assign values
          to. Tuple of (slice string, [loopvars]).
        :arg first: if multiple value tags are being executed by the calling parent, this
          should be true only for the *first* <value> tag processed in the list.
        """
        if self.parent.variable is not None:
            if position in ["before", "assign"]:
                self._codes[position](lines, spacer, slices, first)
            else:
                self._codes[position](lines, spacer, slices)
        else:
            raise ValueError("Trying to assign a value to an unknown variable {}.".format(self.iid))

    def _code_before(self, lines, spacer, slices, first):
        """Calls _code_setvar() if we *are* in repeat mode."""
        if self.repeats:
            self._code_setvar(lines, spacer, slices, first)

    def _code_assign(self, lines, spacer, slices, first):
        """Calls _code_setvar() if we are *not* in repeat mode."""
        if not self.repeats:
            self._code_setvar(lines, spacer, slices, first)

    def _check_autoclass(self):
        """Makes sure that the autoclass for implementing this variable assignment
        exists in the parent test specification.
        """
        varkey = self.parent.name.lower()
        acdict = self.parent.testspec.autoclasses
        if varkey not in acdict:
            acdict[varkey] = AutoClasser(self.parent.variable, self.folder, "",
                                         self.parent.testspec.cases, self.parent.group.coderoot)
        else:
            acdict[varkey].varfile = ""
        acdict[varkey].set_folder(self.folder)
        return acdict[varkey]
            
    def _code_setvar(self, lines, spacer, slices, first):
        """Appends code for assigning the value of the parent variable using
        this value specification.

        :arg slices: for array value assignments, the specific indices to assign values
          to.
        :arg first: if multiple value tags are being executed by the calling parent, this
          should be true only for the *first* <value> tag processed in the list.
        """
        if self.autoclass:
            aclass = self._check_autoclass()
            if self.testsource is not None:
                testpath = self.testsource.filevar
            else:
                testpath = None
            aclass.xcode(lines, "assign", spacer, False, first, testpath)
        elif self.filename is None and self.testsource is None:
            #We handle these each individually when no file assignment is present.
            #The file assigner may handle the embedded and function calls at the right
            #spot when parts are being assigned.
            if self.constant is not None:
                self._code_setvar_value(lines, spacer, self.constant, slices)
            elif self.function is not None:
                self._code_setvar_value(lines, spacer, self.function, slices)
            elif self.embedded is not None:
                self._code_embedded(lines, spacer, slices)
        else:
            self._code_file(lines, spacer, slices, first)

    def _code_setvar_value(self, lines, spacer, value, slices=None):
        """Sets the value of the variable primitively to the specified value."""
        if self.member is None:
            varname = self.parent.name
        else:
            varname = "{}%{}".format(self.parent.name, self.member)
            
        if slices is None:
            lines.append("{}{} = {}".format(spacer, varname, value))
        else:
            lines.append("{}{}({}) = {}".format(spacer, varname, slices[0], value))
    
    def _code_embedded(self, lines, spacer, slices=None, varname=None):
        """Appends code for calling an embedded method in a derived type,
        optionally including all its dependencies.

        :arg varname: in the mode where individual files are being used to set the values
          of different parts of an array of embedded types, 'varname' is the name of the 
          temporary variable created to hold the contents of the file. It will be deallocated
          after the embedded subroutine is called.
        """
        if not self.prereqs and self._derived_type is not None:
            #This is a simple exercise in calling the embedded method. We just 
            #need to track down the parameters list.
            target = self._derived_type.executables[self.embedded].target
            if self.paramlist is None:
                #This should be easy, but the compiler automatically passes in the
                #derived type instance as the first parameter in the method; since
                #that still shows up in the executable parameter list, we need to
                #filter it out based on its kind. We assume that it would be first in
                #the list and have the same kind as the derived type.
                orig_params = list(target.ordered_parameters)
                if orig_params[0].kind == self._derived_type.name:
                    del orig_params[0]
                params = ", ".join([p.name for p in orig_params])
            else:
                params = self.paramlist

            if type(target).__name__  == "Subroutine":
                call = "call "
            else:
                raise ValueError("Embedded type initialization routines must be"
                                 " subroutines, functions aren't supported.")

            if slices is None:
                lines.append("{}{}{}%{}({})".format(spacer, call, self.parent.name, 
                                                    self.embedded, params))
            else:
                i = 1
                while "$" in params:
                    #The developer is using some existing array and wants the loop variable
                    #names substituted for some of the paramaters.
                    params = params.replace("${}".format(i), slices[1][i])
                    i += 1

                #Handle the embedded via file possibility.
                if varname is not None:
                    params = params.replace("@file", varname)

                lines.append("{}{}{}({})%{}({})".format(spacer, call, self.parent.name, 
                                                        slices[0], self.embedded, params))
        #if it does have pre-reqs, they will be handled by the method writer and we
        #don't have to worry about it.

    def _code_file(self, lines, spacer, slices, first):
        """Appends code for initializing the value of a variable that is a
        scalar, vector or 2D matrix from a file.

        :arg first: specifies that this file assignment is the first of a list
          of <value> tags that are siblings in the same context.
        """
        if self.filename is None and self.testsource is None:
            return
        
        flines = []
        #If we have slices being set from file values, we will have a set
        #of files that match a wildcard pattern. They will all have been copied
        #into the test directory. However, we need to replace the wildcard in
        #the filename at *runtime* with the current value of the loop variable.
        if "*" in self.xname and self.suffix != False:
            if self.suffix is None:
                lslice = slices[1]
            else:
                lslice = []
                for i in range(len(slices[1])):
                    if "${}".format(i+1) in self.suffix:
                        lslice.append(slices[1][i])
                        
            indices = "(/ {} /)".format(', '.join(lslice))
            flines.append("call fpy_period_join_indices(" +
                          "{}_pslist, {}, {})".format(self.iid, indices, len(lslice)))
            if self.xname[-1] != "*":
                #This makes sure that if the wildcard is in the middle of the filename
                #we still get it right.
                pre, post = self.xname.split("*")
                rtname = '"{}"//{}_pslist//"{}"'.format(pre, self.iid, post)
            else:
                rtname = '"{}"//{}_pslist'.format(self.xname[:len(self.xname)-1], self.iid)
        else:
            if self.testsource is not None:
                rtname = self.testsource.filevar
            else:
                rtname = "'{}'".format(self.xname)

        #There is also the case where we want to use a single file, but with ragged data
        #lengths on each line.
        if self.ragged:
            ragvar = "{}_rag".format(self.iid)
            if slices is None:
                slices = (ragvar, [ragvar])
            else:
                #Make sure that we have at least one free dimension available for the ragged
                #list in the file.
                if slices[0][-1] == ":":
                    slices[0][-1] = ragvar
                    slices[1].append(ragvar)
                else:
                    raise ValueError("The 'ragged' option can only be used when there is "
                                     "a spare dimension on the array to assign the ragged "
                                     "file values to.")

        if slices is not None and self.member is None:
            varname = "{}_fvar".format(self.iid)
        elif self.member is not None:
            if slices is not None:
                varname = "{}({})%{}".format(self.parent.name, ','.join(slices[1]), self.member)
            else:
                varname = "{}%{}".format(self.parent.name, self.member)
        else:
            varname = self.parent.name

        #We are working with a vector or scalar. Check the dimensionality of
        #the actual variable and see if it needs to be allocated.
        if self.ragged:
            #For the ragged option, this is the only place that we handle it.
            flines.append("call fpy_linevalue_count({}, ".format(rtname) +
                          "'{0}'".format(self.commentchar) + 
                          ", {0}_nlines, {0}_nvalues)".format(self.iid))
            flines.append("allocate({}({}_nlines))".format(self.parent.name, self.iid))
            flines.append("call fpy_linevalue_count_all({}, ".format(rtname) +
                          "'{0}'".format(self.commentchar) + 
                          ", {0}_nlines, {0}_ragvals)".format(self.iid))
            flines.append("open(fpy_newunit({}_funit), ".format(self.iid) + 
                          "file={})".format(rtname))
            flines.append("do {0}_rag=1, {0}_nlines".format(self.iid))
            flines.append("  allocate({0}_fvar({0}_ragvals({0}_rag)))".format(self.iid))
            flines.append("  read({}_funit, *) {}".format(self.iid, varname))
                
        #Now handle the case where the input file fills a 2D variable.
        if not self.ragged:
            fmtstr = "call fpy_read{}({}, '{}', {})"
            if self.member is None:
                modifiers = self.parent.global_attr("modifiers", [])
                D = self.D
                dimensions = self.parent.global_attr("dimensions", "")
            elif type(self.parent.variable).__name__  == "ValueElement":
                #We need to get the code element of the *member* variable and look
                #at its modifiers instead of the parent.
                custype = self.parent.variable.customtype
                if custype is not None and self.member.lower() in custype.members:
                    memvar = custype.members[self.member.lower()]
                    modifiers = memvar.modifiers
                    D = memvar.D
                    dimensions = memvar.dimension
                else:
                    raise ValueError("The member '{}' is not part of user-derived".format(self.member) +
                                     " type '{}'.".format(self.parent.variable.kind))
            else:
                raise ValueError("The variable '{}' is not a user-derived type.".format(self.parent.name))

            if ("pointer" in modifiers and "fvar" not in varname):
                flines.append(fmtstr.format("_p", rtname, self.commentchar, varname))
            elif ("allocatable" not in modifiers and "fvar" not in varname and D > 0 and
                  not re.match("[A-Za-z]+", dimensions)):
                flines.append(fmtstr.format("_f", rtname, self.commentchar, varname))
            else:
                flines.append(fmtstr.format("", rtname, self.commentchar, varname))

        #Handle the mixture of embed/function and filename case.
        if slices is not None:
            if self.embedded is not None:
                self._code_embedded(flines, "", slices, varname)
            if self.function is not None:
                filefun = self.function.replace("@file", varname)
                self._code_setvar_value(lines, spacer, filefun, slices)

        if self.ragged:
            flines.append("  deallocate({0}_fvar)".format(self.iid))
            flines.append("end do")
            flines.append("close({}_funit)".format(self.iid))
        if self.D == 2 and slices is not None:
            flines.append("deallocate({0}_fvar)".format(self.iid))

        #Deallocate the variable for concatenating the loop ids to form the file name
        #if "*" in self.xname:
        #    flines.append("deallocate({0}_pslist)".format(self.iid))
        lines.extend([ spacer + l for l in flines])

    def _code_after(self, lines, spacer, slices=None, first=True):
        """Appends code for deallocating a variable that was assigned a value.
        This is useful so that we can reallocate it again in repeat mode.
        """
        if (self.repeats and self.parent.allocate and 
            ("allocatable" in self.parent.variable.modifiers or
             "pointer" in self.parent.variable.modifiers)):
            lines.append("{}deallocate({})".format(spacer, self.parent.name))

    def _code_init(self, lines, spacer, slices=None, first=True):
        """Appends code to initialize any variables we need for assignment
        operations later on.
        """
        #The test source initialization needs to be coded first since the autoclass
        #could *also* be specified at the same time. In that case we need the variable
        #with the full folder name to be available.
        if self.testsource:
            ts = self.testsource
            lines.append(spacer + "!Read the full path to the test-source output data to use.")
            lines.append("{}call fpy_read('{}', '#', {})".format(spacer, ts.filename, ts.filevar))

        if self.autoclass:
            #We need to analyze the contents of the folder that the variable
            #is getting its value from and write them to the .fortpy.analysis
            #file if it doesn't already exist. The fortran program will suck
            #that file's contents to see the sizes for allocating the arrays.
            aclass = self._check_autoclass()
            aclass.xcode(lines, "init", spacer, False)

    def _code_vars(self, lines, spacer, slices=None, first=True):
        """Appends lines to declare any variables we need for file read
        operations later on.
        """
        if self.testsource:
            vtext = "{}character(250) :: {} ! test-source path variable"
            lines.append(vtext.format(spacer, self.testsource.varname))
        #It is possible to have testsource and autoclass at the *same* time.
        if self.autoclass:
            aclass = self._check_autoclass()
            aclass.xcode(lines, "vars", spacer, False, first)
        elif self.filename is not None:
            if ("*" in self.xname or self.ragged) and self.member is None:
                if self.dtype is None:
                    raise ValueError("Wildcard file names and ragged array inputs both require "
                                     "attribute 'dtype' to be specified when attribute 'member' "
                                     "is not present.")
                if self.kind is None:
                    skind = ""
                else:
                    skind = "({})".format(self.kind)

                fdim = [":"]*self.D
                
                lines.append("{}!Vars for initializing variable {} ".format(spacer, self.parent.name) +
                         "from file {}".format(self.filename))
                lines.append("{}{}{}, allocatable".format(spacer, self.dtype, skind) + 
                             " :: {}_fvar({})".format(self.iid, ','.join(fdim)))

            if "*" in self.xname:
                #This hard-coded 100 assumes that if we have 7 dimensional array and each array
                #dimension has the maximum size of 4.2 billion for an integer, we would have
                #(10 characters (4.2 billion) + 1 (period))*7 = 77
                lines.append("{}character(100) :: {}_pslist".format(spacer, self.iid))
            if self.ragged:
                lines.append("{0}integer :: {1}_nlines, {1}_nvalues, {1}_funit".format(spacer, self.iid))
                lines.append("{}integer :: {}_rag".format(spacer, self.iid))
                lines.append("{}integer, allocatable :: {}_ragvals(:)".format(spacer, self.iid))

    @property
    def D(self):
        """Returns the dimensionality of the variable having its value changed."""
        #We have to do delayed assignment because we don't have a method finder available in the
        #parent assignment until a specific test specification is being setup.
        if "filedim" in self.xml.attrib:
            D = int(self.xml.attrib["filedim"])
            if D > 2:
                raise ValueError("Only 2D arrays are handled automatically via file assignment.")
            return D
        else:
            return self.parent.variable.D

    def _parse_xml(self):
        """Extracts the relevant attributes from the <value> tag."""
        if "identifier" in self.xml.attrib:
            self.identifier = self.xml.attrib["identifier"]
        else:
            #We use this for automatic assignment of value identifiers when there
            #is only one <value> tag and it is unambiguous. The ambiguity is checked
            #by the parent <assignment> tag.
            self.identifier = "default"
        if "member" in self.xml.attrib:
            self.member = self.xml.attrib["member"]

        if "folder" in self.xml.attrib:
            self.folder = self.xml.attrib["folder"]
        if "file" in self.xml.attrib:
            self.filename = self.xml.attrib["file"]
        if "constant" in self.xml.attrib:
            self.constant = self.xml.attrib["constant"]
        if "embedded" in self.xml.attrib:
            self.embedded = self.xml.attrib["embedded"]
        if "function" in self.xml.attrib:
            self.function = self.xml.attrib["function"]
        if "commentchar" in self.xml.attrib:
            self.commentchar = self.xml.attrib["commentchar"]
        if "paramlist" in self.xml.attrib:
            self.paramlist = self.xml.attrib["paramlist"]
        if "rename" in self.xml.attrib:
            self.rename = self.xml.attrib["rename"]
        if "repeats" in self.xml.attrib:
            self.repeats = self.xml.attrib["repeats"].lower() == "true"
        else:
            self.repeats = False
        if "prereqs" in self.xml.attrib and self.xml.attrib["prereqs"] != "false":
            self.prereqs = self.xml.attrib["prereqs"]

        if "ragged" in self.xml.attrib:
            self.ragged = self.xml.attrib["ragged"] == "true"
        if "dtype" in self.xml.attrib:
            self.dtype = self.xml.attrib["dtype"]
        if "kind" in self.xml.attrib:
            self.kind = self.xml.attrib["kind"]

        if "suffix" in self.xml.attrib:
            if "$" in self.xml.attrib["suffix"]:
                self.suffix = self.xml.attrib["suffix"].split(".")
            elif self.xml.attrib["suffix"].lower() == "false":
                self.suffix = False

        if "autoclass" in self.xml.attrib:
            self.autoclass = self.xml.attrib["autoclass"].lower() == "true"
        if "testsource" in self.xml.attrib:
            self.testsource = TestSource(self.xml.attrib["testsource"].lower(), self.parent.testspec)
            #If a test-source is specified, we need to make sure that the relevant
            #attributes from the source <target> tag get copied over to this <value>.
            self.autoclass = self.testsource.target.autoclass

class Condition(object):
    """Represents a single if, elseif or else block to execute."""
    def __init__(self, xmltag, parent):
        """Initializes using an <if>, <elseif> or <else> tag and an
        AssignmentConditional instance."""
        self.xml = xmltag
        self.parent = parent

        self.tag = self.xml.tag
        
        if "condition" in self.xml.attrib:
            self.condition = self.xml.attrib["condition"]
        elif self.tag != "else":
            raise ValueError("'condition' is a required attribute for <if> and"
                             " <elseif> tags.")

        if "value" in self.xml.attrib:
            self.value = re.split(",\s*", self.xml.attrib["value"].lower())
        else:
            raise ValueError("'value' is a required attribute of <if>, <elseif>"
                             " and <else> tags.")

    @property
    def repeats(self):
        """Returns true if any of the values this conditional references have
        the repeat attribute set to true.
        """
        for v in self.value:
            if v in self.parent.values and self.parent.values[v].repeats:
                return True
        else:
            return False  
        
    def code(self, lines, spacer):
        """Appends the code for this condition and its variable assignment."""
        if self.tag in ["if", "elseif"]:
            lines.append("{}{} ({}) then".format(spacer, self.tag, self.condition))
        else:
            lines.append("else")

        #Append the value assignment. To do this we have to look it up in the
        #grand-parents list of possible value assignments.
        for v in self.value:
            if v in self.parent.values:
                valobj = self.parent.values[v]
                valobj.code(lines, "assign", spacer + "  ")
            else:
                raise ValueError("Could not find value '{}' for condition.".format(v))
        
class AssignmentConditional(object):
    """Represents a series of logical tests to perform, each of which results
    in a different value being assigned to the variable."""
    def __init__(self, xmltag, parent):
        """Initializes using a <conditionals> tag and Assignment instance."""
        self.xml = xmltag
        self.parent = parent
        self.conditions = []

        for child in self.xml:
            if child.tag in ["if", "elseif", "else"]:
                self.conditions.append(Condition(child, self))

        self._repeats = None

    @property
    def values(self):
        """Returns the list of possible value objects from the parent
        Assignment object."""
        return self.parent.values

    @property
    def repeats(self):
        """Returns true if any of the values used by this condition are
        repeatable. In that case, the entire block needs to be repeated."""
        if self._repeats is None:
            for c in self.conditions:
                if c.repeats:
                    self._repeats = True
                    break
            else:
                self._repeats = False                

        return self._repeats

    def code(self, lines, position, spacer):
        """Appends the lines to form the entire conditional code block of
        variable assignments."""
        if (len(self.conditions) > 0 and 
            ((position == "before" and self.repeats) or
            (position == "assign" and not self.repeats))):

            #Make sure the first condition is an if and not something else
            if self.conditions[0].tag != "if":
                raise ValueError("The first condition in a logic block must be 'if'.")
            
            for c in self.conditions:
                c.code(lines, spacer)
            lines.append("{}end if".format(spacer))

class Part(object):
    """Represents a part assignment for a <part> tag to set specific dimensions on
    an array-valued variable.

    :attr identifier: the unique identifier for the <part>.
    :attr start: the index of the array to start on.
    :attr value: the value to assign to parts of the array specified by this part.
    :attr limits: a comma-separated list of integer indices to assign the value to
      or a range of min-max integer dimension specifiers. The loop generated by this
      part tag will run over these values.
    :attr parts: a dict of parts nested in this one.
    :attr depth: the depth of this part in a nested part specification.
    :attr assignment: the parent instance of Assignment for this part, even if its
      immediate parent is another Part instance.
    :attr isloop: for a limit specificaton that is a comma-separated list, we don't
      use a loop, only a repeated list of assignments.
    :attr loopid: the name of the loop variable for this part.
    :attr slice_list: list of the dimension specifcations for setting the array part that this
      object is supposed to set.
    :attr repeats: when operating in constant mode, should the part assignment happen before
      each call to the main function being unit tested.
    """
    def __init__(self, xml, parent):
        self.xml = xml
        self.parent = parent

        self.identifier = None
        self.start = 0
        self.value = None
        self.limits = None
        self.parts = {}
        self.isloop = False
        self.assignment = parent if isinstance(parent, Assignment) else parent.assignment
        self.depth = (1 if isinstance(parent, Assignment) else parent.depth + 1) + self.start
        self.loopid = "{}_fpy_{}".format(self.assignment.name, self.depth)
        self.repeats = False

        #We can't determine the remaining attributes until the XML has been parsed.
        self._parse_xml()
        self.slice_list = ([":"]*parent.variable.D if isinstance(parent, Assignment) 
                           else list(parent.slice_list))
        self.slice_list[self.depth-1] = self.loopid if self.isloop else "{}"
        self.slices = ','.join(self.slice_list)

        self._top_parts = None #Lazy assignment for property.

    @property
    def top_parts(self):
        """Returns the parent <part> instances from a set of nested <part> tags.
        """
        if self._top_parts is None:
            part = self
            self._top_parts = [self]
            while isinstance(part.parent, Part):
                part = part.parent
                self._top_parts.insert(0, part)
        return self._top_parts
    
    def _code_vars(self, lines, spacer):
        """Appends lines to declare any variables we need for file read
        operations later on.
        """
        if not self.isloop:
            return

        if self.depth == 1:
            lines.append("{}!V: {} array-part assignment loop".format(spacer, self.assignment.name))
        lines.append("{}integer :: {}".format(spacer, self.loopid))

    def code(self, lines, position, spacer):
        """Adds the code to make this part specification (and loop) work.
        """
        if position == "vars":
            self._code_vars(lines, spacer)
        elif position in ["before", "assign"]:
            if ((position == "before" and self.repeats) or
                (position == "assign" and not self.repeats)):
                self._code_set(lines, spacer, position)

    def _code_set(self, lines, spacer, position):
        """Adds the code to set up the variable assignment with this loop at position
        'before' or 'assign'.
        """
        #Before we can assign a value, we need to make sure the value identifier is valid
        #if it exists.
        for valueid in self.value:
            if valueid is not None and valueid.lower() not in self.assignment.values:
                raise ValueError("{} in part {} is not ".format(valueid, self.identifier) + 
                                 "a valid identifier.")
            
        #Before we attempt to assign values, we need to make sure that the variable is allocated
        #For non-file allocations, the Assignment instance handles allocation; for file types,
        #usually the AssignmentValue instance handles it, except when <part> tags are used.
        #However, if the tag has multiple <value> references, we only need to allocate *once* for
        #the variable referenced by the <assignment> tag.
        hasfile = any([self.assignment.values[v].filename is not None for v in self.value])
        if hasfile:
            if self.assignment.allocate:
                if self.assignment.allocate == True:
                    if self.assignment.variable.D == 0:
                        lines.append("{}allocate({})".format(spacer, self.assignment.name))
                    else:
                        raise ValueError('Using allocate="true" for a <part> tag is not valid '
                                         'for multi-dimensional arrays; specify the dimensionality '
                                         'explicitly: e.g. allocate="size(N, 1)" or allocate="10,5".')
                else:
                    lines.append("{}allocate({}({}))".format(spacer, self.assignment.name,
                                                             self.assignment.allocate))            

        if self.isloop:
            lines.append("{}do {}={}, {}".format(spacer, self.loopid, self.limits[0], self.limits[1]))
            i = 0
            for valueid in self.value:
                i += 1
                self.assignment.values[valueid].code(lines, position, spacer + '  ', 
                                                     (self.slices, self.slice_list), first=i==1)
        else:
            #We need to repeat the assignment for each of the part values
            for ipart in self.limits:
                slices = self.slices.format(ipart)
                i = 0
                for valueid in self.value:
                    i += 1
                    self.assignment.values[valueid].code(lines, position, spacer + '  ', 
                                                         (slices, self.slice_list), first=i==1)

        if self.isloop:
            lines.append("{}end do".format(spacer))
        lines.append("")

    def _loopvar_replace(self, limit):
        """Replaces the $i value in the specified limit with the relevant loop variable.
        """
        #The loop variables have to be one this part's parent <part> tags. If it has
        #no parents, or limit is not a string, we don't have to do anything.
        if not isinstance(limit, str) or isinstance(self.parent, Assignment):
            return limit
        else:
            for i in range(1, len(self.top_parts)):
                rstr = "${}".format(i)
                limit = limit.replace(rstr, self.top_parts[i-1].loopid)
            return limit

    def _parse_xml(self):
        if "identifier" in self.xml.attrib:
            self.identifier = self.xml.attrib["identifier"]
        else:
            raise ValueError("The 'identifier' attribute is required for <part> tags.")

        if "start" in self.xml.attrib:
            self.start = int(self.xml.attrib["start"])-1
        if "value" in self.xml.attrib:
            self.value = re.split(",\s*", self.xml.attrib["value"].lower())
        if "limits" in self.xml.attrib:
            if ":" in self.xml.attrib["limits"]:
                self.isloop = True
                self.limits = self.xml.attrib["limits"].split(":")
                if self.limits[0] != "":
                    self.limits[0] = self._loopvar_replace(self.limits[0])
                else:
                    self.limits[0] = 1
                if self.limits[1] != "":
                    self.limits[1] = self._loopvar_replace(self.limits[1])
                else:
                    self.limits[1] = "size({}, {})".format(self.assignment.name, self.depth)
            else:
                self.isloop = False
                self.limits = list(map(int, self.xml.attrib["limits"].split()))
        else:
            self.isloop = True
            self.limits = (1, "size({}, {})".format(self.assignment.name, self.depth))

        if "repeats" in self.xml.attrib:
            self.repeats = self.xml.attrib["repeats"] == "true"

        for kid in self.xml:
            if kid.tag == "part":
                p = Part(kid, self)
                self.parts[p.identifier] = p

class Assignment(object):
    """Represents an instance value assignment as part of a pre-req chain.

    :attr element: the DocElement instance for the <assignment> tag that this
      assignment object represents.
    :attr parent: the MethodFinder instance who owns this assignment.
    :attr methods: a list that will never get used; needed for compatibility
      with method.MethodWriter._order_dependencies recursion.
    :attr group: the parent TestingGroup instance, even if the immediate parent
      is not a testing group.
    :attr name: the name of the variable whose value will be set.
    :attr conditionals: a list of 'if' blocks that conditionally set the value
      of the variable at run-time.
    :attr values: a dict of AssignmentValues capable of setting the variable
      value from a constant, function, file or by calling an embedded method 
      if the variable is an instance of a derived type.
    :attr allocate: when true, the variable having its value set will be
      allocated by this assignment object rather than some other method called
      previously. Default is True.
    """
    def __init__(self, element, parent):
        self.element = element
        self.parent = parent
        """An instance of TestSpecification or TestingGroup in which this <assignment>
        was defined."""
        self.methods = []
        if isinstance(self.parent, TestingGroup):
            self.group = self.parent
        else:
            self.group = self.parent.testgroup

        self.name = None
        self.conditionals = []
        self.values = {}
        self.parts = {}
        self.allocate = False
        self.value = None
        self.writekey = None
        """Establishes a unique key for the driver writer so that multiple assignments
        of the same variable don't interfere with each other."""
        self.position = "before"
        """If a test specification includes assignments local to the test then *global*
        <assignments> attached to the testing group are added either before/after the
        local assignments. This value specifies that position relative to the local ones."""

        self._variable = None

        #Instances can have their values changed by:
        # - direct assignment to a constant value or function call.
        # - direct assignment as part of a set of logical tests.
        # - from the contents of a file.
        # - by calling an embedded method within a derived type.
        self._parse_xml()

    @property
    def autoclass(self):
        """Returns True if any of the assignment values in an autoclass.
        """
        return any(v.autoclass for v in self.values.values())
        
    @property
    def testspec(self):
        """The currently active test specification for which this Assignment instance
        is being coded.
        """
        if self.group.finder is not None:
            return self.group.finder.test
        
    @property
    def parser(self):
        """Returns this Assignment's parent's CodeParser instance."""
        return self.group.element.module.parent
    
    def global_attr(self, key, default=None):
        """Retuns the value of attribute with the specified key from the GlobalDeclaration
        instance for this variable being assigned.
        """
        g = self.globaldecl
        if g is not None and key in g.attributes:
            return g.attributes[key]
        else:
            return default

    @property
    def globaldecl(self):
        """Returns the GlobalDeclaration instance for the current variable, which may have
        different parameters and modifiers than the actual code element.
        """
        #We need to give preference to the test specification if it is available.
        vars = None
        if self.group.finder is not None and self.group.finder.test is not None:
            vars = self.group.finder.test.variables
        if vars is None:
            vars = self.group.variables

        if self.name.lower() in vars:
            return vars[self.name.lower()]
        else:
            return None

    @property
    def variable(self):
        """Returns the ValueElement instance that this assignment refers to."""
        if self._variable is None:
            #Search the global declarations for the variable and then track down
            #the code element that it represents
            if self.name.lower() in self.group.element.parameters:
                self._variable = self.group.element.parameters[self.name.lower()]
            if self._variable is None and self.globaldecl is not None:
                module = self.group.element.module
                self._variable = self.globaldecl.value_elem(module)

        return self._variable
    
    @property
    def attributes(self):
        """Provides one-level-up access to the XML elements attributes collection."""
        if isinstance(self.element, DocElement):
            return self.element.attributes
        else:
            return self.element.attrib

    @property
    def allocatable(self):
        """Returns true if the variable assigned by this object should be allocated."""
        #We use the global declaration's attributes when deciding how to treat the 
        #variable in the unit testing application. This is because the developer can
        #override some of the variable's behavior with <global> tags and we need to
        #honor those changes.
        return (self.allocate and 
                ("allocatable" in self.global_attr("modifiers", "") or
                 "pointer" in self.global_attr("modifiers", "") or
                 (self.variable.D > 0 and ":" in self.variable.dimension) or
                 #If the dimension is a complicated set of functions or variable names, then
                 #we can't initialize it except with ':' and allocatable. That gets added by
                 #the GlobalDeclaration.definition(), so add this check here.
                 (re.match("[\w]+", self.variable.dimension) and self.variable.D > 0)))
    
    def check_prereqs(self, finder):
        """Checks all the values this assignment may use to make sure they reference
        a valid derived type. If an AssignmentValue specifies to check pre-reqs, 
        they are added to the pre-req chain.

        :arg finder: the instance of MethodFinder for which the pre-reqs are being
          checked and added.
        """
        for v in self.values:
            self.values[v].check_prereqs(finder)

    def copy(self, coderoot, testroot, case, compiler, stagedir=None):
        """Copies the input files needed for this assignment to set a variable.

        :arg coderoot: the full path to folder with the code files.
        :arg testroot: the full path the folder where the test is running.
        :arg case: the case id for multi-case testing.
        :arg compiler: the name of the compiler used to make the unit test executable.
        :arg stagedir: the path to the staging directory to use (override).
        """
        #Just call copy on all the child value objects.
        for v in self.values:
            self.values[v].copy(coderoot, testroot, case, compiler, stagedir)        
        
    def _code_setvar_allocate(self, lines, spacer):
        """Allocates the variables before a general value setting if they need to be
        allocated.
        """
        #This only works if the value they specified includes a specific allocate dimension
        #or we can easily determine what it needs to be.
        if self.allocate is None:
            return

        variable = self.variable if self.globaldecl is None else self.globaldecl
        if (variable.dimension is not None and
            self.allocatable and variable.D >= 1):
            lines.append("{}allocate({}({}))".format(spacer, self.name, self.allocate))

        if ("pointer" in self.global_attr("modifiers", "") and
            ("class" in self.global_attr("type", "") or "type" in self.global_attr("type", ""))):
            if self.allocate == True or not self.allocate:
                lines.append("{}allocate({})".format(spacer, self.name))
            else:
                lines.append("{}allocate({}({}))".format(spacer, self.name, self.allocate))

    def code(self, lines, position, spacer):
        """Appends the code for this assignment to function correctly at the
        specified position in the code file.

        :arg lines: the list of strings that form the code file.
        :arg position: one of ['vars', 'init', 'assign', 'before', 'after'] corresponding
          to variable declaration, initialization, assignment and cleanup.
        """
        if position in ["init", "vars", "after"]:
            #Recognizing that even though the main <assignment> may reference only a single
            #<value> or <part> tag, because of nested parts we probably still need to do
            #the initializations.
            i = 0
            for v in self.values:
                self.values[v].code(lines, position, spacer, False, i==0)
                i+=1
            for p in self.parts:
                self.parts[p].code(lines, position, spacer)
        elif position in ["assign", "before"]:
            #Usually, before any kind of assignment can take place, the variable may need
            #to be allocated. Here we try allocation for the variable if it applies. For
            #file variables, the allocation requires special values extracted at runtime
            #from the input files.
            for v in self.values:
                value = self.values[v]
                if (value.filename is None and (
                    (value.repeats and position=="before") or 
                    (not value.repeats and position=="assign"))):
                    self._code_setvar_allocate(lines, spacer)
                    #We only need to allocate the variable once per assignment.
                    break

            #We also need to handle the case where the developer only wants the variable
            #to be allocated and have nothing else done to it...
            if (self.value is None and self.allocate and position=="assign"):
                self._code_setvar_allocate(lines, spacer)

            #We just need to check whether this assignment had an explicit
            #value or if it has a list of conditionals.
            if self.value is None and len(self.conditionals) > 0:
                for c in self.conditionals:
                    c.code(lines, position, spacer)
            elif self.value is not None:
                #We allow multiple <values> to be assigned in the context of the current
                #coding. They go in the order they showed up in the list.
                i = 0
                for v in self.value:
                    i += 1
                    if v in self.values:
                        self.values[v].code(lines, position, spacer, first=i==1)
                    elif v in self.parts:
                        self.parts[v].code(lines, position, spacer)
                    else:
                        raise ValueError("Value identifier {} not found.".format(v))
        
    def _parse_xml(self):
        """Parses attributes and child tags from the XML element."""
        if "name" in self.attributes:
            if self.attributes["name"] == "[default]":
                self.name = self.group.element.name + "_fpy"
            else:
                self.name = self.attributes["name"]
        else:
            raise ValueError("'name' is a required attribute for <assignment> tags.")

        if "value" in self.attributes:
            self.value = re.split(",\s*", self.attributes["value"].lower())

        def autoval(attr, value):
            """Creates an automatic AssignmentValue instance for the specified
            attribute and value.
            """
            val = AssignmentValue(None, self)
            val.identifier = attr
            setattr(val, attr, value)
            self.values[attr] = val
            if self.value is None:
                self.value = [attr]
            else:
                self.value.append(attr)
            
        if "constant" in self.attributes:
            #We want to manually create an AssignmentValue object for the constant and
            #add it to the child list.
            autoval("constant", self.attributes["constant"])
        if "testsource" in self.attributes:
            tsource = TestSource(self.attributes["testsource"].lower(), self.testspec)
            autoval("testsource", tsource)
            self.values["testsource"].autoclass = tsource.target.autoclass
            
        if "allocate" in self.attributes:
            if self.attributes["allocate"] in ["false", "true"]:
                self.allocate = self.attributes["allocate"] == "true"
            else:
                self.allocate = self.attributes["allocate"]
        if "position" in self.attributes:
            self.position = self.attributes["position"]        

        kids = self.element.xml if isinstance(self.element, DocElement) else self.element
        for child in kids:
            if child.tag == "value":
                val = AssignmentValue(child, self)
                #Having automation of value identifiers disabled the check to make
                #sure that every <value> tag has a unique one. Perform that check now.
                if "default" in self.values and val.identifier == "default":
                    raise ValueError("'identifier' is a required attribute of <value> tags "
                                     "unless there is only one of them in the <assignment>.")
                self.values[val.identifier.lower()] = val
            elif child.tag == "part":
                p = Part(child, self)
                self.parts[p.identifier.lower()] = p
            elif child.tag == "conditionals":
                self.conditionals.append(AssignmentConditional(child, self))
        
        #Construct the key for identifying this assignment operation uniquely.
        self.writekey = self.name
        if len(self.parts) > 0:
            self.writekey += "({})".format(list(self.parts.keys())[0])
        if self.value is not None:
            self.writekey += "_{}".format('|'.join(self.value))

        #If we only have a single value without any fancy conditionals etc.
        #then the value identifier and attribute are unnecessary.
        if (self.value is None and len(self.values) == 1 and
            len(self.parts) == 0 and len(self.conditionals) == 0 and
            "default" in self.values):
            self.value = ["default"]
         
class TestInput(object):
    """Represents information for a single input source as part of a unit test.

    :attr xml: the <input> tag that this object was created from.
    :attr folder: the path to the folder that contains the input file relative
      to the code root.
    :attr filename: the name of the file to use as input.
    :attr rename: the new name to use for the input file for actual execution.
    :attr line: a FileLine instance describing the values that appear in a single
      line of the input file. Use for constant-input mode.
    """
    def __init__(self, xmltag, testspec):
        """Constructor.

        :arg testspec: the TestSpecification parent instance that has access
          to code element information for test-source inputs.
        """
        self.xml = xmltag
        self.testspec = testspec
        """TestSpecification parent instance that has access
          to code element information for test-source inputs."""
        self.folder = None
        self.filename = None
        self.rename = None
        self.line = None
        self.position = "before"
        self.testsource = None
        """Specifies the module.executable.parameter name whose output from a previous
        unit test run will be used as input for this value.
        """

        self._parse_xml()
        
    @property
    def constant(self):
        """Returns true if operating in constant-input mode."""
        return (self.line is not None)

    @property
    def iid(self):
        """Returns the identifier for this input file when running in constant-
        input mode."""
        if self.line is not None:
            return self.line.identifier
        else:
            return ""

    @property
    def fileunit(self):
        """Returns the name used for the variable that holds the file unit
        for reading the file in constant-input mode."""
        #Come up with a name for the file unit that will be unique; since the
        #line template has to have an identifier, use that.
        return "fpy_file_{}".format(self.iid)

    @property
    def repeats(self):
        """Returns the variable that holds the number of times to repeat 
        the read operation and execute the method.
        """
        return "fpy_repeat_{}".format(self.iid)

    @property
    def iostat_read(self):
        """Returns the name of the variable that holds read IOSTAT for the
        read operations when running in constant-input mode."""
        return "fpy_iostat_{}".format(self.iid)

    @property
    def runfile(self):
        """Returns the name of the file that will be in the testing folder
        at run-time.
        """
        #We don't need to take cases into account because the file must be
        #renamed when running in case mode.
        if self.rename is not None:
            return self.rename
        else:
            return self.filename

    def code_vars(self, lines, spacer):
        """Appends the fortran code for declaring any variables need to run 
        the constant-input mode.
        """
        if self.constant:
            lines.append("{}!Constant-input mode variables for {}".format(spacer, self.iid))
            lines.append("{}integer :: {}, {}, {} = 0".format(spacer, self.fileunit,
                                                            self.iostat_read,
                                                            self.repeats))

    def code_init(self, lines, spacer):
        """Appends the initialization code required for opening the input file
        when running in constant-input mode."""
        if self.constant:
            dlines = []
            #We need to open the file with a new file unit, which we can do using
            #the fpy_newunit() function in the fortpy module.
            dlines.append("!Auto-generated for {}; determine number of constants.".format(self.iid))
            dlines.append("open(fpy_newunit({}), file='{}')".format(self.fileunit,
                                                                   self.runfile))
            dlines.append("do")
            dlines.append("  read({}, *, IOSTAT={})".format(self.fileunit, self.iostat_read))
            dlines.append("  if ({} == 0) then".format(self.iostat_read))
            dlines.append("    {0} = {0} + 1".format(self.repeats))
            dlines.append("  else")
            dlines.append("    exit")
            dlines.append("  end if")
            dlines.append("end do")
            dlines.append("rewind({})".format(self.fileunit))
            dlines.append("")

            lines.extend([spacer + l for l in dlines])

    def code_read(self, lines, spacer):
        """Appends the fortran code for reading the values from the file into
        the global variables for a single iteration.
        """
        if self.constant:
            #Use the line parser information to put together a read format string
            varlist = ", ".join(self.line.unique_names)
            lines.append("{}read({}, *) {}".format(spacer, self.fileunit, varlist))
        
    def code_finalize(self, lines, spacer):
        """Appends the code to finalize the constant-input operation."""
        if self.constant:
            lines.append("{}close({})".format(spacer, self.fileunit))

    def copy(self, coderoot, testroot, case="", compiler=None, stagedir=None):
        """Copies the input file from the specified code root directory to
        the folder where the test is being performed. Alternatively, the file is
        just symlinked if that is the global setting.

        :arg coderoot: the full path to the folder that houses all the code files.
        :arg testroot: the full path to the folder that the parent test is being
          performed in.
        :arg case: if a specific case of the same test is being performed, the
          case identifier to use for string formatting.
        :arg compiler: the name of the compiler used to make the unit test executable.
        :arg stagedir: the path to the staging directory to use (override).
        """
        if self.rename is not None:
            target = path.join(testroot, self.rename)
        else:
            target = path.join(testroot, self.filename.format(case))

        if self.testsource is None:
            from fortpy.tramp import coderelpath
            from fortpy.testing.compilers import replace
            relpath = coderelpath(coderoot, self.folder)
            source = replace(path.join(relpath, self.filename.format(case)), compiler)
            source = replace(source, compiler, True)

            from fortpy.utility import symlink
            symlink(source, target)
        else:
            self.testsource.copy(coderoot, testroot, case, compiler, stagedir)

    def _parse_xml(self):
        """Extracts the input file attributes from the xmltag."""
        if "folder" in self.xml.attrib:
            self.folder = self.xml.attrib["folder"]
        if "file" in self.xml.attrib:
            self.filename = self.xml.attrib["file"]
        if "rename" in self.xml.attrib:
            self.rename = self.xml.attrib["rename"]
        if "position" in self.xml.attrib:
            self.position = self.xml.attrib["position"]
        if "testsource" in self.xml.attrib:
            self.testsource = TestSource(self.xml.attrib["testsource"].lower(), self.testspec)
            
        self._parse_lines()

    def _parse_lines(self):
        """Checks for the existence of a <line> tag in the input specification.
        If it exists, the line attribute is initialized with a FileLine object.
        """
        for child in self.xml:
            if child.tag == "line":
                self.line = FileLine(child, None)
                break

class TestOutput(object):
    """Represents an output file *compare* operation that needs to be performed
    by the tester.

    :attr xml: the XML <output> tag that the object was created from.
    :attr identifier: the unique identifier for this output specification.
    :attr folder: for file-mode, the folder that contains 'file'.
    :attr filename: the name of the model file to compare the output to.
    :attr template: the name of the XML file that has template information
      for comparing the file.
    :attr mode: the comparison mode for the FileComparer. Default='default'.
    :attr value: if the output value is a constant, its value.
    :attr tolerance: for value comparisons, as long as the values are within
      this percent, the test still passes. E.g. if tolerance=0.8, then a value
      of 0.92 is considered equivalent to 1.
    """
    def __init__(self, xmltag):
        """Initializes the output comparison using an <output> tag."""
        self.xml = xmltag
        self.identifier = None
        self.folder = None
        self.filename = None
        self.template = None
        self.mode = "default"
        self.value = None
        self.tolerance = None
        self.position = "before"
        self.autoclass = False
        """Specifies whether the output is running in auto-class mode where entire
        arrays of user-defined objects can be compared using a single tag.
        """
        self.actolerance = None
        """Specifies the individual tolerance limit for any single file comparison
        in the auto-class comparision.
        """

        self._parse_attributes()

    @property
    def filemode(self):
        """Returns a value indicating whether the output comparison is with
        a file (as opposed to a hard-coded value)."""
        return (self.value == None)

    def abspath(self, coderoot, caseid="", compiler=None):
        """Returns the absolute path to the output file referenced by this
        output comparison specifier.

        :arg coderoot: the full path to the directory storing the code files.
        :arg case: if a specific case of the same test is being performed, the
          case identifier to use for string formatting.
        """
        #The caseid is the "testid.case". We only want the 'case' part of it
        from fortpy.tramp import coderelpath
        from fortpy.testing.compilers import replace
        relpath = coderelpath(coderoot, self.folder)
        if caseid != "":
            case = caseid.split(".")[-1]
            if self.autoclass:
                formpath = relpath.format(case)
            else:
                formpath = path.join(relpath, self.filename.format(case))
        else:
            if self.autoclass:
                formpath = relpath
            else:
                formpath = path.join(relpath, self.filename)

        #Finally, replace the compiler attributes [c] and [f] for the file path.
        if compiler is not None:
            cpath = replace(formpath, compiler)
            return replace(formpath, compiler, True)
        else:
            return formpath

    def _parse_attributes(self):
        """Extracts output comparsion related attributes from the xml tag."""
        if "identifier" in self.xml.attrib:
            self.identifier = self.xml.attrib["identifier"]
        else:
            raise ValueError("'identifier' is a required attribute for an <output> tag.")
        if "folder" in self.xml.attrib:
            self.folder = self.xml.attrib["folder"]
        if "file" in self.xml.attrib:
            self.filename = self.xml.attrib["file"]
        if "template" in self.xml.attrib:
            self.template = self.xml.attrib["template"]
        if "mode" in self.xml.attrib:
            self.mode = self.xml.attrib["mode"]
        if "value" in self.xml.attrib:
            self.value = self.xml.attrib["value"]
        if "tolerance" in self.xml.attrib:
            self.tolerance = float(self.xml.attrib["tolerance"])
        else:
            self.tolerance = 1.
        if "actolerance" in self.xml.attrib:
            self.actolerance = float(self.xml.attrib["actolerance"])
        else:
            self.actolerance = 1.

        if "position" in self.xml.attrib:
            self.position = self.xml.attrib["position"]
        if "autoclass" in self.xml.attrib:
            self.autoclass = self.xml.attrib["autoclass"].lower()=="true"    
        
class TestTarget(object):
    """Represents a specification to save the state of a variable at the end
    of the program run *or* each time the main method being tested is run.

    :attr xml: the xml <target> tag that the object was created from.
    :attr name: the name of the variable that will be saved.
    :attr varfile: if fortpy is handling the save to file, the name of the file
      that fortpy should save the variable under.
    :attr when: specifies how often the variable should be saved. Possible values
      are ['begin', 'each', 'end'] corresponding to the start of the program,
      each time the method gets run and right before the variables are cleaned up
      by the program.
    :attr generator: the name of a subroutine to call that will save the variable
      to a file. The only parameter passed to the subroutine is the instance of the
      variable to save.
    :attr compareto: the identifier of an <output> tag that the variable's values will
      be compared to after the program has run.
    :attr testspec: the parent test specification that this target belongs to.
    :attr member: the name of a member in a derived type whose value should be saved.
    """
    def __init__(self, xmltag, testspec):
        """Initializes the target specification with a <target> tag."""
        self.xml = xmltag
        self.name = None
        self.varfile = None
        self.when = None
        self.generator = None
        self.compareto = None
        self.testspec = testspec
        self.member = None
        self.position = "before"
        self.autoclass = False
        """Specifies whether the output is running in auto-class mode where entire
        arrays of user-defined objects can be compared using a single tag.
        """

        self._code = None
        """The code to *save* the variable referenced by this target to file."""
        self._parse_attributes()

    def code(self, lines, variables, position, spacer, first=True):
        """Appends the code required to save the value of this target to
        file if the specified position is appropriate.
        """
        #If the target is trying to get a file to be compared to another
        #model output file, then it must start with a dot. In that case
        #we don't generate any code since we assume that the file is being
        #created by the method being unit tested.
        if self.name[0] == '.':
            return
        
        if self._code is None:
            self._process(variables, spacer)

        if position in self.when:
            lines.append(self._code)
        elif position == "vars" and self.autoclass:
            varkey = self._check_autoclass(self.testspec.executable)
            self.testspec.autoclasses[varkey].xcode(lines, "vars", spacer, first=first)

    def init(self, testroot):
        """Creates any directories that this test target needs to function correctly.
        """
        if self.autoclass:
            from os import mkdir
            vardir = path.join(testroot, self.varfile)
            if not path.isdir(vardir):
                mkdir(vardir)

    def clean(self, testroot):
        """Removes the output file for saving the variable's value from the
        testing folder."""
        if self.varfile is not None:
            target = path.join(testroot, self.varfile)
            if self.autoclass:
                from shutil import rmtree
                if path.isdir(target):
                    rmtree(target)
            else:
                target = path.join(testroot, self.varfile)
                if path.isfile(target):
                    remove(target)

    def _process(self, variables, spacer):
        """Processes the code to generate the output of the unit test.
        """
        if self.autoclass:
            self._process_autoclass(variables, spacer)
        else:
            self._process_simple(variables, spacer)

    def _check_autoclass(self, variable):
        """Makes sure that the autoclass for implementing this variable assignment
        exists in the parent test specification.
        """
        varkey = self.name.lower()
        acdict = self.testspec.autoclasses
        if varkey not in acdict:
            acdict[varkey] = AutoClasser(variable, ".", self.varfile, self.testspec.cases,
                                         self.testspec.testgroup.coderoot)
        else:
            acdict[varkey].folder = "."
            acdict[varkey].varfile = self.varfile

        return varkey
            
    def _process_autoclass(self, variables, spacer):
        """Processes the test target using an AutoClasser().
        """
        #See the comments in _process_simple about these logic trees.
        if self.name == "[default]":
            varkey = self._check_autoclass(self.testspec.executable)
        elif self.name in variables:
            glob = variables[self.name]
            if isinstance(glob, GlobalDeclaration):
                finder = self.testspec.testgroup.finder
                variable = glob.value_elem(finder.module)
            elif type(glob).__name__ == "ValueElement":
                variable = glob
            varkey = self._check_autoclass(variable)
        else:
            raise ValueError("Output target auto-class '{}' has not".format(self.name) +
                             " been initialized with a <global> tag or regular=true parameter doc tag.")

        icode = []
        self.testspec.autoclasses[varkey].xcode(icode, "save", spacer)
        self._code = '\n'.join(icode)
            
    def _process_simple(self, variables, spacer):
        """Processes a single outcome involving a variable value comparison.

        :arg variables: a list of the GlobalDeclarations made for the unit test.
        """
        #The framework allows the developer to specify a subroutine to create the
        #output file for comparison. If one was specified, just use that.
        if self.generator is not None:
            if self.varfile is not None:
                self._code = "{}call {}({}, '{}')".format(spacer, self.generator, 
                                                          self.name, self.varfile)
            else:
                self._code = "{}call {}({})".format(spacer, self.generator, self.name)
        else:
            #We need to use the fortpy module or the variables own test_output()
            #to create a file that can be compared later by python. This means
            #that we need to know the type and kind of the variable whose value
            #needs to be compared
            if self.name == "[default]":
                #We need to save the value that was generated by the main *function* being
                #unit tested. The variable defined in the fortpy program is the function
                #name with a suffix of "_fpy".
                self._code = "{}call pysave({}_fpy, '{}')".format(spacer, self.testspec.executable.name,
                                                                  self.varfile)
            elif self.name in variables:
                glob = variables[self.name]
                dtype = glob.attributes["type"]
                if dtype == "class" or dtype == "type":
                    if self.member is not None:
                        self._code = ("{}call pysave({}%{}".format(spacer, self.name, self.member) + 
                                      ", '{}')".format(self.varfile))
                    else:
                        self._code = "{}call {}%test_output('{}')".format(spacer, self.name, self.varfile)
                else:
                    #pysave is an interface in the fortpy module that can accept variables
                    #of different types and write them to file in a deterministic way
                    self._code = "{}call pysave({}, '{}')".format(spacer, self.name, self.varfile)
            else:
                raise ValueError("Output target var {} has not been initialized with ".format(self.name) +
                                 "a <global> tag or regular=true parameter doc tag.")

    def _parse_attributes(self):
        """Extracts implemented attributes from the test target."""
        if "name" in self.xml.attrib:
            self.name = self.xml.attrib["name"].lower()
        else:
            raise ValueError("'name' is a required attribute for a <target> tag.")
        if "varfile" in self.xml.attrib:
            self.varfile = self.xml.attrib["varfile"]
        else:
            self.varfile = "{}".format(self.xml.attrib["name"].replace("%", "."))
        if "when" in self.xml.attrib:
            self.when = re.split(",\s*", self.xml.attrib["when"])
        else:
            self.when = ["end"]
        if "generator" in self.xml.attrib:
            self.generator = self.xml.attrib["generator"]
        if "compareto" in self.xml.attrib:
            self.compareto = self.xml.attrib["compareto"]
        if "member" in self.xml.attrib:
            self.member = self.xml.attrib["member"]
        if "position" in self.xml.attrib:
            self.position = self.xml.attrib["position"]
        if "autoclass" in self.xml.attrib:
            self.autoclass = self.xml.attrib["autoclass"].lower()=="true"

class TestSpecification(object):
    """Represents a single test that needs to be performed. Each test compiles
    a new executable with the same name as the test identifier; however, the
    tests all use the same base dependencies from other modules.

    :attr identifier: a unique id for the test; the unit test executable will
      be created with this name.
    :attr description: included in the summary of the auto-generated program
      *.f90 file.
    :attr cases: a list of cases to run that use the same test specification
      except for a slight change in the names of input/model files.
    :attr runtime: an integer list [min, max] range for the time this test 
      should take to run.
    :attr unit: either 'm' or 'h'; the units for 'runtime'.
    :attr timed: when true, the execution time of the main method being unit tested
      will be accumulated and saved.
    :attr targets: a list of TestTarget instances that specify which variables
      and files need to be verified after the unit test runs.
    :attr inputs: a list of TestInput instances that specify how to initialize
      the values of the variables or that are required by other methods being
      run as pre-reqs of the main unit testing method.
    :attr output: a dict of TestOutput instances; keys are the output identifiers
      the values are the object instances. Represent sources for model data to
      compare the unit tests' generated output to.
    :attr testgroup: the testing group that the specification belongs to.
    """
    def __init__(self, xmltag, testgroup):
        """Initializes the test specification using the <test> tag that is
        a child of the <outcome> tag.
        """
        self.xml = xmltag
        self.identifier = None
        self.description = None
        self.cases = None
        self.runtime = None
        self.unit = None
        self.timed = False
        self.execute = True
        self.runchecks = True

        self.targets = []
        self.inputs = []
        self.outputs = {}
        self.methods = []
        """A set of *local* assignments and pre-reqs that apply to this test."""
        self.variables = {}
        """A set of *local* variable declarations for this test."""
        self._variable_order = []
        self.autoclasses = {}
        """A dict of AutoClasser() instances for generating read/write blocks
        for user-derived type variables that are complex. We only need to write
        the variable declarations once per test specification, which is why we
        put them here.
        """

        self.testgroup = testgroup

        #The name of the variable that holds the number of times to repeat
        #the execution of the method being tested.
        self._repeats = None
        #If any of the inputs runs in constant mode, this variable with be True.
        self._constant = None

        #The identifier is necessary to proceed; all the other xml attributes
        #and tags are parsed once the parent TestingGroup instance has parsed
        #*all* of its children.
        if "identifier" in self.xml.attrib:
            self.identifier = self.xml.attrib["identifier"]
        else:
            #Let's auto-assign an identifier.
            root = "std"
            for child in testgroup.group.xml.iterfind("test"):
                if "identifier" in child.attrib and child.attrib["identifier"] == root:
                    root += "I"
            self.identifier = root

    @property 
    def executable(self):
        """Returns the instance of Executable code element that this specification
        will test."""
        return self.testgroup.element

    @property
    def testable(self):
        """Returns True if this test specification includes at least one output
        specification making it testable."""
        return len(self.outputs) > 0

    @property
    def autoclass(self):
        """Returns True if any of the assignments in the test specification require
        auto-class support.
        """
        reads = any([a.autoclass for a in self.methods if isinstance(a, Assignment)])
        writes = any([o.autoclass for o in self.outputs.values()])
        return reads or writes
    
    @property
    def constant(self):
        """Returns a value indicating whether any of the inputs in this test
        run in constant-input mode."""
        if self._constant is None:
            for i in self.inputs:
                if i.constant:
                    self._constant = True
                    break
            else:
                self._constant = False

        return self._constant

    @property
    def repeats(self):
        """Returns the name of a variable to use for the do loop when the
        method execution is repeated for each line in constant-input file."""
        #Find the first input specification that has operates in constant mode
        if self._repeats is None and self.constant:
            for i in self.inputs:
                if i.constant:
                    self._repeats = i.repeats
                    break

        return self._repeats

    def clean(self, testroot):
        """Removes all variable target output files so that the test folder
        is clean for a new run."""
        for t in self.targets:
            t.clean(testroot)

        #We also need to remove the timing of any previous tests that ran.
        if self.timed:
            timepath = path.join(testroot, "fpytiming.out")
            if path.isfile(timepath):
                remove(timepath)

    def code_validate(self, lines, spacer):
        """Appends the fortran code to validate the length of the constant-valued
        input files if there is more than one of them.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run.
        """
        #Basically, we just want to generate a bunch of if statements and check
        #that the number of entries in each constant-valued input file is the same.
        if self.constant:
            first = None
            for i in self.inputs:
                if i.constant:
                    if first is None:
                        first = i
                    else:
                        lines.extend(self._code_validate_single(first, i, spacer))
        
    def _code_validate_single(self, first, second, spacer):
        """Generates an if/stop code block for checking the number of entries in
        constant-valued input files.
        
        :arg first: an instance of TestInput to compare with 'second'.
        :arg second: an instance of TestInput to compare with 'first'.
        """
        result = []
        result.append("{}if ({} .ne. {}) then".format(spacer, first.repeats, second.repeats))
        result.append("{}stop('{} and {} must have same number of".format(spacer, first.iid, second.iid) + 
                      " entries in their input files')")
        result.append("{}end if".format(spacer))

        return result

    def code_vars(self, lines, spacer):
        """Appends all the code required by all the input file specifications
        in order to run the unit test to the list.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run."""
        for i in self.inputs:
            i.code_vars(lines, spacer)
            
        count=0
        for t in self.targets:
            t.code(lines, self.variables, "vars", spacer, count==0)
            count+=1

        #We need access to the test case locally for reading value with auto-class.
        for a in self.methods:
            if isinstance(a, Assignment) and a.autoclass:
                lines.append("{}character(10) :: fpy_case".format(spacer))
                break

    def code_init(self, lines, spacer):
        """Appends all code required to initialize all the constant-mode input
        file specifications and save starting values of targets.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run."""
        for i in self.inputs:
            i.code_init(lines, spacer)

        for t in self.targets:
            t.code(lines, self.variables, "begin", spacer)

        for a in self.methods:
            if isinstance(a, Assignment) and a.autoclass:
                lines.append("{}call fpy_read('fpy_case', '#', fpy_case)".format(spacer))
                break

    def code_before(self, lines, spacer):
        """Appends all code required to read the next values from a file into
        their variables *before* the main test function is called again.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run."""
        #Before the unit test method executes, we want to get the next lot
        #of constant variable values read from the next lines in their
        #respective files.
        for i in self.inputs:
            i.code_read(lines, spacer)

    def code_after(self, lines, spacer):
        """Appends the code required to save the value of any targets for this
        test to their files when their save frequency is set to 'each'.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run."""
        for t in self.targets:
            t.code(lines, self.variables, "each", spacer)

    def code_finalize(self, lines, spacer):
        """Appends the code for cleaning up open file handles and saving the
        final values of the targets to their files.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run.
        """
        for t in self.targets:
            t.code(lines, self.variables, "end", spacer)

        for i in self.inputs:
            i.code_finalize(lines, spacer)

    def parse(self):
        """Parses the child-tags and attributes of this test specification once
        its parent testing group is finished parsing.
        """
        def group_list(collection, position, target, docopy=False):
            """Appends all the elements in the collection that match the position
            specification to the given target.
            """
            if docopy:
                from copy import deepcopy
            for item in collection:
                if item.position == position:
                    citem = item if not docopy else deepcopy(item)
                    target.append(citem)

        def group_dict(collection, position, target, ordering=None, tordering=None, docopy=False):
            """Updates the 'target' dictionary to include all items in the 'collection'
            that match the specified position.
            """
            if docopy:
                from copy import deepcopy
            if ordering is not None:
                keys = ordering
            else:
                keys = collection
                
            for item in keys:
                if collection[item].position == position:
                    citem = item if not docopy else deepcopy(item)
                    if tordering is not None:                        
                        tordering.append(citem)
                    target[item] = collection[citem]
            
        #We need to add all the pre-reqs and assignments from the parent group to
        #this one since they were specified globally for *all* tests.
        group_list(self.testgroup.methods, "before", self.methods)
        group_list(self.testgroup.inputs, "before", self.inputs, docopy=True)
        group_list(self.testgroup.targets, "before", self.targets)
        group_dict(self.testgroup.outputs, "before", self.outputs)
        group_dict(self.testgroup.variables, "before", self.variables,
                   self.testgroup.variable_order, self._variable_order)
                        
        self._parse_attributes()
        self._parse_children()

        #Now handle the methods from the group that should be added *after* the local
        #assignments and pre-reqs are done.
        group_list(self.testgroup.methods, "after", self.methods)
        group_list(self.testgroup.inputs, "after", self.inputs, docopy=True)
        group_list(self.testgroup.targets, "after", self.targets)
        group_dict(self.testgroup.outputs, "after", self.outputs)
        group_dict(self.testgroup.variables, "after", self.variables,
                   self.testgroup.variable_order, self._variable_order)

        #If there is only a single target and a single output, then match them up
        #even if there is no "compareto" and "identifier" on the <target> and <output>
        #tags.
        if len(self.targets) == 1 and len(self.outputs) == 1:
            if self.targets[0].compareto is None:
                self.targets[0].compareto = self.outputs.keys()[0]

        #Because of the testsource feature, the <input> tags have a reference to the
        #TestSpecification. We copied any from the TestingGroup above, but the Group
        #doesn't initialize the testspec (since it isn't one). We do that here.
        for xinput in self.inputs:
            if xinput.testspec is None:
                xinput.testspec = self
                
    def _parse_children(self):
        """Gets all child tags from the <test> and parses them for relevant
        content to set targets, inputs and output attributes."""
        for child in self.xml:
            if child.tag == "input":
                self.inputs.append(TestInput(child, self))
            elif child.tag == "output":
                outvar = TestOutput(child)
                self.outputs[outvar.identifier] = outvar
            elif child.tag == "target":
                self.targets.append(TestTarget(child, self))
            elif child.tag == "prereq" and "method" in child.attrib:
                self.methods.append(TestPreReq(child, self))
            elif child.tag == "assignment":
                self.methods.append(Assignment(child, self))
            elif (child.tag == "global" and "name" in child.attrib):
                TestingGroup.global_add(self.variables, self._variable_order,
                                        child.attrib["name"].lower(), child)
                
    def _parse_attributes(self):
        """Gets test-level attributes from the xml tag if they exist."""
        if "description" in self.xml.attrib:
            self.description = self.xml.attrib["description"]
        if "runchecks" in self.xml.attrib:
            self.runchecks = self.xml.attrib["runchecks"] == "true"
        if "execute" in self.xml.attrib:
            self.execute = self.xml.attrib["execute"] == "true"
            if not self.execute:
                #If we don't run the executable, we obviously can't check outputs!
                self.runchecks = False
        if "cases" in self.xml.attrib:
            self.cases = _expand_cases(self.xml.attrib["cases"])
        if "runtime" in self.xml.attrib:
            runtime = list(self.xml.attrib["runtime"])
            self.unit = runtime.pop()
            self.runtime = [ int(t) for t in "".join(runtime).split("-") ]
        if "timed" in self.xml.attrib:
            self.timed = self.xml.attrib["timed"] == "true"
        else:
            self.timed = True

class TestPreReq(object):
    """Represents a subroutine or function that must run before the main
    executable being unit-tested can run.
    """
    def __init__(self, tag, parent):
        self.xml = tag
        self.parent = parent
        """The testing group or test specification instance that owns this tag."""
        self.method = None
        """The 'module.executable' identifier for the pre-req method to run."""
        self.position = "before"
        """For pre-reqs specified in a testing group, the user must decide whether
        the global methods run before/after the local test specification methods.
        """
        self.paramlist = None
        """Overrides the default behavior of the parameter list."""
        self.chainto = None
        """Specifies the ID of the test specification that should be used when
        chaining the pre-reqs.
        """

        self._parse_xml()

    def _parse_xml(self):
        if "method" in self.xml.attrib:
            self.method = self.xml.attrib["method"]
        else:
            raise ValueError("'method' is a required attribute of <prereq> tags.")
        if "chainto" in self.xml.attrib:
            self.chainto = self.xml.attrib["chainto"]
        
        if "position" in self.xml.attrib:
            self.position = self.xml.attrib["position"]
        if "paramlist" in self.xml.attrib:
            self.paramlist = self.xml.attrib["paramlist"]    
            
class TestingGroup(object):
    """Represents a unit testing documentation group with information on
    how to automate a unit test for a specific executable.

    :attr tests: a dict of TestSpecifications for running multiple tests
      that use the same variables, mappings and other parameters.
    :attr mappings: a dict with keys as the parameter names from the method
      being tested and values being the variable names to pass in.
    :attr assignments: a dict of assignments for changing the values of the
      variables.
    :attr children: a list of the DocElement instances who belong to this
      testing group.
    :attr methods: a dict of variable assignments and prereq methods for changing 
      variable values or executing pre-req subroutines. 
      The actual assignments are made as part of the executable chaining
      to make sure that values get changed in the correct order. 
    :attr finder: the MethodFinder instance that this group belongs to.
    """
    def __init__(self, group, element):
        """Parses the specified DocGroup to extract the unit testing info.
        
        :arg group: an instance of a DocGroup with purpose="testing".
        :arg element: the *code* element that the doc group belongs to.
        """
        self.group = group
        self.element = element
        self.tests = {}
        self.mappings = {}
        self.variables = {}
        self.children = []
        self.methods = []
        self.assignments = {}
        """A dict of test-group level variable names that have their values
        changed by a test-group level assignment. Keys are lowered variable names
        values are dicts of variable attributes.
        """
        self.inputs = []
        """A list of input files that need to be copied to the test execution directories
        for *all* the test specifications.
        """
        self.outputs = {}
        """Dict of comparison outputs that should be available for *all* test specifications
        in the unit test.
        """
        self.targets = []
        """List of the variables that should have their values saved to file and compared
        with model output for *all* the test specifications.
        """
        self.staging = None
        """The path of the staging directory relative to the code directory."""
        
        self.finder = None
        """The currently active MethodFinder instance being used by the MethodWriter
        to create the executable driver for a specific testid."""
        self.variable_order = []
        """The order in which the global declarations appear in the group."""
        self.coderoot = None
        """The absolute path to the code folder where the module file resides
        for the unit test being run."""
        self._codes = {}
        self._find_children()
        self._parse_xml()
        self._init_codes()

    @property
    def method_fullname(self):
        """Returns the full module.method name of the executable that this
        group decorates.
        """
        return "{}.{}".format(self.element.module.name, self.element.name)

    @property
    def name(self):
        """Returns the name of the underlying DocGroup."""
        if self.group.name is not None:
            return self.group.name
        else:
            return "unittests"
    
    def set_finder(self, finder):
        """Sets the currently active MethodFinder instance being used by the MethodWriter
        to create the executable driver for a specific testid."""
        self.finder = finder
    
    def _init_codes(self):
        """Adds a dictionary of function pointers for appending relevant
        code to the final executable for each test in the outcome."""
        for key in self.tests:
            mapper = {
                "vars": self.tests[key].code_vars,
                "init": self.tests[key].code_init,
                "validate": self.tests[key].code_validate,
                "before": self.tests[key].code_before,
                "after": self.tests[key].code_after,
                "final": self.tests[key].code_finalize
            }
            self._codes[key] = mapper

    def code(self, testid, section, lines, spacer):
        """Appends the relevant code for the specified test and section of
        the final executable file.

        :arg testid: the identifier of the test to append code for.
        :arg section: one of ['vars', 'init', 'validate', 'before', 'after',
          'final']. Refers to section of the PROGRAM *.f90 file being created.
        :arg lines: a list of the strings collected so far for the contents
          of the *.f90 program file.
        """
        if testid in self._codes and section in self._codes[testid]:
            self._codes[testid][section](lines, spacer)

    def _find_children(self):
        """Finds all the DocElement instances in the parent code element 
        that belong to this testing group."""
        for docel in self.element.docstring:
            if docel.group is not None:
                if ((isinstance(docel.group, str) and docel.group == self.group.name) or
                    (isinstance(docel.group, DocGroup) and docel.group.name == self.group.name)):
                    self.children.append(docel)

    def _parse_xml(self):
        """Parses the XML structure into the relevant dictionaries and
        properties for the testing group.
        """
        #Try to determine the staging directory from the group definition.
        if self.staging is None and "staging" in self.group.xml.attrib:
            self.staging = self.group.xml.attrib["staging"]

        for child in self.children:
            if child.doctype == "test":
                test = TestSpecification(child.xml, self)
                self.tests[test.identifier] = test
            elif child.doctype == "prereq" and "method" in child.attributes:
                self.methods.append(TestPreReq(child.xml, self))
            elif child.doctype == "assignment":
                #This is a variable assignment. Sometimes variable values are changed
                #between calls to functions. The order in which prereq/instance elements
                #appear defines when these assignments take place. We will treat a var
                #value assignment as a method to simplify the implementation.
                self.methods.append(Assignment(child, self))
                #If the variable gets its value changed multiple times, the first time
                #should have allocation if applicable, otherwise it wouldn't work.
                self.assignments[child.attributes["name"].lower()] = child.attributes
            elif child.doctype == "input":
                self.inputs.append(TestInput(child.xml, None))
            elif child.doctype == "output":
                outvar = TestOutput(child.xml)
                self.outputs[outvar.identifier] = outvar
            elif child.doctype == "target":
                self.targets.append(TestTarget(child.xml, self))
            elif child.doctype == "global":
                self._parse_childvar(child)

        self._parse_mappings()
        if type(self.element).__name__ in ["Subroutine", "Function"]:
            self._get_param_globals()

        #Now that the testing group has all its children parsed, we can
        #let the tests parse. This is because some of the global tags in the
        #testing group will apply to every test.
        for test in self.tests:
            self.tests[test].parse()

    def _parse_childvar(self, child):
        """Adds a <global> tag to the variables dict."""
        if (child.doctype == "global" and "name" in child.attributes):
            #Handle the special case of an override definition for the functions
            #output variable that holds its return value.
            if child.attributes["name"].lower() == "[default]":
                child.attributes["name"] = self.element.name.lower() + "_fpy"
            TestingGroup.global_add(self.variables, self.variable_order,
                                    child.attributes["name"].lower(), child)

    def _get_param_globals(self):
        """Extracts global declarations from the parameter regular="true" attribute
        of paramater decorations."""
        for name in self.element.parameters:
            param = self.element.parameters[name]
            #If the executable we are handling parameters for is an embedded procedure
            #in a custom type, then we want to initialize the variable instance for
            #the type, even if there is no reglar="true"; developers will seldom put
            #a <parameter> tag on the first argument to an embedded procedure because
            #it is always the class instance.
            if self.element.is_type_target and param.kind == self.element.is_type_target.name:
                #is_type_target is the instance of the CustomType that owns it.
                glob = self._global_from_param(name)
                if glob is not None:
                    TestingGroup.global_add(self.variables, self.variable_order, name, glob)
            else:
                #Check the docstrings for a regular="true" attribute.
                for doc in param.docstring:
                    if doc.doctype == "parameter" and \
                       "regular" in doc.attributes and doc.attributes["regular"].lower() == "true":
                        glob = self._global_from_param(name)
                        if glob is not None:
                            TestingGroup.global_add(self.variables, self.variable_order, name, glob)

    def _global_from_param(self, name):
        """Creates a global DocElement from the specified executable's parameter "name"."""
        if name in self.element.parameters:
            param = self.element.parameters[name]
            result = DocElement(None, None)
            result.doctype = "AUTOPARAM"
            result.attributes["name"] = param.name
            result.attributes["type"] = param.dtype
            result.attributes["D"] = param.D
            self._global_clean_param(result, "kind", param.kind)
            self._global_clean_param(result, "modifiers", ", ".join(param.modifiers))
            self._global_clean_param(result, "dimensions", param.dimension)
            self._global_clean_param(result, "default", param.default)
            
            #if the variable is a deferred shape array and we have been told to allocate
            #it during the initializing phase then we must add "allocatable" as a modifier
            #if it isn't already allocatable or a pointer.
            if name in self.assignments and param.D > 0 and ":" in param.dimension:
                allocate = ("allocate" in self.assignments[name] and
                            self.assignments[name]["allocate"] != "false")
                if ("modifiers" in result.attributes and not
                    ("allocatable" in param.modifiers or
                     "pointer" in param.modifiers)):
                    result.attributes["modifiers"] += ", allocatable"
                elif "modifiers" not in result.attributes:
                    result.attributes["modifiers"] = "allocatable"
                    
            #For regular="true" parameters of the class variable this/self, we need to
            #explicitly add the pointer modifier so that the allocation works.
            needsadj = ((param.dtype == "class" or param.dtype == "type") and
                        ("pointer" not in result.attributes["modifiers"] and
                         "allocatable" not in result.attributes["modifiers"]))
            if needsadj:
                if result.attributes["modifiers"] == "":
                    result.attributes["modifiers"] = "pointer"
                else:
                    result.attributes["modifiers"] += ", pointer"

            return result
        else:
            return None

    def _global_clean_param(self, result, name, value):
        """Adds the specified parameter name and value to the result if the value is not None."""
        if value is not None:
            result.attributes[name] = value

    @staticmethod
    def global_add(variables, order, name, doc):
        """Adds a global declaration for the specified DocElement if it isn't already
        represented in the list."""
        if name not in variables:
            variables[name] = GlobalDeclaration(doc)
            order.append(name)
        elif (hasattr(variables[name].element, "doctype") and
              variables[name].element.doctype == "AUTOPARAM"):
            #We can override the existing variable declaration because it was from a
            #parent test group and a global tag takes precedence. AUTOPARAMS are only
            #generated from regular="true" tags, so the [default] <globals> would never
            #be handled in this section.
            variables[name] = GlobalDeclaration(doc)
            order.append(name)            
        else:
            #We need to make sure that it is unique compared to the existing
            #one. If it isn't, stop the execution. If it is, we don't need
            #to add it again.
            existing = variables[name]
            if not existing.ignore and not existing.compare(doc):
                msg.err("variables in the calling namespace have the same name," + \
                        " but different types: \n{}{}".format(existing.element, doc))
                exit(1)

    def _parse_mappings(self):
        """Searches the group for <mapping> tags and adds them to the 
        mapping dict."""
        for child in self.children:
            if (child.doctype == "mapping" and "name" in child.attributes
                and "target" in child.attributes):
                mappings[child.attributes["target"]] = child.attributes["name"]
        
class MethodFinder(object):
    """Class for recursively finding methods for executing pre-req methods
    that takes recursive dependencies into account.

    :arg identifier: the "module.executable" that this method finder represents
    :arg parser: the code parser instance for inter-module access.
    :arg element: the DocElement that specified this method as a pre-req.
    :arg testid: the identifier of the test specification that this method finder is
      being constructed for.
    :arg fatal_if_missing: specified whether the framework chokes if it can't
      find a method that is referenced as a pre-req.
    :attr executable: the code element instance of the executable found from
      the module that the method belongs to.
    :attr main: when true, this instance is for the *main* method being unit tested.
    :attr basic: when true, the Finder only initializes the executable attribute for
      the given identifier, but doesn't recursively find pre-req and other dependency
      methods.
    """
    def __init__(self, identifier, parser, element, testid, fatal_if_missing = True,
                 main=False, basic=False):
        self.identifiers = identifier.lower().split(".")
        self.name = identifier
        self.writekey = self.name
        """Key to uniquely identify this method finder from any other that exists because
        of multiple calls to the same pre-requisite methods."""
        self.methods = []
        self.element = element
        self.testid = testid
        self.main = main
        self.test = None
        """The TestSpecification instance for the unit test that this method finder is
        operating under."""

        self._parser = parser
        self._module = None
        self._fatal = fatal_if_missing

        self.executable = self._find_executable()
        #If we have a result, get a pointer to the tests so we don't have to go one level
        #deep all the time.
        if self.executable is not None:
            self.group = self.executable.test_group
        else:
            self.group = None

        if "repeats" in self.attributes:
            self.repeats = self.attributes["repeats"] == "true"
        elif element is None:
            #The main method is always repeatable; it is only the pre-reqs that can
            #be selected as repeatable or not.
            self.repeats = True
        else:
            self.repeats = False

        if basic or ("terminate" in self.attributes and self.attributes["terminate"] == "true"):
            #'basic' is important so that we can efficiently find the executables with this
            #class but not have the overhead of the recursive find. If they manually specified
            #a termination, we should honor that now as well.
            return

        if self.group is None:
            msg.warn("Executable {} has no testing group; aborting pre-req search".format(identifier))
            return

        #Get the reference to the specific test that this finder is running for.
        if testid is not None:
            if self.testid in self.group.tests:
                self.test = self.group.tests[self.testid]
            else:
                raise ValueError("Unit test {} cannot be found in docstrings.".format(testid))
        else:
            #By default, if there isn't a testid specified, we can just use the
            #first test in the testing group.
            if len(self.group.tests) > 0:
                self.testid = list(self.group.tests.keys())[0]
                self.test = self.group.tests[self.testid]
            else:
                raise ValueError("Cannot auto-detect test identifier for {}.".format(identifier))

        #Recursively get all the pre-requisite for this method. However if the "terminate"
        #attribute is in the doc element, then don't. The main method being unit tested
        #is at the top of the recursive tree of MethodFinders and it is initialized with
        #element==None; all the lower branches in the recursion tree are initialized with
        #their respective doc element instances.
        if element is None or \
           not ("terminate" in self.attributes 
                and self.attributes["terminate"] == "true"):
            self._get_prereqs()

    @property
    def terminate(self):
        """Specifies whether any pre-reqs chained to this method should be executed or not.
        """
        return ("terminate" in self.attributes and self.attributes["terminate"] == "true")        
            
    @property
    def parser(self):
        """Returns the CodeParser instance for inter-modular access."""
        return self._parser
            
    @property
    def variables(self):
        """Returns a list of all the global variables declared by the testing
        group of this method if they exist."""
        #We give preference to test specification variable lists.
        if self.test is not None:
            return self.test.variables
        elif self.group is not None:
            return self.group.variables
        else:
            return []
        
    @property
    def attributes(self):
        """Returns the dictionary of attributes from the underlying xml tag."""
        if self.element is not None:
            if isinstance(self.element, DocElement):
                return self.element.attributes
            else:
                return self.element.attrib
        else:
            return {}                  

    def add_prereq(self, identify, tag, chainto):
        """Adds a MethodFinder instance to the current one for the specified 
        'module.executable' identity.

        :arg identify: a string of "modulename.executablename".
        :arg tag: the DocElement instance that resulted in this prereq being added.
        :arg chainto: the test identifier specified by the pre-req that should be used
          for chaining.
        """
        #The only reason this method exists is to prevent an import loop problem
        #It gets used by the AssignmentValue._code_embedded().
        self.methods.append(MethodFinder(identify, self._parser, tag, chainto, self._fatal))

    def _get_prereqs(self):
        """Compiles a list of subroutines that need to run. Handles recursive calls."""
        #Make sure we have tests to run off of; some of the routines that are dependencies
        #won't necessarily have any docstrings.
        for method in self.test.methods:
            if isinstance(method, TestPreReq):
                self.methods.append(MethodFinder(method.method, self._parser, method.xml,
                                                 method.chainto, self._fatal))
            elif isinstance(method, Assignment):
                #This is a variable assignment. Sometimes variable values are changed
                #between calls to functions. The order in which prereq/instance elements
                #appear defines when these assignments take place. We will treat a var
                #value assignment as a method to simplify the implementation.
                self.methods.append(method)

        #Now that we have the method tree setup; we need to check each of the Assignments
        #to get their pre-reqs straightened out.
        for m in self.methods:
            if isinstance(m, Assignment):
                m.check_prereqs(self)
        
    def timed(self, testid):
        """Returns true if this method will produce code to time its execution."""
        if self.group is not None and self.main:
            testspec = self.group.tests[testid]
            return testspec.timed
        else:
            return False        

    def code(self, lines, position, spacer, testid):
        """Appends the code to execute *only* this method's main executable as
        part of the main program.

        :arg lines: the list of strings that form the contents of the *.f90 file
          being created.
        """
        #We only need to worry about vars if this is the *main* method that is
        #being unit tested.
        if self.main and position == "vars":
            self._code_vars(lines, spacer, testid)

        if position == "call":
            self._code_call(lines, spacer, testid)

        if self.timed(testid) and position == "final":
            lines.append("{}call pysave(fpy_elapsed, 'fpytiming.out')".format(spacer))

    def _code_vars(self, lines, spacer, testid):
        """Appends code to initialize the variables that accept the return values
        if the method being called is a function. Also appends the variables needed
        to time the method's execution if specified by the test.
        """
        #This method only gets called if we are the main executable being tested.
        if type(self.executable).__name__ == "Function":
            if (self.executable.name + "_fpy").lower() not in self.variables:
                lines.append("{}{}".format(spacer, self.executable.definition("_fpy")))

        if self.timed(testid):
            lines.append("{}real(fdp) :: fpy_start, fpy_end, fpy_elapsed = 0".format(spacer))

    def _code_call(self, lines, spacer, testid):
        """Appends the code to call the executable that this method finder represents."""
        #If we are timing the executable, we want to get the time before and after the
        #execution of just this method and then add it to the elapsed time.
        lines.append("")
        if self.timed(testid):
            lines.append("{}call cpu_time(fpy_start)".format(spacer))

        #This is either a call to a subroutine or a function. We need to make
        #sure that we handle any mapping tags or call tags in the docstrings
        if type(self.executable).__name__ == "Subroutine":
            prefix = "call " 
        else:
            #For a function, we still need to save the value somewhere so we
            #can compare it.
            pntr = ">" if "pointer" in [m.lower() for m in self.executable.modifiers] else ""
            prefix = "{}_fpy ={} ".format(self.executable.name, pntr)

        spacing = len(list(prefix)) + len(list(self.executable.name)) + len(spacer)

        #Unfortunately, the fpy_auxiliary module freaks out if we call a public
        #method directly which is *also* a module procedure for the type (since
        #the .mod files vary because of the extra public subroutine on the type,
        #the compiler dies with Abort trap 6 error. So, we need to call embedded
        #methods using the %-syntax.
        if not self.executable.is_type_target:
            callname = self.executable.name
            xtype = None
        else:
            xtype = self.executable.is_type_target
            for emname, emexec in xtype.executables.items():
                if emexec.target is self.executable:
                    callname = emname
                    break
            else:
                callname = self.executable.name
        
        if "paramlist" in self.attributes:
            #The developer has decided on an alternate call signature from
            #the one that gets auto-generated.
            specified = re.split(",\s*", self.attributes["paramlist"])
            if xtype is not None and len(specified) == len(self.executable.ordered_parameters):
                msg.warn("Explicit parameter list for embedded executable should *exclude* the "
                         "reference to 'self' in the first argument.")
            cleaned = present_params(specified, spacing, 90)
            lines.append("{}{}{}({})".format(spacer, prefix, callname, cleaned))
        else:
            #We can construct the actual list of parameters to use in the call
            #and substitute mappings where appropriate.
            calllist = []
            for ip, param in enumerate(self.executable.ordered_parameters):
                #Because of ignorable parameters, if the parameter is optional, explicitly
                #specify its name.
                if "optional" in param.modifiers:
                    optstr = "{}=".format(param.name)
                else:
                    optstr = ""
                    
                if self.group is not None and param in self.group.mappings:
                    #The first parameter name will also be the name of the variable that gets
                    #created and called.
                    if xtype is not None and ip == 0:
                        callname = "{}%{}".format(self.group.mappings[param], callname)
                    else:
                        calllist.append(optstr + self.group.mappings[param])
                else:
                    var = None
                    pname = param.name.lower()
                    #The test specification takes precedence over the testing group for variables.
                    if self.test is not None and pname in self.test.variables:
                        var = self.test.variables[pname]
                    if var is None and self.group is not None and pname in self.group.variables:
                        var = self.group.variables[pname]
                    
                    if xtype is not None and ip == 0:
                        callname = "{}%{}".format(param.name, callname)
                    else:
                        if var is not None and not var.ignore:
                            calllist.append(optstr + param.name)
                        elif var is None:
                            calllist.append(optstr + param.name)

            lines.append("{}{}{}({})".format(spacer, prefix, callname,
                                              present_params(calllist, spacing, 90)))        

        if self.timed(testid):
            lines.append("{}call cpu_time(fpy_end)".format(spacer))
            lines.append("{}fpy_elapsed = fpy_elapsed + fpy_end - fpy_start".format(spacer))
                
    def _find_executable(self):
        """Finds the executable code element from the code parser instance."""
        #First, we need to determine which module the identifier points to.
        #See if the module is in the code parser
        result = None

        if self.module is not None:
            #We want to use the code parsers tree search to get the executable
            #If the identifier refers to a type method, we need to use the
            #type search
            if "%" in self.identifiers[1]:
                tbase = self.identifiers[1].split("%")[0].strip()
                element = self._parser.type_search(tbase, self.identifiers[1], self._module)
                if isinstance(element, Executable):
                    result = element                   
            else:
                #This is just a standard exectable inside the module
                if self.identifiers[1] in self._module.executables:
                    result = self._module.executables[self.identifiers[1]]
                elif self._fatal:
                    raise ValueError("FATAL: a specified executable was not found in" + 
                                     "the module: {}".format(self.identifiers[1]))

        return result

    @property
    def module(self):
        """Returns the module instance for this method."""
        if self._module is None:
            if not self.identifiers[0] in self._parser.modules:
                #Try a dependency search for the specific module
                self._parser.load_dependency(self.identifiers[0], True, True, False)
                #If the module was found, we can set the value of the local variable
                #otherwise we will end the program by default
                if self.identifiers[0] in self._parser.modules:
                    self._module = self._parser.modules[self.identifiers[0]]
                elif self._fatal:
                    raise ValueError("FATAL: the module for a pre-requisite method could " + 
                                     "not be located: {}".format(self.identifiers[0]))

            if self.identifiers[0] in self._parser.modules:
                self._module = self._parser.modules[self.identifiers[0]]
        
        return self._module        

