import re
from . import msg
from fortpy.testing.elements import TestingGroup

#This module has all the classes for holding the structure of a fortran
#code file and its docstrings.
class CodeElement(object):
    """Represents a code element (e.g. parameter, module, subroutine) etc."""
    
    def __init__(self, name, modifiers, parent):
        self.name = name
        self.modifiers = modifiers
        #Docstring is a list of DocElement intstances describing the code
        self.docstring = []
        #groups is a list of DocGroup() that contain the groupings of the DocElements
        #within the code element.
        self.groups = {}
        self.parent = parent
        #The start and end character indices of the code element definition in the file
        self.start = 0
        self.end = 0

        self._tests = None
        self._testgroup = None
        self._module = None
        self._full_name = None
        self._summary = None

        #If the modifiers passed in is None, set it to an empty list
        #sometimes when there is no regex match on the modifiers we get None
        if self.modifiers is None:
            self.modifiers = []
        else:
            self.clean_mods(self.modifiers)

    @property
    def lname(self):
        """Returns the lowered name of the code element.
        """
        return self.name.lower()
            
    def overwrite_docs(self, doc):
        """Adds the specified DocElement to the docstring list. However, if an
        element with the same xml tag and pointsto value already exists, it
        will be overwritten."""
        for i in range(len(self.docstring)):
            if (self.docstring[i].doctype == doc.doctype and
                self.docstring[i].pointsto == doc.pointsto):
                del self.docstring[i]
                break

        self.docstring.append(doc)

    def __getstate__(self):
        """Cleans up the object so that it can be pickled without any pointer
        issues. Saves the full names of parents etc. so that they can be 
        restored once the unpickling is completed."""
        odict = self.__dict__.copy() # copy the dict since we change it
        del odict['_tests']
        del odict['_testgroup']
        del odict['_module']
        del odict['parent']
        del odict['_summary']
        return odict

    def __setstate__(self, dict):
        self._tests = None
        self._testgroup = None
        self._module = None
        self._summary = None
        self.parent = None
        self.__dict__.update(dict)

    def unpickle_docs(self):
        """Sets the pointers for the docstrings that have groups."""
        for doc in self.docstring:
            if (doc.parent_name is not None and 
                doc.parent_name in self.groups):
                doc.group = self.groups[doc.parent_name]                

    def unpickle(self, parent):
        """Sets the parent pointer references for the type executable."""
        self.parent = parent
        self.unpickle_docs()

    @property
    def embedded(self):
        """Value indicates whether this type declaration is embedded in an executable
        rather than the module, which is the natural default."""
        return not isinstance(self.parent, Module)

    @property
    def absstart(self):
        """Returns the absolute start of the element by including docstrings
        outside of the element definition if applicable."""
        if hasattr(self, "docstart") and self.docstart > 0:
            return self.docstart
        else:
            return self.start

    @property
    def module(self):
        """Returns the module that this code element belongs to."""
        if self._module is None:
            root = self
            while self._module is None and root is not None:
                if isinstance(root, Module):
                    self._module = root
                else:
                    root = root.parent

        return self._module

    @property
    def summary(self):
        """Returns the docstring summary for the code element if it exists."""
        if self._summary is None:
            self._summary = "No summary for element."
            for doc in self.docstring:
                if doc.doctype == "summary":
                    self._summary = doc.contents
                    break

            #If a parameter, member or local tag has dimensions or other children,
            #then the inner-text is not the right thing to use; find a grand-child
            #summary tag instead.
            if self._summary == "No summary for element." and len(self.docstring) > 0:
                summary = self.doc_children("summary")
                if len(summary) > 0:
                    self._summary = summary[0].contents
                else:
                    self._summary = self.docstring[0].contents

        return self._summary

    @property
    def full_name(self):
        """Returns the full name of this element by visiting every
        non-None parent in its ancestor chain."""
        if self._full_name is None:
            ancestors = [ self.name ]
            current = self.parent
            while current is not None and type(current).__name__ != "CodeParser":
                ancestors.append(current.name)
                current = current.parent

            self._full_name = ".".join(reversed(ancestors))
        
        return self._full_name

    @property
    def tests(self):
        """Returns a the contents of a group with purpose="testing" if it exists."""
        if self._tests is None:
            self._tests = []
            docgrp = self.test_group
            if docgrp is not None:
                for docel in self.docstring:
                    if docel.group == docgrp.name:
                        self._tests.append(docel)

        return self._tests

    @property 
    def test_group(self):
        """Returns the doc group with purpose="testing" if it exists."""
        if self._testgroup is None:
            for gkey in self.groups:
                docgrp = self.groups[gkey]
                if "purpose" in docgrp.attributes and \
                   docgrp.attributes["purpose"].lower() == "testing":
                    self._testgroup = TestingGroup(docgrp, self)

        return self._testgroup

    @property
    def has_docstring(self):
        """Specifies whether this code element has a docstring."""
        return type(self.docstring) != type(None)        

    def doc_children(self, doctype, limiters=[]):
        """Finds all grand-children of this element's docstrings that match
        the specified doctype. If 'limiters' is specified, only docstrings
        with those doctypes are searched.
        """
        result = []
        for doc in self.docstring:
            if len(limiters) == 0 or doc.doctype in limiters:
                result.extend(doc.children(doctype))

        return result

    def warn(self, collection):
        """Checks this code element for documentation related problems."""
        if not self.has_docstring():
            collection.append("WARNING: no docstring on code element {}".format(self.name))

    def clean_mods(self, modifiers):
        """Cleans the modifiers to remove empty entries."""
        if "" in modifiers and isinstance(modifiers, list):
            modifiers.remove("")

class ValueElement(CodeElement):
    """Represents a code element that can hold a value."""

    def __init__(self, name, modifiers, dtype, kind, default, dimension, parent, D=0):
        super(ValueElement, self).__init__(name, modifiers, parent)
        self.dtype = dtype
        self.kind = self._clean_args(kind)
        self.default = default
        self.dimension = dimension
        self.D = D
        """Returns the integer number of dimensions that this variable
        is declared as having."""

        self._direction = None
        self._ctypes_parameter = None
        self._kind_module = None
        """This is the instance of fortpy.elements.Module that has the definition
        for the kind of this value element."""

    def __str__(self):
        return self.definition()

    # def __eq__(self, other):
    #     if not isinstance(other, ValueElement):
    #         return False
    #     mods = ["allocatable", "pointer"]
    #     return (self.name.lower() == other.name.lower() and
    #             self.kind.lower() == other.kind.lower() and
    #             self.dtype.lower() == other.dtype.lower() and
    #             self.D == other.D and
    #             all([m in other.modifiers for m in self.modifiers if m in mods]))              

    def matched(self, other):
        """Returns True if the two ValueElement instances differ only by name,
        default value or some other inconsequential modifier.
        """
        mods = ["allocatable", "pointer"]
        return (self.kind.lower() == other.kind.lower() and
                self.dtype.lower() == other.dtype.lower() and
                self.D == other.D and
                all([m in other.modifiers for m in self.modifiers if m in mods]))              
    
    @property
    def strtype(self):
        """Returns a string representing the type and kind of this value element."""
        if self.kind is not None:
            return "{}({})".format(self.dtype, self.kind)
        else:
            return self.dtype

    def _get_ctypes_name(self, index=None):
        """Returns a formatted name for the ctypes Fortran wrapper module.

        :arg index: the index of the array to return a name for. None for the new parameter,
          otherwise the index for the integer variable name of that dimension.
        """
        if index is None:
            if ("allocatable" not in self.modifiers and "pointer" not in self.modifiers
                and self.dtype != "logical"):
                #The fortan logical has type 4 by default, whereas c_bool only has 1
                return self.name
            else:
                return "{}_c".format(self.name)
        else:
            return "{0}_{1:d}_c".format(self.name, index)
        
    def ctypes_parameter(self):
        """Returns the parameter list for this ValueElement adjusted for interoperability
        with the ctypes module.
        """
        if self._ctypes_parameter is None:
            #Essentially, we just need to check if we are an array that doesn't have explicitly
            #defined bounds. Assumed-shape arrays have to be 'pointer' or 'allocatable'. However,
            #the deffered/assumed shape arrays always use ':' as the array dimension.
            if self.dimension is not None and ":" in self.dimension:
                result = [self._get_ctypes_name()]
                result.extend([self._get_ctypes_name(i+1) for i in range(self.D)])
                if self.direction == "(inout)" and ("allocatable" in self.modifiers or
                                                    "pointer" in self.modifiers):
                    result.append("{}_o".format(self.name))
                self._ctypes_parameter = result
            elif self.dtype == "logical":
                self._ctypes_parameter = [self._get_ctypes_name()]
            elif hasattr(self, "parameters"):
                #This is the handler for the function return types that are simple
                self._ctypes_parameter = [self.name + "_f"]
            else:
                #Nothing special, just return the standard name.
                self._ctypes_parameter = [self.name]

        return self._ctypes_parameter
        
    @property
    def direction(self):
        """Returns the direction of the variable if it is a parameter. Possible values
        are ["": no intent, "(in)", "(out)", "(inout)"].
        """
        if self._direction is None:
            if hasattr(self, "parameters"):
                #This is actually a function that inherited from ValueElement, it is always
                #intent(out) for our purposes.
                self._direction = "(out)"
            else:
                intent = [m for m in self.modifiers if "intent" in m]
                if len(intent) > 0:
                    self._direction = intent[0].replace("intent", "").strip()
                else:
                    self._direction = ""

        return self._direction

    def definition(self, suffix = "", local=False, ctype=None, optionals=True,
                   customdim=None, modifiers=None):
        """Returns the fortran code string that would define this value element.

        :arg suffix: an optional suffix to append to the name of the variable.
          Useful for re-using definitions with new names.
        :arg local: when True, the parameter definition is re-cast as a local
          variable definition that has the "intent" and "optional" modifiers removed.
        :arg ctype: if a ctype should be used as the data type of the variable
          instead of the original type, specify that string here.
        :arg optionals: removes the "optional" modifier from the definition before
          generating it.
        :arg customdim: if the dimension string needs to be changed, specify the
          new one here.
        :arg modifiers: specify an additional list of modifiers to add to the
          variable definition.
        """
        kind = "({})".format(self.kind) if self.kind is not None else ""   
        cleanmods = [m for m in self.modifiers if m != "" and m != " "
                     and not (local and ("intent" in m or m == "optional"))
                     and not (not optionals and m == "optional")]
        if modifiers is not None:
            cleanmods.extend(modifiers)
        if len(cleanmods) > 0:
            mods = ", " + ", ".join(cleanmods) + " " 
        else:
            mods = " "
        if customdim is not None:
            dimension = "({})".format(customdim)
        else:
            dimension = "({})".format(self.dimension) if self.dimension is not None else ""

        if self.default is None:
            default = ""
        else:
            if ">" in self.default: #We have a pointer, don't add an extra space.
                default = " ={}".format(self.default) if self.default is not None else ""
            else:
                default = " = {}".format(self.default) if self.default is not None else ""
        name = "{}{}".format(self.name, suffix)
        stype = self.dtype if ctype is None else ctype
        return "{}{}{}:: {}{}{}".format(stype, kind, mods, name, dimension, default)

    def _clean_args(self, arg):
        """Removes any leading and trailing () from arguments."""
        if arg is not None:
            return re.sub("[()]", "", arg)
        else:
            return None

    @property
    def py_initval(self):
        """Returns the python value that would initialize a python variable representing
        this parameter in the standard way.
        """
        valdict = {
            "real": "0.0",
            "integer": "0",
            "complex": "0+0j",
            "character": "''",
            "logical": "False"
        }
        if self.dtype in valdict:
            return valdict[self.dtype]
        else:
            return ""
        
    @property
    def argtypes(self):
        """Returns the ctypes argtypes for use with the method.argtypes assignment for
        an executable loaded from a shared library.
        """
        if self.dimension is not None:
            result = []
            if "in" in self.direction:
                #The only complication here is that the 'known' dimensionality could actually
                #be a function like "size" that needs information about other variables.
                #If we choose to ignore known shapes, we lose the error checking for the passed
                #in variables from the python side.
                if self.direction == "(inout)" and ":" not in self.dimension:
                    wstr = ", writeable"
                else:
                    wstr = ""

                if ":" in self.dimension or "size" in self.dimension:
                    template = 'ndpointer(dtype={}, ndim={}, flags="F{}")'
                    result.append(template.format(self.pytype, self.D, wstr))
                else:
                    template = 'ndpointer(dtype={}, ndim={}, shape=({}), flags="F{}")'
                    sdim = self.dimension + ("" if self.D > 1 else ",")
                    result.append(template.format(self.pytype, self.D, sdim, wstr))
            elif self.direction == "(out)":
                result.append("c_void_p")

            if self.D > 0 and ":" in self.dimension:
                result.extend(["c_int_p" for i in range(self.D)])
            if (self.direction == "(inout)" and ":" in self.dimension and
                ("allocatable" in self.modifiers or "pointer" in self.modifiers)):
                result.append("c_void_p")
            return result
        else:
            ctype = self.ctype
            if ctype is not None:
                return ["{}_p".format(ctype.lower())]
        
    @property
    def pytype(self):
        """Returns the python type of the underlying parameter or variable.
        """
        lookup = {
            "logical": "bool",
            "real": "float",
            "integer": "int",
            "complex": "complex",
            "character": "str"
        }
        if self.dtype in lookup:
            return lookup[self.dtype]
        else:
            return None
        
    @property
    def ctype(self):
        """Returns the name of the c_type from iso_c_binding to use when declaring
        the output parameter for interaction with python ctypes.
        """
        if self.dtype == "logical":
            return "C_BOOL"
        elif self.dtype == "complex":
            #We don't actually know what the precision of the complex numbers is because
            #it is defined by the developer when they construct the number with CMPLX()
            #We just return double to be safe; it is a widening conversion, so there
            #shouldn't be any issues.
            return "C_DOUBLE_COMPLEX"
        elif self.dtype == "character":
            return "C_CHAR"
        elif self.dtype in ["integer", "real"]:
            if self.kind is None:
                if self.dtype == "integer":
                    return "C_INT"
                else:
                    return "C_FLOAT"
                
            if self._kind_module is None and self.kind is not None:
                self.dependency()
            if self._kind_module is None and self.kind is not None:
                raise ValueError("Can't find the c-type for {}".format(self.definition()))
            elif self._kind_module is not None:
                #We look up the parameter in the kind module to find out its
                #precision etc.
                import re
                default = self._kind_module.members[self.kind].default
                vals = default.split("(")[1].replace(")", "")
                ints = list(map(int, re.split(",\s*", vals)))
                if self.dtype == "integer" and len(ints) == 1:
                    if ints[0] <= 15:
                        return "C_SHORT"
                    elif ints[0] <= 31:
                        return "C_INT"
                    elif ints[0] <= 63:
                        return "C_LONG"
                elif self.dtype == "real" and len(ints) == 2:
                    if ints[0] <= 24 and ints[1] < 127:
                        return "C_FLOAT"
                    elif ints[0] <= 53 and ints[1] < 1023:
                        return "C_DOUBLE"
                    
    def dependency(self):
        """Returns the module name from the code parser that this parameter/variable
        needs to compile successfully.
        """
        from re import match
        if self._kind_module is not None:
            #We have already done this search at some point.
            return self._kind_module.name

        if (self.kind is not None and "len=" not in self.kind.lower()
            and not match("\d", self.kind[0])):
            #Find the module that declares this kind as:
            # 1) derived type
            # 2) parameter/member
            #The good news is that in order for the module being unit tested to compile,
            #it must have a 'use' statement for the derived type. We can just search from
            #that module with a tree find. Also, at the end of the day, it will have to
            #be declared as public, in order to be used in other modules.
            found, foundmod = self.module.parent.tree_find(self.kind.lower(),
                                                           self.module,
                                                           "publics")
            self._kind_module = foundmod #At worst this will be None again
            if found is not None:
                return foundmod.name

    @property
    def customtype(self):
        """If this variable is a user-derivedy type, return the CustomType instance that
        is its kind.
        """
        result = None
        if self.is_custom:
            #Look for the module that declares this variable's kind in its public list.
            self.dependency()
            if self._kind_module is not None:
                if self.kind.lower() in self._kind_module.types:
                    result = self._kind_module.types[self.kind.lower()]

        return result
            
    @property
    def is_custom(self):
        """Returns a value indicating whether this value element is of a derived type."""
        return self.dtype == "class" or self.dtype == "type"

class Dependency(object):
    """Represents a call to a function or subroutine from within another
    thus making one executable dependent on the others."""
    def __init__(self, name, argslist, isSubroutine, parent):
        self._name = name
        self.argslist = self.clean(argslist)
        self.parent = parent
        self.isSubroutine = isSubroutine
        
        self._module = None
        self._target = None

    def __str__(self):
        if self.isSubroutine:
            call = "call "
        else:
            call = ""
        return "\t{}{}({})".format(call, self._name, ", ".join(self.argslist))

    def unpickle(self, parent):
        """Sets the parent pointer references for the type executable."""
        self.parent = parent

    @property
    def argnames(self):
        """Returns a list of the variable names being passed to the dependency
        as arguments.
        """
        return [re.split("[(]", a)[0].lower() for a in self.argslist]
        
    @property
    def operator(self):
        """Returns true if this dependency is to an operator interface."""
        return self._name[0] == self._name[-1] and self._name[0] == "."

    @property
    def name(self):
        """Returns the lower case name of the dependency."""
        return self._name.lower()

    @property
    def fullname(self):
        """Returns the original name of the dependency as found in the code."""
        return self._name
        
    @property
    def external_name(self):
        """Returns the modulename.executable string that uniquely identifies
        the executable that this dependency points to."""
        target = self.target
        if target is not None:
            return "{}.{}".format(target.name.lower(), self.name)
        else:
            return "{}.{}".format(self.module.name.lower(), self.name)

    @property
    def target(self):
        """Returns the executable code element that this dependency points to
        if it can be found.
        """
        if self._target is None:
            if '%' in self.name:
                parts = self.name.split('%')
                base = self.module.parent.type_search(parts[0], self.name, self.module)
                if base is not None:
                    self._target = base.target
                else:
                    self._target = False
            else:
                found, foundmod = self.module.parent.tree_find(self.name, self.module, "executables")
                if found is not None:
                    self._target = found
                else:
                    self._target = False

        if isinstance(self._target, Executable):
            return self._target
        else:
            return None
    
    @property
    def module(self):
        """Returns the module that this code element belongs to."""
        if self._module is None:
            root = self
            while self._module is None and root is not None:
                if isinstance(root, Module):
                    self._module = root
                else:
                    root = root.parent

        return self._module

    def clean(self, argslist):
        """Cleans the argslist."""
        result = []
        for arg in argslist:
            if type(arg) == type([]):
                if len(result) > 0:
                    result[-1] = result[-1] + "(*{})".format(len(self.clean(arg)))
                elif "/" not in arg[0]:
                    msg.warn("argument to function call unrecognized. {}".format(arg))
            else:
                cleaner = re.sub("[:,]+", "", arg).strip()
                if len(cleaner) > 0:
                    result.append(cleaner)
        
        return result

class Decoratable(object):
    """Represents a class that can have an *external* docstring attached."""
    def __init__(self):
        #The start and end characters for the docstring that decorates this code element
        self.docstart = 0
        self.docend = 0

    def find_section(self, charindex):
        """Returns a value indicating whether the specified character index
        is owned by the current object."""
        #All objects instances of decorable also inherit from CodeElement,
        #so we should have no problem accessing the start and end attributes.
        result = None
        if hasattr(self, "start") and hasattr(self, "end"):
            #The 8 seems arbitrary, but it is the length of type::b\n for a
            #really short type declaration with one character name.
            if charindex > self.docend and charindex - self.start < 8:
                result = "signature"
            elif charindex >= self.start and charindex <= self.end:
                result = "body"

        if (result is None and charindex >= self.docstart 
            and charindex <= self.docend):
            result = "docstring"
            
        return result
    
class Executable(ValueElement, Decoratable):
    """Represents a function or subroutine that can be executed with parameters."""
    def __init__(self, name, modifiers, dtype, kind, default, dimension, parent):
        super(Executable, self).__init__(name, modifiers, dtype, kind, 
                                         default, dimension, parent)
        Decoratable.__init__(self)
        self.members = {}
        self.dependencies = {}
        #Initialize dicts for the embedded types and executables.
        self.types = {}
        self.executables = {}

        #The order in which the parameters are presented to the function
        self.paramorder = []
        #The string between the end of the signature and the start of the end
        #token for this executable.
        self.contents = None
        #When an instance is add from just a signature and doesn't have an 
        #end_token, this is set to True
        self.incomplete = False
        self._parameters = {}
        self._assignments = []
        #Lazy assignment for this variable; see the property with the similar name.
        self._is_type_target = None

    def search_dependencies(self):
        """Returns a list of modules that this executable needs in order to run
        properly. This includes special kind declarations for precision or derived
        types, but not dependency executable calls.
        """
        #It is understood that this executable's module is obviously required. Just
        #add any additional modules from the parameters.
        result = [p.dependency() for p in self.ordered_parameters]
        result.extend([v.dependency() for k, v in list(self.members.items())])
        for ekey, anexec in list(self.executables.items()):
            result.extend(anexec.search_dependencies())
        return [m for m in result if m is not None and m != self.module.name]
        
    @property
    def primitive(self):
        """Returns True if this executable only has parameters with standard types
        (i.e. no user-derived types in the parameter list).
        """
        #Just look through the types of each parameter to see if any have
        #'type' or 'class'.
        return not any([p.is_custom for p in self.ordered_parameters])
        
    def __getstate__(self):
        """Cleans up the object so that it can be pickled without any pointer
        issues.
        """
        odict = self.__dict__.copy() # copy the dict since we change it
        del odict["_is_type_target"]
        return odict

    def __setstate__(self, dict):
        self._is_type_target = None
        self.__dict__.update(dict)

    def unpickle(self, parent):
        """Sets the parent pointer references for the module *and* all of its
        child classes that also have pointer references."""
        self.parent = parent
        self._unpickle_collection(self.members)
        self._unpickle_collection(self.dependencies)
        self._unpickle_collection(self.types)
        self._unpickle_collection(self.executables)
        self._unpickle_collection(self._parameters)
        self.unpickle_docs()
        
    def _unpickle_collection(self, collection):
        """Unpickles all members of the specified dictionary."""
        for mkey in collection:
            if isinstance(collection[mkey], list):
                for item in collection[mkey]:
                    item.unpickle(self)
            else:
                collection[mkey].unpickle(self)

    def rt_update(self, statement, linenum, mode, xparser):
        """Uses the specified line parser to parse the given line.

        :arg statement: a string of lines that are part of a single statement.
        :arg linenum: the line number of the first line in the list relative to
          the entire module contents.
        arg mode: either 'insert', 'replace' or 'delete'
        :arg xparser: an instance of the executable parser from the real
          time update module's line parser.
        """
        section = self.find_section(self.module.charindex(linenum, 1))

        if section == "body":
            xparser.parse_line(statement, self, mode)
        elif section == "signature":
            if mode == "insert":
                xparser.parse_signature(statement, self)
        #NOTE: docstrings are handled at a higher level by the line parser
        #because of the XML dependence.

    def update_name(self, name):
        """Changes the name of this executable and the reference to it in the
        parent module."""
        if name != self.name:
            self.parent.executables[name] = self
            del self.parent.executables[self.name]
            self.name = name

    @property
    def is_type_target(self):
        """Returns the CustomType instance if this executable is an embedded 
        procedure in a custom type declaration; else False.
        """
        if self._is_type_target is None:
            #All we need to do is search through the custom types in the parent
            #module and see if any of their executables points to this method.
            self._is_type_target = False
            for tkey in self.module.types:
                custype = self.module.types[tkey]
                for execkey in custype.executables:
                    if custype.executables[execkey].target is self:
                        self._is_type_target = custype
                        break
                if self._is_type_target:
                    break

        return self._is_type_target

    @property
    def refstring(self):
        """The string from which this executable was extracted with regex. It is the
        section after the CONTAINS statement in the module."""
        if self.parent is not None:
            return self.parent.contains
        else:
            return ""

    def changed(self, symbol, checked = None):
        """Returns true if the specified symbol has it's value changed
        by this executable or *any* of its dependencies."""
        #Initialize the dictionary if we are the first executable to call this.
        myname = "{}.{}".format(self.module.name, self.name).lower()
        if checked is None:
            checked = []

        if myname not in checked:
            matches = self._get_assignments_in(self._parameters, symbol)
            if (matches is not None and
                (isinstance(matches, bool) and matches == True) or
                (isinstance(matches, list) and len(matches) > 0)):
                return myname
            else:
                #Make sure we don't check any executable twice.
                checked.append(myname)
                #We need to check the dependencies of this executable and
                #see if any of them modify the parameter.
                for dependlist in self.dependencies:
                    for dependency in self.dependencies[dependlist]:
                        if dependency.operator:
                            continue

                        if symbol in dependency.argnames:
                            pindex = dependency.argnames.index(symbol)
                            iexec = dependency.target
                            if iexec is not None:
                                pname = iexec.ordered_parameters[pindex].name.lower()
                                if iexec.changed(pname, checked) != "":
                                    return iexec.full_name
                            else:
                                checked.append(dependency.external_name)
                        else:
                            checked.append(dependency.external_name)
        else:
            return None

    def local_assignments(self):
        """Returns a list of local variable code elements whose values change in
        this executable."""
        return self._get_assignments_in(self.members)

    def external_assignments(self):        
        """Returns a list of parameter code elements whose values change in
        this executable."""
        return self._get_assignments_in(self._parameters)

    def _get_assignments_in(self, filterlist, symbol = ""):
        """Returns a list of code elements whose names are in the specified object.

        :arg filterlist: the list of symbols to check agains the assignments.
        :arg symbol: when specified, return true if that symbol has its value
          changed via an assignment."""
        if symbol != "":
            lsymbol = symbol
            for assign in self._assignments:
                target = assign.split("%")[0].lower()
                if target == lsymbol:
                    return True
        else:
            result = []
            for assign in self._assignments:
                target = assign.split("%")[0].lower()
                if target in filterlist:
                    result.append(assign)
            return result

    @property
    def assignments(self):
        """Returns a list of the names of all the objects whose values change
        in this executable."""
        return self._assignments

    @property
    def ordered_parameters(self):
        """Returns a list of the ordered parameters."""
        return [ self._parameters[k] for k in self.paramorder]

    @property
    def parameters(self):
        """Returns the dictionary of parameters in this exectuable."""
        return self._parameters
    
    def get_parameter(self, index):
        """Returns the ValueElement corresponding to the parameter
        at the specified index."""
        result = None
        if index < len(self.paramorder):
            key = self.paramorder[index]
            if key in self._parameters:
                result = self._parameters[key]

        return result

    def add_parameter(self, parameter):
        """Adds the specified parameter value to the list."""
        if parameter.name.lower() not in self.paramorder:
            self.paramorder.append(parameter.name.lower())
        self._parameters[parameter.name.lower()] = parameter

    def remove_parameter(self, parameter_name):
        """Removes the specified parameter from the list."""
        if parameter_name in self.paramorder:
            index = self.paramorder.index(parameter_name)
            del self.paramorder[index]

        if parameter_name in self._parameters:
            del self._parameters[parameter_name]

    def parameters_as_string(self):
        """Returns a comma-separated list of the parameters in the executable definition."""
        params = ", ".join([ p.name for p in self.ordered_parameters ])
        return params

    def add_assignment(self, value):
        """Adds the name of a variable/parameter whose value is changed by
        this exectuable."""
        if not value in self._assignments:
            self._assignments.append(value)

    def add_dependency(self, value):
        """Adds the specified executable dependency to the list for this executable."""
        if value.name in self.dependencies:
            self.dependencies[value.name.lower()].append(value)
        else:
            self.dependencies[value.name.lower()] = [ value ]
            
class Function(Executable):
    """Represents a function in a program or module."""    
    def __init__(self, name, modifiers, dtype, kind, parent):
        super(Function, self).__init__(name, modifiers, dtype, kind, None, None, parent)
        self.update_dtype()
        
    def __str__(self):
        params = self.parameters_as_string()
        
        depend = "{} dependencies ".format(len(list(self.dependencies.keys())))
        if len(list(self.dependencies.keys())) == 0:
            depend = ""
        assign = "{} assignments".format(len(self.external_assignments()))
        if len(self.external_assignments()) == 0:
            assign = ""
        if depend != "" or assign != "":
            info = "\n\t  - {}{}".format(depend, assign)
        else:
            info = ""

        return "{}FUNCTION {}({}){}".format(self.returns, self.name, 
                                                    params, info)

    @property
    def signature(self):
        """Returns the signature definition for the function."""
        return "{}FUNCTION {}({})".format(self.returns, self.name,
                                          self.parameters_as_string())
        
    @property
    def end_token(self):
        """Gets the end [code type] token for this instance."""
        return "end function"

    def update(self, name, modifiers, dtype, kind):
        """Updates the attributes for the function instance, handles name changes
        in the parent module as well."""
        self.update_name(name)
        self.modifiers = modifiers
        self.dtype = dtype
        self.kind = kind
        self.update_dtype()

    def update_dtype(self, resvar=None):
        """Updates the dtype attribute of the function. This is required because
        fortran functions can have their types declared either as a modifier on
        the function *or* as a member inside the function.

        :arg resvar: the name of the variable declared using the result(var)
          construct after the function signature.
        """
        if self.dtype is None:
            #search the members of this function for one that has the same name
            #as the function. If it gets found, overwrite the dtype, kind and
            #modifiers attributes so the rest of the code works.
            for m in self.members:
                if m == self.name.lower() or m == resvar:
                    member = self.members[m]
                    self.dtype = member.dtype
                    self.modifiers = member.modifiers
                    self.kind = member.kind
                    self.default = member.default
                    self.dimension = member.dimension
                    del self.members[m]
                    break

    @property
    def returns(self):
        """Gets a string showing the return type and modifiers for the
        function in a nice display format."""
        kind = "({}) ".format(self.kind) if self.kind is not None else ""      
        mods = ", ".join(self.modifiers) + " "
        dtype = self.dtype if self.dtype is not None else ""
        return "{}{}{}".format(dtype, kind, mods)

class Subroutine(Executable):
    """Represents a function in a program or module."""   
    def __init__(self, name, modifiers, parent):
        super(Subroutine, self).__init__(name, modifiers, None, None, None, None, parent)
    
    def __str__(self):
        params = self.parameters_as_string()
        mods = ", ".join(self.modifiers)

        depend = "{} dependencies ".format(len(list(self.dependencies.keys())))
        if len(list(self.dependencies.keys())) == 0:
            depend = ""
        assign = "{} assignments".format(len(self.external_assignments()))
        if len(self.external_assignments()) == 0:
            assign = ""
        if depend != "" or assign != "":
            info = "\n\t  - {}{}".format(depend, assign)
        else:
            info = ""

        return "{} SUBROUTINE {}({}){}".format(mods, self.name, params, info)

    @property
    def signature(self):
        """Returns the signature definition for the subroutine."""
        mods = ", ".join(self.modifiers)
        return "{} SUBROUTINE {}({})".format(mods, self.name,
                                             self.parameters_as_string())

    @property
    def end_token(self):
        """Gets the end [code type] token for this instance."""
        return "end subroutine"

    def update(self, name, modifiers):
        """Updates the attributes for the subroutine instance, handles name changes
        in the parent module as well."""
        self.update_name(name)
        self.modifiers = modifiers

class TypeExecutable(CodeElement):
    """Represents a function or subroutine declared in a type that can be executed."""
    
    def __init__(self, name, modifiers, parent, pointsto = None):
        super(TypeExecutable, self).__init__(name, modifiers, parent)
        self.pointsto = pointsto

    def __str__(self):
        mods = ", ".join(self.modifiers)
        pointsto = " => {}".format(self.pointsto) if self.pointsto is not None else ""
        return "{} {}{}".format(mods, self.name, pointsto)

    def parseline(self, line, lineparser):
        """Uses the specified line parser to parse the given line."""
        lineparser.tparser.parseline(self, line)

    def unpickle(self, parent):
        """Sets the parent pointer references for the type executable."""
        self.parent = parent
        self.unpickle_docs()

    @property
    def target(self):
        """Returns the code element that is the actual executable that this type
        executable points to."""
        if self.pointsto is not None:
            #It is in the format of module.executable.
            xinst = self.module.parent.get_executable(self.pointsto.lower())
            return xinst
        else:
            #The executable it points to is the same as its name.
            fullname = "{}.{}".format(self.module.name, self.name)
            return self.module.parent.get_executable(fullname.lower())

class CustomType(CodeElement, Decoratable):
    """Represents a custom defined type in a fortran module."""
    
    def __init__(self, name, modifiers, members, parent):
        super(CustomType, self).__init__(name, modifiers, parent)
        Decoratable.__init__(self)
        #A list of ValueElements() that were declared in the body of the type.
        self.members = members
        #A list of Executable() declared within the contains section of the type.        
        self.executables = {}
        #When an instance is add from just a signature and doesn't have an 
        #end_token, this is set to True
        self.incomplete = False

    def __str__(self):
        execs = "\n\t  - ".join([ x.__str__() for x in self.executables ])
        mods = ", ".join(self.modifiers)
        allexecs = "\n\t  - {}".format(execs) if len(self.executables) > 0 else ""
        mems = "\n\t - ".join([x.__str__() for x in self.members ])
        return "TYPE {} ({}){}\nMEMBERS\n\t{}".format(self.name, mods, allexecs, mems)

    @property
    def fixedvar(self):
        """Returns the name of a member in this type that is non-custom
        so that it would terminate the auto-class variable context chain.
        """
        possible = [m for m in self.members.values() if not m.is_custom]
        #If any of the possible variables is not allocatable or pointer, it will always
        #have a value and we can just use that.
        sufficient = [m for m in possible if "allocatable" not in m.modifiers
                      and "pointer" not in m.modifiers]
        if len(sufficient) > 0:
            return [sufficient[0].name]
        else:
            return [m.name for m in possible]
    
    @property
    def recursive(self):
        """When True, this CustomType has at least one member that is of the same
        type as itself.
        """
        for m in self.members.values():
            if m.kind is not None and m.kind.lower() == self.name.lower():
                return True
        else:
            return False
    
    @property
    def signature(self):
        """Returns the signature definition for the derived type."""
        mods = ", ".join(self.modifiers)
        return "TYPE {} ({})".format(self.name, mods)

    @property
    def end_token(self):
        """Gets the end [code type] token for this instance."""
        return "end type"

    def update_name(self, name):
        """Updates the name of the custom type in this instance and its
        parent reference."""
        if name != self.name:
            self.parent.types[name] = self
            del self.parent.types[self.name]
            self.name = name

    def rt_update(self, statement, linenum, mode, tparser):
        """Uses the specified line parser to parse the given line.

        :arg statement: a string of lines that are part of a single statement.
        :arg linenum: the line number of the first line in the statement relative to
          the entire module contents.
        arg mode: either 'insert', 'replace' or 'delete'
        :arg tparser: an instance of the type parser from the real
          time update module's line parser.
        """
        section = self.find_section(self.module.charindex(linenum, 1))
        if section == "body":
            tparser.parse_line(statement, self, mode)
        elif section == "signature":
            if mode == "insert":
                tparser.parse_signature(statement, self)
        #NOTE: docstrings are handled at a higher level by the line parser
        #because of the XML dependence.

    def unpickle(self, parent):
        """Sets the parent pointer references for the type *and* all of its
        child classes that also have pointer references."""
        self.parent = parent
        self._unpickle_collection(self.members)
        self._unpickle_collection(self.executables)
        self.unpickle_docs()
        
    def _unpickle_collection(self, collection):
        """Unpickles all members of the specified dictionary."""
        for mkey in collection:
            collection[mkey].unpickle(self)

    @property
    def refstring(self):
        """Returns a reference to the string from which this custom type was parsed."""
        if this.parent is not None:
            return this.parent.contents
        else:
            return ""

class Interface(CodeElement):
    """Represetns a fortran assignment, operator or generic interface."""
    def __init__(self, name, modifiers, parent, symbol=None):
        if name in ["operator", "assignment"]:
            if symbol is not None:
                name = "{}|{}".format(name, symbol.lower())
            else:
                raise ValueError("Operator and assignment interfaces require a symbol "
                                 "to identify them in the code.")

        super(Interface, self).__init__(name, modifiers, parent)

        self.procedures = []
        self.symbol = symbol
        self.operator = (name == "operator")
        self.assignment = (name == "assignment")

        #This is a dict of the actual CodeElement instances of the module procedures
        #referenced in the generic interface.
        self._targets = None
        #This is the first module procedure in the list for which we have a valid
        #code element reference.
        self._first = None

    @property
    def parameters(self):
        """Returns the list of parameters from the first valid reference of an
        embedded module procedure in the interface.
        """
        if self.first:
            return self.first.parameters
        else:
            return []

    @property
    def ordered_parameters(self):
        """Returns the list ordered parameters from the first valid reference of an
        embedded module procedure in the interface.
        """
        if self.first:
            return self.first.ordered_parameters
        else:
            return []

    @property
    def subroutine(self):
        """Returns true if the embedded module procedures in this interface are subroutines.
        """
        return isinstance(self.first, Subroutine)

    def changed(self, name):
        """Returns true if the parameter with the specified name has its value changed by
        the *first* module procedure in the interface.

        :arg name: the name of the parameter to check changed status for.
        """
        if self.first:
            return self.first.changed(name)
        else:
            return False

    def get_parameter(self, index):
        """Gets the list of parameters at the specified index in the calling argument
        list for each of the module procedures in the interface.

        :arg index: the 0-based index of the parameter in the argument list.
        """
        result = []
        for target in self.targets:
            if target is not None:
                result.append(target.get_parameter(index))

        return result
        
    def find_docstring(self):
        """Sets the docstring of this interface using the docstrings of its embedded
        module procedures that it is an interface for.
        """
        #The docstrings for the interface are compiled from the separate docstrings
        #of the module procedures it is an interface for. We choose the first of these
        #procedures by convention and copy its docstrings over.
        #Just use the very first embedded method that isn't None and has a docstring.
        if self.first and len(self.docstring) == 0:
            self.docstring = self.first.docstring

    def describe(self):
        """Returns a home-grown description that includes a summary of the calling interface
        for all the module procedures embedded in this interface.
        """
        #Interfaces are tricky because we can have an arbitrary number of embedded procedures
        #that each have different calling interfaces and return types. We are trying to 
        #summarize that information in a single interface. Interfaces can't mix executable
        #types; they are either all subroutines or all functions. 
        self.find_docstring()
        if not self.subroutine:
            #We have a return type to worry about, compile it from all the targets
            ret_type = []
            for target in self.targets:
                if target.returns not in ret_type:
                    ret_type.append(target.returns)

            ret_text = ', '.join(ret_type)
        else:
            ret_text = ""

        #The list of types it accepts is constructed the same way for functions and subroutines.
        act_type = []
        for target in self.targets:
            param = target.ordered_parameters[0]
            if param.strtype not in act_type:
                act_type.append(param.strtype)

        act_text = ', '.join(act_type)

        if self.subroutine:
            return "SUMMARY: {} | ACCEPTS: {}".format(self.summary, act_text)
        else:
            return "SUMMARY: {} | RETURNS: {} | ACCEPTS: {}".format(self.summary, ret_text, act_text)

    def __getstate__(self):
        """Cleans up the object so that it can be pickled without any pointer
        issues. Saves the full names of parents etc. so that they can be 
        restored once the unpickling is completed."""
        odict = self.__dict__.copy() # copy the dict since we change it
        del odict['_targets']
        del odict['_first']
        return odict

    def __setstate__(self, dict):
        self._targets = None
        self._first = None
        self.__dict__.update(dict)

    @property
    def first(self):
        """Returns the first module procedure embedded in the interface that has a valid
        instance of a CodeElement.
        """
        if self._first is None:
            for target in self.targets:
                if target is not None:
                    self._first = target
                    break
            else:
                self._first = False

        return self._first                

    @property
    def targets(self):
        """Provides an ordered list of CodeElement instances for each of the embedded
        module procedures in the generic interface.
        """
        if self._targets is None:
            self._targets = {}
            for key in self.procedures:
                element = self.module.parent.get_executable(key)
                if element is not None:
                    self._targets[key] = element
                else:
                    self._targets[key] = None
                    msg.warn("unable to locate module procedure {}".format(key) +
                             " in interface {}".format(self.name))

        return [self._targets[key] for key in self.procedures]

class Module(CodeElement, Decoratable):
    """Represents a fortran module."""
    
    def __init__(self, name, modifiers, dependencies, publics, contents, parent):
        super(Module, self).__init__(name, modifiers, parent)
        Decoratable.__init__(self)
        #Dependencies is a list of strings in the format module.member that
        #this module requires to operate correctly.
        self.dependencies = dependencies
        #The string contents between the module and end module keywords.
        self.contents = contents
        #A list of methods declared inside the module that were made public
        #using the public keyword in the body of the module (vs. as a 
        #modifier on the member itself).
        self.publics = publics
        self.interfaces = {}
        self.members = {}
        """A list of ValueElements() that were declared in the body of the module."""
        self.types = {}
        """A list of CustomType() declared in the module using fortran type...end type
        """
        self.executables = {}
        """A list of Executable() declared within the contains section of the module.
        """
        self.predocs = {}
        """The dictionary of docstrings extracted from the preamble section of the
        module's contents."""
        self.contains = ""
        """The section in the module after CONTAINS keyword."""
        self.preamble = ""
        """The original string that contains all the members and types before CONTAINS.
        """
        self.refstring = ""
        """The string (file contents) from which the module was parsed."""
        self.filepath = None
        """The path to the library where this module was parsed from."""
        self.change_time = None
        """The datetime that the file was last modified."""
        self.changed = False
        """keeps track of whether the module has had its refstring updated
        via a real-time update since the sequencer last analyzed it."""
        self.precompile = False
        """Does this module require pre-compilation; i.e. does it have 
        pre-processor directives."""
        self.public_linenum = 0
        """The number of the line that contains the first 'public' keyword declaration."""

        #Lines and character counts for finding where matches fit in the file
        self._lines = []
        self._chars = []
        self._contains_index = None
        #Lazy calculation of the adjusted compilation file path.
        self._compile_path = None

    def search_dependencies(self):
        """Returns a list of other modules that this module and its members and executables
        depend on in order to run correctly. This differs from the self.dependencies attribute
        that just hase explicit module names from the 'use' clauses.
        """
        result = [self.module.name]
        #First we look at the explicit use references from this module and all
        #its dependencies until the chain terminates.
        stack = self.needs
        while len(stack) > 0:
            module = stack.pop()
            if module in result:
                continue
            
            self.parent.load_dependency(module, True, True, False)
            if module in self.parent.modules:
                for dep in self.parent.modules[module].needs:
                    modname = dep.split(".")[0]
                    if modname not in result:
                        result.append(modname)
                    if modname not in stack:
                        stack.append(modname)

        #Add any incidentals from the automatic construction of code. These can be from
        #executables that use special precision types without importing them or from
        #derived types. Same applies to the local members of this module.
        for ekey, anexec in list(self.executables.items()):
            for dep in anexec.search_dependencies():
                if dep is not None and dep not in result:
                    result.append(dep)

        for member in list(self.members.values()):
            dep = member.dependency()
            if dep is not None and dep not in result:
                result.append(dep)
                        
        return result

    @property
    def codefolder(self):
        """Returns the full path to the code folder (library) that this module is in.
        """
        from os import path
        return path.dirname(self.filepath)
    
    @property
    def compile_path(self):
        """Returns the file path to this module taking the pre-processing into account.
        If the module requires pre-processing, the extension is reported as F90; otherwise,
        the regular self.filepath is returned.
        """
        if not self.precompile:
            return self.filepath
        else:
            if self._compile_path is None:
                segs = self.filepath.split('.')
                segs.pop()
                segs.append("F90")
                self._compile_path = '.'.join(segs)

            return self._compile_path

    @property
    def xmlpath(self):
        """Returns the full path to this module's accompanying XML file."""
        if self.filepath is None:
            return None
        else:
            return self.parent.get_xmldoc_path(self.filepath)

    def all_to_public(self):
        """Sets all members, types and executables in this module as public as
        long as it doesn't already have the 'private' modifier.
        """
        if "private" not in self.modifiers:
            def public_collection(attribute):
                for key in self.collection(attribute):
                    if key not in self.publics:
                        self.publics[key.lower()] = 1
            public_collection("members")
            public_collection("types")
            public_collection("executables")
        
    def rt_update(self, statement, linenum, mode, modulep, lineparser):
        """Uses the specified line parser to parse the given statement.

        :arg statement: a string of lines that are part of a single statement.
        :arg linenum: the line number of the first line in the statement relative to
          the entire module contents.
        :arg mode: either 'insert', 'replace' or 'delete'
        :arg modulep: an instance of ModuleParser for handling the changes
          to the instance of the module.
        :arg lineparser: a line parser instance for dealing with new instances
          of types or executables that are add using only a signature.
        """
        #Most of the module is body, since that includes everything inside
        #of the module ... end module keywords. Since docstrings are handled
        #at a higher level by line parser, we only need to deal with the body
        #Changes of module name are so rare that we aren't going to bother with them.
        section = self.find_section(self.module.charindex(linenum, 1))
        if section == "body":
            modulep.rt_update(statement, self, mode, linenum, lineparser)

        #NOTE: docstrings are handled at a higher level by the line parser
        #because of the XML dependence.

    def unpickle(self, parent):
        """Sets the parent pointer references for the module *and* all of its
        child classes that also have pointer references."""
        self.parent = parent
        self._unpickle_collection(self.members)
        self._unpickle_collection(self.executables)
        self._unpickle_collection(self.types)
        
    def _unpickle_collection(self, collection):
        """Unpickles all members of the specified dictionary."""
        for mkey in collection:
            collection[mkey].unpickle(self)

    def __str__(self):
        output = []
        #Run statistics on the lines so it displays properly
        self.linenum(1)
        preprocess = ", preprocess" if self.precompile else ""
        output.append("MODULE {} ({} lines{})\n\n".format(self.name, len(self._lines), preprocess))
        uses = "\n\t".join(self.sorted_collection("dependencies"))
        output.append("USES:\n\t{}\n\n".format(uses))

        if len(self.types) > 0:
            types = "\n\t".join([ t[1].__str__() for t in list(self.types.items()) ])
            output.append("TYPES:\n\t{}\n\n".format(types))

        if len(self.executables) > 0:
            functions = "\n\t".join([ x[1].__str__() for x in list(self.functions().items()) ])
            subroutines = "\n\t".join([ x[1].__str__() for x in list(self.subroutines().items()) ])
            output.append("EXECUTABLES:\n\t{}\n\n\t{}\n\n".format(functions, subroutines))

        if len(self.members) > 0:
            members = "\n\t".join([ x[1].__str__() for x in list(self.members.items()) ])
            output.append("MEMBERS:\n\t{}\n\n".format(members))

        return "".join(output)

    def set_public_start(self, start):
        """Sets the line number of the first declaration of public members for the module.
        
        :arg start: the *character* number of the first member marked as public.
        """
        if start > 0:
            self.public_linenum = self.linenum(start)

    def get_dependency_element(self, symbol):
        """Checks if the specified symbol is the name of one of the methods
        that this module depends on. If it is, search for the actual code
        element and return it."""
        for depend in self.dependencies:
            if "." in depend:
                #We know the module name and the executable name, easy
                if depend.split(".")[1] == symbol.lower():
                    found = self.parent.get_executable(depend)
                    break
            else:
                #We only have the name of an executable or interface, we have to search
                #the whole module for the element.
                fullname = "{}.{}".format(depend, symbol)
                if self.parent.get_executable(fullname) is not None:
                    found = self.parent.get_executable(fullname)
                    break
                if self.parent.get_interface(fullname) is not None:
                    found = self.parent.get_interface(fullname)
                    break
        else:
            return None

        return found

    def completions(self, symbol, attribute, recursive = False):
        """Finds all possible symbol completions of the given symbol that belong
        to this module and its dependencies.

        :arg symbol: the code symbol that needs to be completed.
        :arg attribute: one of ['dependencies', 'publics', 'members', 
          'types', 'executables'] for specifying which collections to search.
        :arg result: the possible completions collected so far in the search.
        """
        possible = []
        for ekey in self.collection(attribute):
            if symbol in ekey:
                possible.append(ekey)

        #Try this out on all the dependencies as well to find all the possible
        #completions.
        if recursive:
            for depkey in self.dependencies:
                #Completions need to be fast. If the module for the parent isn't already
                #loaded, we just ignore the completions it might have.
                if depkey in self.parent.modules:
                    possible.extend(self.parent.modules[depkey].completions(symbol, attribute))
            
        return possible

    @property
    def end_token(self):
        """Gets the end [code type] token for this instance."""
        return "end module"

    @property
    def contains_index(self):
        """Returns the *line number* that has the CONTAINS keyword separating the
        member and type definitions from the subroutines and functions."""
        if self._contains_index is None:
            max_t = 0
            for tkey in self.types:
                if self.types[tkey].end > max_t and not self.types[tkey].embedded:
                    max_t = self.types[tkey].end

            #Now we have a good first guess. Continue to iterate the next few lines
            #of the the refstring until we find a solid "CONTAINS" keyword. If there
            #are no types in the module, then max_t will be zero and we don't have
            #the danger of running into a contains keyword as part of a type. In that
            #case we can just keep going until we find it.
            i = 0
            start = self.linenum(max_t)[0]
            max_i = 10 if max_t > 0 else len(self._lines)

            while (self._contains_index is None and i < max_i
                   and start + i < len(self._lines)):
                if "contains" in self._lines[start + i].lower():
                    self._contains_index = start + i
                i += 1

            if self._contains_index is None:
                #There must not be a CONTAINS keyword in the module
                self._contains_index = len(self._lines)-1

        return self._contains_index

    @property
    def needs(self):
        """Returns a unique list of module names that this module depends on."""
        result = []
        for dep in self.dependencies:
            module = dep.split(".")[0].lower()
            if module not in result:
                result.append(module)

        return result

    def absolute_charindex(self, string, start, end):
        """Gets the absolute start and end indices for a regex match
        with reference to the original module file."""
        search = string[start:end]
        abs_start = self.refstring.index(search)
        return abs_start, (end - start) + abs_start

    def functions(self):
        """Returns a dictionary of all the functions in the module."""
        return self._filter_execs(False)
        
    def subroutines(self):        
        """Returns a dictionary of all the subroutines in the module."""
        return self._filter_execs(True)

    def _filter_execs(self, isSubroutine):
        """Filters the executables in the dictionary by their type."""
        result = {}
        for key in self.executables:
            if (isinstance(self.executables[key], Subroutine) and isSubroutine) or \
               (isinstance(self.executables[key], Function) and not isSubroutine):
                result[key] = self.executables[key]
        
        return result                

    def warn(self, collection):
        """Checks the module for documentation and best-practice warnings."""
        super(CodeElement, self).warn(collection)
        if not "implicit none" in self.modifiers:
            collection.append("WARNING: implicit none not set in {}".format(self.name))

    def type_search(self, basetype, symbolstr):
        """Recursively traverses the module trees looking for the final
        code element in a sequence of %-separated symbols.

        :arg basetype: the type name of the first element in the symbol string.
        :arg symblstr: a %-separated list of symbols, e.g. this%sym%sym2%go.
        """
        return self.parent.type_search(basetype, symbolstr, self)

    def sorted_collection(self, attribute):
        """Returns the names of all elements in a collection sorted."""
        return sorted(self.collection(attribute))
    
    def collection(self, attribute):
        """Returns the collection corresponding the attribute name."""
        return {
            "dependencies": self.dependencies,
            "publics": self.publics,
            "members": self.members,
            "types": self.types,
            "executables": self.executables,
            "interfaces": self.interfaces
        }[attribute]

    def update_refstring(self, string):
        """Updates the refstring that represents the original string that
        the module was parsed from. Also updates any properties or attributes
        that are derived from the refstring."""
        self.refstring = string
        self._lines = []
        self._contains_index = None
        self.changed = True

        #The only other references that become out of date are the contains
        #and preamble attributes which are determined by the parsers.
        #Assuming we did everything right with the rt update, we should be
        #able to just use the new contains index to update those.
        icontains = self.contains_index
        ichar = self.charindex(icontains, 0)
        self.preamble = string[:ichar]
        self.contains = string[ichar + 9:]

    def update_elements(self, line, column, charcount, docdelta=0):
        """Updates all the element instances that are children of this module
        to have new start and end charindex values based on an operation that
        was performed on the module source code.

        :arg line: the line number of the *start* of the operation.
        :arg column: the column number of the start of the operation.
        :arg charcount: the total character length change from the operation.
        :arg docdelta: the character length of changes made to types/execs
          that are children of the module whose docstrings external to their
          definitions were changed.
        """
        target = self.charindex(line, column) + charcount

        #We are looking for all the instances whose *start* attribute lies
        #after this target. Then we update them all by that amount.
        #However, we need to be careful because of docdelta. If an elements
        #docstring contains the target, we don't want to update it.
        if line < self.contains_index:
            for t in self.types:
                if self._update_char_check(self.types[t], target, docdelta):
                    self._element_charfix(self.types[t], charcount)

            for m in self.members:
                if self.members[m].start > target:
                    self.members[m].start += charcount
                    self.members[m].end += charcount

            self._contains_index = None
        else:
            for iexec in self.executables:
                if self._update_char_check(self.executables[iexec], target, docdelta):
                    self._element_charfix(self.executables[iexec], charcount)

    def _update_char_check(self, element, target, docdelta):
        """Checks whether the specified element should have its character indices
        updated as part of real-time updating."""
        if docdelta != 0:
            if (element.docstart <= target and 
                element.docend >= target - docdelta):
                return True
            else:
                return element.absstart > target
        else:
            return element.absstart > target

    def _element_charfix(self, element, charcount):
        """Updates the start and end attributes by charcount for the element."""
        element.start += charcount
        element.docstart += charcount
        element.end += charcount
        element.docend += charcount

    def get_element(self, line, column):
        """Gets the instance of the element who owns the specified line
        and column."""
        ichar = self.charindex(line, column)
        icontains = self.contains_index
        result = None

        if line < icontains:
            #We only need to search through the types and members.
            maxstart = 0
            tempresult = None
            for t in self.types:
                if ichar >= self.types[t].absstart:
                    if self.types[t].absstart > maxstart:
                        maxstart = self.types[t].absstart
                        tempresult = self.types[t]

            #This takes the possibility of an incomplete type into account
            if (tempresult is not None and (ichar <= tempresult.end or 
                                            tempresult.incomplete)):
                result = tempresult

            if not result:
                #Members only span a single line usually and don't get added
                #without an end token.
                for m in self.members:
                    if (ichar >= self.members[m].start and 
                        ichar <= self.members[m].end):
                        result = self.members[m]
                        break
        else:
            #We only need to search through the executables
            tempresult = None
            maxstart = 0
            for iexec in self.executables:
                if (ichar >= self.executables[iexec].absstart):
                    if self.executables[iexec].absstart > maxstart:
                        maxstart = self.executables[iexec].absstart
                        tempresult = self.executables[iexec]

            if tempresult is not None and (ichar <= tempresult.end or
                                           tempresult.incomplete):
                result = tempresult

        if result is None:
            #We must be in the text of the module, return the module
            return self
        else:
            return result

    def update_embedded(self, attribute):
        """Updates the elements in the module 'result' that have character indices
        that are a subset of another element. These correspond to cases where the
        type or subroutine is declared within another subroutine.

        :attr attribute: the name of the collection to update over.
        """
        #The parser doesn't handle embeddings deeper than two levels.
        coll = self.collection(attribute)
        keys = list(coll.keys())
        for key in keys:
            element = coll[key]
            new_parent = self.find_embedded_parent(element)
            if new_parent is not None:
                #Update the parent of the embedded element, add it to the collection
                #of the parent element and then delete it from the module's collection.
                element.parent = new_parent
                if attribute == "types":
                    new_parent.types[key] = element
                else:
                    new_parent.executables[key] = element
                del coll[key]

    def find_embedded_parent(self, element):
        """Finds the parent (if any) of the embedded element by seeing
        if the start and end indices of the element are a subset of 
        any type or executable in this module.
        """
        result = None
        for t in self.types:
            if (element.start > self.types[t].start and
                element.end < self.types[t].end):
                result = self.types[t]
                break
        else:
            for iexec in self.executables:
                if (element.start > self.executables[iexec].start and
                    element.end < self.executables[iexec].end):
                    result = self.executables[iexec]
                    break

        return result

    def charindex(self, line, column):
        """Gets the absolute character index of the line and column
        in the continuous string."""
        #Make sure that we have chars and lines to work from if this
        #gets called before linenum() does.
        if len(self._lines) == 0:
            self.linenum(1)

        if line < len(self._chars):
            return self._chars[line - 1] + column
        else:
            return len(self.refstring)

    def linenum(self, index):
        """Gets the line number of the character at the specified index.
        If the result is unknown, -1 is returned."""
        if len(self._lines) == 0 and self.refstring != "":
            self._lines = self.refstring.split("\n")
            #Add one for the \n that we split on for each line
            self._chars = [ len(x) + 1 for x in self._lines ]
            #Remove the last line break since it doesn't exist
            self._chars[-1] -= 1

            #Now we want to add up the number of characters in each line
            #as the lines progress so that it is easy to search for the
            #line of a single character index
            total = 0
            for i in range(len(self._chars)):
                total += self._chars[i]
                self._chars[i] = total

        if len(self._lines) > 0:
            #Now, find the first index where the character value is >= index
            result = -1
            i = 0
            while result == -1 and i < len(self._chars):
                if index <= self._chars[i]:
                    result = [ i, self._chars[i] - index]
                i += 1

            return result
        else:
            return [ -1, -1 ]
