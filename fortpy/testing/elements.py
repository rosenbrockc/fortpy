from .. import msg
from fortpy.testing.templates import FileLine
from fortpy.docelements import DocElement
from os import path, remove
from shutil import copyfile
import re

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
           (self.attributes["kind"].lower() == element.attributes["kind"].lower()):
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
        """Returns the fortran declaration that defines a global variable."""
        if self.ignore:
            return "  ! Variable '{}' was set to be ignored.".format(self.attributes["name"])
        
        result = []
        if "type" not in self.attributes or "name" not in self.attributes:
            msg.err("required variable for execution missing some attributes." + 
                    " {}".format(self.attributes))
            exit(1)

        result.append(self.attributes["type"])
        if "kind" in self.attributes and self.attributes["kind"] is not None:
            result.append("({})".format(self.attributes["kind"]))

        if "modifiers" in self.attributes:
            mods = re.split(",\s*", self.attributes["modifiers"])
            smods = self._clean_mods(mods)
        else:
            smods = ""

        if smods != "":
            result.append(", " + smods)
        result.append(" :: ")
        result.append(self.attributes["name"])

        if "dimensions" in self.attributes:
            result.append("({})".format(self.attributes["dimensions"]))

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
        else:
            return None

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
 
    def copy(self, coderoot, testroot, case):
        """Copies the input files needed for this value to set a variable.

        :arg coderoot: the full path to folder with the code files.
        :arg testroot: the full path the folder where the test is running.
        :arg case: the case id for multi-case testing.
        """
        #We only want to do the copy if we are assigning a value from a file.
        if self.folder is not None:
            source = path.join(coderoot, self.folder[2:], self.filename.format(case))
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
                        copyfile(filename, target)
            else:
                if self.rename is not None:
                    target = path.join(testroot, self.rename)
                else:
                    target = path.join(testroot, self.filename.format(case))
                copyfile(source, target)

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
            self._derived_type, typemod = self.parent.parser.tree_find(target, module, "types")

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
 
    def code(self, lines, position, spacer, slices=None):
        """Appends the code lines to initialize the parent variable.

        :arg lines: the list of strings to append to.
        :arg position: one of ['vars', 'init', 'assign', 'before', 'after'] indicating 
          where in the fortran program the code will be appended.
        :arg slices: for array value assignments, the specific indices to assign values
          to. Tuple of (slice string, [loopvars]).
        """
        if self.parent.variable is not None:
            self._codes[position](lines, spacer, slices)
        else:
            raise ValueError("Trying to assign a value to an unknown variable {}.".format(self.iid))

    def _code_before(self, lines, spacer, slices):
        """Calls _code_setvar() if we *are* in repeat mode."""
        if self.repeats:
            self._code_setvar(lines, spacer, slices)

    def _code_assign(self, lines, spacer, slices):
        """Calls _code_setvar() if we are *not* in repeat mode."""
        if not self.repeats:
            self._code_setvar(lines, spacer, slices)

    def _code_setvar(self, lines, spacer, slices):
        """Appends code for assigning the value of the parent variable using
        this value specification.

        :arg slices: for array value assignments, the specific indices to assign values
          to.
        """
        if self.filename is None:
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
            self._code_file(lines, spacer, slices)

    def _code_setvar_value(self, lines, spacer, value, slices=None):
        """Sets the value of the variable primitively to the specified value."""
        if slices is None:
            lines.append("{}{} = {}".format(spacer, self.parent.name, value))
        else:
            lines.append("{}{}({}) = {}".format(spacer, self.parent.name, 
                                                slices[0], value))
    
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

                #Handle the embedded via file possibility.
                if varname is not None:
                    params = params.replace("@file", varname)

                lines.append("{}{}{}({})%{}({})".format(spacer, call, self.parent.name, 
                                                        slices[0], self.embedded, params))
        #if it does have pre-reqs, they will be handled by the method writer and we
        #don't have to worry about it.

    def _code_file(self, lines, spacer, slices):
        """Appends code for initializing the value of a variable that is a
        scalar, vector or 2D matrix from a file.
        """
        if self.filename is not None:
            flines = []
            #If we have slices being set from file values, we will have a set
            #of files that match a wildcard pattern. They will all have been copied
            #into the test directory. However, we need to replace the wildcard in
            #the filename at *runtime* with the current value of the loop variable.
            if "*" in self.xname:
                indices = "(/ {} /)".format(', '.join(slices[1]))
                flines.append("call fpy_period_join_indices(" +
                              "{}_pslist, {}, {})".format(self.iid, indices, len(slices[1])))
                rtlen = "len({}_pslist, 1) + {}".format(self.iid, len(self.xname)-1)
                rtname = '"{}"//{}_pslist'.format(self.xname[:len(self.xname)-1], self.iid)
            else:
                rtname = "'{}'".format(self.xname)
                rtlen = len(self.xname)

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

            if slices is not None:
                varname = "{}_fvar".format(self.iid)
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
            if (self.D in [0, 1, 2]) and not self.ragged:
                fmtstr = "call fpy_read{}({}, '{}', {})"
                if ("pointer" in self.parent.global_attr("modifiers", []) and
                    "fvar" not in varname):
                    flines.append(fmtstr.format("_p", rtname, self.commentchar, varname))
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

    def _code_after(self, lines, spacer, slices=None):
        """Appends code for deallocating a variable that was assigned a value.
        This is useful so that we can reallocate it again in repeat mode.
        """
        if (self.repeats and self.parent.allocate and 
            ("allocatable" in self.parent.variable.modifiers or
             "pointer" in self.parent.variable.modifiers)):
            lines.append("{}deallocate({})".format(spacer, self.parent.name))

    def _code_init(self, lines, spacer, slices=None):
        """Appends code to initialize any variables we need for assignment
        operations later on.
        """

    def _code_vars(self, lines, spacer, slices=None):
        """Appends lines to declare any variables we need for file read
        operations later on.
        """
        if self.filename is not None:
            if "*" in self.xname or self.ragged:
                if self.dtype is None:
                    raise ValueError("Wildcard file names and ragged array inputs both require "
                                     "attribute 'dtype' to be specified.")
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
            raise ValueError("'identifier' is a required attribute of <value> tags.")

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
            self.value = self.xml.attrib["value"]
        else:
            raise ValueError("'value' is a required attribute of <if>, <elseif>"
                             " and <else> tags.")

    def code(self, lines, spacer):
        """Appends the code for this condition and its variable assignment."""
        if self.tag in ["if", "elseif"]:
            lines.append("{}{} ({}) then".format(spacer, self.tag, self.condition))
        else:
            lines.append("else")

        #Append the value assignment. To do this we have to look it up in the
        #grand-parents list of possible value assignments.
        if self.value in self.parent.values:
            valobj = self.parent.values[self.value]
            valobj.code(lines, "assign", spacer + "  ")
        else:
            raise ValueError("Could not find value '{}' for condition.".format(self.value))
        
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
                if c.value in self.values and self.values[c.value].repeats:
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
        """Adds the code to make this part specification (and loop) work."""
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
        if self.value is not None and self.value.lower() not in self.assignment.values:
            raise ValueError("{} in part {} is not ".format(self.value, self.identifier) + 
                             "a valid identifier.")
            
        #Before we attempt to assign values, we need to make sure that the variable is allocated
        #For non-file allocations, the Assignment instance handles allocation; for file types,
        #usually the AssignmentValue instance handles it, except when <part> tags are used.
        if self.assignment.values[self.value].filename is not None:
            if self.assignment.allocate:
                if self.assignment.allocate == True:
                    lines.append("{}allocate({})".format(spacer, self.assignment.name))
                else:
                    lines.append("{}allocate({}({}))".format(spacer, self.assignment.name,
                                                             self.assignment.allocate))            

        if self.isloop:
            lines.append("{}do {}={}, {}".format(spacer, self.loopid, self.limits[0], self.limits[1]))
            self.assignment.values[self.value].code(lines, position, spacer + '  ', 
                                                    (self.slices, self.slice_list))
        else:
            #We need to repeat the assignment for each of the part values
            for ipart in self.limits:
                slices = self.slices.format(ipart)
                self.assignment.values[self.value].code(lines, position, spacer + '  ', 
                                                        (slices, self.slice_list))

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
            self.value = self.xml.attrib["value"]
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
                self.limits = map(int, self.xml.attrib["limits"].split())
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
                self._variable = self.globaldecl

        return self._variable

    # def _recurse_find_variable(name, finder):
    #     """Finds the ValueElement instance for the variable of the specified name
    #     if it exists. Found by searching the parameter list of all methods that get
    #     called as part of the pre-req chain."""
    #     result = None

    #     if isinstance(finder, MethodFinder) and finder.executable is not None:
    #         if name in finder.executable.parameters:
    #             result = finder.executable.parameters[name]
    #         else:
    #             for nested in finder.methods:
    #                 result = MethodFinder.recurse_find_variable(name, nested)
    #                 if result is not None:
    #                     break
    #     return result                        
    
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
                 (self.variable.D > 0 and ":" in self.variable.dimension)))
    
    def check_prereqs(self, finder):
        """Checks all the values this assignment may use to make sure they reference
        a valid derived type. If an AssignmentValue specifies to check pre-reqs, 
        they are added to the pre-req chain.

        :arg finder: the instance of MethodFinder for which the pre-reqs are being
          checked and added.
        """
        for v in self.values:
            self.values[v].check_prereqs(finder)

    def copy(self, coderoot, testroot, case):
        """Copies the input files needed for this assignment to set a variable.

        :arg coderoot: the full path to folder with the code files.
        :arg testroot: the full path the folder where the test is running.
        :arg case: the case id for multi-case testing.
        """
        #Just call copy on all the child value objects.
        for v in self.values:
            self.values[v].copy(coderoot, testroot, case)        
        
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
            self.allocatable and variable.D == 1):
            lines.append("{}allocate({}({}))".format(spacer, self.name, self.allocate))
                    
        if (variable.dimension is not None and
            variable.D == 2 and self.allocatable):                  
            allocstr = "{2}allocate({0}({1}))"
            lines.append(allocstr.format(self.name, self.allocate, spacer))

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
            for v in self.values:
                self.values[v].code(lines, position, spacer)
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
                if self.value in self.values:
                    self.values[self.value].code(lines, position, spacer)
                elif self.value in self.parts:
                    self.parts[self.value].code(lines, position, spacer)
                else:
                    raise ValueError("Value identifier {} not found.".format(self.value))
        
    def _parse_xml(self):
        """Parses attributes and child tags from the XML element."""
        if "name" in self.attributes:
            self.name = self.attributes["name"]
        else:
            raise ValueError("'name' is a required attribute for <assignment> tags.")

        if "value" in self.attributes:
            self.value = self.attributes["value"]
        if "constant" in self.attributes:
            #We want to manually create an AssignmentValue object for the constant and
            #add it to the child list.
            val = AssignmentValue(None, self)
            val.identifier = "constant"
            val.constant = self.attributes["constant"]
            self.values["constant"] = val
            self.value = "constant"
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
                self.values[val.identifier] = val
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
            self.writekey += "_{}".format(self.value)
         
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
    def __init__(self, xmltag):
        self.xml = xmltag
        self.folder = None
        self.filename = None
        self.rename = None
        self.line = None
        self.position = "before"

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

    def copy(self, coderoot, testroot, case=""):
        """Copies the input file from the specified code root directory to
        the folder where the test is being performed.

        :arg coderoot: the full path to the folder that houses all the code files.
        :arg testroot: the full path to the folder that the parent test is being
          performed in.
        :arg case: if a specific case of the same test is being performed, the
          case identifier to use for string formatting.
        """
        source = path.join(coderoot, self.folder[2:], self.filename.format(case))
        if self.rename is not None:
            target = path.join(testroot, self.rename)
        else:
            target = path.join(testroot, self.filename.format(case))
        copyfile(source, target)

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
        
        self._parse_lines()

    def _parse_lines(self):
        """Checks for the existence of a <line> tag in the input specification.
        If it exists, the line attribute is initialized with a FileLine object.
        """
        for child in self.xml:
            if child.tag == "line":
                self.line = FileLine(child)
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

        self._parse_attributes()

    @property
    def filemode(self):
        """Returns a value indicating whether the output comparison is with
        a file (as opposed to a hard-coded value)."""
        return (self.value == None)

    def abspath(self, coderoot, caseid=""):
        """Returns the absolute path to the output file referenced by this
        output comparison specifier.

        :arg coderoot: the full path to the directory storing the code files.
        :arg case: if a specific case of the same test is being performed, the
          case identifier to use for string formatting.
        """
        #The caseid is the "testid.case". We only want the 'case' part of it
        if caseid != "":
            case = caseid.split(".")[-1]
            return path.join(coderoot, self.folder[2:], self.filename.format(case))
        else:
            return path.join(coderoot, self.folder[2:], self.filename)

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
            self.tolerance = 1
        if "position" in self.xml.attrib:
            self.position = self.xml.attrib["position"]
        
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

        self._code = None
        self._parse_attributes()

    def code(self, lines, variables, position, spacer):
        """Appends the code required to save the value of this target to
        file if the specified position is appropriate.
        """
        #If the target is trying to get a file to be compared to another
        #model output file, then it must start with a dot. In that case
        #we don't generate any code since we assume that the file is being
        #created by the method being unit tested.
        if position in self.when and self.name[0] != ".":
            if self._code is None:
                self._process(variables, spacer)
                
            lines.append(self._code)

    def clean(self, testroot):
        """Removes the output file for saving the variable's value from the
        testing folder."""
        if self.varfile is not None:
            target = path.join(testroot, self.varfile)
            if path.isfile(target):
                remove(target)

    def _process(self, variables, spacer):
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
            if self.name in variables:
                glob = variables[self.name]
                dtype = glob.attributes["type"]
                if dtype == "class" or dtype == "type":
                    if self.member is not None:
                        self._code = ("{}call pysave({}%{}".format(spacer, self.name, self.member) + 
                                      ", '{}', {})".format(self.varfile, len(self.varfile)))
                    else:
                        self._code = "{}call {}%test_output('{}')".format(spacer, self.name, self.varfile)
                else:
                    #pysave is an interface in the fortpy module that can accept variables
                    #of different types and write them to file in a deterministic way
                    self._code = "{}call pysave({}, '{}', {})".format(spacer, self.name, self.varfile,
                                                                      len(self.varfile))
            elif self.name == "[default]":
                #We need to save the value that was generated by the main *function* being
                #unit tested. The variable defined in the fortpy program is the function
                #name with a suffix of "_fpy".
                self._code = "{}call pysave({}_fpy, '{}', {})".format(spacer, self.testspec.executable.name,
                                                                      self.varfile, len(self.varfile))
            else:
                raise ValueError("Variable {} has not been initialized with ".format(self.name) +
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
            self.varfile = "{}.fortpy".format(self.name.replace("%", "."))
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

        self.targets = []
        self.inputs = []
        self.outputs = {}
        self.methods = []
        """A set of *local* assignments and pre-reqs that apply to this test."""
        self.variables = {}
        """A set of *local* variable declarations for this test."""
        self._variable_order = []

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
            raise ValueError("'identifier' is a required attribute for <test> tags.")

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
        def group_list(collection, position, target):
            """Appends all the elements in the collection that match the position
            specification to the given target.
            """
            for item in collection:
                if item.position == position:
                    target.append(item)

        def group_dict(collection, position, target, ordering=None, tordering=None):
            """Updates the 'target' dictionary to include all items in the 'collection'
            that match the specified position.
            """
            if ordering is not None:
                keys = ordering
            else:
                keys = collection
                
            for item in keys:
                if collection[item].position == position:
                    if tordering is not None:
                        tordering.append(item)
                    target[item] = collection[item]
            
        #We need to add all the pre-reqs and assignments from the parent group to
        #this one since they were specified globally for *all* tests.
        group_list(self.testgroup.methods, "before", self.methods)
        group_list(self.testgroup.inputs, "before", self.inputs)
        group_list(self.testgroup.targets, "before", self.targets)
        group_dict(self.testgroup.outputs, "before", self.outputs)
        group_dict(self.testgroup.variables, "before", self.variables,
                   self.testgroup.variable_order, self._variable_order)
                        
        self._parse_attributes()
        self._parse_children()

        #Now handle the methods from the group that should be added *after* the local
        #assignments and pre-reqs are done.
        group_list(self.testgroup.methods, "after", self.methods)
        group_list(self.testgroup.inputs, "after", self.inputs)
        group_list(self.testgroup.targets, "after", self.targets)
        group_dict(self.testgroup.outputs, "after", self.outputs)
        group_dict(self.testgroup.variables, "after", self.variables,
                   self.testgroup.variable_order, self._variable_order)
            
    def _parse_children(self):
        """Gets all child tags from the <test> and parses them for relevant
        content to set targets, inputs and output attributes."""
        for child in self.xml:
            if child.tag == "input":
                self.inputs.append(TestInput(child))
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
        if "cases" in self.xml.attrib:
            # For specifying the cases in a unit test, we should allow ranges like
            # "standard.cr[1-12]" so that the developer doesn't need to enter each
            # of the cases separately. We should still allow a comma-separated list
            # of cases, but each must allow the shorthand notation.
            rawcases = re.split(",\s*", self.xml.attrib["cases"])
            if "[" in self.xml.attrib["cases"]:
                self.cases = []
                rxcase = r"(?P<prefix>[^[]*)\[(?P<range>\d+-\d+)](?P<suffix>.*)"
                recase = re.compile(rxcase)
                for craw in rawcases:
                    m = recase.match(craw)
                    if m:
                        prefix = m.group("prefix")
                        vals = map(int, m.group("range").split("-"))
                        suffix = m.group("suffix")
                        if prefix is None:
                            prefix = ""
                        if suffix is None:
                            suffix = ""
                        for v in range(vals[0], vals[1]+1):
                            self.cases.append("{}{}{}".format(prefix, v, suffix))
                    else:
                        self.cases.append(craw)
            else:
                self.cases = rawcases
        if "runtime" in self.xml.attrib:
            runtime = list(self.xml.attrib["runtime"])
            self.unit = runtime.pop()
            self.runtime = [ int(t) for t in "".join(runtime).split("-") ]
        if "timed" in self.xml.attrib:
            self.timed = self.xml.attrib["timed"] == "true"

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
        
        self.finder = None
        """The currently active MethodFinder instance being used by the MethodWriter
        to create the executable driver for a specific testid."""
        self.variable_order = []
        """The order in which the global declarations appear in the group."""
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
        return self.group.name

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
                    docel.group.name == self.group.name):
                    self.children.append(docel)

    def _parse_xml(self):
        """Parses the XML structure into the relevant dictionaries and
        properties for the testing group.
        """
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
                self.inputs.append(TestInput(child.xml))
            elif child.doctype == "output":
                outvar = TestOutput(child.xml)
                self.outputs[outvar.identifier] = outvar
            elif child.doctype == "target":
                self.targets.append(TestTarget(child.xml, self))

        self._parse_mappings()
        self._parse_variables()

        #Now that the testing group has all its children parsed, we can
        #let the tests parse. This is because some of the global tags in the
        #testing group will apply to every test.
        for test in self.tests:
            self.tests[test].parse()

    def _parse_variables(self):
        """Searches for <global> tags and adds them to the variables dict."""
        for child in self.children:
            if (child.doctype == "global" and "name" in child.attributes):
                TestingGroup.global_add(self.variables, self.variable_order,
                                       child.attributes["name"].lower(), child)

        if type(self.element).__name__ in ["Subroutine", "Function"]:
            self._get_param_globals()

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
                       "regular" in doc.attributes and doc.attributes["regular"] == "true":
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
            self._global_clean_param(result, "kind", param.kind)
            self._global_clean_param(result, "modifiers", ", ".join(param.modifiers))
            self._global_clean_param(result, "dimensions", param.dimension)
            self._global_clean_param(result, "default", param.default)
            
            #if the variable is a deferred shape array and we have been told to allocate
            #it during the initializing phase then we must add "allocatable" as a modifier
            #if it isn't already allocatable or a pointer.
            if name in self.assignments and param.D > 0:
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
        if self.group is not None:
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
            lines.append("{}call pysave(fpy_elapsed, 'fpytiming.out', 13)".format(spacer))

    def _code_vars(self, lines, spacer, testid):
        """Appends code to initialize the variables that accept the return values
        if the method being called is a function. Also appends the variables needed
        to time the method's execution if specified by the test.
        """
        #This method only gets called if we are the main executable being tested.
        if type(self.executable).__name__ == "Function":
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
            prefix = "{}_fpy = ".format(self.executable.name)

        spacing = len(list(prefix)) + len(list(self.executable.name)) + len(spacer)

        if "paramlist" in self.attributes:
            #The developer has decided on an alternate call signature from
            #the one that gets auto-generated.
            specified = re.split(",\s*", self.attributes["paramlist"])
            cleaned = self._present_params(specified, spacing)
            lines.append("{}{}{}({})".format(spacer, prefix, self.executable.name, cleaned))
        else:
            #We can construct the actual list of parameters to use in the call
            #and substitute mappings where appropriate.
            calllist = []
            for param in self.executable.ordered_parameters:
                #Because of ignorable parameters, if the parameter is optional, explicitly
                #specify its name.
                if "optional" in param.modifiers:
                    optstr = "{}=".format(param.name)
                else:
                    optstr = ""
                    
                if self.group is not None and param in self.group.mappings:
                    calllist.append(optstr + self.group.mappings[param])
                else:
                    var = None
                    if self.group is not None and param.name in self.group.variables:
                        var = self.group.variables[param.name]
                    if var is None and self.test is not None and param.name in self.test.variables:
                        var = self.test.variables[param.name]
                    if var is not None and not var.ignore:
                        calllist.append(optstr + param.name)
                    elif var is None:
                        calllist.append(optstr + param.name)

            lines.append("{}{}{}({})".format(spacer, prefix, self.executable.name,
                                              self._present_params(calllist, spacing)))        

        if self.timed(testid):
            lines.append("{}call cpu_time(fpy_end)".format(spacer))
            lines.append("{}fpy_elapsed = fpy_elapsed + fpy_end - fpy_start".format(spacer))

    def _present_params(self, paramlist, spacing = 0):
        """Creates the (paramlist) for a method call formatted nicely for calls
        with lots of parameters."""
        #The +2 is spacing is for the tab indent at the start of the line.
        #The +3 is for indent and the extra parenthesis at the start of the call.
        line = []
        length = 0
        result = []
        for param in paramlist:
            extra = len(list(param))
            if length + extra + 2 + spacing > 90:
                result.append(", ".join(line) + ", &")
                line = [ param ]
                length = extra + 2
            else:
                line.append(param)
                length += extra + 2

        #Add on the remaining bits of the line
        result.append(", ".join(line))

        return "\n{}".format(" ".join([ "" for i in range(spacing + 3)])).join(result)
                
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

