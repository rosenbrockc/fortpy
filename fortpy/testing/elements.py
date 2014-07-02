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

    def compare(self, element):
        """Determines whether the specified element defines the same variable
        type and name as this one.

        The two are considered equal if they have the same:
         - name
         - type
         - dimensionality
         - allocatable/pointer modifiers
        ."""
        if self.attributes["name"].lower() != element.attributes["name"].lower():
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
        return self.element.attributes

    def definition(self):
        """Returns the fortran declaration that defines a global variable."""
        result = []
        if "type" not in self.attributes or "name" not in self.attributes:
            print("FATAL: required variable for execution missing some attributes. {}".format(
                self.attributes))
            exit(1)

        result.append(self.attributes["type"])
        if "kind" in self.attributes and self.attributes["kind"] is not None:
            result.append("({})".format(self.attributes["kind"]))

        if "modifiers" in self.attributes:
            mods = self.attributes["modifiers"].split(",")
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
    :attr prereq: if true, and embedde method on a derived type will be treated
      along with all its prereqs as part of the execution chain.
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
        self.prereq = None
        self.paramlist = None

        self._derived_type = None
        self._codes = {
            "vars": self._code_vars,
            "init": self._code_init,
            "assign": self._code_assign,
            "after": self._code_after,
            "before": self._code_before
        }

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

    @property
    def allocatable(self):
        """Returns true if the variable assigned by this object should be allocated."""
        return (self.parent.allocate and 
                ("allocatable" in self.parent.variable.modifiers or
                 "pointer" in self.parent.variable.modifiers))        

    def copy(self, coderoot, testroot, case):
        """Copies the input files needed for this value to set a variable.

        :arg coderoot: the full path to folder with the code files.
        :arg testroot: the full path the folder where the test is running.
        :arg case: the case id for multi-case testing.
        """
        #We only want to do the copy if we are assigning a value from a file.
        if self.folder is not None:
            source = path.join(coderoot, self.folder[2:], self.filename.format(case))
            if self.rename is not None:
                target = path.join(testroot, self.rename)
            else:
                target = path.join(testroot, self.filename.format(case))
            copyfile(source, target)

    def check_prereqs(self):
        """Checks whether this value element requires an embedded method to be
        added to the prereq chain of the method writer."""
        if self.embedded is not None:
            #The embedded method could use a pointer so that the method we are
            #calling is *not* the name of the actual method. We need to find the
            #instance of the TypeExecutable to locate its actual target.
            var = self.parent.variable
            target = var.kind
            self._derived_type = self.parent.parser.tree_find(target, var.module, "types")

            if derived_type is None:
                raise ValueError("The type for embedded method {} cannot be found.".format(self.embedded))

            if self.prereq:
                if derived_type.pointsto is not None:
                    key = "{}.{}".format(derived_type.module, derived_type.pointsto)
                else:
                    key = "{}.{}".format(derived_type.module, derived_type.name)

                self.parent.add_prereq(key, self.parent.element)
 
    def code(self, lines, position, spacer):
        """Appends the code lines to initialize the parent variable.

        :arg lines: the list of strings to append to.
        :arg position: one of ['vars', 'init', 'assign', 'before', 'after'] indicating 
          where in the fortran program the code will be appended.
        """
        if self.parent.variable is not None:
            self._codes[position](lines, spacer)
        else:
            raise ValueError("Trying to assign a value to an unknown variable {}.".format(self.iid))

    def _code_before(self, lines, spacer):
        """Calls _code_setvar() if we *are* in repeat mode."""
        if self.repeats:
            self._code_setvar(lines, spacer)

    def _code_assign(self, lines, spacer):
        """Calls _code_setvar() if we are *not* in repeat mode."""
        if not self.repeats:
            self._code_setvar(lines, spacer)

    def _code_setvar(self, lines, spacer):
        """Appends code for assigning the value of the parent variable using
        this value specification."""
        if self.constant is not None:
            lines.append("{}{} = {}".format(spacer, self.parent.name, self.constant))
        elif self.function is not None:
            lines.append("{}{} = {}".format(spacer, self.parent.name, self.function))
        elif self.embedded is not None:
            self._code_embedded(lines, spacer)
        else:
            self._code_file(lines, spacer)
            
    def _code_embedded(self, lines, spacer):
        """Appends code for calling an embedded method in a derived type,
        optionally including all its dependencies."""
        if not self.prereqs and self._derived_type is not None:
            #This is a simple exercise in calling the embedded method. We just 
            #need to track down the parameters list.
            if self.paramlist is None:
                target = self._derived_type.target
                params = ", ".join(target.ordered_parameters)
            else:
                params = self.paramlist
            lines.append("{}{}%{}({})".format(spacer, self.parent.name, 
                                          self.embedded, paramlist))
        #if it does have pre-reqs, they will be handled by the method writer and we
        #don't have to worry about it.

    def _code_file(self, lines, spacer):
        """Appends code for initializing the value of a variable that is a
        scalar, vector or 2D matrix from a file.
        """
        if self.filename is not None:
            flines = []
            flines.append("open(fpy_newunit({}_funit), ".format(self.iid) + 
                          "file='{}')".format(self.xname))

            #We are working with a vector or scalar. Check the dimensionality of
            #the actual variable and see if it needs to be allocated.
            if (self.parent.variable.dimension is not None and
                self.allocatable and self.parent.variable.D == 1):
                flines.append("allocate({}({}_nvalues))".format(self.parent.name, self.iid))
                flines.append("read({}_funit, *) {}".format(self.iid, self.parent.name))
                    
            #Now handle the case where the input file fills a 2D variable.
            if (self.parent.variable.dimension is not None and
                self.parent.variable.D == 2 and self.allocatable):                  
                allocstr = "allocate({0}({1}_nlines, {1}_nvalues))"
                flines.append(allocstr.format(self.parent.name, self.iid))

                fmtstr = "read({0}_funit,*)({1}(nrow,:), nrow =1, {0}_nlines)"
                flines.append(fmtstr.format(self.iid, self.parent.name))

            flines.append("close({}_funit)\n".format(self.iid))
            lines.extend([ spacer + l for l in flines])

    def _code_after(self, lines, spacer):
        """Appends code for deallocating a variable that was assigned a value.
        This is useful so that we can reallocate it again in repeat mode.
        """
        if (self.repeats and self.parent.allocate and 
            ("allocatable" in self.parent.variable.modifiers or
             "pointer" in self.parent.variable.modifiers)):
            lines.append("{}deallocate({})".format(spacer, self.parent.name))

    def _code_init(self, lines, spacer):
        """Appends code to initialize any variables we need for file read
        operations later on.
        """
        if self.filename is not None:
            lines.append("{}!Line/value counting for {}.".format(spacer, self.xname))
            lines.append("{}call fpy_linevalue_count('{}', ".format(spacer, self.xname) + 
                         "{2}, '{0}', {1}_nlines, {1}_nvalues)".format(self.commentchar, self.iid,
                                                                      len(self.xname)))

    def _code_vars(self, lines, spacer):
        """Appends lines to declare any variables we need for file read
        operations later on.
        """
        if self.filename is not None:
            lines.append("{}!Vars for initializing variable {} ".format(spacer, self.parent.name) +
                         "from file {}".format(self.filename))
            lines.append("{0}integer :: {1}_nlines, {1}_nvalues, {1}_funit".format(spacer, self.iid))

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
        else:
            self.commentchar = '#'
        if "paramlist" in self.xml.attrib:
            self.paramlist = self.xml.attrib["paramlist"]
        if "rename" in self.xml.attrib:
            self.rename = self.xml.attrib["rename"]
        else:
            self.commentchar = "#"
        if "repeats" in self.xml.attrib:
            self.repeats = self.xml.attrib["repeats"].lower() == "true"
        else:
            self.repeats = False
        if "prereq" in self.xml.attrib:
            self.prereq = self.xml.attrib["prereq"].lower() == "true"
        else:
            self.prereq = False

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

class Assignment(object):
    """Represents an instance value assignment as part of a pre-req chain.

    :attr element: the DocElement instance for the <assignment> tag that this
      assignment object represents.
    :attr parent: the MethodFinder instance who owns this assignment.
    :attr methods: a list that will never get used; needed for compatibility
      with method.MethodWriter._order_dependencies recursion.
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

        self.name = None
        self.conditionals = []
        self.values = {}
        self.allocate = None
        self.value = None

        self._variable = None

        #Instances can have their values changed by:
        # - direct assignment to a constant value or function call.
        # - direct assignment as part of a set of logical tests.
        # - from the contents of a file.
        # - by calling an embedded method within a derived type.
        self._parse_xml()

    @property
    def variable(self):
        """Returns the ValueElement instance that this assignment refers to."""
        if self._variable is None:
            #Search the global declarations for the variable and then track down
            #the code element that it represents
            self._variable = MethodFinder.recurse_find_variable(self.name.lower(), self.parent)

        return self._variable

    @property
    def attributes(self):
        """Provides one-level-up access to the XML elements attributes collection."""
        return self.element.attributes

    def check_prereqs(self):
        """Checks all the values this assignment may use to make sure they reference
        a valid derived type. If an AssignmentValue specifies to check pre-reqs, 
        they are added to the pre-req chain."""
        for v in self.values:
            self.values[v].check_prereqs()

    def copy(self, coderoot, testroot, case):
        """Copies the input files needed for this assignment to set a variable.

        :arg coderoot: the full path to folder with the code files.
        :arg testroot: the full path the folder where the test is running.
        :arg case: the case id for multi-case testing.
        """
        #Just call copy on all the child value objects.
        for v in self.values:
            self.values[v].copy(coderoot, testroot, case)
        
    def code(self, lines, position, spacer):
        """Appends the code for this assignment to function correctly at the
        specified position in the code file.

        :arg lines: the list of strings that form the code file.
        :arg position: one of ['vars', 'init', 'assign', 'before', 'after'] corresponding
          to variable declaration, initialization, assignment and cleanup.
        """
        if position in ["init", "vars", "after"]:
            for v in self.values:
                self.values[v].code(lines, position, spacer)
        elif position in ["assign", "before"]:
            #We just need to check whether this assignment had an explicit
            #value or if it has a list of conditionals.
            if self.value is None:
                for c in self.conditionals:
                    c.code(lines, position, spacer)
            else:
                if self.value in self.values:
                    self.values[self.value].code(lines, position, spacer)
        
    def _parse_xml(self):
        """Parses attributes and child tags from the XML element."""
        if "name" in self.element.attributes:
            self.name = self.element.attributes["name"]
        else:
            raise ValueError("'name' is a required attribute for <assignment> tags.")

        if "value" in self.attributes:
            self.value = self.attributes["value"]
        if "allocate" in self.attributes:
            self.allocate = self.attributes["allocate"] == "true"
        else:
            self.allocate = True

        for child in self.element.xml:
            if child.tag == "value":
                val = AssignmentValue(child, self)
                self.values[val.identifier] = val
            elif child.tag == "conditionals":
                self.conditionals.append(AssignmentConditional(child, self))
         
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

    def _parse_attributes(self):
        """Extracts implemented attributes from the test target."""
        if "name" in self.xml.attrib:
            self.name = self.xml.attrib["name"]
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

        self.testgroup = testgroup

        #The name of the variable that holds the number of times to repeat
        #the execution of the method being tested.
        self._repeats = None
        #If any of the inputs runs in constant mode, this variable with be True.
        self._constant = None

        self._parse_attributes()
        self._parse_children()

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

    def code_validate(self, lines, variables, spacer):
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

    def code_vars(self, lines, variables, spacer):
        """Appends all the code required by all the input file specifications
        in order to run the unit test to the list.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run."""
        for i in self.inputs:
            i.code_vars(lines, spacer)

    def code_init(self, lines, variables, spacer):
        """Appends all code required to initialize all the constant-mode input
        file specifications and save starting values of targets.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run."""
        for i in self.inputs:
            i.code_init(lines, spacer)

        for t in self.targets:
            t.code(lines, variables, "begin", spacer)

    def code_before(self, lines, variables, spacer):
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

    def code_after(self, lines, variables, spacer):
        """Appends the code required to save the value of any targets for this
        test to their files when their save frequency is set to 'each'.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run."""
        for t in self.targets:
            t.code(lines, variables, "each", spacer)

    def code_finalize(self, lines, variables, spacer):
        """Appends the code for cleaning up open file handles and saving the
        final values of the targets to their files.

        :arg lines: the list of strings that will form the contents of the
          *.f90 file being created for the executable.
        :arg variables: the dictionary of all global variables declared in
          the program to be used by any methods being run.
        """
        for t in self.targets:
            t.code(lines, variables, "end", spacer)

        for i in self.inputs:
            i.code_finalize(lines, spacer)

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
        
    def _parse_attributes(self):
        """Gets test-level attributes from the xml tag if they exist."""
        if "identifier" in self.xml.attrib:
            self.identifier = self.xml.attrib["identifier"]
        else:
            raise ValueError("'identifier' is a required attribute for <test> tags.")
        if "description" in self.xml.attrib:
            self.description = self.xml.attrib["description"]
        if "cases" in self.xml.attrib:
            self.cases = re.split(",\s*", self.xml.attrib["cases"])
        if "runtime" in self.xml.attrib:
            runtime = list(self.xml.attrib["runtime"])
            self.unit = runtime.pop()
            self.runtime = [ int(t) for t in "".join(runtime).split("-") ]
        if "timed" in self.xml.attrib:
            self.timed = self.xml.attrib["timed"] == "true"

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

        #The order in which the global declarations appear in the group.
        self._variable_order = []
        self._codes = {}
        self._find_children()
        self._parse_xml()
        self._init_codes()

    @property
    def name(self):
        """Returns the name of the underlying DocGroup."""
        return self.group.name

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
            self._codes[testid][section](lines, self.variables, spacer)

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

        self._parse_mappings()
        self._parse_variables()

    def _parse_variables(self):
        """Searches for <global> tags and adds them to the variables dict."""
        for child in self.children:
            if (child.doctype == "global" and "name" in child.attributes):
                self._global_add(child.attributes["name"], child)

        if type(self.element).__name__ in ["Subroutine", "Function"]:
            self._get_param_globals()

    def _get_param_globals(self):
        """Extracts global declarations from the parameter regular="true" attribute
        of paramater decorations."""
        for name in self.element.parameters:
            param = self.element.parameters[name]
            for doc in param.docstring:
                if doc.doctype == "parameter" and \
                   "regular" in doc.attributes and doc.attributes["regular"] == "true":
                    glob = self._global_from_param(name)
                    if glob is not None:
                        self._global_add(name, glob)

    def _global_from_param(self, name):
        """Creates a global DocElement from the specified executable's parameter "name"."""
        if name in self.element.parameters:
            param = self.element.parameters[name]
            result = DocElement(None, None)
            result.doctype = "AUTOPARAM"
            result.attributes["name"] = name
            result.attributes["type"] = param.dtype
            self._global_clean_param(result, "kind", param.kind)
            self._global_clean_param(result, "modifiers", ", ".join(param.modifiers))
            self._global_clean_param(result, "dimensions", param.dimension)
            self._global_clean_param(result, "default", param.default)
            return result
        else:
            return None

    def _global_clean_param(self, result, name, value):
        """Adds the specified parameter name and value to the result if the value is not None."""
        if value is not None:
            result.attributes[name] = value

    def _global_add(self, name, doc):
        """Adds a global declaration for the specified DocElement if it isn't already
        represented in the list."""
        if name not in self.variables:
            self.variables[name] = GlobalDeclaration(doc)
            self._variable_order.append(name)
        else:
            #We need to make sure that it is unique compared to the existing
            #one. If it isn't, stop the execution. If it is, we don't need
            #to add it again.
            existing = self.variables[name]
            if not existing.compare(doc):
                print("FATAL: variables in the calling namespace have the same name," + \
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
    :arg fatal_if_missing: specified whether the framework chokes if it can't
      find a method that is referenced as a pre-req.
    :attr executable: the code element instance of the executable found from
      the module that the method belongs to.
    :attr main: when true, this instance is for the *main* method being unit tested.
    """
    def __init__(self, identifier, parser, element, fatal_if_missing = True, main=False):
        self.identifiers = identifier.split(".")
        self.name = identifier
        self.methods = []
        self.element = element
        self.main = main

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
    def tags(self):
        """Returns a list of the unit testing related child tags of the testing
        group associated with the executable in this MethodFinder."""
        if self.group is not None:
            return self.group.children
        else:
            return []
        
    @property
    def attributes(self):
        """Returns the dictionary of attributes from the underlying DocElement."""
        if self.element is not None:
            return self.element.attributes
        else:
            return {}                  
                
    @staticmethod
    def recurse_find_variable(name, finder):
        """Finds the ValueElement instance for the variable of the specified name
        if it exists. Found by searching the parameter list of all methods that get
        called as part of the pre-req chain."""
        result = None

        if isinstance(finder, MethodFinder) and finder.executable is not None:
            if name in finder.executable.parameters:
                result = finder.executable.parameters[name]
            else:
                for nested in finder.methods:
                    result = MethodFinder.recurse_find_variable(name, nested)
                    if result is not None:
                        break
        return result                        

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
            lines.append("{}real(dp) :: fpy_start, fpy_end, fpy_elapsed = 0".format(spacer))

    def _code_call(self, lines, spacer, testid):
        """Appends the code to call the executable that this method finder represents."""
        #If we are timing the executable, we want to get the time before and after the
        #execution of just this method and then add it to the elapsed time.
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
                if param in self.group.mappings:
                    calllist.append(self.group.mappings[param])
                else:
                    calllist.append(param.name)
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

    def add_prereq(self, identify, tag):
        """Adds a MethodFinder instance to the current one for the specified 
        'module.executable' identity.

        :arg identify: a string of "modulename.executablename".
        :arg tag: the DocElement instance that resulted in this prereq being added.
        """
        #The only reason this method exists is to prevent an import loop problem
        #It gets used by the AssignmentValue._code_embedded().
        self.methods.append(MethodFinder(identify, self._parser, tag, self._fatal))

    def _get_prereqs(self):
        """Compiles a list of subroutines that need to run. Handles recursive calls."""
        #Make sure we have tests to run off of; some of the routines that are dependencies
        #won't necessarily have any docstrings.
        for tag in self.tags:
            if tag.doctype == "prereq" and "method" in tag.attributes:
                identify = tag.attributes["method"]
                self.methods.append(MethodFinder(identify, self._parser, tag, self._fatal))
            elif tag.doctype == "assignment":
                #This is a variable assignment. Sometimes variable values are changed
                #between calls to functions. The order in which prereq/instance elements
                #appear defines when these assignments take place. We will treat a var
                #value assignment as a method to simplify the implementation.
                self.methods.append(Assignment(tag, self))

        #Now that we have the method tree setup; we need to check each of the Assignments
        #to get their pre-reqs straightened out.
        for m in self.methods:
            if isinstance(m, Assignment):
                m.check_prereqs()

    @property
    def module(self):
        """Returns the module name for this method."""
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

