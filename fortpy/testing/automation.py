"""Enables unit tests to be automated by generating the XML files
that the testing framework soaks up and runs. We store the choices that
the human made in a separate, serialized class instance and use it
to generate the testing XML file.
"""
from fortpy import msg
def has_outputs(wizard, restrict=False):
    """Returns True if the currently active executable in the wizard has
    all of its model outputs existing.
    """
    #This needs to happen across all the test specifications, also it
    #runs only if "runchecks" was not specified explicitly.
    if wizard.xauto.test_group is None:
        return False
    if not wizard.tauto.runchecks:
        return True
    
    testid = None if not restrict else wizard.tauto.identifier
    caseid = None if not restrict else wizard.caseauto
    codefolder = wizard.xauto.module.codefolder
    def oexists(outvar, abspath):
        from os import path
        return ((not outvar.autoclass and path.isfile(abspath)) or
                (outvar.autoclass and path.isdir(abspath)))

    exists = True
    for tid, test in wizard.xauto.test_group.tests.items():
        if not exists:
            break
        if testid is not None and tid != testid:
            continue    
        for target in test.targets:
            if target.compareto in test.outputs:
                outvar = test.outputs[target.compareto]
                if test.cases is None or len(test.cases) == 0:
                    abspaths = [outvar.abspath(codefolder, "")]
                else:
                    abspaths = [outvar.abspath(codefolder, cid) for cid in test.cases
                                if caseid is None or cid == caseid]
                exists = all([oexists[apath] for apath in abspaths])
                if not exists:
                    break

    return exists

def run(wizard, rerun=False, restrict=False):
    """Runs the unit tests for the currently active executable with output
    checking disabled. This creates the output for the unit tests so that
    the user can confirm that it is correct.

    :arg restrict: when True, only the *active* testid and caseid in the
      wizard are re-run.
    """
    #TODO: we still need to get lower-level filtering built in for the script
    #running so that individual tests and test cases can be specified.
    
    #We only want to run the test if we are missing output or it is an
    #explicit re-run.
    if not rerun and has_outputs(wizard, restrict):
        return

    #We set up the args dictionary and then import the script for running
    #tests.
    from fortpy import settings
    from fortpy.scripts.runtests import do_testing
    testid = None if not restrict else wizard.tauto.identifier
    caseid = None if not restrict else wizard.caseauto
    settings.use_test_cache = True
    args = {"quiet": True, "compiler": wizard.compiler, "codedir": wizard.xauto.codefolder,
            "rerun": wizard.xauto.name, "stagedir": None, "verbose": False,
            "templates": None, "nodebug": False, "fortpy": None, "profile": False,
            "strict": True, "testid": testid, "cases": [caseid] if caseid is not None else None}
    do_testing(args)
    
class SimpleCompleter(object):
    
    def __init__(self, options, completer=None):
        self.options = sorted(options)
        self.completer = completer
        return

    def complete(self, text, state):
        response = None
        if state == 0:
            # This is the first time for this text, so build a match list.
            if self.completer is not None:
                self.matches = self.completer(text, text, 0, len(text))
            else:
                if text:
                    self.matches = [s for s in self.options
                                    if s and s.startswith(text)]
                else:
                    self.matches = self.options[:]
                    
            #Be sure to append '?' as a global choice that is always possible.
            self.matches.append('?')
            
        # Return the state'th item from the match list,
        # if we have that many. 
        try:
            response = self.matches[state]
        except IndexError:
            response = None
        return response

#Get an input requester that is independent of python version.
try: input = raw_input
except NameError: pass

def _help_general(tag, attribute):
    """Presents help for the specified tag and attribute to the user.
    """
    

def _complete_general(prompt="$ ", options=None, completer=None, leader=None,
                      regex=None, cast=None, current=None):
    """Prompts the user for input from the specified options.
    """
    import readline
    readline.parse_and_bind("tab: complete")
    if completer is not None:
        readline.set_completer(SimpleCompleter(None, completer).complete)
    elif options is not None:
        readline.set_completer(SimpleCompleter(options).complete)
    if leader is not None:
        msg.okay(leader)
    if current is not None:
        msg.gen("The current value is set as '{}'.".format(current))
        
    choice = input(prompt)
    if choice.strip() == '?':
        return choice.strip()
    if choice.strip() == '':
        return None    
    
    if regex is not None:
        import re
        if not re.match(regex, choice):
            msg.warn("'{}' is not a valid choice.".format(choice))
            return None

    if cast is not None:
        try:
            result = cast(choice.strip())
        except ValueError:
            return None
    else:
        result = choice

    return result

def _prompt_general(header, skeys, sources=None, args=None):
    msg.okay(header)
    for i, k in enumerate(skeys):
        msg.std("{}. {}".format(i, k))
    msg.blank(1,0)
    choice = input("Selection: ")
    if re.match("[\d]+", choice.strip()):
        value = int(choice)
        if sources is None:
            return value
        else:
            if value in sources:
                sources[value](*args)
            else:
                msg.warn("choice '{}' not in list of possibilities.".format(choice))           
    else:
        msg.warn("invalid (non-numeric) choice: '{}'".format(choice))           

def _prompt_attributes(xml, attribs):
    """Prompts the user for selections for each attribute in the specified dictionary.
    Returns a list of attributes and the values selected by the user.

    :arg attribs: {"attribute": {"options": [],
                                 "completer": function(text, line, istart, iend),
                                 "leader": "Text to display to the user."}
                  }. The options and completer are xor optional; leader and prompt
      are also optional. Leader is a question to write before soliciting input,
      prompt is the text to display before the input is typed.
    """
    result = {}
    for aname, adict in attribs.items():
        if aname in xml.attrib:
            #Update the current value if one exists. The user might find it useful.
            adict["current"] = xml.attrib[aname]

        choice = _complete_general(**adict)
        if choice == '?':
            _help_general(xml.tag, aname)
            _complete_general(**adict)
        elif choice is not None:
            result[aname] = choice

    return result
        
def _has_input(wizard, parameter):
    """Determines if the specified parameter has input defined for it
    already.
    """
    if wizard.tauto is None:
        return False
    
    from fortpy.testing.elements import Assignment
    for method in wizard.tauto.methods:
        if isinstance(method, Assignment):
            if (method.name.lower() == parameter.name.lower() and
                method.value is not None):
                return True
    else:
        return False

def _input_testsource(wizard, parameter):
    """Prompts the user to set up a hook between existing unit test output
    and this parameter's input for the active test.
    """
    pass

def _input_existing(wizard, parameter):
    """Prompts the user for a file name to use as the source for the input.
    """
    pass
    
def _input_function(wizard, parameter):
    """Prompts the user for a combination of numpy functions to build the
    input data for the parameter.
    """
    pass

def _input_constant(wizard, parameter):
    """Prompts the user for a constant value to set for the parameter's input.
    """
    pass

def prompt_input(wizard, parameter):
    """Prompts for the location of the input for the specified parameter.
    """
    if _has_input(wizard, parameter):
        msg.gen("This parameter has input specified already.")
        reset = wizard.input("Would you like to add another? ").lower()[0] == 'y'
        if not reset:
            return
        else:
            msg.blank(1,0)

    skeys = ["Unit Test Output", "Existing File", "Known Function", "Constant"]
    sources = {
        0: _input_testsource,
        1: _input_existing,
        2: _input_function,
        3: _input_constant,
    }
    _prompt_general("Choose an input source:", skeys, sources, (wizard, parameter))

def _find_target(wizard, parameter):
    """Checks whether the specified parameter has a target specified for
    saving its output.
    """
    if wizard.tauto is None:
        return False

    for target in wizard.tauto:
        if target.name.lower() == parameter.name.lower():
            return target
    else:
        return None
    
def _prompt_target(wizard, parameter, target):
    """Asks the user for the details for the parameter <target> tag.
    """
    #Get a list of the attributes we need to complete on.
    skeys = ["Automate.", "Call another subroutine."]
    choice = _prompt_general("How should we save the variable?", skeys)

    #Before we start querying the user, we should get hold of the <target> we
    #are modifying.
    if target is None:
        from xml.etree.ElementTree import Element
        ltarget = Element("target", {"name": parameter.name})
    else:
        ltarget = target.xml
        
    attribs = {}    
    if choice == 1:
        attribs["generator"] = {"completer": wizard._complete_executables,
                                "leader": "Which subroutine should be called?"}
    if wizard.tauto.repeats:
        attribs["when"] = {"options": ["begin", "each", "end"],
                           "leader": "How often should we save the variable?"}

    hasmems = False
    if parameter.is_custom:
        custype = parameter.customtype
        if custype is not None:
            attribs["member"] = {"options": ["*"] + custype.members.keys(),
                                 "leader": "Which member should we save; use * for all (default)?"}
            hasmems = True
            
    selattrs = _prompt_attributes("target", attribs)
    if hasmems and ("member" not in selattrs or selattrs["member"] == "*"):
        selattrs["autoclass"] = "true"

    for k, v in selattrs.items():
        ltarget.set(k, v)
            
    if target is None:
        wizard.tauto.targets.append(TestTarget(ltarget, wizard.tauto))
        #else: we were editing the target that is already in the TestSpecifications list.

def _examine_output(wizard, parameter, target):
    """Starts the program in $EDITOR so the user can examine the variable's output
    for correctness.
    """
    pass

def _print_outpath(wizard, parameter, target):
    """Prints the path to the model output for the parameter so the user
    can work with it.
    """
    pass

def _get_casepath(wizard):
    """Returns the full path to the execution directory for the current
    test and case.
    """
    from fortpy.testing.compilers import replace
    xstage = path.join(wizard.stagedir, wizard.xauto.full_name)
    target = replace(xstage + ".[c]", wizard.compiler)
    tfolder = "{}{}".format(wizard.tauto.identifier, '.' + wizard.caseauto if wizard.caseauto is None else "")
    return path.join(target, "tests", tfolder)
    
def _start_debug(wizard, parameter, target):
    """Starts a debug session for the unit test executable of the current test
    specification and case identifier.
    """
    from fortpy.utility import which
    if which("gdb"):
        #Give the user the option of setting some breakpoints. Basically, list the
        #1) main program entry, 2) list of pre-req calls before the executable,
        #3) executable itself.
        brks = ["Main Unit Test Entry Point.", wizard.xauto.signature]
        xnames = ["main", wizard.xauto.name]
        from fortpy.testing.elements import Executable
        for method in wizard.tauto.methods:
            if not isinstance(method, Executable):
                #This is an executable that needs to be called before the main one.
                brks.append(method.signature)
                xnames.append(method.name)

        bchoice = _prompt_general("Which methods would you like breakpoints for?", brks)
        #Start a gdb session for the compiled unit test executable of the active
        #test and test case.
        xs = []
        for x in bchoice:
            if x < len(xnames):
                xs.append("-ex '{}'".format(xnames[x]))

        from os import system, path
        from fortpy.testing.compilers import compile_general          
        xstage = path.join(wizard.stagedir, wizard.xauto.full_name)
        code, success, target = compile_general(xstage, wizard.compiler, wizard.tauto.identifier,
                                                True, quiet=True)
        testpath = _get_casepath(wizard)
        cmd = "cd {}; gdb {} -ex 'run' ../../{}.x"
        system(cmd.format(testpath, ' '.join(xs), wizard.tauto.identifier))
    else:
        msg.err("The GNU debugger 'gdb' is not available on the system.")

def _set_correct(wizard, parameter, target):
    """Sets the current model output for the target and parameter as correct
    so that it will be used for future unit test comparisons.
    """
    casepath = _get_casepath(wizard)
    #We don't worry about updating the XML tag. That is done by a calling
    #procedure. The job here is to copy model output to the global output's
    #folder specified by the user for the current module, and then return the name
    #of the file in that folder.

    #As such, we need to choose a naming convention for the auto-assigned model
    #outputs. The safest thing to do is use the
    #module.executable.parameter.testid.caseid as the file name. But that is
    #*really* long and if the test folder is already buried pretty deep, it can
    #become unwieldly... On the other hand, if multiple unit tests want to compare
    #to the same model output, then we don't want to duplicate that, especially
    #considering the size of some of the files in scientific computing.
    from os import path
    varfile = path.join(casepath, target.varfile)
    newname = ""
    dst = path.join(wizard.folders[wizard.xauto.module], newname)
    if target.autoclass:
        if path.isdir(varfile):
            from fortpy.utility import copytree
            copytree(varfile, dst)
    else:
        if path.isfile(varfile):
            from fortpy.utility import copyfile
            copyfile(varfile, dst)

    return (newname, None)

def _set_existing(wizard, parameter, target):
    """Prompts the user to select an existing file to use as the model output
    for the parameter.
    """

def _prompt_model(wizard, parameter, target):
    """Prompts the user for details about the model output for comparison.
    """
    #This is essentially the prompt to setup the <output> tags. Look for the existing
    #<output> in the test specification.
    if target.compareto is None or target.compareto is not in wizard.tauto.outputs:
        from xml.etree.ElementTree import Element
        output = Element("output", {"identifier": "{}.auto".format(target.name)})
        newoutput = True
    else:
        output = wizard.tauto.outputs[target.compareto].xml
        newoutput = False

    skeys = ["Automate.", "Set a constant value.", "Don't compare to model output."]
    choice = _prompt_general("Choose the source of the model output:", skeys)
    attribs = {
        "tolerance": {"leader": "Enter % accuracy required for comparisons: 0.0 to 1.0 (default)?",
                      "cast": float}
    }
    #Keep track of whether we actually need to set model output for this parameter.
    skip = choice == 2
    
    if choice == 0:
        #The biggest problem we have here is that we need to specify the file to use
        #for the model output. We could have the user search for one, but really the
        #way this works is that we need to compile the test, run it without checks and
        #then present the user with the output to verify for each test case. Since the
        #automator presents the input parameters first and then allows the targets to
        #be established without model outputs first, we should be able to compile and
        #run the test (if it hasn't been done already).
        rdict = {
            1: _examine_output,
            2: _print_outpath,
            3: _start_debug,
            4: _set_correct,
            5: _set_existing
        }
        rkeys = ["Re-run the tests to re-create output(s).",
                 "Examine the variables' output(s) for correctness.",
                 "Print the location of the model output(s).",
                 "Start the debugger for the unit test program.",
                 "Set the variable output(s) as correct.",
                 "Specify an existing file as model output.",
                 "Exit the correction loop."]
        varfile = None
        correct = True
        runonce = False

        while correct:
            msg.blank()
            if has_outputs(wizard, True):
                rchoice = _prompt_general("The model output for the active test case exists.\n"
                                          "What would you like to do?", rkeys)
                if rchoice in rdict:
                    varfile = rdict[rchoice](wizard, parameter, target)
                    if rchoice == 4:
                        correct = False
                elif rchoice == 0:
                    run(wizard, True, True)
                elif rchoice == 5:
                    correct = False
            else:
                #First run the tests to generate the output, then present it to the user
                #so that they can check it is right.
                if not runonce:
                    msg.info("Running the unit test to generate output for examination.")
                    run(wizard, False, True)
                    runonce = True
                else:
                    msg.warn("the model outputs weren't generated by running the unit test.\n"
                             "Check error messages.")
                    correct = False
            
        if varfile is not None:
            output.set("file", varfile)
                
        if "autoclass" in target.xml.attrib and target.xml.attrib["autoclass"] == "true":
            output.set("autoclass", "true")
            if "tolerance" in selattrs:
                output.set("actolerance", selattrs["tolerance"])
    elif choice == 1:
        attribs["value"] = {"leader": "Enter the correct value; can be any valid python code."}

    if skip:
        #We need to remove an output if one already exists; otherwise do nothing.
        if target.compareto is not None:
            if target.compareto in wizard.tauto.outputs:
                del wizard.tauto.outputs[target.compareto]
            target.compareto = None
    else:
        #Prompts for the last, sundry attributes on the <output> tag.
        selattrs = _prompt_attributes("output", attribs)
        for k, v in selattrs.items():
            output.set(k, v)
            
        if newoutput:
            target.compareto = output.attrib["identifier"]
            wizard.tauto.outputs[target.compareto] = TestOutput(output)

def prompt_target(wizard, parameter, auto=True):
    """Prompts the user whether to save the output parameter's values to file.
    """
    target = _find_target(wizard, parameter)
    if not auto:
        #Give the user the option of overriding previous choices.
        if target is not None:
            msg.gen("This parameter is being saved already.")
            retarg = wizard.input("Would you like to re-specify the details?")[0].lower() == 'y'
            msg.blank(1,0)
        else:
            retarg = True
    else:
        #We only care about getting new things that we haven't got already
        retarg = target is None

    if retarg:
        _prompt_target(wizard, parameter, target)    
            
def prompt_output(wizard, parameter, auto=True):
    """Prompts for the setting up of model output for the specified parameter.
    """
    target = _find_target(wizard, parameter)
    if not auto:
        #Give the user the option of overriding previous choices.
        if target.compareto is not None:
            msg.gen("This parameter is already being compared to model output.")
            remodel = wizard.input("Would you like to re-specify the details?")[0].lower() == 'y'
            msg.blank(1,0)
        else:
            remodel = True
    else:
        #We only care about getting new things that we haven't got already
        remodel = target.compareto is None

    if remodel:
        _prompt_model(wizard, parameter, target)
        
class Hook(object):
    """Represents the input/output selection of the developer for a single
    parameter in an executable's signature."""
    def __init__(self, parameter, assignments, allocates):
        """Initializes the specified parameter for automated testing using
        the input/output settings specified.
        """
        
    
class Wizard(object):
    """Represents the choices that the developer made to hook the unit
    tests up together. A wizard only applies to a single code library.
    """
    def __init__(self, coderoot):
        """Loads the wizard and code parser for automating the tests.

        :arg coderoot: the full path to the directory housing all the code files.
        """
        self.xtests = {}
        """A dictionary of test settings indexed by module.executable name.
        Values are dictionaries keyed by parameter name with values being
        instances of Hook: one for each parameter in the executable's
        signature."""
        self.stagedir = None
        """The directory to stage the unit tests in relative to the code root.
        """
        self.folder = None
        """The folder in which the unit test input and model output files will
        be stored for this module's unit tests.
        """
        self.coderoot = coderoot
        """The full path to the directory housing all the code files."""
        self.globxml = {}
        """A dictionary of the <globals> tag indexed for each module in the
        code directory; key -> module name; value -> dict of <globals> relevant tags.
        """
        self.parser = self._load_parser()
        """An instance of fortpy.code.CodeParser for interacting with the code elements.
        """
        
    def __getstate__(self):
        """Cleans up the object so that it can be pickled without any pointer
        issues.
        """
        odict = self.__dict__.copy() # copy the dict since we change it
        del odict["parser"]
        return odict

    def __setstate__(self, dict):
        self.__dict__.update(dict)
        self.parser = self._load_parser()

    def _load_parser(self):
        """Loads the code parser instance for interacting with the code elements
        that comprise the modules in the directory.
        """
        
    def save(self):
        """Saves the state of the wizard to binary in the code root as .fpytests
        """
        from os import path
        from fortpy import msg
        fullpath = path.join(self.coderoot, ".fpytests")
        with open(fullpath, 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        msg.okay("Saved wizard data to {}".format(fullpath))

    @staticmethod
    def load(cls):
        """Loads the state of the wizard from pickle in the code root as .fpytests
        """
        from os import path
        fullpath = path.join(self.coderoot, ".fpytests")
        with open(fullpath, 'rb') as f:
            try:
                gc.disable()
                result = pickle.load(f)
            finally:
                gc.enable()

        return result
