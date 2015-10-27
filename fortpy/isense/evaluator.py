from fortpy.elements import Function, Subroutine, CustomType, ValueElement
from fortpy.elements import Module, Executable, Interface
from fortpy.docelements import DocElement
from . import cache
from .classes import Completion

class Evaluator(object):
    """Uses the user context and code parsers to perform code
    completion and definition lookups.
    """
    def __init__(self, context, pos):
        """Create a new evaluator for the given context.
        
        :arg context: an instance of UserContext for the given source
          and position.
        """
        self.context = context
        self.line, self.column = pos

        #Initialize the lazy variables
        self._names = None
        self._possible = None

    @property
    def element(self):
        """Finds the instance of the parsed element in the code parser
        using the name given in the context."""
        return self.context.element
            
    @property
    def names(self):
        """Returns a list of possible completions for the symbol under the
        cursor in the current user context."""
        #This is where the context information is extremely useful in 
        #limiting the extent of the search space we need to examine.
        if self._names is None:
            attribute = self.get_attribute()
            if self.context.module is not None:
                symbol = self.context.symbol
                fullsymbol = self.context.full_symbol
                self._names = self._complete_el(symbol, attribute, fullsymbol)
            else:
                self._names = []

        return self._names
        
    def bracket_complete(self):
        """Returns a function call signature for completion whenever
        a bracket '(' is pressed."""
        #The important thing to keep track of here is that '(' can be
        #pushed in the following places:
        # - in a subroutine/function definition
        # - when calling a function or subroutine, this includes calls within
        #   the argument list of a function, e.g. function(a, fnb(c,d), e)
        # - when specifying the dimensions of an array for element selection
        # - inside of if/do/while blocks.

        #If there is a \s immediately preceding the bracket, then short_symbol
        #will be null string and we are most likely doing arithmetic of some
        #sort or a block statement.
        if self.context.short_symbol == "":
            return {}

        line = self.context.current_line.lower()[:self.context.pos[1]-1]
        if "subroutine" in line or "function" in line:
            #We are in the definition of the subroutine or function. They
            #are choosing the parameter names so we don't offer any suggestions.
            return {}
            
        #Sometimes people do if() with the condition immediately after
        #the keyword, this regex will catch that.
        symbol = self.context.short_symbol.lower()
        if symbol in ["if", "do", "while", "elseif", "case", "associate"]:
            return {}
        
        #All that should be left now are dimensions and legitimate function
        #calls.
        fullsymbol = self.context.short_full_symbol.lower()
        return self._bracket_complete_sig(symbol, fullsymbol)

    def _bracket_complete_sig(self, symbol, fullsymbol):
        """Returns the call signature and docstring for the executable
        immediately preceding a bracket '(' that was typed."""
        if symbol != fullsymbol:
            #We have a sym%sym%... chain and the completion just needs to
            #be the signature of the member method.
            target, targmod = self._get_chain_parent_symbol(symbol, fullsymbol)

            if symbol in target.executables:
                child = target.executables[symbol]
                return self._compile_signature(child.target, child.name)
            elif symbol in target.members:
                #We are dealing with a dimension request on an array that
                #is a member of the type.
                child = target.members[symbol]
                return self._bracket_dim_suggest(child)
            else:                
                return {}
        else:
            #We must be dealing with a regular executable or builtin fxn
            #or a regular variable dimension.
            iexec = self._bracket_exact_exec(symbol)
            if iexec is not None:
                #It is indeed a function we are completing for.
                return self._compile_signature(iexec, iexec.name)
            else:
                #We need to look at local and global variables to find the
                #variable declaration and dimensionality.
                ivar = self._bracket_exact_var(symbol)
                return self._bracket_dim_suggest(ivar)
            
    def _bracket_dim_suggest(self, variable):
        """Returns a dictionary of documentation for helping complete the
        dimensions of a variable."""
        if variable is not None:
            #Look for <dimension> descriptors that are children of the variable
            #in its docstrings.
            dims = variable.doc_children("dimension", ["member", "parameter", "local"])
            descript = str(variable)
            if len(dims) > 0:
                descript += " | " + " ".join([DocElement.format_dimension(d) for d in dims])

            return dict(
                params=[variable.dimension],
                index=0,
                call_name=variable.name,
                description=descript,
            )        
        else:
            return []

    def get_definition(self):
        """Checks variable and executable code elements based on the current
        context for a code element whose name matches context.exact_match
        perfectly.
        """
        #Check the variables first, then the functions.
        match = self._bracket_exact_var(self.context.exact_match)
        if match is None:
           match = self._bracket_exact_exec(self.context.exact_match)

        return match

    def _bracket_exact_var(self, symbol):
        """Checks local first and then module global variables for an exact
        match to the specified symbol name."""
        if isinstance(self.element, Executable):
            if symbol in self.element.parameters:
                return self.element.parameters[symbol]            
            if symbol in self.element.members:
                return self.element.members[symbol]

        if symbol in self.element.module.members:
            return self.element.module.members[symbol]

        return None

    def _bracket_exact_exec(self, symbol):
        """Checks builtin, local and global executable collections for the
        specified symbol and returns it as soon as it is found."""
        if symbol in self.context.module.executables:
            return self.context.module.executables[symbol]

        if symbol in self.context.module.interfaces:
            return self.context.module.interfaces[symbol]

        if symbol in cache.builtin:
            return cache.builtin[symbol]

        #Loop through all the dependencies of the current module and see
        #if one of them is the method we are looking for.
        return self.context.module.get_dependency_element(symbol)

    def _compile_signature(self, iexec, call_name):
        """Compiles the signature for the specified executable and returns
        as a dictionary."""
        if iexec is not None:
            summary = iexec.summary
            if isinstance(iexec, Function):
                summary = iexec.returns + "| " + iexec.summary
            elif isinstance(iexec, Subroutine) and len(iexec.modifiers) > 0:
                summary = ", ".join(iexec.modifiers) + " | " + iexec.summary
            elif isinstance(iexec, Interface):
                summary = iexec.describe()
            else:
                summary = iexec.summary

            #Add the name of the module who owns the method. Useful in case the
            #same executable is defined in multiple modules, but only one is
            #referenced in the current context.
            if iexec.parent is not None:
                summary += " | MODULE: {}".format(iexec.module.name)
            else:
                summary += " | BUILTIN"

            return dict(
                params=[p.name for p in iexec.ordered_parameters],
                index=0,
                call_name=call_name,
                description=summary,
            )        
        else:
            return []

    def in_function_call(self):
        """This function is called whenever the cursor/buffer goes idle for
        a second. The real workhorse of the intellisense. Decides what kind
        of intellisense is needed and returns the relevant response."""
        #See if we are calling a function inside a module or other function.
        result = []
        if (self.context.el_section == "body" and
            self.context.el_call in [ "sub", "fun" ]):
            #Do a signature completion for the function/subroutine call.
            result = self.signature()

        if result == []:
            return self.complete()

        return result

    def signature(self):
        """Gets completion or call signature information for the current cursor."""
        #We can't really do anything sensible without the name of the function
        #whose signature we are completing.
        iexec, execmod = self.context.parser.tree_find(self.context.el_name, 
                                              self.context.module, "executables")
        if iexec is None:
            #Look in the interfaces next using a tree find as well
            iexec, execmod = self.context.parser.tree_find(self.context.el_name, self.context.module,
                                                  "interfaces")
        if iexec is None:
            return []

        return self._signature_index(iexec)

    def _signature_index(self, iexec):
        """Determines where in the call signature the cursor is to decide which
        parameter needs to have its information returned for the intellisense.
        """
        #Find out where in the signature the cursor is at the moment.
        call_index = self.context.call_arg_index
        if call_index is not None:
            #We found the index of the parameter whose docstring we want
            #to return.
            param = iexec.get_parameter(call_index)
            paramlist = [ p.name for p in iexec.ordered_parameters ]
            paramlist[call_index] = "*{}*".format(paramlist[call_index])

            if not isinstance(param, list) and param is not None:
                #Write a nice description that includes the parameter type and
                #intent as well as dimension.
                summary = "{} | {}".format(str(param), param.summary)
                #We also want to determine if this parameter has its value changed
                #by the function we are completing on.
                changedby = iexec.changed(param.name)
                if changedby is not None:
                    summary += " *MODIFIED*"
            elif isinstance(param, list) and len(param) > 0:
                act_type = []
                for iparam in param:
                    if iparam is not None and iparam.strtype not in act_type:
                        act_type.append(iparam.strtype)

                act_text = ', '.join(act_type)
                summary = "SUMMARY: {} | ACCEPTS: {}".format(param[0].summary, act_text)

                if iexec.changed(param[0].name):
                    summary += " | *MODIFIED*"
            else:
                summary = "No matching variable definition."

            #Add the name of the module who owns the executable.
            summary += " | MODULE: {}".format(iexec.module.name)

            return dict(
                params=paramlist,
                index=call_index,
                call_name=self.context.el_name,
                description=summary,
            )
        else:
            return result

    def complete(self):
        """Gets a list of completion objects for the symbol under the cursor."""
        if self._possible is None:
            self._possible = []
            for possible in self.names:
                c = Completion(self.context, self.names[possible], len(self.context.symbol))
                self._possible.append(c)

        return self._possible

    def _symbol_in(self, symbol, name):
        """Checks whether the specified symbol is part of the name for completion."""
        lsymbol = symbol.lower()
        lname = name.lower()
        return lsymbol == lname[:len(symbol)] or "_" + lsymbol in lname
        
    def _complete_el(self, symbol, attribute, fullsymbol):
        """Suggests a list of completions based on the el_* attributes
        of the user_context."""
        if symbol != fullsymbol:
            #We have a sym%sym%... chain and the completion just needs to
            #be a member variable or method of the type being referenced.
            return self._complete_type_chain(symbol, fullsymbol)

        if self.context.el_section == "params":
            #They are in the process of defining a new executable and are
            #picking the names themselves, return normal word complete.
            return self._complete_word(symbol, attribute)
        elif self.context.el_section == "body":
            if self.context.el_call in ["sub", "fun"]:
                return self._complete_sig(symbol, attribute)
            else:
                return self._complete_word(symbol, attribute)
        else:
            return self._complete_word(symbol, attribute)

    def _get_chain_parent_symbol(self, symbol, fullsymbol):
        """Gets the code element object for the parent of the specified
        symbol in the fullsymbol chain."""
        #We are only interested in the type of the variable immediately preceding our symbol
        #in the chain so we can list its members.
        chain = fullsymbol.split("%")

        #We assume that if symbol != fullsymbol, we have at least a % at the end that 
        #tricked the symbol regex.
        if len(chain) < 2:
            return ([], None)
        previous = chain[-2].lower()

        #Now we need to use the name of the variable to find the actual type name
        target_name = ""
        if previous in self.element.members:
            target_name = self.element.members[previous].kind
        #The contextual element could be a module, in which case it has no parameters
        if hasattr(self.element, "parameters") and previous in self.element.parameters:
            target_name = self.element.parameters[previous].kind

        if target_name == "":
            return (None, None)

        return self.context.parser.tree_find(target_name, self.context.module, "types")

    def _complete_type_chain(self, symbol, fullsymbol):
        """Suggests completion for the end of a type chain."""
        target, targmod = self._get_chain_parent_symbol(symbol, fullsymbol)
        if target is None:
            return {}

        result = {}
        #We might know what kind of symbol to limit the completion by depending on whether
        #it was preceded by a "call " for example. Check the context's el_call
        if symbol != "":
            if self.context.el_call != "sub":
                for mkey in target.members:
                    if self._symbol_in(symbol, mkey):
                        result[mkey] = target.members[mkey]

            for ekey in target.executables:
                if (self._symbol_in(symbol, ekey)):
                    if self.context.el_call == "sub":
                        if (isinstance(target.executables[ekey], Subroutine)):
                            result[ekey] = target.executables[ekey]
                    else:
                        if (isinstance(target.executables[ekey], Function)):
                            result[ekey] = target.executables[ekey]
        else:
            if self.context.el_call != "sub":
                result.update(target.members)
                subdict = {k: target.executables[k] for k in target.executables
                           if isinstance(target.executables[k].target, Function)}
                result.update(subdict)
            else:
                subdict = {k: target.executables[k] for k in target.executables
                           if isinstance(target.executables[k].target, Subroutine)}
                result.update(subdict)

        return result

    def _complete_sig(self, symbol, attribute):
        """Suggests completion for calling a function or subroutine."""
        #Return a list of valid parameters for the function being called
        fncall = self.context.el_name
        iexec, execmod = self.context.parser.tree_find(fncall, self.context.module, "executables")

        if iexec is None:
            #Try the interfaces as a possible executable to complete.
            iexec, execmod = self.context.parser.tree_find(fncall, self.context.module, "interfaces")

        if iexec is not None:
            if symbol == "":
                return iexec.parameters
            else:
                result = {}
                for ikey in iexec.parameters:
                    if self._symbol_in(symbol, ikey):
                        result[ikey] = iexec.parameters[ikey]
                return result
        else:
            return self._complete_word(symbol, attribute)

    def _complete_word(self, symbol, attribute):
        """Suggests context completions based exclusively on the word
        preceding the cursor."""
        #The cursor is after a %(,\s and the user is looking for a list
        #of possibilities that is a bit smarter that regular AC.
        if self.context.el_call in ["sub", "fun", "assign", "arith"]:
            if symbol == "":
                #The only possibilities are local vars, global vars or functions
                #presented in that order of likelihood. 
                return self._complete_values()
            else:
                #It is also possible that subroutines are being called, but that
                #the full name hasn't been entered yet.
                return self._complete_values(symbol)
        else:
            return self.context.module.completions(symbol, attribute, True)        

    def _generic_filter_execs(self, module):
        """Filters the specified dict of executables to include only those that are not referenced
        in a derived type or an interface.
        
        :arg module: the module whose executables should be filtered.
        """
        interface_execs = []
        for ikey in module.interfaces:
            interface_execs.extend(module.interfaces[ikey].procedures)

        return {k: module.executables[k] for k in module.executables if 
                ("{}.{}".format(module.name, k) not in interface_execs and 
                 not module.executables[k].is_type_target)}

    def _complete_values(self, symbol = ""):
        """Compiles a list of possible symbols that can hold a value in
        place. These consist of local vars, global vars, and functions."""
        result = {}
        #Also add the subroutines from the module and its dependencies.
        moddict = self._generic_filter_execs(self.context.module)
        self._cond_update(result, moddict, symbol)
        self._cond_update(result, self.context.module.interfaces, symbol)
        for depend in self.context.module.dependencies:
            if depend in self.context.module.parent.modules:
                #We don't want to display executables that are part of an interface, or that are embedded in
                #a derived type, since those will be called through the type or interface
                filtdict = self._generic_filter_execs(self.context.module.parent.modules[depend])
                self._cond_update(result, filtdict, symbol)
                self._cond_update(result, self.context.module.parent.modules[depend].interfaces, symbol)

        #Add all the local vars if we are in an executable
        if (isinstance(self.context.element, Function) or 
            isinstance(self.context.element, Subroutine)):
            self._cond_update(result, self.element.members, symbol)
        #Next add the global variables from the module
        if self.context.module is not None:
            self._cond_update(result, self.context.module.members, symbol)
            #Next add user defined functions to the mix
            for execkey in self.context.module.executables:
                iexec = self.context.module.executables[execkey]
                if isinstance(iexec, Function) and self._symbol_in(symbol, iexec.name):
                    result[iexec.name] = iexec

        #Finally add the builtin functions to the mix. We need to add support
        #for these in a separate file so we have their call signatures.
        if symbol == "":
            #Use the abbreviated list of most common fortran builtins
            self._cond_update(result, cache.common_builtin, symbol)
        else:
            #we can use the full list as there will probably not be that
            #many left over.
            self._cond_update(result, cache.builtin, symbol)           
        
        return result

    def _cond_update(self, first, second, symbol, maxadd = -1):
        """Overwrites the keys and values in the first dictionary with those
        of the second as long as the symbol matches part of the key.

        :arg maxadd: specifies a limit on the number of entries from the second
          dict that will be added to the first."""
        if symbol != "":
            added = 0
            for key in second:
                if self._symbol_in(symbol, key) and (maxadd == -1 or added <= maxadd):
                    first[key] = second[key]
                    added += 1
        else:
            first.update(second)
        
    def get_attribute(self):
        """Gets the appropriate module attribute name for a collection 
        corresponding to the context's element type."""
        attributes = ['dependencies', 'publics', 'members', 
                          'types', 'executables']
        #Find the correct attribute based on the type of the context
        if self.context.el_type in [Function, Subroutine]:
            attribute = attributes[4]
        elif self.context.el_type == CustomType:
            attribute = attributes[3]
        else:
            attribute = attributes[2]

        return attribute
