from .. import msg
import re
from difflib import SequenceMatcher
from fortpy.elements import Executable, CustomType, Module
import xml.etree.ElementTree as ET
       
class LineParser(object):
    """Parses individual lines or blocks of code to retrieve code element
    representations from them. Exposes low level regex methods for type
    testing of strings.

    :attr parent: the ModuleUpdater instance that owns this line parser.
    """
    def __init__(self, parent):
        self.parent = parent
        self.additions = []

        self._current_op = None
        #Lines and character counts for finding where matches fit in the file
        self._lines = []
        self._chars = []

    @property
    def refstring(self):
        """Returns the source code of the buffer for the current operation,
        or if there is no operation, an empty string."""
        if self._current_op is None:
            return ""
        else:
            return self._current_op.context.refstring   

    def parse(self, op, mode=None):
        """Handles the parsing of the specified operation using the current
        statements from the buffer/cache. If 'comment' is not None, the 
        statement in comment is parsed instead.

        :arg mode: override the 'insert', 'replace' or 'delete' specifier.
        """
        #We need to handle docstrings separately because they are XML and slight
        #insertions/changes can't be handled using only one line. We have to get
        #hold of the entire XML block and parse it all at once again.
        if op.state[0] is None:
            linenum, statement, length = op.curcached
        else:
            linenum, statement, length = op.curbuffer

        #If there is nothing to run, it is pointless carrying on down the
        #stack of functions before this gets realized.
        if statement != "":
            self._current_op = op
            #We may need to override the mode for the 'replace' operation.
            nmode = mode if mode is not None else op.mode

            #We are dealing with pure code, decide what type the owner is and run
            #their method for real time updating.
            if isinstance(op.element, Executable):
                op.element.rt_update(statement, linenum, nmode, 
                                     op.context.parser.modulep.xparser)
            elif isinstance(op.element, CustomType):
                op.element.rt_update(statement, linenum, nmode, 
                                     op.context.parser.modulep.tparser)
            elif isinstance(op.element, Module):
                op.element.rt_update(statement, linenum, nmode, 
                                     op.context.parser.modulep, self)

        if len(self.additions) > 0:
            #Some new instances were created from a single line signature.
            #probably any additional statements in the list being executed
            #include items that should belong to the new instance.
            for add in self.additions:
                instance, module = add
                self.parent.update_instance_extent(instance, module, op)

            self.additions = []       

    def absolute_charindex(self, string, start, end):
        """Finds the absolute character index of the specified regex match
        start and end indices in the *buffer* refstring."""
        search = string[start:end]
        abs_start = self.refstring.index(search)
        return abs_start, (end - start) + abs_start

    def charindex(self, line, char, context):
        """Determines the absolute character index for the specified line
        and char using the *buffer's* code string."""
        #Make sure that we have chars and lines to work from
        if len(context.bufferstr) > 0 and len(self._chars) == 0:
            #Add one for the \n that we split on for each line
            self._chars = [ len(x) + 1 for x in context.bufferstr ]
            #Remove the last line break since it doesn't exist
            self._chars[-1] -= 1

            #Now we want to add up the number of characters in each line
            #as the lines progress so that it is easy to search for the
            #line of a single character index
            total = 0
            for i in range(len(self._chars)):
                total += self._chars[i]
                self._chars[i] = total
        
        return self._chars[line - 1] + char

class Operation(object):
    """Represents an insert, delete or replace operation for turning the cached
    version of a module into it's version in the emacs buffer.

    :attr context: the context of the buffer source code.
    :attr parser: the line parser instance to use for handling the operations.
    :attr mode: either 'replace', 'insert' or 'delete'.
    :attr icached: the [start,end] line index in the cached source code.
    :attr ibuffer: the [start,end] line index in the buffer source code.
    :attr index: the index of this operation in the module updater's operation list.
    :attr buffered: the list of complete fortran statements contained in the
      reference lists from the buffer source code.
    :attr cached: the list of complete fortran statements contained in the
      reference lists from the cached source code.
    :attr state: a tuple of (buffer index, cache index) for the statement that
      is currently being processed by the operation.
    """
    def __init__(self, context, parser, operation, index):
        self.context = context
        self.parser = parser
        self.docparser = self.context.parser.modulep.docparser
        self.mode = operation[0]
        self.icached = operation[1:3]
        self.ibuffer = operation[3:5]
        self.index = index

        self.buffered = self._get_buffered()
        self.cached = self._get_cached()

        #State holds a list of the current (buffer index, cache index) that is 
        #being parsed by the line parser.
        self.state = None
        self.bar_extent = False

        self._element = None
        self._docelement = None
        #This variable will hold a list of the lines that actually participated
        #in a real-time docstring update (buffer lines).
        self._doclines = None
        #We need to keep track of where we used to be so that we know
        #how to update the positions of the rest of instances in the module.
        self.length = 0
        #Doc delta tracks changes in length that were made to an element by a
        #docstring before of it's definition signature. Since the element was
        #already updated, we don't want it to get hit again by the module
        #updating it's child elements.
        self.docdelta = 0

    def __str__(self):
        if self.state is not None and self.state[0] is None:
            line, statement, charindex = self.curcached
        else:
            line, statement, charindex = self.curbuffer

        a = "{} ({},{})".format(self.mode, 
                                "-".join([str(c) for c in self.icached]),
                                "-".join([str(b) for b in self.ibuffer]))
        s = "{}({}): {}".format(line, charindex, statement)
        e = self.element.name
        return a + '\n' + s + "\n" + e

    def set_element(self, element):
        """Overrides the element instance that this operation will modify."""
        self._element = element

    @property
    def start(self):
        """Returns the line, column of the first statement in the cached 
        code that was affected."""
        return (self.icached[0], 0)

    @property
    def curbuffer(self):
        """Returns the current buffer statement for updating using the cached 
        execution state for the operation."""
        if self.state is not None:
            return self.buffered[self.state[0]]
        else:
            return (None, None, None)

    @property
    def curcached(self):
        """Returns the current cached statement for updating using the cached 
        execution state for the operation."""
        if self.state is not None:
            return self.cached[self.state[1]]
        else:
            return (None, None, None)

    @property
    def curlength(self):
        """Returns the character length of the statement currently being run."""
        if self.state[0] is None:
            return self.curcached[2]
        else:
            return self.curbuffer[2]

    @property
    def element(self):
        """Returns the instance of the element who owns the first line
        number for the operation in the cached source code."""
        #We assume here that the entire operation is associated with a single
        #code element. Since the sequence matcher groups operations by contiguous
        #lines of code to change, this is a safe assumption.
        if self._element is None:
            line = self.icached[0]
            #If we are inserting a new line, the location at the start of the line
            #that used to be there interferes with the element finder.
            if self.mode == "insert":
                line -= 1
            self._element = self.context.module.get_element(line, 0)

        return self._element

    @property
    def docelement(self):
        """Returns the instance of the element whose body owns the docstring
        in the current operation.
        """
        #This is needed since the decorating documentation
        #for types and executables is in the body of the module, but when they
        #get edited, the edit belongs to the type/executable because the character
        #falls within the absstart and end attributes.
        if self._docelement is None:
            if isinstance(self.element, Module):
                self._docelement = self.element
            else:
                ichar = self.element.module.charindex(self.icached[0], 1)
                if (ichar > self.element.docstart and ichar <= self.element.docend):
                    self._docelement = self.element.parent
                else:
                    self._docelement = self.element

        return self._docelement
        
    def handle(self):
        """Handles the real time update of some code from the cached representation
        of the module.
        """
        #If we have more statements in the buffer than the cached, it doesn't matter, 
        #we just run the first few replacements of the cache concurrently and do 
        #what's left over from the buffer.
        #REVIEW
        if self.mode == "insert": #then self.icached[0] == self.icached[1]:
            #We are inserting the lines from the buffer into the cached version
            for ib in range(len(self.buffered)):
                self.state = (ib, None)
                self.parser.parse(self)
                self._update_extent()

        elif self.mode == "delete": #then self.ibuffer[0] == self.ibuffer[1]:
            #We are deleting lines from the cached version
            for ic in range(len(self.cached)):
                self.state = (None, ic)
                self.parser.parse(self)
                self._update_extent()
        
        else: # mode == 'replace'
            #Need lines from both the buffer and the cached version
            #First we run all the statements in cached as deletions
            for ic in range(len(self.cached)):
                self.state = (None, ic)
                self.parser.parse(self, "delete")
                self._update_extent()
            #Then run all the buffer statements as insertions.
            for ib in range(len(self.buffered)):
                self.state = (ib, None)
                self.parser.parse(self, "insert")
                self._update_extent()

        self._handle_docstrings()
        
    def _handle_docstrings(self):
        """Searches through the lines affected by this operation to find
        blocks of adjacent docstrings to parse for the current element.
        """
        #Docstrings have to be continuous sets of lines that start with !!
        #When they change in any way (i.e. any of the three modes), we 
        #have to reparse the entire block because it has XML dependencies
        #Because of that, the cached version of the docstrings is actually
        #pointless and we only need to focus on the buffered.
        blocks = self._docstring_getblocks()

        if len(blocks) == 0:
            return

        xmldict = self._docstring_parse(blocks)
        delta = 0

        if isinstance(self.docelement, Module):
            delta += self.docparser.rt_update_module(xmldict, self.docelement)
        else:
            #We just need to handle the type and executable internal defs.
            if self.docelement.name in xmldict:
                docs = self.docparser.to_doc(xmldict[self.docelement.name][0],
                                             self.docelement.name)
                self.docparser.process_memberdocs(docs, self.docelement, False)
            #Also update the docstrings for any embedded types or executables.
            if isinstance(self.docelement, Executable):
                delta += self.docparser.process_embedded(xmldict, 
                                                         self.docelement, False)

        #Finally, we need to handle the overall character length change
        #that this update caused to the element first and then for the
        #operation as a whole for updating the module and its children.
        buffertot = sum([len(self.context.bufferstr[i]) for i in self._doclines])
        cachedtot = 0
    
        for i in range(self.icached[0],self.icached[1]):
            if self.docparser.RE_DOCS.match(self.context.cachedstr[i]):
                cachedtot += len(self.context.cachedstr[i])

        self.length = buffertot - cachedtot

        if delta == 0:
            #The update must have been to members variables of the module or the
            #executables/types. The element who owns the members is going to get
            #missed when the module updates its children.
            self.docelement.end += self.length
        else:
            #All the individual elements have been updated already, so just
            #set the length change for this operation.
            self.docdelta = delta
            
    def _docstring_parse(self, blocks):
        """Parses the XML from the specified blocks of docstrings."""
        result = {}
        for block, docline, doclength, key in blocks:
            doctext = "<doc>{}</doc>".format(" ".join(block))
            try:
                docs = ET.XML(doctext)
                docstart = self.parser.charindex(docline, 0, self.context)
                if not key in result:
                    result[key] = [list(docs), docstart, docstart + doclength]
                else:
                    #If there are docblocks separated by whitespace in the
                    #same element we can't easily keep track of the start and
                    #end character indices anymore.
                    result[key][0].extend(list(docs))
            except ET.ParseError:
                msg.warn(doctext)

        return result

    def _docstring_getblocks(self):
        """Gets the longest continuous block of docstrings from the buffer
        code string if any of those lines are docstring lines.
        """
        #If there are no lines to look at, we have nothing to do here.
        if self.ibuffer[0] == self.ibuffer[1]:
            return []

        lines = self.context.bufferstr[self.ibuffer[0]:self.ibuffer[1]]
        docblock = []
        result = []
        self._doclines = []

        #We need to keep track of the line number for the start of the
        #documentation strings.
        docline = 0
        doclength = 0

        first = self.docparser.RE_DOCS.match(lines[0])
        if first is not None:
            docblock.append(first.group("docstring"))
            docline = self.ibuffer[0]
            self._doclines.append(docline)
            doclength += len(lines[0]) + 1 # + 1 for \n removed by split.

            #We need to search backwards in the main buffer string for
            #additional tags to add to the block
            i = self.ibuffer[0] - 1
            while i > 0:
                current = self.context.bufferstr[i]
                docmatch = self.docparser.RE_DOCS.match(current)
                if docmatch is not None:
                    docblock.append(docmatch.group("docstring"))
                    docline = i
                    doclength += len(current) + 1
                else:
                    break
                i -= 1

        #Reverse the docblock list since we were going backwards and appending.
        if len(docblock) > 0:
            docblock.reverse()

        #Now handle the lines following the first line. Also handle the
        #possibility of multiple, separate blocks that are still valid XML.
        #We have to keep going until we have exceed the operational changes
        #or found the decorating element.
        i = self.ibuffer[0] + 1
        while (i < len(self.context.bufferstr) and 
               (i < self.ibuffer[1] or len(docblock) > 0)):
            line = self.context.bufferstr[i]
            docmatch = self.docparser.RE_DOCS.match(line)
            if docmatch is not None:
                docblock.append(docmatch.group("docstring"))
                doclength += len(line)
                if docline == 0:
                    docline = i
                #Only track actual documentation lines that are within the 
                #operations list of lines.
                if i < self.ibuffer[1]:
                    self._doclines.append(i)

            elif len(docblock) > 0:
                key = self._docstring_key(line)
                result.append((docblock, docline, doclength, key))
                docblock = []
                docline = 0
                doclength = 0

            #We need to exit the loop if we have exceeded the length of
            #the operational changes
            if len(docblock) == 0 and i > self.ibuffer[1]:
                break
            i += 1

        return result

    def _docstring_key(self, line):
        """Returns the key to use for the docblock immediately preceding
        the specified line."""
        decormatch = self.docparser.RE_DECOR.match(line)
        if decormatch is not None:
            key = "{}.{}".format(self.docelement.name, decormatch.group("name"))
        else:
            key = self.element.name

        return key

    def _update_extent(self):
        """Updates the extent of the element being altered by this operation
        to include the code that has changed."""
        #For new instances, their length is being updated by the module
        #updater and will include *all* statements, so we don't want to
        #keep changing the endpoints.
        if self.bar_extent:
            return
            
        original = self.element.end
        if self.mode == "insert":
            #The end value needs to increase by the length of the current
            #statement being executed.
            self.element.end += self.curlength
        elif self.mode == "delete":
            #Reduce end by statement length
            self.element.end -= self.curlength
        elif self.mode == "replace":
            #Check whether we are currently doing the delete portion of the
            #replacement or the insert portion.
            if self.state[0] is None:
                self.element.end -= self.curlength
            else:
                self.element.end += self.curlength

        #Keep track of the total effect of all statements in this operation
        #so that it is easy to update the module once they are all done.
        self.length += self.element.end - original
                
    def _get_buffered(self):
        """Gets a list of the statements that are new for the real time update."""
        lines = self.context.bufferstr[self.ibuffer[0]:self.ibuffer[1]]
        return self._get_statements(lines, self.ibuffer[0])

    def _get_cached(self):
        """Gets a list of statements that the operation will affect during the real
        time update."""
        lines = self.context.cachedstr[self.icached[0]:self.icached[1]]
        return self._get_statements(lines, self.icached[0])

    def _get_statements(self, lines, start):
        """Returns a list of complete Fortran statements for the specified lines by
        dealing with comments and line continuations. Returns a list of tuples of
        the form (linenum, statement, orginal character length).

        :arg start: the start index in the overall document for these lines.
        """
        result = []
        statement = []
        nocomment = [l.split("!")[0] for l in lines]
        length = 0

        for i in range(len(nocomment)):
            line = nocomment[i].strip()
            linenum = start + i
            length += len(lines[i]) + 1

            if len(line) == 0 or line[-1] != "&":
                statement.append(line)
                result.append((linenum-len(statement)+1, 
                               " ".join(statement), length))
                statement = []
                length = 0
            else:
                #Append the line without the continuation character.
                statement.append(line[:len(line)-1])

        return result        

class ModuleUpdater(object):
    """Updates the representations of the fortran modules in memory using
    new source code supplied from the emacs buffer."""
    def __init__(self):
        self.parser = LineParser(self)
        self.matcher = SequenceMatcher()
        self._operations = []
        self.unset = True

    def update(self, context):
        """Updates all references in the cached representation of the module
        with their latest code from the specified source code string that
        was extracted from the emacs buffer.

        :arg context: the buffer context information.
        """
        #Get a list of all the operations that need to be performed and then
        #execute them.
        self._operations = self._get_operations(context)

        for i in range(len(self._operations)):
            self._operations[i].handle()
            self.update_extent(self._operations[i])

        #Last of all, we update the string content of the cached version of
        #the module to have the latest source code.
        if len(self._operations) > 0:
            context.module.update_refstring(context.refstring)

    def update_extent(self, operation):
        """Updates the extent of the *module* and its elements using the
        specified operation which has already been executed."""
        #Operations represent continuous lines of code.
        #The last operation to be executed is still the current one.
        line, col = operation.start
        operation.context.module.update_elements(line, col, operation.length,
                                                 operation.docdelta)

    def _get_operations(self, context):
        """Returns a list of operations that need to be performed to turn the
        cached source code into the one in the buffer."""
        #Most of the time, the real-time update is going to fire with
        #incomplete statements that don't result in any changes being made
        #to the module instances. The SequenceMatches caches hashes for the
        #second argument. Logically, we want to turn the cached version into
        #the buffer version; however, the buffer is the string that keeps 
        #changing.

        #in order to optimize the real-time update, we *switch* the two strings
        #when we set the sequence and then fix the references on the operations
        #after the fact.
        if context.module.changed or self.unset:
            self.matcher.set_seq2(context.cachedstr)
            self.unset = False
            #Set the changed flag back to false now that the sequencer has
            #reloaded it.
            context.module.changed = False
        self.matcher.set_seq1(context.bufferstr)

        opcodes = self.matcher.get_opcodes()
        result = []

        #Index i keeps track of how many operations were actually added because
        #they constituted a change we need to take care of.
        i = 0
        for code in opcodes:
            if code[0] != "equal":
                #Replacements don't have a mode change. All of the operations
                #switch the order of the line indices for the two strings.
                if code[0] == "insert":
                    newcode = ("delete", code[3], code[4], code[1], code[2])
                elif code[0] == "delete":
                    newcode = ("insert", code[3], code[4], code[1], code[2])
                else:
                    newcode = ("replace", code[3], code[4], code[1], code[2])

                op = Operation(context, self.parser, newcode, i)
                result.append(op)
                i += 1

        return result

    def update_instance_extent(self, instance, module, operation):
        """Updates a new instance that was added to a module to be complete
        if the end token is present in any remaining, overlapping operations.
        """
        #Essentially, we want to look in the rest of the statements that are
        #part of the current operation to see how many more of them pertain 
        #to the new instance that was added.

        #New signatures only result in instances being added if mode is "insert"
        #or "replace". In both cases, the important code is in the buffered
        #statements, *not* the cached version. Iterate the remaining statements
        #in the buffer and look for the end_token for the instance. If we don't
        #find it, check for overlap between the operations' index specifiers.
        instance.end -= operation.curlength
        end_token = instance.end_token
        (ibuffer, length) = self._find_end_token(end_token, operation)
        cum_length = length
        
        opstack = [operation]
        while ibuffer is None and opstack[-1].index + 1 < len(self._operations):
            #We didn't find a natural termination to the new instance. Look for
            #overlap in the operations
            noperation = self._operations[opstack[-1].index + 1]
            #We only want to check the next operation if it is a neighbor
            #in line numbers in the buffer.
            if noperation.ibuffer[0] - opstack[-1].ibuffer[1] == 1:
                (ibuffer, length) = self._find_end_token(end_token, noperation)
                cum_length += length
                opstack.append(noperation)
            else:
                break
            
        if ibuffer is not None:
            instance.incomplete = False
            instance.end += cum_length
            for op in opstack:
                op.bar_extent = True
                op.set_element(instance)
        else:
            #We set the element for the current operation to be the new instance
            #for the rest of statements in its set.
            operation.set_element(instance)
            
    def _find_end_token(self, end_token, operation):
        """Looks for a statement in the operation's list that matches the specified
        end token. Returns the index of the statement in the operation that matches.
        """
        ibuffer, icache = operation.state
        length = operation.buffered[ibuffer][2]
        result = None

        for i in range(len(operation.buffered) - ibuffer - 1):
            linenum, statement, charlength = operation.buffered[i + ibuffer + 1]
            length += charlength
            
            if end_token in statement.lower():
                #We already have the absolute char index for the start of the
                #instance; we just need to update the end.
                result = (i + ibuffer + 1, length)
                break        

        #If we didn't find a terminating statement, the full length of the operation
        #is a good estimate for the extent of the instance.
        if result is None:
            result = (None, length)

        return result
