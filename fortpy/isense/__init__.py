from . import cache
import fortpy.debug as debug
from . import context
from . import evaluator
import re
import os
import fortpy.settings as settings
from .classes import BaseDefinition

class Script(object):
    """
    A Script is the base for completions, goto or whatever you want to do with
    |fortpy|. This is modified from the original |jedi| Script class.

    You can either use the ``source`` parameter or ``path`` to read a file.
    Usually you're going to want to use both of them (in an editor).

    :param source: The source code of the current file, separated by newlines.
    :type source: str
    :param line: The line to perform actions on (starting with 1).
    :type line: int
    :param col: The column of the cursor (starting with 0).
    :type col: int
    :param path: The path of the file in the file system, or ``''`` if
        it hasn't been saved yet.
    :type path: str or None
    """
    def __init__(self, source=None, line=None, column=None, path=None):
        self._orig_path = path
        self.path = None if path is None else os.path.abspath(path)

        #Here we need to add support for tramp over SSH for the file that
        #is being edited. We will have cached versions of the other modules
        #in the global code cache, but we need to update the file being
        #edited in that cache so that any other buffers with reference to
        #this module being edited have the latest information.
        if source is None:
            source = cache.parser().tramp.read(path)

        #This section makes sure that the argument we have from the buffer
        #for the line number actually makes sense for the source code we have.
        lines = source.splitlines() or ['']
        if source and source[-1] == '\n':
            lines.append('')
        line = max(len(lines), 1) if line is None else line
        if not (0 < line <= len(lines)):
            raise ValueError('`line` parameter is not in a valid range.')

        #Same thing as above, except we are checking validity of column.
        line_len = len(lines[line - 1])
        column = line_len if column is None else column
        if not (0 <= column <= line_len):
            raise ValueError('`column` parameter is not in a valid range.')
        self._pos = line-1, column

        #Clear the time dependent caches. These are things we cached related
        #to completion but *not* to parsing. This is because multiple calls
        #can be made for completion of a symbol while on the same line of code.
        cache.clear_caches()
        #This is the timing module for the parsing etc. to find bottlenecks
        #and help optimize the code. It sets the timer to zero.
        debug.reset_time()

        #This is the actual parsing of the source file using the fortpy parsers
        #NOTE: to avoid errors in pre-optimization, I will leave this like so
        #We can use some cache control on the UserContext later to speed things up. TODO
        if cache.parser().tramp.is_ssh(path):
            parser_key = "ssh"
        else:
            parser_key = "default"
        self._user_context = context.UserContext(source, self._pos,
                                                 self._orig_path, parser_key)
        self._evaluator = evaluator.Evaluator(self._user_context, self._pos)

        #Time how long all that parsing took.
        debug.speed('init')

    @property
    def context(self):
        """Gets the current cursor context from the file."""
        return self._user_context

    @property
    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, repr(self._orig_path))

    def bracket_complete(self):
        """Returns a function call signature for completion whenever
        a bracket '(' is pressed."""
        return self._evaluator.bracket_complete()

    def in_function_call(self):
        """Return either completion information or a call signature for
        the function definition that we are on currently."""
        #The last thing to do before we can form completions etc. is perform
        #a real-time update of the in-memory versions of the modules.
        if settings.real_time_update:
            cache.rt.update(self._user_context)

        return self._evaluator.in_function_call()

    def completions(self):
        """
        Return :class:`classes.Completion` objects. Those objects contain
        information about the completions, more than just names.

        :return: Completion objects, sorted by name.
        :rtype: list of :class:`classes.Completion`
        """
        debug.speed('completions start')
        comps = self._evaluator.complete()
        debug.speed('completions end')

        return sorted(comps, key=lambda x: (x.name.lower()))

    def goto_definitions(self):
        """
        Return the definition of a the symbol under the cursor via exact match.
        Goes to that definition with a buffer.
        """
        element = self._evaluator.get_definition()
        if element is not None:
            return BaseDefinition(self._user_context, element)
        else:
            return None

    def usages(self, additional_module_paths=()):
        """
        Return :class:`classes.Definition` objects, which contain all
        names that point to the definition of the name under the cursor. This
        is very useful for refactoring (renaming), or to show all usages of a
        variable.

        .. todo:: Implement additional_module_paths

        :rtype: list of :class:`classes.Definition`
        """

    def call_signatures(self):
        """
        Return the function object of the call you're currently in.

        E.g. if the cursor is here::

            abs(# <-- cursor is here

        This would return the ``abs`` function. On the other hand::

            abs()# <-- cursor is here

        This would return ``None``.

        :rtype: list of :class:`classes.CallSignature`
        """
    
