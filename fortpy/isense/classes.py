from . import cache
import re
from fortpy.elements import Function, Subroutine, ValueElement, Executable, TypeExecutable
from fortpy.printing import docs

class BaseDefinition(object):
    def __init__(self, context, element):
        self._context = context     
        self._element = element

    @property
    def name(self):
        """Returns the name of the underlying code element."""
        return self._element.name

    @property
    def type(self):
        """
        The type of the definition. One of the classes in elements.py.
        """
        return type(self._element).__name__

    @property
    def module_path(self):
        """The module path."""
        return self._element.module.filepath

    @property
    def module_name(self):
        """
        The module name
        """
        return str(self._element.module.name)

    @property
    def line(self):
        """The line where the definition occurs (starting with 0)."""
        return self._element.module.linenum(self._element.start)

    @property
    def column(self):
        """The column where the definition occurs (starting with 0)."""
        line = self.line
        return self._element.start - self._element.module.charindex(line, 0)

    @property
    def docstring(self):
        """Returns the <summary> tag of the code element represented by this
        completion.
        """
        if hasattr(self._element, "summary"):
            return self._element.summary
        else:
            return "Undocumented Symbol"

    @property
    def fulldoc(self):
        """Returns the full docstring for the code element."""
        return docs.format(self._element)

    @property
    def description(self):
        """Returns the full docstring information for the element suggested 
        as a completion."""
        result = ""
        if isinstance(self._element, ValueElement):
            if self._element.kind is not None:
                result = "{}({}) | {}".format(self._element.dtype, self._element.kind,
                                          self._element.summary)
            else:
                result = "{} | {}".format(self._element.dtype,
                                          self._element.summary)
        elif isinstance(self._element, Executable):
            result = "({})".format(self._element.parameters_as_string())
        elif isinstance(self._element, str):
            result = "Intrinsic Fortran Symbol"
        elif isinstance(self._element, TypeExecutable):
            result = self._type_description()

        #Clean off any line breaks from the XML and excessive whitespace.
        cleaned = re.sub("\s+", " ", result.replace("\n", " "))
        return cleaned

    @property
    def full_name(self):
        """Dot-separated path of this object."""
        return self._context.element.full_name

    @property
    def params(self):
        """
        Raises an ``AttributeError``if the definition is not callable.
        Otherwise returns a list of `ValueElement` that represents the params.
        """
        if self.context.el_type in [Function, Subroutine]:
            return self.evaluator.element.parameters

    def parent(self):
        return self.evaluator.element.parent

    def __repr__(self):
        return "<%s %s>" % (type(self).__name__, self.description)

class Completion(BaseDefinition):
    """
    `Completion` objects are returned from :meth:`isense.Script.completions`. They
    provide additional information about a completion.
    """
    def __init__(self, context, element, like_name_length):
        super(Completion, self).__init__(context, element)
        self._like_name_length = like_name_length

    def _complete(self, like_name):
        append = ''
        if settings.add_bracket_after_function \
                and self.context.el_type in [Function, Subroutine]:
            append = '('

        #Get the first name in the evaluator list since it is closest in
        #context to the cursor.
        if like_name:
            name = self._name[self._like_name_length:]
        return name + append

    @property
    def complete(self):
        """
        Return the rest of the word, e.g. completing ``isinstance``::

            subrou# <-- Cursor is here

        would return the string 'tine'. It also adds additional stuff, depending
        on your `settings.py`.
        """
        return self._complete(True)

    @property
    def name(self):
        """
        Similar to :attr:`complete`, but return the whole word, for
        example::

            subrout

        would return `subroutine`.
        """
        if hasattr(self._element, "name"):
            return self._element.name
        else:
            return str(self._element)

    def __repr__(self):
        return '<%s: %s>' % (type(self).__name__, self.name)

    def _type_description(self):
        """Gets the completion description for a TypeExecutable."""
        #This is a little tricker because the docstring is housed
        #inside of the module that contains the actual executable.
        #These TypeExecutables are just pointers.
        iexec = self._element.target
        if iexec is not None:
            result = "method() | " + iexec.summary
        else:
            result = "Type Method: points to executable in module."    
        return result
        
    def follow_definition(self):
        """g
        Return the original definitions. I strongly recommend not using it for
        your completions, because it might slow down |fortpy|. If you want to
        read only a few objects (<=20), it might be useful, especially to get
        the original docstrings.
        """
        return "original definition pending"
