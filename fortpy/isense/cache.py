from fortpy.code import CodeParser
from . import builtin
import time
import re
import pyparsing
from . import rtupdate

#Time cache has expiration on the items inside that get used only
#temporarily during code completion
_time_caches = []
#This is our instance of the CodeParser. It handles the parsing
#of all the Fortran modules and has its own caching builtin.
_parsers = { "default": CodeParser() }
nester = pyparsing.nestedExpr("(",")")
#Get a generic module updater for doing real-time updates on
#source code sent from the emacs buffer.

def parser(key = "default"):
    """Returns the parser for the given key, (e.g. 'ssh')"""
    #Make sure we have a parser for that key. If we don't, then set
    #one up if we know what parameters to use; otherwise return the
    #default parser.
    if key not in _parsers:
        if key == "ssh":
            _parsers["ssh"] = CodeParser(True, False)
        else:
            key = "default"

    return _parsers[key]

def clear_caches(delete_all=False):
    """Fortpy caches many things, that should be completed after each completion
    finishes.

    :param delete_all: Deletes also the cache that is normally not deleted,
        like parser cache, which is important for faster parsing.
    """
    global _time_caches

    if delete_all:
        _time_caches = []
        _parser = { "default": CodeParser() }
    else:
        # normally just kill the expired entries, not all
        for tc in _time_caches:
            # check time_cache for expired entries
            for key, (t, value) in list(tc.items()):
                if t < time.time():
                    # delete expired entries
                    del tc[key]

def time_cache(time_add_setting):
    """ This decorator works as follows: Call it with a setting and after that
    use the function with a callable that returns the key.
    But: This function is only called if the key is not available. After a
    certain amount of time (`time_add_setting`) the cache is invalid.
    """
    def _temp(key_func):
        dct = {}
        _time_caches.append(dct)

        def wrapper(optional_callable, *args, **kwargs):
            key = key_func(*args, **kwargs)
            value = None
            if key in dct:
                expiry, value = dct[key]
                if expiry > time.time():
                    return value
            value = optional_callable()
            time_add = getattr(settings, time_add_setting)
            if key is not None:
                dct[key] = time.time() + time_add, value
            return value
        return wrapper
    return _temp

@time_cache("call_signatures_validity")
def cache_call_signatures(source, user_pos, stmt):
    """This function calculates the cache key."""
    index = user_pos[0] - 1
    lines = source.splitlines() or ['']
    if source and source[-1] == '\n':
        lines.append('')

    before_cursor = lines[index][:user_pos[1]]
    other_lines = lines[stmt.start_pos[0]:index]
    whole = '\n'.join(other_lines + [before_cursor])
    before_bracket = re.match(r'.*\(', whole, re.DOTALL)

    module_path = stmt.get_parent_until().path
    return None if module_path is None else (module_path, before_bracket, stmt.start_pos)

rt = rtupdate.ModuleUpdater()

#Setup and compile a bunch of regular expressions that we use
#every time the user context needs to be determined.
RE_COMMENTS = re.compile("\n\s*!")
_RX_MODULE = r"(\n|^)\s*module\s+(?P<name>[a-z0-9_]+)"
RE_MODULE = re.compile(_RX_MODULE, re.I)

_RX_TYPE = (r"\s+type(?P<modifiers>,\s+(public|private))?(\s+::)?"
                         "\s+(?P<name>[A-Za-z0-9_]+)")
RE_TYPE = re.compile(_RX_TYPE)
_RX_EXEC = r"\n\s*((?P<type>character|real|type|logical|integer)?" + \
                r"(?P<kind>\([a-z0-9_]+\))?)?(,(?P<modifiers>[^\n]+?))?\s*" + \
                r"(?P<codetype>subroutine|function)\s+(?P<name>[^(]+)" + \
                r"\s*\((?P<parameters>[^)]+)\)"
RE_EXEC = re.compile(_RX_EXEC, re.I)
RE_MEMBERS = parser().modulep.vparser.RE_MEMBERS
_RX_DEPEND = r"^\s*(?P<sub>call\s+)?(?P<exec>[a-z0-9_%]+)\s*(?P<args>\([^\n]+)$"
RE_DEPEND = re.compile(_RX_DEPEND, re.M | re. I)
_RX_CURSOR = r"(?P<symbol>[^\s=,%(]+)"
RE_CURSOR = re.compile(_RX_CURSOR, re.I)
_RX_FULL_CURSOR = r"(?P<symbol>[^\s=,(]+)"
RE_FULL_CURSOR = re.compile(_RX_FULL_CURSOR, re.I)

#A list of all the function names that are 'builtin' in Fortran
builtin = builtin.load(_parsers["default"].modulep.docparser, _parsers["default"].serialize)

_common_builtin = [ "OPEN", "CLOSE", "PRESENT", "ABS", "DIM", "ALLOCATE", "SQRT"]
common_builtin = {}
for fun in _common_builtin:
    common_builtin[fun.lower()] = fun
