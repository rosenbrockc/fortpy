"""
This module contains variables with global |fortpy| settings. To change the
behavior of |fortpy|, change the variables defined in :mod:`fortpy.settings`.

Plugins should expose an interface so that the user can adjust the
configuration.

Example usage::

    from fortpy import settings
    settings.case_insensitive_completion = True

Completion output
~~~~~~~~~~~~~~~~~

.. autodata:: case_insensitive_completion
.. autodata:: add_bracket_after_function
.. autodata:: no_completion_duplicates

Filesystem cache
~~~~~~~~~~~~~~~~

.. autodata:: cache_directory
.. autodata:: use_filesystem_cache

Parser
~~~~~~

.. autodata:: real_time_update

Dynamic stuff
~~~~~~~~~~~~~

.. autodata:: dynamic_params

Recursions
~~~~~~~~~~

First off there is a global limit :data:`max_executions`. This limit 
is important, to set a maximum amount of time, the completion may use.

The default values are based on experiments while completing the |fortpy| library
itself (inception!). These settings make the completion
definitely worse in some cases. But a completion should also be fast.

.. autodata:: max_until_execution_unique
.. autodata:: max_function_recursion_level
.. autodata:: max_executions
.. autodata:: scale_call_signatures

Caching
~~~~~~~

.. autodata:: call_signatures_validity

Testing
~~~~~~~

.. autodata:: unit_testing_mode
.. autodata:: use_test_cache

"""
import os
import platform

# ----------------
# completion output settings
# ----------------

case_insensitive_completion = True
"""
The completion is by default case insensitive.
"""

add_bracket_after_function = False
"""
Adds an opening bracket after a function, because that's normal behaviour.
Removed it again, because in VIM that is not very practical.
"""

no_completion_duplicates = True
"""
If set, completions with the same name don't appear in the output anymore,
but are in the `same_name_completions` attribute.
"""

# ----------------
# Filesystem cache
# ----------------

use_filesystem_cache = True
"""
Use filesystem cache to save once parsed files with pickle.
"""

if platform.system().lower() == 'windows':
    _cache_directory = os.path.join(os.getenv('APPDATA') or '~', 'Fortpy',
                                    'Fortpy')
elif platform.system().lower() == 'darwin':
    _cache_directory = os.path.join('~', 'Library', 'Caches', 'Fortpy')
    _temp_directory = os.path.join('~', 'var', 'tmp', 'Fortpy')
else:
    _cache_directory = os.path.join(os.getenv('XDG_CACHE_HOME') or '~/.cache',
                                    'fortpy')
cache_directory = os.path.expanduser(_cache_directory)
"""
The path where all the caches can be found.

On Linux, this defaults to ``~/.cache/fortpy/``, on OS X to
``~/Library/Caches/Fortpy/`` and on Windows to ``%APPDATA%\\Fortpy\\Fortpy\\``.
On Linux, if environment variable ``$XDG_CACHE_HOME`` is set,
``$XDG_CACHE_HOME/fortpy`` is used instead of the default one.
"""

# ----------------
# parser
# ----------------

real_time_update = False
"""
Use the fast parser. This means that reparsing is only being done if
something has been changed e.g. to a function. If this happens, only the
function is being reparsed.
"""

# ----------------
# dynamic stuff
# -----
dynamic_params = True
"""
A dynamic param completion, finds the callees of the function, which define
the params of a function.
"""

# ----------------
# recursions
# ----------------

max_until_execution_unique = 50
"""
This limit is probably the most important one, because if this limit is
exceeded, functions can only be one time executed. So new functions will be
executed, complex recursions with the same functions again and again, are
ignored.
"""

max_function_recursion_level = 5
"""
`max_function_recursion_level` is more about whether the recursions are
stopped in deepth or in width. The ratio beetween this and
`max_until_execution_unique` is important here. It stops a recursion (after
the number of function calls in the recursion), if it was already used
earlier.
"""

max_executions = 250
"""
A maximum amount of time, the completion may use.
"""

scale_call_signatures = 0.1
"""
Because call_signatures is normally used on every single key hit, it has
to be faster than a normal completion. This is the factor that is used to
scale `max_executions` and `max_until_execution_unique`:
"""

# ----------------
# caching validity (time)
# ----------------
call_signatures_validity = 3.0
"""
Finding function calls might be slow (0.1-0.5s). This is not acceptible for
normal writing. Therefore cache it for a short time.
"""

unit_testing_mode = False
"""Specifies that the parsers should run in unit test mode, which means that
the config.xml file is not loaded by the code parser. It makes the parsers
function more like they would on a regular system that hasn't been configured yet.
"""

use_test_cache = False
"""Use the cache in Fortpy_Testing so that we don't have competition with the
live caches during development when the serializer version may be incremented.
"""
