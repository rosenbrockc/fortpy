"""This module is copied into package directories that are created by ftypes
for interoperability with Fortran code. It has a single public class that 
makes it more convenient to interact with the results of the ctypes calls. It
also provides a method to find the static type for a given symbol name.
"""
import weakref
from numpy import require, array, asanyarray, dtype as _dtype
from numpy.ctypeslib import ndpointer, as_array as farray
from ctypes import (c_int, c_double, POINTER, c_bool, c_float,
                    c_short, c_long, byref, c_void_p, addressof)
c_int_p = POINTER(c_int)
c_double_p = POINTER(c_double)
c_float_p = POINTER(c_float)
c_bool_p = POINTER(c_bool)
c_short_p = POINTER(c_short)
c_long_p = POINTER(c_long)

libs = {}
"""The dict of shared libraries loaded so far by this module. Indexed by the
lower-case path to the *.so."""
compilers = {}
"""Dict of compilers found for each shared library. Indexed by lower-case path
to the shared library. Value is either 'gfortran' or 'ifort'.
"""
fresults = {}
"""Dict of FtypeResult instances indexed by the 'module.executable' key of the
wrapper subroutine that generated the output. Used to manage shared memory
being overwritten from repeated calls to the same executables.
"""

def as_array(pointer, shape):
    """Returns a contiguous (c-ordered) array for the specified pointer
    returned from a fortran allocatable array.

    :arg shape: a tuple of the integer size of each dimension.
    """
    return farray(pointer, shape).T

def static_symbol(module, method, lib, bindc=False):
    """Returns the symbol for the specified *fortran* module and subroutine/function
    that has been compiled into a shared library.

    :arg lib: the full path to the shared library *.so file.
    """
    from os import path
    from numpy.ctypeslib import load_library
    libpath = path.abspath(lib)
    libkey = libpath.lower()

    if libkey not in compilers:
        compilers[libkey] = detect_compiler(libpath)
        if compilers[libkey] == False:
            raise ValueError("Couldn't auto-detect the compiler for {}".format(libpath))
    compiler = compilers[libkey]

    if libkey not in libs:
        libs[libkey] = load_library(libpath, "")
        if not libs[libkey]:
            raise ValueError("Couldn't auto-load the shared library with ctypes.")
    slib = libs[libkey]
    
    identifier = "{}.{}".format(module, method)
    if bindc:
        symbol = method
    elif compiler == "gfortran":
        symbol = "__" + identifier.lower().replace(".", "_MOD_")
    else:
        #We assume the only other option is ifort.
        symbol = identifier.lower().replace(".", "_mp_") + "_"

    if hasattr(slib, symbol):
        return getattr(slib, symbol)
    else:
        return None

def detect_compiler(libpath):
    """Determines the compiler used to compile the specified shared library by
    using the system utilities.

    :arg libpath: the full path to the shared library *.so file.
    """
    from os import waitpid, path
    from subprocess import Popen, PIPE
    command = "nm {0}".format(path.abspath(libpath))
    child = Popen(command, shell=True, executable="/bin/bash", stdout=PIPE)
    # Need to do this so that we are sure the process is done before moving on
    waitpid(child.pid, 0)
    contents = child.stdout.readlines()
    i = 0
    found = False
    while i < len(contents) and found == False:
        if "_MOD_" in contents[i]:
            found = "gfortran"
        elif "_mp_" in contents[i]:
            found = "ifort"
        i += 1

    return found

class Ftype(object):
    """Represents an output data set from python for interacting with a fortran shared
    library.
    """
    def __init__(self, pointer, indices, libpath):
        """
        :arg pointer: the c-pointer to the memory address.
        :arg indices: the integer size of each dimension in the pointer array.
        """
        self.pointer = pointer
        self.indices = indices
        self.libpath = libpath
        self.deallocated = False
        """Specifies whether this c-pointer has already deallocated the Fortran
        memory that it is referencing."""
 
    def clean(self):
        """Deallocates the fortran-managed memory that this ctype references.
        """
        if not self.deallocated:
            #Release/deallocate the pointer in fortran.
            method = self._deallocator()
            if method is not None:
                dealloc = static_symbol("ftypes_dealloc", method, self.libpath, True)
                if dealloc is None:
                    return
                arrtype = ndpointer(dtype=int, ndim=1, shape=(len(self.indices),), flags="F")
                dealloc.argtypes = [c_void_p, c_int_p, arrtype]
                nindices = require(array([i.value for i in self.indices]), int, "F")
                dealloc(byref(self.pointer), c_int(len(self.indices)), nindices)

        self.deallocated = True                
            
    def _deallocator(self):
        """Returns the name of the subroutine in ftypes_dealloc.f90 that can
        deallocate the array for this Ftype's pointer.

        :arg ctype: the string c-type of the variable.
        """
        lookup = {
            "c_bool": "logical",
            "c_double": "double",
            "c_double_complex": "complex",
            "c_char": "char",
            "c_int": "int",
            "c_float": "float",
            "c_short": "short",
            "c_long": "long"            
        }
        ctype = type(self.pointer).__name__.replace("LP_", "").lower()
        if ctype in lookup:
            return "dealloc_{0}_{1:d}d".format(lookup[ctype], len(self.indices))
        else:
            return None
        
class FtypesResult(object):
    """Represents the result from executing a fortran subroutine or function
    using the ctypes interface.
    """
    def __init__(self, module, name, utype):
        """
        :arg name: the name of the subroutine/function in the *original* code.
        :arg utype: the underlying type of the executable. Either 'subroutine'
          or 'function'. Affects the behavior of calling an instance.
        :arg result: a dictionary of the parameter values that were of intent
          "out" or "inout" and their values on exit.
        """
        self.identifier = "{}.{}".format(module, name).lower()
        self.utype = utype
        self.result = {}
        self._finalizers = {}
        """A dictionary of Ftype instances with c-pointer information for finalizing
        the fortran arrays.
        """
        
        def on_die(kref):
            for key in list(self._finalizers.keys()):
                self._finalizers[key].clean()            
        self._del_ref = weakref.ref(self, on_die)
        
        if self.identifier in fresults:
            #Cleanup the previous result obtained from this method (this involves
            #reallocating and copying the values to a new array). Change the active
            #result (i.e. with active pointers to fortran-managed memory) to
            #be this new one.
            fresults[self.identifier].cleanup()
        fresults[self.identifier] = self

    def cleanup(self):
        """Cleans up this result so that all its pointers reference memory controlled
        *outside* of the shared library loaded with ctypes.
        """
        #First we *copy* the arrays that we currently have pointers to. This is not
        #the optimal solution; however, because of limitations in ctypes, we don't
        #know anything better at the moment.
        for key in self.result:
            self.result[key] = self.result[key].copy()

        #Now deallocate the pointers managed by Fortran in the shared library so that
        #any subsequent calls to the executable that generated this result creates
        #new instances in memory.
        for key in list(self._finalizers.keys()):
            self._finalizers[key].clean()
        
    def add(self, varname, result, pointer=None):
        """Adds the specified python-typed result and an optional Ftype pointer
        to use when cleaning up this object.

        :arg result: a python-typed representation of the result.
        :arg pointer: an instance of Ftype with pointer information for deallocating
          the c-pointer.
        """
        self.result[varname] = result
        setattr(self, varname, result)
        if pointer is not None:
            self._finalizers[varname] = pointer
        
    def __call__(self, key=None):
        """Returns the value of the result from calling the subroutine/function
        for the given key.
        """
        if self.utype == "function":
            return self.result["{}_o"]
        elif key is not None:
            return self.result[key]
        elif len(self.result) == 1:
            return list(self.result.values())[0]
        else:
            raise ValueError("Can't call result without a key for subroutines.")
