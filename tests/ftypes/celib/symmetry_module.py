"""Auto-generated python module for interaction with Fortran shared library
via ctypes. Generated for module symmetry_module.
"""
from ctypes import *
from ftypes import *
from numpy.ctypeslib import load_library, ndpointer
from numpy import require
from os import path

def make_primitive(avecs, atomtype, atom_pos, lattcoords, eps_):
    """No XML documentation summary.
    """
    libpath = path.join(path.dirname(__file__), "ftypes.celib.so")
    method = static_symbol("symmetry_module_c", "make_primitive_c", libpath, True)
    method.argtypes = [ndpointer(dtype=float, ndim=2, shape=(3, 3), flags="F, writeable"), 
                       ndpointer(dtype=int, ndim=1, flags="F"), c_int_p, c_void_p, 
                       ndpointer(dtype=float, ndim=2, flags="F"), c_int_p, c_int_p, 
                       c_void_p, c_bool_p, c_double_p, c_void_p, c_int_p]

    avecs_a = require(avecs, float, "F, writeable")
    atomtype_a = require(atomtype, int, "F")
    atomtype_0 = c_int(atomtype_a.shape[0])
    atomtype_o = POINTER(c_int)()

    atom_pos_a = require(atom_pos, float, "F")
    atom_pos_0 = c_int(atom_pos_a.shape[0])
    atom_pos_1 = c_int(atom_pos_a.shape[1])
    atom_pos_o = POINTER(c_double)()

    lattcoords_c = c_bool(lattcoords)
    eps__c = c_double(eps_)
    removed__o = POINTER(c_int)()
    removed__0 = c_int(0)

    method(avecs_a, atomtype_a, byref(atomtype_0), byref(atomtype_o), atom_pos_a, 
           byref(atom_pos_0), byref(atom_pos_1), byref(atom_pos_o), byref(lattcoords_c), 
           byref(eps__c), byref(removed__o), byref(removed__0))

    result = FtypesResult("symmetry_module", "make_primitive", "Subroutine")
    result.add("avecs", avecs_a)
    atomtype_ft = Ftype(atomtype_o, [atomtype_0], libpath)
    result.add("atomtype", as_array(atomtype_o, (atomtype_0.value,)), atomtype_ft)
    atom_pos_ft = Ftype(atom_pos_o, [atom_pos_0, atom_pos_1], libpath)
    result.add("atom_pos", as_array(atom_pos_o, (atom_pos_1.value, atom_pos_0.value)), atom_pos_ft)
    if not all([removed__0.value==0]):
        removed__ft = Ftype(removed__o, [removed__0], libpath)
        result.add("removed_", as_array(removed__o, (removed__0.value,)), removed__ft)
    else:
        result.add("removed_", None)

    return result

def does_mapping_exist(v, this_type, atom_pos, atomtype, eps):
    """No XML documentation summary.
    """
    libpath = path.join(path.dirname(__file__), "ftypes.celib.so")
    method = static_symbol("symmetry_module_c", "does_mapping_exist_c", libpath, True)
    method.argtypes = [ndpointer(dtype=float, ndim=1, shape=(3,), flags="F"), c_int_p, 
                       ndpointer(dtype=float, ndim=2, flags="F"), c_int_p, c_int_p, 
                       ndpointer(dtype=int, ndim=1, flags="F"), c_int_p, c_bool_p, 
                       c_double_p]

    v_a = require(v, float, "F")
    this_type_c = c_int(this_type)
    atom_pos_a = require(atom_pos, float, "F")
    atom_pos_0 = c_int(atom_pos_a.shape[0])
    atom_pos_1 = c_int(atom_pos_a.shape[1])
    atomtype_a = require(atomtype, int, "F")
    atomtype_0 = c_int(atomtype_a.shape[0])
    mapped_c = c_bool(False)
    eps_c = c_double(eps)
    method(v_a, byref(this_type_c), atom_pos_a, byref(atom_pos_0), byref(atom_pos_1), 
           atomtype_a, byref(atomtype_0), byref(mapped_c), byref(eps_c))

    result = FtypesResult("symmetry_module", "does_mapping_exist", "Subroutine")
    result.add("mapped", mapped_c.value)

    return result

def bring_into_cell(v, cart_to_latt, latt_to_cart, eps):
    """No XML documentation summary.
    """
    libpath = path.join(path.dirname(__file__), "ftypes.celib.so")
    method = static_symbol("symmetry_module_c", "bring_into_cell_c", libpath, True)
    method.argtypes = [ndpointer(dtype=float, ndim=1, shape=(3,), flags="F, writeable"), 
                       ndpointer(dtype=float, ndim=2, shape=(3, 3), flags="F"), 
                       ndpointer(dtype=float, ndim=2, shape=(3, 3), flags="F"), c_double_p]

    v_a = require(v, float, "F, writeable")
    cart_to_latt_a = require(cart_to_latt, float, "F")
    latt_to_cart_a = require(latt_to_cart, float, "F")
    eps_c = c_double(eps)
    method(v_a, cart_to_latt_a, latt_to_cart_a, byref(eps_c))

    result = FtypesResult("symmetry_module", "bring_into_cell", "Subroutine")
    result.add("v", v_a)

    return result

def get_spacegroup(avecs, atomtype, input_pos, lattcoords, eps_):
    """No XML documentation summary.
    """
    libpath = path.join(path.dirname(__file__), "ftypes.celib.so")
    method = static_symbol("symmetry_module_c", "get_spacegroup_c", libpath, True)
    method.argtypes = [ndpointer(dtype=float, ndim=2, shape=(3, 3), flags="F"), 
                       ndpointer(dtype=int, ndim=1, flags="F"), c_int_p, 
                       ndpointer(dtype=float, ndim=2, flags="F"), c_int_p, c_int_p, 
                       c_void_p, c_int_p, c_int_p, c_int_p, c_void_p, c_int_p, c_int_p, 
                       c_bool_p, c_double_p]

    avecs_a = require(avecs, float, "F")
    atomtype_a = require(atomtype, int, "F")
    atomtype_0 = c_int(atomtype_a.shape[0])
    input_pos_a = require(input_pos, float, "F")
    input_pos_0 = c_int(input_pos_a.shape[0])
    input_pos_1 = c_int(input_pos_a.shape[1])
    sg_op_o = POINTER(c_double)()
    sg_op_0 = c_int(0)
    sg_op_1 = c_int(0)
    sg_op_2 = c_int(0)

    sg_fract_o = POINTER(c_double)()
    sg_fract_0 = c_int(0)
    sg_fract_1 = c_int(0)

    lattcoords_c = c_bool(lattcoords)
    eps__c = c_double(eps_)
    method(avecs_a, atomtype_a, byref(atomtype_0), input_pos_a, byref(input_pos_0), 
           byref(input_pos_1), byref(sg_op_o), byref(sg_op_0), byref(sg_op_1), 
           byref(sg_op_2), byref(sg_fract_o), byref(sg_fract_0), byref(sg_fract_1), 
           byref(lattcoords_c), byref(eps__c))

    result = FtypesResult("symmetry_module", "get_spacegroup", "Subroutine")
    result.add("atomtype", atomtype_a)
    if not all([sg_op_0.value==0, sg_op_1.value==0, sg_op_2.value==0]):
        sg_op_ft = Ftype(sg_op_o, [sg_op_0, sg_op_1, sg_op_2], libpath)
        result.add("sg_op", as_array(sg_op_o, (sg_op_2.value, sg_op_1.value, sg_op_0.value)), sg_op_ft)
    else:
        result.add("sg_op", None)
    if not all([sg_fract_0.value==0, sg_fract_1.value==0]):
        sg_fract_ft = Ftype(sg_fract_o, [sg_fract_0, sg_fract_1], libpath)
        result.add("sg_fract", as_array(sg_fract_o, (sg_fract_1.value, sg_fract_0.value)), sg_fract_ft)
    else:
        result.add("sg_fract", None)

    return result
