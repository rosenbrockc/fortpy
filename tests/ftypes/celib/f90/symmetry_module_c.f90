!!<summary>Auto-generated Fortran module for interaction with ctypes
!!through python. Generated for module symmetry_module.</summary>
MODULE symmetry_module_c
  use symmetry_module
  use ISO_C_BINDING
  use num_types
  implicit none
CONTAINS
  subroutine make_primitive_c(aVecs, atomType_c, atomType_1_c, atomType_o, atom_pos_c, &
                              atom_pos_1_c, atom_pos_2_c, atom_pos_o, lattCoords_c, eps_, &
                              removed__c, removed__1_c) BIND(C)
    real(dp), intent(inout) :: aVecs(3, 3)
    integer, intent(inout) :: atomType_1_c
    integer, intent(in) :: atomType_c(atomType_1_c)
    type(C_PTR), intent(inout) :: atomType_o
    
    integer, intent(inout) :: atom_pos_1_c, atom_pos_2_c
    real(dp), intent(in) :: atom_pos_c(atom_pos_1_c, atom_pos_2_c)
    type(C_PTR), intent(inout) :: atom_pos_o
    
    logical(C_BOOL), intent(in) :: lattCoords_c
    real(dp), intent(in) :: eps_
    integer, intent(inout) :: removed__1_c
    type(C_PTR), intent(inout) :: removed__c
    
    integer, pointer :: atomType(:)
    integer(C_INT), allocatable, target, save :: atomType_t(:)
    
    real(dp), pointer :: atom_pos(:, :)
    real(C_DOUBLE), allocatable, target, save :: atom_pos_t(:, :)
    
    logical :: lattCoords
    integer, pointer :: removed_(:)
    integer(C_INT), allocatable, target, save :: removed__t(:)
    
    allocate(atomType(atomType_1_c))
    atomType = atomType_c

    allocate(atom_pos(atom_pos_1_c, atom_pos_2_c))
    atom_pos = atom_pos_c

    lattCoords = lattCoords_c
    call make_primitive(avecs, atomtype, atom_pos, lattcoords, eps_, removed_)

    atomType_t = atomType
    atomType_o = C_LOC(atomType_t)
    atomType_1_c = size(atomType, 1)

    atom_pos_t = atom_pos
    atom_pos_o = C_LOC(atom_pos_t)
    atom_pos_1_c = size(atom_pos, 1)
    atom_pos_2_c = size(atom_pos, 2)

    if (associated(removed_)) then
      removed__t = removed_
      removed__c = C_LOC(removed__t)
      removed__1_c = size(removed_, 1)
    end if

  end subroutine make_primitive_c

  subroutine does_mapping_exist_c(v, this_type, atom_pos, atom_pos_1_c, atom_pos_2_c, &
                                  atomType, atomType_1_c, mapped_c, eps) BIND(C)
    real(dp), intent(in) :: v(3)
    integer, intent(in) :: this_type
    integer, intent(in) :: atom_pos_1_c, atom_pos_2_c
    real(dp), intent(in) :: atom_pos(atom_pos_1_c, atom_pos_2_c)
    integer, intent(in) :: atomType_1_c
    integer, intent(in) :: atomType(atomType_1_c)
    logical(C_BOOL), intent(out) :: mapped_c
    real(dp), intent(in) :: eps
    logical :: mapped
    call does_mapping_exist(v, this_type, atom_pos, atomtype, mapped, eps)

    mapped_c = mapped
  end subroutine does_mapping_exist_c

  subroutine bring_into_cell_c(v, cart_to_latt, latt_to_cart, eps) BIND(C)
    real(dp), intent(inout) :: v(3)
    real(dp), intent(in) :: cart_to_latt(3, 3)
    real(dp), intent(in) :: latt_to_cart(3, 3)
    real(dp), intent(in) :: eps
    
    call bring_into_cell(v, cart_to_latt, latt_to_cart, eps)

  end subroutine bring_into_cell_c

  subroutine get_spacegroup_c(aVecs, atomType, atomType_1_c, input_pos_c, input_pos_1_c, &
                              input_pos_2_c, sg_op_c, sg_op_1_c, sg_op_2_c, sg_op_3_c, &
                              sg_fract_c, sg_fract_1_c, sg_fract_2_c, lattcoords_c, eps_) BIND(C)
    real(dp), intent(in) :: aVecs(3, 3)
    integer, intent(inout) :: atomType_1_c
    integer, intent(inout) :: atomType(atomType_1_c)
    integer, intent(in) :: input_pos_1_c, input_pos_2_c
    real(dp), intent(in) :: input_pos_c(input_pos_1_c, input_pos_2_c)
    integer, intent(inout) :: sg_op_1_c, sg_op_2_c, sg_op_3_c
    type(C_PTR), intent(inout) :: sg_op_c
    
    integer, intent(inout) :: sg_fract_1_c, sg_fract_2_c
    type(C_PTR), intent(inout) :: sg_fract_c
    
    logical(C_BOOL), intent(in) :: lattcoords_c
    real(dp), intent(in) :: eps_
    real(dp), pointer :: input_pos(:, :)
    real(dp), pointer :: sg_op(:, :, :)
    real(C_DOUBLE), allocatable, target, save :: sg_op_t(:, :, :)
    
    real(dp), pointer :: sg_fract(:, :)
    real(C_DOUBLE), allocatable, target, save :: sg_fract_t(:, :)
    
    logical :: lattcoords
    allocate(input_pos(input_pos_1_c, input_pos_2_c))
    input_pos = input_pos_c

    lattcoords = lattcoords_c
    call get_spacegroup(avecs, atomtype, input_pos, sg_op, sg_fract, lattcoords, eps_)

    sg_op_t = sg_op
    sg_op_c = C_LOC(sg_op_t)
    sg_op_1_c = size(sg_op, 1)
    sg_op_2_c = size(sg_op, 2)
    sg_op_3_c = size(sg_op, 3)

    sg_fract_t = sg_fract
    sg_fract_c = C_LOC(sg_fract_t)
    sg_fract_1_c = size(sg_fract, 1)
    sg_fract_2_c = size(sg_fract, 2)

  end subroutine get_spacegroup_c

END MODULE symmetry_module_c