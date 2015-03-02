!!<summary>Module with subroutines for deallocating pointer targets created
!!by Fortran with the 'save' attribute, which persist in memory until a ctypes
!!request releases them.</summary>
module ftypes_dealloc
  use ISO_C_BINDING
  implicit none
CONTAINS
  subroutine dealloc_char_1d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    character(C_CHAR), pointer :: fptr(:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_char_1d

  subroutine dealloc_char_2d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    character(C_CHAR), pointer :: fptr(:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_char_2d

  subroutine dealloc_char_3d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    character(C_CHAR), pointer :: fptr(:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_char_3d

  subroutine dealloc_char_4d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    character(C_CHAR), pointer :: fptr(:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_char_4d

  subroutine dealloc_char_5d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    character(C_CHAR), pointer :: fptr(:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_char_5d

  subroutine dealloc_char_6d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    character(C_CHAR), pointer :: fptr(:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_char_6d

  subroutine dealloc_char_7d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    character(C_CHAR), pointer :: fptr(:,:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_char_7d

  subroutine dealloc_complex_1d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    complex(C_DOUBLE_COMPLEX), pointer :: fptr(:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_complex_1d

  subroutine dealloc_complex_2d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    complex(C_DOUBLE_COMPLEX), pointer :: fptr(:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_complex_2d

  subroutine dealloc_complex_3d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    complex(C_DOUBLE_COMPLEX), pointer :: fptr(:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_complex_3d

  subroutine dealloc_complex_4d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    complex(C_DOUBLE_COMPLEX), pointer :: fptr(:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_complex_4d

  subroutine dealloc_complex_5d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    complex(C_DOUBLE_COMPLEX), pointer :: fptr(:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_complex_5d

  subroutine dealloc_complex_6d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    complex(C_DOUBLE_COMPLEX), pointer :: fptr(:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_complex_6d

  subroutine dealloc_complex_7d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    complex(C_DOUBLE_COMPLEX), pointer :: fptr(:,:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_complex_7d

  subroutine dealloc_long_1d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_LONG), pointer :: fptr(:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_long_1d

  subroutine dealloc_long_2d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_LONG), pointer :: fptr(:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_long_2d

  subroutine dealloc_long_3d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_LONG), pointer :: fptr(:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_long_3d

  subroutine dealloc_long_4d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_LONG), pointer :: fptr(:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_long_4d

  subroutine dealloc_long_5d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_LONG), pointer :: fptr(:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_long_5d

  subroutine dealloc_long_6d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_LONG), pointer :: fptr(:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_long_6d

  subroutine dealloc_long_7d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_LONG), pointer :: fptr(:,:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_long_7d

  subroutine dealloc_int_1d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_INT), pointer :: fptr(:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_int_1d

  subroutine dealloc_int_2d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_INT), pointer :: fptr(:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_int_2d

  subroutine dealloc_int_3d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_INT), pointer :: fptr(:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_int_3d

  subroutine dealloc_int_4d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_INT), pointer :: fptr(:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_int_4d

  subroutine dealloc_int_5d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_INT), pointer :: fptr(:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_int_5d

  subroutine dealloc_int_6d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_INT), pointer :: fptr(:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_int_6d

  subroutine dealloc_int_7d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_INT), pointer :: fptr(:,:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_int_7d

  subroutine dealloc_double_1d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_DOUBLE), pointer :: fptr(:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_double_1d

  subroutine dealloc_double_2d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_DOUBLE), pointer :: fptr(:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_double_2d

  subroutine dealloc_double_3d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_DOUBLE), pointer :: fptr(:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_double_3d

  subroutine dealloc_double_4d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_DOUBLE), pointer :: fptr(:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_double_4d

  subroutine dealloc_double_5d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_DOUBLE), pointer :: fptr(:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_double_5d

  subroutine dealloc_double_6d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_DOUBLE), pointer :: fptr(:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_double_6d

  subroutine dealloc_double_7d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_DOUBLE), pointer :: fptr(:,:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_double_7d

  subroutine dealloc_float_1d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_FLOAT), pointer :: fptr(:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_float_1d

  subroutine dealloc_float_2d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_FLOAT), pointer :: fptr(:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_float_2d

  subroutine dealloc_float_3d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_FLOAT), pointer :: fptr(:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_float_3d

  subroutine dealloc_float_4d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_FLOAT), pointer :: fptr(:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_float_4d

  subroutine dealloc_float_5d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_FLOAT), pointer :: fptr(:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_float_5d

  subroutine dealloc_float_6d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_FLOAT), pointer :: fptr(:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_float_6d

  subroutine dealloc_float_7d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    real(C_FLOAT), pointer :: fptr(:,:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_float_7d

  subroutine dealloc_logical_1d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    logical(C_BOOL), pointer :: fptr(:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_logical_1d

  subroutine dealloc_logical_2d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    logical(C_BOOL), pointer :: fptr(:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_logical_2d

  subroutine dealloc_logical_3d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    logical(C_BOOL), pointer :: fptr(:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_logical_3d

  subroutine dealloc_logical_4d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    logical(C_BOOL), pointer :: fptr(:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_logical_4d

  subroutine dealloc_logical_5d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    logical(C_BOOL), pointer :: fptr(:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_logical_5d

  subroutine dealloc_logical_6d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    logical(C_BOOL), pointer :: fptr(:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_logical_6d

  subroutine dealloc_logical_7d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    logical(C_BOOL), pointer :: fptr(:,:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_logical_7d

  subroutine dealloc_short_1d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_SHORT), pointer :: fptr(:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_short_1d

  subroutine dealloc_short_2d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_SHORT), pointer :: fptr(:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_short_2d

  subroutine dealloc_short_3d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_SHORT), pointer :: fptr(:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_short_3d

  subroutine dealloc_short_4d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_SHORT), pointer :: fptr(:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_short_4d

  subroutine dealloc_short_5d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_SHORT), pointer :: fptr(:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_short_5d

  subroutine dealloc_short_6d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_SHORT), pointer :: fptr(:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_short_6d

  subroutine dealloc_short_7d(handle, n, indices) BIND(C)
    TYPE(C_PTR), intent(in) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: indices(n)
    integer(C_SHORT), pointer :: fptr(:,:,:,:,:,:,:)
    call C_F_POINTER(handle, fptr, indices)
    if (associated(fptr)) deallocate(fptr)
  end subroutine dealloc_short_7d
end module ftypes_dealloc
