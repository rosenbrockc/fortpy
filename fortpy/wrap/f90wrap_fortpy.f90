! Module fortpy defined in file fortpy.f90

subroutine f90wrap_fpy_set_verbosity(v)
    use fortpy, only: fpy_set_verbosity
    implicit none
    
    integer, intent(in) :: v
    call fpy_set_verbosity(v=v)
end subroutine f90wrap_fpy_set_verbosity

subroutine f90wrap_fpy_period_join_indices(pslist, indices, n, n0)
    use fortpy, only: fpy_period_join_indices
    implicit none
    
    character, intent(out) :: pslist
    integer, intent(in), dimension(n0) :: indices
    integer :: n
    integer :: n0
    !f2py intent(hide), depend(indices) :: n0 = shape(indices,0)
    call fpy_period_join_indices(pslist=pslist, indices=indices, n=n)
end subroutine f90wrap_fpy_period_join_indices

subroutine f90wrap_random_real(variable, min_bn, max_bn)
    use fortpy, only: random_real
    implicit none
    
    real(8), intent(out) :: variable
    real(8), intent(in) :: min_bn
    real(8), intent(in) :: max_bn
    call random_real(variable=variable, min=min_bn, max=max_bn)
end subroutine f90wrap_random_real

subroutine f90wrap_random_integer(variable, min_bn, max_bn)
    use fortpy, only: random_integer
    implicit none
    
    integer, intent(out) :: variable
    integer, intent(in) :: min_bn
    integer, intent(in) :: max_bn
    call random_integer(variable=variable, min=min_bn, max=max_bn)
end subroutine f90wrap_random_integer

subroutine f90wrap_random_init(seed, n0)
    use fortpy, only: random_init
    implicit none
    
    integer, intent(in), optional, dimension(n0) :: seed
    integer :: n0
    !f2py intent(hide), depend(seed) :: n0 = shape(seed,0)
    call random_init(seed=seed)
end subroutine f90wrap_random_init

subroutine f90wrap_char_escape_word(word, escaped)
    use fortpy, only: char_escape_word
    implicit none
    
    character, intent(in) :: word
    character, intent(inout) :: escaped
    call char_escape_word(word=word, escaped=escaped)
end subroutine f90wrap_char_escape_word

subroutine f90wrap_char_write_trimmed(variable, n0)
    use fortpy, only: char_write_trimmed
    implicit none
    
    character, intent(in), dimension(n0) :: variable
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call char_write_trimmed(variable=variable)
end subroutine f90wrap_char_write_trimmed

subroutine f90wrap_file_open(filename, n, template_name)
    use fortpy, only: file_open
    implicit none
    
    character, intent(in) :: filename
    integer, intent(in) :: n
    character, intent(in) :: template_name
    call file_open(filename=filename, n=n, template_name=template_name)
end subroutine f90wrap_file_open

subroutine f90wrap_file_close
    use fortpy, only: file_close
    implicit none
    
    call file_close()
end subroutine f90wrap_file_close

subroutine f90wrap_fpy_value_count(line, length, ret_fpy_value_count, ischar)
    use fortpy, only: fpy_value_count
    implicit none
    
    character, intent(in) :: line
    integer, intent(in) :: length
    integer, intent(out) :: ret_fpy_value_count
    logical, optional, intent(in) :: ischar
    ret_fpy_value_count = fpy_value_count(line=line, length=length, ischar=ischar)
end subroutine f90wrap_fpy_value_count

subroutine f90wrap_fpy_linevalue_count(filename, commentchar, nlines, nvalues, &
    ischar)
    use fortpy, only: fpy_linevalue_count
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer, intent(out) :: nlines
    integer, intent(out) :: nvalues
    logical, optional :: ischar
    call fpy_linevalue_count(filename=filename, commentchar=commentchar, &
        nlines=nlines, nvalues=nvalues, ischar=ischar)
end subroutine f90wrap_fpy_linevalue_count

subroutine f90wrap_fpy_newunit(ret_fpy_newunit, unit)
    use fortpy, only: fpy_newunit
    implicit none
    
    integer, intent(out) :: ret_fpy_newunit
    integer, intent(out), optional :: unit
    ret_fpy_newunit = fpy_newunit(unit=unit)
end subroutine f90wrap_fpy_newunit

subroutine f90wrap_pysave_realsp_0d(variable, filename)
    use fortpy, only: pysave
    implicit none
    
    real(4), intent(in) :: variable
    character, intent(in) :: filename
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realsp_0d

subroutine f90wrap_pysave_realsp_1d(variable, filename, n0)
    use fortpy, only: pysave
    implicit none
    
    real(4), intent(in), dimension(n0) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realsp_1d

subroutine f90wrap_pysave_realsp_2d(variable, filename, n0, n1)
    use fortpy, only: pysave
    implicit none
    
    real(4), intent(in), dimension(n0,n1) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realsp_2d

subroutine f90wrap_pysave_realsp_3d(variable, filename, n0, n1, n2)
    use fortpy, only: pysave
    implicit none
    
    real(4), intent(in), dimension(n0,n1,n2) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realsp_3d

subroutine f90wrap_pysave_realsp_4d(variable, filename, n0, n1, n2, n3)
    use fortpy, only: pysave
    implicit none
    
    real(4), intent(in), dimension(n0,n1,n2,n3) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realsp_4d

subroutine f90wrap_pysave_realsp_5d(variable, filename, n0, n1, n2, n3, n4)
    use fortpy, only: pysave
    implicit none
    
    real(4), intent(in), dimension(n0,n1,n2,n3,n4) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realsp_5d

subroutine f90wrap_pysave_realsp_6d(variable, filename, n0, n1, n2, n3, n4, n5)
    use fortpy, only: pysave
    implicit none
    
    real(4), intent(in), dimension(n0,n1,n2,n3,n4,n5) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realsp_6d

subroutine f90wrap_pysave_realsp_7d(variable, filename, n0, n1, n2, n3, n4, n5, &
    n6)
    use fortpy, only: pysave
    implicit none
    
    real(4), intent(in), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realsp_7d

subroutine f90wrap_pysave_realdp_0d(variable, filename)
    use fortpy, only: pysave
    implicit none
    
    real(8), intent(in) :: variable
    character, intent(in) :: filename
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realdp_0d

subroutine f90wrap_pysave_realdp_1d(variable, filename, n0)
    use fortpy, only: pysave
    implicit none
    
    real(8), intent(in), dimension(n0) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realdp_1d

subroutine f90wrap_pysave_realdp_2d(variable, filename, n0, n1)
    use fortpy, only: pysave
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realdp_2d

subroutine f90wrap_pysave_realdp_3d(variable, filename, n0, n1, n2)
    use fortpy, only: pysave
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realdp_3d

subroutine f90wrap_pysave_realdp_4d(variable, filename, n0, n1, n2, n3)
    use fortpy, only: pysave
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realdp_4d

subroutine f90wrap_pysave_realdp_5d(variable, filename, n0, n1, n2, n3, n4)
    use fortpy, only: pysave
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3,n4) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realdp_5d

subroutine f90wrap_pysave_realdp_6d(variable, filename, n0, n1, n2, n3, n4, n5)
    use fortpy, only: pysave
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3,n4,n5) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realdp_6d

subroutine f90wrap_pysave_realdp_7d(variable, filename, n0, n1, n2, n3, n4, n5, &
    n6)
    use fortpy, only: pysave
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_realdp_7d

subroutine f90wrap_pysave_integer_0d(variable, filename)
    use fortpy, only: pysave
    implicit none
    
    integer, intent(in) :: variable
    character, intent(in) :: filename
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integer_0d

subroutine f90wrap_pysave_integer_1d(variable, filename, n0)
    use fortpy, only: pysave
    implicit none
    
    integer, intent(in), dimension(n0) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integer_1d

subroutine f90wrap_pysave_integer_2d(variable, filename, n0, n1)
    use fortpy, only: pysave
    implicit none
    
    integer, intent(in), dimension(n0,n1) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integer_2d

subroutine f90wrap_pysave_integer_3d(variable, filename, n0, n1, n2)
    use fortpy, only: pysave
    implicit none
    
    integer, intent(in), dimension(n0,n1,n2) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integer_3d

subroutine f90wrap_pysave_integer_4d(variable, filename, n0, n1, n2, n3)
    use fortpy, only: pysave
    implicit none
    
    integer, intent(in), dimension(n0,n1,n2,n3) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integer_4d

subroutine f90wrap_pysave_integer_5d(variable, filename, n0, n1, n2, n3, n4)
    use fortpy, only: pysave
    implicit none
    
    integer, intent(in), dimension(n0,n1,n2,n3,n4) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integer_5d

subroutine f90wrap_pysave_integer_6d(variable, filename, n0, n1, n2, n3, n4, n5)
    use fortpy, only: pysave
    implicit none
    
    integer, intent(in), dimension(n0,n1,n2,n3,n4,n5) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integer_6d

subroutine f90wrap_pysave_integer_7d(variable, filename, n0, n1, n2, n3, n4, n5, &
    n6)
    use fortpy, only: pysave
    implicit none
    
    integer, intent(in), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integer_7d

subroutine f90wrap_pysave_integersp_0d(variable, filename)
    use fortpy, only: pysave
    implicit none
    
    integer(4), intent(in) :: variable
    character, intent(in) :: filename
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integersp_0d

subroutine f90wrap_pysave_integersp_1d(variable, filename, n0)
    use fortpy, only: pysave
    implicit none
    
    integer(4), intent(in), dimension(n0) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integersp_1d

subroutine f90wrap_pysave_integersp_2d(variable, filename, n0, n1)
    use fortpy, only: pysave
    implicit none
    
    integer(4), intent(in), dimension(n0,n1) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integersp_2d

subroutine f90wrap_pysave_integersp_3d(variable, filename, n0, n1, n2)
    use fortpy, only: pysave
    implicit none
    
    integer(4), intent(in), dimension(n0,n1,n2) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integersp_3d

subroutine f90wrap_pysave_integersp_4d(variable, filename, n0, n1, n2, n3)
    use fortpy, only: pysave
    implicit none
    
    integer(4), intent(in), dimension(n0,n1,n2,n3) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integersp_4d

subroutine f90wrap_pysave_integersp_5d(variable, filename, n0, n1, n2, n3, n4)
    use fortpy, only: pysave
    implicit none
    
    integer(4), intent(in), dimension(n0,n1,n2,n3,n4) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integersp_5d

subroutine f90wrap_pysave_integersp_6d(variable, filename, n0, n1, n2, n3, n4, &
    n5)
    use fortpy, only: pysave
    implicit none
    
    integer(4), intent(in), dimension(n0,n1,n2,n3,n4,n5) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integersp_6d

subroutine f90wrap_pysave_integersp_7d(variable, filename, n0, n1, n2, n3, n4, &
    n5, n6)
    use fortpy, only: pysave
    implicit none
    
    integer(4), intent(in), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integersp_7d

subroutine f90wrap_pysave_integerdp_0d(variable, filename)
    use fortpy, only: pysave
    implicit none
    
    integer(8), intent(in) :: variable
    character, intent(in) :: filename
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integerdp_0d

subroutine f90wrap_pysave_integerdp_1d(variable, filename, n0)
    use fortpy, only: pysave
    implicit none
    
    integer(8), intent(in), dimension(n0) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integerdp_1d

subroutine f90wrap_pysave_integerdp_2d(variable, filename, n0, n1)
    use fortpy, only: pysave
    implicit none
    
    integer(8), intent(in), dimension(n0,n1) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integerdp_2d

subroutine f90wrap_pysave_integerdp_3d(variable, filename, n0, n1, n2)
    use fortpy, only: pysave
    implicit none
    
    integer(8), intent(in), dimension(n0,n1,n2) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integerdp_3d

subroutine f90wrap_pysave_integerdp_4d(variable, filename, n0, n1, n2, n3)
    use fortpy, only: pysave
    implicit none
    
    integer(8), intent(in), dimension(n0,n1,n2,n3) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integerdp_4d

subroutine f90wrap_pysave_integerdp_5d(variable, filename, n0, n1, n2, n3, n4)
    use fortpy, only: pysave
    implicit none
    
    integer(8), intent(in), dimension(n0,n1,n2,n3,n4) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integerdp_5d

subroutine f90wrap_pysave_integerdp_6d(variable, filename, n0, n1, n2, n3, n4, &
    n5)
    use fortpy, only: pysave
    implicit none
    
    integer(8), intent(in), dimension(n0,n1,n2,n3,n4,n5) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integerdp_6d

subroutine f90wrap_pysave_integerdp_7d(variable, filename, n0, n1, n2, n3, n4, &
    n5, n6)
    use fortpy, only: pysave
    implicit none
    
    integer(8), intent(in), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_integerdp_7d

subroutine f90wrap_pysave_complexsp_0d(variable, filename)
    use fortpy, only: pysave
    implicit none
    
    complex(4), intent(in) :: variable
    character, intent(in) :: filename
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexsp_0d

subroutine f90wrap_pysave_complexsp_1d(variable, filename, n0)
    use fortpy, only: pysave
    implicit none
    
    complex(4), intent(in), dimension(n0) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexsp_1d

subroutine f90wrap_pysave_complexsp_2d(variable, filename, n0, n1)
    use fortpy, only: pysave
    implicit none
    
    complex(4), intent(in), dimension(n0,n1) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexsp_2d

subroutine f90wrap_pysave_complexsp_3d(variable, filename, n0, n1, n2)
    use fortpy, only: pysave
    implicit none
    
    complex(4), intent(in), dimension(n0,n1,n2) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexsp_3d

subroutine f90wrap_pysave_complexsp_4d(variable, filename, n0, n1, n2, n3)
    use fortpy, only: pysave
    implicit none
    
    complex(4), intent(in), dimension(n0,n1,n2,n3) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexsp_4d

subroutine f90wrap_pysave_complexsp_5d(variable, filename, n0, n1, n2, n3, n4)
    use fortpy, only: pysave
    implicit none
    
    complex(4), intent(in), dimension(n0,n1,n2,n3,n4) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexsp_5d

subroutine f90wrap_pysave_complexsp_6d(variable, filename, n0, n1, n2, n3, n4, &
    n5)
    use fortpy, only: pysave
    implicit none
    
    complex(4), intent(in), dimension(n0,n1,n2,n3,n4,n5) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexsp_6d

subroutine f90wrap_pysave_complexsp_7d(variable, filename, n0, n1, n2, n3, n4, &
    n5, n6)
    use fortpy, only: pysave
    implicit none
    
    complex(4), intent(in), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexsp_7d

subroutine f90wrap_pysave_complexdp_0d(variable, filename)
    use fortpy, only: pysave
    implicit none
    
    complex(8), intent(in) :: variable
    character, intent(in) :: filename
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexdp_0d

subroutine f90wrap_pysave_complexdp_1d(variable, filename, n0)
    use fortpy, only: pysave
    implicit none
    
    complex(8), intent(in), dimension(n0) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexdp_1d

subroutine f90wrap_pysave_complexdp_2d(variable, filename, n0, n1)
    use fortpy, only: pysave
    implicit none
    
    complex(8), intent(in), dimension(n0,n1) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexdp_2d

subroutine f90wrap_pysave_complexdp_3d(variable, filename, n0, n1, n2)
    use fortpy, only: pysave
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexdp_3d

subroutine f90wrap_pysave_complexdp_4d(variable, filename, n0, n1, n2, n3)
    use fortpy, only: pysave
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2,n3) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexdp_4d

subroutine f90wrap_pysave_complexdp_5d(variable, filename, n0, n1, n2, n3, n4)
    use fortpy, only: pysave
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2,n3,n4) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexdp_5d

subroutine f90wrap_pysave_complexdp_6d(variable, filename, n0, n1, n2, n3, n4, &
    n5)
    use fortpy, only: pysave
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2,n3,n4,n5) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexdp_6d

subroutine f90wrap_pysave_complexdp_7d(variable, filename, n0, n1, n2, n3, n4, &
    n5, n6)
    use fortpy, only: pysave
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_complexdp_7d

subroutine f90wrap_pysave_character_0d(variable, filename)
    use fortpy, only: pysave
    implicit none
    
    character, intent(in) :: variable
    character, intent(in) :: filename
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_character_0d

subroutine f90wrap_pysave_character_1d(variable, filename, n0)
    use fortpy, only: pysave
    implicit none
    
    character, intent(in), dimension(n0) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_character_1d

subroutine f90wrap_pysave_character_2d(variable, filename, n0, n1)
    use fortpy, only: pysave
    implicit none
    
    character, intent(in), dimension(n0,n1) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_character_2d

subroutine f90wrap_pysave_character_3d(variable, filename, n0, n1, n2)
    use fortpy, only: pysave
    implicit none
    
    character, intent(in), dimension(n0,n1,n2) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_character_3d

subroutine f90wrap_pysave_character_4d(variable, filename, n0, n1, n2, n3)
    use fortpy, only: pysave
    implicit none
    
    character, intent(in), dimension(n0,n1,n2,n3) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_character_4d

subroutine f90wrap_pysave_character_5d(variable, filename, n0, n1, n2, n3, n4)
    use fortpy, only: pysave
    implicit none
    
    character, intent(in), dimension(n0,n1,n2,n3,n4) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_character_5d

subroutine f90wrap_pysave_character_6d(variable, filename, n0, n1, n2, n3, n4, &
    n5)
    use fortpy, only: pysave
    implicit none
    
    character, intent(in), dimension(n0,n1,n2,n3,n4,n5) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_character_6d

subroutine f90wrap_pysave_character_7d(variable, filename, n0, n1, n2, n3, n4, &
    n5, n6)
    use fortpy, only: pysave
    implicit none
    
    character, intent(in), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_character_7d

subroutine f90wrap_pysave_logical_0d(variable, filename)
    use fortpy, only: pysave
    implicit none
    
    logical, intent(in) :: variable
    character, intent(in) :: filename
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_logical_0d

subroutine f90wrap_pysave_logical_1d(variable, filename, n0)
    use fortpy, only: pysave
    implicit none
    
    logical, intent(in), dimension(n0) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_logical_1d

subroutine f90wrap_pysave_logical_2d(variable, filename, n0, n1)
    use fortpy, only: pysave
    implicit none
    
    logical, intent(in), dimension(n0,n1) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_logical_2d

subroutine f90wrap_pysave_logical_3d(variable, filename, n0, n1, n2)
    use fortpy, only: pysave
    implicit none
    
    logical, intent(in), dimension(n0,n1,n2) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_logical_3d

subroutine f90wrap_pysave_logical_4d(variable, filename, n0, n1, n2, n3)
    use fortpy, only: pysave
    implicit none
    
    logical, intent(in), dimension(n0,n1,n2,n3) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_logical_4d

subroutine f90wrap_pysave_logical_5d(variable, filename, n0, n1, n2, n3, n4)
    use fortpy, only: pysave
    implicit none
    
    logical, intent(in), dimension(n0,n1,n2,n3,n4) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_logical_5d

subroutine f90wrap_pysave_logical_6d(variable, filename, n0, n1, n2, n3, n4, n5)
    use fortpy, only: pysave
    implicit none
    
    logical, intent(in), dimension(n0,n1,n2,n3,n4,n5) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_logical_6d

subroutine f90wrap_pysave_logical_7d(variable, filename, n0, n1, n2, n3, n4, n5, &
    n6)
    use fortpy, only: pysave
    implicit none
    
    logical, intent(in), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    character, intent(in) :: filename
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call pysave(variable=variable, filename=filename)
end subroutine f90wrap_pysave_logical_7d

subroutine f90wrap_fpy_read_realsp_0d(filename, commentchar, variable, success_, &
    strict_)
    use fortpy, only: fpy_read
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(4), intent(inout) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    call fpy_read(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realsp_0d

subroutine f90wrap_fpy_read_realdp_0d(filename, commentchar, variable, success_, &
    strict_)
    use fortpy, only: fpy_read
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(8), intent(inout) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    call fpy_read(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realdp_0d

subroutine f90wrap_fpy_read_integer_0d(filename, commentchar, variable, &
    success_, strict_)
    use fortpy, only: fpy_read
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer, intent(inout) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    call fpy_read(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integer_0d

subroutine f90wrap_fpy_read_integersp_0d(filename, commentchar, variable, &
    success_, strict_)
    use fortpy, only: fpy_read
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(4), intent(inout) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    call fpy_read(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integersp_0d

subroutine f90wrap_fpy_read_integerdp_0d(filename, commentchar, variable, &
    success_, strict_)
    use fortpy, only: fpy_read
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(8), intent(inout) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    call fpy_read(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integerdp_0d

subroutine f90wrap_fpy_read_complexsp_0d(filename, commentchar, variable, &
    success_, strict_)
    use fortpy, only: fpy_read
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(4), intent(inout) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    call fpy_read(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexsp_0d

subroutine f90wrap_fpy_read_complexdp_0d(filename, commentchar, variable, &
    success_, strict_)
    use fortpy, only: fpy_read
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(8), intent(inout) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    call fpy_read(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexdp_0d

subroutine f90wrap_fpy_read_character_0d(filename, commentchar, variable, &
    success_, strict_)
    use fortpy, only: fpy_read
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    character, intent(inout) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    call fpy_read(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_character_0d

subroutine f90wrap_fpy_read_logical_0d(filename, commentchar, variable, &
    success_, strict_)
    use fortpy, only: fpy_read
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    logical, intent(inout) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    call fpy_read(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_logical_0d

subroutine f90wrap_fpy_read_realsp_1df(filename, commentchar, variable, &
    success_, strict_, n0)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(4), intent(inout), dimension(n0) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realsp_1df

subroutine f90wrap_fpy_read_realsp_2df(filename, commentchar, variable, &
    success_, strict_, n0, n1)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(4), intent(inout), dimension(n0,n1) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realsp_2df

subroutine f90wrap_fpy_read_realsp_3df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(4), intent(inout), dimension(n0,n1,n2) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realsp_3df

subroutine f90wrap_fpy_read_realsp_4df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(4), intent(inout), dimension(n0,n1,n2,n3) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realsp_4df

subroutine f90wrap_fpy_read_realsp_5df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(4), intent(inout), dimension(n0,n1,n2,n3,n4) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realsp_5df

subroutine f90wrap_fpy_read_realsp_6df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(4), intent(inout), dimension(n0,n1,n2,n3,n4,n5) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realsp_6df

subroutine f90wrap_fpy_read_realsp_7df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5, n6)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(4), intent(inout), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realsp_7df

subroutine f90wrap_fpy_read_realdp_1df(filename, commentchar, variable, &
    success_, strict_, n0)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(8), intent(inout), dimension(n0) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realdp_1df

subroutine f90wrap_fpy_read_realdp_2df(filename, commentchar, variable, &
    success_, strict_, n0, n1)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(8), intent(inout), dimension(n0,n1) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realdp_2df

subroutine f90wrap_fpy_read_realdp_3df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(8), intent(inout), dimension(n0,n1,n2) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realdp_3df

subroutine f90wrap_fpy_read_realdp_4df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(8), intent(inout), dimension(n0,n1,n2,n3) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realdp_4df

subroutine f90wrap_fpy_read_realdp_5df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(8), intent(inout), dimension(n0,n1,n2,n3,n4) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realdp_5df

subroutine f90wrap_fpy_read_realdp_6df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(8), intent(inout), dimension(n0,n1,n2,n3,n4,n5) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realdp_6df

subroutine f90wrap_fpy_read_realdp_7df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5, n6)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    real(8), intent(inout), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_realdp_7df

subroutine f90wrap_fpy_read_integer_1df(filename, commentchar, variable, &
    success_, strict_, n0)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer, intent(inout), dimension(n0) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integer_1df

subroutine f90wrap_fpy_read_integer_2df(filename, commentchar, variable, &
    success_, strict_, n0, n1)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer, intent(inout), dimension(n0,n1) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integer_2df

subroutine f90wrap_fpy_read_integer_3df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer, intent(inout), dimension(n0,n1,n2) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integer_3df

subroutine f90wrap_fpy_read_integer_4df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer, intent(inout), dimension(n0,n1,n2,n3) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integer_4df

subroutine f90wrap_fpy_read_integer_5df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer, intent(inout), dimension(n0,n1,n2,n3,n4) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integer_5df

subroutine f90wrap_fpy_read_integer_6df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer, intent(inout), dimension(n0,n1,n2,n3,n4,n5) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integer_6df

subroutine f90wrap_fpy_read_integer_7df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5, n6)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer, intent(inout), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integer_7df

subroutine f90wrap_fpy_read_integersp_1df(filename, commentchar, variable, &
    success_, strict_, n0)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(4), intent(inout), dimension(n0) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integersp_1df

subroutine f90wrap_fpy_read_integersp_2df(filename, commentchar, variable, &
    success_, strict_, n0, n1)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(4), intent(inout), dimension(n0,n1) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integersp_2df

subroutine f90wrap_fpy_read_integersp_3df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(4), intent(inout), dimension(n0,n1,n2) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integersp_3df

subroutine f90wrap_fpy_read_integersp_4df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(4), intent(inout), dimension(n0,n1,n2,n3) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integersp_4df

subroutine f90wrap_fpy_read_integersp_5df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(4), intent(inout), dimension(n0,n1,n2,n3,n4) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integersp_5df

subroutine f90wrap_fpy_read_integersp_6df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(4), intent(inout), dimension(n0,n1,n2,n3,n4,n5) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integersp_6df

subroutine f90wrap_fpy_read_integersp_7df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5, n6)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(4), intent(inout), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integersp_7df

subroutine f90wrap_fpy_read_integerdp_1df(filename, commentchar, variable, &
    success_, strict_, n0)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(8), intent(inout), dimension(n0) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integerdp_1df

subroutine f90wrap_fpy_read_integerdp_2df(filename, commentchar, variable, &
    success_, strict_, n0, n1)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(8), intent(inout), dimension(n0,n1) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integerdp_2df

subroutine f90wrap_fpy_read_integerdp_3df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(8), intent(inout), dimension(n0,n1,n2) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integerdp_3df

subroutine f90wrap_fpy_read_integerdp_4df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(8), intent(inout), dimension(n0,n1,n2,n3) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integerdp_4df

subroutine f90wrap_fpy_read_integerdp_5df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(8), intent(inout), dimension(n0,n1,n2,n3,n4) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integerdp_5df

subroutine f90wrap_fpy_read_integerdp_6df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(8), intent(inout), dimension(n0,n1,n2,n3,n4,n5) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integerdp_6df

subroutine f90wrap_fpy_read_integerdp_7df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5, n6)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    integer(8), intent(inout), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_integerdp_7df

subroutine f90wrap_fpy_read_complexsp_1df(filename, commentchar, variable, &
    success_, strict_, n0)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(4), intent(inout), dimension(n0) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexsp_1df

subroutine f90wrap_fpy_read_complexsp_2df(filename, commentchar, variable, &
    success_, strict_, n0, n1)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(4), intent(inout), dimension(n0,n1) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexsp_2df

subroutine f90wrap_fpy_read_complexsp_3df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(4), intent(inout), dimension(n0,n1,n2) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexsp_3df

subroutine f90wrap_fpy_read_complexsp_4df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(4), intent(inout), dimension(n0,n1,n2,n3) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexsp_4df

subroutine f90wrap_fpy_read_complexsp_5df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(4), intent(inout), dimension(n0,n1,n2,n3,n4) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexsp_5df

subroutine f90wrap_fpy_read_complexsp_6df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(4), intent(inout), dimension(n0,n1,n2,n3,n4,n5) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexsp_6df

subroutine f90wrap_fpy_read_complexsp_7df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5, n6)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(4), intent(inout), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexsp_7df

subroutine f90wrap_fpy_read_complexdp_1df(filename, commentchar, variable, &
    success_, strict_, n0)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(8), intent(inout), dimension(n0) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexdp_1df

subroutine f90wrap_fpy_read_complexdp_2df(filename, commentchar, variable, &
    success_, strict_, n0, n1)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(8), intent(inout), dimension(n0,n1) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexdp_2df

subroutine f90wrap_fpy_read_complexdp_3df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(8), intent(inout), dimension(n0,n1,n2) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexdp_3df

subroutine f90wrap_fpy_read_complexdp_4df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(8), intent(inout), dimension(n0,n1,n2,n3) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexdp_4df

subroutine f90wrap_fpy_read_complexdp_5df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(8), intent(inout), dimension(n0,n1,n2,n3,n4) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexdp_5df

subroutine f90wrap_fpy_read_complexdp_6df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(8), intent(inout), dimension(n0,n1,n2,n3,n4,n5) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexdp_6df

subroutine f90wrap_fpy_read_complexdp_7df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5, n6)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    complex(8), intent(inout), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_complexdp_7df

subroutine f90wrap_fpy_read_character_1df(filename, commentchar, variable, &
    success_, strict_, n0)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    character, intent(inout), dimension(n0) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_character_1df

subroutine f90wrap_fpy_read_character_2df(filename, commentchar, variable, &
    success_, strict_, n0, n1)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    character, intent(inout), dimension(n0,n1) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_character_2df

subroutine f90wrap_fpy_read_character_3df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    character, intent(inout), dimension(n0,n1,n2) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_character_3df

subroutine f90wrap_fpy_read_character_4df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    character, intent(inout), dimension(n0,n1,n2,n3) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_character_4df

subroutine f90wrap_fpy_read_character_5df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    character, intent(inout), dimension(n0,n1,n2,n3,n4) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_character_5df

subroutine f90wrap_fpy_read_character_6df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    character, intent(inout), dimension(n0,n1,n2,n3,n4,n5) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_character_6df

subroutine f90wrap_fpy_read_character_7df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5, n6)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    character, intent(inout), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_character_7df

subroutine f90wrap_fpy_read_logical_1df(filename, commentchar, variable, &
    success_, strict_, n0)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    logical, intent(inout), dimension(n0) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_logical_1df

subroutine f90wrap_fpy_read_logical_2df(filename, commentchar, variable, &
    success_, strict_, n0, n1)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    logical, intent(inout), dimension(n0,n1) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_logical_2df

subroutine f90wrap_fpy_read_logical_3df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    logical, intent(inout), dimension(n0,n1,n2) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_logical_3df

subroutine f90wrap_fpy_read_logical_4df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    logical, intent(inout), dimension(n0,n1,n2,n3) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_logical_4df

subroutine f90wrap_fpy_read_logical_5df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    logical, intent(inout), dimension(n0,n1,n2,n3,n4) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_logical_5df

subroutine f90wrap_fpy_read_logical_6df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    logical, intent(inout), dimension(n0,n1,n2,n3,n4,n5) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_logical_6df

subroutine f90wrap_fpy_read_logical_7df(filename, commentchar, variable, &
    success_, strict_, n0, n1, n2, n3, n4, n5, n6)
    use fortpy, only: fpy_read_f
    implicit none
    
    character, intent(in) :: filename
    character, intent(in) :: commentchar
    logical, intent(inout), dimension(n0,n1,n2,n3,n4,n5,n6) :: variable
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer :: n0
    !f2py intent(hide), depend(variable) :: n0 = shape(variable,0)
    integer :: n1
    !f2py intent(hide), depend(variable) :: n1 = shape(variable,1)
    integer :: n2
    !f2py intent(hide), depend(variable) :: n2 = shape(variable,2)
    integer :: n3
    !f2py intent(hide), depend(variable) :: n3 = shape(variable,3)
    integer :: n4
    !f2py intent(hide), depend(variable) :: n4 = shape(variable,4)
    integer :: n5
    !f2py intent(hide), depend(variable) :: n5 = shape(variable,5)
    integer :: n6
    !f2py intent(hide), depend(variable) :: n6 = shape(variable,6)
    call fpy_read_f(filename=filename, commentchar=commentchar, variable=variable, &
        success_=success_, strict_=strict_)
end subroutine f90wrap_fpy_read_logical_7df

subroutine f90wrap_fortpy__get__fileunit(f90wrap_fileunit)
    use fortpy, only: fortpy_fileunit => fileunit
    implicit none
    integer, intent(out) :: f90wrap_fileunit
    
    f90wrap_fileunit = fortpy_fileunit
end subroutine f90wrap_fortpy__get__fileunit

subroutine f90wrap_fortpy__set__fileunit(f90wrap_fileunit)
    use fortpy, only: fortpy_fileunit => fileunit
    implicit none
    integer, intent(in) :: f90wrap_fileunit
    
    fortpy_fileunit = f90wrap_fileunit
end subroutine f90wrap_fortpy__set__fileunit

subroutine f90wrap_fortpy__get__fpy_verbose(f90wrap_fpy_verbose)
    use fortpy, only: fortpy_fpy_verbose => fpy_verbose
    implicit none
    integer, intent(out) :: f90wrap_fpy_verbose
    
    f90wrap_fpy_verbose = fortpy_fpy_verbose
end subroutine f90wrap_fortpy__get__fpy_verbose

subroutine f90wrap_fortpy__set__fpy_verbose(f90wrap_fpy_verbose)
    use fortpy, only: fortpy_fpy_verbose => fpy_verbose
    implicit none
    integer, intent(in) :: f90wrap_fpy_verbose
    
    fortpy_fpy_verbose = f90wrap_fpy_verbose
end subroutine f90wrap_fortpy__set__fpy_verbose

subroutine f90wrap_fortpy__get__seeded(f90wrap_seeded)
    use fortpy, only: fortpy_seeded => seeded
    implicit none
    logical, intent(out) :: f90wrap_seeded
    
    f90wrap_seeded = fortpy_seeded
end subroutine f90wrap_fortpy__get__seeded

subroutine f90wrap_fortpy__set__seeded(f90wrap_seeded)
    use fortpy, only: fortpy_seeded => seeded
    implicit none
    logical, intent(in) :: f90wrap_seeded
    
    fortpy_seeded = f90wrap_seeded
end subroutine f90wrap_fortpy__set__seeded

subroutine f90wrap_fortpy__get__fdp(f90wrap_fdp)
    use fortpy, only: fortpy_fdp => fdp
    implicit none
    integer, intent(out) :: f90wrap_fdp
    
    f90wrap_fdp = fortpy_fdp
end subroutine f90wrap_fortpy__get__fdp

subroutine f90wrap_fortpy__get__fsp(f90wrap_fsp)
    use fortpy, only: fortpy_fsp => fsp
    implicit none
    integer, intent(out) :: f90wrap_fsp
    
    f90wrap_fsp = fortpy_fsp
end subroutine f90wrap_fortpy__get__fsp

subroutine f90wrap_fortpy__get__fsi(f90wrap_fsi)
    use fortpy, only: fortpy_fsi => fsi
    implicit none
    integer, intent(out) :: f90wrap_fsi
    
    f90wrap_fsi = fortpy_fsi
end subroutine f90wrap_fortpy__get__fsi

subroutine f90wrap_fortpy__get__fli(f90wrap_fli)
    use fortpy, only: fortpy_fli => fli
    implicit none
    integer, intent(out) :: f90wrap_fli
    
    f90wrap_fli = fortpy_fli
end subroutine f90wrap_fortpy__get__fli

! End of module fortpy defined in file fortpy.f90

