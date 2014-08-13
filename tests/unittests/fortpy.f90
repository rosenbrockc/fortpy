!!<summary>Provides an interface for saving the values of multiple variable
!!types using a single call. Used as part of the FORTPY unit testing framework.</summary>
module fortpy
  implicit none
  private
  public pysave, dp, sp, si, li

  !!<member name="fileunit">I/O unit for the file to write the output to.</member>
  integer :: fileunit
  !!<member name="seeded">Specifies whether the random number generator has
  !!already been seeded for this program execution.</member>
  logical :: seeded

  integer, parameter:: dp=selected_real_kind(15,307)
  integer, parameter:: sp=selected_real_kind(6,37)
  integer, parameter:: si=selected_int_kind(1) ! very short integer -10..10 range
  integer, parameter:: li=selected_int_kind(18) ! Big integer -10^18..10^18 range

  !!<summary>Provides a single call interface for integer and real data types
  !!for single values, 1D and 2D arrays. Also handles outputting true/false for
  !!logical type variables.</summary>
  interface pysave
     module procedure pysave_integer, pysave_integer_1d, pysave_integer_2d, &
          pysave_real, pysave_real_1d, pysave_real_2d, pysave_logical
  end interface pysave

  interface randomfpy
     module procedure random_real, random_integer, random_real_1d, random_real_2d, &
          random_integer_1d, random_integer_2d
  end interface randomfpy

  interface fpyin_range
     module procedure in_range_integer, in_range_integer_1d, in_range_integer_2d, &
          in_range_real, in_range_real_1d, in_range_real_2d
  end interface fpyin_range
contains  
  logical function in_range_real_2d(variable, min, max)
    real(dp), intent(in) :: variable(:,:)
    real(dp), intent(in) :: min, max
    in_range_real_2d = any(variable .ge. min .and. variable .le. max)
    return
  end function in_range_real_2d

  logical function in_range_real_1d(variable, min, max)
    real(dp), intent(in) :: variable(:)
    real(dp), intent(in) :: min, max
    in_range_real_1d = any(variable .ge. min .and. variable .le. max)
    return
  end function in_range_real_1d

  logical function in_range_real(variable, min, max)
    real(dp), intent(in) :: variable
    real(dp), intent(in) :: min, max
    in_range_real = (variable .ge. min .and. variable .le. max)
    return
  end function in_range_real

  logical function in_range_integer_2d(variable, min, max)
    integer, intent(in) :: variable(:,:)
    integer, intent(in) :: min, max
    in_range_integer_2d = any(variable .ge. min .and. variable .le. max)
    return
  end function in_range_integer_2d

  logical function in_range_integer_1d(variable, min, max)
    integer, intent(in) :: variable(:)
    integer, intent(in) :: min, max
    in_range_integer_1d = any(variable .ge. min .and. variable .le. max)
    return
  end function in_range_integer_1d

  logical function in_range_integer(variable, min, max)
    integer, intent(in) :: variable
    integer, intent(in) :: min, max
    in_range_integer = (variable .ge. min .and. variable .le. max)
    return
  end function in_range_integer

  subroutine random_real_2d(variable, min, max)
    real(dp), intent(out) :: variable(:,:)
    real(dp), intent(in) :: min, max
    integer :: r, c, i, j

    r = size(variable, 1)
    c = size(variable, 2)
    do i = 1, r
       do j = 1, c
          call random_real(variable(i, j), min, max)
       end do
    end do
  end subroutine random_real_2d

  subroutine random_real_1d(variable, min, max)
    real(dp), intent(out) :: variable(:)
    real(dp), intent(in) :: min, max
    integer :: r, i

    r = size(variable, 1)
    do i = 1, r
       call random_real(variable(i), min, max)
    end do
  end subroutine random_real_1d

   subroutine random_real(variable, min, max)
    real(dp), intent(out) :: variable
    real(dp), intent(in) :: min, max
    real(dp) :: x
    !Since we only get real numbers in [0, 1], we have to set the value to be 
    !in the range manually.
    call random_number(x)    
    variable = min + x * (max - min)
  end subroutine random_real

 subroutine random_integer_2d(variable, min, max)
    integer, intent(out) :: variable(:,:)
    integer, intent(in) :: min, max
    integer :: r, c, i, j

    r = size(variable, 1)
    c = size(variable, 2)
    do i = 1, r
       do j = 1, c
          call random_integer(variable(i, j), min, max)
       end do
    end do
  end subroutine random_integer_2d

  subroutine random_integer_1d(variable, min, max)
    integer, intent(out) :: variable(:)
    integer, intent(in) :: min, max
    integer :: r, i

    r = size(variable, 1)
    do i = 1, r
       call random_integer(variable(i), min, max)
    end do
  end subroutine random_integer_1d

  subroutine random_integer(variable, min, max)
    integer, intent(out) :: variable
    integer, intent(in) :: min, max
    real(dp) :: x
    !Since we only get real numbers, we have to set the value to be in the range
    !and then turn it to an integer.
    call random_number(x)    
    variable = ceiling(min + x * (max - min))
  end subroutine random_integer

  !!<summary>Initializes the seed of the random number generator.</summary>
  subroutine random_init(seed)
    integer, intent(in), optional, dimension(1) :: seed
    if (.not. seeded) then
       if (present(seed)) then
          call random_seed(PUT=seed)
       else
          call random_seed()
       end if
       seeded = .true.
    end if
  end subroutine random_init

  subroutine pysave_integer(variable, filename)
    integer, intent(in) :: variable
    character(80), intent(in) :: filename

    call file_open(filename, 'integer')
    write(fileunit, '(i12)') variable

    call file_close()
  end subroutine pysave_integer

  subroutine pysave_integer_1d(variable, filename)
    integer, intent(in) :: variable(:)
    character(80), intent(in) :: filename
    integer :: c, i

    c = size(variable, 1)
    call file_open(filename, 'integer')

    do i = 1, c
       write(fileunit, '(i12)') variable(c)
    end do

    call file_close()
  end subroutine pysave_integer_1d

  subroutine pysave_integer_2d(variable, filename)
    integer, intent(in) :: variable(:,:)
    character(80), intent(in) :: filename
    integer :: r, c, i

    r = size(variable, 1)
    c = size(variable, 2)

    call file_open(filename, 'integer')
    do i = 1, r
       write(fileunit, '(<c>i12)') variable(r, :)
    end do

    call file_close()
  end subroutine pysave_integer_2d

  subroutine pysave_real(variable, filename)
    real(dp), intent(in) :: variable
    character(80), intent(in) :: filename

    call file_open(filename, 'float')
    write(fileunit, '(f12.7)') variable

    call file_close()
  end subroutine pysave_real

  subroutine pysave_real_1d(variable, filename)
    real(dp), intent(in) :: variable(:)
    character(80), intent(in) :: filename
    integer :: c, i

    c = size(variable, 1)
    call file_open(filename, 'float')
    do i = 1, c
       write(fileunit, '(f12.7)') variable(c)
    end do
    call file_close()
  end subroutine pysave_real_1d

  subroutine pysave_real_2d(variable, filename)
    real(dp), intent(in) :: variable(:,:)
    character(80), intent(in) :: filename
    integer :: r, c, i

    r = size(variable, 1)
    c = size(variable, 2)

    call file_open(filename, 'float')
    do i = 1, r
       write(fileunit, '(<c>f12.7)') variable(r, :)
    end do

    call file_close()
  end subroutine pysave_real_2d

  subroutine pysave_logical(variable, filename)
    logical, intent(in) :: variable
    character(80), intent(in) :: filename

    call file_open(filename, 'logical')
    if (variable) then
       write(fileunit, *) 'true'
    else
       write(fileunit, *) 'false'
    end if
    
    call file_close()
  end subroutine pysave_logical

  subroutine file_close()
    !This is just a one line routine, but if we decide to add cleanup later
    !it makes it easier to do it.
    close(fileunit)
  end subroutine file_close
 
  subroutine file_open(filename, template_name)
    character(80), intent(in) :: filename    
    character(*), intent(in) :: template_name
    integer :: ioerr
    open(newunit(fileunit), file=filename, status='replace', iostat=ioerr)
    if (ioerr /= 0) then
       print *, "ERROR opening file for pysave in fortpy", ioerr
    else
       write(fileunit, *) '# <fortpy version="1" template="', template_name, '"></fortpy>'
    end if
  end subroutine file_open

  !<summary>Returns lowest i/o unit number not in use.</summary>
  !<parameter name="unit">Out parameter that will contain the lowest i/o number.</parameter>
  integer function newunit(unit) result(n)
    integer, intent(out), optional :: unit
    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent
    do n = nmin, nmax
       inquire(unit=n, opened=inuse)
       if (.not. inuse) then
          if (present(unit)) unit=n
          return
       end if
    end do
    stop "newunit ERROR: available unit not found."
  end function newunit
end module fortpy
