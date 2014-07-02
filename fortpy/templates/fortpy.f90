!!<summary>Provides an interface for saving the values of multiple variable
!!types using a single call. Used as part of the FORTPY unit testing framework.</summary>
module fortpy
  implicit none
  private
  public pysave, dp, sp, si, li, fpy_linevalue_count, fpy_newunit, pysave_integer_li, &
       fpy_value_count

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
          pysave_integer_li, pysave_integer_1d_li, pysave_integer_2d_li, &
          pysave_real, pysave_real_1d, pysave_real_2d, pysave_logical
  end interface pysave

  !!<summary>Generates a 0-2 dimensional array of integer or real values.</summary>
  interface randomfpy
     module procedure random_real, random_integer, random_real_1d, random_real_2d, &
          random_integer_1d, random_integer_2d
  end interface randomfpy

  !!<summary>Determines whether the specified value lies within the range for
  !!0-2 dimensional integer and real values.</summary>
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
    integer, intent(in), optional :: seed(:)
    if (.not. seeded) then
       if (present(seed)) then
          call random_seed(PUT=seed)
       else
          call random_seed()
       end if
       seeded = .true.
    end if
  end subroutine random_init

  subroutine pysave_integer(variable, filename, n)
    integer, intent(in) :: variable
    integer, intent(in) :: n
    character(n), intent(in) :: filename

    call file_open(filename, n, 'integer')
    write(fileunit, '(i12)') variable

    call file_close()
  end subroutine pysave_integer

  subroutine pysave_integer_1d(variable, filename, n)
    integer, intent(in) :: variable(:)
    integer, intent(in) :: n
    character(n), intent(in) :: filename
    integer :: c, i

    c = size(variable, 1)
    call file_open(filename, n, 'integer')

    do i = 1, c
       write(fileunit, '(i12)') variable(c)
    end do

    call file_close()
  end subroutine pysave_integer_1d

  subroutine pysave_integer_2d(variable, filename, n)
    integer, intent(in) :: variable(:,:)
    integer, intent(in) :: n
    character(n), intent(in) :: filename
    character(20) :: FMT
    integer :: r, c, i

    r = size(variable, 1)
    c = size(variable, 2)
    write(FMT, *)

    call file_open(filename, n, 'integer')
    do i = 1, r
       write(fileunit, '('// adjustl(FMT) // 'i12)') variable(r, :)
    end do

    call file_close()
  end subroutine pysave_integer_2d

  subroutine pysave_integer_li(variable, filename, n)
    integer(li), intent(in) :: variable
    integer, intent(in) :: n
    character(n), intent(in) :: filename

    call file_open(filename, n, 'integer')
    write(fileunit, '(i25)') variable

    call file_close()
  end subroutine pysave_integer_li

  subroutine pysave_integer_1d_li(variable, filename, n)
    integer(li), intent(in) :: variable(:)
    integer, intent(in) :: n
    character(n), intent(in) :: filename
    integer :: c, i

    c = size(variable, 1)
    call file_open(filename, n, 'integer')

    do i = 1, c
       write(fileunit, '(i25)') variable(c)
    end do

    call file_close()
  end subroutine pysave_integer_1d_li

  subroutine pysave_integer_2d_li(variable, filename, n)
    integer(li), intent(in) :: variable(:,:)
    integer, intent(in) :: n
    character(n), intent(in) :: filename
    character(20) :: FMT
    integer :: r, c, i

    r = size(variable, 1)
    c = size(variable, 2)
    write(FMT, *)

    call file_open(filename, n, 'integer')
    do i = 1, r
       write(fileunit, '('// adjustl(FMT) // 'i25)') variable(r, :)
    end do

    call file_close()
  end subroutine pysave_integer_2d_li

  subroutine pysave_real(variable, filename, n)
    real(dp), intent(in) :: variable
    integer, intent(in) :: n
    character(n), intent(in) :: filename

    call file_open(filename, n, 'float')
    write(fileunit, '(f12.7)') variable

    call file_close()
  end subroutine pysave_real

  subroutine pysave_real_1d(variable, filename, n)
    real(dp), intent(in) :: variable(:)
    integer, intent(in) :: n
    character(n), intent(in) :: filename
    integer :: c, i

    c = size(variable, 1)
    call file_open(filename, n, 'float')
    do i = 1, c
       write(fileunit, '(f12.7)') variable(c)
    end do
    call file_close()
  end subroutine pysave_real_1d

  subroutine pysave_real_2d(variable, filename, n)
    real(dp), intent(in) :: variable(:,:)
    integer, intent(in) :: n
    character(n), intent(in) :: filename
    character(20) :: FMT
    integer :: r, c, i

    r = size(variable, 1)
    c = size(variable, 2)
    write(FMT, *)

    call file_open(filename, n, 'float')
    do i = 1, r
       write(fileunit, '(' // adjustl(FMT) // 'f12.7)') variable(r, :)
    end do

    call file_close()
  end subroutine pysave_real_2d

  subroutine pysave_logical(variable, filename, n)
    logical, intent(in) :: variable
    integer, intent(in) :: n
    character(n), intent(in) :: filename

    call file_open(filename, n, 'logical')
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
 
  subroutine file_open(filename, n, template_name)
    integer, intent(in) :: n
    character(n), intent(in) :: filename    
    character(*), intent(in) :: template_name
    integer :: ioerr
    logical :: exist

    inquire(file=filename, exist=exist)
    if (exist) then
       open(fpy_newunit(fileunit), file=filename, status="old", position="append",  &
            action="write", iostat=ioerr)
    else
       open(fpy_newunit(fileunit), file=filename, status="new", action="write", iostat=ioerr)
       if (ioerr == 0) write(fileunit, *) '# <fortpy version="1" template="', template_name, '"></fortpy>'
    end if

    if (ioerr /= 0) then
       print *, "ERROR opening file for pysave in fortpy", ioerr
    end if
  end subroutine file_open

  !!<summary>Returns the number of values in the specified line assuming
  !!that they are separated by spaces or tabs.</summary>
  !!<parameter name="length">The number of characters in line.</parameter>
  !!<parameter name="line">The string for the line to count values in.</parameter>
  integer function fpy_value_count(line, length)
    integer, intent(in) :: length
    character(length), intent(in) :: line
    character(2) :: whitespace
    integer           :: success, i, indx, prev = 1, beginning = 1
    real              :: value

    !Initialize the whitespace array. We will cycle through all the characters
    !in the specified line looking for whitespace. Each time we find it, if the
    !character immediately preceding it was not whitespace, we have a value.
    whitespace = '  ' // char(9)
    fpy_value_count = 0

    do i = 1, length
       !indx will be zero if the current character is not a whitespace character.
       indx = index(whitespace, line(i:i))
       if (indx > 0 .and. prev == 0) then
          !We found the first whitespace after the end of a value we want.
          fpy_value_count = fpy_value_count + 1
       end if

       prev = indx
    end do

    !If the last value on the line ends right before \n, then we wouldn't have
    !picked it up; add an extra one.
    if (indx == 0) fpy_value_count = fpy_value_count + 1
  end function fpy_value_count

  !!<summary>Returns the number of lines in the file that aren't comments and
  !!the number of whitespace-separated values on the first non-comment line.</summary>
  !!<parameter name="filename">The name of the file to pass to open.</parameter>
  !!<parameter name="n">The number of characters in 'filename'.</parameter>
  !!<parameter name="commentchar">A single character which, when present at the start
  !!of a line designates it as a comment.</parameter>
  subroutine fpy_linevalue_count(filename, n, commentchar, nlines, nvalues)
    character(n), intent(in) :: filename
    integer, intent(in) :: n
    character(1), intent(in) :: commentchar
    integer, intent(out) :: nlines, nvalues
    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit, i
    character(500) :: line

    !Initialize the value for the result; if we get an error during the read, we
    !end the loop. It can be caused by badly formatted data or the EOF marker.
    nlines = 0
    nvalues = 0

    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
       do
          read(funit, "(A)", iostat=ioerr) line
          if (ioerr == 0) then
             cleaned = trim(adjustl(line))
             if (len(cleaned) .gt. 0 .and. cleaned(1:1) /= commentchar) then
                nlines = nlines + 1
                !We only need to get the number of values present in a line once.
                !We restrict the file structure to have rectangular arrays.
                if (nvalues == 0) then
                   nvalues = fpy_value_count(cleaned, len(cleaned))
                end if
             end if
          else
             exit
          end if
       end do
    end if
    close(funit)
  end subroutine fpy_linevalue_count

  !!<summary>Returns lowest i/o unit number not in use.</summary>
  !!<parameter name="unit">Out parameter that will contain the lowest i/o number.</parameter>
  integer function fpy_newunit(unit)
    integer, intent(out), optional :: unit
    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent
    integer :: n

    do n = nmin, nmax
       inquire(unit=n, opened=inuse)
       if (.not. inuse) then
          if (present(unit)) unit=n
          fpy_newunit = n
          return
       end if
    end do
    stop "newunit ERROR: available unit not found."
  end function fpy_newunit
end module fortpy
