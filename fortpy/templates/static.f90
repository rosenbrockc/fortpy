!!<fortpy codeversion="__version__" />
!!<summary>Provides an interface for saving the values of multiple variable
!!types using a single call. Used as part of the FORTPY unit testing framework.</summary>
module fortpy
  implicit none
  private
  public pysave, fdp, fsp, fsi, fli, fpy_linevalue_count, fpy_newunit, fpy_read, &
       fpy_value_count, fpy_period_join_indices, fpy_linevalue_count_all, fpy_read_p, &
       fpy_read_f, fpy_vararray, autoclass_analyze, fpy_set_verbosity, fpy_verbose

  !!<member name="fileunit">I/O unit for the file to write the output to.</member>
  integer :: fileunit
  !!<member>Output/debugging verbosity for auxiliary reads/writes.</member>
  integer :: fpy_verbose
  
  !!<member name="seeded">Specifies whether the random number generator has
  !!already been seeded for this program execution.</member>
  logical :: seeded

  integer, parameter:: fdp=selected_real_kind(15,307)
  integer, parameter:: fsp=selected_real_kind(6,37)
  integer, parameter:: fsi=selected_int_kind(1) ! very short integer -10..10 range
  integer, parameter:: fli=selected_int_kind(18) ! Big integer -10^18..10^18 range

  !!<summary>Provides a single call interface for integer and real data types
  !!for single values, 1D and 2D arrays. Also handles outputting true/false for
  !!logical type variables.</summary>
  interface pysave
     module procedure __pysave__
  end interface pysave

  !!<summary>Reads values from a data file into a variable.</summary>
  interface fpy_read
     module procedure __fpy_read__
  end interface fpy_read

  !!<summary>Interface for fpy_read for variables with fixed dimensions (i.e. they
  !!don't have the 'pointer' or 'allocatable' attributes.</summary>
  interface fpy_read_f
     module procedure __fpy_read_f__
  end interface fpy_read_f
  
  !!<summary>Provides an interface for fpy_read for pointer-valued variables that
  !!need to be allocated before reading data.</summary>
  interface fpy_read_p
     module procedure __fpy_read_p__
  end interface fpy_read_p

  !!<summary>Provides a structure for variable length arrays that need to have their
  !!cartesian product taken.</summary>
  !!<usage>
  !!type(fpy_vararray) single
  !!allocate(single)
  !!single%init(array)
  !!</usage>
  type fpy_vararray
     !!<member name="items">The array of items to take cartesian product over.</member>
     !!<member name="length">The number of items in @CREF[this.items].</member>
     integer, pointer :: items(:)
     integer :: length
  contains
    procedure, public :: init => vararray_init
  end type fpy_vararray
contains
  !!<summary>Sets the global verbosity for the auxiliary reads/writes.</summary>
  subroutine fpy_set_verbosity(v)
    integer, intent(in) :: v
    fpy_verbose = v
    if (v > 0) write (*, *) "Fortpy F90 verbosity set to", v, "."
  end subroutine fpy_set_verbosity

  !!<summary>Analyzes the specified .fortpy.analysis file to determine the
  !!actual dimensionality of the data being read in via auto-class.</summary>
  subroutine autoclass_analyze(filename, analysis)
    character(len=*), intent(in) :: filename
    class(fpy_vararray), allocatable, intent(out) :: analysis(:)

    !!<local name="line">The integer dimensionality specified by a single line of the file.</local>
    !!<local name="ragvals">The number of values on each line of the file.</local>
    integer, allocatable :: line(:), ragvals(:)
    integer :: nlines, nvalues, funit, i

    call fpy_linevalue_count_all(filename, '#', nlines, ragvals)
    allocate(analysis(nlines))

    open(fpy_newunit(funit), file=filename)
    do i=1, nlines
      allocate(line(ragvals(i)))
      read(funit, *) line
      call analysis(i)%init(line, alloc=.true.)
      deallocate(line)
    end do
    close(funit)
  end subroutine autoclass_analyze
  
  !!<summary>Initializes the array items and length property.</summary>
  subroutine vararray_init(self, array, length, alloc)
    class(fpy_vararray) :: self
    integer, target, optional, intent(in) :: array(:)
    integer, optional, intent(in) :: length
    logical, optional, intent(in) :: alloc

    logical :: nalloc
    !We need to see if we are *copying* the array, or just referencing it.
    if (present(alloc)) then
       nalloc = alloc
    else
       nalloc = .false.
    end if

    if (present(array)) then
       if (nalloc) then
          allocate(self%items(size(array, 1)))
          self%items = array
       else
          self%items => array
       end if
       self%length = size(self%items, 1)
    else
       allocate(self%items(length))
       self%length = length
    end if
  end subroutine vararray_init

  __fxpy_read__
  __fxpy_read_f__
  __fxpy_read_p__
  
  !!<summary>Joins an array of integer indices as a period-separated string for
  !!concatenating to a file name.</summary>
  !!<parameter name="pslist">The period-separated result.</parameter>
  !!<parameter name="indices">An array of integer indices to join together.</parameter>
  !!<parameter name="n">The number of entries in the 'indices' array.</parameter>
  subroutine fpy_period_join_indices(pslist, indices, n)
    character(100), intent(out) :: pslist
    integer :: n
    integer, intent(in) :: indices(n)

    character(n*10) :: tempstr
    character(40) :: buffer
    integer :: i

    tempstr = ""
    pslist = ""
    do i=1, n
       write(buffer, '(I10)') indices(i)
       buffer = adjustl(buffer)
       if (i .eq. n) then
          tempstr = trim(pslist) // trim(buffer)
          pslist = tempstr
       else
          tempstr = trim(pslist) // trim(buffer) // "."
          pslist = tempstr
       end if
    end do
    pslist = trim(tempstr)
  end subroutine fpy_period_join_indices

  subroutine random_real(variable, min, max)
    real(fdp), intent(out) :: variable
    real(fdp), intent(in) :: min, max
    real(fdp) :: x
    !Since we only get real numbers in [0, 1], we have to set the value to be 
    !in the range manually.
    call random_number(x)    
    variable = min + x * (max - min)
  end subroutine random_real

  subroutine random_integer(variable, min, max)
    integer, intent(out) :: variable
    integer, intent(in) :: min, max
    real(fdp) :: x
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

  !!<summary>Escapes all the apostrophes in each word and surrounds it
  !!with a set of apostrophes.</summary>
  subroutine char_escape_word(word, escaped)
    character(len=*), intent(in) :: word
    character(1000), intent(inout) :: escaped
    integer :: ichar, ibuff

    escaped = "'"
    ibuff = 2
    do ichar=1, len(trim(adjustl(word)))
       if (word(ichar:ichar) .eq. "'") then
          escaped(ibuff:ibuff+1) = "''"
          ibuff = ibuff + 2
       else
          escaped(ibuff:ibuff) = word(ichar:ichar)
          ibuff = ibuff + 1
       end if
    end do
    escaped(ibuff:ibuff) = "'"
  end subroutine char_escape_word

  !!<summary>Writes a vector array of *escaped* strings as a single line
  !!to the currently open file unit.</summary>
  subroutine char_write_trimmed(variable)
    character(len=*), intent(in) :: variable(:)
    !!<local name="fmt">The compiled format string that includes an 'A#' directive
    !!for each string in the variable array where # is the length of the string.</local>
    !!<local name="word">Holds the length of the current array element as a string.</local>
    character(len=:), allocatable :: fmt
    character(10) :: word
    integer :: i
    character(1000) :: escaped(size(variable, 1))

    !First we need to escape all the quotes in each word in the array
    !so that they can be read in correctly by Fortran later.
    do i=1, size(variable, 1)
       call char_escape_word(variable(i), escaped(i))
    end do

    fmt = "("
    do i=1, size(escaped, 1)
       write (word, "(i10)") (len(trim(adjustl(escaped(i)))) + 2)
       fmt = fmt // "A" // trim(adjustl(word))
       if (i .lt. size(escaped, 1)) then
          fmt = fmt // ","
       end if
    end do
    fmt = fmt // ")"
    write(fileunit, fmt) escaped
  end subroutine char_write_trimmed
  
  __xpysave__  
  subroutine file_open(filename, n, template_name)
    integer, intent(in) :: n
    character(n), intent(in) :: filename    
    character(*), intent(in) :: template_name
    character(len=:), allocatable :: delim
    integer :: ioerr
    logical :: exist

    if (template_name .eq. "character") then
       delim = "APOSTROPHE"
    else
       delim = "NONE"
    end if
    
    inquire(file=filename, exist=exist)
    if (exist) then
       open(fpy_newunit(fileunit), file=filename, status="old", position="append",  &
            action="write", iostat=ioerr)
    else
       open(fpy_newunit(fileunit), file=filename, status="new", action="write", iostat=ioerr)
       if (ioerr == 0) write(fileunit, "(A)") '# <fortpy version="1" template="' // template_name // '"></fortpy>'
    end if

    if (ioerr /= 0) then
       print *, "ERROR opening file ", filename, " for pysave in fortpy", ioerr
    end if
  end subroutine file_open

  subroutine file_close()
    !This is just a one line routine, but if we decide to add cleanup later
    !it makes it easier to do it.
    close(fileunit)
  end subroutine file_close

  !!<summary>Returns the number of values in the specified line assuming
  !!that they are separated by spaces or tabs.</summary>
  !!<parameter name="length">The number of characters in line.</parameter>
  !!<parameter name="line">The string for the line to count values in.</parameter>
  !!<parameter name="ischar">When true, the type of data being read in is character
  !!so that whitespace is not the separator, but rather isolated apostrophes or quotes.</parameter>
  integer function fpy_value_count(line, length, ischar)
    integer, intent(in) :: length
    character(length), intent(in) :: line
    logical, optional, intent(in) :: ischar
    
    character(2) :: whitespace
    integer :: i, ichar, indx, prev, cindx(2)
    logical :: ischar_, isquote(2), fquote, fquote_set, qopen
    character :: cchar(2), quote='"', apost="'"

    if (present(ischar)) then
       ischar_ = ischar
    else
       ischar_ = .false.
    end if

    !Initialize the whitespace array. We will cycle through all the characters
    !in the specified line looking for whitespace. Each time we find it, if the
    !character immediately preceding it was not whitespace, we have a value.
    whitespace = ' ' // char(9)
    fquote_set = .false.
    qopen = .false.

    fpy_value_count = 0
    prev = -1
    ichar = 1

    do i = 1, length
       !indx will be zero if the current character is not a whitespace character.
       cchar(ichar) = line(i:i)
       if (ischar_) then
          !We need to identify whether the current character is a quote or apostrophe.
          !We also need to remember which one it is to handle the escaping of '' or "".
          if (cchar(ichar) .eq. apost) then
             isquote(ichar) = .false.
             cindx(ichar) = 1
             if (.not. fquote_set) then
                fquote = .false.
                fquote_set = .true.
             end if
          else
             if (cchar(ichar) .eq. quote) then
                isquote(ichar) = .true.
                cindx(ichar) = 1
                if (.not. fquote_set) then
                   fquote = .true.
                   fquote_set = .true.
                end if
             end if
          end if
          
          !Next, we can analyze whether this character and the last are both the same
          !character. If they are, it is as good as a non-match. If they aren't we
          !check to see if either is a quote and whether they match the first type of
          !quote in the line.
          if (cindx(1) .gt. 0) then
             !First, make sure that the previous character isn't the same since that would
             !be an escaped quote.
             if (ichar .gt. 1) then
                if (cchar(ichar) .eq. cchar(ichar-1)) then
                   !Reset the sequencer so that triple quotes count as an escaped quote
                   !followed by a non-escaped one.
                   indx = 0
                else
                   !If the previous one was a quote and this one is not, then we may have
                   !a match; in that case, make sure the quote matches the first one in
                   !the line.
                   if (cindx(1) .eq. 1 .and. cindx(2) .eq. 0) then
                      if (isquote(1) .eqv. fquote) then
                         indx = 1
                         qopen = (.not. qopen)
                      else
                         indx = 0
                      end if
                   end if
                end if
                ichar = 1
                cindx = 0
                cchar = ""
             else
                !we can't make decisions using an isolated character because of the
                !double character escaping standard in Fortran.
                ichar = ichar + 1
             end if
          else
             !Since neither was a quote, it is the equivalent of matching a non-whitespace
             !character type for the numeric counters.
             indx = 0
             cindx(1) = 0
             cchar(1) = ""
          end if
       else
          indx = index(whitespace, cchar(ichar))
       end if
       
       if (indx > 0 .and. prev == 0) then
          !We found the first whitespace/quote after the end of a value we want.
          if (ischar_) then
             if (qopen) fpy_value_count = fpy_value_count + 1
          else
             fpy_value_count = fpy_value_count + 1
          end if
       end if

       prev = indx
    end do

    if (ischar_) then
       !Check if the last character was a quote and of the same variety as the first
       !quote character on the line.
       if (ichar .eq. 1 .and. cindx(1) .eq. 1 .and. isquote(1) .eqv. fquote .and. qopen) &
            fpy_value_count = fpy_value_count + 1
    else
       !If the last value on the line ends right before \n, then we wouldn't have
       !picked it up; add an extra one.
       if (indx == 0) fpy_value_count = fpy_value_count + 1
    end if
  end function fpy_value_count

  !!<summary>Gets the lengths of ragged-array structured lines for each line
  !!in the specified data file.</summary>
  subroutine fpy_linevalue_count_all(filename, commentchar, nlines, nvalues, ischar)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    integer, intent(out) :: nlines
    integer, allocatable, intent(out) :: nvalues(:)
    logical, optional :: ischar

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit, i, firstnval
    character(250000) :: line
    logical :: ischar_

    if (present(ischar)) then
       ischar_ = ischar
    else
       ischar_ = .false.
    end if

    !Initialize the value for the result; if we get an error during the read, we
    !end the loop. It can be caused by badly formatted data or the EOF marker.
    call fpy_linevalue_count(filename, commentchar, nlines, firstnval)
    allocate(nvalues(nlines))
    nvalues = 0
    i = 0

    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
       do
          read(funit, "(A)", iostat=ioerr) line
          if (ioerr == 0) then
             cleaned = trim(adjustl(line))
             if (len(cleaned) .gt. 0) then
                if (cleaned(1:1) /= commentchar) then
                   i = i + 1
                   nvalues(i) = fpy_value_count(cleaned, len(cleaned), ischar_)
                end if
             end if
          else
             exit
          end if
       end do
    end if
    close(funit)
  end subroutine fpy_linevalue_count_all

  !!<summary>Returns the number of lines in the file that aren't comments and
  !!the number of whitespace-separated values on the first non-comment line.</summary>
  !!<parameter name="filename">The name of the file to pass to open.</parameter>
  !!<parameter name="n">The number of characters in 'filename'.</parameter>
  !!<parameter name="commentchar">A single character which, when present at the start
  !!of a line designates it as a comment.</parameter>
  subroutine fpy_linevalue_count(filename, commentchar, nlines, nvalues, ischar)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    integer, intent(out) :: nlines, nvalues
    logical, optional :: ischar
    
    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    character(250000) :: line
    logical :: ischar_, exists

    if (present(ischar)) then
       ischar_ = ischar
    else
       ischar_ = .false.
    end if

    !Initialize the value for the result; if we get an error during the read, we
    !end the loop. It can be caused by badly formatted data or the EOF marker.
    nlines = 0
    nvalues = 0
    inquire(file=filename, exist=exists)
    if (.not. exists) then
       write(*,*) "The file ", filename, " does not exist."
    end if
    
    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
       do
          read(funit, "(A)", iostat=ioerr) line
          if (ioerr == 0) then
             cleaned = trim(adjustl(line))
             if (abs(len(cleaned) - len(line)) < 10) write (*,*) "Number of characters in line ", len(cleaned), &
                  & " likely exceeds the hard-coded limit of ", len(line), " in file '", filename, "'."

             if (len(cleaned) .gt. 0) then
                if (cleaned(1:1) /= commentchar) then
                   nlines = nlines + 1
                   !We only need to get the number of values present in a line once.
                   !We restrict the file structure to have rectangular arrays.
                   if (nvalues == 0) then
                      nvalues = fpy_value_count(cleaned, len(cleaned), ischar_)
                   end if
                end if
             end if
          else
             if (ioerr .ne. -1) then
                write(*,*) "IO error (", ioerr, ") counting file lines and values. Found ", &
                     & nlines, " and ", nvalues, " in '", filename, "'."
             end if
             exit
          end if
       end do
    end if
    if (fpy_verbose > 0) write (*,*) "Found ", nlines, " lines and ", nvalues, " values in '", filename, "'."
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
