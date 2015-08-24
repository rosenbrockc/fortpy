!!<fortpy version="1.7.0" codeversion="1.7.0" />
!!<summary>Auto-generated auxiliary module exposing interfaces to save
!!instances of user-derived types.</summary>
module fpy_auxiliary
  use fortpy
  use iso_c_binding, only: c_ptr, c_associated, c_loc
  use enumeration_types
  use autoclass
  implicit none
  private
  
  public auxsave, auxread

  !!<summary>Provides a single call interface for saving user-derived data types
  !!for single values and arrays.</summary>
  interface auxsave
     module procedure auxsave_rotpermlist1d, auxsave_rotpermlist0d, auxsave_challenge0d, &
                       auxsave_oplist1d, auxsave_oplist0d, auxsave_simple1d, &
                       auxsave_simple0d
  end interface auxsave

  !!<summary>Provides a single call interface for reading user-derived data types
  !!for single values and arrays.</summary>
  interface auxread
     module procedure auxread_rotpermlist1d, auxread_rotpermlist1d_p, &
                       auxread_rotpermlist0d, auxread_challenge0d, auxread_oplist1d, &
                       auxread_oplist1d_p, auxread_oplist0d, auxread_simple1d, &
                       auxread_simple1d_p, auxread_simple0d
  end interface auxread

  !!<summary>Stores the prefix path and address of a pointer-type user-derived type
  !!variable that is being saved.</summary>
  type, public :: fpy_address
     !!<member name="prefix">The file prefix to the instances top-level variable.</member>
     !!<member name="address">64-bit adress of the pointer in memory.</member>
     character(len=:), allocatable :: prefix
     type(c_ptr) :: address
  end type fpy_address
contains
  subroutine auxsave_RotPermList1d_p(variable, pfolder, nested_)
    character(len=*), intent(in) :: pfolder
    type(RotPermList), pointer, intent(in) :: variable(:)
    logical, optional, intent(in) :: nested_

    character(len=:), allocatable :: folder
    integer :: fpy0
    character(100) :: pslist

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('enumeration_types.RotPermList', folder//'.fpy.type')
    end if

    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      call auxsave_RotPermList0d(variable(fpy0), folder//'-'//pslist, .true.)
    end do
  end subroutine auxsave_RotPermList1d_p

  subroutine auxsave_RotPermList1d(variable, folder)
    character(len=*), intent(in) :: folder
    type(RotPermList), allocatable, target, intent(in) :: variable(:)
    type(RotPermList), pointer :: lvar(:)
        
    lvar => variable
    call auxsave_RotPermList1d_p(lvar, folder)
  end subroutine auxsave_RotPermList1d

  subroutine auxsave_RotPermList0d(variable, pfolder, nested_)
    character(len=*), intent(in) :: pfolder
    type(RotPermList), intent(in) :: variable
    logical, optional, intent(in) :: nested_

    character(len=:), allocatable :: folder
    !V: dperms auto-class support variables

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('enumeration_types.RotPermList', folder//'.fpy.type')
    end if

    if (associated(variable%RotIndx)) then
      call pysave(variable%RotIndx, folder//'-rotindx')
    end if
    if (associated(variable%perm)) then
      call pysave(variable%perm, folder//'-perm')
    end if
    if (associated(variable%v)) then
      call pysave(variable%v, folder//'-v')
    end if
  end subroutine auxsave_RotPermList0d

  subroutine auxsave_challenge0d(variable, pfolder, nested_)
    character(len=*), intent(in) :: pfolder
    type(challenge), intent(in) :: variable
    logical, optional, intent(in) :: nested_

    character(len=:), allocatable :: folder
    !V: c auto-class support variables
    integer :: variable_R_ac1
    character(100) :: variable_R_acps2

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('autoclass.challenge', folder//'.fpy.type')
    end if

    do variable_R_ac1=1, size(variable%R, 1)
      call fpy_period_join_indices(variable_R_acps2, (/ variable_R_ac1 /), 1)
      if (allocated(variable%R)) then
        call auxsave(variable%R(variable_R_ac1), folder//'-r-'//trim(adjustl(variable_R_acps2)), .true.)
      end if
    end do
    call pysave(variable%F, folder//'-f')
  end subroutine auxsave_challenge0d

  subroutine auxsave_opList1d_p(variable, pfolder, nested_)
    character(len=*), intent(in) :: pfolder
    type(opList), pointer, intent(in) :: variable(:)
    logical, optional, intent(in) :: nested_

    character(len=:), allocatable :: folder
    integer :: fpy0
    character(100) :: pslist

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('enumeration_types.opList', folder//'.fpy.type')
    end if

    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      call auxsave_opList0d(variable(fpy0), folder//'-'//pslist, .true.)
    end do
  end subroutine auxsave_opList1d_p

  subroutine auxsave_opList1d(variable, folder)
    character(len=*), intent(in) :: folder
    type(opList), allocatable, target, intent(in) :: variable(:)
    type(opList), pointer :: lvar(:)
        
    lvar => variable
    call auxsave_opList1d_p(lvar, folder)
  end subroutine auxsave_opList1d

  subroutine auxsave_opList0d(variable, pfolder, nested_)
    character(len=*), intent(in) :: pfolder
    type(opList), intent(in) :: variable
    logical, optional, intent(in) :: nested_

    character(len=:), allocatable :: folder
    !V: tmpOp auto-class support variables

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('enumeration_types.opList', folder//'.fpy.type')
    end if

    if (associated(variable%shift)) then
      call pysave(variable%shift, folder//'-shift')
    end if
    if (associated(variable%rot)) then
      call pysave(variable%rot, folder//'-rot')
    end if
  end subroutine auxsave_opList0d

  subroutine auxsave_simple1d_p(variable, pfolder, nested_)
    character(len=*), intent(in) :: pfolder
    type(simple), pointer, intent(in) :: variable(:)
    logical, optional, intent(in) :: nested_

    character(len=:), allocatable :: folder
    integer :: fpy0
    character(100) :: pslist

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('autoclass.simple', folder//'.fpy.type')
    end if

    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      call auxsave_simple0d(variable(fpy0), folder//'-'//pslist, .true.)
    end do
  end subroutine auxsave_simple1d_p

  subroutine auxsave_simple1d(variable, folder)
    character(len=*), intent(in) :: folder
    type(simple), allocatable, target, intent(in) :: variable(:)
    type(simple), pointer :: lvar(:)
        
    lvar => variable
    call auxsave_simple1d_p(lvar, folder)
  end subroutine auxsave_simple1d

  subroutine auxsave_simple0d(variable, pfolder, nested_)
    character(len=*), intent(in) :: pfolder
    type(simple), intent(in) :: variable
    logical, optional, intent(in) :: nested_

    character(len=:), allocatable :: folder
    !V: R auto-class support variables

    if (present(nested_)) then
      folder = pfolder
    else
      call system("mkdir -p "//pfolder)
      folder = pfolder//'_'
      call pysave('autoclass.simple', folder//'.fpy.type')
    end if

    if (allocated(variable%A)) then
      call pysave(variable%A, folder//'-a')
    end if
    if (associated(variable%B)) then
      call pysave(variable%B, folder//'-b')
    end if
  end subroutine auxsave_simple0d


  subroutine auxread_RotPermList1d_p(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(RotPermList), pointer, intent(inout) :: variable(:)
    type(RotPermList), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), optional, intent(in) :: rstack(:)
    
    type(RotPermList), pointer :: pointer_stack(:)
    integer :: fpy0
    character(100) :: pslist
    type(fpy_address), allocatable :: stack(:)
    integer, allocatable :: drange(:)
    type(RotPermList), pointer :: tempvar

    if (present(multi_stack) .and. present(rstack)) then
      stack = rstack
      pointer_stack => multi_stack
    else
      call fpy_read_address(stack, folder)
      allocate(pointer_stack(size(stack)))
    end if

    call fpy_get_drange(drange, folder//'-', (/ 'RotIndx', 'perm   ', 'v      ' /))
    if (.not. allocated(drange)) return
    allocate(variable(drange(1)))
    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      tempvar => variable(fpy0)
      call auxread_RotPermList0d(tempvar, folder//'-'//trim(adjustl(pslist)), pointer_stack, stack)
    end do
  end subroutine auxread_RotPermList1d_p

  subroutine auxread_RotPermList1d(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(RotPermList), allocatable, target, intent(inout) :: variable(:)
    type(RotPermList), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), optional, intent(in) :: rstack(:)
        
    type(RotPermList), pointer :: lvar(:)
    lvar => variable
    call auxread_RotPermList1d_p(lvar, folder, multi_stack, rstack)
    variable = lvar
  end subroutine auxread_RotPermList1d

  subroutine auxread_RotPermList0d(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(RotPermList), intent(inout) :: variable
    type(RotPermList), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), allocatable, intent(in), optional :: rstack(:)

    character(len=:), allocatable :: lfolder

    if (present(multi_stack) .and. present(rstack)) then
       lfolder = folder
    else
       lfolder = folder//'_'
    end if

    call fpy_read_p(lfolder//'-rotindx', '#', variable%rotindx)
    call fpy_read_p(lfolder//'-perm', '#', variable%perm)
    call fpy_read_p(lfolder//'-v', '#', variable%v)
  end subroutine auxread_RotPermList0d

  subroutine auxread_challenge0d(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(challenge), intent(inout) :: variable
    type(challenge), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), allocatable, intent(in), optional :: rstack(:)

    character(len=:), allocatable :: lfolder

    if (present(multi_stack) .and. present(rstack)) then
       lfolder = folder
    else
       lfolder = folder//'_'
    end if

    allocate(variable%r(0))
    call auxread_simple1d(variable%r, lfolder//'-r')
    call fpy_read_f(lfolder//'-f', '#', variable%f)
  end subroutine auxread_challenge0d

  subroutine auxread_opList1d_p(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(opList), pointer, intent(inout) :: variable(:)
    type(opList), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), optional, intent(in) :: rstack(:)
    
    type(opList), pointer :: pointer_stack(:)
    integer :: fpy0
    character(100) :: pslist
    type(fpy_address), allocatable :: stack(:)
    integer, allocatable :: drange(:)
    type(opList), pointer :: tempvar

    if (present(multi_stack) .and. present(rstack)) then
      stack = rstack
      pointer_stack => multi_stack
    else
      call fpy_read_address(stack, folder)
      allocate(pointer_stack(size(stack)))
    end if

    call fpy_get_drange(drange, folder//'-', (/ 'shift', 'rot  ' /))
    if (.not. allocated(drange)) return
    allocate(variable(drange(1)))
    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      tempvar => variable(fpy0)
      call auxread_opList0d(tempvar, folder//'-'//trim(adjustl(pslist)), pointer_stack, stack)
    end do
  end subroutine auxread_opList1d_p

  subroutine auxread_opList1d(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(opList), allocatable, target, intent(inout) :: variable(:)
    type(opList), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), optional, intent(in) :: rstack(:)
        
    type(opList), pointer :: lvar(:)
    lvar => variable
    call auxread_opList1d_p(lvar, folder, multi_stack, rstack)
    variable = lvar
  end subroutine auxread_opList1d

  subroutine auxread_opList0d(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(opList), intent(inout) :: variable
    type(opList), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), allocatable, intent(in), optional :: rstack(:)

    character(len=:), allocatable :: lfolder

    if (present(multi_stack) .and. present(rstack)) then
       lfolder = folder
    else
       lfolder = folder//'_'
    end if

    call fpy_read_p(lfolder//'-shift', '#', variable%shift)
    call fpy_read_p(lfolder//'-rot', '#', variable%rot)
  end subroutine auxread_opList0d

  subroutine auxread_simple1d_p(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(simple), pointer, intent(inout) :: variable(:)
    type(simple), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), optional, intent(in) :: rstack(:)
    
    type(simple), pointer :: pointer_stack(:)
    integer :: fpy0
    character(100) :: pslist
    type(fpy_address), allocatable :: stack(:)
    integer, allocatable :: drange(:)
    type(simple), pointer :: tempvar

    if (present(multi_stack) .and. present(rstack)) then
      stack = rstack
      pointer_stack => multi_stack
    else
      call fpy_read_address(stack, folder)
      allocate(pointer_stack(size(stack)))
    end if

    call fpy_get_drange(drange, folder//'-', (/ 'A', 'B' /))
    if (.not. allocated(drange)) return
    allocate(variable(drange(1)))
    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      tempvar => variable(fpy0)
      call auxread_simple0d(tempvar, folder//'-'//trim(adjustl(pslist)), pointer_stack, stack)
    end do
  end subroutine auxread_simple1d_p

  subroutine auxread_simple1d(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(simple), allocatable, target, intent(inout) :: variable(:)
    type(simple), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), optional, intent(in) :: rstack(:)
        
    type(simple), pointer :: lvar(:)
    lvar => variable
    call auxread_simple1d_p(lvar, folder, multi_stack, rstack)
    variable = lvar
  end subroutine auxread_simple1d

  subroutine auxread_simple0d(variable, folder, multi_stack, rstack)
    character(len=*), intent(in) :: folder
    type(simple), intent(inout) :: variable
    type(simple), pointer, optional, intent(in) :: multi_stack(:)
    type(fpy_address), allocatable, intent(in), optional :: rstack(:)

    character(len=:), allocatable :: lfolder

    if (present(multi_stack) .and. present(rstack)) then
       lfolder = folder
    else
       lfolder = folder//'_'
    end if

    call fpy_read(lfolder//'-a', '#', variable%a)
    call fpy_read_p(lfolder//'-b', '#', variable%b)
  end subroutine auxread_simple0d


  !!<summary>Adds the specified value to the end of the integer-valued list.</summary>
  !!<parameter name="list">The integer-valued list to append to.</parameter>
  !!<parameter name="appendage">The extra value to append to the list.</parameter>
  subroutine append_address(list, appendage)
    type(fpy_address), allocatable, intent(inout) :: list(:)
    type(fpy_address), intent(in) :: appendage
    type(fpy_address), allocatable :: templist(:)

    allocate(templist(size(list,1)+1))
    templist(1:size(list,1)) = list
    templist(size(list,1)+1) = appendage

    call move_alloc(templist, list)
  end subroutine append_address

  !!<summary>Returns the integer index of the specified address in the stack.</summary>
  !!<parameter name="stack">Array of 64-bit addresses to pointers.</parameter>
  !!<parameter name="avalue">Address to search for in the stack.</parameter>
  integer function address_loc(stack, avalue, qprefix)
    type(fpy_address), intent(in) :: stack(:), avalue
    logical, optional, intent(in) :: qprefix
    integer :: i
    address_loc = -1
    do i=1, size(stack, 1)
       if ((.not.present(qprefix)) .and. c_associated(stack(i)%address, avalue%address)) then
          address_loc = i
          exit
       else if (present(qprefix) .and. stack(i)%prefix .eq. avalue%prefix) then
          address_loc = i
          exit
       end if       
    end do
  end function address_loc

  !!<parameter name="drange">The dimensionality of the variable at the current prefix context.</parameter>
  !!<parameter name="prefix">The file context (ending in -) of the variable to check data range on.</parameter>
  !!<parameter name="members">The name of a non-array variable that terminates the variable chain.</parameter>
  subroutine fpy_get_drange(drange, prefix, members)
    integer, allocatable, intent(out) :: drange(:)
    character(len=*), intent(in) :: prefix, members(:)

    integer :: D, i, j, k
    logical :: exists, memexists
    character(len=:), allocatable :: catstr
    character(50) :: istr
    
    !First, check the dimensionality of the range by checking for a member at the first position.
    exists = .true.
    memexists = .false.
    i = 0
    D = 0
    catstr = ''

    do while (i .lt. 7)
       if (len(catstr) .gt. 0) then
          catstr = catstr//'.1'
       else
          catstr = '1'
       end if

       memexists = .false.
       j = 1
       do while (j .le. size(members) .and. .not. memexists)
          inquire(file=prefix//catstr//'-'//trim(members(j)), exist=memexists)
          j = j+1
       end do
       
       i = i + 1
       if (memexists) then
          D = i
       end if
    end do
    if (D .eq. 0) return
    allocate(drange(D))

    do k=1, D
       i = 1
       exists = .true.
       do while(exists)
          write (istr, *) i
          catstr = ''
          
          do j=1, k
             if (len(catstr) .gt. 0) catstr = catstr//'.'
             
             if (j .lt. k) then
                catstr = catstr//'1'
             else
                catstr = catstr//trim(adjustl(istr))
             end if
          end do
          do j=k, D-1
             catstr = catstr//'.1'
          end do

          memexists = .false.
          j = 1
          do while (j .le. size(members) .and. .not. memexists)
             inquire(file=prefix//catstr//'-'//trim(members(j)), exist=memexists)
             j = j+1
          end do

          exists = memexists
          if (exists) then
             i = i+1
          else
             i = i-1
          end if
       end do
       drange(k) = i
    end do    
  end subroutine fpy_get_drange
  
  !!<summary>Returns the integer index of the address in @CREF[param.stack] that
  !!has the same prefix.</summary>
  integer function fpy_address_index(stack, prefix)
    type(fpy_address), intent(in) :: stack(:)
    character(len=*), intent(in) :: prefix
    integer :: i, ind
    character(len=:), allocatable :: tprefix

    fpy_address_index = 0
    ind = index(prefix, '/', .true.)
    if (ind .gt. 0) then
       tprefix = prefix(ind+1:)
    else
       tprefix = prefix
    end if

    do i=1, size(stack)
       if (stack(i)%prefix .eq. tprefix) then
          fpy_address_index = i
          exit
       end if
    end do
  end function fpy_address_index
  
  !!<summary>Saves the list of prefixes in order for deserializing later.</summary>
  !!<parameter name="stack">The stack to save to file.</parameter>
  !!<parameter name="folder">The path to the variable's saving folder.</parameter>
  subroutine fpy_save_addresses(stack, folder)
    type(fpy_address), intent(in) :: stack(:)
    character(len=*), intent(in) :: folder
    !!<local name="savelist">The list of prefixes to save, in order.</local>
    character(100) :: savelist(size(stack, 1))
    integer :: i

    do i=1, size(stack, 1)
       savelist(i) = stack(i)%prefix
    end do
    call pysave(savelist, folder//'.fpy.address')
  end subroutine fpy_save_addresses

  !!<summary>Restores the stack of fpy_addresses so that a variable's folder can be
  !!deserialized back via aux_read.</summary>
  subroutine fpy_read_address(stack, folder)
    type(fpy_address), intent(out), allocatable :: stack(:)
    character(len=*), intent(in) :: folder
    character(100), allocatable :: savelist(:)
    integer :: i
    logical :: fpy_success
    
    call fpy_read(folder//'.fpy.address', '#', savelist, fpy_success)
    if (.not. fpy_success) then
       allocate(stack(0))
    else
       allocate(stack(size(savelist)))
       do i=1, size(savelist)
          stack(i)%prefix = savelist(i)
       end do
    end if
  end subroutine fpy_read_address
  
  !!<summary>Returns the 64-bit integer address of the pointer at @CREF[param.cloc].</summary>
  !!<parameter name="ploc">The object returned by calling the @CREF[loc] interface on
  !!the variable.</parameter>
  !!<parameter name="prefix">The file location prefix for the variable being located.</parameter>
  type(fpy_address) function fpy_get_address(ploc, prefix)
    type(c_ptr), intent(in) :: ploc
    character(len=*), intent(in) :: prefix
    integer :: ind

    fpy_get_address%address = ploc
    ind = index(prefix, '/', .true.)
    if (ind .gt. 0) then
       fpy_get_address%prefix = prefix(ind+1:)
    else
       fpy_get_address%prefix = "_"
    end if       
  end function fpy_get_address    
end module fpy_auxiliary
