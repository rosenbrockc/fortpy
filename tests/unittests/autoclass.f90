!!<summary>Tests the autoclass unit-test generation and comparision
!!checking for nested user-derived types.</summary>
module autoclass
  use num_types
  implicit none

  type simple
     integer, allocatable :: A(:,:)
     real(dp), pointer :: B(:)
  end type simple

  type challenge
     type(simple), allocatable :: R(:)
     integer :: F(2,2)
  end type challenge
contains
  !!<summary>We need something to hook the testing framework to.</summary>
  !!<parameter name="c" regular="true">The challenging user-type to test with.</parameter>
  subroutine do_nothing(c)
    type(challenge), intent(in) :: c
    integer :: local = 0
    local = local + 1
  end subroutine do_nothing
end module autoclass
