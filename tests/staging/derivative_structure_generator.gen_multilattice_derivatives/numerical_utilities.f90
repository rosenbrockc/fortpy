MODULE numerical_utilities
use num_types
implicit none
private
public equal

! Overloaded procedure for comparing real types
   INTERFACE equal
      MODULE PROCEDURE equal_scalar, &
                       equal_rank1,  &
                       equal_rank2,  &
                       equal_rank3,  &
                       equal_rank1_rank0, &
                       equal_rank2_rank0, &
                       equal_rank2_real_int, &
                       equal_rank1_real_int, &
                       equal_scalar_real_int

   END INTERFACE

CONTAINS
!****************************************************************************************
! This function takes two real entities and compares them to see if they are equal
! within some tolerance. This prevents failed comparisons due to numbers that are "equal"
! but differ due to small differences arising from finite precision.
function equal_rank1(a, b, tolerance)
logical equal_rank1
real(dp) a(:), b(:), tolerance

equal_rank1 = .false.
if(all(abs(a - b) < tolerance)) equal_rank1 = .true.
end function equal_rank1

function equal_rank2(a, b, tolerance)
logical equal_rank2
real(dp) a(:,:), b(:,:), tolerance

equal_rank2 = .false.
if(all(abs(a - b) < tolerance)) equal_rank2 = .true.
end function equal_rank2

function equal_rank3(a, b, tolerance)
logical equal_rank3
real(dp) a(:,:,:), b(:,:,:), tolerance

equal_rank3 = .false.
if(all(abs(a - b) < tolerance)) equal_rank3 = .true.
end function equal_rank3

function equal_scalar(a, b, tolerance)
logical equal_scalar
real(dp) a, b, tolerance
equal_scalar = .false.
if(abs(a - b) < tolerance) equal_scalar = .true.
end function equal_scalar

function equal_scalar_real_int(a, b, tolerance)
logical equal_scalar_real_int
real(dp) a, tolerance
integer b

equal_scalar_real_int = .false.
if(abs(a - b) < tolerance) equal_scalar_real_int = .true.
end function equal_scalar_real_int

function equal_rank1_rank0(a, b, tolerance)
logical equal_rank1_rank0
real(dp) a(:), b, tolerance

equal_rank1_rank0 = .false.
if(all(abs(a - b) < tolerance)) equal_rank1_rank0 = .true.
end function equal_rank1_rank0

function equal_rank2_rank0(a, b, tolerance)
logical equal_rank2_rank0
real(dp) :: a(:,:), b, tolerance

equal_rank2_rank0 = .false.
if(all(abs(a - b) < tolerance)) equal_rank2_rank0 = .true.
end function equal_rank2_rank0

function equal_rank2_real_int(a, b, tolerance)
logical equal_rank2_real_int
real(dp) :: a(:,:), tolerance
integer :: b(:,:)
equal_rank2_real_int = .false.
if(all(abs(a - b) < tolerance)) equal_rank2_real_int = .true.
end function equal_rank2_real_int

function equal_rank1_real_int(a, b, tolerance)
logical equal_rank1_real_int
real(dp) :: a(:), tolerance
integer :: b(:)
equal_rank1_real_int = .false.
if(all(abs(a - b) < tolerance)) equal_rank1_real_int = .true.
end function equal_rank1_real_int

END MODULE numerical_utilities
