MODULE rational_mathematics
  use num_types
  use vector_matrix_utilities
  use utilities_module, only: ralloc
implicit none
private
public gcd, SmithNormalForm, HermiteNormalForm, is_a_rational_in_range,&
       get_rationals_in_range

! Overloaded procedure for computing the greatest common denominator
INTERFACE gcd
   MODULE PROCEDURE gcd_2ints, gcd_rank1, gcd_3ints, gcd_4ints
END INTERFACE

CONTAINS
!*******************************************************************************
! This routine takes an integer 3x3 matrix and computes its Smith Normal Form
subroutine SmithNormalForm(H,A,M,B)
integer, intent(in) :: H(3,3) ! Input matrix
! Smith normal form (M), left (A) & right (B) transforms
integer, intent(out), dimension(3,3) :: M, A, B 

integer i, minm, maxidx, minidx, multiple, j, nondividx(2), check(9)
logical Ldiv(2,2)
integer itCnt


if(determinant(H)<1) stop "SmithNormalForm routine failed because the input matrix had a negative determinant"
A = 0; B = 0; M = H ! M starts out as H, the input matrix
forall(i=1:3); A(i,i) = 1; B(i,i) = 1; end forall ! A & B = identity

j=1
itCnt = 0

do ! Keep doing steps (1) and (2) until all elements in j-th row and column are
   ! zero except the one on the diagonal. When j == 1, if the (1,1) element doesn't
   ! divide every other non-zero element left in the matrix at the end, then add the
   ! offending row to row 1 and try again.
   !(1) Use row operations to zero out the column
   ! Divide the row with the smallest value into the largest

   itCnt = itCnt + 1
   if (itCnt>=100) stop "ERROR bad programming in SmithNormalForm"

   do while (count(M(:,j)/=0) > 1) ! Keep going until only 1 non-zero element in column j
      call get_minmax_indices(M(:,j),minidx,maxidx)
      minm = M(minidx,j)
      ! Subtract a multiple of the row containing the smallest element from
      ! the row containing the largest element
      multiple = M(maxidx,j)/minm
      M(maxidx,:) = M(maxidx,:) - multiple*M(minidx,:)
      A(maxidx,:) = A(maxidx,:) - multiple*A(minidx,:)
      if (any(matmul(matmul(A,H),B)/=M)) stop "ROW: Transformation matrices didn't work"
   enddo ! End of step 1


   if (M(j,j) == 0) then; call swap_row(A,M,j)
   endif
   if (any(matmul(matmul(A,H),B)/=M)) stop "ROWSWAP: Transformation matrices didn't work"
   
   !(2) Use column operations to zero out first row
   ! Divide the colum with the smallest value into the largest
   do while (count(M(j,:)/=0) > 1) ! Keep going until only 1 non-zero element in row 1
      !call printAMB(A,M,B)
      call get_minmax_indices(M(j,:),minidx,maxidx)
      minm = M(j,minidx)
      ! Subtract a multiple of the column containing the smallest element from
      ! the row containing the largest element
      multiple = M(j,maxidx)/minm ! Factor to multiply by
      M(:,maxidx) = M(:,maxidx) - multiple*M(:,minidx)
      B(:,maxidx) = B(:,maxidx) - multiple*B(:,minidx)
      if (any(matmul(matmul(A,H),B)/=M)) stop "COLS: Transformation matrices didn't work"
   enddo ! End of step 2


   if (M(j,j)<0) then ! Change signs
      M(:,j) = -M(:,j); B(:,j) = -B(:,j)
   elseif (M(j,j) == 0) then;call swap_column(M,B,j)
   endif
   if (count(M(j,:)/=0) > 1 .or. count(M(:,j)/=0) > 1) cycle

   if (any(matmul(matmul(A,H),B)/=M)) stop "COLSWAP: Transformation matrices didn't work"

   Ldiv = mod(M(2:,2:),M(1,1)) == 0
   if (j==1 .and. any(Ldiv .eqv. .false.)) then! Add the offending row to row 1 
!      nondividx = maxloc(mod(M(2:,2:),M(1,1))) ! Find one of the elements that isn't 
      nondividx = maxloc(abs( mod(M(2:,2:),M(1,1)) )) ! Find one of the elements that isn't 
      M(1,:) = M(1,:) + M(nondividx(1)+1,:)    ! divided by the diagonal element, and 
      A(1,:) = A(1,:) + A(nondividx(1)+1,:)    ! add the row it's in to row 1
      cycle ! Go back to the top of the outer do loop
   endif
   if (j==2) then 
      if (mod(M(3,3),M(2,2))/=0) then
         M(2,:) = M(2,:) + M(3,:)
         A(2,:) = A(2,:) + A(3,:)
         cycle
      endif
   else
      j = 2;cycle;endif ! Start row/col 2
   ! Try again if the matrix still isn't diagonal
   if (j == 2 .and. (M(3,2)/=0 .or. M(2,3)/=0)) then; cycle; endif
   exit ! We should be done if we hit this point
enddo
if (M(3,3)<0) then ! Change signs
   M(:,3) = -M(:,3); B(:,3) = -B(:,3);
endif
if (any(matmul(matmul(A,H),B)/=M)) stop "END: Transformation matrices didn't work"
check = reshape(M,(/9/))
if (any(check((/2,3,4,6,7,8/))/=0)) stop "Not diagonal"
if (mod(M(2,2),M(1,1))/=0 .or. mod(M(3,3),M(2,2))/=0) stop "SNF conditions not met"
ENDSUBROUTINE SmithNormalForm

!***************************************************************************************************
! Find the Hermite normal form of a a given integer matrix. Similar to the SNF finder above but a
! a little simpler. This routine is not very elegant, just brute force. Don't be disappointed.
subroutine HermiteNormalForm(S,H,B)
integer, intent(in) :: S(3,3) ! Input matrix
integer, intent(out), dimension(3,3) :: H, B ! HNF output and transformation matrix

integer i, minm, maxidx, minidx, multiple, j,  check(9), tempcol(3)

if (determinant(S) == 0) stop "Singular matrix passed to HNF routine"
B = 0; H = S ! H starts out as S, the input matrix
forall(i=1:3); B(i,i) = 1; end forall ! B = identity

do ! Keep doing column operations until all elements in row 1 are
   ! zero except the one on the diagonal. 
   ! Divide the column with the smallest value into the largest
   do while (count(H(1,:)/=0) > 1) ! Keep going until only zeros beyond first element
      call get_minmax_indices(H(1,:),minidx,maxidx)
      minm = H(1,minidx)
      ! Subtract a multiple of the column containing the smallest element from
      ! the row containing the largest element
      multiple = H(1,maxidx)/minm ! Factor to multiply by
      H(:,maxidx) = H(:,maxidx) - multiple*H(:,minidx)
      B(:,maxidx) = B(:,maxidx) - multiple*B(:,minidx)
      if (any(matmul(S,B)/=H)) stop "COLS: Transformation matrices didn't work"
   enddo ! End of step row 1
   if (H(1,1) == 0) call swap_column(H,B,1) ! swap columns if (1,1) is zero
   if (H(1,1)<0)then; H(:,1) = -H(:,1); B(:,1) = -B(:,1) ! Change sign if (1,1) is negative
   endif
   if (count(H(1,:)/=0) > 1) stop "Didn't zero out the rest of the row"
   if (any(matmul(S,B)/=H)) stop "COLSWAP: Transformation matrices didn't work"
   exit
enddo
! Now work on element H(2,3)
do
   do while (H(2,3)/=0)
      if (H(2,2) == 0) then
         tempcol = H(:,2); H(:,2) = H(:,3); H(:,3) = tempcol 
         tempcol = B(:,2); B(:,2) = B(:,3); B(:,3) = tempcol 
         if (H(2,3) == 0) exit
      endif
      if (abs(H(2,3))<abs(H(2,2))) then; maxidx = 2; minidx = 3
      else; maxidx = 3; minidx = 2; endif
      multiple = H(2,maxidx)/H(2,minidx)
      H(:,maxidx) = H(:,maxidx) - multiple*H(:,minidx)
      B(:,maxidx) = B(:,maxidx) - multiple*B(:,minidx)
      if (any(matmul(S,B)/=H)) stop "COLS: Transformation matrices didn't work"
   enddo
   if (H(2,2) == 0) then
      tempcol = H(:,2); H(:,2) = H(:,3); H(:,3) = tempcol 
   endif
   if (H(2,2)<0) then ! Change signs
      H(:,2) = -H(:,2); B(:,2) = -B(:,2); endif
   if (H(2,3)/=0) stop "Didn't zero out last element"

   if (any(matmul(S,B)/=H)) stop "COLSWAP: Transformation matrices didn't work"
   exit
enddo
if (H(3,3)<0) then ! Change signs
   H(:,3) = -H(:,3); B(:,3) = -B(:,3);
endif
check = reshape(H,(/9/))
if (any(check((/4,7,8/))/=0)) stop "Not lower triangular"
if (any(matmul(S,B)/=H)) stop "END PART1: Transformation matrices didn't work"
 
! Now that the matrix is in lower triangular form, make sure the lower off-diagonal 
! elements are non-negative but less than the diagonal elements
do while (H(2,2) <= H(2,1) .or. H(2,1)<0)
      if (H(2,2) <= H(2,1)) then; multiple = 1
      else; multiple = -1; endif
      H(:,1) = H(:,1) - multiple*H(:,2)
      B(:,1) = B(:,1) - multiple*B(:,2)
enddo
do j = 1,2
   do while (H(3,3) <= H(3,j) .or. H(3,j)<0)
         if (H(3,3) <= H(3,j)) then; multiple = 1
         else; multiple = -1; endif
         H(:,j) = H(:,j) - multiple*H(:,3)
         B(:,j) = B(:,j) - multiple*B(:,3)
   enddo
enddo
if (any(matmul(S,B)/=H)) stop "END: Transformation matrices didn't work"
check = reshape(H,(/9/))
if (any(check((/4,7,8/))/=0)) stop "Not lower triangular"
if (any(check((/2,3,6/))<0)) stop "Negative elements in lower triangle"
if (check(2) > check(5) .or. check(3) > check(9) .or. check(6) > check(9)) stop "Lower triangular elements bigger than diagonal"
ENDSUBROUTINE HermiteNormalForm

!***************************************************************************************************
! Support routines for the Smith and Hermite normal form finders. 
subroutine swap_row(A,M,k) ! Swap rows of M (and A)
integer, intent(inout) :: M(3,3), A(3,3)
integer, intent(in) :: k
integer tmpRow(3), maxidx(1)
 
maxidx = maxloc(abs(M(k:,k)))+k-1  ! find index of the non-zero element in col k
tmpRow = A(k,:); A(k,:) = A(maxidx(1),:); A(maxidx(1),:) = tmpRow
tmpRow = M(k,:); M(k,:) = M(maxidx(1),:); M(maxidx(1),:) = tmpRow
endsubroutine swap_row

subroutine swap_column(M,B,k) ! Swap columns of M (and B)
integer, intent(inout) :: M(3,3), B(3,3)
integer, intent(in) :: k
integer tmpCol(3), maxidx(1)
 
maxidx = maxloc(abs(M(k,k:)))+k-1 ! find index of the non-zero element in row k
tmpCol = B(:,k); B(:,k) = B(:,maxidx(1)); B(:,maxidx(1)) = tmpCol
tmpCol = M(:,k); M(:,k) = M(:,maxidx(1)); M(:,maxidx(1)) = tmpCol
endsubroutine swap_column

subroutine printAMB(A,M,B)
integer, intent(in), dimension(3,3) :: A,M,B
integer i
do i = 1,3
   write(*,'(3(3x,3i4))') A(i,:),M(i,:),B(i,:);enddo;print *
endsubroutine printAMB

subroutine printHB(H,B)
integer, intent(in), dimension(3,3) :: H,B
integer i
do i = 1,3
   write(*,'(2(3x,3i4))') H(i,:),B(i,:);enddo;print *
endsubroutine printHB

subroutine get_minmax_indices(invec,min,max)
integer, intent(in) :: invec(3)
integer, intent(out) :: min, max

integer :: tmpmin(1), tmpmax(1), vec(3)
vec = abs(invec)
tmpmin = minloc(vec,vec>0)
! Search from the right for the max so it will be different from the minimum
! even if the min and max are the same value
tmpmax = 4 - maxloc(vec(3:1:-1),vec(3:1:-1)>0)
min = tmpmin(1)
max = tmpmax(1)
endsubroutine get_minmax_indices

!*****************************************************************************
! This function finds the greatest common denominator of several integers
! ****************************************************************************
! 
! This case works for two integers, given as separate arguments
function gcd_2ints(x1, x2) result(divisor)
integer, intent(in) :: x1, x2
integer divisor

integer a, b
a = abs(x1); b = abs(x2) ! Make sure inputs are positive
if (b>a) call swap(a,b)  

do ! Keep dividing a by b, until one of them is zero
   if (b>a) call swap(a,b) ! Keep the bigger number in a's place
   if (b == 0) exit ! we're done when b == 0
   a = mod(a,b) ! Divide a by b and keep only the remainder
enddo
divisor = a

contains
  subroutine swap(x,y) ! Swap two values
    integer x,y,tmp
    tmp = x; x = y; y = tmp
  endsubroutine swap
end function gcd_2ints

! This case works on a list of integers (a rank-1 array)
function gcd_rank1(x) result(divisor)
integer, intent(in) :: x(:)
integer divisor

integer a(size(x)), N, indx(1), big2

N = size(x); a = abs(x)
if (any(a<0)) stop "GCD requires non-negative integers"
do ! Divide the biggest number by the second biggest until 
   ! the second biggest is zero
   indx = maxloc(a)  ! Find the location of the biggest number
   if (all(a == a(indx(1)))) then ! check if all numbers are the same
      big2 = a(indx(1))
   else   ! The "real" around 'a' is a workaround for a problem in the Absoft compiler
      big2 = maxval(real(a),mask=a < a(indx(1))) ! Find the size of the 2nd biggest number
   endif
   if (big2 == 0) exit
   a(indx(1)) = mod(a(indx(1)),big2)
enddo
divisor = a(indx(1)) ! The divisor is the number left when every other member==0
endfunction gcd_rank1

! This case works on 3 integers, not in an array
function gcd_3ints(x1,x2,x3)
integer, intent(in) :: x1,x2,x3
integer gcd_3ints
gcd_3ints = gcd_rank1((/x1,x2,x3/))
end function gcd_3ints

! This case works on 4 integers, not in an array
function gcd_4ints(x1,x2,x3,x4)
integer, intent(in) :: x1,x2,x3,x4
integer gcd_4ints
gcd_4ints = gcd_rank1((/x1,x2,x3,x4/))
end function gcd_4ints

! this function checks to see if there is a rational number in a range
! (assume a 1-D range. This isn't general for ternaries, etc.)
function is_a_rational_in_range(range,n,eps_)
integer range(3)
integer n
real(dp), optional :: eps_
logical is_a_rational_in_range
integer j2
real(dp) eps,left,right

is_a_rational_in_range = .false.

if (present(eps_)) then
   eps = eps_
else
   eps = 1e-12_dp
endif
left  = real(range(1),dp)/range(3) - eps
right = real(range(2),dp)/range(3) + eps
do j2 = 0, n ! numerators
   if      (real(j2,dp)/n > left &
      .and. real(j2,dp)/n < right) then ! there is a rational in the range
      is_a_rational_in_range = .true.
      return
   endif
enddo
! If the loops end, then the if condition was never met and
! there isn't a rational number in the range
endfunction

! This routine generates all of the rational numbers of a given
! numerator that lie within a specified range.
! (assume a 1-D range. This isn't general for ternaries, etc.)
subroutine get_rationals_in_range(range,n,numerators,eps_)
integer range(3)
integer n
integer,   pointer :: numerators(:)
real(dp), optional :: eps_
integer j, cR
real(dp) eps,left,right

nullify(numerators)

if (present(eps_)) then
   eps = eps_
else
   eps = 1e-12_dp
endif
left  = real(range(1),dp)/range(3) - eps
right = real(range(2),dp)/range(3) + eps
allocate(numerators(n+1))
numerators = 0

cR = 0
do j = 0, n ! numerators
   if      (real(j,dp)/n > left &
      .and. real(j,dp)/n < right) then ! there is a rational in the range
      cR = cR + 1
      numerators(cR) = j
   endif
enddo
if (cR==0) stop "get_rationals_in_range didn't find any in the range"
numerators => ralloc(numerators,cR)
end subroutine get_rationals_in_range

END MODULE rational_mathematics
