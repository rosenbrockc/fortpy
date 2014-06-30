!!<summary>Auto-generated unit test for derivative_structure_generator.gen_multilattice_derivatives
!!using FORTPY. Generated on 2014-06-27 12:49:15.713785.
!!Original enumerations tests of basic structures</summary>
PROGRAM UNITTEST_gen_multilattice_derivatives
  use derivative_structure_generator
  use io_utils, only: read_input
  use fortpy

  integer :: nmax
  logical :: full
  integer :: LatDim
  real(dp), pointer :: dFull(:,:) => null()
  integer, pointer :: cRange(:,:)
  integer :: ndfull
  integer, pointer :: equivalencies(:)
  real(dp) :: parlv(3,3)
  character(1) :: plattyp
  character(80) :: title
  integer :: k
  real(dp) :: eps
  logical :: conc_check
  character(80) :: fname = 'struct_enum.in'
  integer, pointer :: labelFull(:,:)
  integer, pointer :: digitFull(:)
  integer :: nmin


  call read_input(title, LatDim, parLV, nDFull, dFull, k, equivalencies, nMin, nMax, eps, &
                   full, labelFull, digitFull, fname, cRange, conc_check)
  if (LatDim==3) then
    platTyp = 'b'
  elseif (LatDim==2) then
    platTyp = 's'
  end if
  call gen_multilattice_derivatives(title, parLV, nDFull, dFull, k, nMin, nMax, pLatTyp, &
                                     eps, full, labelFull, digitFull, equivalencies, &
                                     conc_check, cRange)


END PROGRAM UNITTEST_gen_multilattice_derivatives
