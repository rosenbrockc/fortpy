!!<summary>Auto-generated unit test for derivative_structure_generator.gen_multilattice_derivatives
!!using FORTPY. Generated on 2014-08-13 15:52:33.711581.
!!Original enumerations tests of basic structures</summary>
PROGRAM UNITTEST_gen_multilattice_derivatives
  use io_utils, only: read_input
  use derivative_structure_generator
  use fortpy

  character(1) :: plattyp
  integer :: LatDim
  integer, pointer :: labelFull(:,:)
  real(dp) :: eps
  integer :: nmax
  character(80) :: title
  character(80) :: fname = 'struct_enum.in'
  real(dp), pointer :: dFull(:,:) => null()
  integer :: ndfull
  real(dp) :: parlv(3,3)
  integer, pointer :: cRange(:,:)
  integer, pointer :: equivalencies(:)
  logical :: conc_check
  integer :: nmin
  logical :: full
  integer :: k
  integer, pointer :: digitFull(:)


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