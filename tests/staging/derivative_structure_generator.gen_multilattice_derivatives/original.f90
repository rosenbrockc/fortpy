!!<summary>Auto-generated unit test for derivative_structure_generator.gen_multilattice_derivatives
!!using FORTPY. Generated on 2014-08-07 16:28:23.158460.
!!Original enumerations tests of basic structures</summary>
PROGRAM UNITTEST_gen_multilattice_derivatives
  use io_utils, only: read_input
  use derivative_structure_generator
  use fortpy

  real(dp) :: eps
  logical :: conc_check
  integer :: k
  integer, pointer :: labelFull(:,:)
  integer :: LatDim
  real(dp) :: parlv(3,3)
  character(1) :: plattyp
  character(80) :: fname = 'struct_enum.in'
  integer :: nmax
  integer, pointer :: equivalencies(:)
  integer, pointer :: digitFull(:)
  integer :: nmin
  integer :: ndfull
  logical :: full
  character(80) :: title
  real(dp), pointer :: dFull(:,:) => null()
  integer, pointer :: cRange(:,:)


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