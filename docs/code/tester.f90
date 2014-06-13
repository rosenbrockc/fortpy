!!<summary>Auto-generated unit test for derivative_structure_generator.gen_multilattice_derivatives
!!using FORTPY. Generated on 2013-09-12 17:32:25.818542.</summary>
PROGRAM UNITTEST_gen_multilattice_derivatives
  use derivative_structure_generator
  use io_utils, only: read_input
  use fortpy

  integer :: nmax
  logical :: full
  character(1) :: plattyp
  character(80) :: title
  integer :: LatDim
  integer :: k
  real(dp) :: eps
  real(dp), pointer :: dFull(:,:) => null()
  integer, pointer :: equivalencies(:)
  logical :: conc_check
  integer, pointer :: cRange(:,:)
  integer :: nmin
  character(80) :: fname = 'struct_enum.in'
  integer :: ndfull
  integer, pointer :: labelFull(:,:)
  integer, pointer :: digitFull(:)
  real(dp) :: parlv(3,3)

  call read_input(title, LatDim, parLV, nDFull, dFull, k, equivalencies, nMin, nMax, eps, &
                  full, labelFull, digitFull, fname, cRange, conc_check)
  if (LatDim==3) then
    platTyp = 'b'
  elseif (LatDim==2) then
    platTyp = 's'
  endif
  call gen_multilattice_derivatives(title, parLV, nDFull, dFull, k, nMin, nMax, pLatTyp, &
                                    eps, full, labelFull, digitFull, equivalencies, &
                                    conc_check, cRange)
  
END PROGRAM UNITTEST_gen_multilattice_derivatives
