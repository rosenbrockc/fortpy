!*****************************************************************************  
! this module calculates the spacegroup of a crystal, the
! pointgroup of a bravais lattice, and 
! checks whether or not a given unit cell is primitive
!
! Gus Hart
! UC Davis
! Started 10/98
! 
!     Finite precision error fixed 8/99
!
! Started using the module again for the UNCLE code
! GLWH BYU Dec. 2006
!
! Added error flag passback on get_transformations
!
! Modified SG routine to include fractional shifts for non-primitive lattices
! (10/26/2009)
MODULE symmetry_module
use num_types
use numerical_utilities, only: equal
use vector_matrix_utilities, only: matrix_inverse, norm, cross_product, volume
implicit none
private
public get_spaceGroup, get_spaceGroup_atomTypes, rm_3d_operations, make_primitive, get_lattice_pointGroup,&
      check_spaceGroup, does_mapping_exist, get_transformations,&
      bring_into_cell, find_site_equivalencies

CONTAINS
!*******************************************************************************
! This routine takes a crystal structure (basis vectors and basis atoms &
! positions) and returns the point operators and fractional translations of the
! space group. The routine assumes that the given crystal structure is already
! primitive. To reduce a  non-primitive structure to a primitive one, use the
! function "make_primitive" in this same module.
!
! No assumptions are made about the orientation of the input vectors. The
! positions of the basis atoms may be given in lattice coordinates or in
! cartesian coordinates. 
!
! The main steps are:
!
! (1) Check inputs and generate matrices for converting vectors (atom positions
! an lattice points) from (to) lattice coordinates to (from) Cartesian
! coordinates 

! (2) Convert atom positions from lattice coordinates, if necessary
! (3) Translate all atoms into primitive unit cell
! (4) Find the point operators of the given lattice
!     i) Generate all triplets of lattice points that preserve the length of
!        the original lattice vectors.
!    ii) Eliminate triplets that aren't primitive (cell volume is changed)
!   iii) Compute the transformation that takes the original lattice vectors to
!        the new triplet of points.
!    iv) Check that the transformation is orthogonal. If so, it is part of the
!        point group of the lattice.          
! (5) Find which of the point operators are part of the space group and compute
!     the corresponding fractional translations.
subroutine get_spaceGroup(aVecs, atomType, input_pos,  sg_op, sg_fract, lattcoords, eps_)

real(dp), intent(in):: aVecs(3,3)       ! Real space primitive lattice vectors
integer, intent(inout):: atomType(:)    ! Integers representing type of each basis atom 
real(dp), intent(in), pointer:: input_pos(:,:)      ! Positions of the basis atoms
real(dp), intent(out), pointer:: sg_op(:,:,:)        ! Rotations in the space group
real(dp), intent(out), pointer:: sg_fract(:,:)       ! Translations of the space group
logical, intent(in) :: lattcoords    ! If .true., atom positions are assumed to be in lattice
                      ! coordinates. Otherwise, they are treated as cartesian
real(dp), intent(in), optional:: eps_   ! "epsilon" for checking equivalence in 
                                        ! floating point arithmetic 

real(dp) atom_pos(3,size(input_pos,2)) ! Copy of input_pos
real(dp) latt_to_cart(3,3) ! Transformation for lattice coords to cartesian coords
real(dp) cart_to_latt(3,3) ! Transformation for cartesian coords to lattice coords
real(dp), pointer:: lattpg_op(:,:,:) ! The lattice point group operations
real(dp) temp_sgops (3,3,48*size(atomType)) ! Temporary matrix to store space group operations 
real(dp) temp_sgfracts(3,48*size(atomType)) ! Temporary matrix to store fractional translations
                            ! NB: the previous two allocations originally read 48*4, where *4 was a rough estimate by GH.
                            !     We probably should have something like 48*nD (nD=number of d-vectors in the parent), but right
                            !     now I (tk) don't see nD passed into this routine. Because I don't want to change all interfaces,
                            !     I resort to size(atomType) which yields the number of "atoms" or "lattice sites" within a structure.
                            !     This number is probably far too high, but certainly not too low.
real(dp) fract(3)           ! A fractional translation to check
real(dp) v(3)               ! A translated vector to check mapping for
real(dp) v2(3)              ! A second translated vector check
real(dp) eps                ! "epsilon" for checking equivalence in floating point arithmetic 

integer i                   ! Generic loop counter
integer iop                 ! Loop counter for loops over point group operations
integer jAtom, kAtom        ! Loop counters for loops over atoms 
integer this_type           ! Type of current atom being checked
integer sgop_count          ! Counter for number of space group operations
integer nAtoms              ! Number of atoms in the basis

logical err                 ! Used to check for coplanar vectors in the basis
logical mapped              ! True when an operation has mapped an atom to another


if(.not. present(eps_)) then
   eps = 1e-10_dp
else
   eps =  eps_
endif

! Get number of atoms in the basis
nAtoms = size(atomType)

! Save original original atomic input positions
atom_pos = input_pos

! A vector can be represented as either a set of cartesian coordi-
! nates or as a linear combination of primitive lattice vectors
! Get transformation matrices to take us back and forth
!!err=.false.
call get_transformations(aVecs, latt_to_cart, cart_to_latt,err)
!!write(*,*) err
if(err) stop "ERROR in get_spaceGroup: primitive vectors appear to be coplanar"  

! Convert the position of the basis atoms from lattice coordinates, if necessary
if(lattcoords .eqv. .true.) then
   do i = 1, nAtoms
      atom_pos(:,i) = matmul(latt_to_cart, atom_pos(:,i))
   enddo
endif

! Bring all the basis atoms into the unit cell
do i = 1, nAtoms
   call bring_into_cell(atom_pos(:,i), cart_to_latt, latt_to_cart, eps)
enddo

! Now, find the point group for the given lattice
call get_lattice_pointGroup(aVecs, lattpg_op, eps)

! **** Find the elements of the space group ****
!
! Count the elements
sgop_count = 0
! Apply each of the point operators in the point group to the crystal
do iop = 1, size(lattpg_op,3)
! rotate atom 1 and store its position in the vector v
   v = matmul(lattpg_op(:,:,iop), atom_pos(:,1))
! Loop over all possible fractional translations
   do jAtom = 1, nAtoms
      if(atomType(jAtom) /= atomType(1)) cycle
      fract = atom_pos(:,jAtom) - v
      call bring_into_cell(fract, cart_to_latt, latt_to_cart, eps)
! Is each atom of every type mapped by this rotation + translation?
      do kAtom = 1, nAtoms
         this_type = atomType(kAtom)
! Rotate and translate each atom        
         v2 = matmul(lattpg_op(:,:,iop),atom_pos(:,kAtom))
         v2 = v2 + fract
         call bring_into_cell(v2, cart_to_latt, latt_to_cart, eps)
! Try to map this rotated atom onto another the same type
         call does_mapping_exist(v2, this_type, atom_pos, atomType, mapped, eps)
         if(.not. mapped) exit ! no mapping for this atom
      enddo
! if all atoms mapped successfully then count this element (rotation + translation)
      if(mapped) then
         sgop_count = sgop_count + 1 ! Count the number of elements
         temp_sgfracts(:,sgop_count) = fract ! Save the translational part
         temp_sgops(:,:,sgop_count) = lattpg_op(:,:,iop) ! Store the rotational part
         !exit !exit loop over fractional translations and try next op
         ! By removing the preceding exit, we include fractional translations
         ! for non-primitive lattices. (GLWH 10/26/2009) 
      endif
   enddo
enddo

! Now that we know how many space group operations there are, store them in a matrix of
! the appropriate size as well as the fractional translations
allocate(sg_op(3,3,sgop_count),sg_fract(3,sgop_count))
sg_op = temp_sgops(:,:,1:sgop_count)
sg_fract = temp_sgfracts(:,1:sgop_count)

!LN do iOp = 1, sgop_count
!LN    do i = 1, 3
!LN       write(*,'("Op#: ",i2,3x,3(f8.4,1x))') iOp,sg_op(i,:,iOp)
!LN    end do
!LN    write(*,'("Shift#: ",3(f8.4,1x),/)') sg_fract(:,iOp)
!LN end do
!LNprint *, "-------------------------------"


end subroutine get_spaceGroup

!****************************************************************************************
! Determine the atom types "aTyp" used for get_spaceGroup.
! --> IN:
!             label  : the possible labels for a site (e.g. 0/1 or 1/2/3)
!             digit  : the number of possible labels (e.g.  2,     3)
! <-- OUT:
!              aTyp  : the atom "types"
subroutine get_spaceGroup_atomTypes(label,digit,ce,aTyp)
integer, intent(in)  :: label(:,:)  ! { digit#, dvector# }, possible labels (occupations) for each dvector 
integer, intent(in)  :: digit(:) ! { dvector# }, number of labels for each dvector
integer, intent(in)  :: ce(:)    ! { dvector# }, to which CE does a dvector belong to? (refers to CE couplings)
integer, pointer     :: aTyp(:)  ! { dvector# } intent(out)                                                                                    
integer iD1, iD2, nD, maxTyp
integer digit1, digit2
integer, pointer :: uqLabel(:,:) ! { digit#, dvector# }: unique labels for each d-vector (see comments below)
integer :: maxLabel

! "label" holds the possible occupations of the dvectors. If there are coupled lattice systems
!      Lattice = Lattice1 (x) Lattice2
! the dvectors of Lattice1 and Lattice2 are to be treated separately.
! "label", however, is not unique since the occupations on both Lattice1 and Lattice2 are
! denoted by, say, 0/1. On Lattice1, "0" may represent A, "1" B, while on lattice2, "0" may 
! represent C, and "1" D.
! --> we have to make "label" unique.
nD = size(digit)
maxLabel=maxval(label) 
allocate(uqLabel(size(label,1),size(label,2)))
do iD1=1,nD
  do digit1=1,digit(iD1)
    uqLabel(digit1,iD1) = label(digit1,iD1)  +  (ce(iD1)-1) * maxLabel  ! make label unique
  end do
end do


nD=size(digit)
if (associated(aTyp)) then; aTyp => null(); endif  ! Need strange construction here to make
if (.not. associated(aTyp)) allocate(aTyp(nD))     ! sure that an associated pointer is handled correctly

! now assign each dvector a type called "aTyp". All dvectors that have the same possible occupations are
! assigned the same aTyp. 

iD1_loop: do iD1=1,nD
  maxTyp=0
  do iD2=1,iD1-1
    digit1=digit(iD1) ! number of occupations on dvector 1
    digit2=digit(iD2) ! number of occupations on dvector 2
    if (digit1==digit2 &    ! if dvector 1 and dvector 2 have the same number of possible occupations...
        .and. all(uqLabel(:digit1,iD1)==uqLabel(:digit1,iD2))) then  ! ... and if they have the same possible occupations ...
      aTyp(iD1)=aTyp(iD2)                                            ! ... then the aTyp is identical.
      cycle iD1_loop ! found a type for iD1, go on with the next                                                                  
    else
      maxTyp = max(maxTyp,aTyp(iD2)) ! generate a new type for iD1
    endif
  enddo
  ! here: iD2=iD1-1                                                                                                               
  maxTyp = maxTyp+1
  aTyp(iD1) = maxTyp
enddo iD1_loop
end subroutine get_spaceGroup_atomTypes

!****************************************************************************************
subroutine rm_3d_operations(aVecs,sgrots,sgshifts,eps)
real(dp), intent(in):: aVecs(3,3)  
real(dp), pointer:: sgrots(:,:,:), sgshifts(:,:)  ! intent(inout)
real(dp), intent(in) :: eps
real(dp), allocatable :: tSGrots(:,:,:),tSGshifts(:,:)
integer :: nRot, iRot, status, i

if (.not.equal(aVecs(2:3,1),0._dp,eps) .or. &
    .not.equal(aVecs(1,2:3),0._dp,eps) )    &
   stop "Error in rm_3d_operations: only allowed for primitive vectors x00,0xx,0xx"

nRot = size(sgrots,3)
allocate(tSGrots(3,3,nRot),tSGshifts(3,nRot),STAT=status)
if(status/=0) stop "Failed to allocate memory in rm_3d_opertions: tSGrots"

irot = 0
do i = 1, nRot
  if (equal(sgrots(2:3,1,i),0._dp,eps) .and. equal(sgrots(1,2:3,i),0._dp,eps) .and.&
      & equal(abs(sgrots(1,1,i)),1._dp,eps)) then ! this operation is "2D"                                                        
         irot = irot + 1
         tSGrots(:,:,irot) = sgrots(:,:,i)
         tSGshifts(:,irot) = sgshifts(:,i)
  endif
enddo
nRot = irot
deallocate(sgrots,sgshifts)
allocate(sgrots(3,3,nRot),sgshifts(3,nRot),STAT=status)
if(status/=0) stop "Allocation of sgrots failed in rm_3d_operations"

sgrots=tSGrots(:,:,1:nRot)
sgshifts=tSGshifts(:,1:nRot)
end subroutine rm_3d_operations

!****************************************************************************************
! If the given lattice vectors and basis do not form a primitive unit cell, reduce the
! vectors to a set of primitive vectors and reduce the number of basis atoms, if
! necessary. Also, all the atom positions are moved inside the new unit cell
subroutine make_primitive(aVecs, atomType, atom_pos, lattCoords, eps_, removed_)
real(dp), intent(inout):: aVecs(3,3)   ! Primitive real space lattice vectors
integer, intent(inout), pointer:: atomType(:)         ! intent(inout): Atom types represented as integers
real(dp), intent(inout), pointer:: atom_pos(:,:)      ! intent(inout): Positions of the basis atoms
logical, intent(in):: lattCoords       ! .true. if positions are in lattice coordinates
real(dp), intent(in), optional   :: eps_        ! eps
integer, intent(out), pointer, optional       :: removed_(:) ! the indices of the atoms that have been removed
                                                ! if the cell is not primitive. 
integer nFracts         ! Number of possible fractional translations for this structure
integer iAtom, jAtom    ! Loop counters for looping over atoms
integer nAtoms          ! Number of atoms in the input basis
integer this_type       ! Type of atom currently being checked for mapping
integer i, j, k         ! Loop variables for looping over potentially new lattice vectors
integer iFract          ! Loop variable for loop over possible fractional translations
integer removed(size(atomType))     ! Temporary storage for removed_
integer tempType(size(atomType))    ! Temporary storage for use during reallocation
real(dp) temp_pos(3,size(atomType)) ! Temporary storage for use during reallocation
real(dp) fract(3)       ! Fractional translation
real(dp) v(3)           ! Position of current atom after a fractional translation
real(dp) cart_to_latt(3,3)  ! Transforms vectors from cartesian coords to lattice coords
real(dp) latt_to_cart(3,3)  ! Transforms vectors from lattice coords to cartesian coords
real(dp) eps            ! "epsilon" for checking equivalence in floating point arithmetic
logical mapped          ! Set to true by "does_mapping_exist" if the atom is mapped
logical err             ! Used in matrix_inverse to check for coplanar vectors

real(dp), allocatable:: fracts(:,:) ! Array of possible fractional translations
real(dp), allocatable:: lattice_point(:,:) ! Stores extra lattice points


nAtoms = size(atomType)

! optional arguments
if(.not. present(eps_)) then; eps = 1e-10_dp
                        else; eps = eps_
endif


call get_transformations(aVecs, latt_to_cart, cart_to_latt)

! Convert the position of the basis atoms from lattice coordinates, if necessary
if(lattcoords .eqv. .true.) then
   do i = 1, nAtoms
      atom_pos(:,i) = matmul(latt_to_cart, atom_pos(:,i))
   enddo
endif

! Bring all the basis atoms into the unit cell
do i = 1, nAtoms
   call bring_into_cell(atom_pos(:,i), cart_to_latt, latt_to_cart, eps)
enddo

! Number of possible fractional translations can be no larger than the number of
! translations that exist between atom 1 and all other atoms of the same type
nFracts = count(atomType == atomType(1)) - 1
!write (*,'("nFracts:",i3)') nFracts
allocate(fracts(3,nFracts),lattice_point(3,3 + nFracts))

! Count the number of fractional translations that bring the crystal 
! back onto itself. In other words, find any lattice points that are
! inside the cell
nFracts = 0 !counter for the number of fractional translations found

! If there is a lattice vector inside the unit cell, it will bring 
! any atom of any type in the basis back onto itself. Such a fract-
! ional translation must exist for EACH type of atom, so it is suf-
! ficient to check that all of the possible translations for a single
! type (might as well be first type) map every atom of each type onto
! AT LEAST one other atom of the same type
do iAtom = 2, nAtoms
   if(atomType(iAtom) /= atomType(1)) cycle
! Get fractional translation and bring it into cell (if necessary)
   fract = atom_pos(:,iAtom) - atom_pos(:,1)
   call bring_into_cell(fract, cart_to_latt, latt_to_cart, eps)
! Try this fractional translation on ALL atoms 
   do jAtom = 1, nAtoms
      this_type = atomType(jAtom)
! v contains the coordinates of the atom after the translation
      v = fract + atom_pos(:,jAtom)
      call bring_into_cell(v, cart_to_latt, latt_to_cart, eps)
! Check to see that this translation takes each atom to another
      call does_mapping_exist(v, this_type, atom_pos, atomType, mapped, eps)
! If any atom is not mapped successfully, then this fractional trans-
! lation is not a lattice point. Try the next possible translation
      if(.not. mapped) exit
   enddo
! If this loop ends successfully (mapped = .true.) then this translation
! takes all atoms to another of the same type. (It's a lattice point.)
! Count it and save it in an array
   if(mapped) then
      nFracts = nFracts + 1
      fracts(:,nFracts) = fract
   endif
enddo
!write (*,'("nFracts:",i3)') nFracts

! If the lattice is not primitive, then extra lattice points were
! found inside the given unit cell.
if (nFracts > 0) then
!   write(*,*) "Given unit cell is not primitive"
! collect all lattice points (i.e. potential new primitive vectors)
   lattice_point(:,1:nFracts) = fracts(:,1:nFracts)
   lattice_point(:,nFracts + 1:nFracts + 3) = aVecs
! total number of lattice points (including the original primitive vectors)
   nFracts = 3 + nFracts
! Take all the possible triplets of points and check to see if one
! of the triplets forms a set of primitive basis vectors. A triplet
! of points that constitutes a new set of primitive basis vectors
! will have the property that the coefficients of ALL the lattice 
! points will be integer values when the coefficients are given
! in (the new) lattice coordintates
l1:do i = 1, nFracts - 2
      do j = i + 1, nFracts - 1
         do k = j + 1, nFracts
! form the trial lattice vectors
            aVecs(:,1) = lattice_point(:,i)
            aVecs(:,2) = lattice_point(:,j)
            aVecs(:,3) = lattice_point(:,k)
! If the points are coplanar, inverse will be singular. If so, try the next triplet
            call matrix_inverse(aVecs, cart_to_latt, err, eps)
            if(err) cycle        
            do iFract = 1, nFracts
! put points in lattice coordinates using the new basis vectors
               v = matmul(cart_to_latt, lattice_point(:,iFract))
! Are all components integers? If so, new vectors found so exit all
! loops. If not try the next triplet of points
               mapped = .true.
               if(.not. equal(v, anint((v),dp), eps)) &
               then; mapped = .false.; exit; endif
            enddo
            if(mapped) exit l1 ! found new vectors so exit
         enddo
      enddo
   enddo l1
   if(.not. mapped) stop "Error in get_symmetry: New basis vectors not found"
! Bring all the atoms inside the new unit cell
   do iAtom = 1, nAtoms
      call bring_into_cell(atom_pos(:,iAtom), cart_to_latt, aVecs, eps)
   enddo
! Check to see if there were redundant atoms in the basis and remove them if so
   nAtoms = 0
   do iAtom = 1, size(atomType)
      v = atom_pos(:,iAtom)
      this_type = atomType(iAtom)
      mapped = .false.
      do j = iAtom + 1, size(atomType)
         if(atomType(j) == this_type .and. equal(v, atom_pos(:,j), eps)) then
           mapped = .true.
           removed(iAtom)=iAtom  ! store which atom we remove
         endif
      enddo   
      if(.not. mapped) then
         nAtoms = nAtoms + 1
         tempType(nAtoms) = this_type
         temp_pos(:,nAtoms) = v
      endif
   enddo
   deallocate(atomType, atom_pos)
   allocate(atomType(nAtoms), atom_pos(3,nAtoms))
   atomType = tempType(:nAtoms)
   atom_pos = temp_pos(:,:nAtoms)  
   
   if (present(removed_)) then
     allocate(removed_( count(removed/=0) ))
     removed_ = pack(removed,removed/=0)
   endif
endif

! Convert the positions of the basis atoms back to cartesion coordinates, if necessary
if(lattcoords .eqv. .true.) then
   do i = 1, nAtoms
      atom_pos(:,i) = matmul(cart_to_latt, atom_pos(:,i))
   enddo
endif
deallocate(fracts,lattice_point) 
end subroutine make_primitive

!****************************************************************************************
! this routine returns only the point group of the lattice rather than the space group
! of the given crystal structure
subroutine get_lattice_pointGroup(aVecs, lattpg_op, eps_)
real(dp), intent(in):: aVecs(3,3)
real(dp), pointer:: lattpg_op(:,:,:)
real(dp), intent(in), optional:: eps_   ! "epsilon" for checking equivalence in 
                                        ! floating point arithmetic

real(dp) temp_op(3,3,48)         ! Temporary storage for the point group operations
real(dp) new_vectors(3,3)        ! Possibly rotated set of primitive vectors
real(dp) inverse_aVecs(3,3)      ! Inverse aVecs matrix (used to check rotations)
real(dp) this_vector(3)          ! Vector used in testing
real(dp) test_matrix(3,3), rotation_matrix(3,3)
real(dp) cell_volume             ! Volume of the given unit cell
real(dp) max_norm                ! Maximum norm of the given lattice vectors
real(dp), allocatable:: Rvecs(:,:), Rlengths(:)
real(dp) norm_avecs(3)           ! Norms of the given lattice vectors
real(dp) length                  ! Length of currenct vector being checked
real(dp), parameter:: Identity(3,3) = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
real(dp) eps                     ! "epsilon" for checking equivalence in floating point arithmetic

integer n1, n2, n3              ! Upper limit for loops over R vectors
integer i, j, k      ! Loop variables
integer num_Rs       ! Number of R vectors in the sphere that is searched
integer num_ops      ! Used to count the number of point operations found

if(.not. present(eps_)) then; eps = 1e-10_dp; else; eps = eps_;endif
call matrix_inverse(aVecs, inverse_aVecs)
! Store the norms of the three lattice vectors
do i = 1, 3;norm_avecs(i) = norm(aVecs(:,i));enddo

   ! Decide how many lattice points to look in each direction to get all the 
   ! points in a sphere that contains all of the longest _primitive_ vectors
   cell_volume = abs(dot_product(aVecs(:,1),cross_product(aVecs(:,2),aVecs(:,3))))
   max_norm = max(norm(aVecs(:,1)),norm(aVecs(:,2)),norm(aVecs(:,3)))
   n1 = ceiling(max_norm*norm(cross_product(aVecs(:,2),aVecs(:,3))/cell_volume)+eps)
   n2 = ceiling(max_norm*norm(cross_product(aVecs(:,3),aVecs(:,1))/cell_volume)+eps)
   n3 = ceiling(max_norm*norm(cross_product(aVecs(:,1),aVecs(:,2))/cell_volume)+eps)

   ! Count the number of R vectors in the sphere before allocating arrays
   num_Rs = 0
   do i = -n1, n1
   do j = -n2, n2
   do k = -n3, n3
      this_vector = i*aVecs(:,1) + j*aVecs(:,2) + k*aVecs(:,3)
      length = norm(this_vector)
      if(length > max_norm + eps) cycle ! Vector is too long
      num_Rs = num_Rs + 1 ! Count number of R vectors in sphere
   enddo
   enddo
   enddo
   allocate(Rvecs(3, num_Rs), Rlengths(num_Rs))

   ! Store the R vectors that lie within the sphere
   num_Rs = 0
   do i = -n1, n1
   do j = -n2, n2
   do k = -n3, n3
      this_vector = i*aVecs(:,1) + j*aVecs(:,2) + k*aVecs(:,3)
      length = norm(this_vector)
      if(length > max_norm + eps) cycle ! This vector is outside sphere
      num_Rs = num_Rs + 1      
      Rvecs(:,num_Rs) = this_vector 
      Rlengths(num_Rs) = length
   enddo
   enddo
   enddo

   ! Try all R vector triplets in the sphere and see which ones are valid 
   ! rotations of the original basis vectors.
   ! 
   ! The length of all vectors must be preserved under a unitary transformation so skip any
   ! trial vectors that aren't the same length as the original. We also skip any set of 
   ! vectors that have the right lengths but do not form a parallelpiped that has the same
   ! volume as the original set. Also, note that the we skip sets of vectors that contain 
   ! the same vector more than once (i.e. the indices i, j, k must be unique).
   num_ops = 0
   do i = 1, num_Rs
      if(abs(Rlengths(i) - norm_avecs(1)) > eps) cycle
      do j = 1, num_Rs
         if(abs(Rlengths(j) - norm_avecs(2)) > eps) cycle
         if(j == i) cycle
         do k = 1, num_Rs
            if(abs(Rlengths(k) - norm_avecs(3)) > eps) cycle
         
            if(k == i .or. k == j) cycle
            if(abs(cell_volume - abs(volume(Rvecs(:,i),Rvecs(:,j),Rvecs(:,k)))) > eps) cycle

   ! Form the new set of "rotated" basis vectors
            new_vectors = reshape((/Rvecs(:,i),Rvecs(:,j),Rvecs(:,k)/),(/3,3/))
   ! If the transformation matrix that takes the original set to the new set is
   ! an orthogonal matrix then this rotation is a point symmetry of the lattice.
            rotation_matrix = matmul(new_vectors,inverse_aVecs)
   ! Check orthogonality of rotation matrix by [R][R]^T = [1]
            test_matrix = matmul(rotation_matrix,transpose(rotation_matrix))
            if(equal(test_matrix, Identity, eps)) then ! Found valid rotation
               num_ops = num_ops + 1 ! Count number of rotations
               temp_op(:,:,num_ops) = rotation_matrix
               !write(10,'(i5,3i4)') num_ops, i, j, k
            endif
         enddo
      enddo
   enddo

allocate(lattpg_op(3,3,num_ops))
lattpg_op = temp_op(:,:,1:num_ops)

deallocate(Rvecs, Rlengths)
end subroutine get_lattice_pointGroup

!****************************************************************************************
! This routine generates the matrices for converting vectors from lattice coordinates to
! cartesion coordinates and vice versa
subroutine get_transformations(aVecs, prim_to_cart, cart_to_prim,errflag)
real(dp) :: aVecs(3,3), cart_to_prim(3,3), prim_to_cart(3,3)
logical, optional :: errflag

prim_to_cart = aVecs
if (present(errflag))then
call matrix_inverse(prim_to_cart, cart_to_prim,errflag)
else;call matrix_inverse(prim_to_cart, cart_to_prim)
endif
end subroutine get_transformations

!********************************************************************
! This subroutine translates a point back into the unit cell
! In lattice coordinates, the coefficients of the point must all be
! less than one and at least zero
subroutine bring_into_cell(v, cart_to_latt, latt_to_cart, eps)
real(dp), intent(inout) ::  v(3)
real(dp), intent(in)    ::  cart_to_latt(3,3), latt_to_cart(3,3), eps
integer c, maxc


! Put the representation of the point into lattice coordinates
v = matmul(cart_to_latt, v)

! counter to catch compiler bug
c = 0
maxc = max(ceiling(abs(maxval(v))),ceiling(abs(minval(v)))) *2

! If a component >= 1, translate by subtracting a lattice vector
! If a component < 0, translate by adding a lattice vector
do while(any(v >= 1.0_dp - eps) .or. any(v < 0.0_dp - eps)) 
   c = c +1
   v = merge(v, v - 1.0_dp, v <  1.0_dp - eps) 
   v = merge(v, v + 1.0_dp, v >= 0.0_dp - eps)
   if (c>maxc) stop "ERROR: loop does not end in bring_into_cell. Probably compiler bug."
enddo

! Put the point back into cartesion coordinate representation
v = matmul(latt_to_cart, v)
end subroutine bring_into_cell

!******************************************************************************
! checks to see if a mapping exists between the vector v and
! the position of any of the atoms of type "this_type".
! If a mapping exists, then the logical "mapped" is 
! returned .true., otherwise .false.
subroutine does_mapping_exist(v, this_type, atom_pos, atomType, mapped, eps)

real(dp), intent(in):: v(3)            ! Position to check mapping for
integer, intent(in):: this_type        ! Type of atom that is being checked
real(dp), intent(in):: atom_pos(:,:)   ! Positions of basis atoms
integer, intent(in):: atomType(:)      ! Types of atoms in the basis
logical, intent(out):: mapped          ! Set to .true. if a mapping exists
real(dp), intent(in):: eps             ! Epsilon for checking equivalence

integer i   ! Loop over atoms
real(dp) this_position(3) ! Position of atom to be checked
       
mapped = .false.
do i = 1, size(atomType)
   if(atomType(i) == this_type) then
      ! if the coordinates are the same, 
      ! their difference will be zero for every component
      this_position = atom_pos(:,i)
      if(equal(v, this_position, eps)) &
      then; mapped = .true.; exit; endif
   endif
enddo
end subroutine does_mapping_exist

! **********************************************************
subroutine check_spaceGroup(SGop, SGfrac, eps_)
real(dp), pointer :: SGop(:,:,:)
real(dp), pointer :: SGfrac(:,:)
real(dp), optional, intent(in):: eps_

integer i, j, k, Nops
real(dp) eps
real(dp) :: testop(3,3)
logical exists

if(.not. present(eps_)) then; eps = 1e-10_dp; else; eps = eps_;endif

open(11, file="sym_check.out", status="unknown")
! Are the operations unique? (Necessary but *insufficient* condition)
Nops = size(SGop,3)
do i = 1,Nops
   do j = i+1, Nops
      if (equal(SGop(:,:,i),SGop(:,:,j),eps)) then
         write(11,*) "Error: SG Operations that were found are not unique"
         stop
      endif
   enddo
enddo
write(11,*) "The SG operations are unique"

! Do the operations form a group? (Again, this is a necessary but not
! sufficient condition to double check)
write(11,'(3x)',advance="no")
do i = 1, Nops; write(11,'(i3)',advance="no") i; enddo; write(11,*)
do i = 1, Nops
   write(11,'(i3)',advance="no") i   
   do k = 1, i-1; write(11,'(3x)',advance='no');enddo
   do j = i, Nops
      exists = .false. 
      ! Is the product of SG_i x SG_j in the set?
      testop = matmul(SGop(:,:,i),SGop(:,:,j))
      do k = 1, Nops
         if (equal(testop,SGop(:,:,k),eps)) then
            exists = .true.
            write(11,'(i3)',advance="no") k
            exit
         endif
      enddo
      if (.not. exists) then
         write(11,*) "The set of SG ops doesn't form a group"
         stop
      endif
   enddo
   write(11,*)
enddo
write(11,*) "The SG ops form a group"
close(11)

end subroutine check_spaceGroup

!***************************************************************************************************
SUBROUTINE find_site_equivalencies(pBas, siteLabel, pLV, BasEq, nSites, eps_)
real(dp), pointer :: pBas(:,:)
real(dp), intent(in) :: pLV(3,3)
integer, pointer :: BasEq(:), siteLabel(:)
real(dp), intent(in), optional :: eps_
integer, intent(out) :: nSites

real(dp), pointer :: tBas(:,:)
real(dp) :: invLV(3,3), v(3)
real(dp), pointer :: sg_rot(:,:,:), sg_fract(:,:)
integer :: iOps, iBas, jBas
integer :: nOps, nBas
logical :: equivalent
real(dp) :: eps

if(.not. present(eps_)) then
   eps = 1e-10_dp
else
   eps =  eps_
endif

allocate(tBas(size(pBas,1), size(pBas,2)))
tBas = pBas

call matrix_inverse(pLV, invLV)
nBas = size(tBas,2)
allocate(BasEq(nBas))
BasEq = 0
BasEq(1) = 1  !independent of basis size, the first always corresponds to itself
nSites = 1
call get_spaceGroup(pLV, sitelabel, tBas,  sg_rot, sg_fract, .false.,eps)
nOps = size(sg_rot,3)

do iBas=2, nBas
   call bring_into_cell(tBas(:,iBas),invLV,pLV,eps)
   equivalent = .false.
   do jBas=1, iBas-1
      do iOps = 1,nOps
         v = matmul(sg_rot(:,:,iOps),tBas(:,jBas))+sg_fract(:,iOps)
         call bring_into_cell(v,invLV,pLV,eps)
         if (equal(v,tBas(:,iBas),eps)) then
            equivalent = .true.
            exit
         endif
      end do
      if (equivalent) then
         BasEq(iBas) = BasEq(jBas)
         exit
      endif
   end do
   if (.not.equivalent) then
      nSites = nSites + 1
      BasEq(iBas) = nSites
   endif
end do
if (any(BasEq==0)) stop "There was a bug in find_site_equivalencies"

ENDSUBROUTINE find_site_equivalencies

END MODULE symmetry_module
