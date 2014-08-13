MODULE setup_cluster_expansion
use ce_types
use num_types
use symmetry_module
use vector_matrix_utilities
use figure_enumeration
use numerical_utilities
use io_utilities
use combinatorics
use structure_utilities

implicit none

private
public generate_eqvClusters,  check_input_structures, setup_all_correlations, setup_reference_energy, setup_damping

CONTAINS
!***************************************************************************************************
! - Just get the correlations for structures in the list
! - Also calculate correlations for reference structures and get their energy
recursive SUBROUTINE setup_all_correlations(inpStr,clusters,CE,isRef)

type(CE_variables), intent(inout) :: CE    ! intent(INOUT)
type(crystal),      pointer       :: inpStr(:)  ! intent(IN)
type(figRep),       pointer       :: clusters(:)
logical, optional                 :: isRef

type(crystal), pointer :: currstr, refstr(:)
logical isReference
integer iRef, nRef, c

integer, pointer :: label(:)
integer :: iStr, iFg, iRep, L(3,3), D(3,3), iBas
real(dp) invA(3,3)
integer :: nStr, nCl, status, iS, nS
integer(1), pointer :: gTab(:,:,:,:)

if (present(isRef)) then; isReference=isRef; else; isReference=.false.; endif

! Empty the missing_analysis.out file
open(32,file='missing_analysis.out',status='unknown')
close(32,status='delete')

!At first we want to check, how many inequivalent sites we have within our unit cell
!in order to see how many inequivalent onsitecorrelations there are.
!call get_site_labels(label,CE), replaced by
call get_spaceGroup_atomTypes(CE%label, CE%digit,CE%CEBas, label)
call find_site_equivalencies(CE%pBas, label, CE%pLV, CE%BasEq, CE%nBasIneq, CE%eps)
           ! --> celib/symmetry_module.f90
write(*,'("Site equivalencies: ",100i3)') CE%BasEq

!Now we can proceed to the actual calculation of correlations
nStr = size(inpStr)
nCl = size(clusters)
write (*,'(A, I4)') "Number of input structures in setup_all_correlations: ", nStr
 print *, "calculating correlations for structure:"
!________________________________________________________________________________
do iStr = 1, nStr
  if (inpStr(iStr)%isHelpStructure) cycle
  write (*, '(I4,A,I4,5x,A30)') iStr,"/",nStr, adjustl(inpStr(iStr)%title)
  call create_spin_lookup_table(inpStr(iStr),gTab,invA,L,D, determineLD=.true.) ! -> structure_utilities
  call get_gTab_correlations(CE, gTab, L, D, invA, clusters, inpStr(iStr)%PI, mustGetG=.true.)
  deallocate(gTab)
  do iFg = 0, nCl-1
    do iRep = 1, size(clusters(iFg)%Rep)
      deallocate(clusters(iFg)%Rep(iRep)%gPt)
    enddo
  enddo
enddo
!________________________________________________________________________________

!!--------------------------------------------------------------------------------
!! PI output:
open(33,file="inpStr_correlations.out",iostat=status)
if(status/=0) stop "Couldn't open file: inpStr_correlations.out"
open(34,file="PI_matrix.out",iostat=status)
if(status/=0) stop "Couldn't open file: PI_matrix.out"
do iStr = 1, nStr
   write(33,'(A80)') inpStr(iStr)%title
   do iFg = 0, nCl-1
      write(33,'("Clusters #: ",i3,1x,"Corr: ",f7.3)') iFg, &
           inpStr(iStr)%PI(iFg)
      write(34,'(f9.6,1x)',advance="no") inpStr(iStr)%PI(iFg)
   enddo
   write(34,*)
enddo
close(33)
close(34)

!--------------------------------------------------------------------------------
! References:
if (.not. isReference) then
  nRef=CE%nreferences
  do iRef=1,nRef
    !if (associated(refStr)) deallocate(refStr)
    allocate(refStr(nStr))
    do iStr=1,nStr
      refStr(iStr) = inpStr(iStr)%reference(iRef)
    enddo
    call setup_all_correlations(refStr,CE%reference(iRef)%eqvClusters,CE%reference(iRef),isRef=.true.)
    do iStr=1,nStr
      if (inpStr(iStr)%isHelpStructure) cycle
!      refstr(iStr)%ce_energy = calculate_energy(CE%reference(iRef), refStr(iStr), CE%reference(iRef)%J)
      refstr(iStr)%ce_energy = calculate_energy( CE%reference(iRef)%J , refStr(iStr)%PI )
      ! this is the reference energy PER ATOM of the reference structure
      ! ==> need to rescale, in order to fit to an extended input structure:      
      refstr(iStr)%ce_energy = refstr(iStr)%ce_energy * refstr(iStr)%nAtoms * 1.0_dp / inpStr(iStr)%nAtoms

      inpStr(iStr)%reference(iRef)=refStr(iStr)
    enddo
  enddo
endif


ENDSUBROUTINE setup_all_correlations



!***************************************************************************************************
! - This routine takes in a list of input structures and a parent lattice from the lat.in file
!   and makes sure the structures are consistent with the parent lattice and that there are no
!   duplicate structures in the list.
! - Determine the right concentration of the structure (e.g. for surfaces you don't count the bulk
!   like atoms in the structure)
! - Build reference structure and check it too
recursive SUBROUTINE check_input_structures(str, CE, thisRef)
use compare_structures, only: is_derivative
use enumeration_utilities, only: get_HNF_of_derivative_structure, get_gspace_representation, compare_two_gstructures
type(crystal), pointer :: str(:), currstr, refstr(:)
type(crystal)          :: mirrorstr
type(CE_variables), intent(in) :: CE
integer, optional :: thisRef
integer :: iRef, nRef
logical :: isReference

real(dp) :: conVec(3), pVol, currD(3), currPos(3), invLV(3,3)
integer sVol
integer iStr, jStr, iAt, jAt, iuq, dlabel, c, iC, nC, iD, nD
integer nStr, nPar, nAt, x, y, z
integer, allocatable :: uqlist(:)
logical checkPoint, equivalent, found, mapped, err, performedReconstruction

integer idx1, idx2, idxoff, ik, iat1, iat2, i
integer, pointer :: gettype(:), rankSpace(:), missingD(:)

character(80) :: header(4), filename

! structure comparison:
integer :: SNF(3,3), L(3,3), LatDim
logical :: match
logical, allocatable :: mask(:)

call get_SpinRank_mapping(CE%CErank, gettype, rankSpace)


nC = CE%ncouplings

if (present(thisRef)) then; isReference=.true.; else; isReference=.false.; endif
if (.not. isReference) then
  write(*,'(A)') "Checking input structures"
endif

nStr = size(str)
if (nStr==0) stop "ERROR: check_input_structures: No input structures" 
nPar = size(CE%pBas,2)

do iStr = 1, nStr

   currStr => str(iStr)
   
   if (currStr%isHelpStructure) cycle  ! do not check mere help structures

   if (.not. isReference) then
   !LN  write(*,'(A,1x,I5)') "Checking structure #", iStr
   else
     write(*,'(A,1x,I5,1x,A,1x,I5)') "Checking Structure #", iStr, " Reference #", thisRef
   endif


   ! Is each structure a derivative superlattice of the parent lattice?
   if (.not. is_derivative(CE%pLV,currStr%LV,CE%eps)) then
     if (.not. isReference) then
       write(*,'("ERROR: Structure #: ",i4," is not a derivative superlattice.")') iStr
       write(*,'(A,A)') "ERROR:    struc title: ", trim(currStr%title)
       stop
     else
       write(*,'(A,I5,A)') "ERROR: the reference #",iRef," structure is not a derivative superlattice."
       stop
     endif
   endif

   nAt     = size(currStr%spin) ! total number of atoms
   nD      = CE%nBas  ! number of dvectors
   sVol    = nAt / nD ! volume of this structure in terms of unit cells
   

   currStr%xnTyp = 0

   allocate(currStr%pBVlab(nAt))  ; currStr%pBVlab = -1

   ! Do all the interior points of the input structure correspond to points of the parent?
   do iAt = 1, nAt-currStr%missingN

      ! get the dset label of that atom
      dlabel = get_dset_label( CE%pLV, currStr%pos(:,iAt), CE%pBas(:,:), CE%eps )
      currStr%pBVlab(iAt)=dlabel
      ! this atom of this structure does not lie on a lattice point:
      if(dlabel<0) then
        if (.not. isReference) then
          write(*,'(A,I4,A)') "ERROR: in input structure #", iStr
        else
          write(*,'(A,I4,A)') "ERROR: the reference #",iRef," structure"
        endif
        write(*,'(A,A80)')    "ERROR: name:  ",  adjustl(currStr%title)
        write(*,'(A,I4,A)')   "ERROR:   atom ", iAt, " is not on the parent lattice."
        write(*,'(A,3F5.2)')  "ERROR:   cartesian coordinates of atom: ", currStr%pos(:,iAt)
        stop
      endif

      ! 
      ! if we are here: this atom lies on a lattice point, all is ok (up to now at least...)
      !

   enddo ! loop over all atoms in this structure

   ! now reconstruct missing atoms:
   ! (note: this is a brute force method, not really elegant, but works...)
   performedReconstruction=.false.
   if (currStr%missingN>0) then
      performedReconstruction=.true.
      if (currStr%missingN>50) then
        write(*,'(A)')   "    ... reconstructing atoms: this can take a while for >50 atoms"
      else
        write(*,'(A)')   "    ... reconstructing atoms"
      endif
      
      allocate(missingD( currStr%missingN )) ! allocating an array that is going to held the missing dvectors.
                                             ! If a dvector is missed more than once, it appears more than once!
      missingD = 0
      
      idx1=1
      do iD=1,nD ! loop through the dset
        idx2 = idx1 + (sVol - count(currStr%pBVlab==iD)) -1
                      ! <----------------------------> !
                              ! this expression tells you how often a specific dvector is missing. It should be
                              ! there the same number of times as we have unit cells (sVol = structure volume).
        if (idx2>=idx1) then  ! If it is missing, we ...
          missingD(idx1:idx2) = iD   ! ... add the dvector the correct number of times to the array
          idx1 = idx2+1
        endif
      enddo
      iAt = nAt-currStr%missingN
      do c=1,currStr%missingN  ! loop through all missing atoms. The following is not really elegant, but works.

        if (currStr%missingN > 50) write(*,'(I5,1x)',advance='no') c
   
        shiftloop: do x=0,sVol; do y=0,sVol; do z=0,sVol;
          currD = CE%pBas(:, missingD(c)) ! take the dvector...
          currD = currD + x*CE%pLV(:,1) + y*CE%pLV(:,2) + z*CE%pLV(:,3) ! ... and shift it by parent lattice vectors
          currStr%pos(:,iAt+c)  = currD         ! add the current dvector coordinate to the structure
          currStr%pBVlab(iAt+c) = missingD(c)   ! and set the label correctly
          call check_for_duplicate_atoms(currStr, iat1, iat2, CE%eps, limit=iAt+c)   ! check if added basis point
                                                ! collides with an already existing point. If yes, the two integers
                                                ! iat1 and iat2 provide the ID numbers of the atoms that collide
          if (iat1==0 .and. iat2==0) then  ! if no collision: both iat1 and iat2 == 0 and ...
            exit shiftloop ! ... we can go on with the next missing atom.
          endif
        enddo; enddo; enddo shiftloop;

        ! failsafe check:
        if (x==sVol .and. y==sVol .and. z==sVol) stop "ERROR: bad programming in the reconstruction loop"
      enddo
   
   endif ! reconstructing
   if (currStr%missingN > 50) write(*,*)



   ! If this CE has some dvectors with fixed occupation we store the respective atoms
   ! of this structure in a separate variable
   allocate(currStr%fixedOcc( sVol * size(CE%fixedOcc) )) ! there will be "number of unit cell for this struc" 
                                                          ! times "fixed dvector occupation"
                                                          ! number of entries
   if (size(currStr%fixedOcc)>0) then  ! if there are some fixed occupation atoms:
     allocate(mask(nAt))               ! - create a temporary mask
     mask = .false.                    !   and initialise it
     where ( &
                (/   (any(CE%fixedOcc==currStr%pBVlab(iAt)) , iAt=1,nAt)    /) &   
           ) mask=.true.               ! - mask is set true, wherever a pBVlab coincides with a fixed occupation
                                       !   dvector
                                       ! - now pack the indices of the atoms according to mask. In fixedOcc there are
                                       !   now the indices of those atoms that have fixed occupation according to the
                                       !   current CE.
     currStr%fixedOcc = pack(  (/ (iAt,iAt=1,nAt) /), mask )  
     deallocate(mask)
   endif

   ! Determine atom types and spins for the structure
   if (.not. isReference) then ! for references this is done during the buildup of the reference structure itself
     do iAt=1,nAt
       currStr%atyp = 1              ! setup: all atoms have the same type...
       currStr%spin = gettype(1)     ! ... and the respective spin
       do iC = 1, CE%ncouplings
         do ik = 2, CE%CErank        ! now: for all other types...
           
           idxoff = sum(currStr%nTyp(:,:iC-1))         ! all the atoms that were already counted
           idx1   = sum(currStr%nTyp(1:ik-1,iC))+1  + idxoff 
           idx2   = sum(currStr%nTyp(1:ik,iC))      + idxoff           
           currStr%atyp( idx1 : idx2 ) = ik            ! ... set the type
           currStr%spin( idx1 : idx2 ) = gettype(ik)   ! ... and spin of those atoms
         enddo
       enddo
     enddo
   endif


   do iAt=1,nAt

      ! check if the occupation of every site is ok with the possible occupations in lat.in
      ! CE%label(:,iD) has the possible occupations labels (rank space) for dset member iD.
      ! Thus, the following checks if any one of the possible labels fits to the structure's occupation.
      ! Careful: %label = 0,1,2,3...,k-1 like in enumlib
      !          %aTyp  = 1,2,3,...k
      dlabel = currStr%pBVlab(iAt)
      if (.not. any(CE%label(:,dlabel)==currStr%ATyp(iAt)-1)) then
        if (.not. isReference) then
          write(*,'(A,I4,A)')   "ERROR: in input structure #", iStr
        else
          write(*,'(A,I4,A)')   "ERROR: in reference #",iRef," structure"
        endif
        write(*,'(A,A80)')      "ERROR: name:  ", currStr%title
        write(*,'(A)')          "ERROR:   the structure's occupation does not fit the occupation in lat.in."
        write(*,'(A,I3,A,I2)')  "ERROR:   structures.in: structure atom        ", iAt,  " has type  ", currStr%ATyp(iAt)-1
        write(*,'(A,I3,A,10(I2,1x))') &
                                "ERROR:   lat.in       : equivalent basis atom ", dlabel," has types ", &
                                CE%label(:CE%digit(iAt),dlabel)
        stop 
      endif
            
      ! update xnTyp, which stores the total number of atoms of a certain kind, but only
      ! counts x-relevant (i.e. concentration-relevant) structures
      if (CE%xRelevant(dlabel)) then ! if this is a concentration-relevant dset member
        currStr%xnTyp(currStr%ATyp(iAt),CE%CEBas(dlabel))=currStr%xnTyp(currStr%ATyp(iAt),CE%CEBas(dlabel)) + 1
      endif
   enddo

!OLD>   !____________________________________________________________
!OLD>   ! is this a surface CE with symmetric slab?
!OLD>   if (CE%surface%isSurf) then
!OLD>      if (CE%surface%isSymmetric) then
!OLD>
!OLD>        if (.not. structure_is_symmetric(currStr%pos, currStr%spin, CE%surface, CE%eps)) then
!OLD>          write(*,'(A,I4,A)')     "ERROR: in surface input structure #", iStr
!OLD>          write(*,'(A,A80)')      "ERROR: name:  ", currStr%title
!OLD>          write(*,'(A)')          "ERROR:     the structure is not symmetric"
!OLD>          stop
!OLD>        endif
!OLD>
!OLD>      endif ! symmetric slab
!OLD>    endif ! surface CE
!OLD>   !____________________________________________________________


   ! Store parent lattice information for each input structure
   currStr%pLV = CE%pLV 
   allocate(currStr%intPt(3,size(CE%pBas,2)))
   currStr%intPt = CE%pBas
   currStr%eps = CE%eps

   ! update concentration
   !    calculate the concentration from xnTyp, which contains the number of concentration-relevant
   !    atoms in this structure. So, %x contains as many components as there are atom types.
   do iC=1,nC
     currStr%x(:,iC) = currStr%xnTyp(:,iC) * 1.d0 / sum(currStr%xnTyp(:,iC))
   end do
   !    for compatibility reasons: %conc is the quasi-binary concentration.
   currStr%conc = currStr%x(1,1)  

enddo ! loop over all structures, iStr


if (.not. isReference .and. performedReconstruction) then
   filename  = "structures.out"
   header(1) = "# Postprocessed structures.in file"
   header(2) = "#    This file is only written if some of the original atoms were marked"
   header(3) = "#    missing by the use of -n in structures.in, where n is the number of"
   header(4) = "#    missing atoms"
   call write_InpStrucList(filename,str,header)
endif


!--------------------------------------------------------------------------------
! (in routine: check_input_structures)
! References:
! Now setup the reference structure for the current input structures, check if
! they are ok with the lattice definition (etc.)
if (.not. isReference) then
  nRef=CE%nreferences
!  if (associated(refStr)) deallocate(refStr)
  allocate(refStr(nStr))
  do iRef=1,nRef ! setup all references...
    do iStr=1,nStr ! ... for all input structures
      if (str(iStr)%reference(iRef)%automaticReference) then
        write(*,'(A,1x,I5,1x,A,1x,I5)') "Setting up automatic reference #", iRef, " for structure #", iStr
        call setup_reference_structure(str(iStr),refStr(iStr),CE,iRef)
        allocate(refStr(iStr)%nTyp(CE%CErank,CE%ncouplings))   ! don't want to do these allocs in the
        allocate(refStr(iStr)%xnTyp(CE%CErank,CE%ncouplings))  ! setup_reference_structure routine, since it is
        allocate(refStr(iStr)%x(CE%CErank,CE%ncouplings))      ! also going to be used during gss
      else
        write(*,'(A,1x,I5,1x,A,1x,I5)') "Setting up user-input reference #", iRef, " for structure #", iStr
        refStr(iStr) = str(iStr)%reference(iRef)  
      endif
    enddo ! end: setup 1 reference (iRef) for all input structures

    ! check the reference structures for reference iRef:
    call check_input_structures(refStr, CE%reference(iRef),iRef)

    do iStr=1,nStr
      str(iStr)%reference(iRef)=refStr(iStr)
    enddo
  enddo ! iRef
  deallocate(refStr)
endif



! Check that each input structure is unique from every other
if (.not. isReference) then

   ! The comparison---which, in prior versions, was primarily based on symmetry operation on the real 
   ! space lattice---is done in gspace, in the format that also the struct_enum.out file uses.
   ! For this, we...
   ! ... first need to convert the (input) structures to gspace format (part A)
   ! ... then simply compare the gspace formats of all structures      (part B)

   ! set the dimension of the lattice
   if(CE%surface%isSurf) then; LatDim = 2;
                         else; LatDim = 3;
   endif;

   open(100,file='duplicate_structures.out',status='replace')     ! write the duplicate structures to a file for the user
   write(100,'(A)') "# This file lists the duplicate structures in structures.in"
   write(100,'(A)') "#"
   write(100,'(A)') "# index of structure 1 | index of structure 2 | titles of the structures"
   write(100,'(A)') "# in structures.in     | in structures.in     |"
   write(100,'(A)') "#--------------------------------------------------------------------------------"

   write(*,*)
   write(*,'("Checking the uniqueness of the input structures")') 
   write(*,'("- ",i6," structures in the input list.")') nStr

   ! part A
   write(*,'(A)',advance='no'), "- get HNF of structure # "
   do iStr=1,nStr
!     write(*,'(I4,1x)',advance='no') iStr
     call get_HNF_of_derivative_structure(str(iStr)%title,   & ! wer bist du? luzi
                                          str(iStr)%LV,      & ! in
                                          str(iStr)%pos,     & ! inout
                                          str(iStr)%aTyp-1,  & ! in
                                          CE%pLV,            & ! in
                                          CE%pBas,           & ! in
                                          str(iStr)%HNFlist, & ! if the structure was not primitive, HNFList==0 (let's change that later...)
                                          SNF,               & ! out
                                          L,                 & ! out
                                          CE%eps,            &
                                          fixedOcc=str(iStr)%fixedOcc)
     if (.not. all(str(iStr)%HNFlist==0)) &
        call get_gspace_representation(CE%pLV,       & ! in
                                    CE%pBas,           & ! in
                                    str(iStr)%LV,      & ! in
                                    str(iStr)%pos,     &
                                    str(iStr)%aTyp-1,  & ! in
                                    str(iStr)%HNFlist, & ! in
                                    LatDim,            &
                                    str(iStr)%pLabel,  & ! out
                                    CE%eps)

  enddo
   ! part B
   write(*,'(/,A)'), "- comparing #"
   do iStr=1,nStr
   write(*,'(5x,I4)') iStr   ! for complex structures the ifort hangs here if I don't write something
                             ! out, w/o advance='no'
   do jStr=iStr+1,nStr
     match=.false.
     if (all(str(iStr)%HNFlist==0) .or. all(str(jStr)%HNFlist==0)) cycle   ! <--- we should remove that and make it better at some point!
     call compare_two_gstructures(LatDim,CE%pLV,CE%pBas,CE%eps, &
          str(iStr)%HNFlist(:,:,1), str(iStr)%pLabel(1,:), &
          str(jStr)%HNFlist(:,:,:), str(jStr)%pLabel(:,:), &
          match) 
     ! write(*,'(5x,L)') match
     if (match) then ! The two structures match each other
       write(*,*) ""
       write(*,'("WARNING: in structures.in: structures ",i4," and ",i4," appear to be the same.")') iStr, jStr
       write(*,'("WARNING:    see file duplicate_structures.out for details")')
       write(*,*) ""
       write(100,'(I5,1x,I5,5x,A1,A,A,A,A1)') iStr, jStr, """", trim(adjustl(str(iStr)%title)),""" <--> """,trim(adjustl(str(jStr)%title)),""""
       !stop
     endif
   enddo
   enddo
   
   !old comparison>do iStr = 1, nStr-1
   !old comparison>   write (*,'("Checking uniqueness of Str. #",I4,"  ",A50)') iStr, adjustl(str(iStr)%title)
   !old comparison>   do jStr = iStr+1, nStr
   !old comparison>!      write(*,'("Comparing strs. #: ",2(i3,1x))') iStr, jStr
   !old comparison>
   !old comparison>      equivalent=.false.  ! tk: actually, shouldn't need this, but somethings going wrong in 
   !old comparison>                          !     compare_arbitrary_structures, where it is set =.false. if
   !old comparison>                          !     this optional argument is present. But it fails quite often. 
   !old comparison>                          !     Don't ask ME about it...
   !old comparison>      call compare_arbitrary_structures(str(iStr)%LV,str(iStr)%spin,str(iStr)%pos,&
   !old comparison>                                        str(jStr)%LV,str(jStr)%spin,str(jStr)%pos,&
   !old comparison>                                          CE%eps,mapped,identical=equivalent)
   !old comparison>      if (equivalent) then ! The structures are equivalent *and* have exactly the same concentration
   !old comparison>         write(*,'("WARNING: *** In the input list, structures ",i4," and ",i4," appear to be the same.")') iStr, jStr
   !old comparison>         write(*,'("WARNING:    Structure title: ",a80)') adjustl(str(iStr)%title)
   !old comparison>         write(*,'("WARNING:    Structure title: ",a80)') adjustl(str(jStr)%title)
   !old comparison>         write(100,'(I5,1x,I5,5x,A1,A,A,A,A1)') iStr, jStr, """", trim(adjustl(str(iStr)%title)),""" <--> """,trim(adjustl(str(jStr)%title)),""""
   !old comparison>         !stop
   !old comparison>      endif
   !old comparison>   enddo
   !old comparison>enddo

   close(100)
endif ! is no Reference

! Check that each structure doesn't have two atoms on the same site
do iStr = 1, nStr ! Loop over each structures
  call check_for_duplicate_atoms(str(iStr),iAt,jAt,CE%eps)
  
  if (iAt/=0 .or. jAt/=0) then
    if (.not. isReference) then
      write(*,'(A,I4,A)') "ERROR: in input structure #", iStr
    else
      write(*,'(A,I4,A,I4,A)') "ERROR: in reference #",iRef," structure #",istr
    endif
    write(*,'("ERROR:   Structure has duplicate atoms on one site.")') 
    write(*,'("ERROR:   Duplicate atoms are numbers ",i3," and ",i3)') iAt, jAt
    write(*,'("ERROR:   Structure title: ",a80)') adjustl(str(iStr)%title)
    stop
  endif

enddo

! Last check is to see if each structure has the correct volume/atom. If it does, and 
! it passed the checks above then all the atomic sites in the cell must be occupied.
pVol = abs(determinant(CE%pLV)) ! Volume of the parent cell
if (nD==1) then
do iStr = 1, nStr ! Loop over each structure
   ! Check the volume/atom ratio of each structures with the volume of the parent cell
   ! This will not work as coded if the parent cell has more than one site/cell (that's why
   ! there's the if statement around the do-loop)
   if (.not.(equal(nPar*abs(determinant(str(iStr)%LV))/size(str(iStr)%spin,1),pVol,CE%eps))) then
      write(*,'("There is a problem with the volume/atom ratio input structure ",i4)') iStr
      write(*,'("Structure title: ",a80)') adjustl(str(iStr)%title)
      stop
   endif
enddo
endif

! Last check. Is the number of types in the input structure consistent with the CE rank from lat.in?
! This routine won't win any elegance awards but it seems to do the trick.
allocate(uqlist(CE%CErank))
uqlist = -99  ! set it to something very unlikely
do iStr = 1, nStr ! Loop over each structure
   nAt = size(str(iStr)%pos,2)
   do iAt = 1, nAt ! Loop over each atom in the current structure
      found = .false.
      do iuq = 1, count(uqlist/=-99) ! Loop over every known (so far) atom type
         if (iuq>CE%CErank) then
            write(*,'("More atom types in input structures than in lat.in file")')
            write(*,'("Structure title: ",a80)') adjustl(str(iStr)%title)
            stop
         endif
!         write(*,'("iuq",i3,"  rank",i3,"  str",i3,"  iAt",i3," uqlist",i3,"  spin",i3)') &
!                       iuq,CE%CErank,istr,iAt,uqlist(iuq),str(iStr)%spin(iAt)
         if (uqlist(iuq)==str(iStr)%spin(iAt)) found = .true.
      enddo ! End loop over unique types found so far
      if (.not. found) uqlist(iuq) = str(iStr)%spin(iAt)
   enddo
enddo
! I don't think there are any other possibilities of bad input if all these tests pass...
ENDSUBROUTINE check_input_structures

!*************************************************************************************************
SUBROUTINE generate_eqvClusters(clusters, eqvClusters, CE, MPI_rank,iRef)
type(figRep), pointer :: eqvClusters(:)  ! intent(out): This is the list of symmetrically-equivalent figures
type(figure), intent(in) :: clusters(0:) ! { cluster# }
type(CE_variables), intent(in) :: CE
integer, intent(in) :: MPI_rank
integer, optional :: iRef

real(dp) :: invLV(3,3), diffVec(3), lattCoordCheck(3)
integer, pointer :: sitelabel(:)
integer nCl, nV, nstore, nOps, nS, i, j   ! counter for figures, vertices, stored figs
integer iCl, iV, isym, istore, iBas, iS, jV ! loop index for figures, vertices, symops
type(figure), pointer :: tempStore(:) ! Need these to temporarily store the symmetry "brothers"
type(figure) :: rotCluster

real(dp), allocatable :: diff(:,:), rotPts(:,:)
real(dp), pointer :: rot(:,:,:), shift(:,:)
real(dp) :: tVec(3) ! temporary copy of vertex # 1 (to shift into parent cell 0)
real(dp) :: shiftVec(3) ! amount to shift each vertex to move vertex #1 int parent cell 0)
logical fixed, diffcheck, err, stored
integer iOps, dupc, ip, ic, jp, sRank
integer, allocatable :: permute(:), ps(:)
integer, pointer :: sTemp(:,:)
integer :: ierror

logical pr

if (MPI_rank == 0) then; pr=.true.; else; pr=.false.; endif

if (pr) then
!LN  write (*,'(A)',advance='no') "Generating equivalent clusters." 
  if (present(iRef)) then; write(*,'(A,I4)') " Reference ", iref
                     else; endif
endif

! We also need the symmetry of the current parent lattice so read in the parent
! lattice from the lat.in file and calculate the symmetries.
!call read_lattdef(title,cetype,CErank,LV,aBas,aBasCoord,Ninterior,Rmax,eps)

! call get_site_labels(sitelabel,CE), replaced by
call get_spaceGroup_atomTypes(CE%label, CE%digit,CE%CEBas, sitelabel)
!call make_primitive(CE%pLV, sitelabel, CE%pBas, .false., CE%eps, removed)
call get_spaceGroup(CE%pLV, sitelabel, CE%pBas, rot, shift, .false.,CE%eps)
if (CE%surface%isSurf) then
  call rm_3d_operations(CE%pLV, rot, shift, CE%eps)
!  if (.not. CE%surface%isSymmetric) stop "Rethink rm_3d_operations for non symmetric slabs"
endif

call matrix_inverse(CE%pLV,invLV,err)
if (err) stop "ERROR: generate_eqvClusters: input vectors are coplanar"

nOps = size(rot,3)
nCl = size(clusters)    ! total number of clusters, as read in from a J files, including the constant term (cluster #0)

! We need these "dummy" allocations to avoid problems for iCl == 1
! iteration below---can't deallocate if it hasn't be allocated first
allocate(diff(1,1),rotCluster%vertex(1,1)) ! 
allocate(permute(1),ps(1)) 
allocate(rotPts(1,1))


allocate(eqvClusters(0:nCl-1))   ! for each cluster there will be equivalent clusters

! setup the constant term, cluster #0
allocate(eqvClusters(0)%Rep(1))   ! has only 1 representativ
allocate(eqvClusters(0)%Rep(1)%vertex(3,0))
allocate(eqvClusters(0)%Rep(1)%label(0))
allocate(eqvClusters(0)%Rep(1)%s(0))
eqvClusters(0)%Rep(1) = clusters(0)
eqvClusters(0)%count  = 1
eqvClusters(0)%damping= clusters(0)%damping
eqvClusters(0)%nV     = clusters(0)%nV

! setup all other clusters
do iCl = 1, nCl-1
   nV = size(clusters(iCl)%vertex,2)

   deallocate(diff)
   deallocate(rotPts)
   deallocate(rotCluster%vertex)
   deallocate(permute)
   allocate(diff(3,nV), rotPts(3,nV))
   allocate(rotCluster%vertex(3,nV), rotCluster%s(nV), permute(nV), rotCluster%label(nV))
   nstore = 0
   allocate(tempStore(nOps))
   do iOps=1,nOps;allocate(tempStore(iOps)%vertex(3,nV),tempStore(iOps)%label(nV),tempStore(iOps)%s(nV)); enddo

   do isym = 1,nOps  ! Loop over all the symmetry operations
      stored            = .false.
      rotCluster%s     = clusters(iCl)%s
      rotCluster%vertex = matmul(rot(:,:,isym),clusters(iCl)%vertex) ! Rotate the figure
      do iV = 1, nV  ! Add the non-symmorphic shift, if it exists
         rotCluster%vertex(:,iV) = rotCluster%vertex(:,iV)+shift(:,isym)
      enddo
      rotPts = rotCluster%vertex
      call sort_figure_vertices(rotCluster,CE%eps,sortS=.true.)  ! sort the vertices and with them their s-vectors

      ! Move the figure so that the first vertex is inside parent cell 0
      tVec = rotCluster%vertex(:,1) ! Store first vertex in temp variable
      call bring_into_cell(tVec,invLV,CE%pLV,CE%eps) ! Move it inside parent cell 0
      shiftVec = rotCluster%vertex(:,1) - tVec ! Find the needed shift
      do iV = 1, nV ! Shift each vertex
        rotCluster%vertex(:,iV) = rotCluster%vertex(:,iV) - shiftVec 
      enddo
     
      do istore = 1,nstore ! Loop over all the instances of _this_ figure that we have stored already
                           ! old obsolete comment: (therefore, we don't have to check the s-vectors, since they are already known to be unique
                           ! for THIS figure)
                           ! new comment (03/05/12 tk): the old comment is *not* true; in fact it is *wrong*. We *have to* check if the s-vector
                           !                            is equal or not, because it is sorted together with the vertices. (see above).
                           !                            If we don't do this, a bcc-NN cluster with s=(1,2) only gets 4 representatives, just
                           !                            like the bcc-NN cluster with s=(1,1) or s=(2,2). However, there are 8 representatives!
         diff = rotCluster%vertex - tempStore(istore)%vertex 
         ! Check that each vector in "diff" is identical (necessary but insufficient condition)
         diffcheck = .true.
         do iV = 2, size(diff,2)
            if (.not. equal(diff(:,1),diff(:,iV),CE%eps)) then
               diffcheck = .false.; exit
            endif
         enddo
         if (.not. diffcheck) cycle
         lattCoordCheck = matmul(invLV,diff(:,1))
         if (equal(lattCoordCheck, (/0._dp,0._dp,0._dp/),CE%eps)) then
         if (all(rotCluster%s==tempStore(istore)%s)) then  ! also check for the s-vectors (see comment above)
            stored = .true.
            exit
         endif
         endif

      enddo ! End loop over already-stored instances of this figure

      if (.not. stored) then
         nstore = nstore + 1
         ! until we know how many symmetrically-equivalent figs there
         ! are, we can't allocate the number of Reps in eqvFig 
         tempStore(nstore)%vertex = rotCluster%vertex
         tempStore(nstore)%s      = rotCluster%s
         tempStore(nstore)%avgR   = clusters(iCl)%avgR
         tempStore(nstore)%damping= clusters(iCl)%damping
         tempStore(nstore)%nV     = nV

         do iV = 1, nV
            do iBas = 1, CE%nBas
               diffVec = tempStore(nstore)%vertex(:,iV) - CE%pBas(:,iBas)
               lattCoordCheck = matmul(invLV,diffVec)
               if (equal(lattCoordCheck, nint(lattCoordCheck),CE%eps)) then
                  tempStore(nstore)%label(iV) = iBas
               endif
            enddo
         enddo
         if (any(tempStore(nstore)%label==0)) then
            print *,"Didn't find a label for a figure: generate_eqvClusters"
            write(*,'("Clusters #:",i2)') iCl
            write(*,'("nstore: ",i2)') nstore
            write(*,'("# of vertices: ",i2)') tempStore(nstore)%nV
            write(*,'("vertices: ",3(f7.3,1x))') (tempStore(nstore)%vertex(:,iV),iV=1,nV)
            
            write(*,'("Labels: ",20(i1,1x))') tempStore(nstore)%label
            stop "ERROR: Didn't find a label for a figure generate_eqvClusters"
         endif

      endif
   enddo ! End of loop over the symmetry operators
   ! Now that we know how many symmetry brothers there are, allocate eqvClusters(iCl)%Rep and store them
   allocate(eqvClusters(iCl)%Rep(nstore))
   eqvClusters(iCl)%damping = clusters(iCl)%damping     ! damping is also stored in the top-level 
                                                        ! of eqvClusters and not only in the level %Rep for later convenience
   eqvClusters(iCl)%nV      = clusters(iCl)%nV          ! ditto
   eqvClusters(iCl)%Rep(1:nstore) = tempStore(1:nstore)
   eqvClusters(iCl)%count = nstore

   deallocate(tempStore)
enddo ! iCl, End of loop over the list of symmetrically-inequivalent figures

ENDSUBROUTINE generate_eqvClusters


!********************************************************************************
subroutine setup_damping(clusters,GA)
type(figure), intent(inout)    :: clusters(0:)   
type(GA_variables), intent(in) :: GA

integer iCl, nCl
integer nV

nCl = size(clusters)
do iCl = 0, nCl - 1
  nV = clusters(iCl)%nV
  clusters(iCl)%damping = GA%dampingC(nV) * clusters(iCl)%avgR**GA%dampingLambda(nV)
enddo
end subroutine setup_damping


!********************************************************************************
! --> IN: 
!            str  :  input structures
!            CE   :  CE data
! <-- OUT:
!            str  :  updated input structures
!                    This routine adds up all the reference energies of the structures
!                    and puts it into str%ce_refEnergy
!                    Furthermore, it sets the fit energy str%fitEnergy.
!
subroutine setup_reference_energy(str,CE)
type(crystal), pointer         :: str(:)   ! intent(INOUT)
type(CE_variables), intent(in) :: CE

integer :: iStr, nStr, iRef, nRef

nRef = CE%nreferences
nStr = size(str)

str(:)%ce_refEnergy = 0.0_dp
str(:)%fitEnergy    = str(:)%energy

if (nRef>0) then
  print *
  write(*,*) "Adding up all reference energies..."
  do iStr=1,nStr
    str(iStr)%ce_refEnergy = sum(str(iStr)%reference(:)%ce_energy)
    str(iStr)%fitEnergy    = str(iStr)%energy - str(iStr)%ce_refEnergy
  end do

  ! Output reference energy information
  print *
  write(*,*) "Reference energy information:"
  print *
  write (*, '(A31, A13, A13, A36)') "structure title |","Edft/atom |","ERef/atom |", "single References E/atom (m=manual)"
  write (*, '(A)')         "---------------------------------------------------------------------------------------------"
  do iStr=1,nStr
    write(*,'(A30,1x,F11.5,2x,F11.5,A3)',advance='no') str(iStr)%title,str(iStr)%energy,str(iStr)%ce_refenergy,"| "
    do iRef=1,nRef
      write(*,'(F10.5,1x)',advance='no') str(iStr)%reference(iRef)%ce_energy
      if (str(iStr)%reference(iRef)%automaticReference) then; write(*,'(A1)',advance='no') " "; else; write(*,'(A1)',advance='no') "m"; endif
    enddo
    write(*,*)
  enddo
endif
end subroutine setup_reference_energy



subroutine check_for_duplicate_atoms(str,iAt,jAt,eps,limit)
use compare_structures

type(crystal), intent(in) :: str
integer, intent(out)      :: iAt, jAt
real(dp), intent(in)      :: eps
integer, intent(in), optional :: limit

integer nAt
real(dp) :: conVec(3)

nAt = size(str%pos,2)
if (present(limit)) nAt = limit

do iAt = 1, nAt ! Loop over each atom
  do jAt = iAt+1, nAt ! Loop over all other atoms ahead of iAt in the list
  ! Check that the i-th and j-th atoms aren't connected by a lattice vector
    conVec = str%pos(:,iAt) - str%pos(:,jAt)
    if (is_lattice_point(str%LV, conVec, eps)) then
      return
    endif
  enddo
enddo

iAt = 0
jAt = 0

end subroutine check_for_duplicate_atoms

ENDMODULE setup_cluster_expansion

