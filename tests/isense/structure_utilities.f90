MODULE structure_utilities
use num_types
use ce_types
use enumeration_types, only: maxLabLength
use compare_structures
use numerical_utilities, only : equal
use vector_matrix_utilities, only : matrix_inverse
use utilities_module, only: ralloc
use sort
implicit none

!------------------------------------------------------------------------------------------
! Definition of a "bad" global variable, that is really "good". We need the
! Chebychev's everywhere, so a lookup table is really helpful and boosts performance. 
!
! chebychev_lookup
! The lookup table for the Chebychev values. This lookup table is always filled, even
! in the useChebychevs=F case, where it represents the identity function.
real(dp), pointer  :: chebychev_lookup(:,:)   
!------------------------------------------------------------------------------------------


public create_spin_lookup_table, prepare_input_structures, &
     & calculate_energy, get_SpinRank_mapping,&
     & get_gTab_correlations, &
     & chebychev, setup_chebychev, find_g, find_x, mod_g3, create_structure, planes_from_structure, mirror_at_plane, get_dset_label, &
     & setup_reference_structure, get_GSs, get_lowEstrucs, calc_structure_DHf, &
     & write_GSList, &
     & get_pureElements, &
     & get_ilabeling, get_ilabeling_x, set_intersection, &
     & map_labeling_position_2_dvector_index, &
     & is_dvector_removed_in_reference_ce
private

CONTAINS





SUBROUTINE planes_from_structure(str, nAt, useRSpaceFile, RSpaceFile, norm_vec, basis_vec)

type rspaceformat
real(dp)    :: pos(3)
integer(si) :: spin
end type rspaceformat

type(crystal), intent(in) :: str
real(dp), intent(in) :: norm_vec(:,:), basis_vec(:,:)
integer :: nPlane, iPlane

integer, intent(in) :: nAt
integer :: iAt, i
integer, allocatable :: nInPlane(:)
integer, allocatable :: j(:), k(:), l(:), pfh(:), dfh(:)
real(dp) :: x(3), prod
character(12), allocatable :: pname(:), dname(:)    ! Pov filename, Data filename
real(dp), allocatable :: avgR(:,:)

logical, intent(in) :: useRSpaceFile
character(30)       :: RSpaceFile
type(rspaceformat) :: rspacedata
integer :: datalength


!nAt = size(str%pos,2)

nPlane = size(norm_vec,2)
if (nPlane /= size(basis_vec,2)) stop "Bad programming"
allocate(pname(nPlane),dname(nPlane),j(nPlane),k(nPlane),l(nPlane),pfh(nPlane),dfh(nPlane))
allocate(nInPlane(nPlane),avgR(3,nPlane))
do iPlane=1,nPlane
  j(iPlane) = (iPlane/100)+1
  k(iPlane) = ((iPlane - 100 * (iPlane/100))/10)+1
  l(iPlane) = iPlane -100*(iPlane/100)-10 * ((iPlane - 100 * (iPlane/100))/10) + 1
   
  pname(iPlane) =  "plane" //"0123456789" (j(iPlane):j(iPlane))//"0123456789" (k(iPlane):k(iPlane))//"0123456789" (l(iPlane):l(iPlane))// ".pov"
  dname(iPlane) =  "plane" //"0123456789" (j(iPlane):j(iPlane))//"0123456789" (k(iPlane):k(iPlane))//"0123456789" (l(iPlane):l(iPlane))// ".dat"
  pfh(iPlane)    = 260 + iPlane -1
  dfh(iPlane)    = 460 + iPlane -1
  open(pfh(iPlane), file=pname(iPlane))
  open(dfh(iPlane), file=dname(iPlane))
  write (pfh(iPlane), '(A)') '#include "colors.inc"'
  write (pfh(iPlane), '(A)') '#include "finish.inc"'
  write (pfh(iPlane), '(A)') '#include "textures.inc"'
  write (pfh(iPlane), '(A)') '#include "stones.inc"'
  write (pfh(iPlane), '(A)') '#include "woods.inc"'
  write (pfh(iPlane), '(A)') "camera { orthographic right 1*x angle"
  write (pfh(iPlane), '(A,3(F10.4,A))') "  location <",basis_vec(1,iPlane)+10*norm_vec(1,iPlane),",", &
       basis_vec(2,iPlane)+10*norm_vec(2,iPlane),",",basis_vec(2,iPlane)+10*norm_vec(3,iPlane)  ,">"
  write (pfh(iPlane), '(A)') "  look_at <0.0,0.0,0.0>   rotate <0,0,0>}"
  write (pfh(iPlane), '(A)') "global_settings { ambient_light rgb <1, 1, 1>}"
  write (pfh(iPlane), '(A,3(F10.4,A))') "light_source{ <",20*norm_vec(1,iPlane),",",20*norm_vec(2,iPlane),",",20*norm_vec(3,iPlane),"> color White}"
  write (pfh(iPlane), '(A)') "background { color White }"
  write (pfh(iPlane), '(A)') "#declare atom0=sphere {<0, 0, 0>, 0.5 texture {pigment {color rgb <0.9,0.9,0.9>} finish { phong 0.75 }}}" ! was Al
  write (pfh(iPlane), '(A)') "#declare atom1=sphere {<0, 0, 0>, 0.2 texture {pigment {color rgbt <0.5,0.5,0.5,0.9>} finish {phong 0.75 }}}" ! was Vac
  write (pfh(iPlane), '(A)') "#declare atom2=sphere {<0, 0, 0>, 0.5 texture {pigment {color rgb <0.2,0.2,0.2>} finish { phong 0.75 }}}" ! was Ni
  write (pfh(iPlane), '(A)') "#declare MCplane = union {"


  write (dfh(iPlane), "(A)") "<Number of Atoms>"
  write (dfh(iPlane), '(A)') "Cu Pt | x, y, z"

  avgR  (:,iPlane) = 0.0
  nInPlane(iPlane) = 0

enddo


if (.not. useRSpaceFile) then 
  do iAt = 1, nAt
    do iPlane=1,nPlane
     x = str%pos(:,iAt) - basis_vec(:,iPlane)
     prod = 0
     do i=1, 3
        prod = prod + x(i)*norm_vec(i,iPlane)
     enddo
     if (equal(prod, 0.0_dp, str%eps)) then
        avgR(:,iPlane) = avgR(:,iPlane) + str%pos(:,iAt)
        nInPlane(iPlane) = nInPlane(iPlane) + 1
        if (nInPlane(iPlane) <= 10) then ! the first 10 atoms
          write(*,'(A,I4,A,I4,A,I4,A,3(F7.3,1x))') "Plane: ",iPlane," #atom: ", nInPlane(iPlane), " spin: ", str%spin(iAt), " pos: ", str%pos(:,iAt)
        endif
!        write (*, '(3F16.4, I3)') str%pos(:,iAt), str%spin(iAt)

        ! write to data file
        write(dfh(iPlane),'(I2,5x,3(F10.5,1x))')  str%spin(iAt), str%pos(:,iAt)

        ! write to povray file
        select case (str%spin(iAt))
        case(-1) 
           write (pfh(iPlane), '(A,3(F10.3,A) )') "object{ atom0 translate <", str%pos(1,iAt),&
                & ","  , str%pos(2,iAt),","  , str%pos(3,iAt),">}"
        case(0)
           write (pfh(iPlane), '(A,3(F10.3,A) )') "object{ atom1 translate <", str%pos(1,iAt),&
                & ","  , str%pos(2,iAt),","  , str%pos(3,iAt),">}"
        case(1)
           write (pfh(iPlane), '(A,3(F10.3,A) )') "object{ atom2 translate <", str%pos(1,iAt),&
                & ","  , str%pos(2,iAt),","  , str%pos(3,iAt),">}"
        case default
           stop "found wrong spin!"
        endselect

     endif
    enddo ! plane loop
  enddo ! atom loop
else ! read structure data from file
  inquire(iolength=datalength) rspacedata
  open(100,file=RSpaceFile,access='direct',form='unformatted',status='old',recl=datalength)
  do iAt = 1, nAt
    if (mod(iAt,10000000)==0) write(*,'(A,F5.1,A)') "Finding plane: ", iAt*100.d0/nAt, "% done"
    do iPlane=1,nPlane
     read(100,rec=iAt) rspacedata
     x = rspacedata%pos(:)- basis_vec(:,iPlane)
     prod = 0
     do i=1, 3
        prod = prod + x(i)*norm_vec(i,iPlane)
     enddo
     if (equal(prod, 0.0_dp, str%eps)) then
        avgR(:,iPlane) = avgR(:,iPlane) + rspacedata%pos(:)
        nInPlane(iPlane) = nInPlane(iPlane) + 1
        if (nInPlane(iPlane) <= 10) then ! the first 10 atoms
          write(*,'(A,I4,A,I4,A,I4,A,3(F7.3,1x))') "Plane: ",iPlane," #atom: ", nInPlane(iPlane), " spin: ", rspacedata%spin, " pos: ", rspacedata%pos(:) 
        endif
 !       write (*, '(3F16.4, I3)') rspacedata%pos(:), rspacedata%spin
        select case (rspacedata%spin)
        case(-1) 
           write (pfh(iPlane), '(A,3(F10.3,A) )') "object{ atom0 translate <", rspacedata%pos(1),&
                & ","  , rspacedata%pos(2),","  , rspacedata%pos(3),">}"
        case(0)
           write (pfh(iPlane), '(A,3(F10.3,A) )') "object{ atom1 translate <", rspacedata%pos(1),&
                & ","  , rspacedata%pos(2),","  , rspacedata%pos(3),">}"
        case(1)
           write (pfh(iPlane), '(A,3(F10.3,A) )') "object{ atom2 translate <", rspacedata%pos(1),&
                & ","  , rspacedata%pos(2),","  , rspacedata%pos(3),">}"
        case default
           stop "found wrong spin!"
        endselect
     endif
    enddo ! plane loop
   enddo ! atom loop
   close(100)
endif

do iPlane=1,nPlane
  avgR(:,iPlane) = (-1)* avgR(:,iPlane) / nInPlane(iPlane)
  write (pfh(iPlane), '(A)') "}"
!write (26, '(A,3(F10.4,A))') "#declare shiftpos = <",avgR(1),",",avgR(2),",",avgR(3),">;"
  write (pfh(iPlane), '(A,3(F10.4,A))') "object { MCplane translate <",avgR(1,iPlane),",",avgR(2,iPlane),",",avgR(3,iPlane),">  scale 1.0 } "
  close(pfh(iPlane))
  close(dfh(iPlane))
enddo

ENDSUBROUTINE planes_from_structure











!*******************************************************************************
! Purpose: Calculate the correlation by conversion of vertices to g-space
!
!
! IMPORTANT:
! MAKE SURE THAT YOU'VE SETUP THE CHEBYCHEV LOOKUP TABLE BEFORE CALCULATING CORRELATIONS
!
!

!<parameter name="gTab">The table of spin values in g-space for the entire monte-carlo cell. Indices: { p, q, r, iPBas }.</parameter>
!<parameter name="L">The left matrix that was used to compute the Smith Normal Form of the parent lattice vectors @CREF[param.D].</parameter>
!<parameter name="D">The Smith Normal Form of the parent lattice vectors.</parameter>
!<parameter name="invA">A is the column form of the original parent lattice vectors. invA is its inverse.</parameter>
SUBROUTINE get_gTab_correlations(CE, gTab, L, D, invA, cluster, PI, mustGetG, atG3)
type(CE_variables), intent(in) :: CE
integer(kind=1), pointer       :: gTab(:,:,:,:)     ! { p, q, r, iPBas }
integer, intent(in)            :: L(3,3), D(3,3)    ! L can be a dummy if mustGetG=.false.
real(dp), intent(in)           :: invA(3,3)         ! invA can be a dummy if mustGetG=.false.
type(figRep), intent(inout)    :: cluster(0:)
real(dp), pointer              :: PI(:)     ! { clusters# }, intent(OUT)
logical, intent(in)            :: mustGetG  ! must this routine calculate the g values of the clusters? (if not, L and invA can be a dummies)
integer, intent(in), optional  :: atG3(3)   ! at which gpoint do we want to evaluate the correlation? if not given: all points are
                                            ! taken into account
integer p, q, r

integer iCl, iV, iRep, iPbas, i, j, counter
integer nAt, nCl, nV
integer g3(3), n(3), g4(4)
real(dp) :: eps!, Sconc, li, thisli ! epsilon, 2x-1, likelihood index, cluster-specific li
real(dp) :: vertex(3)

eps = CE%eps
nAt = size(gTab)
nCl = size(cluster)


allocate(PI(0:nCl-1)) ! constant + number of clusters (including s-vectors)
PI    = 0._dp
!................................................................................
! set correlation for constant, cluster #0
PI(0) = 1._dp        
if (mustGetG) allocate( cluster(0)%Rep(1)%gPt(4,0) )  ! this needs to be allocated for ALL clusters (cf. below)

!................................................................................
! calculate correlations for all other clusters
do iCl = 1, nCl-1 ! Loop over each cluster
   do iRep = 1, size(cluster(iCl)%Rep) ! loop over each representation of current cluster
      nV = cluster(iCl)%Rep(1)%nV ! abbreviation for the number of vertices, which is the same for all
                                  ! representatives.
                                  ! NB: Don't put this statement outside the iRep loop since Rep(1) may
                                  !     not be defined, and the loop will be skipped.

      ! Set up the 4-vectors for the cluster vertices 
      if (mustGetG) then ! get the g-coordinates of the clusters. 
                         ! Will do this only for "normal" structures. In the gss case (i.e. for many many stuctures with
                         ! the same HNF) those will be calculated beforehand, so we don't have to spend time doing it here.

         allocate(cluster(iCl)%Rep(iRep)%gPt(4,nV))
         do iV = 1, nV ! loop over each vertex
            iPbas  = cluster(iCl)%Rep(iRep)%label(iV)
            vertex = cluster(iCl)%Rep(iRep)%vertex(:,iV)
            call find_g(vertex,invA,L,D,CE%pBas,iPbas,eps,g3)
            g4 = (/g3, iPbas/)
            cluster(iCl)%rep(iRep)%gPt(:,iV) = g4
         enddo
      end if
      ! Now calculate correlations for this cluster by looping over all the
      ! entries of the g table (Thanks again, Rod)
      if (present(atG3)) then ! evaluate correlations at a certain point
         p = atG3(1)
         q = atG3(2)
         r = atG3(3)
         PI(iCl) = PI(iCl) + spin_product(p,q,r,cluster(iCl)%Rep(iRep),nV,gTab,D)
       
      else ! evaluate correlations for the whole structures, i.e. all points
         do p = 0, D(1,1)-1
            do q = 0, D(2,2)-1
               do r = 0, D(3,3)-1
                  PI(iCl) = PI(iCl) + spin_product(p,q,r,cluster(iCl)%Rep(iRep),nV,gTab,D)       
               enddo ! end r loop
            enddo ! end q loop
         enddo ! end p loop
       endif
    enddo ! End do over equivalent reps of the cluster
enddo ! End loop clusters

! normalize all clusters, but NOT the constant
forall(iCl=1:nCl-1)
   PI(iCl) = PI(iCl) / real( cluster(iCl)%count*D(1,1)*D(2,2)*D(3,3) , dp ) ! normalize
end forall

!open(33, file="correlations.out", position="append")
!do iV = 0,(size(PI)-1)
!   write(33,'("Fig #",i4," Corr= ",f10.6)') iV, PI(iV)
!enddo
contains


!................................................................................
pure function spin_product(p,q,r,clusterRep,nV,gTab,D)
integer, intent(in)       :: p, q, r
type(figure), intent(in)  :: clusterRep
integer, intent(in)       :: nV
integer(kind=1), pointer  :: gTab(:,:,:,:)  ! { p, q, r, iPBas }
integer, intent(in)       :: D(3,3)

real(dp) spin_product

integer iV, g3(3), g4(4)

spin_product = 1._dp

do iV = 1, nV ! Loop over vertices
   g3 = (/p,q,r/)+clusterRep%gPt(1:3,iV)
   ! Bounce every vertex back into structure cell 0
   g3 = mod(g3,(/D(1,1),D(2,2),D(3,3)/))
   g4 = (/g3, clusterRep%gPt(4,iV)/)
   spin_product = spin_product * chebychev(  gTab(g4(1),g4(2),g4(3),g4(4))   ,  clusterRep%s(iV) )
enddo
end function spin_product
!................................................................................


END SUBROUTINE get_gTab_correlations
!****************************************************************************************************




!***************************************************************************************************
! This routine takes the CE rank (given by k) and maps it onto the spin variables, -k/2..k/2,
! skipping 0 where necessary (for even rank CEs), and maps it onto rank space, 1..k, as well
SUBROUTINE get_SpinRank_mapping(k, spinSpace, rankSpace)
integer, intent(in) :: k
integer, pointer :: spinSpace(:), rankSpace(:)
integer ik, i ! Counters

allocate(spinSpace(k))
allocate(rankSpace(-(k-mod(k,2))/2:(k-mod(k,2))/2))  !Handle k=odd and k=even cases together
do ik = 0, k-1
   i =  ik-k/2 ! convert 0..k-1 label to spin variable -k/2..k/2
   if(mod(k,2)==0) then
      if (.not. i<0) i = i + 1 ! skip 0 as a spin if k is even
   endif
   spinSpace(ik+1) = i 
   rankSpace(i) = iK + 1 ! Add one to make list go 1..k
enddo
ENDSUBROUTINE get_SpinRank_mapping

!***************************************************************************************************
pure FUNCTION calculate_energy(J,PI)
real(dp) calculate_energy
real(dp), intent(in)      :: PI(0:)
real(dp), intent(in)      :: J(0:)


!energy = effI%empty

calculate_energy = dot_product(J(0:),PI(0:))

ENDFUNCTION calculate_energy

!***************************************************************************************************
subroutine get_ilabeling(labeling,ilabeling)
character(maxLabLength), intent(in) :: labeling
integer, pointer          :: ilabeling(:)

integer nl, il

nl=len_trim(labeling)
allocate(ilabeling(nl))
forall (il=1:nl); ilabeling(il)=ichar(labeling(il:il))-48; end forall

end subroutine get_ilabeling

subroutine get_ilabeling_x(nUC,nD,k,nC,dC,ilabeling,x,nTyp,xRel)
integer, intent(in) :: nUC, nD ! number of unit cells in the labeling, number of parent dset members
integer, intent(in) :: k ! CE rank
integer, intent(in) :: nC ! number of CE couplings
integer, intent(in) :: dC(:) ! { dset# }: d-vector to coupling map (i.e. which d-vector belongs to which CE)
integer, intent(in)  :: ilabeling(:) ! integer representation of labeling
real(dp), pointer :: x(:,:) ! out: concentration for this ilabeling for each CE
logical, intent(in) :: xRel(:) ! { dset# }: count only labels that are relevant for concentration
integer, pointer :: nTyp(:,:)

integer :: xRelSitesPerUC(nC) ! { coupling# }: number of x relevent sites per unit cell
integer iEl, iD, iC
integer :: labelpos2dsetlabel(nUC*nD)


allocate(x(k,nC)); x=0._dp
allocate(nTyp(k,nC)); nTyp= 0

! this array tells us the dset number to which a position in the labeling corresponds to,
! e.g. in a labeling xxyyzz for nUC=2 (total of 2 unit cells) and nD=3 (3 dset members), 
! all x belong to dset 1, y to 2, z to 3.
labelpos2dsetlabel =  map_labeling_position_2_dvector_index(nUC=nUC,nD=nD)

forall (iC=1:nC)
  xRelSitesPerUC(iC) = count(xRel .and. dC==iC)
  forall (iEl=1:k)
    nTyp(iEl,iC) = sum((/ (  count(  pack(ilabeling,labelpos2dsetlabel==iD)==iEl-1  ) , iD=1,nD ) /) , xRel .and. dC==iC)
  end forall
  x(:,iC) = nTyp(:,iC)/real(nUC*xRelSitesPerUC(iC),dp)
end forall


!forall (iC=1:nC)
!  nTyp(:,iC) = x(:,iC) * xRelSitesPerUC(iC) * nUC
!end forall

end subroutine get_ilabeling_x


subroutine get_ilabeling_spin(k,ilabeling,spin)
integer, intent(in) :: k
integer, pointer :: ilabeling(:) ! in
integer, pointer :: spin(:) ! out

integer nl, il

nl = size(ilabeling)
allocate(spin(nl))

forall (il=1:nl); spin(il) = ilabeling(il) - k/2; end forall
! if even CE: spin /= 0! ==> "spin" gets adjusted
if (mod(k,2)==0) then; where(spin>=0); spin = spin+1; end where; endif

end subroutine get_ilabeling_spin


!***************************************************************************************************
! Given a point, find out which interior point of the *parent lattice* it is
FUNCTION get_pIP(point, pLV, pBas, eps)
real(dp), intent(in) :: point(3) ! the point/vertex we're mapping into g table
real(dp), intent(in) :: pLV(3,3), pBas(:,:) ! parent lattice vectors and site basis vectors
integer get_pIP
real(dp) :: eps

integer nPar 
integer iPar
logical checkpoint

checkPoint = .false.
nPar = size(pBas,2)
do iPar = 1, nPar ! Try shifting by each parent basis vector to take care of the
   ! general case where the parent isn't a primitive cell
   if (is_lattice_point(pLV,point-pBas(:,iPar),eps)) then
      checkPoint = .true.; exit
   endif
enddo
if (.not. checkPoint) then
   write(*,'("The given point/vertex is not on the parent lattice.")') 
   write(*,'("Error in get_parent_interior_point, at point:",3f7.3)') point
   stop
endif
get_pIP = iPar
ENDFUNCTION get_pIP

!***************************************************************************************************
! Find the coordinates, in the orthorhombic lattice, of a given point. (That is, fofx will be the
! coordinates of the given point, given in multiples of the parent lattice.)
! (See Rod's write-up of Nov. 10)
SUBROUTINE find_g(point,invA,L,D,pBas,iPbas,eps,gpoint)
real(dp), intent(in) :: point(3) ! the point/vertex we're mapping into g table
integer, intent(in), dimension(3,3) ::  L, D ! left transformation and Smith Normal Form
integer, intent(out) :: gpoint(3)
real(dp), intent(in) :: eps, invA(3,3)
real(dp), pointer :: pBas(:,:) ! (input only)
integer,intent(in) :: iPbas ! the ordinal number of the site in the parent cell

real(dp) :: x(3), tx(3)
integer ic, i, z(3)

x = point - pBas(:,iPbas) ! Get the coordinates of the cell rather than the point
tx = matmul(invA,x)       ! Get the lattice coordinates of the cell
if (.not. equal(nint(tx)-tx, 0._dp,eps)) then
   write(*,'("ERROR in find_g: Atom/vertex is not on a lattice point")')
   write(*,'(A,3(F6.3,1x),A,i2)')   "  real space point / dset label ", point(:),"/",iPbas
   write(*,'(A,3(F6.3,1x))')        "  real space lattice point:     ", x(:)
   write(*,'(A,3(F6.3,1x))')        "  real space z point = A^-1 r:  ", tx
   write(*,'(A,3(F6.3,1x),A,F6.3)') "  deviation from int / allowed: ", nint(tx)-tx,"/", eps
   write(*,*) " A^-1: ", invA
!   write(*,'("Diff: ",3(g12.6,1x)," EPS: ",g9.4)') nint(tx)-tx, eps
!   write( *, '(3(f6.3,1x),i2)') point(:), iPbas
!   print *,point - pBas(:,iPbas)
   stop 
end if
z = matmul(L,nint(tx)) ! Get the coordinates of the cell in coordinates of the "orthorhombic" rep. of supercell

call mod_g3(z, D, gpoint)  ! Move the coords of the cell inside the first supercell

END SUBROUTINE find_g  



!********************************************************************************
! The main goal is to map a given g3 gpoint back to the first supercell. However, a simple "mod"
! does not work due to negative numbers. In that case, we have to add a superlattice
! vector for all g3 coordinates to be >0.
! 
! Stand-alone routine since 26/7/10, tk
! Prior, the code was included in the find_g subroutine.
!
subroutine mod_g3(g3in, D, g3out)
integer, intent(in)  :: g3in(3)
integer, intent(in), dimension(3,3) ::  D ! Smith Normal Form
integer, intent(out) :: g3out(3)

integer i, ic

g3out = mod(g3in,(/D(1,1),D(2,2),D(3,3)/)) ! Move the coords of the cell inside the first supercell (except if negative)
do i = 1,3 ! Loop over x,y,z  ! It's possible that the coordinates of the cell are negative. If so, add one sLV.
   ic = 0;
   do while (g3out(i)<0) ! Moving the point to the inside of a single orthorhombic tile
      ic = ic + 1;      ! This should only have to happen once. If not, there is a math mistake in the algorithm
      if (ic>1) then;
         write(*,*) "The shaft is ray-shielded so you'll have to use proton torpedoes."
         write(*,*) "Please email gus.hart@gmail.com."
         stop
      endif
      g3out(i) = g3out(i) + D(i,i) ! Add one supercell lattice vector to the current coordinate of the cell.
   enddo
enddo

end subroutine mod_g3



!***************************************************************************************************
! Find the realspace coordinates, given a certain g_vector of our lookup table
! This is the inverse routine of find_g
SUBROUTINE find_x(g,iPbas,A,invL,pBas,r)
use vector_matrix_utilities, only : matrix_inverse

integer, intent(in) :: g(3)
integer, intent(in) :: iPbas
real(dp), pointer :: pBas(:,:) ! (input only)
real(dp), intent(in) :: A(3,3), invL(3,3)
real(dp), intent(out) :: r(3)

real(dp) :: tx(3), x(3)!, invL(3,3), Lreal(3,3)

tx = matmul(invL, real(g))
x = matmul(A, tx)

r = x + pbas(:,iPbas)

ENDSUBROUTINE find_x

!***************************************************************************************************
! This routine takes in a g_table and gives back the coordinates of the structure.
! inverse routine of create_spin_lookup_table
SUBROUTINE create_structure(gTab, invL, str, nAt, write2file, rspacefile)
use vector_matrix_utilities, only : matrix_inverse

type rspaceformat
real(dp)    :: pos(3)
integer(si) :: spin
end type rspaceformat

integer(1), pointer :: gTab(:,:,:,:)
type(crystal), intent(inout) :: str
real(dp), intent(in) :: invL(3,3)
integer, intent(out) :: nAt

integer :: ix, iy, iz, ib, iAt
integer :: nx, ny, nz, nb
integer :: istat

logical, intent(out)  :: write2file
logical :: useOldFile
character(30), intent(inout) :: rspacefile
type(rspaceformat) :: rspacedata
integer :: datalength

nx = size(gTab,1)
ny = size(gTab,2)
nz = size(gTab,3)
nb = size(gTab,4)
nAt = nx * ny * nz * nb
allocate(str%pos(3,nAt), str%spin(nAt),stat=istat)
if (istat /= 0) then
  allocate(str%pos(3,1), str%spin(1),stat=istat)
  if (istat /= 0) stop "ERROR: in memory allocation"
  write2file = .true.
  inquire(iolength=datalength) rspacedata
  open(100,file=rspacefile,access='direct',form='unformatted',status='new',recl=datalength,iostat=istat)
  if (istat /= 0) then
    print *, "WARNING: will not override file ", rspacefile
    print *, "WARNING: the pov-file we create may be out of date"
    useOldFile = .true.
  else
    useOldFile = .false.
  endif
else
  write2file = .false.
endif 

iAt = 1

if (.not. write2file) then
  ! write to RAM
  do ib = 1, nb
    do iz = 0, nz-1
      do iy = 0, ny-1
        do ix = 0, nx-1
          str%spin(iAt) = gTab(ix,iy,iz,ib)
          call find_x((/ix,iy,iz/), ib, str%pLV, invL, str%intPt, str%pos(:,iAt))
!            write (*, '(A, I3, 3F8.4, I3)')"Atom #", iAt, str%pos(:,iAt), str%spin(iAt)
          iAt = iAt + 1
        enddo
      enddo
    enddo
  enddo
else ! write to file
  if (.not. useOldFile) then ! if an old real space file is present, we don't have to generate a new one
  do ib = 1, nb
    do iz = 0, nz-1
      do iy = 0, ny-1
        do ix = 0, nx-1
          if (mod(iAt,10000000)==0) write(*,'(A,F5.1,A)') "Creating real space structure: ", iAt*100.d0/nAt, "% done."
          rspacedata%spin= gTab(ix,iy,iz,ib)
          call find_x((/ix,iy,iz/), ib, str%pLV, invL, str%intPt, rspacedata%pos(:))
          write(100,rec=iAt) rspacedata
          iAt = iAt + 1 
        enddo
      enddo
    enddo
  enddo
  close(100) ! rspacefile
  endif
endif 
print *
ENDSUBROUTINE create_structure

!***************************************************************************************************
! Find the parent cell lattice coordinates of each atom in the basis of the superstructure.
! The "g-table" array contains the atom type (spin) of each atom in the superlattice.
SUBROUTINE create_spin_lookup_table(structure,gStruct,invA,L,D, determineLD)
use vector_matrix_utilities, only : matrix_inverse
use rational_mathematics, only : SmithNormalForm
type(crystal), intent(in)                    :: structure
integer(1),    pointer                       :: gStruct(:,:,:,:)
real(dp),      intent(out)                   :: invA(3,3)
integer,       dimension(3,3), intent(inout) :: D, L ! Smith Normal form (diagonal) and left transformation
logical,       intent(in)                    :: determineLD  

integer iAt, iPBV, nAt, nBas
real(dp) :: eps
integer, dimension(3,3) :: R ! Right transformations matrix for SNF 
integer, dimension(3,3) ::  N ! lattice coordinate rep. of the superlattice
real(dp) ::  tN(3,3)  ! Inverse of parent LVs, temporary representative of matrix N
logical err
integer, dimension(3) :: gpoint ! Lattice coordinates of x, mapping of z to orthorhombic cell

nAt = size(structure%spin) ! Number of atoms in the superstructure basis
nBas = size(structure%intPt,2)  ! Number of atoms in the basis of the parent lattice
eps = structure%eps  ! Finite precision checking parameter


call matrix_inverse(structure%pLV,invA,err)
if (err) stop "ERROR: Coplanar vectors in 'create_atomic_spin_lookup_table'"

if (determineLD) then ! do we really have to determine L and D? If yes, do so.
  tN = matmul(invA,structure%LV) ! N = A^-1*B; N is the HNF of the superlattice
  
  if (.not. equal(tN-nint(tN),0._dp,eps)) stop "ERROR: N matrix isn't integer in 'create_spin_lookup_table'"
  N = nint(tN)
  call SmithNormalForm(N , L,D,R) 
                    ! in | out
endif 

allocate(gStruct(0:D(1,1)-1,0:D(2,2)-1,0:D(3,3)-1,nBas)) ! Set up the table lattice coords of the basis atoms
gStruct = 0

do iAt = 1, nAt ! Loop over each basis atom and place it in the table
   iPBV = structure%pBVlab(iAt) ! the dset label of that atom
!   print *,"site labels",structure%pBVlab(:)
!    write(*,'("atom: ",3(f7.3,1x),2x,I3)') structure%pos(:,iAt), structure%pbvlab(iAt)
   call find_g(structure%pos(:,iAt),invA,L,D,structure%intPt,iPBV,eps,gpoint)
!   write(*,'(1(3i2,2x))') gpoint
   gStruct(gpoint(1),gpoint(2),gpoint(3),iPBV) = structure%spin(iAt) 
enddo

ENDSUBROUTINE create_spin_lookup_table



!*****************************************************************************
! - sort input structures by energy
! - calculate formation enthalpies
! - find dft groundstate line
SUBROUTINE prepare_input_structures(structures, CE, GA_data)
type(crystal), pointer :: structures(:)
type(crystal), pointer :: tempstructures(:)
type(crystal) :: tstr
type(CE_variables), intent(inout) :: CE
type(GA_variables), intent(in) :: GA_data

integer iStr, jStr, c, i, j, k
integer nStr, cend, nrdftgst

!helper variables to contain the switched data
real(dp) :: LV(3,3)   ! Lattice vectors
real(dp), allocatable :: pos(:,:) ! position of the basis atoms
integer,  allocatable :: spin(:)   ! occupation variables (-1,+1,etc...)
integer,  allocatable :: pBVlab(:) ! parent basis vector label
integer,  allocatable :: nTyp(:)   ! The number of each type of atom
integer,  allocatable :: dftgstlist(:) ! List containing the indices of the groundstate structures
integer               :: weight
character(80) :: title        ! Name of the structure
real(dp) :: pLV(3,3) ! parent lattice vectors
real(dp), allocatable :: intPt (:,:) ! Interior points of the parent lattice
real(dp) :: a0     ! Lattice constant (not really necessary but...)
real(dp) :: energy, quasibinform, form ! Energy of the structure
real(dp) :: eps ! epsilon for finite precision checking
real(dp) :: conc ! binary concentration, in ternary case: quasibinary conc.
real(dp), allocatable :: x(:), x1(:), x2(:)
real(dp) :: con(3) ! ternary concentrations
integer :: stringLength, maxStringLength

! needed for enthalpies and gst:
real(dp) :: bulkenergya, bulkenergyb, bulkenergyc, conca, concb, gstformen, nTotalAtom
integer :: firstgst, lastgst, nrlowestenergy, nAt
integer, allocatable :: lowestenergy(:), templowestenergy(:)
logical :: foundgst, founda, foundb, foundc
integer, allocatable :: adsList(:)
character(100) :: formatstring
character(80)  :: filename
character(80)  :: header(1)

! references
integer :: nref, iref

integer :: ix, nx, ic, nc
real(dp), allocatable :: pureConditions(:,:)

real(dp), pointer :: tmp(:)

nStr = size(structures)

!again those dummy allocations ...
allocate(pos(1,1))
allocate(spin(1))
allocate(pBVlab(1))
allocate(nTyp(1))
allocate(intPt(1,1))
allocate(x(1))

!!GH no longer used fitstrucnumb = nStr - CE%leftoutStr


!do iStr = 1, nStr
!   write (*, '(A20, 4F8.3, F13.6, I3)') structures(iStr)%title, structures(iStr)%x, &
!                                  structures(iStr)%energy, size(structures(iStr)%pos,2)
!end do

!!DEBUG
!maxStringLength=40
!do iStr = 1, nStr
!   stringLength=len(trim(structures(iStr)%title))
!   if (stringLength>maxStringLength) then
!     write (*, '(I3,1x,A15,A3,A22, F8.3,4x, 10(F8.3,1x))') istr,structures(iStr)%title(1:15),"...",structures(iStr)%title(stringLength-21:stringLength), structures(iStr)%energy, structures(iStr)%x
!   else
!     write (*, '(I3,1x,A40, F8.3,4x, 10(F8.3,1x))',advance='no') istr,structures(iStr)%title, structures(iStr)%energy, structures(iStr)%x
!write(*,'(3(F3.1,1x),1x,I2)'), structures(iStr)%pos(:,1),structures(iStr)%pBVlab(1) !DEBUG
!   endif
!end do

print *, " "
print *, "sorting structures by concentration"
print *, " "

nx = CE%CErank
nc = CE%ncouplings
allocate(x1(nx*nc),x2(nx*nc))

do ic=1,nc
do ix=1,nx-1
  do iStr=nStr, 1, -1
    do jStr=2, iStr
      if (ix==1 .and. ic==1) then
        if (structures(jStr-1)%x(ix,ic) > structures(jStr)%x(ix,ic)) then
          tStr = structures(jStr)
          structures(jStr)=structures(jStr-1)
          structures(jStr-1) = tStr
        end if
      else ! ix > 1
        x1 = reshape(structures(jStr-1)%x(:,:),(/1/))
        x2 = reshape(structures(jStr)  %x(:,:),(/1/))

        if ( all( abs(x1(1:(ic-1)*nx+ix-1)-x2(1:(ic-1)*nx+ix-1)) < CE%xeps) ) then ! only sort if all previous concentrations are equal
        if (structures(jStr-1)%x(ix,ic) > structures(jStr)%x(ix,ic)) then
          tStr = structures(jStr)
          structures(jStr)=structures(jStr-1)
          structures(jStr-1) = tStr
        end if
        endif
      endif ! ix check
    end do
  end do
enddo
enddo


!do iStr = 1, nStr
!   write (*, '(A20, 4F8.3, F13.6, I3)') structures(iStr)%title, structures(iStr)%x, &
!                                  structures(iStr)%energy, size(structures(iStr)%pos,2)
!!   do jStr=1,  size(structures(iStr)%pos,2)
!!      write(*, '(3F8.4)') structures(iStr)%pos(:,jStr)
!!   end do
!end do

!!DEBUG
!maxStringLength=40
!do iStr = 1, nStr
!   stringLength=len(trim(structures(iStr)%title))
!   if (stringLength>maxStringLength) then
!     write (*, '(I3,1x,A15,A3,A22, F8.3,4x, 10(F8.3,1x))') istr,structures(iStr)%title(1:15),"...",structures(iStr)%title(stringLength-21:stringLength), structures(iStr)%energy, structures(iStr)%x
!   else
!     write (*, '(I3,1x,A40, F8.3,4x, 10(F8.3,1x))',advance='no') istr,structures(iStr)%title, structures(iStr)%energy, structures(iStr)%x
!write(*,'(3(F3.1,1x),1x,I2)'), structures(iStr)%pos(:,1),structures(iStr)%pBVlab(1) !DEBUG
!   endif
!end do

print *, " "
print *, "sorting structures by energy ..."
print *, " "

do iStr=nStr, 1, -1
  do jStr=2, iStr
    x1 = reshape(structures(jStr-1)%x(:,:),(/1/))
    x2 = reshape(structures(jStr)  %x(:,:),(/1/))

    if ( all( abs(x1-x2) < CE%xeps) ) then ! only sort if all concentrations are equal
      if (structures(jStr-1)%energy < structures(jStr)%energy) then
        tStr = structures(jStr)
        structures(jStr)=structures(jStr-1)
        structures(jStr-1) = tStr
      end if
    endif
  end do
end do

maxStringLength=40
do iStr = 1, nStr
   stringLength=len(trim(structures(iStr)%title))
   if (stringLength>maxStringLength) then
     write (*, '(I3,1x,A15,A3,A22, F8.3,4x, 10(F8.3,1x))') istr,structures(iStr)%title(1:15),"...",structures(iStr)%title(stringLength-21:stringLength), structures(iStr)%energy, structures(iStr)%x
   else
     write (*, '(I3,1x,A40, F8.3,4x, 10(F8.3,1x))') istr,structures(iStr)%title, structures(iStr)%energy, structures(iStr)%x
   endif
end do

print *, " "
print *, "calculating formation energies"
print *, " "

allocate(CE%pureList(CE%CErank))

if (CE%ncouplings == 1) then

  ! get the pure elements of CE 1
  call get_pureElements(1,structures,CE%CErank,CE%ncouplings,CE%pureList)
  print *
  print *, "Pure Elements are structures #", CE%pureList

  forall(iStr=1:nStr)
    structures(iStr)%formationenergy = calc_structure_DHf(structures,iStr, &
                                                          pureEl=CE%pureList, isCE=.false.)
  end forall

elseif (CE%ncouplings == 2 .and. CE%adsorbate%isAds) then

  ! get the pure elements of CE 1 under the following concentration conditions:
  allocate(pureConditions(3,1))  ! 1 condition only:
  pureConditions(1,1) = CE%adsorbate%freeSiteRank ! this concentration ... (i.e. the concentration of the
                                                  ! free adsorbate sites)
  pureConditions(2,1) = 2        ! ... of CE 2 ...
  pureConditions(3,1) = 1        ! ... must be 1
  
  call get_pureElements(1,structures,CE%CErank,CE%ncouplings,CE%pureList,pureConditions,CE%xeps)

  print *, "Pure Elements are structures #", CE%pureList

  ! setup the DHf calculation:
  ! will need the single adsorbate energies within a structure list. So generate
  ! a temporary list, "tempstructures". This will be equal to "structures" for all structures
  ! up to nStr. Then, we append the adsorbate data:
  allocate(adsList(CE%CErank))
  allocate(tempstructures(nStr+CE%CErank))  ! there are additional CErank single adsorbates
  tempstructures(:nStr) = structures(:nStr) ! up to :nStr, "tempstructures" is the same as structures
  tempstructures(nStr+1:)%isHelpStructure = .true.  ! all additional structures are only for helping out
  j = nStr
  do i=1,CE%CErank ! fill the appended structures with information
    j = j+1

    ! - concentrations
    allocate(tempstructures(j)%x(CE%CErank,CE%ncouplings))
    allocate(tempstructures(j)%nTyp(CE%CErank,CE%ncouplings))
    write(tempstructures(j)%title,'(A,I1)') "Single Adsorbate #",i
    tempstructures(j)%x         = 0._dp
    tempstructures(j)%nTyp      = 0
    tempstructures(j)%x(i,2)    = 1._dp
    tempstructures(j)%nTyp(i,2) = 1

    ! - energy (declare the single adsorbate energies)
    tempstructures(j)%energy          = CE%adsorbate%singleEnergy(i)
    tempstructures(j)%formationenergy = 0._dp

    ! - references
    allocate(tempstructures(j)%reference(CE%nreferences))  ! also: the references for the help structures
    do iref=1,CE%nreferences                                  !   do not really exist, so also mark them as
      tempstructures(j)%reference(iref)%isHelpStructure = .true.   !   help structure
      tempstructures(j)%reference(iref)%automaticReference = .false.
    enddo
  enddo
  ! - and their place in the list
  adsList = (/ (i, i=nStr+1,nStr+CE%CErank) /)

  ! calculate the formation energy of structure iStr of tempstructures (=structure(iStr)!)
  forall(iStr=1:nStr)
    structures(iStr)%formationenergy = calc_structure_DHf(tempstructures,iStr, &
                                                          pureEl=CE%pureList, isCE=.false.,adsEl=adsList)
  end forall

else
  stop "ERROR: Formation enthalpies only for non-coupled CEs or coupled adsorbate system CE"
endif



! Tobias Kerscher ???
!!lastgst = fitstrucnumb     !! already sorted by energy => last one is lowest in energy
!!iStr = 1                      !! at given concentration
!!
!!do while (abs(structures(iStr)%conc-structures(1)%conc)<epsilon(structures(iStr)%conc))
!!   iStr = iStr + 1
!!   if (iStr > nStr) exit
!!end do
!!
!!iStr = iStr - 1
!!firstgst = iStr
!!bulkenergya = structures(firstgst)%energy
!!bulkenergyb = structures(lastgst)%energy
!!!print *, "firstgst: ", firstgst, bulkenergya
!!!print *, "lastgst: ", lastgst, bulkenergyb
!!conca = structures(firstgst)%conc
!!concb = structures(lastgst)%conc


founda = .false.
foundb = .false.
foundc = .false.


!!GH (not useful) print *, " "
!!GH (not useful) if (CE%leftoutStr == 0) then
!!GH (not useful)    print *, "all structures will be included in the CE fit!"
!!GH (not useful) else
!!GH (not useful)    print *, "structures, which are not included in actual CE-fit"
!!GH (not useful)    print *, "Structure Title   /  concentration  /  energy      /    formation enthalpy "
!!GH (not useful)    write (*, '(A31, A8, A13, A12)') "structure title |","x(1) |","Edft/atom |","DHf/atom"
!!GH (not useful)    write (*, *) "----------------------------------------------------------------------"
!!GH (not useful)    do iStr = fitstrucnumb+1, nStr
!!GH (not useful)       write (*, '(A30, F8.3, F13.6, F12.5)') adjustl(structures(iStr)%title), structures(iStr)%conc, &
!!GH (not useful)            structures(iStr)%energy, structures(iStr)%formationenergy
!!GH (not useful)    end do
!!GH (not useful) end if
!!GH (not useful) print *, " "






!! Get groundstate line from DFT-values
print *, "getting DFT ground state line..."

if (nStr<2) stop "ERROR: need at least 2 input structures"

!OLD GSL>allocate(dftgstlist(fitstrucnumb))
!OLD GSL>nrdftgst = 1
!OLD GSL>dftgstlist = 0
!OLD GSL>dftgstlist(1) = firstgst


!!GH --not useful to leave out some structures-- fitstrucnumb = nStr - CE%leftoutStr
!!GHcall get_GSs(structures(:nStr-CE%leftoutStr),CE%CErank,CE%ncouplings,CE%xeps,CE%GSlist)
call get_GSs(structures,CE%CErank,CE%ncouplings,CE%xeps,CE%GSlist)

! the following assigments are for legacy reasons:
nrdftgst = size(CE%GSlist)
allocate(dftgstlist(nrdftgst))
dftgstlist = CE%GSlist


print *, ""
print *, "ground state line (DFT): (see also gsl_dft.out)"
print *

filename  = "gsl_dft.out"
header(1) = "# DFT ground state line"
call write_GSList(filename,structures(dftgstlist),header,CE,stdout=.true.)

allocate (CE%dftgstlist(nrdftgst))
do i=1, nrdftgst
   CE%dftgstlist(i) = dftgstlist(i)
end do

END SUBROUTINE prepare_input_structures




!***********************************************************************************
!sascha>! Chebychev section
!sascha>!
!sascha>! Function that calculates the chebychev polynome rank(i) of a given spin sigma
!sascha>! Just stole this from my old ternary code without further update
!sascha>FUNCTION chebychev_function(rank, sigma, i)
!sascha>implicit none
!sascha>integer,     intent(IN) :: rank ! rank of CE
!sascha>integer,     intent(IN) :: i
!sascha>integer(si), intent(IN) :: sigma
!sascha>real(dp) :: chebychev_function
!sascha>
!sascha>!if (useChebychevs) then
!sascha>
!sascha>Select case (rank)
!sascha>!-------------------------------------------------------
!sascha>case(2)   ! binary
!sascha>!-------------------------------------------------------
!sascha>  chebychev_function = real(sigma,dp)
!sascha>
!sascha>!-------------------------------------------------------
!sascha>case(3)   ! ternary
!sascha>!-------------------------------------------------------
!sascha>  Select case (i)
!sascha>  case(1)
!sascha>     chebychev_function = sqrt(3./2.)*real(sigma,dp)
!sascha>  case(2)
!sascha>     chebychev_function = sqrt(2.)*(1-(3./2.)*real(sigma,dp)**2)       
!sascha>  case default
!sascha>     print *, "wrong Chebychev-Polynome specified in input - check s-matrix: ", i
!sascha>     STOP
!sascha>  end Select
!sascha>
!sascha>!-------------------------------------------------------
!sascha>case default
!sascha>!-------------------------------------------------------
!sascha>  write(*,'(A)') "ERROR: the code was not compiled for a rank ", rank, "CE!"
!sascha>  stop
!sascha>
!sascha>end select
!sascha>
!sascha>!!else
!sascha>!
!sascha>!  if (i/=1) stop "ERROR: The code was compiled for binaries only."
!sascha>!  print *
!sascha>!  write(*,'(A)') "WARNING: You are using the Id function for the spins."
!sascha>!  print *
!sascha>!  chebychev_function = real(sigma,dp)
!sascha>
!sascha>!endif
!sascha>
!sascha>END FUNCTION chebychev_function

!sascha_deprecated>Select case (i)
!sascha_deprecated>case(1)
!sascha_deprecated>   chebychev_function = real(sigma,dp)
!sascha_deprecated>case(2)
!sascha_deprecated>   chebychev_function = 2*real(sigma,dp)**2 - 1
!sascha_deprecated>case(3)
!sascha_deprecated>   chebychev_function = 4*real(sigma,dp)**3 - 3*real(sigma,dp) 
!sascha_deprecated>case default
!sascha_deprecated>   print *, "wrong Chebychev-Polynome specified in input - check s-matrix: ", i
!sascha_deprecated>   STOP
!sascha_deprecated>end Select


!I consider this safe for up to 4-comp. systems until somebody proves me otherwise :)
!-sascha Oct/2012
FUNCTION gram_polynomial(sigma, i, rank)
implicit none
integer(si), intent(IN) :: sigma
integer,  intent(IN) :: i
integer, intent(IN) :: rank
real(dp) :: gram_polynomial, gram_polynomial_n_1, gram_polynomial_n_2 !recursion variables
real(dp) :: alpha_n_1, alpha_n_2
real(dp) :: x
integer :: n

!map spins to gramian points, to enforce spacing over [-1,+1]
x= map_spin_to_gramian_point(sigma, rank)

!reasoning behind this remap: larger spins managed to blow up the magnitude of the old chebychevs 
!starting at 4-naries (with spins of +2/-2). Gram-Polynomials in the smaller spin-range remedy this.
!-sascha Oct/2012 

gram_polynomial_n_2= 0 
gram_polynomial_n_1= 1
alpha_n_2 = 1

print*, "recursing Gram-polynomial of order ", i ,"in rank", rank, "!"

do n=1,i !starting a recursion (see Journal of Approximation Theory, 94, 128-143 (1998)
   alpha_n_1 = (real(rank,dp)/n) * sqrt( (n**2 - 0.25) / (rank**2-n**2) ) 
   gram_polynomial=2.0*x*gram_polynomial_n_1-alpha_n_1/alpha_n_2*gram_polynomial_n_2
   gram_polynomial_n_2=gram_polynomial_n_1 !it's a three-way recursion, we need to save the TWO last iterants
   gram_polynomial_n_1=gram_polynomial
enddo

END FUNCTION gram_polynomial


!subroutine to map spins from -rank/2,+rank/2 into the interval [-1,+1]
!-sascha Oct/2012
FUNCTION map_spin_to_gramian_point(spin,rank)
integer(si), intent(in) :: spin
integer, intent(in) :: rank
real(dp) :: map_spin_to_gramian_point

map_spin_to_gramian_point= -1.0 + (2.0*spin - 1.0)/rank

END FUNCTION map_spin_to_gramian_point

!subroutine to setup the lookup table
!I haven't changed a whole lot here; old spins remain. The remap into the Gramian interval is done internally.
!For all practical purposes, the casual user would still think uncle operates with spins between -rank/2,+rank/2
!Because I am a lazy bastard I sneaked in the counter-variable, so I don't have to reverse engineer the spin-counter from the spin.
!-sascha Oct/2012
subroutine setup_chebychev(rank,MPIr)
integer, intent(in) :: rank ! ce rank
integer, intent(in) :: MPIr
integer(si) iSpin
integer     iS
integer(si) :: counter
allocate(chebychev_lookup(-rank/2:+rank/2,1:rank-1))

if (MPIr==0) then
  print *
  write(*,'(A,I2,A)') "Setting up basis functions for a rank ",rank," cluster expansion:"
endif

counter = 1
do iSpin=-rank/2,+rank/2
if (mod(rank,2)==0 .and. iSpin==0) cycle   ! in a 2, 4, 6-nary CE, skip the spin==0
do iS=1,rank-1
  !sascha: I messed around here, just in case some1 notices something weird: that's prolly my fault.
  chebychev_lookup(iSpin,iS) = gram_polynomial(counter,iS, rank) 
  if (MPIr==0) write(*,'(5x,A,I2,A,I2,A,F10.6)') "chebychev(",iSpin,",",iS,")=",chebychev_lookup(iSpin,iS)
enddo
counter = counter +1

enddo
print *

end subroutine setup_chebychev

! Function to look into the lookup table
! (the use of this function is to make chebychev_lookup "universally" visible)
pure function chebychev(sigma,i)
integer, intent(in) :: i
integer(si), intent(in) :: sigma
real(dp) :: chebychev

chebychev = chebychev_lookup(sigma,i)
!It's still called chebychev, because I didn't want to change too much of the layout.
!Anyways, the new basis set works the same way. -sascha

end function chebychev
! END CHEBYCHEV SECTION
!***********************************************************************************


!***********************************************************************************
! Subroutine to mirror a point at a plane
!
! --> IN:
!              point  : 3d point that will be mirrored
!         unitnormal  : the UNIT normal of the mirror plane
!           distance  : the distance of the mirror plane from the origin
! <-- OUT:     
!        mirrorpoint  : the mirrored point
!
subroutine mirror_at_plane(point,unitnormal,distance,mirrorpoint)
real(dp), intent(in) :: point(3)
real(dp), intent(in) :: unitnormal(3)
real(dp), intent(in) :: distance
real(dp), intent(out):: mirrorpoint(3)

real(dp) :: dist2MirrorPlane

! calculate the distance of point to the mirror plane:                                                           
! mind: distance(P) < 0 if (0,0,0) and point P are on the same side of the plane                                           
!       distance(P) > 0 if (0,0,0) and point P are on different sides of the plane                                         
dist2MirrorPlane=dot_product(point,unitnormal)-distance

! mirror the point on the plane:                                                                                     
mirrorpoint= point-2*dist2MirrorPlane*unitnormal

end subroutine mirror_at_plane


!********************************************************************************
! Function that returns the dset-label of a given point
!
! --> IN:
!           pLV  : parent lattice vectors
!         point  : the point we want to know the dset-label of
!          dset  : all the dset points
!           eps  : finite precision check
! <-- OUT:
!        the function returns the dset-label of "point". If the point
!        is NOT part of the lattice (and hence the dset), the value is
!        -1. Take care that your program stops then.
integer function get_dset_label(pLV, point, dset, eps)
real(dp), intent(in) :: pLV(:,:)
real(dp), intent(in) :: point(:)
real(dp), intent(in) :: dset(:,:)
real(dp), intent(in) :: eps

logical :: checkPoint
integer :: id, nd

nd = size(dset,2)
checkPoint = .false.

do id = 1, nd     ! Try shifting by each parent basis vector to take care of the                                            
                  ! general case where the parent cell is a multilattice                                                    
  if (is_lattice_point(pLV,point-dset(:,id),eps)) then
    checkPoint = .true.
    exit
  endif
enddo

if (checkPoint) then
  get_dset_label=id
else
  get_dset_label=-1
endif

end function get_dset_label




!OLD>!********************************************************************************
!OLD>! Setup the site labels for the find_site_equivalencies routine
!OLD>! IN -->
!OLD>!             CE  : CE data, in particular contains
!OLD>!                   - nBas, number of basis atoms (dset size)
!OLD>!                   - surface, surface variables
!OLD>!                     - isSurf, logical switch for surface CE
!OLD>!                     - BasType, type of the basis atoms (dset) "B/S/A"
!OLD>!                     - BulkSpin, spin (i.e. occupation) of the bulk-like atoms
!OLD>! <-- OUT
!OLD>!          label  : site labels
!OLD>!                   Mind: for a surface CE the bulk-like atoms get the fixed occupation (spin) as label
!OLD>!                         for the surface atoms label=10 (arbitrary)
!OLD>!
!OLD>subroutine get_site_labels(label,CE)
!OLD>integer, pointer :: label(:) ! intent: out
!OLD>type(CE_variables), intent(in) :: CE
!OLD>
!OLD>allocate(label(CE%nBas))
!OLD>
!OLD>if (.not. CE%surface%isSurf) then
!OLD>  !____________________
!OLD>  ! Bulk CE
!OLD>  label = 1  ! Spacegroup code needs labels but in our case we only have points,
!OLD>             ! no atoms, so the labels are all the same.
!OLD>else
!OLD>  !____________________
!OLD>  ! Surface CE
!OLD>  where(CE%surface%BasType=='B'); label=CE%surface%BulkSpin; end where
!OLD>  where(CE%surface%BasType=='S'); label=10; end where
!OLD>  if (any(CE%surface%BulkSpin==10)) stop "ERROR: Bad programming. Bulk Occ = 10, whilst label=10 is reserved for surface atoms"
!OLD>  if (any(CE%surface%BasType=='A')) stop "ERROR: Adsorbates not yet implemented"
!OLD>endif
!OLD>end subroutine get_site_labels



!********************************************************************************
! Purpose:
! Takes a structures and constructs a reference structure according to the rules in lat.in
!
! IN:   str    : structure
!       CE     : information about the CE
!      iRef    : the index of the reference
! OUT:  rstr   : reference structure
!
subroutine setup_reference_structure(str,rstr,CE,iRef)
type(crystal), intent(in)      :: str
type(crystal), intent(inout)   :: rstr
type(CE_variables), intent(in) :: CE
integer, intent(in)            :: iRef
character(2) cRef
integer iAt,cAt,j,k

!--------------------------------------------------------------------------------
! Allocate
rstr%nAtoms = size(str%pos,2) - CE%l2r(iRef)%nDel * str%nAtoms * 1.0_dp / CE%nBas

allocate(rstr%pos(3,rstr%nAtoms))
allocate(rstr%spin(rstr%nAtoms))
allocate(rstr%pBVlab(rstr%nAtoms))
allocate(rstr%intPt(3,size(CE%reference(iRef)%pBas,2)))
allocate(rstr%aTyp(rstr%nAtoms))

!--------------------------------------------------------------------------------
! Name structure with
!    Ref(iRef)_str%title
write(cRef,'(I2.2)') iRef
rstr%title = "Ref("//cRef//")_"//str%title

!--------------------------------------------------------------------------------
! Bring the two lattices together
!
! {vec r'} = M . {vec r}
! However, we have to implement it the other way round, since in our input the columns and
! rows of both LV and Mlv are exchanged. So we want to do
! {vec r'} = M^T . {vec r}^T = {vec r}.M;
rstr%LV= matmul(str%LV,CE%l2r(iRef)%Mlv)
!--------------------------------------------------------------------------------
! Bring the atoms to the reference lattice
!
! (1)
! vec r' = M. vec r + O  for all atom positions r.
! Same problem here, have to use M^T in the code.
! (2)
! Remove the atoms that are specified in dsetDel (because they are
! the same as other ones when we construct a reference struc with the above
! lattice vectors)
! (3)
! Set the dset-label right (pBVlab), and the atomic type (aTyp)
cAt = 0

do iAt=1,str%nAtoms
  if (any(CE%l2r(iRef)%dsetDel==str%pBVlab(iAt))) cycle
  cAt = cAt+1
  rstr%pos(:,cAt)  = matmul(transpose(CE%l2r(iRef)%Mds),str%pos(:,iAt)) + CE%l2r(iRef)%Ods(:)
  rstr%spin(cAt)   = str%spin(iAt)
  rstr%pBVlab(cAt) = get_dset_label(CE%reference(iRef)%pLV, rstr%pos(:,cAt), CE%reference(iRef)%pBas, CE%eps)
  rstr%aTyp(cAt)   = str%aTyp(iAt)
enddo


!--------------------------------------------------------------------------------
! Other properties to set
rstr%pLV   = CE%reference(iRef)%pLV
rstr%intPt = CE%reference(iRef)%pBas
rstr%eps   = str%eps
rstr%automaticReference = .true.

end subroutine setup_reference_structure


!********************************************************************************
! Purpose:
! Given a set of structures (containing their concentrations x(:,:) and formation energies)
! and the dimension of the CE (multinary), return the ground states.
!
! IN -->
!              str  : structures
!                k  : rank of the CE          |  the dimension of the E-x-diagram will be
!               nC  : number of CE couplings  |  calculated out of these: dim = (k-1)*nC + 1
!
! <-- OUT
!              list : the indices of the GSs
!              dist : (optional, only together with adjlist) 
!                     for each str: the distance to the convex hull hyperplane
!           adjlist : (optional, only together with dist)
!                     for each str: the list of adjacent GSs, i.e. what GSs span the 
!                     hyperplane above which the str is floating?
!
subroutine get_GSs(str,k,nC,eps,list,dist_,adjlist_)
type(crystal), intent(in)   :: str(:)
integer, intent(in)         :: k  ! rank
integer, intent(in)         :: nC  ! couplings
integer                     :: dim  ! dimension of the E-x diagram
real(dp), intent(in)        :: eps  ! concentration check eps
integer, pointer            :: list(:) ! out
integer, pointer            :: listSortIdx(:)
real(dp), pointer, optional :: dist_(:)
integer, pointer, optional  :: adjlist_(:,:)

real(dp), pointer :: dist(:)
integer, pointer  :: adjlist(:,:)

logical               :: getDist, getAdjList

integer               :: iStr, jStr, nStr, iExPoint, nExPoint, idim, c, nlist, iPlane, nPlane, iV, minPlane(1), iC, iK
real(dp), allocatable :: ExPoint(:), d(:), Edirection(:), r(:,:)
real(dp), pointer     :: N(:,:), O(:)
integer, allocatable  :: IP(:,:)



integer nPositiveDHf
integer istat
logical err
integer rDim, rnP, rnN, rnIP, rnV
character(1) :: rC
integer, allocatable :: rV(:), rIP(:,:)
real(dp), allocatable         :: rP(:,:)
real(dp), allocatable, target :: rN(:,:), rO(:)

character(80) :: line

integer i

if (present(dist_   )) then; getDist=.true.   ; else; getDist=.false.   ; endif
if (present(adjlist_)) then; getAdjList=.true.; else; getAdjList=.false.; endif
if ((getDist .and. .not. getAdjList) .or. (getAdjList .and. .not. getDist)) stop "ERROR: bad programming in get_GSs."

! Determine the dimension of the E-x-diagram:
! For each coupling we need only k-1 dimensions (cf. simple binary expansion: rank=2, but in fact 1 concentration
! suffices; ternary: rank=3, but 2 concentration suffice...), and finally 1 additional dimension for the energy
! coordinate
dim = (k-1)*nC + 1

nStr=size(str)
nExPoint = nStr*dim
allocate(ExPoint(nExPoint))
allocate(Edirection(dim))
!write(*,'(A)') "Convex hull routine:"
!write(*,'(A)') ". setup"

!--------------------------------------------------------------------------------
! Setup the points for finding the convex hull
! A point has the following coordinates:
!   (Dhf, x(1), x(2), ..., x(k-1))
!    ~E  ~x ==> "Ex"-Point
! i.e. 
! Edirection = (1,0,0,...)
!
Edirection = 0; Edirection(1) = 1; ! do not need it here but only for determining the distance
                                   ! of an ExPoint to the convex hyperplane in the E direction.
                                   ! (Have put it here in order not make the direction clear)
c=0
do iStr=1,nStr
  c=c+1
  ExPoint(c)=str(iStr)%formationenergy ! fill ExPoint with the E (energy)
  do iC=1,nC ! go through all couplings
    do iK=1,k-1 ! go through all ranks (can skip the last, since k-1 suffice to determine all concentrations)
      c=c+1
      ExPoint(c) = str(iStr)%x(iK,iC) ! fill ExPoint with the x (concentrations)
    enddo
  enddo
enddo

if (all(str(:)%formationenergy==0.0)) then 
  write(*,'(A)') ". WARNING: all structures have DHf=0."
  write(*,'(A)') ". WARNING: skipping convex hull routine: all structures are ground states."
  allocate(list(nStr))
  list= (/ (iStr,iStr=1,nStr) /)
  return
endif

!--------------------------------------------------------------------------------
! Get convex hull 
! => convex_hull.out
!write(*,'(A)') ". call routine"
call qh_getConvexHull(ExPoint,nStr,dim)

!--------------------------------------------------------------------------------
! Read data from convex_hull.out
! Format of the file as of revision 4:
!
!    # Vertices section
!    <number of vertices>
!       <vertex id 1>
!       ...
!    <dimension>
!    <number of points = number of vertices>
!       <point 1, i.e. concentrations and energy coordinates>
!       ...
!    # Normals and Offsets section (Hyperplanes)
!    <dimension+1>
!    <number of normals = number of offsets = number of hyperplanes of the convex hull>
!        <normal 1> <offset 1> 
!        ...
!    # Incident points
!    <number of incidence lists = number of normals>
!        <list of incident points of plane 1>
!        ...
!
! MIND: the id of the vertices is off by 1
!
!write(*,'(A)') ". read data"
open(100,file="convex_hull.out",iostat=istat)
if (istat/=0) stop "ERROR: finding convex hull or opening file convex_hull.out"
   read(100,*) rC                                                   ! comment
   read(100,*) rnV              ; allocate(rV(rnV))                 ! number of vertices
do i=1, rnV            
   read(100,*) rV(i)                                                ! vertices
enddo
   read(100,*) rC                                                   ! comment
   read(100,*) rDim             ; if (rDim /= dim) stop "ERROR: convex hull"  ! dimension
   read(100,*) rnP              ; allocate(rP(dim,rnP));            ! number of points
do i=1,rnP
   read(100,*) rP(:,i)                                              ! points (reals!!!)
end do
   read(100,*) rC                                                   ! comment
   read(100,*) rDim             ; if (rDim /= dim+1) stop "ERROR: convex hull" ! dimension+1
   read(100,*) rnN              ; allocate(rN(dim,rnN), rO(rnN))    ! number of normals and offsets
do i=1,rnN
   read(100,*) rN(:,i), rO(i)                                       ! normals and offsets
enddo
   read(100,*) rC                                                   ! comment
   read(100,*) rnIP             ; allocate(rIP(dim,rnIP))           ! number of incident points
do i=1,rnIP
   read(100,*) rIP(:,i)                                             ! incident points (integer)
enddo
close(100)

!--------------------------------------------------------------------------------
! Get the struc numbers of the GSL

! setup list
allocate(list(rnV))

! get the id right (the qhull algorithm always starts counting with 0)
list = rV + 1

! ditch positive formation energies
nPositiveDHf = count(str(list)%formationenergy>0.0_dp)
where(str(list)%formationenergy>0.0_dp)
  list = -1
end where
list = pack(list,list>0)
nlist = size(list)-nPositiveDHf
list => ralloc(list,nlist)

! ditch vertices that have a struture beneath them (can happen in particular with 
! scans over the configuration space that are not complete)
do iStr=1,nlist    ! go through all the structures in "list" = the ground states
   do jStr=1,nStr  ! ... and compare them to all other structures.
      if ( equal(str(list(iStr))%x-str(jStr)%x,0._dp,eps) ) then    ! if the concentration matches
         if ( str(jStr)%formationenergy < str(list(iStr))%formationenergy ) then   ! and the other structure has smaller formation energy
            list(iStr)=-1 ! mark the current list entry as non-valid ground state.
            exit
         endif
      endif
   enddo
enddo
nlist = size(list)-count(list<0)
list = pack(list,list>0)
list => ralloc(list,nlist)

!--------------------------------------------------------------------------------
! Get distances to GS hyperplane and adjacent vertices (i.e. for each structure:
! the vertices that span the hyperplane above which the structure floats). 
! Originally, I planned to do the following ONLY if we really want to return
! the distance list from that subroutine, i.e. if (getDist).
! However, there are cases when ExPoints really lie right on a hyperplane, that
! qhull does not return them as vertices of the hyperplane. In order to remedy
! this we determine the distances anyway...
!write(*,'(A)') ". postprocessing"
nPlane = rnN
allocate(dist(nStr))
allocate(d(nPlane))
allocate(adjlist(dim,nStr))
allocate(IP(dim,rnIP))

! In qhull, the hyperplanes are defined by
!    N.r + r0 = 0
! N : normal
! r0: offset
!
N => rN       ! hyperplane normals
O => rO       ! hyperplane offsets
IP = rIP+1    ! incident points (index is off by 1)

allocate(r(dim,nstr))
r = reshape(ExPoint(:),(/dim,nStr/))  ! bring the ExPoints in a dim-dimensional vector format

do iStr=1,nStr
  if (any(list==iStr)) then
    dist(iStr)=0._dp      ! this is a ground state. No further calc needed. The =0 setting
                          ! eschews possible numerical fluctuations.
    adjlist(:,iStr)=iStr
    cycle                 ! go on with next struc
  endif

  ! if the structure is not a ground state: we continue...

  forall(iPlane=1:nPlane)  ! go through all planes and offsets and determine distance
  ! to plane in the energy direction. If the plane is BELOW the 
  ! structure, d > 0. The (unphysical) planes ABOVE the structure
  ! (i.e. the convex hull parts that are not part of the GSL) have 
  ! d < 0.
    d(iPlane) = ( dot_product(N(:,iPlane),r(:,iStr))+O(iPlane) ) / ( dot_product(N(:,iPlane),Edirection(:))  )
  end forall

! The following if-branching and its explanation is not totally true for coupled systems.
!!    if (str(iStr)%formationenergy <= 0._dp) then
  dist(iStr) = minval(d(:),d>0)   ! the minimal distance to all hyperplanes BELOW (in energy direction)
                                  ! the structure is the distance to the convex hull. (I don't want
                                  ! to check for d>=0 here, since there might be ExPoints lying on the
                                  ! unphysical upper hyperplanes (with DHf > 0). Those points will
                                  ! find a physical hyperplane below, if d>0.)
  minPlane = minloc(d(:),d>0)
!!    else
!!      dist(iStr) = minval(d(:),d>=str(iStr)%formationenergy)  ! this construction is a failsafe for 
!!                                      ! structures for positive formation energies. Their distance
!!                                      ! to the GS hyperplane must be larger than their formation energy.
!!                                      ! (It can happen that a structure lies above the (unphysical) "upper"
!!                                      ! convex hull, at least if precision is not high enough. Then the distance
!!                                      ! would be measured to that upper hull)
!!      minPlane = minloc(d(:),d>=str(iStr)%formationenergy) 
!!    endif

  ! Did we find the distance?
  if (minPlane(1)==0) then
      ! no, didn't find a distance. This is most probable the case if the structure is NOT a ground state
      ! in the qhull sense, i.e. it is no vertex of the convex hull, but -- nevertheless -- lies right
      ! on a hyperplane, i.e. for us it IS a ground state.
      ! Admitted, this scenario will probably only occur with perfect input data, but we want to be 
      ! prepared.

    dist(iStr) = minval(d(:),d>=0)   ! the minimal distance to all hyperplanes "BELOW" (in energy direction)
                                     ! the structure is the distance to the convex hull. This time, 
                                     ! we check for d>=0 instead of d>0, i.e. we don't really check if
                                     ! we are ABOVE the hyperplane. If you compare to some lines above 
                                     ! (where we check for d>0), the situation now is different: ExPoints
                                     ! that lie on an unphysical upper hyperplane are already accounted
                                     ! for above, because those necessarily find a lower hyperplane. 
                                     ! Here, we want also to account for those d=0 that sit on a lower
                                     ! hyperplane.
    minPlane = minloc(d(:),d>=0)
      
    if (minPlane(1)==0) then ! now we really should have found something...
        write(*,'(A)') "ERROR: in get_GSs"
        write(*,'(A)') "ERROR:   could not determine distance to convex hull for"
        write(*,'(A,I4)         ') "ERROR:   structure #", iStr
        write(*,'(A,A)          ') "ERROR:      title  : ", adjustl(str(iStr)%title)
        write(*,'(A,10(F8.3,1x))') "ERROR:          x  :", str(iStr)%x
        stop
    endif


    ! if we found a structure right on the hyperplane, we have also to add it to the GS list
    nlist = nlist+1 ! we add one element to the list
    list => ralloc(list,nlist) ! memory shifting
    list(nlist) = iStr ! add the element: iStr      
  endif

  ! Now set adjacent vertices:
  if (existsVertexAtX(str(iStr)%x,iV)) then  ! is there a vertex (GS) at this conc?
    adjList(:,iStr)=list(iV)                 ! "adjacent" vertices are only this vertex
  else
    adjList(:,iStr)=IP(:,minPlane(1))   ! no vertex below that structure's conc. Get info about
                                        ! adjacent vertices from the IP (former read-in rIP)
  endif

enddo ! iStr=1,nStr


!--------------------------------------------------------------------------------
! sort list
allocate(listSortIdx(size(list)))
listSortIdx=(/(i,i=1,size(list))/)
call heapsort(listSortIdx,list)


!--------------------------------------------------------------------------------
! assign optional return values
if (getDist)    then; dist_     => dist;    else; deallocate(dist);    endif
if (getAdjList) then; adjlist_  => adjList; else; deallocate(adjList); endif

!--------------------------------------------------------------------------------
! deallocations
deallocate(rP,rN,rO,rIP,rV)
deallocate(listSortIdx)
deallocate(ExPoint,Edirection)
deallocate(d, r, IP)

contains

logical function existsVertexAtX(x,iV)
integer  :: iV, nV
real(dp) :: x(:,:)
nV=nlist

existsVertexAtX = .false.
do iV=1,nV
  if (equal(str(list(iV))%x,x,eps)) then
    existsVertexAtX = .true.
    exit
  endif
enddo
end function existsVertexAtX

end subroutine get_GSs


!********************************************************************************
! Purpose:
! Given a list of structures it returns a list (with an entry for every structure)
! that tells us what structure is the lowest in energy at that concentration.
! Furthermore, it measures the distance to this structure (for every str)
!
! --> IN:
!              str : structures
! <-- OUT:
!             list : for each structure: the index of the structure that is
!                    lowest in energy at that concentration
!             dist : (optional) for each structure: the distance of that
!                    structure to the structure that is lowest in energy at
!                    that concentration
!           
subroutine get_lowEstrucs(str, eps, list, dist, useFitEnergy)
type(crystal), intent(IN) :: str(:)
real(dp), intent(IN)      :: eps ! concentration check epsilon
integer, pointer          :: list(:) ! intent(out)
real(dp), pointer, optional :: dist(:)
logical, optional           :: useFitEnergy

real(dp), pointer :: d(:) ! internal distance measurement
real(dp)          :: maxd(1)

logical, pointer  :: sameConcStr(:)

logical :: getDist, FitEnergyDist

integer iStr, nStr, jStr


if (present(dist)) then; getDist=.true.; else; getDist=.false.; endif
if (present(useFitEnergy)) then; FitEnergyDist=useFitEnergy; else; FitEnergyDist=.false.; endif

nStr = size(str)

allocate(list(nStr))
if (getDist) then
  allocate(dist(nStr))
  dist=0._dp
endif

allocate(sameConcStr(nStr))
allocate(d(nStr))
do iStr=1,nStr

  ! (1) check what structures have the same concentration as iStr
  do jStr=1,nStr
    if (equal(str(iStr)%x, str(jStr)%x, eps)) then
      sameConcStr(jStr) = .true.
    else
      sameConcStr(jStr) = .false.
    endif
  enddo

  ! (2) for all structures that have the same concentration as iStr:
  !     calculate the difference in (fit)energy to iStr
  !     (tk: I'd like to replace it by a where construct, but confess that I'm not 100% sure
  !     whether something like
  !        where (sameConcStr); d = str(iStr)%energy - str(:)%energy; endwhere
  !     produces the correct result (because of the str(:)) )
  do jStr=1,nStr
    if (sameConcStr(jStr)) then
      if (.not. FitEnergyDist) then; d(jStr) = str(iStr)%energy    - str(jStr)%energy
                               else; d(jStr) = str(iStr)%fitenergy - str(jStr)%fitenergy
                               endif;
    endif
  enddo

  ! (3) for all structures that have the same concentration as iStr:
  !     determine the maximal (positive!) distance => that's the lowest structure at that concentration
  !  a) set the lowest structure into list
  maxd = maxloc(d, (d>=0._dp) .and. (sameConcStr .eqv. .true.) )
  list(iStr) = maxd(1)
  !  b) if desired: set the intent(out) distance
  if (getDist) dist(iStr) = maxval(d, (d>=0._dp) .and. (sameConcStr .eqv. .true.) )

enddo


! old code, that seems to rely on the pre-sorting by energy of the structures
!!do iStr=1,nStr  ! this is probably not the best way to do it, but we only do it once.
!!
!!  list(iStr)=iStr
!!
!!  do jStr=1,nStr
!!    if (jStr==iStr) cycle
!!    if (.not. equal(str(iStr)%x,str(jStr)%x,eps)) cycle
!!    if (str(jStr)%energy < str(iStr)%energy) then   ! we are at 1 concentration 
!!                                                    ! => formationenergy-difference == energy-difference
!!      list(iStr)=jStr
!!      if (getDist) then
!!        if (.not. FitEnergyDist) then; dist(iStr) = str(iStr)%energy    - str(jStr)%energy
!!                                 else; dist(iStr) = str(iStr)%fitenergy - str(jStr)%fitenergy
!!        endif ! distance is measured as energy-difference or as fitenergy-difference
!!      endif ! retrieve distance?
!!    endif ! concentrations are equal
!!  enddo ! jStr
!!enddo ! iStr

end subroutine get_lowEstrucs



!********************************************************************************
! Routine to find the pure elements in a set of structures
!
! IN -->
!               C  : find the pure elements of the CE coupling C
!             str  : the structures to find the pure elements in
!               k  : CE rank
!              nC  : total number of CE couplings (so, C <= nC)
!           xcond  : (optional) condition matrix for the concentrations
!                    { 1,2,3 ; cond# }:  for each condition ( cond# ) supply
!                    the concentration index (1), the CE index (2) and the concentration (3)
!                    that you want to check. If xcond==0, then no condition is assumed.
! <-- OUT
!            list  : the list of pure elements of CE coupling C
!
subroutine get_pureElements(C,str,k,nC,list,xcond,xeps)
integer, intent(in)          :: C
type(crystal), intent(in)    :: str(:)
integer, intent(in)          :: k
integer, intent(in)          :: nC
integer, pointer             :: list(:)   ! intent(out)
!integer, pointer             :: tlist(:)
real(dp), intent(in),optional :: xcond(:,:) ! condition matrix: { value, cond# }
real(dp), intent(in), optional :: xeps

integer iStr, nStr, ix, nx, iC
integer nCond, iCond, nlist
integer, allocatable :: pureEl(:,:,:)
integer, allocatable :: npureEl(:,:)
real(dp), allocatable :: x(:,:,:)  ! { x#, CE#, struc# }

!if (nC>1) write(*,*) "WARNING: get_pureElements is for CE #1 ONLY!"

nStr=size(str)
nx=k

!allocate(x(k,nC,nStr))
!
!
!forall(iStr=1:nStr)
!  x(:,:,iStr) = str(iStr)%x(:,:) 
!end forall

allocate(pureEl(k,nC,nStr))
allocate(npureEl(k,nC))
npureEl = 0
pureEl  = 0
do iC=1,nC
  do ix=1,nx
    do iStr=1,nStr
      if (  str(iStr)%x(ix,iC) >= 1._dp  ) then
        npureEl(ix,iC) = npureEl(ix,iC)+1
        pureEl(ix,iC,npureEl(ix,iC)) = iStr
      endif
    enddo
  enddo
enddo

! this was nice, but not general enough:
!do iC=1,nC
!do ix=1,nx
!  pureEl(ix,iC,:) = maxloc(x(ix,iC,:),x(ix,iC,:)>=1._dp)
!end do
!end do
!


! conditions following:
if (present(xcond)) then
  if (.not. all(xcond==0)) then ! allow for the possibility that all conditions are 0. This should be equal
                                ! to NO condition.
    if (size(xcond,1)/=3) stop "ERROR: in get_pureElements. Dimension 1 of the xcond matrix must be 3"
  
    ncond = size(xcond,2)   ! number of supplied conditions
    do icond=1,ncond       ! check each condition separately
  
      do iC=1,nC
      do ix=1,nx
      nlist = npureEl(ix,iC)
      do iStr=1,nlist
        if (  .not. equal( str(pureEl(ix,iC,iStr))%x(nint(xcond(1,icond)),nint(xcond(2,icond)))  ,  xcond(3,icond),xeps) ) then
          npureEl(ix,iC) = npureEl(ix,iC) - 1
          pureEl(ix,iC,iStr) = 0
        endif
      enddo
      enddo
      enddo
  
    enddo
  endif
endif

! final check and finish
if (any(npureEl(:,C)>1))  stop "ERROR: in get_pureElements. More than 1 pure element found."
if (any(npureEl(:,C)==0)) stop "ERROR: in get_pureElements. No pure element for supplied CE found."

allocate(list(k))
list = pack(pureEl(:,C,:),pureEl(:,C,:)>0)

end subroutine get_pureElements


!********************************************************************************
! Calculate the formation enthalpy
!
! IN:
!         str   : a _set_ of structures
!                 important subvariables are:  %energy _or_ %ce_energy (see flag isCE)
!                                              %x: the concentration of the elements in str
!           i   : the index of the structure in str that we want to calculate DHf of
!      pureEl   : the indices of the structures in str that are "pure" elements.
!                 NOTE: The first entry corresponds to the first pure element, i.e. x(1) = 100%
!        isCE   : (optional) logical flag that determines whether to use %energy (=.false.) or
!                 %ce_energy (=.true.) for the calculation of the formation enthalpy. 
!                 If not given, .false. is assumed (i.e. %energy is used)
!       adsEl   : (optional) the indices of the structures in str that are adsorbate elements.
!                 NOTE: The first entry corresponds to the first adsorbate element
! 
! OUT:
!      calc_structure_DHf  : the function's value is the formation enthalpy (DHf) of str(i)
!
! IMPORTANT:
!
! - calculation of "normal" formation enthalpy uses 
!      str%energy, or str%ce_energy (resp.)
!      str%x(:,1)
! - calculation of adsorbate formation enthalpy uses
!      str%energy, or str%ce_energy (resp.)
!      str%nTyp(:,2)
!
pure function calc_structure_DHf(str,i,pureEl,isCE,adsEl)
real(dp) :: calc_structure_DHf

type(crystal), intent(in) :: str(:)
integer, intent(in) :: i
integer, intent(in) :: pureEl(:)
logical, intent(in), optional :: isCE
integer, intent(in), optional :: adsEl(:)

real(dp) :: DHf
logical :: useCEenergy, calc_adsorbate_DHf


if (present(isCE))  then; useCEenergy=isCE;          else; useCEenergy=.false.;        endif
if (present(adsEl)) then; calc_adsorbate_DHf=.true.; else; calc_adsorbate_DHf=.false.; endif

!--------------------------------------------------------------------------------
! "Normal" formation enthalpy
!
! NB: The pure elements are only defined for the 1st CE, so we can also
! put x(:,1) here
if (.not. useCEenergy) then
  DHf = str(i)%energy - sum(str(pureEl(:))%energy * str(i)%x(:,1) )
else
  DHf = str(i)%ce_energy - sum(str(pureEl(:))%ce_energy * str(i)%x(:,1) )
endif
!--------------------------------------------------------------------------------
! Adsorbate formation enthalpy
if (calc_adsorbate_DHf) then 
  DHf = DHf - sum(str(adsEl(:))%energy * str(i)%nTyp(:,2)) / sum(str(i)%nTyp)    ! all energies are "per site",
                                                             ! so the additional energy of the free adsorbate
                                                             ! must also be scaled by site
endif ! adsorbate DHf?
!--------------------------------------------------------------------------------

calc_structure_DHf = DHf 

end function calc_structure_DHf





!OLD>logical function structure_is_symmetric(pos, spin, surf, eps)
!OLD>real(dp)                , intent(in) :: pos(:,:)
!OLD>integer                 , intent(in) :: spin(:)
!OLD>type(surface_variables) , intent(in) :: surf
!OLD>real(dp)                , intent(in) :: eps
!OLD>
!OLD>real(dp), allocatable :: mpos(:,:) ! mirror pos
!OLD>integer, allocatable  :: mspin(:)  ! mirror spin
!OLD>integer iAt, jAt, nAt
!OLD>character(1)::dummy
!OLD>
!OLD>nAt = size(pos,2)
!OLD>
!OLD>allocate(mpos(3,nAt))
!OLD>allocate(mspin(nAt))
!OLD>
!OLD>
!OLD>! mirror the structure on the symmetry plane to check for symmetry
!OLD>do iAt=1,nAt
!OLD>  call mirror_at_plane( pos(:,iAt), &
!OLD>                        surf%unitNormal, surf%distSymPlane, &
!OLD>                        mpos(:,iAt) )
!OLD>  mspin(iAt)=spin(iAt)
!OLD>enddo
!OLD>
!OLD>! now: compare str and mirrorstr, should be exactly alike
!OLD>do iAt=1,nAt
!OLD>  do jAt=1,nAt
!OLD>    if (equal(mpos(:,iAt),pos(:,jAt),eps)) exit
!OLD>    if (jAt == nAt) then
!OLD>      deallocate(mpos, mspin)
!OLD>      structure_is_symmetric = .false.
!OLD>      return
!OLD>    endif
!OLD>  enddo
!OLD>  if (mspin(iAt) == spin(jAt)) cycle
!OLD>  deallocate(mpos, mspin)
!OLD>  structure_is_symmetric = .false.
!OLD>  return
!OLD>end do
!OLD>
!OLD>deallocate(mpos, mspin)
!OLD>structure_is_symmetric = .true.
!OLD>return
!OLD>
!OLD>end function structure_is_symmetric



subroutine set_intersection(set1,set2,intersection)
integer, intent(in) :: set1(:), set2(:)
integer, pointer    :: intersection(:) ! intent(out)

integer :: n1, n2, ni, iEl

n1 = size(set1)
n2 = size(set2)

ni = 0
do iEl=1, n1
  if (any(set2 == set1(iEl))) ni=ni+1
enddo

allocate(intersection(ni))

ni = 0
do iEl=1, n1
  if (any(set2 == set1(iEl))) then
    ni = ni+1
    intersection(ni) = set1(iEl)
  endif
enddo

end subroutine set_intersection





pure function map_labeling_position_2_dvector_index(nUC,nD)
! function return value:
integer ::  map_labeling_position_2_dvector_index(nUC*nD)
! input values:
integer, intent(in) :: nUC, nD   ! number of unit cells, number of dvectors per unit cell
integer i

! this array tells us the dset number to which a position in the labeling corresponds to,
! e.g. in a labeling xxyyzz for nUC=2 (total of 2 unit cells) and nD=3 (3 dset members), 
! all x belong to dset 1, y to 2, z to 3.
map_labeling_position_2_dvector_index = ( (/ (i,i=1,nUC*nD) /) + (nUC-1) ) / nUC ! integer division

end function map_labeling_position_2_dvector_index



subroutine write_GSList(file,strucs,header,CE,stdout)
character(80), intent(in)  :: file      ! filename
type(crystal), intent(in)  :: strucs(:) ! the GS structures (and ONLY THOSE)
character(80), intent(in)  :: header(:) ! header information of the file
type(CE_Variables), intent(in) :: CE
logical, optional          :: stdout

integer i, n
character(100) :: formatstring

open(100,file=file,status='REPLACE')
do i=1,size(header)
   write(100,'(A)') header(i)
end do

write(100,'(A)') "# x(1:k,1:nC) | DHf / atom [eV] | structure title"
if (present(stdout)) then
   if (stdout) write(  *,'(A)') "  x(1:k,1:nC) | DHf / atom [eV] | structure title"
endif

formatstring="("; write(formatstring(2:4),'(I2)') CE%CErank*CE%ncouplings; 
formatstring = trim(formatstring)//"(F15.10,1x),5x,F16.10,5x,A80)"

write (100, '(A)') "# -------------------------------------------------------------------------------"
if (present(stdout)) then
   if (stdout) write (*, *) "-------------------------------------------------------------------------------"
endif

do i=1, size(strucs)
   write (100,formatstring) strucs(i)%x(:,:), strucs(i)%formationenergy, adjustl(strucs(i)%title)
   if (present(stdout)) then
      if (stdout) write (  *,formatstring) strucs(i)%x(:,:), strucs(i)%formationenergy, adjustl(strucs(i)%title)
   endif
enddo

if (present(stdout)) then
   if (stdout) print *, " "
endif

close(100)

end subroutine write_GSList


! Is a dlabel deleted in a reference ce?
pure function is_dvector_removed_in_reference_ce(dlabel,CEorig,iRef)
logical is_dvector_removed_in_reference_ce
integer, intent(in) :: dlabel            ! the label of the dvector to check
type(CE_variables), intent(in) :: CEorig ! original CE data (not reference!)
integer, intent(in) :: iRef              ! the index of the reference

if ( any(CEorig%l2r(iRef)%dsetDel==dlabel) ) then
  is_dvector_removed_in_reference_ce = .true.
else
  is_dvector_removed_in_reference_ce = .false.
endif
end function is_dvector_removed_in_reference_ce


ENDMODULE structure_utilities
