module io_utilities
use num_types
use ce_types
use vector_matrix_utilities
use numerical_utilities
use structure_utilities, only: get_spinRank_mapping, mirror_at_plane, setup_reference_structure, calc_structure_DHf, &
                               map_labeling_position_2_dvector_index, get_ilabeling
use utilities_module, only: ucase
use sort, only: sort_figure_vertices, heapsort
use symmetry_module, only: bring_into_cell
use enumeration_types, only: maxLabLength
implicit none

integer :: dummy

private ran
public read_lattdef, write_clusters, newunit, &
     write_figure_reps, writeCEresults, co_ca, read_input_structures, &
     read_CEfitting, initializeRND, rnd, initialize_random_number, check_enum_file, read_MCpar, read_MDpar, &
     read_structpred, perform_automaticsearch, &
     read_gss, read_MCevaluation, read_MCcell, write_MCcell, goto_enumStart, skip_lines, concatenate_files, &
     read_adsorbateData, write_multilattice_dvector_occupation, write_InpStrucList, write_MC_correlations, read_CSvars, write_structures_in, read_toy_js_and_noise !, ucase

interface rnd
  module procedure rnd_GA
  module procedure rnd_seed
end interface

interface cmp_arrays
  module procedure cmp_arrays1
  module procedure cmp_arrays2
end interface

CONTAINS
!<summary>Returns lowest i/o unit number not in use.</summary>
integer function newunit(unit) result(n)
  !<parameter>Out parameter that will contain the lowest i/o number.</parameter>
  integer, intent(out), optional :: unit
  logical inuse
  integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
  integer, parameter :: nmax=999  ! may be system-dependent
  do n = nmin, nmax
     inquire(unit=n, opened=inuse)
     if (.not. inuse) then
        if (present(unit)) unit=n
        return
     end if
  end do
  stop "newunit ERROR: available unit not found."
end function newunit

!***************************************************************************************************
! This routine reads MCevaluation.in
SUBROUTINE read_MCevaluation(MC_data)
type(MC_variables), intent(inout) :: MC_data

integer :: i, status
logical :: err
open(25,file="MCevaluation.in",status="old",iostat=status)
if (status/=0) then
   print *, "couldn't open MCevaluation.in ... exiting."
   stop
endif

call co_ca(25, err)
read(25,*) MC_data%fname
call co_ca(25, err)
read(25,*) MC_data%planes
allocate(MC_data%norm_vec(3,MC_data%planes),MC_data%basis_vec(3,MC_data%planes))
do i=1, MC_data%planes
   call co_ca(25, err)
   read(25,*) MC_data%norm_vec(:,i)
   call co_ca(25, err)
   read(25,*) MC_data%basis_vec(:,i)
enddo

ENDSUBROUTINE read_MCevaluation

!***************************************************************************************************
! This routine writes the correlations of an MCcell and several other structures to 'compare_MC_correlations.out'

SUBROUTINE write_MC_correlations(MCstr,structures)

type(crystal) :: MCstr
type(crystal), pointer :: structures(:)

real(dp), pointer :: PIs(:,:)
character(80),pointer :: titles(:)
integer :: i, j

allocate(PIs(size(MCstr%PI),size(structures)+1))
allocate(titles(size(structures)+1))
titles(1) = MCstr%title(:9)
do i=1,size(structures)
   titles(i+1) = structures(i)%title
enddo

do i = 0,(size(MCstr%PI)-1)
   PIs(i+1,1) = MCstr%PI(i)
enddo

do j = 1,(size(structures))
   do i = 0,(size(MCstr%PI)-1)
      PIs(i+1,j+1) = structures(j)%PI(i)
   enddo
enddo
   
open(33, file="compare_MC_correlations.out", status="replace")
write(33,'("                    ",6(A14))')  titles(:)!MCstr%PI(i),structures(:)%PI(i)
   
do i = 1,(size(MCstr%PI))
   write(33,'("Fig #",i4," Corr=    ",30(f10.7,3x))') i, PIs(i,:)!MCstr%PI(i),structures(:)%PI(i)
enddo

close(33)

ENDSUBROUTINE write_MC_correlations


!***************************************************************************************************
! This routine reads a MCcell.???.out
SUBROUTINE read_MCcell(fname, MCcell, scalefactor, MPIr) !
character(30), intent(in) :: fname
integer(kind=1), pointer :: MCcell(:,:,:,:)
integer, intent(in), optional :: scalefactor(3)
integer, optional :: MPIr
integer :: MPI_rank

integer :: oldMCcellSize(3)

logical :: rescale

integer :: ix, iy, ib, iz !, iAt
!integer :: nx, ny, nz, nb, nAt
integer :: ioerr

if (present(MPIr)) then; MPI_rank = MPIr; else; MPI_rank=0; endif

if (MPI_rank == 0) then
  write(*,'(A)'),   ". reading in MCcell"
  write(*,'(A,A)'), "  file: ", fname
endif

if (present(scalefactor)) then
  if (any(scalefactor/=1)) then; rescale=.true.; else; rescale=.false.; endif
else
  rescale=.false.
endif

open(24,file=fname,status="old",iostat=ioerr)
if (ioerr/=0) then
   print *, "couldn't open ", fname, " ... exiting."
   stop
endif

if (.not. rescale) then
   do ib=1, size(MCcell,4)
      do ix=0, size(MCcell,1)-1
         do iy=0, size(MCcell,2)-1
            read (24, *) MCcell(ix, iy, :, ib)
         enddo
      enddo
   enddo
else
  ! Rescaling
  oldMCcellSize = (/size(MCcell,1),size(MCcell,2),size(MCcell,3)/) / scalefactor
  if (MPI_rank==0) then
    print *
    print *, "WARNING: re-scaling MC cell by factors"
    print *, "WARNING:    ", scalefactor
    print *, "WARNING: old cell size was"
    print *, "WARNING:    ", oldMCcellSize
    print *
    write(*,'(A$)'), "    Reading original cell..."
  endif
  do ib=1, size(MCcell,4)
    do ix=0, oldMCcellSize(1)-1
      do iy=0, oldMCcellSize(2)-1
        read (24, *) MCcell(ix,iy,:oldMCcellSize(3)-1,ib)
      enddo
    enddo
  enddo
  if (MPI_rank==0) then
    write(*,'(A)'), " done."
    write(*,'(A$)'), "    Scaling cell..."
  endif
  do ib=1, size(MCcell,4)
    do ix=0,size(MCcell,1)-1
      do iy=0,size(MCcell,2)-1
        do iz=0,size(MCcell,3)-1
          MCcell(ix,iy,iz,ib) = MCcell(mod(ix,oldMCcellSize(1)), mod(iy,oldMCcellSize(2)), mod(iz,oldMCcellSize(3)),ib)
        enddo
      enddo
    enddo
  enddo
  if (MPI_rank==0) write(*,'(A)'), " done."
  if (MPI_rank==0) call write_MCcell(-1,MCcell)
endif


if (MPI_rank==0) write(*,'(A)') "  done."

close(24)

ENDSUBROUTINE read_MCcell



!***************************************************************************************************
subroutine write_MCcell(temperatureStep, site, isTemporary, TimeForWrite)
! Interface
integer                  :: temperatureStep
integer(kind=1), pointer :: site(:,:,:,:)
logical, optional        :: isTemporary
real(dp), optional       :: TimeForWrite

logical :: createTemporary
real(dp) :: walltime

character(20) :: outfilename, tmpoutfilename
integer :: tmpstatus
integer :: j, m,n, b,f,h,i,k

! Timer
real(dp) walltime1,cputime1           ! Time when starting task
real(dp) walltime2,cputime2           ! Time when task has ended


if (present(isTemporary)) then; createTemporary=isTemporary;
                          else; createTemporary=.false.;
                         endif

if (temperatureStep>=0) then
!LN
   f = (temperatureStep/1000000)
   h = ((temperatureStep - 1000000*f)/100000)
   i = ((temperatureStep - 1000000*f - 100000*h)/10000)
   k = ((temperatureStep - 1000000*f - 100000*h - 10000*i)/1000)
   j = ((temperatureStep - 1000000*f - 100000*h - 10000*i - 1000*k)/100)
   m = ((temperatureStep - 1000000*f - 100000*h - 10000*i - 1000*k - 100*j)/10)
   n = temperatureStep - 1000000*f - 100000*h - 10000*i - 1000*k - 100*j - 10*m

!   j = (temperatureStep/100)
!   m = ((temperatureStep - 100*j)/10)
!   n = temperatureStep -100*j -10*m

   f = f+1
   h = h+1
   i = i +1
   k = k + 1
   j=j+1
   m=m+1
   n=n+1


   outfilename =  "MCcell" //"0123456789" (f:f)//"0123456789" (h:h)//"0123456789" (i:i)//"0123456789" (k:k)//"0123456789" (j:j)//"0123456789" (m:m)//"0123456789" (n:n)// ".out"
!   outfilename =  "MCcell" //"0123456789" (j:j)//"0123456789" (m:m)//"0123456789" (n:n)// ".out"
   tmpoutfilename = outfilename(1:13)//".tmp"
else
   outfilename = "MCcell_start.out"
endif

if (.not. createTemporary) then
  write(*,'(3A,$)') "Writing ", adjustl(trim(outfilename)), "... "
  open(24, FILE=outfilename, status='REPLACE')
else
  write(*,'(3A,$)') "        writing ", adjustl(trim(tmpoutfilename)), "... "
  open(24, FILE=tmpoutfilename, status='REPLACE')
endif

call timing(walltime1,cputime1)

do b=1, size(site,4)
  do m=0, size(site,1)-1
    do n=0, size(site,2)-1
      write (24, '(100I3)', advance = 'no') site(m, n, :, b)
      write (24, '(A, 2I3)') " # x,y: ", m, n
    enddo
  enddo
enddo
close(24)
call timing(walltime2,cputime2)
write(*,'(A$)') "done."

walltime=walltime2-walltime1
write(*,'(A,F6.1,A)') " (", walltime, " s)"
if (present(TimeForWrite)) TimeForWrite=walltime

if (.not. createTemporary) then
  open(25, FILE=tmpoutfilename, status='OLD', iostat=tmpstatus)
  if (tmpstatus==0) then
    write(*,'(3A,$)') "Deleting temporary file ", adjustl(trim(tmpoutfilename)), "... "    
    close(25, status='DELETE')
    write(*,'(A)') "done."
  endif
endif


end subroutine write_MCcell


!***************************************************************************************************
SUBROUTINE read_gss(CE) ! title,cetype,CErank,LV,aBas,aBasCoord,Ninterior,cutoff,eps)
type(CE_variables), intent(inout) :: CE

logical :: err
character(1) dummy
character(3) EOF
integer :: status, iline, k, n, ioerr, ik

open(21,file='groundstatesearch.in',status='old',iostat=ioerr)
if(ioerr/=0) stop "ERROR: Couldn't open the file 'groundstatesearch.in'"
call co_ca(21, err)
read(21,*) CE%GSSspan
!!GH we don't really need to do this in a CE context so let's just
!!always make this true. (Maybe we should remove this eventually.)
!! You can get the same functionality that this variable used to
!! provide by specifying a concentration range of zero for one of the
!! components. 
CE%GSS_full = .true.
!call co_ca(21, err)
!read(21,*) CE%GSS_full
call co_ca(21, err)
read(21,*) CE%GSS_concOnly
call co_ca(21, err)
k = CE%CErank
allocate(CE%GSS_concRange(k,3))
CE%GSS_concRange(:,:) = 0
if (CE%GSS_concOnly) then
   if (CE%nCouplings>1) stop "ERROR in groundstatesearch.in: Concentration ranges not available for coupled CEs."
   call co_ca(21, err)
   do iline = 1, k
      call co_ca(21, err)
      read(21,*) CE%GSS_concRange(iline,1:3)
      !write(*,'(" range ",3(i3,1x))') CE%GSS_concRange(iline,1:3)
   enddo
   if (any(CE%GSS_concRange<0)) stop "ERROR: negative input on concentrations in read_input"
   do ik = 1, k
      if (maxval(CE%GSS_concRange(ik,1:2))>CE%GSS_concRange(ik,3)) then
         write(*,'("ERROR: Numerator is larger than denominator.")')
         write(*,'("ERROR: Check the concentration input for element #:",i2)') ik
         stop
      endif
   enddo
else
   dummy = ""
   do while (dummy/="#")
      read(21,*,iostat=status) dummy
      if (status/=0) stop "Problem skipping lines in groundstatesearch.in"
   enddo
endif
call co_ca(21, err)
read(21,*) CE%GSS_MPIcat, CE%GSS_MPIdelete
call co_ca(21, err)

!--------------------------------------------------------------------------------
! End of file :-)
call co_ca(21, err)
read(21, *) EOF
call ucase(EOF)
if (EOF /= "EOF") stop "ERROR in groundstatesearch.in: format changed. Last entry in CEfitting.in must be 'EOF'. Check input file."

close(21)
END SUBROUTINE read_gss

!*************************************************************************************************
! subroutine:  perform_automaticsearch(structpredlist)
! purpose:     This subroutine checks the predicted groundstates and determines automatically
!              the gss structure numbers of those structures which are energetically most favourable
SUBROUTINE perform_automaticsearch(CE,predictions)
! IN:
type(CE_variables), intent(in) :: CE
! OUT:
type(structure_prediction_t) :: predictions


integer :: poscarfh


! binary tree (both for concentrations and formation energies)
! type def:
type conc_energy_tree
  real(dp), pointer :: x(:,:)
!  real(dp) :: conc           ! concentration
  real(dp) :: formenergy     ! formation energy
  integer  :: strucnum   ! structure number
  type (conc_energy_tree), pointer :: cup, cleft, cright ! pointers for the concentration tree
  type (conc_energy_tree), pointer :: eup, eleft, eright ! pointers for the energy tree
end type conc_energy_tree
type (conc_energy_tree), pointer :: cet_root ! common root of the tree
logical :: create_root ! do we still have to create a root

! data read in from gss.out
type dataset
  integer :: strucnum
  real(dp), pointer :: x(:,:)
!  real(dp):: conc
  real(dp):: CE_energy
  real(dp):: CE_formationenergy
  real(dp) :: CE_referenceenergy
  integer :: different_concentrations ! how many different concentrations were found?
  integer :: total_lines ! how many lines did we read in?
  integer, pointer :: dummy(:)
end type dataset
type (dataset) :: data

! data read in from strucpred.in
!integer :: maxnumber_energies    ! how many low energy structures will be printed out?
!real(dp):: conc_start, conc_stop ! for which concentration range?
!real(dp), pointer :: xi(:,:), xf(:,:)  ! concentration ranges: xi = start, xf = stop
!integer :: lsize, i               ! in case of manual input: size of structure list

! file i/o
integer :: ioerr
logical :: err


poscarfh = 21  ! tk moved the check for the POSCAR.out file up a few lines, so that the code no longer
               ! traverses the gss.out if the POSCAR.out does exist.
open(poscarfh,file="POSCAR.out",status='new',iostat=ioerr)
if (ioerr/=0) then
  print *
  stop "ERROR: POSCAR.out already exists. Remove first."
endif

write(*,'(A)') "Traversing gss.out..."

open(20, file='gss.out', status='OLD', iostat=ioerr)
if (ioerr/=0) stop "ERROR: gss.out was not found!"

allocate(data%x(CE%CErank,CE%ncouplings))
allocate(data%dummy(1+CE%CErank*CE%ncouplings))

data%total_lines=0
data%different_concentrations=1
create_root=.true. ! we still need to create a root for the tree

write(*,'(/,A)') "Progress (one dot every 10,000 strucs):"

do 
  ! read in from gss.out
  call co_ca(20,err)

  ! gss.out structure is:
  ! # struc no. | x(1,1) | x(2,1) | x(1,2) | x(2,2) | nUC | nTyp              |       E(CE) |     DHf(CE) |    Eref(CE) |
  read(20,*,iostat=ioerr) data%strucnum, &  ! structure number
                          data%x(:,:), &    ! concentrations
                          data%dummy(:), &  ! nUC and nTyp
                          data%ce_energy, & 
                          data%ce_formationenergy, &
                          data%ce_referenceenergy
  if (ioerr/=0) exit  ! end of file reached?

  ! lines already read in:
  data%total_lines=data%total_lines+1

  ! output
  if (mod(data%total_lines,  10000)==0) write(*,'(A)',advance='no') "."

  if (data%total_lines==1 .or. create_root) then
    ! create root
    if ( check_concentration( data%x, predictions%xi, predictions%xf ) ) then
!    if (data%conc >= conc_start .and. data%conc <= conc_stop) then
      ! but check wheter the actual concentration is in the desired range
      call create_root_node(cet_root,data%x,data%CE_formationenergy,data%strucnum)
      create_root=.false.
    endif
  else
    ! don't need to create root
    if ( check_concentration( data%x, predictions%xi, predictions%xf ) ) then
!    if (data%conc >= conc_start .and. data%conc <= conc_stop) then
      ! check wheter the actual concentration is in the desired range: 
      !   if yes insert a concentration node and count the number of different concentrations already found
      if (create_conc_node(cet_root,data%x,data%CE_formationenergy,data%strucnum)) &
            data%different_concentrations=data%different_concentrations+1
    endif
  endif

end do
close(20)

write(*,*) 

! finished reading from gss.out
! print summary:
write(*,*)
write(*,*) "Read from gss.out finished"
write(*,*) "-- read datasets                 : ", data%total_lines
write(*,*) "-- different concentrations found: ", data%different_concentrations

! allocate for the structure list
allocate(predictions%struclist(data%different_concentrations*predictions%nSPX))
allocate(predictions%energylist(data%different_concentrations*predictions%nSPX))
predictions%struclist = -1  ! failsafe for not used entries

! Now traverse the whole binary tree in order to get a sorted list of concentrations and formation energies.
write(*,*)
write(*,*) "For each concentration:"
write(*,'(a,i5,a)') "determine the", predictions%nSPX," structures that are lowest in energy:"
write(*,*)


write(poscarfh,'(A)') "#--------------------------------------------------------------------------------"
call conc_tree__sort(cet_root,predictions%nSPX,predictions%struclist,predictions%energylist)
write(poscarfh,'(A)') "#--------------------------------------------------------------------------------"
write(poscarfh,'(A)') "#"
write(poscarfh,'(A)') "#================================================================================"

close(poscarfh)

contains

function check_concentration(x,xi,xf)
logical check_concentration
real(dp), intent(in) :: x(:,:)
real(dp), intent(in) :: xi(:,:), xf(:,:)

integer ik, nk
integer iC, nC

check_concentration = .true.

nk = size(x,1)  ! rank 
if (nk /= size(xi,1) .or. nk /= size(xf,1)) stop "ERROR: bad programming in check_concentration (size x,1)"
nC = size(x,2)  ! number of couplings
if (nC /= size(xi,2) .or. nC /= size(xf,2)) stop "ERROR: bad programming in check_concentration (size x,2)"

! check for each concentration of each coupling, if x is within the (xi,xf) range. Things get more
! complicated since maybe xi > xf
do iC=1,nC
do ik=1,nk
  if ( x(ik,iC) >= xi(ik,iC) .and. x(ik,iC) <= xf(ik,iC) ) cycle
  if ( x(ik,iC) <= xi(ik,iC) .and. x(ik,iC) >= xf(ik,iC) ) cycle
  check_concentration = .false.
  return
enddo
enddo
end function check_concentration


! traverse the concentration tree in a sorted manner (ascending concentrations)
recursive subroutine conc_tree__sort(thisnode,maxenergies,structpredlist,energylist)
type (conc_energy_tree), pointer :: thisnode
type (conc_energy_tree), pointer :: left, right
integer maxenergies
integer, pointer :: structpredlist(:)
real(dp), pointer :: energylist(:)

if (associated(thisnode)) then
  left  => thisnode%cleft
  right => thisnode%cright
  
  call conc_tree__sort(left,maxenergies,structpredlist,energylist)

  write(*,'(a,10(f12.6,1x))') "concentrations (1:k,nC)= ", thisnode%x !, thisnode%strucnum
  write(*,'(4a20)') ,"prediction # |","energy # |","formenergy |","struc # |"
  write(21,'("#",1x,a,10f10.4)') "concentrations (1:k,nC) = ", thisnode%x !, thisnode%strucnum
  write(21,'("#",4a20)'),"prediction # |","energy # |","formenergy |","struc # |"

  call energy_tree__sort(thisnode,maxenergies,.true.,structpredlist,energylist)     

  call conc_tree__sort(right,maxenergies,structpredlist,energylist)
endif

end subroutine conc_tree__sort

! traverse the energy tree in a sorted manner (ascending energies)
recursive subroutine energy_tree__sort(thisnode,maxenergies,reset,structpredlist,energylist)
type (conc_energy_tree), pointer :: thisnode
type (conc_energy_tree), pointer :: left, right
real(dp) value
integer strucnum

integer maxenergies
integer, pointer :: structpredlist(:)
real(dp), pointer :: energylist(:)

integer, save :: energy_serial = 0, conc_energy_serial 
logical reset

if (reset) conc_energy_serial=0
if (associated(thisnode)) then
  value = thisnode%formenergy
  strucnum = thisnode%strucnum
  left  => thisnode%eleft
  right => thisnode%eright
  
  call energy_tree__sort(left,maxenergies,.false.,structpredlist,energylist)

  conc_energy_serial=conc_energy_serial+1

  if (conc_energy_serial<=maxenergies)  then
    energy_serial     =energy_serial+1

    write(*,'(i20,i20,f20.5,i20)')  , energy_serial,conc_energy_serial,value, strucnum
    write(21,'("#",i20,i20,f20.5,i20)') , energy_serial,conc_energy_serial,value, strucnum

    structpredlist(energy_serial)=strucnum ! add this struture to the structpredlist
    energylist    (energy_serial)=value    ! add the respective formation energy to the energylist
    
  endif

  call energy_tree__sort(right,maxenergies,.false.,structpredlist,energylist)
endif

end subroutine energy_tree__sort


! create a root for the binary tree
subroutine create_root_node(thisnode,cvalues,evalue,strucnum)
type (conc_energy_tree), pointer :: thisnode
real(dp), intent(in)             :: cvalues(:,:), evalue
integer, intent(in)              :: strucnum

integer :: km, nCm   ! CE rank -1, ncouplings -1

km  = size(cvalues,1)
nCm = size(cvalues,2)

allocate(thisnode)
allocate(thisnode%x(km,nCm))

thisnode%x = cvalues
thisnode%formenergy = evalue
thisnode%strucnum = strucnum

nullify(thisnode%cleft,thisnode%cright,thisnode%cup)
nullify(thisnode%eleft,thisnode%eright,thisnode%eup)

end subroutine create_root_node



! create a node in the concentration tree
logical function create_conc_node(root,cvalues,evalue,strucnum)
type (conc_energy_tree), pointer :: root
type (conc_energy_tree), pointer :: thisnode, prevnode, newnode
real(dp), intent(in)             :: cvalues(:,:), evalue
integer, intent(in)              :: strucnum
logical left, makenewconcnode

integer :: cmp ! compare status of the arrays

integer :: km, nCm   ! CE rank -1, ncouplings -1
km  = size(cvalues,1)
nCm = size(cvalues,2)

thisnode => root
left=.false.
do while (associated(thisnode))

  call cmp_arrays(cvalues,thisnode%x,cmp)

  if (cmp < 0) then
    ! if the read concentration is less than the concentration represented in the node:
    ! add a node left of the current one
    prevnode => thisnode
    thisnode => thisnode%cleft
    left=.true.
    makenewconcnode=.true.
  elseif (cmp > 0) then
    ! if the read concentration is greater than the concentration represented in the node:
    ! add a node right of the current one
    prevnode => thisnode
    thisnode => thisnode%cright
    left=.false.
    makenewconcnode=.true.
  elseif (cmp == 0) then
    ! if the read concentration is equal to the concentration represented in the node:
    ! do not add a new concentration node, but add an energy node. The binary tree for
    ! the energies is attached to the respective concentration node
    makenewconcnode=.false.
    call create_energy_node(thisnode,cvalues,evalue,strucnum)
    exit
  else
    stop "Bad programming (lkjhf)"
  end if
end do

! here we actually add the concentration node:
if (makenewconcnode) then
  allocate(newnode) ! new node
  allocate(newnode%x(km,nCm))
  newnode%x = cvalues ! update values
  newnode%formenergy = evalue
  newnode%strucnum = strucnum
  nullify(newnode%cleft,newnode%cright,newnode%cup)
  nullify(newnode%eleft,newnode%eright,newnode%eup)
  newnode%cup => prevnode ! cup = Concentration UP (upwards pointer)
  
  if (left) then
    prevnode%cleft => newnode ! cleft = Concentration LEFT
  else
    prevnode%cright => newnode ! cright = Concentration RIGHT
  endif
endif

! did we create a node?
create_conc_node = makenewconcnode
return

end function create_conc_node

! create an energy node
! the binary trees for the energies is attached to the respective concentration nodes, 
! thus the "root" of the energy tree is the concentration node
subroutine create_energy_node(root,cvalues,evalue,strucnum)
type (conc_energy_tree), pointer :: root
type (conc_energy_tree), pointer :: thisnode, prevnode, newnode
real(dp), intent(in)             :: cvalues(:,:), evalue
integer                          :: strucnum
logical left, makenewenergynode

integer :: km, nCm   ! CE rank -1, ncouplings -1
km  = size(cvalues,1)
nCm = size(cvalues,2)

thisnode => root
left=.false.
do while (associated(thisnode))
  if (evalue < thisnode%formenergy) then ! smaller formation energy?
    prevnode => thisnode
    thisnode => thisnode%eleft
    left=.true.
    makenewenergynode=.true.
  elseif (evalue >= thisnode%formenergy) then ! greater or equal formation energy?
    prevnode => thisnode
    thisnode => thisnode%eright
    left=.false.
    makenewenergynode=.true.
!  else 
!    makenewenergynode=.false.
!    exit
  end if
end do

! add new node
if (makenewenergynode) then
  allocate(newnode)
  allocate(newnode%x(km,nCm))
  newnode%x = cvalues
  newnode%formenergy = evalue
  newnode%strucnum = strucnum
  nullify(newnode%cleft,newnode%cright,newnode%cup)
  nullify(newnode%eleft,newnode%eright,newnode%eup)
  newnode%eup => prevnode
  
  if (left) then
    prevnode%eleft => newnode
  else
    prevnode%eright => newnode
  endif
endif

end subroutine create_energy_node



ENDSUBROUTINE perform_automaticsearch
!*************************************************************************************************
SUBROUTINE read_MCpar(MC_data, rank, ncouplings, nreferences, eps, xeps, MPI_rank)
use utilities_module, only: ralloc

type(MC_variables) :: MC_data
integer, intent(in) :: rank  ! Rank of the cluster expansion (binary, ternary, etc.)
integer, intent(in) :: ncouplings ! how many CEs are coupled?
integer, intent(in) :: nreferences ! number of CE references
real(dp), intent(in) :: eps   ! finite precision checking
real(dp), intent(in) :: xeps  ! finite precision checking for concentrations
integer, intent(in) :: MPI_rank
real(dp), pointer :: mu(:,:)

logical err
integer :: ioerr, i,j, iC, iMuTemperatures, iMuX

integer :: itmp, nE
integer, pointer :: canonicalSpecies(:)
character(1) :: ctmp
character(30) :: boxsize, c_tmp
character(1000) :: line
real(dp), pointer :: numbers(:)


! substantial settings (had to be moved up to this position here)
MC_data%eps        = eps
MC_data%rank       = rank 
MC_data%ncouplings = ncouplings


open(19, file='MCpar.in', status='OLD', iostat=ioerr)
if (ioerr/=0) stop "ERROR: MCpar.in was not found!"

call get_SpinRank_mapping(rank, MC_data%sp, MC_data%rk)

!--------------------------------------------------
! Read: Cell size
!--------------------------------------------------
call co_ca(19,err)

read(19,'(A1000)') line
call parse_line_for_numbers(line," ",maxEntries=9,readEntries=nE,values=numbers)
if (nE<3) stop "ERROR in MCpar.in: The simulation cell size must be given by 3 integer values."
MC_data%boxSize(1:3) = numbers(1:3)
if (nE>3 .and. nE/=9) stop "ERROR in MCpar.in: You must enter 3 or 9 values to determine the simulation cell."
if (nE==3) then
  MC_data%simOffset(1:3)=0
  MC_data%simSize(1:3)=MC_data%boxSize(1:3)
elseif (nE==9) then
  MC_data%simOffset(1:3)=numbers(4:6)
  MC_data%simSize(1:3)=numbers(7:9)
else
  stop "ERROR in reading MCpar.in. Bad programming."
endif
if ( any(MC_data%simOffset+MC_data%simSize>MC_data%boxSize) ) stop "ERROR in MCpar.in: the restricted simulation cell *must* entirely lie withing the first unit cell of the simulation."
  ! NB: If the "restricted simulation cell" that is defined by simSize and simOffset would not entirely lie within 
  ! the "boxSize" (i.e., the first unit cell of the simulation superlattice), then the site selection process would
  ! need an additional modulo operation. It could be easily included, but I (tk) am a little bit reluctant at the moment
  ! to introduced yet another source of slowing the algorithm down. If it is absolutely needed, maybe we can simply have a
  ! standalone compilation for that.

! convert boxsize to a printable (character) format
write(c_tmp,'(I5)') MC_data%boxsize(1)
boxsize=trim(adjustl(c_tmp))
write(c_tmp,'(I5)') MC_data%boxsize(2)
boxsize=trim(adjustl(boxsize))//" x "//trim(adjustl(c_tmp))
write(c_tmp,'(I5)') MC_data%boxsize(3)
boxsize=trim(adjustl(boxsize))//" x "//trim(adjustl(c_tmp))
MC_data%cboxSize=boxsize

!--------------------------------------------------
! Read: Concentrations and Ensemble
!--------------------------------------------------
allocate(MC_data%x(rank,ncouplings))
do i=1,ncouplings
  call co_ca(19,err)
  read(19,*,iostat=ioerr) MC_data%x(:,i) ! Read in the concentration for each atom type and each CE
  if (ioerr/=0) stop "ERROR in MCpar.in: supply a concentration for each atomic species and each coupling"
  if (.not. equal(sum(MC_data%x(:,i)),1._dp,xeps)) then
    if (MPI_rank==0) then
      write(*,'(/,A,I3,A)')      "WARNING: sum of concentrations of CE# ", i, " does not add up to 1:"
      write(*,'(A,3(F10.5,2x))') "WARNING:   x(:) = ", MC_data%x(:,i)
      write(*,'(A)')             "WARNING: renormalizing the concentrations to"
    endif
    MC_data%x(:,i) = MC_data%x(:,i)/sum(MC_data%x(:,i))
    if (MPI_rank==0) then
      write(*,'(A,3(F10.5,2x))') "WARNING:   x(:) = ", MC_data%x(:,i)
    endif
  endif
enddo

! one more failsafe about the concentrations, should never happen:
if (.not. all( (/ (equal(sum(MC_data%x(:,iC)),1._dp,xeps),iC=1,ncouplings) /) ) ) then
  write(*,'(A)')                "ERROR in MCpar.in: supplied concentration do not add up to 1"
  write(*,'(A,100(F15.13,1x))') "ERROR            : concentrations = ", MC_data%x
  write(*,'(A,F15.13)')         "ERROR            : sum = ", sum(MC_data%x)
  stop
end if


call co_ca(19,err)
read(19,*) ctmp
call ucase(ctmp)
   if (ctmp=='G') then
     if (ncouplings>1) stop "No GC MC for coupled systems yet"
     MC_data%GC  =.true.;   ! this is a gMC run
     MC_data%pGC = .false.
     MC_data%mode = 1
     allocate(canonicalSpecies(0)) ! temporary alloc, needed for the correct setup of the possibleNewSpins
     call read_gMC_parameters ! outsourced the reading of all the mu parameters for the gMC.
                              ! (note: bad subroutine, since it simply uses the variables of this subroutine here)
   elseif (ctmp=='P') then ! partical grand canonical
     if (ncouplings>1) stop "No pGC MC for coupled systems yet"
     MC_data%GC  = .true.
     MC_data%pGC = .true.
     MC_data%mode = 3

     call co_ca(19,err)
     read(19,'(A1000)') line
     call parse_line_for_numbers(line," ",maxEntries=rank,readEntries=nE,values=numbers)
     allocate(MC_data%pGCspecies(nE),MC_data%GCspins(nE))
     allocate(canonicalSpecies(rank))

     MC_data%pGCspecies(:) = int( numbers(:) )  ! the atomic species that will run in GC mode

     canonicalSpecies = (/ (i,i=1,rank) /)      ! the atomic species that will run in C mode: start setup...
     do i=1,nE; if (any(MC_data%pGCspecies(:)==canonicalSpecies(i))) canonicalSpecies(i) = -1; enddo
     canonicalSpecies = pack( canonicalSpecies, canonicalSpecies/=-1 )
     canonicalSpecies => ralloc(canonicalSpecies, rank-nE) ! ... finished

     MC_data%GCspins(:)    = MC_data%sp( MC_data%pGCspecies )
     call read_gMC_parameters ! outsourced the reading of all the mu parameters for the gMC.
                              ! (note: bad subroutine, since it simply uses the variables of this subroutine here)
     deallocate(canonicalSpecies)
   elseif (ctmp=='C') then ! canonical
     MC_data%GC  = .false.
     MC_data%pGC = .false.
     MC_data%mode = 2
   else
     stop "ERROR in MCpar.in: Specify 'G' or 'C' for Grandcanonical or Canonical ensemble"
   endif

!--------------------------------------------------
! Read: temperature parameters
!--------------------------------------------------
call co_ca(19,err)
read(19,*) itmp
MC_data%nTemperatureSchedule=itmp
allocate(MC_data%TemperatureMode(itmp))
allocate(MC_data%startT(itmp), MC_data%endT(itmp), MC_data%changeT(itmp))
allocate(MC_data%changeScheduleT(itmp+1))

do i=1, itmp
  call co_ca(19,err)
  read(19,*) MC_data%startT(i)
  if (i>1) MC_data%endT(i-1) = MC_data%startT(i)
  call co_ca(19,err)
  read(19,*) ctmp, MC_data%changeT(i)
  if (ctmp=="*") then;        MC_data%TemperatureMode(i) =  0;
     elseif (ctmp=="+") then; MC_data%TemperatureMode(i) =  1;
     elseif (ctmp=="-") then; MC_data%TemperatureMode(i) = -1;
     else; stop "ERROR in MCpar.in: Specify '*', '-' or '+' for the temperature mode"
  endif
end do 
call co_ca(19,err)
read(19,*) MC_data%endT(itmp)

MC_data%changeScheduleT(1:itmp) = MC_data%startT(:)
MC_data%changeScheduleT(itmp+1) = MC_data%endT(itmp)

!--------------------------------------------------
! Read: MC parameters
!--------------------------------------------------
call co_ca(19,err)
read(19,*) MC_data%avsteps
call co_ca(19,err)
read(19,*) MC_data%E_conv
   if (MC_data%E_conv < 0) then
     MC_data%maxnSteps = - int(MC_data%E_conv,li)
     MC_data%E_conv    = 0.d0
   else
     MC_data%maxnSteps = -1
   endif

!--------------------------------------------------
! Read: Performance parameters
!--------------------------------------------------
call co_ca(19,err)
read(19,*) itmp
   if (itmp<10 .and. itmp>0) then
     MC_data%nDepLists = itmp
   else
     stop "Error in MCpar.in: number of dependency lists must be >0, <10"
   endif


call co_ca(19,err)
read(19,*) MC_data%depBoxLength
if ( mod(MC_data%boxsize(1),MC_data%depBoxLength(1)) /= 0 .or. &
     mod(MC_data%boxsize(2),MC_data%depBoxLength(2)) /= 0 .or. &
     mod(MC_data%boxsize(3),MC_data%depBoxLength(3)) /= 0 )    &
   stop "Error in MCpar.in: length of dependency box and length of MC cell must be commensurate."

call co_ca(19,err)
read(19,*,iostat=ioerr) MC_data%onlyTestSpeed, MC_data%SIstepsFactor
if (ioerr /= 0) MC_data%SIstepsFactor = 1 

!--------------------------------------------------
! Read: Write out parameters
!--------------------------------------------------
call co_ca(19,err)
read(19,*) MC_data%nTsteps_write_MCcell

!--------------------------------------------------
! Read: continue parameters
!--------------------------------------------------
call co_ca(19,err)
read(19,*) MC_data%continue

if (MC_data%continue) then
  if (nreferences>0) stop "ERROR: no MC continuation for reference-supported CE's yet. Sorry."

  call co_ca(19,err)
  read(19,*) MC_data%continue_fileOrStruc   ! continue from File or from Structure?
  call ucase(MC_data%continue_fileOrStruc)

  if (MC_data%continue_fileOrStruc=='F') then
    ! continue from file. This is equivalent to the situation that had always been there in MCpar.in.
    ! The MC run is continued from a file that had been written (probably) by a previous MC run.
    call co_ca(19,err)
    read(19,*) MC_data%continue_fileName   ! read filename
    MC_data%continue_fileName=adjustl(MC_data%continue_fileName)
  elseif (MC_data%continue_fileOrStruc=='S') then
    ! continue from structure. This is 'new' (May 2012). You can also specify a structure in structures.in-format. 
    allocate(MC_data%continue_structure(1))  ! unfortunately, need a pointer for this 1 struc
    call read_structure_poscar(19,MC_data%continue_structure(1),MC_data%rank,MC_data%nCouplings,.false.,.false.,ioerr,.false.,justReadStructureDetails_=.true.)
    ! mind: pLV of the structure is not yet set!!
  else
    stop "ERROR in MCpar.in: continue parameter 'File or Struc' must be set to either 'F' or 'S'. (MCpar.in changed format May 2nd, 2012)"
  endif

  call co_ca(19,err)
  read(19,*) MC_data%continue_temperatureStep
  call co_ca(19,err)
  read(19,*) MC_data%continue_temperature
  call co_ca(19,err)
  read(19,*) MC_data%continue_energy
  call co_ca(19,err)
  read(19,*) MC_data%continue_energyConverged
  call co_ca(19,err)
  read(19,*,iostat=ioerr) MC_data%continue_scale(:)
  if (ioerr /= 0) MC_data%continue_scale = 1
  if (any(mod(MC_data%boxsize,MC_data%continue_scale)/=0)) &
       stop "ERROR: scale factor for continue-run must be commensurate with MC cell size"
  if (any(MC_data%continue_scale<0)) &
       stop "ERROR: for now, the scale factors must be >0, i.e. up-scaling"
endif



! References:
allocate( MC_data%reference(nreferences))
allocate( MC_data%g2gref   (nreferences))


close(19)


contains

! outsourced routine for the read-in of the mu parameters.
! this is a "bad" subroutine since it uses the variables of the main subroutine.
subroutine read_gMC_parameters
integer k
     ! now set the stage for reading the mu(T,x)-table:
     ! - how many temperatures are supplied?
     ! - how many concentrations are supplied?
     call co_ca(19,err)
     read(19,*) MC_data%nMuTemperatures, MC_data%nMuX

     ! allocate all the parameters:
     allocate(MC_data%mu(lbound(MC_data%rk,1):ubound(MC_data%rk,1),MC_data%nMuTemperatures,MC_data%nMuX))
     allocate(MC_data%muTemperatures(MC_data%nMuTemperatures))
     allocate(MC_data%muX(rank,MC_data%nMuX))
     allocate(MC_data%possibleNewSpins(0:rank-2-size(canonicalSpecies),lbound(MC_data%rk,1):ubound(MC_data%rk,1)))  ! canonical species are excluded in all possibleNewSpins lists
     allocate(MC_data%Dmu(lbound(MC_data%rk,1):ubound(MC_data%rk,1),lbound(MC_data%rk,1):ubound(MC_data%rk,1)))
     allocate(mu(rank,MC_data%nMuX))

     ! read the concentrations (for all ranks!)
     call co_ca(19,err)
     read(19,*) MC_data%muX(:,:)    ! i.e.,   x1A x1B    x2A x2B    ....   where 1,2 is the concentration counter, A,B the different atomic species  

     ! for each temperature:
     ! - read the temperature value
     ! - read the mu(T,x)
     do iMuTemperatures=1,MC_data%nMuTemperatures
       call co_ca(19,err)
       read(19,*) MC_data%muTemperatures( iMuTemperatures ), mu(:,:)

       do iMuX=1,MC_data%nMuX
       do i=1, rank
         MC_data%mu( MC_data%sp(i) , iMuTemperatures, iMuX ) = mu(i,iMuX)
       enddo ! i
       enddo ! iMuX

       do i=1, rank
         MC_data%possibleNewSpins(:,MC_data%sp(i)) = pack(MC_data%sp(:),MC_data%sp/=MC_data%sp(i) .and. .not. (/ (any(canonicalSpecies==k),k=1,rank) /)) 
                        ! the new possible spins for each spin:         <---------------------->        <-------------------------------------------->
                        !                                               should not be the spin type     should not include those spins that belong to 
                        !                                               that the site had before        an atomic species that is simulated canonically
                        !                                               the flip                        (otherwise spins could become "canonical spins"
                        !                                                                               but never the other way round)
       enddo

       
     enddo ! iMuTemperatures

end subroutine read_gMC_parameters

ENDSUBROUTINE read_MCpar

SUBROUTINE read_MDpar(MD_data)

type(MD_variables) :: MD_data
integer :: ioerr, i,j, iC, iMuTemperatures, iMuX

integer :: itmp, nE
integer, pointer :: canonicalSpecies(:)
character(1) :: ctmp
character(30) :: boxsize, c_tmp
character(1000) :: line
real(dp), pointer :: numbers(:)
logical err

open(19, file='MDpar.in', status='OLD', iostat=ioerr)
if (ioerr/=0) stop "ERROR: MDpar.in was not found!"

print *, 'here'
call co_ca(19,err)
read(19,*) MD_data%Temperature
call co_ca(19,err)
read(19,*) MD_data%gaussian_height, MD_data%gaussian_variance
call co_ca(19,err)
read(19,*) MD_data%deltaT, MD_data%deltaI, MD_data%minheight
call co_ca(19,err)
read(19,*) MD_data%num_of_flips_between_gaussian_adds
!call co_ca(19,err)
!read(19,*) MD_data%flatness_criteria
!call co_ca(19,err)
!read(19,*) MD_data%x(:)
!call co_ca(19,err)
!read(19,*) MD_data%boxSize(:)
!call co_ca(19,err)
!read(19,*) MD_data%orientation_vector(:)
call co_ca(19,err)
read(19,*) MD_data%num_of_discrete_cv_values_x
call co_ca(19,err)
read(19,*) MD_data%num_of_discrete_cv_values_y

!MD_data%num_of_discrete_cv_values = MD_data%num_of_discrete_cv_values_x

allocate(MD_data%bias_potential(MD_data%num_of_discrete_cv_values_x,MD_data%num_of_discrete_cv_values_y))
MD_data%bias_potential = 0
MD_data%bias_potential_division_x = 2._dp/(MD_data%num_of_discrete_cv_values_x - 1)  ! Hard coded a 2 in here for the case where the CV is concentration.  For other cases this might not be true.
MD_data%bias_potential_division_y = 2._dp/(MD_data%num_of_discrete_cv_values_y - 1)  ! Hard coded a 2 in here for the case where the CV is concentration.  For other cases this might not be true.
MD_data%bias_potential_division = MD_data%bias_potential_division_x

END SUBROUTINE read_MDpar
!*************************************************************************************************
SUBROUTINE read_clusters(Clusters,clustersFile,GA,eps)
character(30), intent(in) :: clustersFile
type(GA_variables), intent(inout) :: GA
type(figure), pointer :: Clusters(:)
real(dp), intent(in) :: eps
integer :: nF, ifg, nV, iV, ioerr
logical err
character(80) line

! Read in the figures from the clusters.out file
open(10,file=clustersfile,status="old",iostat=ioerr)
if (ioerr /= 0) then
  write(*,'(A,A)') "ERROR: Couldn't open clusters file ", trim(adjustl(clustersfile))
  stop
endif

! First, count how many there are
nF = 1   ! always have the constant cluster
do 
   call co_ca(10,err)
   read(10,*,iostat=ioerr) ! Skip the figure number
   if (ioerr/=0) exit
   call co_ca(10,err)
   read(10,*,iostat=ioerr) nV ! Read in the number of vertices
   call co_ca(10,err) 
   read(10,*,iostat=ioerr) ! Read in the "size" of the figure
   do iv = 1, nV ! Read in the vertex coordinates
     call co_ca(10,err)
     read(10,*,iostat=ioerr)
   enddo
   nF = nF+1
enddo                                                        

!do
!   read(10,'(a80)',iostat=ioerr) line
!   if (ioerr /= 0) exit
!   if (index(line,"Cluster number") >0) nF = nF + 1
!enddo

rewind(10)
GA%nTotClusters = nF

! Now read them in and put them in the clusters variable
allocate(clusters(0:nF-1))   ! 0 = constant cluster, not listed in clusters.out

! setup the constant cluster
allocate(clusters(0)%vertex(3,0))
allocate(clusters(0)%label(0))
allocate(clusters(0)%s(0))
clusters(0)%nV   = 0
clusters(0)%avgR = 0

! setup all the other clusters
do ifg = 1, nF-1
   call co_ca(10,err)
   read(10,*) ! Skip the figure number
   call co_ca(10,err)
   read(10,*) nV ! Read in the number of vertices
   clusters(ifg)%nV = nV
   !if (nV == 2) GA%nrpairs = GA%nrpairs + 1
   call co_ca(10,err) 
   read(10,*) clusters(ifg)%avgR ! Read in the "size" of the figure
   allocate(clusters(ifg)%vertex(3,nV))
   allocate(clusters(ifg)%label(nV))
   allocate(clusters(ifg)%s(nV))
   do iv = 1, nV ! Read in the vertex coordinates
     call co_ca(10,err)
      read(10,*) clusters(ifg)%vertex(:,iv), clusters(ifg)%label(iv), clusters(ifg)%s(iv)
   enddo
   call sort_figure_vertices(clusters(ifg),eps,sortS=.true.) ! We can't assume the vertices are already sorted if 
enddo                                                        ! the figures file is "hand-made" rather than generated
close(10)


ENDSUBROUTINE read_clusters
!******************************************************************************************





!******************************************************************************************
subroutine check_enum_file(CE,istat,MPI_rank)
type(ce_variables), intent(in) :: CE
integer, intent(out) :: istat      ! 0   = ok
                                   ! >0  = enum file exists but is inconsistent
                                   ! <0  = enum file does not exist
integer, intent(in) :: MPI_rank

integer ioerr, nBas, iBas, idx, i
logical err, concCheck, full
real(dp) :: a(3,3), eps
real(dp), pointer :: pBas(:,:)
integer span(2), rank
character(10) title
character(1)  type
character(1000) line
real(dp), pointer :: numbers(:)
integer iE,nE
integer concElement
integer :: concRange(3)

logical :: pr

if (MPI_rank==0) then; pr = .true.; else; pr = .false.; endif

if (pr) then
  print *
  write(*,'(A)') "Checking the consistency of the enumeration file..."
  write(*,'(A)') "Note: If this check fails try to remove struct_enum.out first."
  print *
endif

eps = CE%eps
istat = 0

open(14,file="struct_enum.out",status='old',iostat=ioerr)
if (ioerr/=0) then
   istat = -1
   if (pr) print *, "WARNING: struct_enum.out was not found, will create a new one!"
   close(14)
   return
endif

call co_ca(14,err)
read(14,*) title
call co_ca(14,err)
read(14,*) type
call ucase(type)
if (.not. CE%CEtype == type) then
  print *, "CE type failed"
  istat = istat+1
  close(14)
  return
endif
call co_ca(14,err)
read(14,*) a(:,1)
call co_ca(14,err)
read(14,*) a(:,2)
call co_ca(14,err)
read(14,*) a(:,3)
if (.not. equal(a,CE%pLV, eps)) then
  print *, "pLV failed"
  istat = istat+1
  close(14)
  return
endif
call co_ca(14,err)
read(14,*) nBas
if (nBas /= CE%nBas) then
  print *, "nBas failed"
  istat = istat +1
  close(14)
  return
endif
allocate(pBas(3,nBas))
do iBas=1,nBas
  call co_ca(14,err)
  read(14,*) pBas(:,iBas)
  if (.not. equal(pBas(:,iBas),CE%pBas(:,iBas),eps)) then
    print *, "dset failed"
    write(*,'(A,I3,A,3(F10.5,1x),5x,3(F10.5,1x))') "dvector #",iBas,":", pBas(:,iBas), CE%pBas(:,iBas)
    istat = istat+1
    close(14)
    return
  endif
end do
deallocate(pBas)
call co_ca(14,err)
read(14,'(i2)') rank
call co_ca(14,err)
read(14,*) span
if (.not. all(span==CE%GSSspan)) then
  print *, "span failed"
  istat = istat+1
  close(14)
  return
endif

call co_ca(14,err)
read(14,*) eps
if (.not. eps == CE%eps) then
  print *, "eps failed"
  istat = istat+1
  close(14)
  return
endif

call co_ca(14,err)
read(14,'(A100)') line

call co_ca(14,err)
read(14,*) concCheck

if (concCheck) then
  call co_ca(14,err)
  read(14,'(A100)') line ! Skip the comment line
  do i = 1, rank
     read(14,'(A100)') line
     idx = index(line,"/")
     read(line(idx-2:idx-1),*) concRange(1)
     idx = index(line,"/",back=.true.)
     read(line(idx-2:idx-1),*) concRange(2)
     read(line(idx+1:),*) concRange(3)
  enddo

  if (concCheck .neqv. CE%GSS_concOnly &
    .or. concElement /= CE%GSS_concN &
    .or. .not. all(concRange==CE%GSS_concRange(1,:))) then

    print *, "Concentration check failed"
    istat=istat+1
    close(14)
    return

  endif
endif ! concentration check

call co_ca(14,err)
read(14,*) line ! only ready until first " "
call ucase(line); line=adjustl(trim(line))
if (line(1:1)=='P') then; full=.false.; else; full=.true.; endif
if (full .neqv. CE%GSS_full) then
  print *, "Full labeling check failed"
  istat=istat+1
  close(14)
  return
endif

call co_ca(14,err)
read(14,'(A1000)') line
if (line(1:1)=='E') then ! check only if Equivalency line is present (legacy if)
  line=line(index(line,":")+1:)
  call parse_line_for_numbers(line," ",maxEntries=nBas,readEntries=nE,values=numbers)
  if (nE==0) stop "ERROR: parsing line in struct_enum.out"
  if (nE/=nBas) then
    print *, "number of equivalencies failed"
    istat=istat+1
    close(14)
    return
  endif
  if (any(nint(numbers(:nE))/=CE%equivalencies)) then
    print *, "equivalencies failed"
    istat=istat+1
    close(14)
    return
  endif
endif

close(14)

end subroutine check_enum_file
!***************************************************************************************************


!***************************************************************************************************
SUBROUTINE read_finalECIs(CE, cluster, Jav,Jff)
type(ce_variables) :: CE
type(figure), pointer :: cluster(:)
real(dp), pointer     :: Jav(:), Jff(:)   
real(dp)              :: tmpFinalECI0, tmpFinalECIfullFit0
logical err
character(80) :: filename
character(2)  :: cObs
integer :: iObs, nObs
integer :: rank, isites, Nsites, iV, nV, nFig, figNum, iFig, nSvec

!print *, "read_finalECIs"

nObs = CE%nObservables
if (nObs>1) stop "revise read in of ecis for nObs>1"
! for nObs>1: something like that:
!do iObs=1,nObs
!  write(cObs,'(I2)') iObs
!  filename=trim(CE%J_file)//"."//trim(adjustl(cObs))//".out"
!  call read_J_file(filename, cluster, Jav, Jff)
!enddo

call read_J_file(CE%J_file, cluster, Jav, Jff)

ENDSUBROUTINE read_finalECIs


!***************************************************************************************************
SUBROUTINE read_input_structures(str, CE,filename)
type(crystal), pointer :: str(:)
type(crystal) :: tStr
!type(crystal) :: tStr(1000)
type(CE_variables), intent(inout) :: CE
character(30), intent(in),optional :: filename


real(dp) :: energy
integer iStr, iLV, iL, iAt, c, iTyp, ic, i
integer CErank, indx, ioerr
integer nStr, nAt, nTyp, nSpin, nPos
character(1) :: aBasCoord
character(1000) :: line
logical err, peratom, useweights
integer, pointer :: gettype(:), rankSpace(:)

character(1) :: cTransform
logical :: doTransform
real(dp) :: TransformMatrix(3,3), TransformShift(3)

character(1) :: structFormat
integer :: nSpecialOptions
integer :: iRef, nRef

real(dp), pointer :: numbers(:)


print *,"Reading in the input structures" 

if (present(filename)) then
   open(13,file=filename,status="old",iostat=ioerr)
else
   open(13,file="structures.in",status="old",iostat=ioerr)
end if
if (ioerr/=0) stop "ERROR: Couldn't find 'structures.in' file"
open(14,file="debug_readinstructures.out")
nRef = CE%nReferences

!================================================================================
! Special options for structures.in

nSpecialOptions = 3  ! specify the number here! This is for skipping lines afterwards.

!sascha>! (1)
!sascha>call co_ca(13,err)
!sascha>read(13,*) line ! do you want an expansion for energies only or nonscalar quantities?
!sascha>call ucase(line)
!sascha>if(line=="SCALAR") then; CE%scalar = .true.; CE%nObservables = 1;
!sascha>elseif(line=="MULTIDIMENSIONAL") then 
!sascha>   CE%scalar = .false.
!sascha>   call co_ca(13,err)
!sascha>   nSpecialOptions = nSpecialOptions + 1
!sascha>   read(13,*) CE%nObservables ! read in the number of additional observables, that is _including_ energy
!sascha>   write(*,*) "attempting multidimensional cluster-expansion... user beware"
!sascha>   write(*,*) "structures.in says you are trying to clusterexpand on the free energy and ", CE%nObservables-1, " additional observables"
!sascha>   ! failsafe
!sascha>   if (CE%nObservables<2) stop "ERROR: nObservables must be >1."
!sascha>else; stop "ERROR: in the structures.in file, the first entry should be ""scalar"" or ""multidimensional"""
!sascha>endif
!sascha>

CE%scalar= .true.

!(1) total/peratom
call co_ca(13,err)
read(13,*) line ! are the energies total energies? or given "per atom"?
call ucase(line)
if(line=="PERATOM") then; peratom = .true.
elseif(line=="TOTAL") then; peratom = .false.
else; stop "ERROR: in the structures.in file, the first entry should be ""peratom"" or ""total"""
endif

!(2) weights/noweights
call co_ca(13,err)
read(13,*) line
call ucase(line)
if(line=="WEIGHTS") then; useweights = .true.;
elseif(line=="NOWEIGHTS") then; useweights = .false.
else; stop "ERROR: in the structures.in file, the second entry should be ""weights"" or ""noweights"""
endif

!this flag switches the expected layout of the structures.in file from the VASP format (POSCAR) to the DFTB format (.GEN)

!(3) poscar/gen
call co_ca(13,err)
read(13,*) line
call ucase(line)
if(line=="GEN") then; structFormat='G';
elseif(line=="POSCAR") then; structFormat='P'
else; stop "ERROR: in the structures.in file, the fourth entry should be ""POSCAR"" or ""GEN"""
endif

!================================================================================ 

!depracated>! count how many atom *types* there are
!depracated>call co_ca(13,err)
!depracated>do iL = 1,5; read(13,*); call co_ca(13,err); enddo  ! Skip the first 5 lines
!depracated>read(13,*,iostat=ioerr) cTransform ! transformation wanted?
!depracated>call ucase(cTransform)
!depracated>if (cTransform=='T') then
!depracated>  do iL = 1,4; read(13,*);enddo  ! Skip 4 more lines
!depracated>else
!depracated>  backspace(13)
!depracated>endif
!depracated>cTransform=""
!depracated>
!depracated>read(13,'(a1000)') line ! Read in the number of atoms of each type
!depracated>call parse_line_for_numbers(line," ",maxEntries=CE%CErank*CE%ncouplings,readEntries=c,values=numbers)
!depracated>CErank = c / CE%ncouplings
!depracated>!LNwrite(*,'("structures.in contains structures of CE-rank ",i2)') CErank
!depracated>if (.not.(CErank==CE%CErank)) then 
!depracated>   print *, "ERROR: CE rank given in lat.in and structures.in does not match!"
!depracated>   STOP
!depracated>end if
!depracated>
!depracated>
!depracated>
!depracated>! Create a table (gettype) of the atom label for any k-nary
!depracated>! This is just the list of values that the occupation variables can have
!depracated>call get_SpinRank_mapping(CErank, gettype, rankSpace)
!depracated>!c=0
!depracated>!do iTyp = -CErank/2,CErank/2,+1 !GH We should replace this with a call to get_spinRank_mapping 
!depracated>!   if (mod(CErank,2)==0 .and. iTyp == 0) cycle
!depracated>!   c=c+1
!depracated>!   gettype(c) = iTyp
!depracated>!enddo
!depracated>!err = .false.  !! GH What did we need this for?
!depracated>
!depracated>
!depracated>!--------------------------------------------------------------------------------
!depracated>! Counting structures
!depracated>!--------------------------------------------------------------------------------
!depracated>! start anew
!depracated>rewind(13)
!depracated>! We need to skip the special options
!depracated>do iL=1, nSpecialOptions
!depracated>  call co_ca(13,err)
!depracated>  read(13,*) line
!depracated>enddo
!depracated>
!depracated>
!depracated>write(*,'(A)') "Counting the structures in structures.in..."
!depracated>
!depracated>iStr = 0
!depracated>do ! Read in iStr-th structure   
!depracated>   call co_ca(13,err)
!depracated>   if (err) exit ! Jump out of the loop if we are at the end of the file
!depracated>   iStr = iStr + 1
!depracated> 
!depracated>   !sascha>   write(*,'(I5,2x)',advance='no') iStr
!depracated>   call read_structure(13,tStr,CE%CErank,CE%nCouplings,peratom,useweights,ioerr,isReference=.false.)
!depracated>   !sascha> print*, "done reading ", tStr%title
!depracated>   if (ioerr/=0) exit
!depracated>
!depracated>   allocate(tStr%reference(nRef))
!depracated>   do iRef=1,nRef
!depracated>     call read_structure(13,tStr%reference(iRef),CE%CErank,CE%nCouplings,peratom,useweights,ioerr,isReference=.true.)
!depracated>   enddo
!depracated>enddo ! couting
!depracated>if (ioerr/=0) then 
!depracated>  iStr = iStr - 1 ! There was junk at the end of the file---don't count
!depracated>                  ! it as another input structure.
!depracated>  print *
!depracated>  write(*,'(A)') "WARNING: structures.in contains junk input at the end of file."
!depracated>  write(*,'(A)') "WARNING:   last structure will not be included."
!depracated>  print *
!depracated>endif
!depracated>nStr=iStr
!depracated>write(*,'("read_input_structures found ",i5," input structures.")') nStr
!depracated>CE%totalstructures = nStr
!depracated>!--------------------------------------------------------------------------------
!depracated>! End Counting structures
!depracated>!--------------------------------------------------------------------------------
!depracated>
!depracated>allocate(str(nStr))
!depracated>
!depracated>!--------------------------------------------------------------------------------
!depracated>! Reading structures
!depracated>!--------------------------------------------------------------------------------
!depracated>! start anew
!depracated>rewind(13)
!depracated>! We need to skip the special options
!depracated>do iL=1, nSpecialOptions
!depracated>  call co_ca(13,err)
!depracated>  read(13,*) line
!depracated>enddo
!depracated>do iStr=1,nStr ! Read in iStr-th structure   
!depracated>   call co_ca(13,err)
!depracated>!LN   write(*, '("Reading in structure #",i4, ":  ",$)') iStr
!depracated>
!depracated>   call read_structure(13,str(iStr),CE%CErank,CE%nCouplings,peratom,useweights,ioerr,isReference=.false.)
!depracated>   if (ioerr/=0) exit
!depracated>
!depracated>   allocate(str(iStr)%reference(nRef))
!depracated>   do iRef=1,nRef
!depracated>     call read_structure(13,str(iStr)%reference(iRef),CE%CErank,CE%nCouplings,peratom,useweights,ioerr,isReference=.true.)
!depracated>!     if (str(iStr)%reference(iRef)%automaticReference) call setup_reference_structure(str(iStr),str(iStr)%reference(iRef),CE)
!depracated>   enddo
!depracated>
!depracated>enddo ! reading structures
!depracated>!--------------------------------------------------------------------------------
!depracated>! End Reading structures
!depracated>!--------------------------------------------------------------------------------

!================================================================================
!END_OF_DEPRACTED_CODE


select case (structFormat)
case('P')
! count how many atom *types* there are
call co_ca(13,err)
do iL = 1,5; read(13,*); call co_ca(13,err); enddo  ! Skip the first 5 lines
read(13,*,iostat=ioerr) cTransform ! transformation wanted?
call ucase(cTransform)
if (cTransform=='T') then
  do iL = 1,4; read(13,*);enddo  ! Skip 4 more lines
else
  backspace(13)
endif
cTransform=""

read(13,'(a1000)') line ! Read in the number of atoms of each type
call parse_line_for_numbers(line," ",maxEntries=CE%CErank*CE%ncouplings,readEntries=c,values=numbers)
CErank = c / CE%ncouplings
write(*,'("structures.in contains structures of CE-rank ",i2)') CErank
if (.not.(CErank==CE%CErank)) then 
   print *, "ERROR: CE rank given in lat.in and structures.in does not match!"
   STOP
end if

! Create a table (gettype) of the atom label for any k-nary
! This is just the list of values that the occupation variables can have
call get_SpinRank_mapping(CErank, gettype, rankSpace)
!c=0
!do iTyp = -CErank/2,CErank/2,+1 !GH We should replace this with a call to get_spinRank_mapping 
!   if (mod(CErank,2)==0 .and. iTyp == 0) cycle
!   c=c+1
!   gettype(c) = iTyp
!enddo
!err = .false.  !! GH What did we need this for?


!--------------------------------------------------------------------------------
! Counting structures
!--------------------------------------------------------------------------------
! start anew
rewind(13)
! We need to skip the special options
do iL=1, nSpecialOptions
  call co_ca(13,err)
  read(13,*) line
enddo


write(*,'(A)') "Counting the structures in structures.in..."

iStr = 0
do ! Read in iStr-th structure   
   call co_ca(13,err)
   if (err) exit ! Jump out of the loop if we are at the end of the file
   iStr = iStr + 1
!lance   write(*,'(I5,2x)',advance='no') iStr
   call read_structure_poscar(13,tStr,CE%CErank,CE%nCouplings,peratom,useweights,ioerr,isReference=.false.)

   if (ioerr/=0) exit

   allocate(tStr%reference(nRef))
   do iRef=1,nRef
     call read_structure_poscar(13,tStr%reference(iRef),CE%CErank,CE%nCouplings,peratom,useweights,ioerr,isReference=.true.)
   enddo
enddo ! couting
if (ioerr/=0) then 
  iStr = iStr - 1 ! There was junk at the end of the file---don't count
                  ! it as another input structure.
  print *
  write(*,'(A)') "WARNING: structures.in contains junk input at the end of file."
  write(*,'(A)') "WARNING:   last structure will not be included."
  print *
endif
nStr=iStr
write(*,'("read_input_structures found ",i3," input structures.")') nStr
CE%totalstructures = nStr
!--------------------------------------------------------------------------------
! End Counting structures
!--------------------------------------------------------------------------------

allocate(str(CE%totalstructures))

!--------------------------------------------------------------------------------
! Reading structures
!--------------------------------------------------------------------------------
! start anew
rewind(13)
! We need to skip the special options
do iL=1, nSpecialOptions
  call co_ca(13,err)
  read(13,*) line
enddo
do iStr=1,CE%totalstructures ! Read in iStr-th structure   
   call co_ca(13,err)

   !debug> write(*, '("Reading in structure #",i4, ":  ",$)') iStr

   call read_structure_poscar(13,str(iStr),CE%CErank,CE%nCouplings,peratom,useweights,ioerr,isReference=.false.)
   if (ioerr/=0) exit

   allocate(str(iStr)%reference(nRef))
   do iRef=1,nRef
     call read_structure_poscar(13,str(iStr)%reference(iRef),CE%CErank,CE%nCouplings,peratom,useweights,ioerr,isReference=.true.)
!     if (str(iStr)%reference(iRef)%automaticReference) call setup_reference_structure(str(iStr),str(iStr)%reference(iRef),CE)
   enddo

enddo ! reading structures
!--------------------------------------------------------------------------------
! End Reading structures
!--------------------------------------------------------------------------------

case('G') !this means you are using DFTB+

call get_SpinRank_mapping(CE%CErank, gettype, rankSpace)



write(*,'(A)') "Counting the structures in structures.in..."

iStr = 0
do ! Read in iStr-th structure                                                                                                                     
  call co_ca(13,err)
   if (err) exit ! Jump out of the loop if we are at the end of the file 
   iStr = iStr + 1
   write(*,'(I5,2x)',advance='no') iStr
   call read_structure_gen(13,tStr,CE,peratom,useweights,ioerr,isReference=.false.)
   if (ioerr/=0) exit

   allocate(tStr%reference(nRef))
   do iRef=1,nRef
     call read_structure_gen(13,tStr%reference(iRef),CE,peratom,useweights,ioerr,isReference=.true.)
   enddo
enddo ! counting           

!now iStr = number of structures in structures.in
CE%totalstructures = iStr
write(*,'("read_input_structures found ",i3," input structures.")') CE%totalstructures
allocate(str(CE%totalstructures))

!debug>print*, iStr, " = ", CE%nStructures, " = ", size(str(:))
                                                                                                                 
rewind(13) 

!We need to skip the special options
do iL=1, nSpecialOptions
  call co_ca(13,err)
  read(13,*) line
enddo

do iStr=1,CE%totalstructures ! Read in iStr-th structure                                                                                                                                                  
   call co_ca(13,err)
   !debug> write(*, '("Reading in structure #",i4, ":  ",$)') iStr

   call read_structure_gen(13,str(iStr),CE,peratom,useweights,ioerr,isReference=.false.)
   if (ioerr/=0) exit

   allocate(str(iStr)%reference(nRef))
   do iRef=1,nRef
     call read_structure_gen(13,str(iStr)%reference(iRef),CE,peratom,useweights,ioerr,isReference=.true.)
   enddo
enddo ! reading structures   

!debug>do iStr=1,CE%nStructures ! Read in iStr-th structure                 
!debug>   print*, str(iStr)%nAtoms
!debug>   print*, "x =", str(iStr)%pos(1,:), " y=", str(iStr)%pos(2,:), " z=",str(iStr)%pos(3,:)
!debug>enddo

case default
   print*, "failsave: no valid structFormat"

endselect

close(13)

!--------------------------------------------------------------------------------
! For adsorbates: get single energies
if (CE%adsorbate%isAds) then
  call read_adsorbateData(CE%CErank,CE%adsorbate%singleEnergy,CE%adsorbate%freeSiteRank)
endif

ENDSUBROUTINE read_input_structures



!================================================================================
! Read a single structure (str) from file (fh)
! tk 02/05/2012: changed the interface slightly. We don't need CE% in this routine; CErank and nCouplings are sufficient.
!                In this way, I can also use this subroutine in different contexts without the whole CE% available.
subroutine read_structure_poscar(fh,str,CErank,nCouplings,peratom,useweights,ioerr,isReference,justReadStructureDetails_)
integer, intent(in)  :: fh ! filehandle
type(crystal), intent(out) :: str
integer, intent(in) :: CErank, nCouplings
logical, intent(in)  :: peratom, useweights
integer, intent(out) :: ioerr
logical, intent(in)  :: isReference
logical, intent(in), optional :: justReadStructureDetails_
logical :: justReadStructureDetails
logical              :: automaticReference
integer iStr, iLV, iL, iAt, c, ic, i
integer indx
integer nStr, nAt, nSpin, nPos
integer, allocatable :: nTyp(:,:)
character(1) :: aBasCoord
character(80) :: line
logical err

character(1) :: cTransform
logical :: doTransform
real(dp) :: TransformMatrix(3,3), TransformShift(3)

if (present(justReadStructureDetails_)) then; justReadStructureDetails = justReadStructureDetails_; else; justReadStructureDetails=.false.; endif

if (.not. isReference) then
! This is NO REFERENCE structure
  call co_ca(fh,err)
  read(fh,'(a80)',iostat=ioerr) str%title  !Read in title
  !LNwrite(*, '(A50)') adjustl(str%title)
  call co_ca(fh,err)
  read(fh,*) str%a0  ! Read in lattice constant
else
! This IS a REFERENCE structure
  call co_ca(fh,err)
  read(fh,*) line
  call ucase(line)
  if (line(1:1) == 'A') then
    automaticReference = .true.
    str%automaticReference = .true.
  else
    automaticReference = .false.
    str%automaticReference = .false.
    if (line(1:1) /= 'R') then
      print *
      print *, "ERROR: Expected a reference structure following the input structure"
      print *, "ERROR: Reference structure section must start with 'R'"
      stop
    endif
  endif
endif
if (.not. isReference .or. (isReference .and. .not. automaticReference)) then
! If the structure is a REAL structure (and no reference structure)
! OR
! if it is a REFERENCE structure W/O AUTOMATIC determination of the reference structure:
! ==> read inputs
  do iLV = 1, 3 ! Read in lattice vectors
    call co_ca(fh,err)
    read(fh,*,iostat=ioerr) str%LV(:,iLV)
    write(14,'(3(f8.4,1x))') str%LV(:,iLV)
  enddo
  
  
  allocate(str%nTyp(CErank,nCouplings))  ! number of atoms of each type
  allocate(str%xnTyp(CErank,nCouplings)) ! number of atoms of each type: only x-relevant atoms are counted (x=concentration)
           str%xnTyp = 0         ! init
  allocate(str%x(CErank,nCouplings))     ! concentration of each type
     
  ! do we need to transform the lattice vectors (rotate, shift origin?)
  call co_ca(fh,err)
  read(fh,*,iostat=ioerr) cTransform ! transformation wanted?
  call ucase(cTransform)
  if (cTransform=='T') then
    doTransform=.true.               ! yes
    do i=1,3
      call co_ca(fh,err)
      read(fh,*,iostat=ioerr) TransformMatrix(:,i)  ! read transformation matrix
    enddo
    call co_ca(fh,err)
    read(fh,*,iostat=ioerr) TransformShift(:)       ! read transformation (cartesian) shift of atomic positions
                                                    ! i.e. re-set the origin
    str%LV = matmul(TransformMatrix,str%LV)   ! transform the lattice vectors
  else                               ! no transformation
    doTransform=.false.
    backspace(fh)
  endif
  
  
  call co_ca(fh,err)
  read(fh,*) str%nTyp(:,:) ! Read in number of each type of atom for each CE
  write(14,'("Atoms: ",20(i3,1x))') str%nTyp

  
  if (any(str%nTyp<0)) then
    str%atomsMissing = .true.
    if (count(str%nTyp<0) >1 ) stop "ERROR: At most ONE missing atom sort allowed"
    if (str%nTyp(CErank,nCouplings) >= 0) stop "ERROR: Up to now, only the last type can by missing"

    str%missingTyp   = minloc(str%nTyp)
    str%missingN     = abs(str%nTyp( str%missingTyp(1), str%missingTyp(2) ))
    str%nTyp         = abs(str%nTyp)
  else
    str%atomsMissing = .false.
    str%missingN     = 0
  endif

  call co_ca(fh,err)
  read(fh,*) aBasCoord ! Read in coordinate system (Cartesian or direct)
  call ucase(aBasCoord)
  str%aBasCoord = aBasCoord
  nAt = sum(str%nTyp) ! sum over all types of atoms of all CEs = total number of atoms
  str%nAtoms = nAt
  write(14,'("Number of atoms: ",i3)') str%nAtoms

  
  ! Determine concentration (legacy code, should be removed someday)
  !________________________________________________________________________________
str%conc = -10000 
  ! Remark: the general concentration str%x(:) is set in check_input_structures.


  ! Read atoms
  !________________________________________________________________________________
  allocate(str%pos(3,nAt),str%dBasCoord(3,nAt))
  do iAt = 1, nAt-str%missingN ! Read in each atom
    call co_ca(fh,err)
    read(fh,*,iostat=ioerr) str%pos(:,iAt)
    if (aBasCoord=='D') then !Convert to lattice coordinates
       str%dBasCoord(:,iAt) = str%pos(:,iAt)
       str%pos(:,iAt) = matmul(str%LV,str%pos(:,iAt))
       if (doTransform) str%pos(:,iAt) = str%pos(:,iAt) + TransformShift ! if wanted: transform the atomic
                                                                        ! positions (reset origin)
    else ! cartesian coordinates
       if (doTransform) str%pos(:,iAt) = matmul(TransformMatrix,str%pos(:,iAt)) + TransformShift
                                                              ! if wanted: transform by rotating and shifting 
                                                              ! the atomic position (note: in the direct coordinates the
                                                              ! rotation was already included in the lattice vectors, here
                                                              ! in cartesian coordinates we have to rotate the atoms
                                                              ! separately.
    endif
    if (ioerr/=0) then
      write(*,'("Missing atoms in structure ",i4," in ''read_input_structures''")') iStr
      stop
    endif
  enddo
  
  !________________________________________________________________________________
  
  
  if (.not. justReadStructureDetails) then ! if we're only intereseted in the structure and not in the energy, skip this 
  if (.not. isReference) then ! only read the following information if we are not merely
                              ! reading a reference structure (the energy of which will be
                              ! determined by the code through the use of the reference ECI.
    call co_ca(fh,err)
    read(fh,*) str%energy
    !lance >
!    str%energy = str%energy * 1000
    if (.not. peratom) then ! total energies are given so normalize per atom
      str%energy = str%energy/nAt
    endif
       
    if (useweights) then
      call co_ca(fh,err)
      read(fh,*) str%weight
    else
      str%weight=1
    endif
  endif
  endif

  ! Determine atom types and spins for the structure:
  ! set failsafe values, that'll break the code. The correct values are set in check_input_structures.
  allocate(str%ATyp(nAt))
  allocate(str%spin(nAt))
  str%atyp = -100
  str%spin = -100

endif ! no automatic reference

str%isHelpStructure = .false.   ! whenever we actually read in a structure or build a reference structure,
                                ! the structure is "real" and not only used for helping some other stuff (e.g. 
                                ! convex hull)

end subroutine read_structure_poscar

subroutine read_structure_gen(fh,str,CE,peratom,useweights,ioerr,isReference)
integer, intent(in)  :: fh ! filehandle
type(crystal), intent(out) :: str
type(CE_variables), intent(inout) :: CE
logical, intent(in)  :: peratom, useweights
integer, intent(out) :: ioerr
logical, intent(in)  :: isReference
logical              :: automaticReference
integer iStr, iLV, iL, iAt, c, ic, i, iObs
integer CErank, nC, indx
integer nStr, nAt, nSpin, nPos
integer, allocatable :: nTyp(:,:)
character(1) :: aBasCoord, dummyChar
integer :: dummyInt
character(80) :: line
logical err
integer, allocatable :: atomCounter(:)

character(1) :: cTransform
logical :: doTransform
real(dp) :: TransformMatrix(3,3), TransformShift(3)

CErank = CE%CErank
nC     = CE%ncouplings



call co_ca(fh,err)
read(fh,'(a80)',iostat=ioerr) str%title  !Read in title
write(*, '(A50)') adjustl(str%title)
call co_ca(fh,err)
read(fh,*) str%nAtoms, dummyChar
!print*, str%nAtoms, " ", dummyChar
call co_ca(fh,err)
call ucase(dummyChar)
if(dummyChar/='S') stop "failsave: supercell format needed"

write(14,'("Number of atoms: ",i3)') str%nAtoms

allocate(str%nTyp(CErank,nC))  ! number of atoms of each type
allocate(str%xnTyp(CErank,nC)) ! number of atoms of each type: only x-relevant atoms are counted (x=concentration)
         str%xnTyp = 0         ! init
allocate(str%x(CErank,nC))     ! concentration of each type
allocate(atomCounter(CErank))  ! in gen format we will need to count atoms of each type manually     
allocate(str%pos(3,str%nAtoms))

!debug> print*, size(str%xntyp(:,1)), "  ", CErank, " = ", CE%CErank

atomCounter(:)=0

read(fh,'(a50)') line
call co_ca(fh,err)

do iAt = 1, str%nAtoms
   call co_ca(fh,err)
   read(fh,*,iostat=ioerr) dummyInt, dummyInt, str%pos(:,iAt)
   atomCounter(dummyInt)= atomCounter(dummyInt) + 1

   !we need a failsave here to check whether the structure has the appropriate number of atoms!
   !-sascha 31 May 2012
enddo

str%nTyp(:,1)= atomCounter(:)

!read out null line of .gen format
call co_ca(fh,err)
read(fh,*,iostat=ioerr) str%LV(:,1)

!read in lattice vectors
do iLV = 1, 3 
    read(fh,*,iostat=ioerr) str%LV(:,iLV)
    write(14,'(3(f8.4,1x))') str%LV(:,iLV)
enddo

call co_ca(fh,err)
read(fh,*) str%energy
if (.not. peratom) then ! total energies are given so normalize per atom
   str%energy = str%energy/str%nAtoms 
endif

print*, "nAt=", nAt, str%nAtoms



if (useweights) then
   call co_ca(fh,err)
   read(fh,*) str%weight
else
   str%weight=1
endif
!endif

!Determine atom types and spins for the structure:
!set failsafe values, that'll break the code. The correct values are set in check_input_structures.
allocate(str%ATyp(str%nAtoms))
allocate(str%spin(str%nAtoms))
str%atyp = -100
str%spin = -100

str%conc = -10000 ! let's see what happens....

str%isHelpStructure = .false.   ! whenever we actually read in a structure or build a reference structure,
                                ! the structure is "real" and not only used for helping some other stuff (e.g. 
                                ! convex hull)
end subroutine read_structure_gen


!****************************************************************************************************

subroutine read_adsorbateData(k,energy,free)
integer, intent(in) :: k
real(dp), pointer   :: energy(:)
integer, intent(out):: free

integer i, ioerr
logical err

write(*,'(A)') "Reading adsorbate data..."

open(13,file="adsorbate.in",status="old",iostat=ioerr)
if (ioerr /= 0) stop "ERROR: opening file adsorbate.in"

allocate(energy(k))

call co_ca(13,err)
read(13,*) free

do i=1,k
  call co_ca(13,err)
  read(13,*) energy(i)
enddo

close(13)

if (energy(free)/=0._dp) then
  print*
  write(*,'(A)') "WARNING: The 'non'-adsorbate (vacancy) should have energy 0"
  print *
endif


end subroutine read_adsorbateData





!***************************************************************************************************
SUBROUTINE write_figure_reps(fig,EqvfigFile)
type(figRep), pointer :: fig(:)
character(30), intent(in) :: EqvfigFile
integer nF, nV
integer iFg, iV, iRep

open(12,file=EqvFigFile)
write(12,*) "# List of inequivalent figures and all corresponding"
write(12,*) "# symmetrically equivalent figures"
write(12,*) "#"

Nf = size(fig)
do iFg = 0, nF-1
   write(12,'("# Figure No. : ",i5," Average dist.: ",f10.6, " Damping: ",f10.6)') iFg, fig(iFg)%Rep(1)%avgR, fig(iFg)%damping
   write(12,'("# Degeneracy of figure: ",i5)') fig(iFg)%count
   nV = fig(iFg)%Rep(1)%nV
   do iRep = 1, fig(iFg)%count
      write(12,'("# Number: ",i2)') iRep
      do iV = 1, nV
         write(12,'(3(f10.6,1x),i5,2x,i5)') fig(iFg)%Rep(iRep)%vertex(:,iV), fig(iFg)%Rep(iRep)%label(iV) , fig(iFg)%Rep(iRep)%s(iV)
      enddo
   enddo
   write(12,'("# End of figure number: ",i5,/)') iFg
enddo
close(12)
ENDSUBROUTINE write_figure_reps
!**************************************************************************************************
SUBROUTINE write_clusters(cluster,fh,offset)
type(figure), intent(in) :: cluster(:)
integer, intent(in) :: fh
integer, intent(in) :: offset

integer nCl, iCl, nvert, ivert


if (offset==0) then ! if this is true that's the first time this routine is called
  write(fh,'(A)') "#--------------------------------------------------------------------------------"
  write(fh,'(A)') "# Cluster number: 0"
  write(fh,'(A)') "# This is the empty cluster = constant term for the cluster expansion."
  write(fh,'(A)') "#"
  write(fh,'(A)') "# - The empty cluster is not explicitely listed in this file. "
  write(fh,'(A)') "#   (This comment here only reminds the user of the fact that"
  write(fh,'(A)') "#   the empty cluster does exist in UNCLE.)"
  write(fh,'(A)') "# - The emtpy cluster does not contain any vertices." 
  write(fh,'(A)') "# - The empty cluster is not read in from this file."
  write(fh,'(A)') "# - The emtpy cluster is always used in the fitting procedure."
  write(fh,'(A)') "# - The emtpy cluster counts as a separate cluster for the fitting procedure."
  write(fh,'(A)') "#   (E.g., a simple fit with 1 cluster is a fit with the empty cluster only.)"
end if

nCl = size(cluster)
nvert = cluster(1)%nV   ! we call this routine for every iV, so nvert remains fixed
                    ! for all figures that are written out by this routine
!open(fh,file="figures.out")
do iCl = 1,nCl
   write(fh,'(A)') "#--------------------------------------------------------------------------------"
   write(fh,'("# Cluster number (total/of this kind):",i6," / ",i6,/,i6)') iCl+offset, iCl, iCl
   write(fh,'("# Number of vertices:",/,i4)') nvert
   write(fh,'("# Avg. Distance:",/,f10.6,/,"# Vertices: (x,y,z) | d-vector label | s-vector")') cluster(iCl)%avgR
   do ivert = 1, nvert
      write(fh,'(3(f10.6,2x),5x,i3,5x,i3)') cluster(iCl)%vertex(:,ivert), cluster(iCl)%label(ivert) , cluster(iCl)%s(ivert) 
   enddo
   write(fh,'("#")')
enddo 
!close(fh)
ENDSUBROUTINE write_clusters


SUBROUTINE read_CSvars(CS) ! title,cetype,CErank,LV,aBas,aBasCoord,Ninterior,cutoff,eps)

type(CS_variables), intent(inout) :: CS

integer fh ! filehandle

character(1) ctemp,ctemp2
character(80) csvars_file
logical :: err
integer ioerr

fh = 10
csvars_file = 'CS.in'
open(fh,file=csvars_file,status='old',iostat=ioerr)
if (ioerr /= 0) then
  print *, "ERROR: Opening CS file ", trim(csvars_file)
  stop
endif
!if (pr) write(*,'(A,A)') "Reading lattice definition from ", adjustl(trim(CE%lat_file))
!print *

!-------------- title ---------------------
! Read in the mode ( Add new structures or perform random subsets)
call co_ca(fh, err)
read(fh,*) CS%mode
call ucase(CS%mode)

! Do you want to employ the weighted l-1 norm scheme?
call co_ca(fh, err)
read(fh,*) CS%pltone
call ucase(CS%pltone)

! Loop over values of sigma^2 or just let the hard-coded expression
! for it be used.  I've had excellent luck with using the hard-coded
! expression, but I've kept this functionality in the code in case 
! someone wants to toy around with it.
call co_ca(fh, err)
read(fh,*) CS%loop
call ucase(CS%loop)

! If  you are going to loop over values of sigma^2 then tell me the
! starting point, the ending point and the step size/type.
call co_ca(fh, err)
read(fh,*) CS%musigma_start, CS%musigma_finish

call co_ca(fh,err)
read(fh,*) ctemp, CS%delta_musigma

if (ctemp == '*') then
   CS%musigma_change_mode = 0
elseif (ctemp == '+') then
   CS%musigma_change_mode = 1
else
   CS%musigma_change_mode = -1
endif

call co_ca(fh, err)
read(fh,*) CS%nInpStruct

call co_ca(fh, err)
read(fh,*) CS%nSubsets

call co_ca(fh, err)
read(fh,*) CS%J_cutoff

call co_ca(fh, err)
read(fh,*) CS%penaltyfxn

END SUBROUTINE read_CSvars



SUBROUTINE write_structures_in(structures)

type(crystal) structures(:)

integer ioerr,i,j
integer one
open(7,file='structures.in.toy',iostat=ioerr)

one = 1
!write(7,'(A)') 'scalar'
write(7,'(A)') 'peratom'
write(7,'(A)') 'noweights'
write(7,'(A)') 'poscar'

do i = 1 , size(structures,1)
   write(7,'(A)') '#---------------------------'
   write(7,'(A)') structures(i)%title
   write(7,'(I1)') one
   write(7,'(3F8.3)') structures(i)%LV(:,1)
   write(7,'(3F8.3)') structures(i)%LV(:,2)
   write(7,'(3F8.3)') structures(i)%LV(:,3)
   write(7,'(5I3)') structures(i)%nTyp
   write(7,'(A)') structures(i)%aBasCoord
   if (structures(i)%aBasCoord == 'C') then
      do j = 1, structures(i)%nAtoms
         write(7,'(3F10.6)') structures(i)%pos(:,j)!structures(i)%dBasCoord(:,j)
      end do
   else
      do j = 1, structures(i)%nAtoms
         write(7,'(3F10.6)') structures(i)%dBasCoord(:,j)
      end do
   end if
   write(7,'(A)') '#Energy:'
   write(7,'(F10.5)') structures(i)%energy

end do
close(7)
END SUBROUTINE write_structures_in

SUBROUTINE read_toy_js_and_noise(noise_level,toy_Js,n_ToyCorrs)

!type(CS_variables),intent(inout):: CS
logical :: err
real(dp), pointer:: toy_js(:)
real(dp),intent(inout):: noise_level
integer, intent(in) :: n_ToyCorrs
integer fh, which_J
real(dp) J_value
character(80) toyjs_file
integer ioerr,ios


fh = 10
toyjs_file = 'toyJs.in'

open(fh,file=toyjs_file,status='old',iostat=ioerr)
if (ioerr /= 0) then
  print *, "ERROR: Opening toyJs file ", trim(toyjs_file)
  stop
endif
allocate(toy_Js(n_ToyCorrs))
toy_Js = 0

call co_ca(fh, err)
read(fh,*,iostat=ios) noise_level
call co_ca(fh, err)
do while (.true.)
   read(fh,*,iostat=ios) which_J, J_value
   if (ios /= 0) exit
   toy_Js(which_J) = J_value
enddo


END SUBROUTINE read_toy_js_and_noise
!***************************************************************************************************
recursive SUBROUTINE read_lattdef(CE,MPI_rank,iReference) ! title,cetype,CErank,LV,aBas,aBasCoord,Ninterior,cutoff,eps)
use utilities_module, only: ralloc

type(CE_variables), intent(inout) :: CE
integer, intent(in)               :: MPI_rank
integer, intent(in), optional     :: iReference

logical                           :: isReference
logical pr ! print statments?
integer fh ! filehandle

integer ivec, ipt, ipt1, ipt2, iRef, icoord
logical :: err
real(dp) :: cart2latt(3,3)

integer ioerr
integer iBulkAtoms
integer iTopAtom, iBottomAtom

character(1000) :: line
integer, pointer :: SpinSpace(:), RankSpace(:)

if (present(iReference)) then; isReference=.true.; else; isReference=.false.; endif
if (MPI_rank == 0) then; pr=.true.; else; pr=.false.; endif

fh = 10
if (isReference) fh=fh+iReference
print *, CE%lat_file, "lat file"
open(fh,file=CE%lat_file,status='old',iostat=ioerr)
if (ioerr /= 0) then
  print *, "ERROR: Opening lattice file ", trim(CE%lat_file)
  stop
endif
if (pr) write(*,'(A,A)') "Reading lattice definition from ", adjustl(trim(CE%lat_file))
print *

!-------------- title ---------------------
call co_ca(fh, err)
read(fh,*) CE%title

!-------------- CE type  ------------------
call co_ca(fh, err)
read(fh,*) CE%cetype ! Surface or bulk or adsorbate?
call ucase(CE%cetype)
CE%isBulk          = .false.
CE%surface%isSurf  = .false.
CE%adsorbate%isAds = .false.
CE%ncouplings      = 1
!----- surface
if (CE%cetype=='S') then ! this is a surface CE
  CE%surface%isSurf=.true.
!----- adsorbate
elseif (CE%cetype=='A') then
  CE%CEtype = 'S'            ! will be used for the enumeration, the fact that it is an adsorbate CE (i.e. coupled
                             ! is taken care of by the next variables)
  CE%surface%isSurf =.true.  ! an adsorbate problem is also a surface problem
  CE%adsorbate%isAds=.true.
  CE%ncouplings     =2       ! the CE will be 2 coupled CEs (IMPORTABERT: the rank of each CE must be
                             ! EQUAL at the moment...)
!----- bulk
elseif (CE%cetype=='B') then
  CE%isBulk = .true.
!----- don't know
else
  stop "ERROR: type of CE undetermined. Must be 'B' or 'S' or 'A'"
endif

!-------------- CE rank -------------------
call co_ca(fh, err)
read(fh,*) CE%CErank
call co_ca(fh, err)
read(fh,*) CE%nObservables
!-------------- epsilon -------------------
call co_ca(fh, err)
read(fh,*) CE%eps
!-------------- lattice vectors -----------
call co_ca(fh, err)
call read_pLV(fh,CE%pLV)
!-------------- basis atoms ---------------
call co_ca(fh, err)
read(fh,*) CE%nBas

!!----- surface
!if (CE%surface%isSurf) then
!  call get_SpinRank_mapping(CE%CErank, SpinSpace, RankSpace)
!  allocate(CE%surface%BasType(CE%nBas))
!  allocate(CE%surface%BulkSpin(CE%nBas))
!  allocate(CE%surface%BulkEnumRank(CE%nBas))
!  allocate(CE%surface%BulkList(CE%nBas))
!  CE%BasEnum = .false.  ! the correct enumeration will be determined later
!!----- bulk
!else
!  CE%BasEnum = .true.   ! will enumerate all dset members
!endif

call co_ca(fh, err)
read(fh,*) CE%CoordSys
call ucase(CE%CoordSys)


!________________________________________________________________________________
! Read in: dset
call read_dset(fh,CE%CErank,CE%nBas,CE%pBas,CE%digit,CE%label,CE%xRelevant,CE%fixedOcc,CE%equivalencies)
call transform_dset(CE%CoordSys,CE%pLV,CE%pBas,CE%eps)   ! e.g. direct -> cartesian

! Make sure that all basis atoms are at different locations
call check_for_collisions(CE%pBas)

! Preliminary CEbas setup 
allocate(CE%CEbas(CE%nBas))
CE%CEbas = 1   ! a priori all d-vectors belong to the 1st CE
!---------- end of basis atoms read ----------

!----- surface
if (CE%surface%isSurf) then ! surface CE?

  if (pr) then
    write(*,*)
    write(*,*) "This is a surface CE"
  endif

  call read_surface_normal(fh,CE%pLV,CE%surface%unitNormal)
  call find_topbottom(CE%pLV,CE%pBas,CE%surface%unitNormal,CE%surface%distT,CE%surface%distB)

  if (pr) then
    write(*,'(1x,A,3(F10.3,1x))') "* Surface normal:", CE%surface%unitNormal
    write(*,*) "* Found surfaces at distances from origin:"
    write(*,'(A4,3(A15,A2))')   "|","Bottom Surface", "|", "Top Surface",   "|"
    write(*,'(A4,3(F15.3,A2))') "|",CE%surface%distB, "|", CE%surface%distT,"|"
    print *
  endif

  ! if adsorbate CE: read which dset members are adsorbate sites
  allocate(CE%adsorbate%isAdsd(CE%nBas))
  if (CE%adsorbate%isAds) then 
    call read_adsorbate_sites(fh,CE%nBas,CE%adsorbate%isAdsd)
    where (CE%adsorbate%isAdsd)
      CE%CEbas = 2 ! the adsorbate d-vectors belong to CE 2
                   ! (maybe we could use this variable ONLY, but for the time being I consider it
                   ! more helpful to have the isAdsd too...)
    end where
  else
    CE%adsorbate%isAdsd = .false.
    CE%CEbas            = 1
  endif

endif ! surface CE
!----- bulk
! do nothing

call co_ca(fh, err)
read(fh,*) CE%cutoff

call co_ca(fh,err)
read(fh,*) CE%nreferences

if (isReference .and. CE%nreferences>0) then
  if (pr) write(*,'(/,A)') "WARNING: References only allowed in the master lat.in."
  if (pr) write(*,'(A)')   "WARNING: Rest of lat.in will be ignored."
  CE%nreferences = 0
  allocate(CE%reference(CE%nreferences))
  allocate(CE%l2r(CE%nreferences))
  return ! <------------ return ----------------
endif

allocate(CE%reference(CE%nreferences))
allocate(CE%l2r(CE%nreferences))

do iRef=1,CE%nreferences
  if (pr) then
    write(*,'(A)')    ""
    write(*,'(A)')    "--------------------------------------------------------------------------------"
    write(*,'(A,I5)') "Setting up reference ", iRef
  endif
  call co_ca(fh,err)
  read(fh,*) CE%reference(iRef)%lat_file
  call co_ca(fh,err)
  read(fh,*) CE%reference(iRef)%cluster_file
  call co_ca(fh,err)
  read(fh,*) CE%reference(iRef)%J_file

! Not yet implemented. As for now I use adsorbate.in for the relevant information.  
!  ! in adsorbate mode we need the structures.in of the surface reference
!  if (CE%adsorbate%isAds .and. CE%reference(iRef)%surface%isSurf) then
!    call co_ca(fh,err)
!    read(fh,*) CE%reference(iRef)%str_file
!  endif

  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mlv(:,1)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mlv(:,2)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mlv(:,3)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mds(:,1)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mds(:,2)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mds(:,3)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Ods(:)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%nDel; allocate(CE%l2r(iRef)%dsetDel(CE%l2r(iRef)%nDel))
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%dsetDel(:)
  call read_lattdef(CE%reference(iRef),MPI_rank,iRef)
  if (pr) write(*,'(A)')    "--------------------------------------------------------------------------------"
enddo

!if (CE%surface%isSurf) then
!  if (CE%nReferences == 0) print *, "WARNING: consider using a bulk reference for surface CE"
!  if (CE%nReferences > 1) stop "ERROR: for surface CE only 1 (bulk) reference can be used now"
!  if (any(CE%reference(:)%CEtype /= 'B')) stop "ERROR: for surface CE only bulk references are implemented"
!endif
  
!if (CE%adsorbate%isAds .and. .not. any(CE%reference(:)%surface%isSurf)) & ! adsorbate CE but no surface reference
!     stop "ERROR: an (bi-binary) adsorbate CE needs a surface reference!"

!call co_ca(fh,err)
!read(fh,*) CE%GSSspan

!REMARK: This should be read in someday, but for now I always want this to be true here
!CE%GSS_full = .true.


close(fh)

contains

function ProjectionLength(vec,direction)
  real(dp)  :: ProjectionLength
  real(dp)  :: vec(3)
  real(dp)  :: direction(3)
  
  direction = direction / norm(direction)
  ProjectionLength=abs(dot_product(vec,direction))
  
end function ProjectionLength

subroutine read_pLV(fh,pLV)
  integer, intent(in)   :: fh
  real(dp), intent(out) :: pLV(3,3)
  integer i
  logical err
  
  do i=1,3
    call co_ca(fh, err)
    read(fh,*) pLV(:,i)
  enddo
end subroutine read_pLV



subroutine read_dset(fh,k,nD,d,digit,label,xrel,fixed,eq)
  integer, intent(in) :: fh          ! filehandle
  integer, intent(in) :: k           ! rank
  integer, intent(in) :: nD          ! number of dset points
  real(dp), pointer   :: d(:,:)      ! intent(out): { coord, dset# } dset points
  integer, pointer    :: label(:,:)  ! intent(OUT): { label#, dset# } labels for the enumeration
  integer, pointer    :: digit(:)    ! intent(OUT): { dset# } number of labels for each point
  logical, pointer    :: xrel(:)     ! intent(OUT): { dset# } relevant for concentration?
  integer, pointer    :: fixed(:)    ! intent(OUT): { fixed# } the index of the dvectors that have fixed occ
  integer, pointer    :: eq(:)       ! intent(OUT): { dset# } equivalency list for enumeration

  integer, pointer :: xrelTmp(:)
  integer :: nE
  real(dp), pointer :: numbers(:)

  character(1000) :: line
  logical err
  integer iD, icoord, i
  
  allocate(d(3,nD),label(k,nD),digit(nD),eq(nD))

!--------------------------------------------------------------------------------
! loop over dset  
  do iD=1,nD
    
    call co_ca(fh, err)
    
    read(fh,'(a1000)') line
    call parse_line_for_numbers(line," ",maxEntries=3,readEntries=nE,values=numbers)
    if (nE /= 3) stop "ERROR: parsing line with d-vector coordinates in lat.in"
    do icoord=1,3
      d(icoord,iD) = numbers(icoord)
    enddo
    call parse_line_for_numbers(line,"/",maxEntries=k,readEntries=nE,values=numbers)
    if (nE==0) then ! only the d-vector coordinates were given
       digit(iD) = k 
       label(:,iD) = (/ (i,i=1,k) /) -1
    else ! occupation was supplied too
       do i = 1, k ! Loop over the number of (possible) labels, exit when there are no more /'s
         label(i,iD) = nint(numbers(i))
   
         ! Sanity check on the input for the label (perhaps not sufficient but catches some errors)
         if (label(i,iD) > k-1 .or. label(i,iD) < 0) then
            write(*,'("Incorrect number for a label, ''",i2,"'', on d-vector #",i2)') label(i,iD), iD
            stop
         endif
         if (i==nE) exit
       enddo
       digit(iD) = i ! Store the number of labels that were specified for each d-vector
    endif
  enddo ! iD
!--------------------------------------------------------------------------------

  ! determine which dset members are relevant for concentration
  allocate(xrel(nD)); xrel=.true.
  where (digit==1)
    xrel=.false.
  end where
  ! ... and write the indices of those dvectors that have fixed occupacy:
  allocate(fixed(count(digit==1)))
  fixed = pack( (/ (iD,iD=1,nD) /) ,digit==1 )

  ! read equivalencies (if present, as denoted by "E")
  call co_ca(fh,err)
  read(fh,'(A1000)') line
  if (line(1:1)=='E') then
    call co_ca(fh,err)
    read(fh,*) eq(:)
  else
    eq = (/(i,i=1,nD)/)
    backspace(fh)
  end if

end subroutine read_dset

subroutine transform_dset(coordsys,LV,d,eps)
  character(1), intent(in) :: coordsys
  real(dp), intent(in) :: LV(3,3)
  real(dp), intent(inout) :: d(:,:)
  real(dp), intent(in) :: eps
  
  real(dp) :: LVinv(3,3)
  logical  :: err
  integer ipt
  
  ! convert direct coordinates to cartesian coordinates
  if (CoordSys=='D') d = matmul(LV,d)
 
  ! If the basis atoms are not inside unit cell 0, then shift by lattice vectors until they are
  call matrix_inverse(LV, LVinv, err)
  do ipt = 1, size(d,2)
    call bring_into_cell(d(:,ipt), LVinv, LV, eps)
  enddo
  
end subroutine transform_dset

subroutine check_for_collisions(d)
  real(dp), intent(in) :: d(:,:)
  integer ipt1, ipt2, nD
  
  nD=size(d,2)

  do ipt1=1, nD
    do ipt2=1, nD
      if (ipt1==ipt2) cycle
      if (norm(CE%pBas(:,ipt1)-CE%pBas(:,ipt2))==0) stop "Error: Basis atoms collision!"
    enddo
  enddo
end subroutine check_for_collisions

subroutine read_surface_normal(fh,LV,normal)
  integer, intent(in)   :: fh         ! filehandle
  real(dp), intent(in)  :: LV(3,3)    ! lattice vectors (parent)
  real(dp), intent(out) :: normal(3)  ! normal
  character(1) :: coordsys
  logical :: err
  
  call co_ca(fh,err)
  read(fh,*) coordsys
  call ucase(coordsys)
  
  call co_ca(fh,err)
  read(fh,*) normal
  
  if (coordsys=='D') normal=matmul(LV,normal)

  normal = normal / norm(normal)   ! make unitnormal

  if (.not. all(normal == (/1,0,0/))) stop "ERROR: Surface normal MUST BE 1,0,0"

end subroutine read_surface_normal

subroutine read_topbottom(fh,d,n,distT,distB)
  integer, intent(in)   :: fh             ! filehandle
  real(dp), intent(in)   :: d(:,:)        ! dset
  real(dp), intent(in)  :: n(3)           ! normal
  real(dp), intent(out) :: distT, distB   ! distances to top and bottom layer, resp.

  logical :: err
  integer :: iTop, iBottom

  call co_ca(fh,err)
  read(fh,*) iTop, iBottom

  distT = ProjectionLength(d(:,iTop)   , n)
  distB = ProjectionLength(d(:,iBottom), n)

  
end subroutine read_topbottom


! Find the top and the bottom surface distance
! by Sascha Maisel
subroutine find_topbottom(pLV,d,n,SdistT,SdistB)
real(dp), intent(in) :: pLV(3,3)  ! parent lattice vectors
real(dp), intent(in) :: d(:,:)    ! d-vectors
real(dp), intent(in) :: n(3)      ! surface normal (unit normal!)
real(dp), intent(out) :: SdistT, SdistB  ! top and bottom distance of the surfaces

integer :: iD,nD
real(dp) :: a(3), b(3)
real(dp) :: M(3,3) ! tranformation matrix to planar basis
real(dp), allocatable :: dslash(:,:) !transformed coordinates
real(dp), allocatable :: tmp(:),trueDistance(:), internalDistance(:)
real(dp) :: tmp1, tmp2
integer, pointer :: list(:)
integer :: gapParameter

nD = size(d,2)

tmp1 = 0._dp
tmp2 = 0._dp
gapParameter = nD

allocate (list(nD))
allocate (dslash(3,nD))
allocate (trueDistance(nD), tmp(nD), internalDistance(nD))

! reconstruction of planar basis:
! Take n to be one of the coordinate directions, and construct a and b accordingly to form a 
! full coordinate system. The n-direction is the first coordinate.
a = (/ -(n(3)**2+n(2)**2), -n(1)*n(2), n(1)*n(3) /)
a = a*1.0_dp / norm(a)

b = (/ 0._dp, n(3), -n(2) /)
b = b*1.0_dp / norm(b)

!inversion of transformation matrix
M(1,1:3) = n
M(2,1:3) = a
M(3,1:3) = b

!perform rotation 
forall(iD=1:nD)
  dslash(:,iD) = matmul(M, d(:,iD))
end forall
trueDistance(:) = dslash(1,:)   ! remember: the first coordinate in the new system is the n direction

if (nD==1) then ! only one d-vector: the distances of the surface(s) is quite clear. We have to issue
                ! return, since the following algorithm relies on at least 2 d-vectors
  SdistT = trueDistance(1)
  SdistB = trueDistance(1)
  return
endif

! sort the points according to their distance
list = (/(iD,iD=1,nD)/)
call heapsort(list, trueDistance)  ! trueDistance is now sorted...

!calculate intermediate distances (with respect to the n-direction)
forall(iD=1:nD)
  internalDistance(iD) = abs(trueDistance(iD) - trueDistance(mod(iD, nD) + 1))
end forall

!figure out between which atoms the vacuum is located
do iD=1, nD
   if(internalDistance(gapParameter) < internalDistance(iD)) gapParameter = iD
enddo !after this, we can be sure that the vacuum is between gapParameter and gapParameter+1

if(gapParameter == nD) then !this means: vacuum above and under the system
   sDistT = trueDistance(nD)
   sDistB = trueDistance(1)
else !this means: vacuum inside the system
   sDistT = trueDistance(gapParameter + 1)
   sDistB = trueDistance(gapParameter)
endif
end subroutine find_topbottom


subroutine read_adsorbate_sites(fh,nD,ads)
integer, intent(IN) :: fh      ! filehandle
integer, intent(IN) :: nD      ! number of dset points
logical, pointer    :: ads(:)  ! INTENT(OUT), { dset# }: which of the dvectors are adsorbate sites?

logical err
character(1000)    :: line
integer           :: nE
real(dp), pointer :: numbers(:)

if (.not. associated(ads)) allocate(ads(nD))
ads = .false.  ! a priori: no dvectors are adsorbate sites

call co_ca(fh,err)
read(fh,'(A1000)') line ! read line

call parse_line_for_numbers(line," ",maxEntries=nD,readEntries=nE,values=numbers) ! get the numbers of the adsorbate sites
if (nE==0) stop "ERROR: you're trying an adsorbate CE without any adsorbate"

ads(nint(numbers(:nE))) = .true. ! set "ads" of the adsorbate sites to .true.

end subroutine read_adsorbate_sites



logical function isThereAtom(d,x,n,eps)
! this function finds out if there is an atom in the plane (defined by the surface normal) 
! at point x(:)
real(dp), intent(in) :: d(:,:)  ! d-vectors
real(dp), intent(in) :: x(3)    ! actual position
real(dp), intent(in) :: n(3)    ! surface unit normal
real(dp), intent(in) :: eps     ! finite precision parameter

real(dp)  :: distV(3)  ! distance vector
integer   :: iD, nD

isThereAtom=.false.
nD = size(d,2)

do iD=1,nD
  DistV(:)=x(:)-d(:,iD)
  if (ProjectionLength(DistV,n)<eps) then
    isThereAtom=.true.
    return
  endif
enddo
end function isThereAtom


END SUBROUTINE read_lattdef


!***************************************************************************************************
!!GH These variables are no longer used anywhere in the code
!!SUBROUTINE read_bulkpar(CE)
!!type(CE_variables), intent(inout) :: CE
!!
!!logical :: err
!!
!!open(10,file='bulkpar.in',status='old')
!!call co_ca(10, err)
!!read(10,*) CE%leftoutStr
!!if (CE%leftoutStr/=0) stop "ERROR: please do no longer exclude any structures by the flag in bulkpar.in"
!!call co_ca(10, err)
!!read(10,*) CE%quasibin
!!call co_ca(10, err)
!!read(10,*) CE%concdep
!!call co_ca(10, err)
!!read(10,*) CE%polynome
!!call co_ca(10, err)
!!read(10,*) CE%cstrain
!!call co_ca(10, err)
!!read(10,*) CE%locstrain
!!close(10)
!!
!!ENDSUBROUTINE read_bulkpar




!***************************************************************************************************
SUBROUTINE read_CEfitting(CE, GA)
type(CE_variables), intent(inout) :: CE
type(GA_variables), intent(inout) :: GA

integer :: nDifferentDelta, i, iObs, nObs

integer :: iostat

logical :: err

character(80) line
integer           :: nE
real(dp), pointer :: numbers(:)

character(3) :: EOF


open(13,file='CEfitting.in',status='old',iostat=iostat)
if (iostat/=0) stop "ERROR: could not open file CEfitting.in"

!********************************************************************************
! GENERAL
!********************************************************************************

!--------------------------------------------------------------------------------
! General
call co_ca(13, err)
read(13,*) line
line=adjustl(trim(line))
call ucase(line)
if     (line(1:1)=='O') then; GA%simplefit = .false.
elseif (line(1:1)=='S') then; GA%simplefit = .true.
else;  stop "ERROR in CEfitting.in: Please specify 'O'ptimise or 'S'imple for the kind of fitting"
endif

GA%continue%best = .false.
GA%continue%dead = .false.
GA%continue%all  = .false.
GA%continue%indv = .false.
read(13,*) GA%continue%yes
call co_ca(13, err)
read(13,*) GA%continue%filename
call co_ca(13, err)
read(13,*) line
line=adjustl(trim(line))
call ucase(line)
if     (line(1:1)=='B') then; GA%continue%best = .true.
elseif (line(1:1)=='D') then; GA%continue%dead = .true.
elseif (line(1:1)=='A') then; GA%continue%all  = .true.
else;                         GA%continue%indv = .true.;  read(line,*) GA%continue%iIndv
endif

if (GA%continue%yes) then
  if (GA%simplefit .and. GA%continue%all) &   ! for a simple fit, we assume that the user really selects ONE of the individuals in the file, and not
       stop "ERROR in CEfitting.in: you cannot use 'all' individuals for a simple fit." ! uses ALL individuals.
  if (.not. GA%simplefit .and. (GA%continue%best .or. GA%continue%dead)) & ! for the GA, we assume that the user takes ALL individuals
       stop "ERROR in CEfitting.in: you must use 'all' individuals for continuing the genetic algorithm."
endif

!--------------------------------------------------------------------------------
! Clusters
call co_ca(13, err)
read(13,*) GA%nMaxClusters

GA%dampingC      = 0._dp
GA%dampingLambda = 0._dp
do i=2,6  ! i = number of vertices
  call co_ca(13, err)
  read(13,*) GA%dampingC(i), GA%dampingLambda(i)
enddo

!********************************************************************************
! GARBULSKY-CEDER CONSTRAINTS
!********************************************************************************

!--------------------------------------------------------------------------------
! General
call co_ca(13, err)
read(13,*) nDifferentDelta

nObs = CE%nObservables
allocate(CE%delta(3,nDifferentDelta,nObs))
allocate(CE%gen4DeltaSwitch(nDifferentDelta-1))
call co_ca(13,err)
do iObs=1,nObs; read(13,*) CE%delta(:,1,iObs); enddo
do i=2,nDifferentDelta
  call co_ca(13, err)
  read(13, *) CE%gen4DeltaSwitch(i-1), CE%delta(:,i,1)
  do iObs=2,nObs; read(13, *) CE%delta(:,i,iObs); enddo
enddo


!--------------------------------------------------------------------------------
! Special
call co_ca(13, err)
read(13, *) CE%pureDelta1
call co_ca(13, err)
read(13, *) CE%gstweight

!--------------------------------------------------------------------------------
! Smearing
call co_ca(13, err)
read(13, *) CE%ktfit      !! temperature factor


!********************************************************************************
! GENETIC ALGORITHM
!********************************************************************************

!--------------------------------------------------------------------------------
! Population
call co_ca(13, err)
read(13,*) GA%nPopulations
call co_ca(13, err)
read(13,*) GA%nIndividuals
call co_ca(13, err)
read(13,*) GA%maximumAgeForIndividual
call co_ca(13, err)
read(13,*) GA%diversityLimit
call co_ca(13, err)
read(13,*) GA%diversityEX, GA%diversityMB
!--------------------------------------------------------------------------------
! Children
call co_ca(13, err)
read(13,*) GA%nChildren
!if(mod(GA%nChildren,2)/=0) stop "Number of children should be even"
call co_ca(13, err)
read(13,*) GA%childWithOtherPopulation
call co_ca(13, err)
read(13,*) GA%nReplacements
if (GA%nReplacements > GA%nIndividuals)   stop "ERROR in CEfitting.in: you cannot have more individuals replaced than there are individuals"
if (GA%nReplacements > GA%nChildren) stop "ERROR in CEfitting.in: you cannot have more individuals replaced than there are children."
call co_ca(13, err)
!--------------------------------------------------------------------------------
! Simulation time
call co_ca(13, err)
read(13,*) GA%maxgenerations
!--------------------------------------------------------------------------------
! Cross validation
call co_ca(13, err)
read(13,'(A80)') line
call ucase(line)
if (line(1:1)=='L') then
  GA%CrossValidationLeaveOneOut = .true.
else
  GA%CrossValidationLeaveOneOut = .false.
  call parse_line_for_numbers(line," ",maxEntries=2,readEntries=nE,values=numbers)
  GA%CrossValidationK   = numbers(1)
  GA%CrossValidationRep = numbers(2)
endif

!--------------------------------------------------------------------------------
! Mutation
call co_ca(13, err)
read(13, *) GA%mutationProbability

!--------------------------------------------------------------------------------
! End of file :-)
call co_ca(13, err)
read(13, *) EOF
call ucase(EOF)
if (EOF /= "EOF") stop "ERROR in CEfitting.in: format changed. Last entry in CEfitting.in must be 'EOF'. Check input file."
close(13)
ENDSUBROUTINE read_CEfitting
!***************************************************************************************************


!***************************************************************************************************
subroutine co_ca(unit,error)
!           subroutine was taken from the code of Ralf Drautz
!
!           I: unit to read from
!    
!           co_ca cares about comments and blanks and avoids reading them,
!	        comment lines start with a #
implicit none
character(50) phrase !letter: contains the first letter of every line
logical   com !true if comment line is found
integer  unit, i, ios !unit specifies the unit to read from
logical error

com = .true.; error = .false.
do while ( com .eqv. .true.)
   read(unit,50,iostat=ios) phrase
   if (ios/=0) exit
   !              blank line ?
   if (phrase .ne. ' ') then
      i = index(phrase, '#')
      !                 # not found ?
      if (i .ne. 0) then
         !                    # first letter in line ?
         if (i .ne. 1) then
            phrase = phrase(1:i-1)
            if (phrase .ne. ' ') com= .false.
         endif
      else
         com = .false.
      endif
   endif
end do
50 format(50a)
backspace unit
if (ios/=0) error = .true.
end subroutine co_ca

!***************************************************************************************************
!! start the random number generator 
!! using the system time
SUBROUTINE initializeRND(GA,setseed,printseed)
  !         	intrinsic cpu_time
  intrinsic system_clock
  intrinsic date_and_time
  type(GA_variables), intent(inout) :: GA
  integer, optional :: setseed
  logical, optional :: printseed

  !                integer :: dummy2
  integer, dimension(8) :: dummy3 

  logical :: present_setseed, print_seed

  present_setseed  =.false.
  if(present(setseed))   present_setseed=.true.
  if(present(printseed)) then;  print_seed=printseed; else; print_seed=.true.; endif

  !! version for Absoft Fortran pro
  !!call cpu_time(dummy)
  !!rndseed = -1*int(dummy)
  
  !! version for Lahey Fortran 95
  !!call system_clock(dummy2)
  !!print *,dummy2
  !!rndseed = -1*dummy2

  !!version for pgf90

  ! is seed fixed (debugging)?
  if (GA%fixedseed /= 99) then
    GA%rndseed = -abs(GA%fixedseed)
  else
  ! seed is not fixed
    if (.not. present_setseed) then ! no seed was given => determine it from time
      call date_and_time(VALUES = dummy3)
      GA%rndseed = -1*(dummy3(7)*1000+dummy3(8))
    else ! seed was passed in to function
      GA%rndseed = setseed
    endif
  endif

  if (print_seed) write(*,*), "rnd(): seed = ", GA%rndseed
!  dummy = ran(GA%rndseed)
!  print *, dummy
END SUBROUTINE initializeRND

!***************************************************************************************************
!! wrapper function for the random number generator


REAL FUNCTION rnd_GA(GA)

type(GA_variables), intent(inout) :: GA
real(dp) :: a
        
a = ran(GA%rndseed)   !! call the actual random number generator
rnd_GA = a
end function rnd_GA


REAL FUNCTION rnd_seed(seed)

integer, intent(inout) :: seed
        
rnd_seed = ran(seed)   !! call the actual random number generator

END FUNCTION rnd_seed

!***************************************************************************************************
!! following routine was taken from Numerical recipes for Fortran Volume 2
!! generate a series of random numbers
function ran(idum)
  implicit none
  integer, INTENT(INOUT) :: idum
  real(dp) :: ran

  integer, parameter :: IA=16807, IM=2147483647, IQ=127773, IR=2836
  real(dp), SAVE :: am
  integer, SAVE :: ix=-1, iy=-1,k
  
  if (idum <= 0 .or. iy < 0) then     ! Initialize
     am=nearest(1.0,-1.0)/IM
     iy=ior(ieor(888889999,abs(idum)),1)
     ix=ieor(777755555,abs(idum))
     idum=abs(idum)+1
  end if
  ix=ieor(ix,ishft(ix,13))
  ix=ieor(ix,ishft(ix,-17))
  ix=ieor(ix,ishft(ix,5))
  k=iy/IQ
  iy=IA*(iy-k*IQ)-IR*k
  if (iy < 0) iy=iy+IM
  ran=am*ior(iand(IM,ieor(ix,iy)),1)
end function ran

!***************************************************************************************************
!! initialize the fortran intrinsic random number generator random_number
!! taken from:
!! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
subroutine initialize_random_number(seed,setseed,printseed)
implicit none

INTEGER :: i, n, clock
INTEGER, DIMENSION(:), pointer     :: seed
integer, dimension(:), optional    :: setseed
logical, optional                  :: printseed

if (.not. present(setseed)) then
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
else
  if (.not.associated(seed)) allocate(seed(size(setseed)))
  seed = setseed
endif

CALL RANDOM_SEED(PUT = seed)

if (.not.present(printseed)) printseed=.true.
if (printseed) write(*,*) "random_number(): seed = ", seed

!DEALLOCATE(seed)
end subroutine initialize_random_number

!***************************************************************************************************
SUBROUTINE writeCEresults(CE, GA, structures, CE_res, clusters)
use ce_types
implicit none

type(CE_variables), intent(in) :: CE
type(GA_variables), intent(in) :: GA
type(crystal), intent(in) :: structures(:)
type(t_fitResults), intent(in) :: CE_res(:)  ! { observable# }
!type(FigRep), intent(in) :: figures(:)
type(figure), intent(in) :: clusters(:)

integer :: nrdftgst, fitstrucnumb
integer :: i,j,l, Ncorr
real(dp) :: meanerr, maxerr, acterr, rms
real(dp) :: mixenergy, mixenergy2, mixenergy_dft
real(dp) :: DHf_DFT, DHf_CE
real(dp) :: energy, formen
real(dp) :: bulkenergya,bulkenergyb,bulkenergyc
real(dp) :: conca, concb
real(dp) :: cebulkenergya, cebulkenergyb, cebulkenergyc
real(dp) :: curreci

integer :: iC, ix, iObs, nObs
character :: filestring(40)

! adsorbate stuff
integer, allocatable :: adsList(:)
type(crystal), pointer :: tempstructures(:)
integer :: nStr, iRef, nRef

character(80) :: filename, cObs

nObs = CE%nObservables

do iObs=1,nObs
  !print *, "nOn", CE_res(iObs)%bestIndividual%nOnGenes
  !print *, "On ", CE_res(iObs)%bestIndividual%OnGenes(1:CE_res(1)%bestIndividual%nOnGenes)
  !print *, "cvs", CE_res(iObs)%cvs
  !print *, "Jav", CE_res(iObs)%Jav
  !print *, "Jff", CE_res(iObs)%Jff

  write(cObs,'(I2)') iObs
  filename="J."//trim(adjustl(cObs))//".out"
  call write_J_file(filename,clusters,CE_res(iObs))

  write(cObs,'(I2)') iObs
  filename="CVS."//trim(adjustl(cObs))//".out"
  call write_CVS_file(filename,structures,GA%subset(1:,CE_res(iObs)%bestPopulationIndex),CE_res(iObs)%predictionError(:,:),CE_res(iObs)%CVS)

  write(cObs,'(I2)') iObs
  filename="fittingErrors."//trim(adjustl(cObs))//".out"
  call write_fittingErrors_file(filename,structures,CE_res(iObs)%fittingError(:))

enddo


  
END SUBROUTINE writeCEresults
!***************************************************************************************************


!********************************************************************************
! Takes a character line and parses it for numbers which are then returned.
! (tk: I've taken the basic algorithm for that from Gus's labeling parser)
!
! --> IN:
!             line  :  the character string to be parsed
!        separator  :  the character by which the numbers are separated (e.g. ":", "/", ...)
!       maxEntries  :  the maximal number of numbers to be expected
!
! <-- OUT:
!      readEntries  :  the number of numbers that were read in (<= maxEntries)
!           values  :  the numbers that were read in
!
subroutine parse_line_for_numbers(line,separator,maxEntries,readEntries,values)
character(len=*), intent(inout):: line
character(1), intent(in)      :: separator
integer, intent(in)           :: maxEntries
integer, intent(out)          :: readEntries
real(dp), pointer             :: values(:) ! intent(OUT)                               

integer length
integer iE, ierr

nullify(values)
length = len(line)

!if (associated(values)) deallocate(values)
allocate(values(maxEntries))

readEntries=0
iE = 0

line = adjustl(line) ! Remove preceding blanks                                         
if (trim(line)=="") then ! make sure we didn't get a blank line
  readEntries = 0
  return
endif
             
line(length:length) = separator

do while (iE<maxEntries)
  iE = iE+1
  read(line,*,iostat=ierr) values(iE)
  if (ierr /=0) then
    readEntries=iE-1
    return 
    !stop "ERROR: line parsing for number failed. Maybe format mismatch?"
  endif
  line = adjustl(line(index(line,separator)+1:))  ! remove what we just read in                                            
  if (line=="") exit
enddo

readEntries=iE


end subroutine parse_line_for_numbers



!********************************************************************************
! Sets file pointer to the start of the enumeration section of fh
subroutine goto_enumStart(fh)
integer, intent(in) :: fh ! filehandle
character(5) startFlag

rewind(fh)
startFlag = ""
do while (.not. startFlag=="START")
   read(fh,*) startFlag
   call ucase(startFlag)
enddo
end subroutine goto_enumStart
!********************************************************************************



!********************************************************************************
! Skips nlines lines in fh
subroutine skip_lines(fh,nlines)
integer, intent(in) :: fh ! filehandle
integer, intent(in) :: nlines

character(1) :: dummy
integer iSkip

do iSkip=1, nlines
   read(17,*) dummy
enddo
end subroutine skip_lines
!********************************************************************************



!********************************************************************************
! Reads a character line and stops
subroutine io_debug(fh)
integer, intent(in) :: fh
character(80) :: line

read(fh,*) line
print *, "line read in:"
print *, line
stop "stop requested by io_debug"

end subroutine io_debug
!********************************************************************************



!********************************************************************************

subroutine concatenate_files(infiles,outfile,startComment,delete,verbose)
character(80), intent(in) :: infiles(:), outfile
character(80), intent(in) :: startComment
logical, intent(in)       :: delete
logical, intent(in)       :: verbose

integer :: nFiles, iFiles, fhout, fhin, iostat
character(1000) :: line

character(10000) :: systemcall

! Timer
real(dp) walltime1,cputime1           ! Time when starting task
real(dp) walltime2,cputime2           ! Time when task has ended

nFiles = size(infiles)

if (verbose) then
  print *
  write(*,'(A)') "Concatenating files:"
  do iFiles=1,nFiles
    write(*,'(I13,1x,A,1x,A)') iFiles, ":", infiles(iFiles)
  enddo
  write(*,'(A)') "to single output file:"
  write(*,'(12x,A)') outfile
  write(*,'(A$)') "The input files will "
  if (.not. delete) write(*,'(A$)') "NOT "
  write(*,'(A)') "be DELETED after the process."
  print *
endif

fhout = 100

call timing(walltime1,cputime1)
open(unit=fhout,file=outfile,status='NEW',iostat=iostat)
if (iostat /= 0) then
  write(*,'(A,A)') "WARNING: couldn't open NEW file ", adjustl(outfile)
  write(*,'(A,A)') "WARNING: I'll stop the concatenation process..."
  return
endif
write(fhout,'(A,A)') "# ", startComment

! The following commented lines contain the quick-and-dirty fortran concatenation.
! However, the system call "cat" is MUCH faster; that's why it's implemented below (after the comments)
!
!!if (verbose) write(*,'(A$)') "Processing: "
!!do iFiles=1,nFiles
!!  if (verbose) write(*,'(I3$)') iFiles
!!  fhin = fhout + iFiles
!!  open(unit=fhin,file=infiles(iFiles),status='OLD',iostat=iostat)
!!  if (iostat /= 0) then
!!    print *, "WARNING: couldn't open EXISTING file ", adjustl(infiles(iFiles))
!!    print *, "WARNING: I'll stop the concatenation process..."
!!    return
!!  endif
!!  do 
!!    read(fhin,'(A100)',iostat=iostat) line
!!    if (iostat /= 0) exit
!!    write(fhout,'(A)',iostat=iostat) trim(line)
!!    if (iostat /= 0) stop "ERROR: Write Error"
!!  enddo
!!  if (delete) then; close(fhin,status='DELETE'); else; close(fhin); endif
!!  if (verbose) write(*,'(A$)') ". "
!!enddo
!!if (verbose) print *

close(fhout)

! the system call is MUCH MUCH Faster!
systemcall="cat"
do iFiles=1,nFiles
  systemcall=trim(systemcall) // " " // adjustl(trim(infiles(iFiles)))
enddo
systemcall=trim(systemcall) // ">> gss.out"
call system(systemcall)

if (delete) then
  systemcall="rm -f"
  do iFiles=1,nFiles
    systemcall=trim(systemcall) // " " // adjustl(trim(infiles(iFiles)))
  enddo
  call system(systemcall)
endif


call timing(walltime2,cputime2)

if (verbose) then
  write(*,'(A,1x,F10.2,1x,A)') "Time for concatentation of files: ", walltime2-walltime1, "s"
endif

end subroutine concatenate_files
!********************************************************************************

  

!********************************************************************************
subroutine cmp_arrays1(array1,array2,cmp)
real(dp), intent(in) :: array1(:), array2(:)
integer, intent(out) :: cmp  ! -1: array1 < array2
                             !  0: array1 = array2
                             !  1: array1 > array2

integer i, n

n = size(array1)
if (size(array2) /= n) stop "ERROR in cmp_array: cannot compare 2 arrays of different sizes"

cmp = -99

do i=1, n

  if (array1(i)>array2(i)) then
    cmp = 1
    return
  endif
  
  if (array1(i)<array2(i)) then
    cmp = -1
    return
  endif

  ! loop only continues if array1(i)==array2(i)
enddo

! if we end up here, all array elements are equal
cmp = 0

end subroutine cmp_arrays1
!********************************************************************************



!********************************************************************************
subroutine cmp_arrays2(array1,array2,cmp)
real(dp), intent(in) :: array1(:,:), array2(:,:)
integer, intent(out) :: cmp  ! -1: array1 < array2
                             !  0: array1 = array2
                             !  1: array1 > array2

integer n

n = size(array1)

call cmp_arrays1( reshape(array1,(/n/)), reshape(array2,(/n/)), cmp)

end subroutine cmp_arrays2
!********************************************************************************



!********************************************************************************
subroutine write_multilattice_dvector_occupation(CE,MPI_rank)
type(CE_variables), intent(in) :: CE
integer, intent(in) :: MPI_rank

integer :: HNF(3,3), a, b, c, d, e, f, nUC, strN, hnfN, sizeN, pgOps, diag(3), Hdegen, labdegen,Tdegen
integer :: L(3,3) ! Left transform matrix for SNF
character(maxLabLength) labeling
integer, pointer ::  ilabeling(:)

integer iostat
integer :: iLab, iD, nD
integer, allocatable :: Docc(:), labelpos2dsetlabel(:)

nD = CE%nBas
allocate(Docc(nD))

if (nD>100) stop "ERROR: bad programming in write_multilattice_dvector_occupation. Check write format strings"

open(100,file="struct_enum.out",status='OLD',iostat=iostat)
if (iostat/=0) stop "ERROR: coudln't open file struct_enum.out"

open(200,file="dvector_occ.out",status='NEW',iostat=iostat)
if (iostat/=0) stop "ERROR: couldn't create new file dvector_occ.out"

write(200,'(A)') "# struc num | label | for each dvector: number of occupancies with label"
write(200,'(A)') "#--------------------------------------------------------------------------------"

call goto_enumStart(100)

do
  read(100,*,iostat=iostat) strN, hnfN, Hdegen,labdegen,Tdegen,sizeN, nUC, pgOps, diag, a,b,c,d,e,f, L, labeling

  if (iostat/=0) exit

  call get_ilabeling(labeling,ilabeling)

  allocate(labelpos2dsetlabel(nD*nUC))
  labelpos2dsetlabel =  map_labeling_position_2_dvector_index(nUC=nUC,nD=nD)

  do iLab=0,CE%CErank-1
    do iD=1,nD
      
      Docc(iD) = count(  pack(ilabeling,labelpos2dsetlabel==iD)==iLab ) ! how often does iLab occur in the 
                                            ! labeling, provided that we only look at those labels that
                                            ! correspond to the dvector iD
    enddo

    write(200,'(I15,5x,I4,5x,100(I4,1x))') strN, iLab, Docc(:)

  enddo
  

  deallocate(labelpos2dsetlabel)
  deallocate(ilabeling)
end do

close(200)
close(100)
end subroutine write_multilattice_dvector_occupation
!********************************************************************************




!********************************************************************************
! Purpose: read in options from structpred.in
subroutine read_structpred(CE,predictions)
type(CE_variables), intent(in)            :: CE
type(structure_prediction_t), intent(out) :: predictions

integer :: n, i, ioerr
logical :: err

write(*,'(A)') "Reading data for struture search..."
open(20, file='structpred.in', status='OLD', iostat=ioerr)
if (ioerr/=0) stop "file: structpred.in was not found!"

call co_ca(20,err)
read(20,*) predictions%coordsys
call ucase(predictions%coordsys)
if (predictions%coordsys/='D' .and. predictions%coordsys/='C') stop "ERROR: in structpred.in. First entry must be coordinate system used for output"


call co_ca(20,err)
read(20,*) predictions%compareStrucs


call co_ca(20,err)
read(20,*) n

if (n<0) then
   ! want to do a manual search:
   predictions%automatic = .false.
   predictions%nS        = -n       ! number of structures
else 
   predictions%automatic = .true.
   predictions%nSPX      = n        ! number of structures per concentration (SPX)
endif

allocate(predictions%xi(CE%CErank,CE%ncouplings),predictions%xf(CE%CErank,CE%ncouplings))
do i=1,CE%ncouplings
  call co_ca(20,err)
  read(20,*) predictions%xi(:,i), predictions%xf(:,i)
enddo

if (.not. predictions%automatic) then ! manual input selected
  call co_ca(20,err)
  allocate(predictions%struclist(predictions%nS))
  do i=1, predictions%nS
    read (20,*) predictions%struclist(i) ! read in list
  end do
endif
close(20)
! finished reading from structpred.in

end subroutine read_structpred
!****************************************************************************************************






!****************************************************************************************************
subroutine write_InpStrucList(file,strucs,header)
character(80), intent(in)  :: file      ! filename
type(crystal), intent(in), target  :: strucs(:) ! the structures to write out
character(80), intent(in)  :: header(:) ! header information of the file
type(crystal), pointer     :: currStr
integer i, n, iStr, nStr
integer iRef, nRef
real(dp) :: invLV(3,3)
logical err

nStr=size(strucs)
if (nStr>=1) then
  nRef=size(strucs(1)%reference)
else
  return
endif

open(100,file=file,status='REPLACE')
do i=1,size(header)
   write(100,'(A)') header(i)
end do

write(100,'(A)') "peratom"
write(100,'(A)') "weights"
write(100,'(A)') "POSCAR"
do iStr = 1, nStr
  currStr => strucs(iStr)
  call matrix_inverse(currStr%LV,invLV,err)
  write(100,'(A)') "#================================================================================"
  write(100,'(A)') currStr%title
  write(100,'(A)') "1    # lattice constant for CE"
  do i=1,3; write(100,*) currStr%LV(:,i); enddo
  write(100,*) currStr%nTyp
  write(100,'(A)') "Direct"
  do i=1,currStr%nAtoms
    write(100,*) matmul(invLV,currStr%pos(:,i))
  enddo
  write(100,'(A)')"# Energy:"
  write(100,*) currStr%energy
  write(100,'(A)')"# Weight:"
  write(100,*) currStr%weight
  do iRef=1,nRef
      write(100,'(A)')"# References:"
      if (currStr%reference(iRef)%automaticReference ) then
        write(100,'(A,1x,I3)')"Automatic Reference ", iRef
      else
        write(100,'(A)')"Sorry. Manual references not supported for structures.out."
      endif
  enddo
enddo
close(100)

end subroutine write_InpStrucList
!****************************************************************************************************



!****************************************************************************************************
subroutine write_fittingErrors_file(filename,structures,fittingError)
character(80), intent(in)      :: filename
type(crystal), intent(in)      :: structures(:)
!sascha>real(dp), intent(in)           :: allowedDelta(:)
real(dp), intent(in)           :: fittingError(:)

real(dp), parameter :: xeps = 0.00001

character(1000) formatstring, formatstring_attention
integer iStr, nStr, iX, nX, nC, iC

real(dp), pointer :: currX(:,:), prevX(:,:)
real(dp)          :: currE,      prevE

nStr = size(structures)
nX   = size(structures(1)%x,1)   ! number of concentrations per coupling
nC   = size(structures(1)%x,2)   ! number of CE couplings

allocate(currX(nX,nC))
allocate(prevX(nX,nC))

formatstring="(I4,1x,"
formatstring_attention="(I4,1x,"
do iC=1,nC
do iX=1,nX
  formatstring=trim(adjustl(formatstring))//"f8.5,2x,"
  formatstring_attention=trim(adjustl(formatstring_attention))//"f8.5,2x,"
enddo
enddo !I messed around here 1/Oct/2012 -sascha
formatstring          =trim(adjustl(formatstring))          //"3x,6(f15.10,2x),  1x,3x,1x,  A)"
formatstring_attention=trim(adjustl(formatstring_attention))//"3x,6(f15.10,2x),  1x,A3,1x,  A)"

open(100,file=filename,status='replace')
write(100,'(A)') "# no. | x(1:k,1:nC) | E_DFT | D(E_DFT,E_CE) | E_CE | Delta1 (GC) | DHf_DFT | DHf_CE | Dist to GSL DFT | title"
write(100,'(A)') "#------------------------------------------------------------------------------------------------------------"
do iStr=1,nStr


!print*, formatstring
!debug>print*, iStr
!debug>print*, structures(istr)%x
!debug>print*, structures(iStr)%energy
!debug>!print*, fittingError(iStr)
!debug>print*, structures(iStr)%ce_energy !sascha> allowedDelta(iStr), &                                        
!debug>print*, structures(iStr)%formationenergy 
!debug>print*, "and", structures(iStr)%ce_formationenergy
!debug>print*, "problemski: ",structures(iStr)%distanceToInputGSL
!debug>print*, trim(adjustl(structures(iStr)%title))


   write(100,formatstring) &
        iStr, &
        structures(iStr)%x, &
        structures(iStr)%energy, fittingError(iStr), structures(iStr)%ce_energy, & !sascha> allowedDelta(iStr), &
        structures(iStr)%formationenergy, structures(iStr)%ce_formationenergy, &
        structures(iStr)%distanceToInputGSL, &
        trim(adjustl(structures(iStr)%title))


!sascha  if (iStr==1) then
!sascha    prevX = structures(iStr)%x
!sascha    prevE = structures(iStr)%ce_energy
!sascha  endif
!sascha  currX = structures(iStr)%x
!sascha  currE = structures(iStr)%ce_energy

  ! Note: 
  ! The flags -X- and H-- and HX- that are set by the following loop should be handled much more 
  ! elegantly than what is done here.
  !
  ! "X"    this structure has a larger deviation from the input than allowed by the 1st GC constraint
  ! "H"    this structure has a problem with the energetical hierarchy


!  if (.not. all(currX-prevX<xeps)) then
!    write(100,'(A)') "#--------------------------------------------------------------------------------"
!  endif
!  if ( abs(fittingError(iStr)) <= allowedDelta(iStr) ) then
!    if ( (all(currX-prevX<xeps) .and. (currE <= prevE)) .or. (.not. all(currX-prevX<xeps)) ) then
    ! check if the CE hierarchy represents the DFT hierarchy: if the concentrations match, see if the
    ! energy decreases. This should be the case, since the structures are ordered by their DFT energy.
    ! This is a very rudimentary check, should be made better...
!sascha>      prevE = currE
!sascha>      prevX = currX
!sascha>    else
!sascha>      write(100,formatstring_attention) &
!sascha>         iStr, &
!sascha>         structures(iStr)%x, &
!sascha>         structures(iStr)%energy, fittingError(iStr), structures(iStr)%ce_energy, allowedDelta(iStr), &
!sascha>         structures(iStr)%formationenergy, structures(iStr)%ce_formationenergy, &
!sascha>         structures(iStr)%distanceToInputGSL, &
!sascha>         "H--", trim(adjustl(structures(iStr)%title))
!sascha>      !prevE = currE
!sascha>      !prevX = currX
!sascha>    endif
!sascha>
!sascha>  else
!sascha>
!sascha>    if ( (all(currX-prevX<xeps) .and. (currE <= prevE)) .or. (.not. all(currX-prevX<xeps)) ) then
!sascha>    ! see comments above about the hierarchy
!sascha>      write(100,formatstring_attention) &
!sascha>         iStr, &
!sascha>         structures(iStr)%x, &
!sascha>         structures(iStr)%energy, fittingError(iStr), structures(iStr)%ce_energy, allowedDelta(iStr), &
!sascha>         structures(iStr)%formationenergy, structures(iStr)%ce_formationenergy, &
!sascha>         structures(iStr)%distanceToInputGSL, &
!sascha>         "-X-", trim(adjustl(structures(iStr)%title))
!sascha>      prevE = currE
!sascha>      prevX = currX
!sascha>    else
!sascha>      write(100,formatstring_attention) &
!sascha>         iStr, &
!sascha>         structures(iStr)%x, &
!sascha>         structures(iStr)%energy, fittingError(iStr), structures(iStr)%ce_energy, allowedDelta(iStr), &
!sascha>         structures(iStr)%formationenergy, structures(iStr)%ce_formationenergy, &
!sascha>         structures(iStr)%distanceToInputGSL, &
!sascha>         "HX-", trim(adjustl(structures(iStr)%title))
!sascha>      !prevE = currE
!sascha>      !prevX = currX
!sascha>    endif
!sascha>  endif
enddo

write(100,'(A, F15.5)') "RMS Error: ", sqrt(sum(fittingError(:)**2)/nStr)
close(100)

end subroutine write_fittingErrors_file
!****************************************************************************************************




!****************************************************************************************************
subroutine   write_CVS_file(filename,structures,subset,predictionErrors,cvs)
character(80), intent(in)      :: filename
type(crystal), intent(in)      :: structures(:)            ! { structures# }
type(t_subset), intent(in)     :: subset(:)                ! { subset# }
real(dp), intent(in)           :: predictionErrors(:,:)    ! { structures#, subset# }
real(dp), intent(in)           :: cvs

integer iSS, nSS, iStr, nStr

nSS = size(subset)

open(100,file=filename,status='replace')
do iSS=1,nSS
  nStr=subset(iSS)%nPred
  write(100,'(A)') "--------------------------------------------------------------------------------"
  write(100,'(A,I3)') "Subset #", iSS
  write(100,'(A)') "--------------------------------------------------------------------------------"
  do iStr=1,nStr
    write(100,'(I5,2x,I5,2x,F15.10,10x,A100)') iStr, subset(iSS)%predList(iStr), predictionErrors(iStr,iSS), structures( subset(iSS)%predList(iStr) )%title
  enddo
enddo

write(100,'(A)') ""
write(100,'(A)') "=============================="
write(100,'(A)') "Final CVS:"
write(100,'(F15.10)') cvs
write(100,'(A)') "=============================="
write(100,'(A)') ""

close(100)

end subroutine write_CVS_file
!****************************************************************************************************






!****************************************************************************************************
subroutine write_J_file(filename,clusters,results)
character(80), intent(in)      :: filename
type(figure), intent(in)       :: clusters(0:)     ! { cluster# }
type(t_fitResults), intent(in) :: results

character(80) :: summaryFilename
integer iCl, nCl, iClOrig, iV, dotoutIdx

dotoutIdx       = index(filename,".out")
summaryFilename = trim(adjustl(filename(1:dotoutIdx-1)))//".summary.out"

nCl    = results%bestIndividual%nOnGenes


open(100,file=filename,        status='replace')
open(101,file=summaryFilename, status='replace')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Please note that the formatting for the J-file below is repeating in the routine write_final_Js
!! in the cs_utility.f90 file. If you change the format of the J-file here, make sure that you
!! change it there too.  --GLWH July 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(100,'(A)') "# number of clusters (including constant term)"
write(100,*) nCl

write(100,'(A)') "#"
write(100,'(A)') "#--------------------------------------------------------------------------------"
iCl=0
write(100,*) "constant", results%Jav(iCl), results%Jff(iCl)

do iCl=1,nCl-1   ! go through all clusters that are activated in the genome (constant term separately, see above)
  write(100,'(A)') "#--------------------------------------------------------------------------------"
  iClOrig= results%bestIndividual%onGenes(iCl+1)   ! +1 because the first entry is the constant term (cluster #0)
  write(100,'(A)') "# Cluster number: in this list | original in clusters.out"
  write(100,'(I6,5x,I6)') iCl, iClOrig
  write(100,'(A)') "# J_average | J_fullfit"
  write(100,*) results%Jav(iCl), results%Jff(iCl)
  write(100,'(A)') "# Number of vertices"
  write(100,'(I3)') clusters( iClOrig )%nV
  write(100,'(A)') "# Average Distance"
  write(100,'(f12.8)') clusters( iClOrig )%avgR
  write(100,'(A)') "# Damping"
  write(100,'(f12.8)') clusters( iClOrig )%damping
  write(100,'(A)') "# Vertices: (x,y,z) | d-vector label | s-index"
  do iV=1, clusters( iClOrig )%nV
    write(100,'(3(f12.8,2x),5x,I3,3x,I3)') clusters( iClOrig )%vertex(:,iV), clusters( iClOrig )%label(iV), clusters( iClOrig )%s(iV)
  enddo
enddo

close(100)

write(101,'(A1,1x,5A18)') "#","cluster# |", "#vertices |", "Jav |", "avg dist |", "damping |"
write(101,'(A)')          "#------------------------------------------------------------------------------------------"
do iCl=0,nCl-1
  iClOrig= results%bestIndividual%onGenes(iCl+1) 
  write(101,'( 2(2x,i16), 3(2x,f16.8))') iCl, clusters(iClOrig)%nV, results%Jav(iCl), clusters(iClOrig)%avgR, clusters(iClOrig)%damping
enddo
close(101)

end subroutine write_J_file
!****************************************************************************************************




!***************************************************************************************************
subroutine read_J_file(filename, cluster, Jav, Jff)
character(80), intent(in) :: filename
type(figure), pointer     :: cluster(:)
real(dp), pointer         :: Jav(:), Jff(:)   

real(dp), pointer :: numbers(:)
character(1000) :: line
character(8)    :: dummy
integer iCl, nCl, iV
integer iostatus, nE
logical err

!print *, "read_J_file"

open(100,file=filename, status='old', iostat=iostatus)
if (iostatus/=0) then
  write(*,'(A)') "ERROR: could not open file ", trim(adjustl(filename))
  stop
end if

call co_ca(100,err)
read(100,*) nCl   ! number of clusters

allocate(Jav(0:nCl-1), Jff(0:nCl-1))
allocate(cluster(0:nCl-1))

! setup constant term as cluster #0
allocate( cluster(0)%vertex(3,0) )
allocate( cluster(0)%label(0)    )
allocate( cluster(0)%s(0)        )
cluster(0)%nV = 0
call co_ca(100,err)
read(100,*) dummy, Jav(0), Jff(0)

! setup all other clusters
do iCl=1,nCl-1
  call co_ca(100,err)
  read(100,*) dummy    ! cluster numbers, don't care about that
  call co_ca(100,err)
  read(100,*) Jav(iCl), Jff(iCl)
  call co_ca(100,err)
  read(100,*) cluster(iCl)%nV
  allocate( cluster(iCl)%vertex(3,cluster(iCl)%nV) )
  allocate( cluster(iCl)%label(   cluster(iCl)%nV) )
  allocate( cluster(iCl)%s(       cluster(iCl)%nV) )
  call co_ca(100,err)
  read(100,*) cluster(iCl)%avgR
  call co_ca(100,err)
  read(100,*) cluster(iCl)%damping
  do iV=1,cluster(iCl)%nV
    call co_ca(100,err)
    read(100,*) cluster(iCl)%vertex(:,iV), cluster(iCl)%label(iV), cluster(iCl)%s(iV)
  enddo
enddo

close(100)

print*, "read in clusters!"

end subroutine read_J_file
!***************************************************************************************************



!***************************************************************************************************
subroutine save_population( individuals, dead, iPopulation, iGeneration, iObs )
type(t_individual), intent(in)     :: individuals(:)
type(t_individual), intent(in)     :: dead
integer, intent(in)                :: iPopulation
integer, intent(in)                :: iGeneration
integer, intent(in)                :: iObs

character(80) filename
character(6)  cObs, cGeneration, cPopulation
integer iIndv, nIndv
logical deadExists

character(80) :: format_onGenes, format_genome, c_integer


nIndv = size(individuals)
if (dead%cvs /= -1 .and. .not. isnan(dead%cvs)) then
!Stop "Can't get gfortran 4.2 to compile 'isnan'"
  deadExists=.true.
else
  deadExists=.false.
endif

write(cPopulation,'(I3)') iPopulation
write(cGeneration,'(I6)') iGeneration
write(cObs,'(I2)') iObs

if (iGeneration==0) cGeneration="last"

filename="population.g"//trim(adjustl(cGeneration))//".o"//trim(adjustl(cObs))//".p"//trim(adjustl(cPopulation))//".out"

open(100,file=filename,status='replace')
write(100,'(A)') "# Total number of individuals: alive | dead"
if (deadExists) then; write(100,*) nIndv, 1;
                else; write(100,*) nIndv, 0
endif
write(100,'(A,I5)') "# Best individual is: "
if (deadExists) then
  if (dead%CVS<minval(individuals(:)%CVS)) then
    write(100,*) nIndv+1
  else
    write(100,*) minloc(individuals(:)%CVS)
  endif
else 
  write(100,*) minloc(individuals(:)%CVS)
endif
write(100,'(A)') "#"
write(100,'(A)') "# For each individual:"
write(100,'(A)') "#   - number of individual"
write(100,'(A)') "#   - CVS"
write(100,'(A)') "#   - number of genes that are switched on"
write(100,'(A)') "#   - the indices of the genes that are switched on"
write(100,'(A)') "#   - the full genome"
write(100,'(A)') "#"




do iIndv=1,nIndv
  
  ! setting up format (don't need for gfortran, but ifort)
  write(c_integer,*) individuals(iIndv)%nOnGenes
  format_onGenes = "("//trim(adjustl(c_integer))//"(I5,1x))"
  write(c_integer,*) size(individuals(iIndv)%genome)
  format_genome  = "("//trim(adjustl(c_integer))//"(L1,1x))"

  ! write
  write(100,'(A)') "#--------------------------------------------------------------------------------"
  write(100,             *) iIndv
  write(100,             *) individuals(iIndv)%CVS
  write(100,             *) individuals(iIndv)%nOnGenes
  write(100,format_onGenes) individuals(iIndv)%onGenes( 1:individuals(iIndv)%nOnGenes )
  write(100,format_genome ) individuals(iIndv)%genome
enddo

if (deadExists) then

  ! setting up format (don't need for gfortran, but ifort)
  write(c_integer,*) dead%nOnGenes
  format_onGenes = "("//trim(adjustl(c_integer))//"(I5,1x))"
  write(c_integer,*) size(dead%genome)
  format_genome  = "("//trim(adjustl(c_integer))//"(L1,1x))"

  ! write
  write(100,'(A)') "#--------------------------------------------------------------------------------"
  write(100,'(A)') "# Best dead individual:"
  write(100,             *) nIndv + 1
  write(100,             *) dead%CVS
  write(100,             *) dead%nOnGenes
  write(100,format_onGenes) dead%onGenes( 1:dead%nOnGenes )
  write(100,format_genome ) dead%genome
endif

close(100)

end subroutine save_population
!***************************************************************************************************




!***************************************************************************************************
subroutine read_population( individuals, dead, iPop, iObs, continueParameters )
type(t_individual), intent(out)     :: individuals(:)
type(t_individual), intent(out)     :: dead
integer, intent(in)                 :: iPop, iObs   ! the index of the population and of the observable
type(t_GA_continue), intent(in)     :: continueParameters

type(t_individual), pointer :: tmpIndividual

character(80) filename, line
character(6)  cObs, cPop, cGeneration
integer iIndv, iStoreIndv, nIndv, nIndvAlive, nIndvDead, iIndvBest, dummy, iostat, OBSindex, POPindex
logical err


nIndv = size(individuals)
dead%cvs = -1   ! failsafe if dead invidiual is not specified in file

write(cObs,'(I2)') iObs
write(cPop,'(I2)') iPop

! Determine the correct filename:
filename=continueParameters%filename
filename=adjustl(trim(filename))
! replace $OBS 
OBSindex=index(filename,"$OBS")
filename=filename(1:OBSindex-1)//trim(adjustl(cObs))//filename(OBSindex+4:)
! replace $POP
POPindex=index(filename,"$POP")
if (POPindex /= 0) &   ! there is $POP in the string
     filename=filename(1:POPindex-1)//trim(adjustl(cPop))//filename(POPindex+4:)
! done


open(100,file=filename,status='old',iostat=iostat)
if (iostat/=0) then
  write(*,'(A)') "ERROR: could not open file ", filename
  stop "ERROR: could not open file."
endif

write(*,'(/,A,I3,A,A,A/)') "Reading population #", iPop, " from file """, trim(adjustl(filename)), """"

! Important:
! If a simple fit is selected, then we assume that we only want to read ONE individual, and nIndv==1.
! For the continuation of the GA, we assume that we read ALL individuals.

call co_ca(100, err)
read(100,*) nIndvAlive, nIndvDead
call co_ca(100, err)
read(100,*) iIndvBest

iStoreIndv = 1
do iIndv=1,nIndvAlive + nIndvDead
  if ( (continueParameters%indv .and. iIndv/=continueParameters%iIndv) &
  .or. (continueParameters%best .and. iIndv/=iIndvBest) &
  .or. (continueParameters%dead .and. iIndv<=nIndvAlive) ) then ! didn't find the right individuum
    call co_ca(100, err); read(100,*) line
    call co_ca(100, err); read(100,*) line
    call co_ca(100, err); read(100,*) line
    call co_ca(100, err); read(100,*) line
    call co_ca(100, err); read(100,*) line
    cycle
  endif
 
  if (iIndv<=nIndvAlive .or. continueParameters%best .or. continueParameters%dead) then
    call co_ca(100, err); read(100,*) dummy
    call co_ca(100, err); read(100,*) individuals(iStoreIndv)%CVS
    call co_ca(100, err); read(100,*) individuals(iStoreIndv)%nOnGenes
    call co_ca(100, err); read(100,*) individuals(iStoreIndv)%onGenes( 1:individuals(iStoreIndv)%nOnGenes )
    call co_ca(100, err); read(100,*) individuals(iStoreIndv)%genome
    iStoreIndv = iStoreIndv + 1
  else
    call co_ca(100, err); read(100,*) dummy
    call co_ca(100, err); read(100,*) dead%CVS
    call co_ca(100, err); read(100,*) dead%nOnGenes
    call co_ca(100, err); read(100,*) dead%onGenes( 1:dead%nOnGenes )
    call co_ca(100, err); read(100,*) dead%genome
  endif
enddo



close(100)
end subroutine read_population
!***************************************************************************************************


END MODULE io_utilities

