!
!
! WARNING:
!    Be VERY CAREFUL with the J and PI:
!    The EMPTY CLUSTER has index 0 in J and PI!
!    If you pass J or PI to a routine, either use the POINTER attribute (and hence deferred shape) or
!    specify the range by, say, J(0:) if you want to use INTENT(IN).
!    INTENT(IN) does NOT defer the ubound and lbound, it always assumes 1: shape!!
!
!

#if 0
!================================================================================
!
! Important:
!
! This fortran file contains preprocessor directives.
!    #ifdef __MKL__
!    <MKL code>
!    #endif
! Thus, you have to preprocess this file (which is usually done by the 
! Makefile, so you don't have to worry about that) before compiling it.
!
!================================================================================
#endif


MODULE ce_fit
use ce_types
use io_utilities
use structure_utilities
use numerical_utilities
use utilities_module, only: ralloc
use combinatorics, only: permutation_parity
use enumeration_utilities !added
use setup_cluster_expansion
use sort

implicit none


logical :: debug = .false.

private generate_subsets, make_fitting_sets
public do_clusterexpansion, simple_prediction_of_inputs, get_enumStrcorrs


CONTAINS




!***************************************************************************************************
! This routine takes a list of structures, and their correlations, and a set of ECIs and calculates
! the CE energy for each one.
SUBROUTINE simple_prediction_of_inputs(CE, structures, J)
type(CE_variables), intent(in) :: CE
type(crystal), intent(inout)   :: structures(:)
real(dp), intent(in)           :: J(0:)

integer iStr, nStr
integer iObs, nObs

nStr = size(structures)
nObs = CE%nObservables

do iObs=1,nObs
print *, "J",J
  do iStr = 1, nStr
     structures(iStr)%ce_energy = structures(iStr)%ce_refEnergy + calculate_energy(J,structures(iStr)%PI)
     !write(*, '("PI: ", 3(f9.5,1x))'),structures(iStr)%PI
     !write(*, '("PI"')),structures(iStr)%PI
     !write(*, '("En: ", 3(f9.5,1x))'),structures(iStr)%ce_energy
  enddo

!  if(debug) then
     open(37,file="fitted_energies.out")
     write(37, '("# Conc | DFT-DHf     | CE-DHf      | DeltaE     |Structurename         ")')
     do iStr = 1, nStr
        write(37, '(F7.4, 1x, 3(F12.5,1x),A)') structures(iStr)%conc, &
             structures(iStr)%energy,structures(iStr)%ce_energy, &
             structures(iStr)%energy-structures(iStr)%ce_energy, adjustl(trim(structures(iStr)%title))
     enddo
     close(37)
!  endif

enddo

ENDSUBROUTINE simple_prediction_of_inputs
!***************************************************************************************************





!***************************************************************************************************
! Driver routine for the fitting
SUBROUTINE do_clusterexpansion(CE, GA, structures, eqvClusters, results)

type(CE_variables), intent(inout) :: CE
type(GA_variables), intent(inout) :: GA
type(crystal), intent (inout)     :: structures(:)    ! { structure# }
type(FigRep), intent(in)          :: eqvClusters(0:)  ! { cluster# }
type(t_fitResults), pointer       :: results(:)       ! { observable# }
integer :: iObs

call read_CEfitting(CE, GA)

if (GA%nMaxClusters>size(eqvClusters)) then
  GA%nMaxClusters = size(eqvClusters)
  write(*,'(A)') "WARNING: The cluster pool is smaller than the maximum number of clusters allowed in GApar.in."
end if

allocate(results(CE%nObservables))

if (GA%simplefit) then ! run a simple cluster expansion?
  !GA%nSubsets       = 0
  GA%CrossValidationK   = 0
  GA%CrossValidationRep = 0
  GA%maxGenerations     = 0
  GA%nIndividuals       = 1
  GA%nPopulations       = 1
#ifdef __MKL__
  print *
  write(*,'(A)') "WARNING: You are running a simple fit and have the linear independence check"
  write(*,'(A)') "WARNING:    enabled. You have to make sure that the clusters that you provide"
  write(*,'(A)') "WARNING:    for the simple fit are linear independent. Otherwise, the fit will fail!"
  print *
#endif
endif

call Generate_subsets(GA, CE, structures)  ! generate fitting and prediction subsets

! Loop over all observables and do a separate genetic algorithm.
! (The initial setup was to do this loop inside the GA routine; however, I (tk) think it might
! be better to do a big loop here (Sascha?); then, OMP or even MPI parallelisation should be rather easy...)
! Notes:
! - think about separate input parameters for GA for the different observables
!   Different deltas can already be read in (see io_utilities.f90
! - what about totally different CEfitting.in files for the different observables?
!   That would mean: reading GA(iObs)?
do iObs=1, CE%nObservables

   ! run the genetic algorithm to find the best genes
   ! Note: if simpleFit is activated, we only have 1 subset and
   ! do not enter the generations loop in the GA.
   print *, 'entering genetic algorithm'
   call geneticAlgorithm_for_clusterExpansion(GA, CE, eqvClusters, structures, iObs, results(iObs))  
   
   ! set the ce values for the input structures
   call CE_for_input(structures, results(iObs))  

   ! calculate the formation enthalpy for the input structures
   if (iObs==1) call calculate_DHf_for_input(CE, structures)  

end do ! observables

END SUBROUTINE do_clusterexpansion
!***************************************************************************************************





!***************************************************************************************************
subroutine calculate_DHf_for_input(CE, structures)
type(CE_variables), intent(in)     :: CE
type(crystal), intent(inout)       :: structures(:)
type(crystal), pointer   :: adsStructures(:)
integer iStr, jStr, nStr, ik

nStr = size(structures)


if (CE%ncouplings==1) then
  do iStr=1,nStr
    structures(iStr)%ce_formationenergy = calc_structure_DHf(structures, &
                                                             iStr, &
                                                             CE%pureList, &
                                                             isCE = .true.)
  enddo
elseif (CE%ncouplings == 2 .and. CE%adsorbate%isAds) then
  allocate(adsStructures(CE%CErank))   ! there are CErank additional adsorbate structures
  do ik=1,CE%CErank
    adsStructures(ik)%energy = CE%adsorbate%singleEnergy(ik)
  enddo

  do iStr=1,nStr
    structures(iStr)%ce_formationenergy = calc_structure_DHf( (/ structures(:), adsStructures(:) /), &
                                                             iStr, &
                                                             pureEl= CE%pureList, &
                                                             isCE  = .true.,      &
                                                             adsEl = (/ (jStr, jStr=nStr+1,nStr+CE%CErank) /))

  enddo
else
  stop "ERROR in calculate_DHf_for_input: Formation enthalpies only for non-coupled CEs or coupled adsorbate system CE"
endif

end subroutine calculate_DHf_for_input
!***************************************************************************************************




!***************************************************************************************************
subroutine CE_for_input(structures, results)
type(crystal), intent(inout)  :: structures(:)! { structure# }
type(t_fitResults), intent(in):: results

integer iStr, nStr

nStr = size(structures)

do iStr=1,nStr
  structures(iStr)%ce_energy = dot_product( &
                                 results%Jav(0:) , &
                                 structures(iStr)%PI( results%bestIndividual%onGenes( 1:results%bestIndividual%nOnGenes ) ) &
                               ) + structures(iStr)%CE_refEnergy
enddo


end subroutine CE_for_input
!***************************************************************************************************











!***************************************************************************************************
!! Generate the subsets of structures to be included
SUBROUTINE generate_subsets(GA, CE, str)
use ce_types                      !! general variables
implicit none

type(GA_variables), intent(inout) :: GA 
type(CE_variables), intent(in) :: CE 
type(crystal), intent (inout) :: str(:)

integer :: startp, endp, fitstrucnumb, nrdftgst
integer :: ExDim    ! the dimension of the E-x-diagram
integer :: i, j,l,m
integer :: temp
real(dp):: rnd1
integer :: pos

integer, pointer :: indx(:)
integer, pointer :: GSlist(:), AdjStr(:,:), StrLowE(:), PureList(:)
real(dp), pointer :: DistList(:)
real(dp), pointer :: mixture(:,:)
real(dp), pointer :: deltaList(:,:)
real(dp), allocatable :: pureConditions(:,:)
integer, pointer :: excludeList(:), predList(:), fitList(:)
integer iSS, nSS, iPop, nPop
!! first set check if the endpoints shall always be included


call InitializeRND(GA,printseed=.true.)          !! initialize random number generator
!DEBUG
!call InitializeRND(GA,setseed=-31114,printseed=.true.)
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print *, "WARNING: SEED FIXED"
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!END DEBUG

fitstrucnumb = CE%totalstructures - CE%leftoutStr
nrdftgst = size(CE%dftgstlist)

! If we don't want to include the "endpoints", i.e. the pure elements, in the predictions set, we have
! to know what those endpoints are in the first place. Here: determine the "excludeList" which contains the 
! indices of those structures that are to be excluded from prediction (i.e. they are always in the 
! fitting set).
if (.not. GA%includeend) then
  ! get the pure elements of CE 1 (do we want to generalise this for arbitrary couplings?)
  write(*,'(/,A,/)') "Getting the pure elements..."   ! without that output, I got a SEGFAULT with 
                                                      ! ifort, OMP turned on and seed fixed (see above)
  if (CE%ncouplings==1) then
    allocate(pureConditions(1,1))
    pureConditions = 0
  elseif (CE%ncouplings==2 .and. CE%adsorbate%isAds) then
    allocate(pureConditions(3,1))  ! 1 condition only:
    pureConditions(:,1) = (/CE%adsorbate%freeSiteRank,2,1/) ! concentration fSR of CE 2 must be 0
  else
    stop "ERROR: subset setup only for non-coupled CEs or coupled adsorbate system CE"
  endif

  call get_pureElements(1,str,CE%CErank,CE%ncouplings,excludeList,pureConditions,CE%xeps)
          ! excludeList contains the pure elements, which are to be excluded from the prediction subsets

  indx => null()                  ! sigh, sometimes the "if(associated) deallocate" in heapsort throws an error :-/
  call heapsort(indx,excludeList) ! sorting them is essential for later correcting the indices of the 
                                  ! generated subsets
else
  stop "ERROR: you should not include the pure elements for the prediction"
endif

!================================================================================
! now: setup subsets
!================================================================================
print *
write(*,'(A)') "setting up subset lists:" 

nPop= GA%nPopulations   ! number of populations
if (GA%CrossValidationLeaveOneOut) then   ! setup the K and Rep for the LOO case
  GA%CrossValidationK   = fitStrucNumb-size(excludeList)
  GA%CrossValidationRep = 1
endif
nSS = GA%CrossValidationK * GA%CrossValidationRep   ! total number of subsets
allocate(GA%subset(0:nSS,1:nPop)) ! subsets 1...nSS contain the "regular" subsets, i.e. fitting structures and
                           ! prediction structures.
                           ! subset 0 contains ALL input structures (that are not excluded explicitely from
                           ! the fit). This subset will be used to call clusterexpand_subset in order to
                           ! replace clusterexpand_qp
do iPop=1,nPop; allocate(GA%subset(0,iPop)%predList(0)); enddo


do iPop=1,nPop ! set up different subsets for the different populations!

   print *
   write(*, '(A,I4)') "- population: ", iPop

   ! (1) setup the predList and fitList for all subsets. We don't take into account the structures that
   ! are going to be excluded from the prediction set (see above, pure structures, e.g.). Therefore, the 
   ! resulting lists "predList" and "fitList" don't know about those structures.
   call make_fitting_sets(fitStrucNumb-size(excludeList), &      ! number of structures available (also see below @ (2))
                          GA%CrossValidationK, &                 ! number of subsets of ONE repetition \   multiply for the
                          GA%CrossValidationRep, &               ! repetition                          /   total number of subsets
                          GA%subset(1:,iPop), &                  ! subsets to generate
                          GA)                                        
   ! (2) correct the lists for the indices of the excluded structures
   do iSS=1,nSS
     do i=1,size(excludeList)
       ! for each subset, for each excludeList entry: the indices of the structures need to be adjusted
       where(GA%subset(iSS,iPop)%predList >= excludeList(i))
         GA%subset(iSS,iPop)%predList = GA%subset(iSS,iPop)%predList + 1
       end where
       where(GA%subset(iSS,iPop)%fitList >= excludeList(i))
         GA%subset(iSS,iPop)%fitList = GA%subset(iSS,iPop)%fitList + 1
       end where
     enddo
     
     ! for each subset: the prediction-excluded structures are appended to the fitting list
     i = size(GA%subset(iSS,iPop)%fitList)
     GA%subset(iSS,iPop)%fitList => ralloc( GA%subset(iSS,iPop)%fitList , i+size(excludeList))
     GA%subset(iSS,iPop)%fitList(i+1:) = excludeList(:)
   
     write(*,'(A,I4,A)') "  subset #", iSS, " will predict structures #:"
     write(*,*) GA%subset(iSS,iPop)%predList
   enddo ! iSS
   
   ! setup the special subset "0":
   allocate(GA%subset(0,iPop)%fitList(fitstrucnumb))
   GA%subset(0,iPop)%fitList = (/(i,i=1,fitstrucnumb)/)    ! subset 0 contains ALL input structures
   !================================================================================
   
   !================================================================================
   ! Determine further properties of the structures in the subset
   !
   !! Find the GSs in each subset (so that we don't have to calculate the GSL over and over 
   !! again during the GA) 
   !!   AND
   !! Find the distance of each structure in the subset to the GSL hyperplanes
   !!   AND
   !! Get the mixture of GSs
   !!   AND
   !! Get the low energy structures and the distances to them
   !!   AND
   !! sets the GC deltas for each subset and each structure
   
   
   do iSS=0,nSS ! k=1,n
   
     GA%subset(iSS,iPop)%nFit   = size(GA%subset(iSS,iPop)%fitList)
     GA%subset(iSS,iPop)%nPred  = size(GA%subset(iSS,iPop)%predList)
   
     ! call the get_GSs routine. This routine invokes the qhull library and writes information
     ! to the file "convex_hull.out".
     ! parameters:
     ! -> the structures in the current subset
     ! -> the dimension of the E-x-diagram = rank of CE (always?)
     ! <- a list containing the indices of the groundstates in the subset
     ! <- a list containing the distance of each structure in the subset to the convex hull (GSL)
     ! <- a list containing the indices of the gs-structures in the subset that span the convex hull
     !    below each subset structure ("adjacent" gs-structures)
     call get_GSs(str(GA%subset(iSS,iPop)%fitList(:)),CE%CErank,CE%ncouplings,CE%xeps,GSlist,DistList,AdjStr)
     if (iSS==0) str(:)%distanceToInputGSL = DistList
   
     ! get the pure elements of CE 1 (do we want to generalise this for arbitrary couplings?)
     if (allocated(pureConditions)) deallocate(pureConditions)
     if (CE%ncouplings==1) then
       allocate(pureConditions(1,1))
       pureConditions = 0
     elseif (CE%ncouplings==2 .and. CE%adsorbate%isAds) then
       allocate(pureConditions(3,1))  ! 1 condition only:
       pureConditions(:,1) = (/CE%adsorbate%freeSiteRank,2,1/)! concentration fSR of CE 2 must be 1
     else
       stop "ERROR: subset setup only for non-coupled CEs or coupled adsorbate system CE"
     endif
     call get_pureElements(1,str(GA%subset(iSS,iPop)%fitList(:)),CE%CErank,CE%ncouplings,pureList,pureConditions,CE%xeps)
   
     ExDim = (CE%CErank-1)*CE%nCouplings+1 ! see comment in get_GSs for this dimension
   
     GA%subset(iSS,iPop)%nGS       = size(GSlist)
       allocate(GA%subset(iSS,iPop)%GSList(size(GSList)))
       allocate(GA%subset(iSS,iPop)%pureList(CE%CErank))
       allocate(GA%subset(iSS,iPop)%GSofPlane(ExDim,GA%subset(iSS,iPop)%nFit))
     GA%subset(iSS,iPop)%GSList    = GSlist    ! this is the list of ground states WITHIN the fitlist (so an entry 12 means: on the
                                        ! 12th position IN FITLIST there is a pure element.
     GA%subset(iSS,iPop)%pureList  = pureList  ! this is the list of pure elements WITHIN the fitlist (see comment above)
     GA%subset(iSS,iPop)%GSofPlane = AdjStr    ! dito!
   
   
     call get_GSLmixture(str(GA%subset(iSS,iPop)%fitList(:)), AdjStr, CE%CErank, CE%ncouplings, mixture)
       allocate(GA%subset(iSS,iPop)%GSmixture(ExDim,GA%subset(iSS,iPop)%nFit))
     GA%subset(iSS,iPop)%GSmixture =mixture
   
     ! Get the distance of each structure to the fit-energy convex hull (well, the vertices
     ! of the convex hull remain the same, i.e. the vertices are the vertices of the real energy
     ! convex hull, but the energies and thus the distance to this hull are calculated by the
     ! fit-energies.)
     if (CE%nreferences>0) then
       deallocate(DistList)
       call adjust_GSLdist(str(GA%subset(iSS,iPop)%fitList(:)), AdjStr, mixture, DistList)
     endif
       allocate(GA%subset(iSS,iPop)%GSLdist(size(DistList)))
     GA%subset(iSS,iPop)%GSLdist = DistList
       deallocate(GSlist,pureList,DistList,AdjStr,mixture)
   
     ! Get the true energy distance to structures that are lowest in energy at a specific concentration.
   ! (will do that separately)  ! Then, set the GC deltas according to a Boltzmann-factor
     call get_lowEstrucs(str(GA%subset(iSS,iPop)%fitList(:)), CE%xeps, StrLowE, DistList, useFitEnergy=.false.)
   !  call set_GCdeltas(DistList, CE, 1, deltaList)  ! 1 = current Delta index. At first, we use the first Deltas
       allocate(GA%subset(iSS,iPop)%StrDeltaDist(GA%subset(iSS,iPop)%nFit))
       allocate(GA%subset(iSS,iPop)%Delta(3,GA%subset(iSS,iPop)%nFit,CE%nObservables))
     GA%subset(iSS,iPop)%StrDeltaDist = DistList
   !  GA%subset(iSS,iPop)%Delta        = deltaList  ! { delta#, struc# }
   
   !    deallocate(deltaList)
   
     ! Get the fit energy distance to structures that are lowest in energy at a specific concentration.
     ! Store this information for later use (GC constraints).
     if (CE%nreferences>0) then
       deallocate(StrLowE, DistList)
       call get_lowEstrucs(str(GA%subset(iSS,iPop)%fitList(:)), CE%xeps, StrLowE, DistList, useFitEnergy=.true.)
     endif
       allocate(GA%subset(iSS,iPop)%StrLowE(size(StrLowE)))
       allocate(GA%subset(iSS,iPop)%StrLowEDist(size(DistList)))
     GA%subset(iSS,iPop)%StrLowE     = StrLowE      ! for each struc: index of struc with least energy
     GA%subset(iSS,iPop)%StrLowEdist = DistList     ! for each struc: distance to that low energy struc
       deallocate(StrLowE,DistList)
     
   enddo ! iSS (subsets)

enddo ! iPop (populations)

END SUBROUTINE generate_subsets
!***************************************************************************************************



!***************************************************************************************************
! Purpose:
! If we use a reference energy then the GC constraints use the distance of
! each structure to the fit-energy convex hull (GSL). The "convex hull" itself is comprised
! of the same vertices as the energy-(or formation energy)-convex hull. This means, we
! fix those vertices, even if the fit-energy convex hull would be different. Now, we have
! to determine the fit-energy distance to the GSL, wherefore we employ the mixture list.
!
! --> IN:
!            str : structures
!        adjlist : for each structure: the index of those structures that build the real
!                  convex hull plane, above which the structure is floating
!        mixture : mixture of the adjacent structure's concentrations that
!                  gives the original structure's concentration.
!                  (see subroutine get_GSLmixture for details)
!
! <-- OUT:
!           dist : for each structure: the distance to the convex hull in terms
!                  of fit-energy
!
subroutine adjust_GSLdist(str, adjlist, mixture, dist)
type(crystal), intent(IN) :: str(:)
integer, intent(IN)       :: adjlist(:,:)  ! { idx adj str, struc }
real(dp), intent(IN)      :: mixture(:,:)  ! { idx adj str, struc }
real(dp), pointer         :: dist(:)

integer iStr, nStr
real(dp) :: EfitConvexHull

nStr=size(str)
allocate(dist(nStr))

do iStr=1,nStr
!forall( iStr=1:nStr )
  dist(iStr) = str(iStr)%fitenergy - sum(str(adjlist(:,iStr))%fitenergy * mixture(:,iStr))
!!  print *, "##### ",iStr
!!  print *, "fE",str(iStr)%fitenergy
!!  print *, "adj",adjlist(:,iStr), str(adjlist(:,iStr))%fitenergy
!!  print *, "mix",mixture(:,iStr)
!!  print *, "d", dist(iStr)
!end forall
enddo

end subroutine adjust_GSLdist
!***************************************************************************************************






!***************************************************************************************************
! Purpose :
! - calculate AND set the proper structure-dependent Garbulsky-Ceder constraints for all subsets
! - this includes Boltzman smearing, groundstate weights and an extra GC1 for the pure elements
!
! IN -->
!         deltaPure    :   the 1st delta (total deviation) for the pure elements (as supplied in fitpar.in)
!             delta    :   the GC deltas as they were read in from fitpar.in
!                          1st index: the constraint, i.e. 1 for the 1st delta, 3 for the 3rd
!                          2nd index: what "current index" tells which set of deltas we are
!                                     using at the moment. At first, we start at 1, then---after
!                                     swtching the deltas---this index is 2, and so on.
!                          3rd index: what observable are we looking at? 1st observable is energy.
!           current    :   what set of deltas are we using? (see delta for more comments on that)
!                kT    :   the Boltzmann kT factor that smears the deltas according to the 
!                          structures' distance from the low stuctures.
!          GSweight    :   ground-state weight. Actually, it's rather a _factor_ (not a weight)
!                          by which the deltas are multiplied for the ground states.
!           subsets    :   the subsets (all, starting at 0) that were setup in generat_subsets.
!                          This includes the distances of all structures to the low structures,
!                          as well as the deltas for all structures.
!          
! <-- OUT
!           subsets    :   The deltas for all structures are set.
!
subroutine set_GCdeltas(deltaPure, delta, current, kT, GSweight, subsets)
real(dp), intent(in)          :: delta(:,:,:) ! { constraint#, current idx, observable idx }
real(dp), intent(in)          :: deltaPure
integer, intent(in)           :: current
real(dp), intent(in)          :: kT
type(t_subset), intent(inout) :: subsets(0:) ! you HAVE TO TELL fortran that the subsets start at 0!!
real(dp), intent(in)          :: GSweight

integer iSS, nSS, iFitStr, iDelta


!write(*,'(A$)') "entering set_GCdeltas"

nSS=size(subsets)
!print *, "nss=",nss

do iSS=0,nSS-1  ! subset 0 contains ALL structures, 
                ! all other subsets are TRUE SUB-SETS (mathematically) of subset 0

   forall(iFitStr=1:subsets(iSS)%nFit, iDelta=1:3) ! apply Boltzman correction
      subsets(iSS)%Delta(iDelta, iFitStr, :)= delta(iDelta, current, :) * &
                                              exp( subsets(iSS)%StrDeltaDist(iFitStr) / kT )  ! all deltas are smeared wrt their (energy!)-distance
   end forall 


   ! For first observable (energy), additional restrictions apply:
   ! (A) be extra picky at all groundstates
   forall(iDelta=1:3) 
      subsets(iSS)%Delta(iDelta, subsets(iSS)%GSList, 1)   = delta(iDelta,current,1)*GSweight 
   end forall

   ! (B) be super picky at the pure elements
   subsets(iSS)%Delta(1, subsets(iSS)%pureList, 1) = deltaPure 

enddo
!write(*,*) "done"

end subroutine set_GCdeltas
!***************************************************************************************************




!***************************************************************************************************
! Purpose:
! Given a set of structures and "adjacent" (GS) structures, determine the mixture
! of the GS concentrations that will give the actual structure's concentration.
!
! --> IN: 
!              str : structures (probably a subset from generate_subsets)
!          adjlist : adjacent structures (see get_GSs routine for explanation)
!                k : dimension (CE rank)
!               nC : number of CE couplings
! <-- OUT:
!          mixture : mixture of the adjacent structure's concentrations that
!                    gives the original structure's concentration.
!
!                    Example (binary):
!                    Suppose structure B, conc = xb, is above the GSL. The
!                    left GS is A, conc = xa, the right GS is C, conc = xc.
!                    Then, we search for the mixture factors ("probabilities")
!                    pa and pc for which:
!                              pa xa + pc xc == xb
!                              pa + pc == 1
!
!                    For solutions to this set of equations in arbitrary
!                    dimensions, try the following Mathematica script.
!                     
! Mathematica script to determine the solution:
!
!    For[n = 2, n <= 3, n++,
!     For[komp = 1, komp <= n - 1, komp++,
!      eq1[komp] = \!\(
!    \*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(p[i] 
!           x[i, komp]\)\) == c[komp];
!      ];
!     eq2 = \!\(
!    \*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(p[i]\)\) == 1;
!     eq1Table = Table[eq1[komp], {komp, 1, n - 1}];
!     AllEq = AppendTo[eq1Table, eq2];
!     sol = Solve[AllEq, Table[p[i], {i, 1, n}]];
!     sol = FullSimplify[sol];
!     Print[n, sol]
!     ]
!
!
! !!!!! note for future improvements !!!!!
! This is a routine that determines the weight factors for a linear interpolation.
! tk has programmed something *very* similar for grand-canonical MC with varying 
! mu parameters. There it was solved in a general way, without resorting to
! actually solving those equations. tk thinks we should also use that way of
! determining the weight parameters here. So, have a look at the Monte Carlo
! routines, in particular at linear_interpolation.
!
!
subroutine get_GSLmixture(str,adjlist,k,nC,mixture)
type(crystal), intent(IN) :: str(:)
integer, intent(IN)       :: adjlist(:,:)  ! { index adj str, struc }
integer, intent(IN)       :: k, nC         ! rank of CE, number of CE couplings
real(dp), pointer         :: mixture(:,:)  ! { index adj str, struc }
integer                   :: dim  ! dimension of the GS diagram (see comment at subroutine get_GSs)

real(dp), allocatable :: c(:)
real(dp), allocatable :: x(:,:) ! { #adj str, #conc }
integer :: iStr, nStr, iAdj, nAdj, i
integer, allocatable  :: IX(:) ! independent concentrations
real(dp) :: numerator, denominator
logical :: isAboveVertex

integer, allocatable :: adjacent(:)

! dimension of the GS: diagram
dim = (k-1)*nC + 1

nStr=size(str)
nAdj=dim

allocate(c(k*nC))
allocate(x(dim,k*nC))
allocate(mixture(dim,nStr))
allocate(adjacent(dim))

do iStr=1,nStr
  isAboveVertex = .false.

  c(:) = reshape(str(iStr)%x(:,:),(/1/)) ! put the coupled concentrations of str(iStr) into a single 
                                         ! concentration vector
  
  adjacent = adjlist(:,iStr)  ! store the adjacent structures of strutures iStr (intent IN)...

  if (all(adjacent==adjacent(1))) then
    mixture(:,iStr) = 1._dp / dim
    cycle
  endif

  do iAdj=1,nAdj
    x(iAdj,:)= reshape(str(adjacent(iAdj))%x(:,:),(/1/))  
                                            ! ... and get the concentrations of those adjacent structures.
                                            ! IMPORTANT: I do not like this double-step assignement. However,
                                            ! the _single_ line x(iAdj,:)=str(adjlist(iAdj,iStr))%x(:) threw
                                            ! a SEGFAULT in case of ifort 10.1 and >200 input structures.
                                            ! To circumvent this, you have to touch "adjlist" before the loop.
                                            ! Either you use a dummy variable for that, or you do what I (tk) did here.
  enddo

  do i=1,k*nC
    if (all(x(:,i)==x(1,i))) then; isAboveVertex=.true.; else; isAboveVertex=.false.; exit; endif
  enddo

  if (allocated(IX)) deallocate(IX)

  if (isAboveVertex) then
    mixture(:,iStr) = 0._dp
    mixture(1,iStr) = 1._dp
  else ! .not. isAboveVertex
    select case(dim)
    !--------------------------------------------------------------------------------
    case(2) ! binary CE, no coupling (=> all concentrations have coupling index 1 (the last index))
      allocate(IX(1)); IX = 1

      denominator = x(1,IX(1))-x(2,IX(1))

      numerator   =  c( IX(1) )-x(2,IX(1))
      mixture(1,iStr) = numerator / denominator

      numerator   = -c( IX(1) )+x(1,IX(1))
      mixture(2,iStr) = numerator / denominator
    !--------------------------------------------------------------------------------
    case(3)
      allocate(IX(2)) ! there exist 2 independent concentrations

      select case(nC)
      case(1) ! one ternary CE
        IX(1) = 1; IX(2) = 2;
      case(2) ! two coupled binary CEs
        IX(1) = 1; IX(2) = 3;  ! the independent concentrations are concentration #1 and #3.
                               ! #2 is only the inverse of #1, #4 the inverse of #3
      case default
        print *, "ERROR: get_GSLmixture is not setup for (k,couplings)=",k,nC
        stop
      end select
       
      denominator = x(1, IX(2) ) * (  x(2, IX(1) )-x(3, IX(1) ) ) &
                  + x(2, IX(2) ) *    x(3, IX(1) ) &
                  - x(2, IX(1) ) *    x(3, IX(2) ) &
                  + x(1, IX(1) ) * ( -x(2, IX(2) )+x(3, IX(2) ) )

      numerator =   c( IX(2) ) * ( x(2, IX(1) ) - x(3, IX(1) ) ) &
                + x(2, IX(2) ) * ( -c( IX(1) )  + x(3, IX(1) ) ) &
                + x(3, IX(2) ) * ( c( IX(1) )   - x(2, IX(1) ) ) 
      mixture(1,iStr) = numerator / denominator

      numerator = x(1, IX(2) ) * ( c( IX(1) )   - x(3, IX(1) ) ) &
                +   c( IX(2) ) * ( -x(1, IX(1) )+ x(3, IX(1) ) ) &
                + x(3, IX(2) ) * (-c( IX(1) )   + x(1, IX(1) ) )
      mixture(2,iStr) = numerator / denominator

      numerator =   c( IX(2) ) * ( x(1, IX(1) ) - x(2, IX(1) ) ) &
                + x(1, IX(2) ) * ( -c( IX(1) )  + x(2, IX(1) ) ) &
                + x(2, IX(2) ) * ( c( IX(1) )   - x(1, IX(1) ) )
      mixture(3,iStr) = numerator / denominator

    !--------------------------------------------------------------------------------
    case default
      print*, "ERROR: get_GSLmixture is not setup for dimension ", dim
      mixture(:,:)=0
      !stop
    !--------------------------------------------------------------------------------
    end select
  endif ! isAboveVertex

enddo ! iStr

deallocate(c,x)

end subroutine get_GSLmixture
!***************************************************************************************************




!****************************************************************************************************
! Purpose:
! Make a cluster-expansion fit of the input structures' obersvables. Usually, a genetic algorithm
! is applied to find the best genome. However, the routine is also used for a simple fit, where
! the genetic algorithm itself is basically deactivated by only having 1 generation and a full fit.
!
! Remark: 
! - no lockouts are implemented. Instead, individuals can die after a certain number of generations
!
!
! --> IN:
!                GA : variables for the genetic algorithm, including the subsets
!                     (also intent(out) because of the seed and other small settings)
!                CE : general variables for the cluster expansion, including fitting deltas (Garbulsky-Ceder)
!                     (also intent(out) because some fitting parameters are set right in front of the quadprg routine)
!          clusters : all clusters that are used for the CE (including empty cluster and onsites)
!        structures : input structures
!              iObs : the current observable
! <-- OUT:
!           results : information about the GA results, including the best individual, the CVS and the J's
!
subroutine geneticAlgorithm_for_clusterExpansion(GA, CE, clusters, structures, iObs, results)
! Interface:
type(GA_variables), intent(inout) :: GA
type(CE_variables), intent(inout) :: CE
type(FigRep), intent(in)          :: clusters(0:)       ! { cluster# }
type(crystal), intent(in)         :: structures(:)
integer, intent(in)               :: iObs               ! the current observable
type(t_fitResults)                :: results            ! the final results of the fitting

! Population for the GA:
type(t_individual), pointer :: individuals(:,:)           ! { individual#, population# }           
type(t_individual), pointer :: bestDeadIndividual(:,:)    ! { (1, population#) } sorry, must be a 1d array in order to use all the nice functions...
real(dp)          , pointer :: avgDeviation(:)            ! { population# }
integer                     :: bestIndividualIndex, bestPopulationIndex

! Children generation
real(dp), pointer           :: roulette(:,:)            ! { individual#, pop# }
type(t_GA_parents), pointer :: parents(:,:)             ! { child#, pop# }
type(t_individual), pointer :: children(:,:)            ! { child#, pop# }
integer                     :: iChild, nChild

! population
integer iPop, jPop
integer :: nPop
integer rndPop
type(t_individual), pointer :: tmpIndividual(:)
logical :: virginBirth(3)

! counters
integer iGeneration, nGeneration
integer nStr, iStr, iCl
integer iSS, nSS
integer iIndv, nIndv
integer, allocatable :: currDeltaIdx(:)  ! { pop# }
integer iRepeat
integer generationLastDiversityFail

! temporary
real(dp), pointer :: PI(:,:)    ! storage for the (activated clusters') PIs of all structures
integer failIndividual          ! what indiviual failed? (not used at the moment)
real(dp) :: random
character(10), pointer :: diversityCountermeasure(:)
real :: NaN = 0.                ! want to produce NaN by division 0./0., which gfortran doesn't like if it's not obscured from the compiler

! ETA
real(dp) walltime_curr, cputime_curr
real(dp) walltime_prev, cputime_prev

! flags
logical diversityFail


!----------------------------------------------------------------------------------------------------
! Initialisation

GA%nTotClusters = size(clusters)  ! how many clusters in total? (constant and inequivalent sites are now INCLUDED in the clusters)
GA%nAlwaysOn    = count(clusters(:)%nV<=1)   ! constant term + onsite terms are supposed to be always activated
                                             ! ( Originally, we wanted to give every single cluster free for the GA to switch on or
                                             ! off. But we found out that: 1) it's really good to always have the constant term on
                                             !                             2) in the LaH multilattice case switching off the onsites
                                             !                                proved fatal for the GSS-PREDICTIONS, but NOT for the
                                             !                                CE vs. DFT values. So, the CVS was really good, but the
                                             !                                GSS predictions were terrible.
                                             !                                In monolattice cases, we didn't observe that behaviour. 
                                             !                                I (tk) guess---and this makes sense IMHO---that the 
                                             !                                onsites are really (and only) important if there are
                                             !                                inequivalent onsites!
                                             ! )
nStr             = size(structures(:))       ! number of structures (total)
nSS              = size(GA%subset(:,1))-1    ! subset(0) are ALL structures, i.e. it is not a true subset ==> subtract 1 in order to get only true subsets
nPop             = GA%nPopulations           ! number of populations
nIndv            = GA%nIndividuals           ! number of indiviuals
nGeneration      = GA%maxGenerations         ! number of generations
nChild           = GA%nChildren              ! number of children
GA%nPairClusters = count(clusters(1:)%nV==2) ! number of pair clusters in the cluster pool

generationLastDiversityFail = 0

allocate( individuals(nIndv, nPop) )
allocate( tmpIndividual(1) )
allocate( bestDeadIndividual(1, nPop) )
allocate( children(nChild,nPop)   ) 
allocate( roulette(nIndv, nPop), parents(nChild, nPop) )
allocate( diversityCountermeasure(nPop) )
allocate( avgDeviation(nPop) )
allocate( currDeltaIdx(nPop) )

currDeltaIdx= 1                      ! we start with the first GC-deltas that are supplied

!----------------------------------------------------------------------------------------------------

!====================================================================================================
! 1st generation
iGeneration=1

! setup a temporary individual
call setup_population( tmpIndividual(:),1 , GA%nMaxClusters, GA%nTotClusters, GA%nAlwaysOn, GA, warn=.false.)

! go through all populations and setup the individuals
do iPop=1,nPop
   call set_GCdeltas( CE%pureDelta1 , CE%delta , currDeltaIdx(iPop) , CE%kTfit , CE%gstWeight , GA%subset(0:,iPop) )
   
   ! setup the population:
   if (GA%simpleFit) then ! if we don't really want to run a genetic algorithm, but only make a "simple fit", then we want to
                          ! activate the first GA%nMaxClusters genes
     call setup_population( individuals(:,iPop)          , nIndv , GA%nMaxClusters, GA%nTotClusters, GA%nAlwaysOn, GA, &
                            warn=.false., activateFirstGenes=.true.)
   
   else ! we are running a genetic algorithm. Setup diverse population
   
     ! setup a new population
     call setup_population( individuals(:,iPop)          , nIndv , GA%nMaxClusters, GA%nTotClusters, GA%nAlwaysOn, GA, warn=.true.)
   
     ! correct the pair selection:
     ! the old-UNCLE scheme used to select a certain _number_ of pairs, n_pair, and then activated the first n_pair pairs,
     ! instead of activating a random selection of the total set of pairs. If you want to use the old procedure, uncomment
     ! the next line and all other calls to correct_pairs subroutine!
     call correct_pairs(GA%selectSmallestPairsOnly, GA%nPairClusters, individuals(:,iPop), GA%nAlwaysOn)
   
     ! setup the "dead individual": only 1 individual, with no activated clusters
     call setup_population( bestDeadIndividual(:,iPop)   , 1     , 0              , GA%nTotClusters, 1           , GA, warn=.false.)
   
   endif
   
   ! Do we read in a population?
   ! This will overwrite all individuals that were setup before. (The setup was necessary, however, for the allocation of memory)
   if (GA%continue%yes) call read_population( individuals(:,iPop), bestDeadIndividual(1,iPop), iPop, iObs, GA%continue ) 
   
   ! tell the user
   diversityCountermeasure = ""
   call output__population_info(iPop)
enddo ! iPop

call output__GCdeltas(currDeltaIdx(1))

write(*,'(//,A)') "Making initial fit for all populations. This can take a while..."
!$OMP PARALLEL DO private(failIndividual)
do iPop=1,nPop
   ! make the initial fit
   call make_ceFit_and_calculate_CVS_for_individuals( CE, structures, GA%subset(1:,iPop), individuals(:,iPop), iObs, clusters(0:)%damping, failIndividual)  
   if (GA%continue%yes .and. bestDeadIndividual(1,iPop)%cvs /= -1) &
        call make_ceFit_and_calculate_CVS_for_individuals( CE, structures, GA%subset(1:,iPop), bestDeadIndividual(:,iPop), iObs, clusters(0:)%damping, failIndividual)
!$OMP CRITICAL
   write(*,'(A,I3,A)') "Fit for population ",iPop, " done."
!$OMP END CRITICAL
enddo ! iPop
!$OMP END PARALLEL DO

do iPop=1,nPop
   ! will call the diversity routine only to set the correct value for avgDeviation
   call check_population_diversity( individuals(:,iPop), GA, avgDeviation(iPop), diversityFail ) 
   
   ! tell the user about the first generation
   if (.not. GA%simpleFit) then
     if (iPop==1) call output__generation_info_header
     call output__generation_info(1,iPop)
   endif
enddo ! iPop

!====================================================================================================
! all other generations
do iGeneration=2,nGeneration

!$OMP PARALLEL DO private(iChild, iIndv, random, bestIndividualIndex, virginBirth, failIndividual)
! Note about the OMP parallelisation:
! I (tk) have run into several race conditions that I couldn't resolve, as well as data corruption that must be due
! to some data interaction between the OMP threads that I didn't consider. So, I decided to basically run the whole
! loop in serial mode (via the OMP CRITICAL pragma) and only do the fitting in parallel. Since this is the place
! were most time is spent, this dumb-and-dull parallelisation should already suffice to bring about a substantial
! increase in speed.
   do iPop=1,nPop

      ! Do we want to switch the Garbulsky-Ceder deltas? ..................................................
      if ( any(CE%gen4DeltaSwitch == iGeneration) ) then
!$OMP CRITICAL

        ! even if currDeltaIdx is always the same for all population, we're running into trouble in case
        ! of OMP parallelisation. Therefore, currDeltaIdx has an index that designates it to the 
        ! corresponding population.
        if (iPop==1) then
          call output__switch_GCdeltas(currDeltaIdx(iPop), currDeltaIdx(iPop)+1)    ! tell the user what we're doing
        endif
        currDeltaIdx(iPop) = currDeltaIdx(iPop) + 1
        call set_GCdeltas( CE%pureDelta1 , CE%delta , currDeltaIdx(iPop) , CE%kTfit , CE%gstWeight , GA%subset(0:,iPop) )  ! set the deltas for all subsets
        ! for different deltas, the fit changes. Therefore, re-fit and re-evaluate the fitness of the individuals
        ! before continuing with the new deltas that---maybe---won't return a valid fit
        ! anymore, we're writing a backup file for the previous generation:
        call save_population( individuals(:,iPop), bestDeadIndividual(1,iPop), iPop, iGeneration-1, iObs )  
!$OMP END CRITICAL
        ! make new fit:
        call make_ceFit_and_calculate_CVS_for_individuals( CE, structures, GA%subset(1:,iPop), individuals(:,iPop)       , iObs, clusters(0:)%damping, failIndividual)
        if (bestDeadIndividual(1,iPop)%cvs /= -1) &
           call make_ceFit_and_calculate_CVS_for_individuals( CE, structures, GA%subset(1:,iPop), bestDeadIndividual(:,iPop), iObs, clusters(0:)%damping, failIndividual)
!$OMP BARRIER
      endif
      !....................................................................................................

!$OMP CRITICAL  
        
      ! complete the setup of children. Those pointers will be nullfied after the good children have replaced
      ! some individuals of the parent generation. So we cannot allocate them above together with children(:) itself.
      do iChild=1,nChild
        allocate(children(iChild,iPop)%genome(0:GA%nTotClusters-1))
        allocate(children(iChild,iPop)%onGenes(1:GA%nTotClusters))
      enddo

      ! The CVS represents the fitness of an individual. The smaller the CVS, the fitter the individual.
      ! In the following, we want to choose parents for the next generation according to their fitness: fitter individuals
      ! are more likely to be parents than weaker ones. 
      ! Therefore, we translate the CVS to a roulette wheel: According to their fitness, the single individuals take more or
      ! less room on that wheel.
      call map_CVS_to_roulette( individuals(:,:)%CVS, roulette(:,:) )
   
      ! Select parents:
      ! The roulette wheel is "rotated" and where the roulette wheel "stops" determines the parent. Due to the construction
      ! of the roulette wheel, fitter parents are selected more often.

      call make_parent_list(roulette(:,:), parents(:,iPop), iPop, GA)
   
      ! Beget children
      ! nChild children are created from those parents that were selected in make_parent_list.
      diversityCountermeasure(iPop) = ""
      do iChild=1,nChild
        if (parents(iChild,iPop)%father==0 .and. parents(iChild,iPop)%mother==0) then
          print *, "cvs=",individuals(:,iPop)%cvs
          print *, "roulette=",roulette(:,:)
          stop "ERROR: bad programming while begetting children."
        endif
        call mate_2_parents( individuals( parents(iChild,iPop)%father, parents(iChild,iPop)%fpop ), &
                             individuals( parents(iChild,iPop)%mother, parents(iChild,iPop)%mpop ), &
                             children(iChild,iPop), GA)
        if     (parents(iChild,iPop)%fpop/=iPop .and. parents(iChild,iPop)%mpop/=iPop) then; diversityCountermeasure(iPop) = trim(diversityCountermeasure(iPop))//"O";
        elseif (parents(iChild,iPop)%fpop/=iPop .or.  parents(iChild,iPop)%mpop/=iPop) then; diversityCountermeasure(iPop) = trim(diversityCountermeasure(iPop))//"A"; 
        endif
      enddo
   
      ! Children are mutated. 
      call mutate_individuals( children(:,iPop), GA%mutationProbability, GA ) 
      ! Check if the mating of parents + mutation yielded a valid genome. If not, the child is invalid 
      ! and therefore reSetup ("virgin birth"):
      ! - genome with too many genes switched on
      call check_individuals_for_too_many_onGenes( children(:,iPop), GA, virginBirth(1) )
      ! - genome identical to a dead individual ("Lockout" the best one), any population
      call check_individuals_for_equality( children(:,iPop), GA, bestDeadIndividual(1,:), virginBirth(2) )
      ! - genome identical to any living individual, any population
      call check_individuals_for_equality( children(:,iPop), GA, reshape( individuals(:,:), (/ nIndv*nPop /)), virginBirth(3))
      if (any(virginBirth.eqv..true.)) then; diversityCountermeasure(iPop) = trim(diversityCountermeasure(iPop))//"V"; endif

!$OMP END CRITICAL
      
      ! Make a CE fit for the children and evaluate their CVS
      call make_ceFit_and_calculate_CVS_for_individuals( CE, structures, GA%subset(1:,iPop), children(:,iPop), iObs, clusters(0:)%damping, failIndividual)  
   
!$OMP CRITICAL

      ! make a new generation by replacing the worst parents with the best children
      call make_new_generation( individuals(:,iPop), children(:,iPop), GA, bestDeadIndividual(1,iPop) )
      call save_population( individuals(:,iPop), bestDeadIndividual(1,iPop), iPop, 0, iObs )   ! write a snapshot file for the individuals
                                                                         ! iGeneration=0 means: current generation, will be replaced after the
                                                                         ! the next step. "Real" backup files that are not replaced are only written
                                                                         ! when deltas change, or when population is resetup
   
      ! see if the population is diverse enough. 
      call check_population_diversity( individuals(:,iPop), GA, avgDeviation(iPop), diversityFail ) 
   
      ! manipulate the population if
      ! - the diversity is too low (i.e. diversityFail is true)
      ! - and the last manipulation took place at least 10 generations ago
      if (diversityFail .and. iGeneration >= generationLastDiversityFail + 1) then
        generationLastDiversityFail = iGeneration
   
        ! We try to increase the diversity of the population by a mutation burst.
        ! Instead of increasing the mutation rate for all the subsequent generations, we
        ! simply take all (but the best) CURRENT individuals and mutate them with an increased 
        ! mutation probability
        random =  rnd(GA)
        if (random<GA%diversityEX) then
          diversityCountermeasure(iPop) = "EX" ! EXchange individuals
          call get_index_of_best_individual( individuals(:,iPop) , bestIndividualIndex )
          do iIndv=1,nIndv
            if (iIndv==bestIndividualIndex) cycle
            rndPop = ceiling(rnd(GA)*nPop)
            if (rndPop /= iPop) then
              tmpIndividual(1) = individuals( iIndv, iPop)
              individuals( iIndv, iPop )  = individuals( iIndv, rndPop )
              individuals( iIndv, rndPop )= tmpIndividual(1)
            endif
          end do
        elseif (random<GA%diversityEX+GA%diversityMB) then
          diversityCountermeasure(iPop) = "MB" ! Mutation Burst
          call get_index_of_best_individual( individuals(:,iPop) , bestIndividualIndex )
          call mutate_individuals( individuals(:bestIndividualIndex-1,iPop), 10*GA%mutationProbability, GA ) ! mutation burst part 1 \  leave out (i.e. don't mutate) 
          call mutate_individuals( individuals(bestIndividualIndex+1:,iPop), 10*GA%mutationProbability, GA ) ! mutation burst part 2 /  the best individual
          call check_individuals_for_too_many_onGenes( individuals(:,iPop), GA, virginBirth(1) )
          call check_individuals_for_equality( individuals(:,iPop), GA, bestDeadIndividual(1,:), virginBirth(2) )
          call check_individuals_for_equality( individuals(:,iPop), GA, reshape( individuals(:,:), (/ nIndv*nPop /)), virginBirth(3) )
          call correct_pairs(GA%selectSmallestPairsOnly, GA%nPairClusters, individuals(:,iPop), GA%nAlwaysOn)       ! see comment above the first correct_pairs call
      
          call check_population_diversity( individuals(:,iPop), GA, avgDeviation(iPop), diversityFail ) ! call this routine only to get the avgDeviation right
          call make_ceFit_and_calculate_CVS_for_individuals( CE, structures, GA%subset(1:,iPop), individuals(:,iPop), iObs, clusters(0:)%damping, failIndividual)  
        else
          diversityCountermeasure(iPop) = "RS" ! ReSetup
          call resetup_population( individuals(:,iPop), GA ) 
          call correct_pairs(GA%selectSmallestPairsOnly, GA%nPairClusters, individuals(:,iPop), GA%nAlwaysOn)      ! see comment above the first correct_pairs call
          call make_ceFit_and_calculate_CVS_for_individuals( CE, structures, GA%subset(1:,iPop), individuals(:,iPop), iObs, clusters(0:)%damping, failIndividual)  
        endif
   
        call check_population_diversity( individuals(:,iPop), GA, avgDeviation(iPop), diversityFail ) ! call this routine only to get the avgDeviation right
   
   
      endif ! diversity fail
      call output__generation_info(iGeneration,iPop)
   
      ! nullifiy the children. THIS IS TERRIBLY IMPORTANT, as the "make_new_generation" routine copies "children to individuals" as a whole
      ! and not step by step, e.g. by copying "children%genome to individuals%genome" and so on. Therefore, we have to "detach" the individuals
      ! from their "from location", i.e. the children.
      ! If you don't do that, the individuals are incorrect.
      ! Don't deallocate the children's %genome or %onGenes, both of those continue to live in the individuals!
      do iChild=1,nChild; 
        children(iChild,iPop)%genome  => null()
        children(iChild,iPop)%onGenes => null()
      enddo;   

!$OMP END CRITICAL   

   enddo ! iPop
!$OMP END PARALLEL DO

end do ! Generation loop
! end all other generations
!====================================================================================================


! Having finished the GA loop, we now turn our attention to the best individual.


!----------------------------------------------------------------------------------------------------
! get best individual
!----------------------------------------------------------------------------------------------------
if (GA%simpleFit) then
  ! simple fit
  bestIndividualIndex = 1
  bestPopulationIndex = 1
else
  ! GA
  bestPopulationIndex=1
  do iPop=1,nPop
     call get_index_of_best_individual( individuals(:,iPop), bestIndividualIndex )
     if (.not.isnan( bestDeadIndividual(1,iPop)%CVS ) .and. bestDeadIndividual(1,iPop)%CVS/=-1) then
     if ( (bestDeadIndividual(1,iPop)%CVS < individuals( bestIndividualIndex, iPop )%CVS) & ! if a dead invidiual is better (wrt CVS) than the current best individual...
          .or. &                                                                            ! or...
          ( (bestDeadIndividual(1,iPop)%CVS == individuals( bestIndividualIndex, iPop )%CVS) .and. (bestDeadIndividual(1,iPop)%nOnGenes < individuals( bestIndividualIndex, iPop )%nOnGenes) ) &
          ! ... if a dead individual is equally good but needs fewer clusters than the current best individual, then ...
        ) then  
       individuals( bestIndividualIndex, iPop )  = bestDeadIndividual(1,iPop)    ! ... we simply replace the current best with the dead one
     endif
     if (individuals( bestIndividualIndex, iPop)%CVS < minval( individuals( :, :iPop-1)%CVS ) ) bestPopulationIndex = iPop
     endif
  enddo ! iPop
  call get_index_of_best_individual( individuals(:,bestPopulationIndex), bestIndividualIndex )
endif

results%bestIndividualIndex = bestIndividualIndex
results%bestPopulationIndex = bestPopulationIndex

!----------------------------------------------------------------------------------------------------
! Averaged J: Jav
! get averaged J (subsets) and determine final cvs
!----------------------------------------------------------------------------------------------------

! If you want to have the Jav without damping, uncomment the following lines:
! Jav w/o damping>call make_ceFit_and_calculate_CVS_for_individuals( CE, &
! Jav w/o damping>                                                   structures, &
! Jav w/o damping>                                                   GA%subset(1:,bestPopulationIndex), &
! Jav w/o damping>                                                   individuals(bestIndividualIndex:bestIndividualIndex,bestPopulationIndex), &
! Jav w/o damping>                                                   iObs, &
! Jav w/o damping>                                                   (/ (0.0_dp,iCl=0,GA%nTotClusters-1) /), &
! Jav w/o damping>                                                   failIndividual)  

! Set the best individual:
results%bestIndividual = individuals( bestIndividualIndex, bestPopulationIndex ) 

! allocate mem for the averaged Js and the prediction errors:
allocate(results%Jav( 0:results%bestIndividual%nOnGenes-1 ))
allocate(results%predictionError( maxval(GA%subset(:,results%bestPopulationIndex)%nPred),nSS) ); results%predictionError=0._dp

! average the J's of the subsets...
call average_J(results%bestIndividual%J(0:,:), results%Jav(0:), subsetMask=results%bestIndividual%subsetFitSuccessful )  

! ... and calculate the CVS for those averaged J's
allocate(PI(nStr, 0:results%bestIndividual%nOnGenes-1)) ! put the structures' PIs in PI
do iStr=1,nStr; PI(iStr,0:) = structures(iStr)%PI( results%bestIndividual%onGenes( 1:results%bestIndividual%nOnGenes ) ); enddo
results%CVS = calc_cvs_for_individual( &  ! <-- get CVS 
                                             structures(:)%fitenergy, &                      ! IN:  { observable# }
                                             GA%subset(1:nSS,results%bestPopulationIndex), & ! IN:  { subset# }
                                             (/(.true.,iSS=1,nSS) /), &                      ! IN:  want to take all subsets
                                             spread(results%Jav(0:),dim=2,nCopies=nSS), &    ! IN:  J for each subset, all are the same Jav!
                                             PI(:,0:), &                                     ! IN:  { structure#, on-cluster# } PI
                                             predictionErrors=results%predictionError(:,:) ) ! <-- also get the prediction errors
deallocate(PI)  ! this is strange: I have to deallocate PI here and allocate (+assign values) it again later on
                ! if I want the simple fit to work.

!----------------------------------------------------------------------------------------------------
! full fit J: Jff
! get full fit J (all structures, i.e. subset(0))
!----------------------------------------------------------------------------------------------------
allocate(results%Jff( 0:results%bestIndividual%nOnGenes-1 ))
! make a fit of subset 0
call make_ceFit_and_calculate_CVS_for_individuals( CE, structures, &
                                                   GA%subset(0:0,results%bestPopulationIndex), &
                                                   individuals(results%bestIndividualIndex:results%bestIndividualIndex,results%bestPopulationIndex), &
                                                   iObs, &
                                                   clusters(0:)%damping, &
                                                   failIndividual)        ! we 'misuse' the CVS function: putting in only 1 subset---namely
                                                   ! subset 0---, we get the full fit ECIs "Jff". All other returns from that function are not used.
results%Jff(0:) = individuals(results%bestIndividualIndex, results%bestPopulationIndex)%J(0:,1) ! the J's are the full fit J's. On the rhs, the 2nd index of J 
                                                   ! must be 1, since there was only 1 subset used for the fit.

! failsafe
if (size(individuals(results%bestIndividualIndex, results%bestPopulationIndex)%J,2)>1) &
     stop "ERROR: bad programming, something went wrong with the full fit J's"

!----------------------------------------------------------------------------------------------------
! was this a simple fit run? (i.e. no "real" genetic algorithm, when we simply 'misused' the routine)
!---------------------------------------------------------------------------------------------------- 
if (GA%simpleFit) then
  allocate(results%Jav( 0:results%bestIndividual%nOnGenes-1 ))
  results%Jav(0:) = results%Jff(0:)     ! no Jav exist, because no average is possible. However, since we always use Jav in the following,
                                        ! we set Jav = Jff
  results%CVS     = NaN / NaN           ! CVS is not defined
endif

!----------------------------------------------------------------------------------------------------
! get some fitting statistics
!----------------------------------------------------------------------------------------------------
! * fitting errors: error between dft and ce value for all structures
allocate(PI(nStr, 0:results%bestIndividual%nOnGenes-1)) ! put the structures' PIs in PI
do iStr=1,nStr; PI(iStr,0:) = structures(iStr)%PI( results%bestIndividual%onGenes( 1:results%bestIndividual%nOnGenes ) ); enddo
allocate(results%fittingError(nStr))
call get_errors_of_CE( structures(:)%fitenergy, &    
                       results%Jav(0:), &                  ! IN:  { on-cluster# } effective cluster interactions
                       PI( :,0:), &                        ! IN:  { structure#, on-cluster# }
                       results%fittingError(:) )           ! OUT: { structure# }
! * failsafe
if (all(results%Jav==0.) .and. all(results%Jff==0.)) then
  write(*,'(A,I4)') "ERROR: in observable #", iObs
  stop              "ERROR: fitting routines failed to produce an appropriate fit."
endif
if (any(  abs(results%fittingError(:)) > GA%subset(0,results%bestPopulationIndex)%Delta( 1,:,iObs ) )) then
  write(*,'(A)') "WARNING: some fitting errors are larger than the allowed Delta1. Check appropriate fittingErrors.out file."
endif

! * maximum fitting error
results%fittingErrorMax = maxval(results%fittingError)
! * average fitting error
call average(results%fittingError, results%fittingErrorAv)
! * RMS fitting error
results%fittingErrorRMS = sqrt( sum( (results%fittingError - results%fittingErrorAv)**2 ) / real(nStr,dp) )
deallocate(PI)

contains


!....................................................................................................
! This is only an OUTPUT routine, you have to call set_GCdeltas afterwards
subroutine output__switch_GCdeltas(curr,next)
integer, intent(in) :: curr, next

print *
!switching garbulsky-ceder for all observables (including energy)
  if(iObs == 1) then
    write(*,'(A,I2,A)') "Switching deltas for observable ", iObs, ", i.e. energy"
  else
    write(*,'(A,I2)') "Switching deltas for observable ", iObs
  endif

  write(*,'(A,F12.5,A,F12.5)') "delta1:",CE%delta(1,curr,iObs)," --> ",CE%delta(1,next,iObs)
  write(*,'(A,F12.5,A,F12.5)') "delta2:",CE%delta(2,curr,iObs)," --> ",CE%delta(2,next,iObs)
  write(*,'(A,F12.5,A,F12.5)') "delta3:",CE%delta(3,curr,iObs)," --> ",CE%delta(3,next,iObs)
  print *
  write(*,'(A)') "Update CVS..."
  print *
end subroutine output__switch_GCdeltas
!....................................................................................................
subroutine output__GCdeltas(curr)
integer, intent(in) :: curr

print *
!switching garbulsky-ceder for all observables (including energy)
  if(iObs == 1) then
    write(*,'(A,I2,A)') "deltas for observable ", iObs, ", i.e. energy"
  else
    write(*,'(A,I2)') "deltas for observable ", iObs
  endif

  write(*,'(A,F12.5,A,F12.5)') "delta1:",CE%delta(1,curr,iObs)
  write(*,'(A,F12.5,A,F12.5)') "delta2:",CE%delta(2,curr,iObs)
  write(*,'(A,F12.5,A,F12.5)') "delta3:",CE%delta(3,curr,iObs)
  print *
end subroutine output__GCdeltas
!....................................................................................................
subroutine output__population_info(iP)
integer iP
print *
write(*,'(A,1x,I2)') "Observable #", iObs
write(*,'(A,1x,I2)') "Population #", iP
write(*,'(A)') "--------------------------------------------------------------------------------"
write(*,'(A)',advance='no') "Number of activated genes: "
do iIndv=1,nIndv
  write(*,'(I4,3x)',advance='no') individuals(iIndv,iP)%nOnGenes
enddo
print *
print *
end subroutine output__population_info
!....................................................................................................
subroutine output__generation_info_header

write(*,'(//)')
write(*,'(A11," | ",A23," | ",A21," | ",A17," | ",A8," | ",A3)') &
     "generation","population","best alive","best dead","ETA","div"
write(*,'(A11," | ",A23," | ",A21," | ",A17," | ",8x," | ")') &
     "","#, <CVS>, diversity","CVS (#genes,age)","CVS (#genes)"
write(*,'(98("-"))')

end subroutine output__generation_info_header
!....................................................................................................
subroutine output__generation_info(iG,iP) ! output the information line after each generation
integer iG,iP
integer bestIdxArray(1), bestIdx
real(dp),save :: ETA
character(3),save :: ETAunit
character(8),save :: ETAstring
integer, parameter :: ETAsteps = 10

! Determine ETA:
if (iG==1 .or. any(CE%gen4DeltaSwitch == iG)) then
  ETAstring = "     N/A"
  call timing(walltime_curr,cputime_curr)
else
  if (mod(iG,ETAsteps)==0 .and. iP==1) then
     walltime_prev = walltime_curr   ! save the old values
     cputime_prev  = cputime_curr
     call timing(walltime_curr,cputime_curr)  ! get the new values
     ETA = (walltime_curr - walltime_prev) * (nGeneration-iG) / (1.*ETAsteps) / 60.
     ETAunit = "min"
     if (ETA>60) then; ETA=ETA/60.d0; ETAunit="h";
       if (ETA>24)  then; ETA=ETA/24.d0; ETAunit="d"; endif;
     endif
     write(ETAstring,'(F4.1,1x,A3)') ETA, ETAunit
  endif
endif
! End Determine ETA


bestIdxArray = minloc(individuals(:,iP)%cvs)
bestIdx      = bestIdxArray(1)


!print *, "avg activated=", sum(individuals(:)%nOnGenes)*1.0/size(individuals)

!if (iP==1) then
  write(*,'(I5,A1,I5," | ")',advance='no') iG,"/",nGeneration
!else
!  write(*,'(A5,A1,A5," | ")',advance='no') " ", " ", " "
!endif

write(*,'(I5,1x,F10.5,1x,F6.1," | ",F11.6," (",I3,1x,I3,")"," | ")', advance='no') &
       iP, &
       sum(individuals(:,iP)%cvs)/nindv, avgDeviation(iP), &
       individuals(bestIdx,iP)%cvs, individuals(bestIdx,iP)%nonGenes, individuals(bestIdx,iP)%age

if (bestDeadIndividual(1,iP)%cvs/=-1) then ! this corresponds to: there is a dead individual
  write(*,'(F11.6," (",I3,")"," | ",A8, " | ", A2)') &
       bestDeadIndividual(1,iP)%cvs, bestDeadIndividual(1,iP)%nonGenes, &
       ETAstring, &
       adjustl(diversityCountermeasure(iP))
else
     write(*,'(A11," (",A3,")"," | ",A8," | ",A2)') &
          "N/A","N/A", &
          ETAstring, &
          adjustl(diversityCountermeasure(iP))
endif

!DEBUG>! DEBUG: write out the distribution of pairs/triplets/...
!DEBUG>write(666,'(I7,4x,5(I3,2x))') &
!DEBUG>iGeneration, &
!DEBUG>count( clusters( individuals(bestIdx)%onGenes( 1:individuals(bestIdx)%nonGenes ) )%nV == 2), &
!DEBUG>count( clusters( individuals(bestIdx)%onGenes( 1:individuals(bestIdx)%nonGenes ) )%nV == 3), &
!DEBUG>count( clusters( individuals(bestIdx)%onGenes( 1:individuals(bestIdx)%nonGenes ) )%nV == 4), &
!DEBUG>count( clusters( individuals(bestIdx)%onGenes( 1:individuals(bestIdx)%nonGenes ) )%nV == 5), &
!DEBUG>count( clusters( individuals(bestIdx)%onGenes( 1:individuals(bestIdx)%nonGenes ) )%nV == 6)

end subroutine output__generation_info

end subroutine GeneticAlgorithm_for_clusterExpansion
!****************************************************************************************************





!****************************************************************************************************
! Purpose:
! The population diversity with respect to the best (living) individual is checked. If the individuals
! are (on average) less than a user-specified number of genes apart from each other, the population is setup anew (while the
! best individual is kept, though)
subroutine check_population_diversity( individuals, GA, avgDeviation, resetup )
type(t_individual), intent(in)    :: individuals(:)    ! the individuals
type(GA_variables), intent(inout) :: GA                ! some GA variables
real(dp), intent(out)             :: avgDeviation      ! average deviation of the individuals from the best individual
logical, intent(out)              :: resetup           ! do we have to setup the population anew?

integer            :: bestIndv        ! index of the best individual

integer iIndv, jIndv, nIndv       ! counter

nIndv        = size(individuals)
resetup      = .false.
avgDeviation = 0._dp

! ! calculate the average deviation of the individuals from each other
! do iIndv=1,nIndv
!   do jIndv=iIndv+1,nIndv
!   ! add up how many genes are different
!     avgDeviation = avgDeviation + count( individuals(iIndv)%genome .neqv. individuals(jIndv)%genome )
!   enddo
! enddo
! avgDeviation = avgDeviation / ( nIndv * (nIndv-1) ) * 2._dp


! calculate the average deviation of the individuals from the best individual
call get_index_of_best_individual( individuals , bestIndv ) 
do iIndv=1,nIndv
    ! add up how many genes are different
    avgDeviation = avgDeviation + count( individuals(iIndv)%genome .neqv. individuals(bestIndv)%genome )
enddo
avgDeviation = avgDeviation / (nIndv-1) ! get the average value. (nIndv-1 and not nIndv because the best individual
                                        ! itself should not be taken into account

! if the average deviation is too small then the best individual is stored, the
! population as a whole freshly populated, and the best individual restored.
if (avgDeviation<GA%diversityLimit) then
  resetup=.true.
endif

end subroutine check_population_diversity
!****************************************************************************************************






!****************************************************************************************************
! Purpose:
! Take the average of "quantities" according to "mask" and put it into "avg"
subroutine average( quantities, avg, mask )
real(dp), intent(in)  :: quantities(:)   
real(dp), intent(out) :: avg
logical, intent(in), optional :: mask(:)

integer nQ

if (present(mask)) then
  if (size(mask)/=size(quantities)) stop "ERROR: wrong mask setup for average of array rank 1"
  nQ      = count(mask.eqv..true.)
  avg     = sum( quantities(:), mask=mask ) / real(nQ,dp)
else
  nQ      = size(quantities,1)
  avg     = sum( quantities(:) ) / real(nQ,dp)
endif

end subroutine average
!****************************************************************************************************


!****************************************************************************************************
! Purpose:
! Take the J's for each subset and average over the subsets. Uses the subroutine "average" above.
subroutine average_J(J, Jav, subsetMask)
real(dp), intent(in)      :: J(0:,:)   ! { cluster#, subset# }
real(dp), intent(out)     :: Jav(0:)   ! { cluster# }
logical, intent(in)       :: subsetMask(:)    ! { subset# }

integer iSS, nSS
integer iCl, nCl

nCl = size(J,1)    ! number of clusters
nSS = size(J,2)    ! number of subsets
if (nSS /= size(subsetMask)) then
  print *, "ERROR in average_J:"
  print *, "- size subsetMask: ", size(subsetMask)
  print *, "- subsetMask: ", subsetMask
  print *, "- J: number of clusters: ", size(J,1)
  print *, "- J: number of subsets:  ", size(J,2)
  stop "ERROR: bad programming: wrong mask setup for average_J"
endif

do iCl=0,nCl-1
  call average( J(iCl,:), Jav(iCl), mask=subsetMask)
enddo


end subroutine average_J
!****************************************************************************************************



!****************************************************************************************************
! Purpose:
! Out of a list of individuals ("ind"), return the index "idx" of the individual that has the lowest
! CVS. Obviously, this is a two-liner. It was put into a separate subroutine, because FORTRAN's minloc
! function always returns an array, and we only want to have a number, lest it get's too messed up in
! the main code.
pure subroutine get_index_of_best_individual( ind , idx )
type(t_individual), intent(in)    :: ind(:)     ! { # }
integer           , intent(out)   :: idx

integer iObs, nObs, bestIdx(1)

bestIdx    = minloc(ind(:)%CVS)  
idx        = bestIdx(1)

end subroutine get_index_of_best_individual
!****************************************************************************************************





!****************************************************************************************************
! Purpose:
! The new interface of the quadprg fitting routine.
! Provide the values that are fitted ("observable"), the weight for each value ("weight"), the correlation
! functions ("PI") and the constraint matrix and vector, in order to get the "J" as output.
!
subroutine make_ceFit_with_constraints( observable, weight, PI, clustersDamping, constraint_matrix, constraint_vector, J, CE, success )
use quadprog                      !! Quadratic programming
real(dp),           intent(in)    :: observable(:)      ! { structure# }, will return the value of observable as a function of the structure index 
integer,            intent(in)    :: weight(:)          ! { structure# }: weight of a structure
real(dp),           intent(in)    :: PI(:,0:)           ! { structure#, cluster# }
real(dp),           intent(out)   :: J(0:)              ! { cluster# }
real(dp),           intent(in)    :: clustersDamping(0:)! { cluster# }: for each cluster the corresponding damping (=c*r^lambda)
real(dp),           intent(in)    :: constraint_matrix(:,0:), constraint_vector(:)   ! constraints
type(CE_variables), intent(inout) :: CE
logical,            intent(out)   :: success

integer nConstraints
integer rvalue

integer  :: iStructures, nStructures, iClusters, jClusters, nClusters
real(dp) :: linear(0:size(PI,2)-1), quadratic(0:size(PI,2)-1, 0:size(PI,2)-1)

real(dp) :: weighted_observable(size(observable)), det
real(dp) :: weighted_PI(size(PI,1),0:size(PI,2)-1)

!real(dp) :: avgAbsJ, sigmaAbsJ

#ifdef __MKL__
! MKL call
integer ierr
real(dp) :: quadratic_tmp(0:size(PI,2)-1, 0:size(PI,2)-1)
integer  :: jpvt(size(PI,2))
!real(dp) :: work(16*(3*size(PI,2)+1))
real(dp) :: work(3*size(PI,2))
real(dp) :: tau(size(PI,2))
#endif


! According to the PhD Thesis of Ole, he used the equation   
!     ct x + 1/2 xt Q x =!= min (*)
! where Q ~ quadratic, ct ~ linear
!
! Our minimisation of chi^2 includes a -2 in front of the linear J term
! and +1 in front of the quadratic term. 
! For "quadratic" Ole didn't use a factor 1/2; this basically means he multiplied
! equation (*) by 2 which doesn't harm the minimisation. However, then he would not need
! a - sign for the linear term. He did not comment the function that gets minimised in the 
! quadprg routine... 
!
! By some kind of trial and error and by comparing to Ole's setup of the fitting, we
! conclude that we determine J from
!
! -2 ( observable^T * PI ) * J + J^T * ( PI^T * PI ) * J =!= min
! <---------------------->             <----------->
!      linear factor                  quadratic factor
!
! Beware that both "observable" and "PI" are scaled with sqrt(weight) for each structure (see below).


!write(*,'(A$)') "entering make_ceFit_with_constraints"

nStructures = size(PI,1)   ! number of input structres
nClusters   = size(PI,2)   ! number of fitting variables: clusters
nConstraints= size(constraint_vector) / nStructures  / 2


call seteps(CE, 1.0E-11_dp)         ! sets CE%qp_eps
call setmaxcounts(CE, 10000)       ! sets CE%maxcounts
call setinfinity(CE, 10.0E+10_dp)  ! sets CE%inf


! TK's setup: ................................................................................
!
! (A) incorporate the weights of the structures
! (B) setup linear term
!     - without damping (! damping in the linear term is wrong !)
! (C) setup quadratic term
!     - without any damping
!     - add clusters damping (! damping belongs to the quadratic term !)
!     - add static regularization (a hard-coded constan) for better numerical stability 
!       (add least, that's what tk guesses; Ole had it in his code but never really commented
!       about it (surprise, surprise...))
!
!
! (A)....................
! For the minimisation of chi^2 we can also provide "weights": those (de)emphasise certain input
! values by
!      chi_i^2 = w_i (E_i(DFT) - E_i(CE))^2
! where the index i denotes a certain input observable (e.g. the energy of structure i). w_i is
! the respective weight. In the following, we adsorb w_i into the "( )" by:
forall(iStructures=1:nStructures)
  weighted_observable(iStructures) = observable(iStructures) * sqrt(real(weight(iStructures),dp))
  weighted_PI(iStructures,0:)      = PI(iStructures,0:)      * sqrt(real(weight(iStructures),dp))
end forall

! (B)....................
linear      = - 2._dp * matmul( weighted_observable, weighted_PI )   ! linear term without damping

! (C)....................
! Note: According to the Goldfarb-Idnani Algorithm for solving the quadratic programming
! problem, the matrix representing the quadratic factor must be positive definite. We will
! partially take care of that below. (Maybe this test is also included in quadprog.f90.)
quadratic   = matmul( transpose(weighted_PI) , weighted_PI )         ! quadratic term


#ifdef __MKL__

! check for linear independence of the quadratic term
quadratic_tmp = quadratic
call DGEQPF(nClusters,nClusters,quadratic_tmp(0:,0:),nClusters,jpvt,tau,work,ierr) 
! DGEQPF computes a QR factorization with column pivoting of a real M-by-N matrix 
!       A: A*P = Q*R
!
! MIND: if you want to determine, say, the determinant of A, then
!           det(A) = det(Q)*det(R) = det(R)
!       in the usual QR decomposition. BUT the routine DGEQPF uses some pivoting (P),
!       so some columns may have switched order and therefore the determinant might
!       have changed the sign:
!           det(A) = det(R) * sign(P)
!
                                                                                   
    ! Note: according to an INTEL help page (which might not be applicable to all MKL implementations), the matrix
    !       A is overwritten in the following way:
    !       - m>=n: the elements below the diagonal are overwritten by the details of the unitary (orthogonal) matrix Q, 
    !               and the upper triangle is overwritten by the corresponding elements of the upper triangular matrix R.
    !       - m< n: the strictly lower triangular part is overwritten by the details of the matrix Q, 
    !               and the remaining elements are overwritten by the corresponding elements of the m-by-n upper trapezoidal matrix R.

if(ierr/=0) stop "ERROR: bad programming: problems in call to DGEQPF"
det = permutation_parity(jpvt) * product((/(quadratic_tmp(iClusters,iClusters),iClusters=0,nClusters-1)/))

if ( det < CE%qp_eps ) then  ! if abs(det)<eps then the quadratic term is not linearly independent.
                             ! if det<0 at all it cannot be positive definite.
!   if (det < 0) then
!      write(*,'(///,"Potential failsafe triggered!")')
!      print *,"It should be the case (based on my tests [GH]) that the PI^T*PI matrix must be positive definite, by construction"
!      print *,"If this message triggered, then perhaps something is wrong"
!      stop
!   endif
!if ( equal(det,0._dp,CE%qp_eps) ) then    ! (obsolete now)
   success = .false.
   J = 0._dp
!   write(*,'(A)',advance='no') "!"
   return
endif

#else

! Adding some weight to the diagonal elements. Ole had this in his code.
! It makes the quadratic matrix linearly independent by force: adding a small fraction
! of an identity matrix. (What happens to "positive definite" check I don't know. Mabye
! the algorithm in quadprg fails then.)
forall(iClusters=0:nClusters-1)
  quadratic(iClusters,iClusters) = quadratic(iClusters,iClusters) + 1*CE%qp_eps
end forall

#endif

! Add the damping to the quadratic term
! (this should probably not happen before the linear-independency check above, because the
! damping makes the quadratic matrix linear independent most of the time)
forall(iClusters=0:nClusters-1)
  quadratic(iClusters,iClusters) = quadratic(iClusters,iClusters) + clustersDamping(iClusters)  
end forall


!********************************************************************************
!********************************************************************************
! FITTING

J = 0._dp

!print *, "using for fit:"
!print *, "nclusters=",nclusters
!print *, "nstructures=",nstructures
!write(990,*) "constraint_matrix=",constraint_matrix
!print *, "constraint_vector=",constraint_vector
call quadprg(nClusters,   &
             nStructures * nConstraints * 2, &
             quadratic,   & ! quadratic part of eq.
             linear,      & ! linear part of eq.
             constraint_matrix, & ! origin: call setup_quadprog_constraints
             constraint_vector, & ! dito
             J,     & ! <-- this is actually what we want to determine here!
             rvalue, & 
             CE) 
!print *, "rvalue=",rvalue
!stop "after quadprg"
! FITTING
!********************************************************************************
!********************************************************************************

!the following did not really work out... :-/
! ! J concistency check
! nClusters = size( J(11:) )
! avgAbsJ   = sum( abs( J(11:) ) ) / nClusters
! sigmaAbsJ = sqrt( dot_product( (abs(J(11:))-avgAbsJ),(abs(J(11:))-avgAbsJ) ) ) / (nClusters-1)
! if (any( abs(J(11:)) > avgAbsJ + 2*sigmaAbsJ )) then
! !  write(*,*) J(0:10)
! !  write(*,*) J(11:)
! !  write(*,*) "avg|J|",avgAbsJ
! !  write(*,*) "sigma",sigmaAbsJ
! 
! !!!!
! !!!!
! !!!! DO NOT YET SUBMIT
! !!!! IF THIS WORKS YOU HAVE TO HAVE THE ONSITES AS A VARIABLE
! !!!!
! !!!!
! 
!   rvalue = -2
!   !write(*,'(A)',advance='no') "!"
! endif



! Some notes:
! * rvalue = -1 : the quadratic matrix was not positive definite
! * if the fit did not succeed (rvalue /= 0)
!   - either: all J == 0, because we set them to 0 before the call.
!   - or    : J are extremely large
!   Either way, we set them to 0:
if (rvalue==0) then; success=.true. ;
else;                success=.false.; J = 0._dp;
endif


!write(*,*) "done"

end subroutine make_ceFit_with_constraints
!****************************************************************************************************



!****************************************************************************************************
! Purpose:
! Given the CE J's and the PI's (that belong to the observable, e.g. the energy of a structure), 
! return the errors to the real observables (e.g. determined by DFT)
subroutine get_errors_of_CE( observable, J, PI, errors )
real(dp), intent(in)    :: observable(:)   ! 
real(dp), intent(in)    :: J(0:), PI(:,0:) ! { clusters# }, { structures#, clusters# }
real(dp), intent(out)   :: errors(:)       ! 

integer iO, nO

nO=size(observable)

!print*, "obs: ", observable(7)
!print*, " J ", J(4)
!print*, " PIT ", PI(7,4)

errors = 0._dp
forall(iO=1:nO)
  errors(iO) = - (observable(iO) - dot_product( J(0:) , PI(iO,0:) ) )   ! tk: makes IMHO more sense to have the "-"
                                                                        ! sign. A negative error thus means, that the
                                                                        ! CE's observable is *lower* than expected, and
                                                                        ! vice versa.
end forall

end subroutine get_errors_of_CE
!****************************************************************************************************




!****************************************************************************************************
! gh+tk, adapted by tk+sm
!
! Purpose:
! Setup "nIndv" individuals "indv" and their genome with "nGenes" genes: activate a maximum of
! "maxFitpar" genes. If the optional "activateFirstGenes" is true, then do not choose random
! genes, but activate the first "maxFitpar" ones.
!
SUBROUTINE setup_population(indv, nIndv, maxFitpar, nGenes, nGenesAlwaysOn, GA, warn, activateFirstGenes)
type(t_individual), intent(out)   :: indv(:)    ! intent out. Initialize and passout
integer, intent(in)               :: nIndv      ! number of individuals ("genomes")
integer, intent(in)               :: maxFitpar  ! maximal number of fitparameters turned on
integer, intent(in)               :: nGenes     ! total number of fitparameters ("genes"): including CONSTANT and ONSITES!!
integer, intent(in)               :: nGenesAlwaysOn  ! number of genes that are supposed to be always switched on (i.e.
                                                     ! constant and onsites). This number is always connected to
                                                     ! the FIRST genes in the genome, i.e., the first nGenesAlwysOn genes
                                                     ! are switched on.
type(GA_variables),intent(inout)  :: GA         ! just need this for the random seed
logical, intent(in)               :: warn       ! give a warning message if not all genes were used during the setup?
logical, intent(in), optional     :: activateFirstGenes

integer             :: iIndv, ifp, gene, listL, iGenes
integer             :: list(nGenes)

integer :: nSwitchOn
integer :: iRefresh

!write(*,'(A$)') "entering setup_population"


do iIndv = 1, nIndv ! Allocate the genomes and zero them out
   allocate(indv(iIndv)%genome(0:nGenes-1))
   allocate(indv(iIndv)%OnGenes(1:nGenes))
   indv(iIndv)%genome                     = .false.
   indv(iIndv)%genome(0:nGenesAlwaysOn-1) = .true.  ! always switch those on
   indv(iIndv)%age                        =  1
   indv(iIndv)%cvs                        = -1      ! initial failsafe value
enddo

if (present(activateFirstGenes) .and. activateFirstGenes) then

   ! For each individual, activate the first maxFitpar genes
  do iIndv = 1, nIndv
    indv(iIndv)%genome(0:maxFitpar-1) = .true.
  enddo

else

   iRefresh = 0
   list = (/ (iGenes,iGenes=0,nGenes-1) /)  ! starts at 0 !!!
   listL = nGenes ! Keep track of the length of the unused list
   
   ! For each individual, select random genes from the pool
   do iIndv = 1, nIndv 
     nSwitchOn = maxFitpar-nGenesAlwaysOn  ! maximal number to switch on 
     nSwitchOn = nSwitchOn * ( ( 1 - 0.0 )*rnd(GA) + 0.0 )  ! switch on between 0% and 100% of the maximal number. Originally, GH had intended to
                                                            ! always switch on the maximum number of genes. However, the current form should make
                                                            ! the population a little bit more diverse. We (GH+TK) don't think that we want to have
                                                            ! an additional user input for that; shouldn't matter much.
                                                            ! However, the most important thing---as GH pointed out---is that each genes gets activated
                                                            ! roughly the same number of times. That's why we use the "take_random_element_out_of_list"
                                                            ! function. We will also print an Warning message if some of the genes were never activated
                                                            ! at all (and if the flag warn==.true.)
     do ifp = 1, nSwitchOn
       call take_random_element_out_of_list(list,listL,gene,GA)
       indv(iIndv)%genome(gene) = .true.
       if(listL<1)then ! list of genes is exhausted. Refresh it and keep going
         iRefresh = iRefresh + 1
         list = (/ (iGenes,iGenes=0,nGenes-1) /)
         listL = nGenes
       endif
     enddo
   enddo

   if (iRefresh == 0 .and. warn) then ! maxFitpar == 0 for the setup of the "dead" individual...
     print *
     write(*,'(A)') "WARNING: The number of individuals was too low to activate each gene at least once."
     print *
   endif

endif

do iIndv=1,nIndv  ! don't use FORALL here!
  call set_genome_information(indv(iIndv)%genome,indv(iIndv)%nOnGenes,indv(iIndv)%onGenes,condition=.true.)
end do

!write(*,*) "done"

END SUBROUTINE setup_population
!****************************************************************************************************



!****************************************************************************************************
! This is a helper routine to mimic the behaviour of the old UNCLE GA concerning the pair selection
! during the GA. It is assumed that the individuals indv(:) are already set up (maybe those are children!)
! an the first nGenesAlwaysOn are always activated (e.g., constant and onsite). nPairs is the total 
! number of pairs in the cluster pool, i.e., the total number of pair-genes in the individuals' genome.
!
! The routine counts how many (random) pairs are activated, then deactivates those and activates the
! same number of pairs again, starting from the smallest pair in a consecutive sequence.
!
subroutine correct_pairs(executeRoutine, nPairs, indv, nGenesAlwaysOn)
logical           , intent(in)    :: executeRoutine  ! will we execute this routine at all?
integer           , intent(in)    :: nPairs
type(t_individual), intent(inout) :: indv(:)
integer           , intent(in)    :: nGenesAlwaysOn

integer iIndv, nIndv
integer nPairsActivated
integer firstPairIdx

if (.not. executeRoutine) return ! i.e., we don't want to correct the pair selection


nIndv = size(indv)

firstPairIdx = nGenesAlwaysOn-1 +1   ! that's the index (in the genome) of the first pair
                                     ! offset of -1 because the genome starts at 0
                                     ! additional offset of +1 because the first after the always-on-genes
                                     ! is the first pair. (yes, I know we could simply use nGenesAlwaysOn, but
                                     ! I think it's more clear in that way)

do iIndv=1,nIndv
  ! count how many pairs are activated
  nPairsActivated = count(  indv(iIndv)%genome( firstPairIdx : firstPairIdx + nPairs-1 ) .eqv. .true. )

  ! set all pairs to F
  indv(iIndv)%genome( firstPairIdx : firstPairIdx + nPairs-1 ) = .false.
  
  ! activate the same number of pairs again, but the smallest pairs only
  indv(iIndv)%genome( firstPairIdx : firstPairIdx + nPairsActivated-1 ) = .true.

  ! set info
  call set_genome_information(indv(iIndv)%genome,indv(iIndv)%nOnGenes,indv(iIndv)%onGenes,condition=.true.)
enddo

end subroutine correct_pairs


!***************************************************************************************************
! This routine sets up all the individuals anew and keeps only the best individual for the new
! population. Though, this best individual is aged to its maximum age, so it's probably going to be
! killed in the next generation.
subroutine resetup_population( individuals, GA )
type(t_individual), intent(inout)    :: individuals(:)    ! the individuals
type(GA_variables), intent(inout)    :: GA

integer            :: bestIndv        ! index of the best individual
type(t_individual) :: bestIndividual  ! temporary storage for best individual

integer iIndv, nIndv       ! counter

nIndv        = size(individuals)

call get_index_of_best_individual( individuals , bestIndv ) 
bestIndividual = individuals(bestIndv)
do iIndv=1,nIndv
  if (iIndv==bestIndv) cycle
  deallocate(individuals(iIndv)%onGenes)
  deallocate(individuals(iIndv)%genome)
enddo
call setup_population( individuals(:), nIndv, GA%nMaxClusters, GA%nTotClusters, GA%nAlwaysOn, GA, warn=.false.)
individuals(1) = bestIndividual
individuals(1)%age = GA%maximumAgeForIndividual   ! age the best individual

end subroutine resetup_population
!***************************************************************************************************




!***************************************************************************************************
!
! TAKEN OUT OF THE DEVEL-NEWGAFIT BRANCH of UNCLE
! gh, tk
!
! This routine takes in a list and an index (length of the remaining list) and returns one element  
! of the list and then "shrinks" the list. The routine assumes that the list has no negative values 
! as part of the original data in the list.                            
SUBROUTINE take_random_element_out_of_list(list,length,item,GA)
integer, intent(inout) :: list(:)
integer, intent(inout) :: length
integer, intent(out)  :: item ! the item that was select and removed from the list                  
type(GA_variables) GA ! need this for the seed only                                                 

real(dp) r
integer indx

if (any(list(1:length)==-1)) stop "Bad list in take_element_out_of_list routine in ga_fitting.f90"
if (length==0) stop "ERROR: Cannot choose an element of the list since its length is 0"
r = rnd(GA)
indx = ceiling(r*length) ! Pick the index of item at random from the remaining list                 
item = list(indx) ! Get the value of the list at the selected index                                 
list = (/list(1:indx-1),list(indx+1:size(list)),-1/) ! Update the list of items                     
length = length - 1 ! shrink the list                                                               
if(length<0) stop "Blew it in the take_element_out_of_list routine in ga_fitting.f90"
END SUBROUTINE take_random_element_out_of_list
!***************************************************************************************************






!***************************************************************************************************
!
! TAKEN OUT OF THE DEVEL-NEWGAFIT BRANCH of UNCLE
! gh, tk
!
! CHANGED a little bit:
! Replaced "fitSet" and "validSet" by a mere "subset".
! The substructure of those sets was "%str(:)".
! The corresponding substructure of "subset" is "%fitList" and "%predList", resp.
!
!
! This routine generates a list of fitting (training) sets to be used in a K-fold cross-validation
! approach to making a fit. See "cross validation" at Wikipedia. Actually what we do is make a list
! of validation sets, take those away from the entire data set, to produce the fitting sets. We use
! K-fold cross validation, repeated (with different sets) nRep times. 
SUBROUTINE make_fitting_sets(nStr,k,nRep,subset,GA)
integer, intent(in) :: nStr, k ! number of structures, number of pieces to divide the data into (i.e. number of subsets)
integer, intent(in) :: nRep ! number of times to repeat the k-fold division
type(t_subset), intent(out) :: subset(:)
type(GA_variables) :: GA ! Just need this for the random number seed

integer iRep, offset, ik
integer, allocatable :: validSetMask(:), fitSetIndx(:)
integer, pointer :: indx(:)

integer i

! failsafe
if (k > nStr) stop "ERROR: number of subsets > number of structures that are used for subsets"

allocate(indx(nStr),validSetMask(nStr),fitSetIndx(nStr))
fitSetIndx = (/ (i,i=1,nStr) /)
validSetMask = ceiling(fitSetIndx / (real(nStr,dp)/k) )
do iRep = 1, nRep
   call randomize_list(fitSetIndx,GA)
   offset = (iRep -1)*k
   do ik = 1, k
      allocate(subset(offset+ik)%fitList(count(validSetMask/=ik)),&
               subset(offset+ik)%predList(count(validSetMask==ik)))
      subset(offset+ik)%fitList  = pack(fitSetIndx,validSetMask/=ik)
      subset(offset+ik)%predList = pack(fitSetIndx,validSetMask==ik)
      call heapsort(indx,subset(offset+ik)%fitList)
      call heapsort(indx,subset(offset+ik)%predList)
   enddo
enddo

END SUBROUTINE make_fitting_sets



!***************************************************************************************************
!
! TAKEN OUT OF THE DEVEL-NEWGAFIT BRANCH of UNCLE
! gh, tk
!
! This routine takes a list and randomizes it.
SUBROUTINE randomize_list(list,GA)
integer, intent(inout) :: list(:)
type(GA_variables) GA ! need this for the seed only

real(dp) r
integer indx, iL, nL, item

nL = size(list)
do iL = 1, nL
   r = rnd(GA)
   indx = ceiling(r*(nL - iL + 1)) ! Pick the index of item at random from the remaining list
   item = list(indx) ! Get the value of the list at the selected index
   list = (/list(1:indx-1),list(indx+1:size(list)),item/) ! Update the list of items
enddo
END SUBROUTINE randomize_list
!****************************************************************************************************



!****************************************************************************************************
! Take current individuals and replace the worst ones by their children (don't have to be better, will
! always replace the worst parents)
subroutine make_new_generation( individuals, children, GA, bestDeadIndividual )
type(t_individual), intent(inout) :: individuals(:)         ! current individuals (IN), individuals of new generation (OUT)
type(t_individual), intent(inout) :: children(:)            ! current children. Are going to replace some of the individuals (IN)
type(GA_variables), intent(inout) :: GA
type(t_individual), intent(inout) :: bestDeadIndividual     ! storage for the best individual that has died.

integer nIndv, iIndv, nChild, iChild

integer, pointer  :: indexIndividuals(:), indexChildren(:)
real(dp), pointer :: individualsCVS(:), childrenCVS(:)

real :: NaN = 0.

integer deadIndex

nIndv = size(individuals)   ! number of individuals
nChild= size(children)      ! number of children


allocate(individualsCVS  (1:nIndv))  ! temporary container of the CVS of the individuals
allocate(childrenCVS     (1:nChild)) ! temporary container of the CVS of the children
individualsCVS = individuals(:)%CVS
childrenCVS    = children(:)%CVS


! see if we have to kill some individuals due to their age
deadIndex = -1   ! failsafe starting value
do iIndv=1,nIndv
  if (individuals(iIndv)%age >= GA%maximumAgeForIndividual) then ! is the individual too old?
    ! yes, it is too old. Now compare it to the best individual that has died before and see, if it is better
    if (.not.isnan(individualsCVS(iIndv))) then  ! is my CVS defined at all?
    if (individualsCVS(iIndv) < bestDeadIndividual%cvs .or. bestDeadIndividual%cvs<0 .or. isnan(bestDeadIndividual%cvs)) then
        !<------------------------------------------->!     !<-------------------->!      !<------------------------->!
        ! is my CVS better than the dead one?               ! is there a dead one?        ! is the dead one defined at all?
                                                         ! If no individual has ever died,
                                                         ! this is -1
      bestDeadIndividual = individuals(iIndv) ! if this is the best individuum ever: keep a record
      deadIndex          = iIndv              ! if an individual will die and it's good: keep its index
    endif
    endif

    individualsCVS(iIndv) = maxval(individualsCVS)+1 ! kill it. Actually, this not really "kills" it, but sets it extremely unfavourable, so it gets replaced
                                                     ! by children afterwards. (Note that by this construction the true CVS hierarchy gets mixed up. Also note
                                                     ! that if there are fewer children than "kills" some of the kills will not really be replaced and survive,
                                                     ! with a really bad CVS)
  endif
enddo

if (deadIndex /= -1) then
individualsCVS(deadIndex) = maxval(individualsCVS)+1 ! make sure that at least the best old individual gets REALLY killed,
                                                     ! even if there are many other individuals that we want to kill too.
                                                     ! (Note that by the construction above this might NOT happen if there are
                                                     ! fewer children than "kills". Here, we set the best "dead" one so 
                                                     ! so unfavourable that it gets killed even if there is only 1 child available.
                                                     ! Otherwise we run into memory/pointer problems because the code assumes
                                                     ! that the best dead one is no longer stored in the "living indivuals" list.
endif

! heapsort setup
allocate(indexIndividuals(1:nIndv)) 
allocate(indexChildren   (1:nChild)) 
indexIndividuals= (/ (iIndv , iIndv =1,nIndv ) /) 
indexChildren   = (/ (iChild, iChild=1,nChild) /)

! heapsort
call heapsort(indexIndividuals, individualsCVS)  ! sorts smallest first
call heapsort(indexChildren   , childrenCVS   )

! replace individuals by children:
do iIndv= nIndv - GA%nReplacements + 1, nIndv ! replace worst parents with best children
                                              ! IMPORTANT: see comment below the call to make_new_generation !!!! 
                                              ! You have to detach children by nullifying their contents afterward !!!!
  if ( indexIndividuals(iIndv) /= deadIndex ) then  ! if this is NOT the new best dead individual
                                                    ! (if it is, don't deallocate, because those arrays are now in
                                                    ! use in bestDeadIndividual)
    deallocate(individuals( indexIndividuals(iIndv) )%onGenes             ) ! deallocate OLD arrays
    deallocate(individuals( indexIndividuals(iIndv) )%genome              ) ! deallocate OLD arrays
    deallocate(individuals( indexIndividuals(iIndv) )%J                   ) ! deallocate OLD arrays
    deallocate(individuals( indexIndividuals(iIndv) )%subsetFitSuccessful ) ! deallocate OLD arrays
  endif
  individuals( indexIndividuals(iIndv) ) = children( indexChildren( -nIndv + GA%nReplacements + iIndv) ) ! NEW arrays are assigned
                                                                                                         ! directly from children
end do

! deallocate the other children (those that didn't find their way into the new generation)
do iChild = GA%nReplacements+1, nChild
  deallocate(children( indexChildren( iChild ) )%onGenes             )
  deallocate(children( indexChildren( iChild ) )%genome              )
  deallocate(children( indexChildren( iChild ) )%J                   )
  deallocate(children( indexChildren( iChild ) )%subsetFitSuccessful )
enddo

! failsafe
! If the contents of children are correctly detached (see comment above), then the following failsafe loop that
! updates the individuals' information according to their genome should not be necessary. It doesn't cost that
! much, though, so tk's leaving it in place for now.
do iIndv=1,nIndv
  individuals(iIndv)%age = individuals(iIndv)%age + 1
  call set_genome_information( individuals(iIndv)%genome, individuals(iIndv)%nOnGenes, individuals(iIndv)%onGenes, condition=.true.)
enddo

! failsafe
! If the contents of a dying individual are copied correctly, we shouldn't need the following. But, again, doesn't matter
! to be sure.
call set_genome_information( bestDeadIndividual%genome, bestDeadIndividual%nOnGenes, bestDeadIndividual%onGenes, condition=.true.)

deallocate(indexIndividuals,individualsCVS)
deallocate(indexChildren,childrenCVS)

end subroutine make_new_generation
!****************************************************************************************************




!****************************************************************************************************
! Purpose:
! Take "individuals" and mutate their genomes at random positions, according to the probability
! given
!
subroutine mutate_individuals( individuals, mutationProbability, GA ) 
type(t_individual), intent(inout) :: individuals(:)
real(dp),           intent(in)    :: mutationProbability
type(GA_variables), intent(inout) :: GA ! used for general GA variables and for the random numbers in particular

integer iGene, nTrueGenes, nOnTrueGenes
integer, allocatable :: onTrueGenes(:)
integer iMutate, nMutate
integer, allocatable :: list(:)
integer listL
integer nSwitchOn, nSwitchOff
integer iSwitch
integer iIndv, nIndv

!integer, allocatable :: offGenes(:)
!integer              :: nOffGenes

integer :: HighestGene
integer :: nOnGenes

integer element


nIndv = size(individuals)

do iIndv=1, nIndv

   HighestGene = ubound(individuals(iIndv)%genome,1)
   
   ! Split the mutation process up into 2 parts: first, determine how many genes will be mutated and
   ! then select this amount of random genes and mutate them. In principle, this split would not be
   ! necessary and you could have one single loop, like
   !     do iGene=1,nGene
   !          if (rnd<mutationProb) mutate( gene(iGene) )
   !     end
   ! However, if you have lots of genes (i.e. clusters) then this process would turn on far too
   ! many genes even for small mutatione probabilities. I (tk) therefore regard it as more useful
   ! only to mutated at most nMaxClusters, or even individual%nOnGenes. Both of them provide THE
   ! genome-scale for mutations. The first (nMaxClusters) nevertheless tends to switch on more
   ! and more genes as the algorithm proceeds.
   !
   !
   ! A) Determine how many genes shall be mutated
   nMutate=0
   do iGene=1,individuals(iIndv)%nOnGenes ! OR: do iGene=1,GA%nMaxClusters
     if (rnd(GA)<mutationProbability) nMutate = nMutate + 1
   enddo
   
   ! B) mutate genes
   allocate(list(GA%nAlwaysOn:HighestGene))
   list  = (/ (iGene,iGene=GA%nAlwaysOn,HighestGene) /)
   listL = size(list)
   do iMutate=1,nMutate
     call take_random_element_out_of_list(list,listL,iGene,GA) ! iGene now holds the index of the gene that is going to be mutated
     individuals(iIndv)%genome(iGene) = .not. individuals(iIndv)%genome(iGene)      ! mutate gene
   enddo
   deallocate(list)

   ! correct pairs
   call correct_pairs(GA%selectSmallestPairsOnly, GA%nPairClusters, individuals(iIndv:iIndv), GA%nAlwaysOn)
   
   ! set individual information:
   call set_genome_information( individuals(iIndv)%genome, individuals(iIndv)%nOnGenes, individuals(iIndv)%onGenes, condition=.true.)

enddo ! iIndv

end subroutine mutate_individuals
!****************************************************************************************************




!****************************************************************************************************
! Purpose:
! Check if too many genes are turned on. If yes, the respective individuals is either setup from
! scratch (reSetup==.true.) or some of the onGenes are simply turned off (reSetup==.false.).
! (reSetup==.false. had already been realised in mutate_individual in a prior version of the code.
! It turns out, however, that the code is given to never ever leave the maximum amount of genes
! activated, if it ran into this by mutation. So, I guess (tk), that reSetup==.true. is the better
! alternative. It simply means that if an individual has too many genes turned on, it is not fit
! for survival. We kill it and do a "virign birth", i.e. setup anew.
!
subroutine check_individuals_for_too_many_onGenes( individuals, GA, virginBirth )
type(t_individual), intent(inout) :: individuals(:)
type(GA_variables), intent(inout) :: GA
logical, intent(out)              :: virginBirth

integer iIndv, nIndv
integer nSwitchOff, nOnTrueGenes, nTrueGenes, iSwitch, iGene
integer, allocatable :: onTrueGenes(:)

nIndv = size(individuals)

virginBirth = .false.
do iIndv=1,nIndv

  if (individuals(iIndv)%nOnGenes <= GA%nMaxClusters) cycle ! everything's OK
  
  ! if we are here: too many clusters are activated
  virginBirth = .true.

  ! setup the respective individual from scratch ("virgin birth")
  deallocate(individuals(iIndv)%genome, individuals(iIndv)%onGenes)
  call setup_population( individuals(iIndv:iIndv), 1, GA%nMaxClusters, GA%nTotClusters, GA%nAlwaysOn, GA, warn=.false.)
  call correct_pairs(GA%selectSmallestPairsOnly, GA%nPairClusters, individuals(iIndv:iIndv), GA%nAlwaysOn)
enddo

end subroutine check_individuals_for_too_many_onGenes
!****************************************************************************************************



!****************************************************************************************************
! Purpose:
! Check if the individuals that are given via the interface equal the dead individual.
! If yes, the respective individual of the set individualsA is setup anew ("virgin birth")
!
subroutine check_individuals_for_equality( individualsA, GA, individualsB, virginBirth )
type(t_individual), intent(inout) :: individualsA(:)
type(t_individual), intent(in)    :: individualsB(:)
type(GA_variables), intent(inout) :: GA
logical, intent(out)              :: virginBirth

integer iIndvA, nIndvA
integer iIndvB, nIndvB

nIndvA = size(individualsA)
nIndvB = size(individualsB)

virginBirth = .false. 

do iIndvA=1,nIndvA

  if (individualsA(iIndvA)%cvs == -1 .or. isnan(individualsA(iIndvA)%cvs)) &
       cycle   ! if there is no dead individual (because no individual has died yet) the cvs 
               ! is -1 and we skip this routine

  do iIndvB=1,nIndvB

    if (individualsB(iIndvB)%cvs == -1 .or. isnan(individualsB(iIndvB)%cvs)) &
          cycle   ! if there is no dead individual (because no individual has died yet) the cvs 
                  ! is -1 and we skip this routine

    if (all(individualsA(iIndvA)%genome .eqv. individualsB(iIndvB)%genome)) then
      virginBirth = .true.
      deallocate(individualsA(iIndvA)%genome, individualsA(iIndvA)%onGenes)
      call setup_population(individualsA(iIndvA:iIndvA), 1, GA%nMaxClusters, GA%nTotClusters, GA%nAlwaysOn, GA, warn=.false.)
      call correct_pairs(GA%selectSmallestPairsOnly, GA%nPairClusters, individualsA(iIndvA:iIndvA), GA%nAlwaysOn)
    endif

  enddo 

enddo

end subroutine check_individuals_for_equality
!****************************************************************************************************




!****************************************************************************************************
! Purpose:
! Given a "genome", return the "number" of genes that are switched on or off ("condition"), and return
! a "list" of those genes.
!
! E.g. Genome = ttft
!      condition=t   ==>  number = 3,  list = { 0,1,3 }
!      condition=f   ==>  number = 1,  list = { 2 }
!
subroutine set_genome_information(genome,number,list,condition)
logical, intent(in)  :: genome(0:)
integer, intent(out) :: number
integer, intent(out) :: list(:)
logical, intent(in)  :: condition

integer :: idx( size(list) ), i

idx = (/ (i,i=0,size(genome)-1) /) 
list  = -1

number          = count(       genome(0:) .eqv. condition)
list(1:number)  = pack ( idx , genome(0:) .eqv. condition)

end subroutine set_genome_information
!***************************************************************************************************




!***************************************************************************************************
! This routine takes two parents (selected previously) and creates a new child by mating with
! n-point crossover. 
SUBROUTINE mate_2_parents(Par1, Par2, Child, GA)
type(t_individual), intent(in) :: Par1, Par2
type(t_individual), intent(out):: Child
type(GA_variables),intent(inout) :: GA ! Need this just for the random number generator seed

integer iGene, nGene

if(.not. associated(Child%genome)) stop "Child was not allocated before mate_2_parents in ga_fitting.f90"

nGene = size(Par1%genome) ! Number of parameters in the fit; i.e., number of genes

Child%genome(0:GA%nAlwaysOn-1) = .true.  ! the first ones are always activated
do iGene = GA%nAlwaysOn, nGene-1   ! the gene with index nAlwaysOn is the first that can be OFF (list of genes
                                   ! starts at 0!)
   if (rnd(GA)<.5) then; Child%genome(iGene) = Par1%genome(iGene)
                   else; Child%genome(iGene) = Par2%genome(iGene); endif
end do

! set child information:
Child%age = 0
call set_genome_information(child%genome,child%nOnGenes,child%onGenes,condition=.true.)

END SUBROUTINE mate_2_parents
!***************************************************************************************************




!***************************************************************************************************
! This routine takes a list of cross validation scores and returns a "roulette wheel" for selecting
! parents to be used in the mating. The relative weight of each individual on the roulette wheel is
! evaluated by the weight function (contained in this subroutine). Put your own weight function there.
SUBROUTINE map_CVS_to_roulette(CV, roulette)

real(dp), intent(in)  :: CV(:,:)         ! { individual#, population# }
real(dp), intent(out) :: roulette(:,:)   ! { individual#, population# }

integer   :: nCV           ! Number of CV scores for each population (= number of individuals in each pop)
integer   :: iCV           ! counter
real(dp)  :: w(size(CV,1)) ! temporary storage for the weights (size = number of individuals in each pop)

integer nPop, iPop

nCV  = size(CV,1)       ! = number of individuals in each pop
nPop = size(roulette,2)

do iPop=1,nPop

   ! determine the weights:
   w = weights(CV(:,iPop),nCV)
   ! sum the weights up (cumulative function) to give roulette
   roulette(:,iPop) = (/ (sum(w(1:iCV)),iCV=1,nCV) /)

enddo ! iPop

CONTAINS

function weights(CV,nCV)
  real(dp), intent(in)  :: CV(:)
  integer               :: nCV
  real(dp)              :: weights(nCV) ! function value

  real(dp)              :: Z ! Normalisation constant


  !--------------------------------------------------------------------------------
  ! place your kind of weight function here
  !--------------------------------------------------------------------------------

  ! GH+TK:
  ! * Weight function used here:
  !   (beta*CV_max - CV_i)
  ! * beta factor: 1 < beta
  !   - influences how strongly better parents are selected over less 
  !     fit ones to go into the mating routine
  !   - the smaller beta, the larger the difference


  integer, parameter :: beta = 1.05_dp ! Should be bigger than 1; old ce_fit used 2; gus proposed 1.05
                                       ! mind: if all CV happen to be rather small (1E-6) or equal then having a
                                       ! beta ~ 1 results in all weights == 0, therefore Z = NaN
  real(dp) CVmax  
  
  CVmax   = maxval(CV)
  weights = (beta*CVmax - CV)
  !--------------------------------------------------------------------------------

  ! check for NaN
  if (all(isnan(weights))) then
    weights = 1._dp
    print *, "weights=",weights
  endif
  where ( isnan(weights) )
    weights = 0._dp
  end where

  ! check for all == 0
  if (all(weights==0)) weights = 1._dp

  ! normalise the weights:
  Z = sum(weights)
  weights = weights / Z

end function weights

END SUBROUTINE map_CVS_to_roulette
!***************************************************************************************************





!***************************************************************************************************
! This routine uses the roulette (relative fitness scores of the parents) and selects parents
SUBROUTINE make_parent_list(roulette,parents,currpop,GA)
real(dp), intent(in) :: roulette(:,:)  ! roulette wheel
type(t_GA_parents), intent(out) :: parents(:) ! { child# }
integer, intent(in)  :: currpop
type(GA_variables)   :: GA           ! need rnd seed and nChildren

integer iChild, nChild, loc(1), iPar
integer iPop, nPop

nChild = GA%nChildren
nPop   = size(roulette,2)

!print*, "roulette=",roulette

do iChild = 1, nChild  ! loop over all children that we want to beget
   do iPar = 1, 2      ! each child needs 2 parents

      if (rnd(GA)<GA%childWithOtherPopulation) then ! make child with different population?
        iPop = ceiling(rnd(GA)*nPop) ! yes: choose population
      else  
        iPop = currpop ! no: use current population
      endif
      
      loc = minloc(roulette(:,iPop),rnd(GA)<roulette(:,iPop))  ! draw a  random number, rnd(GA) and compare it to
                                               ! the roulette wheel.
                                               ! e.g.     rnd = 0.6
                                               !      roulette = 0.1, 0.2, 0.4, 0.7, 1
                                               ! ==> choose #3 as parent

      if (iPar==1) then
        parents(iChild)%father = loc(1)
        parents(iChild)%fpop   = iPop
      else
        parents(iChild)%mother = loc(1)
        parents(iChild)%mpop   = iPop
      endif

   enddo
end do

END SUBROUTINE make_parent_list
!***************************************************************************************************





!***************************************************************************************************
! This routine initializes a matrix A and a vector b that represent the inequality constraints on
! the fitting. In our case, these constraints are the Garbulsky-Ceder constraints. Ax <= b (and x in
! this case is just the vector of J's that we are trying to find.) There are three separate G-C
! constraints .
! 
SUBROUTINE setup_quadprog_constraints(Gamma, PI, lowE_vs, lowE_dist, gsl_adj, gsl_mix, gsl_dist, delta, A, b)
real(dp), intent(in) :: Gamma(:)     ! fit values (e.g. energies, bulk moduli, etc...)
real(dp), intent(in) :: PI(:,0:)     ! PIs of fit structures (1st index) regarding fit parameter (2nd index)
integer, intent(in)  :: lowE_vs(:)   ! { structure# } vector subscript identifying which structure is lowest energy
                                     ! at the concentration of the current structure
real(dp), intent(in) :: lowE_dist(:) ! { structure# } for each structure: its distance to the structure that is lowest in energy
                                     ! at the concentration of the current structure
integer,  intent(in) :: gsl_adj(:,:) ! { adjacent#, structure# }, vector subscript identifying which structures are adjacent ground states
real(dp), intent(in) :: gsl_mix(:,:) ! { adjacent#, structure# }, for each structure: mixture of those ground states 
real(dp), intent(in) :: gsl_dist(:)  ! { structure# }, distance to DFT convex hull
real(dp), intent(in) :: delta(:,:)   ! { constraint#, structure# } 

real(dp), pointer    :: A(:,:), b(:) ! see comment below

real(dp), pointer    :: Apart(:,:,:) ! Temporary holder for values for A matrix
real(dp), pointer    :: bpart(:,:)   ! Temporary holder for values for b vector

integer iStr, nStr, iAdj, nAdj

! mind:
!   size(Gamma) = number of input structures for this fit
!   size(PI,2)  = number of fit parameters (formerly aka figures)
!   size(delta) = number of constraints
!
! A = 2*size(delta) "row chunks": 
!
!     /                    \                               
!     |  Apart1            |                                     
!     | . . . . . . . . .  |                                                
!     | -Apart1            |                                     
!     | .................  |                                               
! A = |  Apart2            |                                     
!     | . . . . . . . . .  |                                               
!     | -Apart2            |                                     
!     | .................  |                                               
!     | etc                |                                 
!     \                    /                             
!
! size(delta) = number of abs-constraints
! 2*          = abs has to be resolved in +/-
!
!
!
!
!     /                    \    /         \                     
!     |  bpart1            |    | delta1  |                            
!     | . . . . . . . . .  |    |         |                           
!     | -bpart1            |    | delta1  |                            
!     | .................  |    |         |                                
! b = |  bpart2            | +  | delta2  |
!     | . . . . . . . . .  |    |         |                                
!     | -bpart2            |    | delta2  |                             
!     | .................  |    |         |                                
!     | etc                |    | etc     |                     
!     \                    /    \         /               
!
! mind: delta = delta(sigma)!

integer nGamma, nCstrt ! number of supplied fit structures,
                         ! number of constraints
integer i                ! counter/loop index
integer x,y
integer rows(4)          

real(dp) mixPI, iCl, nCl

!write(*,'(A$)') "entering setup_quadprog_contraints"

          ! vv 2 * number of inputs for each delta, number of fitting parameters
allocate(A(2*size(delta, 1)*size(Gamma),0:size(PI,2)-1 ))
allocate(b(2*size(delta, 1)*size(Gamma)))

allocate(Apart(size(Gamma),0:size(PI,2)-1,size(delta,1)))
allocate(bpart(size(Gamma),               size(delta,1)))

nGamma = size(Gamma)
nCstrt = size(delta,1)
nStr=nGamma
nAdj=size(gsl_mix,1)

A     = 0._dp;
Apart = 0._dp; 

!--------------------------------------------------------------------------------
! Setup your constraints here:

!......................................
! Condition matrix A
Apart(:,0:,1) = PI(:,0:)                  ! 1st GC constraint
Apart(:,0:,2) = PI(:,0:)-PI(lowE_vs,0:)   ! 2nd GC constraint
Apart(:,0:,3) = PI(:,0:)                  ! 3rd GC constraint (dist to GSL) 
do iStr=1,nStr; do iAdj=1,nAdj;          
  Apart(iStr,0:,3) = Apart(iStr,0:,3) - PI( gsl_adj(iAdj,iStr) , 0: ) * gsl_mix( iAdj, iStr )
enddo; enddo; ! no forall possible!
!......................................
! Condition vector b
bpart(:,1) = -Gamma(:)                    ! 1st GC constraint                  (positive! >0)        
bpart(:,2) = lowE_dist(:)                 ! 2nd GC constraint                  (positive! >0)      
bpart(:,3) = gsl_dist(:)                  ! 3rd GC constraint (dist to GSL)    (positive! >0)
!......................................

!write(991,'(A,1000(59(F4.1,1x),/))') "A1=",apart(:,2,1)
!write(992,'(A,1000(59(F4.1,1x),/))') "A2=",apart(:,2,2)
!write(993,'(A,1000(59(F4.1,1x),/))') "A3=",apart(:,2,3)

!--------------------------------------------------------------------------------

! concatenate the "part"-matrices/vectors together 
!   Apart1, Apart2, ... ==> A
!   bpart1, bpart2, ... ==> b
! in order to be ready for the constraints 
!   Ax <= b
! when fitting the x (ECIs)

do i = 1, nCstrt

   ! i represents Apart_i and bpart_i, resp.

   rows(:) = (/ (i-1)*2*nGamma +1 , &    ! row where  Apart_i/ bpart_i starts in A/b ("/" means "or")
                (2*i-1)*nGamma    , &    ! row where  Apart_i/ bpart_i stops  in A/b
                (2*i-1)*nGamma +1 , &    ! row where -Apart_i/-bpart_i starts in A/b
                 i*2*nGamma           /) ! row where -Apart_i/-bpart_i stops  in A/b

   A(rows(1):rows(2),0:)   =  Apart(:,0:,i)
   A(rows(3):rows(4),0:)   = -Apart(:,0:,i)

   b(rows(1):rows(2))     =  bpart(:,i)+delta(i,:)
   b(rows(3):rows(4))     = -bpart(:,i)+delta(i,:)

end do

!write(994,'(A,1000(59(F4.1,1x),/))') "A=",a(:,2)

deallocate(Apart, bpart)

!write(*,*) "done"

END SUBROUTINE setup_quadprog_constraints
!***************************************************************************************************




!****************************************************************************************************
! Purpose:
! Get all the information needed to make a CE fit via Ole's quadprg routine. The errors of the fit
! may be evaluated. More importantly, the errors of the predictions is used to calculate the 
! cross validation score (CVS) for each individual.
!
! --> IN
!               CE  :  variables for the cluster expansion in general. Here: number of observables and some 
!                      general variables (eps) for the fitting
!       structures  :  ALL input structures
!          subsets  :  the subsets used for fitting and prediction
!      individuals  :  the current individuals with their genes (i.e. clusters) switched on or off
!             iObs  :  observable#
!  clustersDamping  :  the damping for each clusters = c * avgR^lambda
!                     
! <-- OUT             
!    failIndividual : the index of the individual that failed to produce any fits
!             "CVS" : the cross validation score for the individuals will be contained in individuals
!
subroutine make_ceFit_and_calculate_CVS_for_individuals(CE, structures, subsets, individuals, iObs, clustersDamping, failIndividual)
! interface:
type(CE_variables), intent(inout) :: CE
type(crystal), intent(in)         :: structures(:)
type(t_subset), intent(in)        :: subsets(:)    
type(t_individual), intent(inout) :: individuals(:) ! { individual# } only for ONE generations, intent(out) is the CVS
integer, intent(in)               :: iObs           ! which observable are we looking at?
integer, intent(out)              :: failIndividual
real(dp), intent(in)              :: clustersDamping(0:) ! { cluster# }, for each cluster the corresponding damping
! additional:
real(dp), pointer           :: fErrors(:,:)        ! { prediction structure#, subset# }: fitting errors
real(dp), pointer           :: J(:)                ! { cluster# }
real(dp), pointer           :: PI(:,:)             ! { structure#, cluster# } PI correlation matrices for the fitting set and for the prediction set
real(dp), pointer           :: matA(:,:), vecB(:)  ! constraint matrices for the quadprog

integer nStr, iStr, firstSS, iSS, nSS, iIndv, nIndv

nStr        = size(structures)  ! number of structures (total)
firstSS     = lbound(subsets,1)      
nSS         = ubound(subsets,1)   
nIndv       = size(individuals)



allocate( fErrors( maxval(subsets(:)%nFit ),nSS) )        ! only used for failsafe

do iIndv=1,nIndv

    ! get PI-matrix (according to genome) for all structures
    allocate(PI(nStr, 0:individuals(iIndv)%nOnGenes-1))
    do iStr=1,nStr
      PI(iStr,0:) = structures(iStr)%PI( individuals(iIndv)%onGenes( 1:individuals(iIndv)%nOnGenes ) ); 
    enddo

    if (associated( individuals(iIndv)%J )) individuals(iIndv)%J => null()
    allocate( individuals(iIndv)%J(0:individuals(iIndv)%nOnGenes-1,nSS) )
    if (associated( individuals(iIndv)%subsetFitSuccessful )) individuals(iIndv)%subsetFitSuccessful => null()
    allocate( individuals(iIndv)%subsetFitSuccessful(nSS) )

    !.................................... SUBSETS ................................................................
    do iSS=firstSS,nSS

      ! failsafe check: is the constant term switched on?
      if (.not.any(individuals(iIndv)%onGenes( 1:individuals(iIndv)%nOnGenes )==0)) then
           print *, iindv,iss,"nOn=",individuals(iindv)%nongenes
           print *, individuals(iindv)%ongenes(1:individuals(iindv)%nongenes )
           stop "ERROR: bad programming: constant term not switched on"
      endif

      ! set constraints
      call setup_quadprog_constraints( structures( subsets(iSS)%fitList )%fitenergy, &  ! IN:  fit the structures in fitList
                                       PI( subsets(iSS)%fitList ,:), &                  ! IN:  { structure#, on-cluster# }
                                       lowE_vs  =subsets(iSS)%StrLowE, &                ! IN:  { structure# }: index of the struc with lowest energy
                                       lowE_dist=subsets(iSS)%StrLowEDist(:), &         ! IN:  { structure# }: distance to that structure
                                       gsl_adj  =subsets(iSS)%GSofPlane, &              ! IN:  { gs#, structure# }: indices of the adjacent ground-states
                                       gsl_mix  =subsets(iSS)%GSmixture, &              ! IN:  { gs#, structure# }
                                       gsl_dist =subsets(iSS)%GSLdist, &                ! IN:  { structure# }: distance to the convex hull
                                       delta    =subsets(iSS)%delta(:,:,iObs), &        ! IN:  { constraint#, fit structure# }
                                       A=matA, b=vecB )                                 ! OUT: constraint matrix and vector

      ! fit
      call make_ceFit_with_constraints( structures( subsets(iSS)%fitList )%fitenergy, &   ! IN:  fit the strucutres in fitList
                                        structures( subsets(iSS)%fitList )%weight, &      ! IN:  their weight
                                        PI( subsets(iSS)%fitList(:),:), &                 ! IN:  { structure#, on-cluster# }
                                        clustersDamping( individuals(iIndv)%onGenes( 1:individuals(iIndv)%nOnGenes )), & 
                                                                                          ! IN:  { on-cluster# }, the damping parameters for the on-clusters
                                        matA, vecB, &                                     ! IN:  constraints
                                        individuals(iIndv)%J(0:,iSS) , &                  ! OUT: { on-cluster# }: the effective cluster interactions
                                        CE, &                                             ! IN:  CE data
                                        individuals(iIndv)%subsetFitSuccessful(iSS) )     ! OUT: successful fit?
      ! failsafe:
      !if (individuals(iIndv)%subsetFitSuccessful(iSS)) then
      !
      !  ! get errors of FIT for this subset and this observable
      !  call get_errors_of_CE( structures( subsets(iSS)%fitList )%fitenergy, &      ! IN:  observable
      !                         individuals(iIndv)%J(0:,iSS), &                      ! IN:  effective cluster interactions
      !                         PI( subsets(iSS)%fitList ,0:), &                     ! IN:  correlations
      !                         fErrors(:,iSS) )                                     ! OUT: { fit-structure# }
      !  ! if this failsafe triggers, it means that seteps() in front of quadprg uses a different epsilon than here
      !  if ( any(abs(fErrors( 1:subsets(iSS)%nFit ,iSS)) > subsets(iSS)%delta(1,:,iObs)+0.00001 )) &
      !     stop "ERROR: bad programming in the fitting: errors > delta1 !"
      !
      ! Remark: We would need something more sophisticated here, since small errors can happen (cf. the fittingErrors.out file,
      ! where those structures that have delta<error are marked with an "-X-"...)
      !
      !endif ! successful fit?

       deallocate(matA,vecB)

    enddo ! iSS
    !................................. END SUBSETS ...................................................................

    ! see (*) below
    !    ! check if at least 1 fit succeeded
    !    if (all(individuals(iIndv)%subsetFitSuccessful.eqv..false.)) then
    !      ! no (subset) fit succeeded for the individual iIndv
    !      failIndividual = iIndv   ! iIndv failed!
    !      return / or: stop "No fit."
    !    end if
    !    ! if we are here, at least 1 subset yielded a fit for individual iIndv
    
    individuals(iIndv)%CVS = calc_cvs_for_individual( structures(:)%fitenergy, &                ! IN:  { observable# }
                                                      subsets(1:), &                            ! IN:  { successful subset# }
!                                                      individuals(iIndv)%subsetFitSuccessful, & ! IN:  which of those are successful?
                                                      (/ (.true.,iSS=firstSS,nSS) /), &         ! take all subsets into account, see (*)
                                                      individuals(iIndv)%J(0:,:), &             ! IN:  { on-cluster#, subset# } J for each subset
                                                      PI(:,0:) )                                ! IN:  { structure#, on-cluster# } PI

    ! Remark (*) by TK:
    ! The question is: do we want to calculate the CVS only by using those subsets that returned a fit? 
    ! If we do use only the "successful" subsets, then the fact, that the current individual was not able to
    ! produce a fit for some of the subsets (the "unsuccessful" ones), does not enter the final CV score. 
    ! In fact, I dare say that in this case, the CVS is always smaller for an individual that produces 
    ! fewer fits. Thus, I don't think that's a good idea: wouldn't we breed individuals that make fewer fits?
    ! Let's take into account ALL subsets, and the CVS of the individual will more truly show its fitness. 
    ! Those subsets that don't return a successful fit, will give a really BAD prediction, therefore also increasing
    ! the CVS for that individual.
    ! (GH agrees)
  
    deallocate(PI)

enddo ! iIndv

failIndividual = 0  ! no individual failed to make a fit

deallocate(ferrors)

end subroutine make_ceFit_and_calculate_CVS_for_individuals
!****************************************************************************************************




!****************************************************************************************************
function calc_cvs_for_individual( observable, subset, subsetMask, J, PI, predictionErrors  )
real(dp)  :: calc_cvs_for_individual
real(dp), intent(in)            :: observable(:)          ! { structure# }
type(t_subset), intent(in)      :: subset(:)              ! { subset# }
logical, intent(in)             :: subsetMask(:)          ! { subset# } which subsets of subset(:) to use?
real(dp), intent(in)            :: J(0:,:)                ! { cluster#, subset# }
real(dp), intent(in)            :: PI(:,0:)               ! { structure#, cluster# }
real(dp), optional, intent(out) :: predictionErrors(:,:)  ! { structure#, subset# }

integer :: nPredictions( size(subset) )
integer iSS, nSS

real(dp), allocatable :: predictionErrors_(:,:)

nSS          = size(subset)      ! number of subsets passed in
nPredictions = subset(:)%nPred   ! copy of the number of predictions in each subset

! get memory for the prediction errors
allocate(predictionErrors_(maxval(subset(:)%nPred),nSS))
         predictionErrors_ = 0._dp

! determine the predictions errors for each subset
do iSS=1,nSS
    ! get PREDICTION errors for this subset and this observable
         call get_errors_of_CE( observable( subset(iSS)%predList ), &     ! IN:  { observable# in prediction list }
                                J(0:,iSS), &                              ! IN:  { on-cluster# } effective cluster interactions
                                PI( subset(iSS)%predList ,0:), &          ! IN:  { structure# in prediction list, on-cluster# }
                                predictionErrors_(:,iSS) )                ! OUT: { prediction-structure# }
enddo   
! fitness of iIndv wrt iObs: 
! use the prediction errors pErrors to evaluate the CVS
calc_cvs_for_individual = calc_cvs_from_errors( predictionErrors_(:,:) , nPredictions, subsetMask ) 
                                                                         ! IN:  { prediction-structure#, subset# }
                                                                         ! IN:  number of predictions for each subset
                                                                         ! OUT: cross validation score for individual iIndv wrt iObs

if (present(predictionErrors)) then      ! *argh*, ifort has a problem if we write that in ONE if-statement
if (size(predictionErrors)>0) then
  predictionErrors = predictionErrors_
  return ! get out of the function and don't deallocate predictionErrors
endif
endif

! deallocate predictionErrors
deallocate(predictionErrors_)

end function calc_cvs_for_individual
!****************************************************************************************************





!****************************************************************************************************
function calc_cvs_from_errors( errors, nerr, subsetMask )
real(dp)                      :: calc_cvs_from_errors
real(dp),       intent(inout) :: errors(:,:)  ! { error#, subset# }: residuals for all structures
integer,        intent(in) :: nerr(:)      ! { subset# }: for each subset: number of predictions
logical, optional, intent(in) :: subsetMask(:)

real(dp) :: CVS2  ! squared CVS

integer iSS, nSS
real :: NaN = 0.  ! want to generate NaN by the fraction 0./0. Although not a problem in ifort, gfortran gives compile error.
                  ! so, have to fool the compiler by "hiding" the 0./0. by the use of NaN/NaN

nSS = size(errors,2)   ! number of subsets

if (all(subsetMask.eqv..false.)) then ! no good subset!
  calc_cvs_from_errors = NaN / NaN
endif

! failsafe: all errors above nerr(iSS) must be zero!
do iSS=1,nSS; if ( any(errors(nerr(iSS)+1:,iSS)/=0) ) stop "ERROR: bad programming in calc_cvs_from_errors"; enddo

if (present(subsetMask)) then                              ! follow the comments here: 1), 2a,b), 3):
  CVS2     = real( sum( &                                         ! 2a) sum the result (an array of dimension 1 = subset#)  ...
                        sum( errors(:,:)**2 , dim=1 ), &   ! 1) sum squared errors along dimension 1 (error#)
                   subsetMask.eqv..true.) , dp) &                 ! 2b) ... where the subset mask equals T
             / real( sum( nerr(:), subsetMask.eqv..true. ) , dp )       ! 3) divide the total sum by the total number of errors (only where mask equals T)
else
  CVS2     = real( sum( errors(:,:)**2 ) , dp) / real( sum( nerr(:) ) , dp )
endif
calc_cvs_from_errors = sqrt( CVS2 )

end function calc_cvs_from_errors
!****************************************************************************************************



!SUBROUTINE calculate_gamma(PI_matrix,lambda)
!
!use vector_matrix_utilities, only : matrix_inverse
!
!real(dp), intent(in):: PI_matrix(:,:)
!integer :: err
!real(dp), allocatable :: mean(:),tester(:), trans(:,:), invtrans(:,:), tracematrix(:,:), test(:,:)
!real(dp), allocatable :: deviation(:,:), deviation_squared(:,:), temporary(:), covar_list(:),covar_matrix(:,:)
!integer :: i, j, k, n,l,  number_of_structures, num_of_PIs, counter = 0, indexone, indextwo
!real(dp) :: running_total, sum, secondterm, firstterm
!real(dp), intent(out):: lambda
!
!num_of_PIs = size(PI_matrix,2)
!
!number_of_structures = size(PI_matrix,1)
!
!allocate(mean(num_of_PIs),tester(num_of_PIs))
!
!allocate(covar_list(num_of_PIs**2),trans(num_of_PIs,num_of_PIs),invtrans(num_of_PIs,num_of_PIs),test(num_of_PIs,num_of_PIs),tracematrix(num_of_PIs,num_of_PIs))
!
!allocate(deviation(num_of_PIs,number_of_structures),deviation_squared(num_of_PIs**2,number_of_structures),temporary(number_of_structures),covar_matrix(num_of_PIs,num_of_PIs))
!
!!calculate mean vector
!do i=1,num_of_PIs
!   running_total=0
!   do j=1,number_of_structures
!      running_total = running_total + PI_matrix(j,i)
!   enddo
!   mean(i) = running_total/real(number_of_structures,dp)
!
!enddo
!
!!calculate the deviation matrix
!do k=1,num_of_PIs
!
!   do j=1,number_of_structures
!      deviation(k,j) = PI_matrix(j,k) - mean(k)
!   enddo
!
!enddo
!
!
!!calculate the deviation squared matrix
!counter = 0
!do n=1,num_of_PIs
!   
!   do k=1,num_of_PIs
!      do j=1,number_of_structures
!         temporary(j) = deviation(n,j) * deviation(k,j) 
!      enddo
!      counter = counter + 1
!      do l=1,size(temporary,1)
!         deviation_squared(counter,l) = temporary(l)
!      enddo
!
!         
!   enddo
!
!enddo
!
!
!
!sum=0.0
!do k=1,size(deviation_squared,1)
!   sum=0
!   do j=1,size(deviation_squared,2)
!      sum = sum + deviation_squared(k,j)
!   enddo
!   covar_list(k)= sum/real(number_of_structures-1)
!enddo
!covar_matrix = reshape(covar_list,shape(covar_matrix))
!
!
!
!trans = matmul(transpose(PI_matrix),PI_matrix)
!
!call FINDInv(trans,invtrans,num_of_PIs,err)
!test = matmul(trans,invtrans)
!
!print *, "INV TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!do k=1,size(test,1)
!   do j=1,size(test,2)
!      write(*,'(f9.6,1x)',advance="no")  test(k,j)
!   enddo
!   write(*,*)
!enddo
!
!
!!call matrix_inverse(trans,invtrans,err)
!!if (err) return "ERROR: Problem with inverse matrix function------>>><<<<<<"
!
!
!tracematrix = matmul(invtrans,covar_matrix)
!
!
!
!
!tester = matmul(invtrans,mean)
!
!secondterm = 0
!do i =1,size(tester,1)
!   secondterm = secondterm + tester(i)*mean(i)
!enddo
!
!firstterm=0
!do i =1,size(tracematrix,1)
!   firstterm = firstterm + tracematrix(i,i)
!enddo
!
!lambda = firstterm + secondterm
!
!end subroutine calculate_gamma


!********************************************************************************
! Purpose: calculate correlations for an enumerated structure in gss.out
subroutine get_enumStrcorrs(CE, str, &
                              strN,hnfN,sizeN,n,pgOps,diag,a,b,c,d,e,f,L,labeling) ! reference


type(CE_variables), target               :: CE
type(crystal)                            :: str

integer, intent(in)            :: a, b, c, d, e, f, n, strN, hnfN,sizeN, pgOps, diag(3)
integer, intent(inout)         :: L(3,3)
character(maxLabLength), intent(inout)   :: labeling

integer :: nBas, HNF(3,3)
real(dp), pointer :: mapping_routine_x(:)
integer :: Dmat(3,3)
integer istat, iFg, iRep
integer(kind=1), pointer :: gTab(:,:,:,:)

real(dp) :: invA(3,3)

integer, pointer :: gIndx(:) ! Label ordinal for each atom (position in the labeling)  

nBas = CE%nBas

str%pLV = CE%pLV
str%eps = CE%eps
allocate(str%pBVlab(n*nBas),str%intPt(3,CE%nBas),str%aTyp(n*nBas))

str%pBVlab =  map_labeling_position_2_dvector_index(nUC=n,nD=nBas)

str%intPt = CE%pBas

HNF = 0
HNF(1,1) = a; HNF(2,1) = b; HNF(2,2) = c
HNF(3,1) = d; HNF(3,2) = e; HNF(3,3) = f


call map_enumStr_to_real_space(CE%CErank, n, HNF, labeling, str%pLV,&
       & CE%pBas, CE%eps, str%LV, str%pos, str%spin, gIndx, mapping_routine_x, L, diag, minkowskiReduce=.true.)




! calculate correlations
call create_spin_lookup_table(str,gTab,invA,L,Dmat, determineLD=.true.)
call get_gTab_correlations(CE, gTab, L, Dmat, invA, CE%eqvClusters, str%PI, mustGetG=.true.)


do iFg = 1, size(CE%eqvClusters)
! Don't forget to subtract 1, eqvClusters is a zero indexed array
   do iRep = 1, size(CE%eqvClusters(iFg-1)%Rep)
      !         if (associated(fig(iFg)%Rep(iRep)%gPt))then
      deallocate(CE%eqvClusters(iFg-1)%Rep(iRep)%gPt)
      !            print *, "YELL! The culprid is: gPt",iFg, iRep
      !         endif
   enddo
enddo

deallocate(gTab)
deallocate(mapping_routine_x)

deallocate(gIndx)
deallocate(str%pBVlab)
deallocate(str%aTyp)
deallocate(str%intPt)
deallocate(str%pos)
deallocate(str%spin)

end subroutine get_enumStrcorrs

!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm           
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
!SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
!  IMPLICIT NONE
!  !Declarations
!  INTEGER, INTENT(IN) :: n
!  INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
!  REAL(dp), INTENT(IN) :: matrix(n,n)  !Input matrix
!  REAL(dp), INTENT(OUT) :: inverse(n,n) !Inverted matrix
!  
!  LOGICAL :: FLAG = .TRUE.
!  INTEGER :: i, j, k, l
!  REAL :: m
!  REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix
!	
!  !Augment input matrix with an identity matrix
!  DO i = 1, n
!     DO j = 1, 2*n
!        IF (j <= n ) THEN
!           augmatrix(i,j) = matrix(i,j)
!        ELSE IF ((i+n) == j) THEN
!           augmatrix(i,j) = 1
!        Else
!           augmatrix(i,j) = 0
!        ENDIF
!     END DO
!  END DO
!	
!  !Reduce augmented matrix to upper traingular form
!  DO k =1, n-1
!     IF (augmatrix(k,k) == 0) THEN
!        FLAG = .FALSE.
!        DO i = k+1, n
!           IF (augmatrix(i,k) /= 0) THEN
!              DO j = 1,2*n
!                 augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
!              END DO
!              FLAG = .TRUE.
!              EXIT
!           ENDIF
!!           IF (FLAG .EQV. .FALSE.) THEN
!!              PRINT*, "Matrix is non - invertible"
!!              inverse = 0
!!              errorflag = -1
!!              return
!!           ENDIF
!        END DO
!        IF (FLAG .EQV. .FALSE.) THEN
!           PRINT*, "Matrix is non - invertible"
!           inverse = 0
!           errorflag = -1
!           return
!        ENDIF
!     ENDIF
!     DO j = k+1, n			
!        m = augmatrix(j,k)/augmatrix(k,k)
!        DO i = k, 2*n
!           augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
!        END DO
!     END DO
!  END DO
!  
!  !Test for invertibility
!  DO i = 1, n
!     IF (augmatrix(i,i) == 0) THEN
!        PRINT*, "Matrix is non - invertible"
!        inverse = 0
!        errorflag = -1
!        return
!     ENDIF
!  END DO
!	
!  !Make diagonal elements as 1
!  DO i = 1 , n
!     m = augmatrix(i,i)
!     DO j = i , (2 * n)				
!        augmatrix(i,j) = (augmatrix(i,j) / m)
!     END DO
!  END DO
!	
!  !Reduced right side half of augmented matrix to identity matrix
!  DO k = n-1, 1, -1
!     DO i =1, k
!        m = augmatrix(i,k+1)
!        DO j = k, (2*n)
!           augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
!        END DO
!     END DO
!  END DO
!	
!  !store answer
!  DO i =1, n
!     DO j = 1, n
!        inverse(i,j) = augmatrix(i,j+n)
!     END DO
!  END DO
!  errorflag = 0
!END SUBROUTINE FINDinv






END MODULE ce_fit

