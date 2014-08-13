MODULE ce_types
use num_types
implicit none

type label_pt ! A structure for defining points that become vertices in a figure
   integer  :: label  ! basis point of the parent lattice for this R
   real(dp) :: pos(3) ! Cartesian coordinates of the point
   real(dp) :: len    ! Distance of the point from the origin
endtype label_pt

type label_pt_list
   integer                      :: thisLabel   
   type(label_pt)     , pointer :: thisRList(:)
   integer                      :: thisRListSize
   type(label_pt_list), pointer :: prev, next
end type label_pt_list

type figure ! Stores a CE figure and associated information
   integer  :: nV                   ! number of vertices
   real(dp) :: avgR                 ! average distance of each point from the center of gravity
   real(dp) :: damping              ! = c/2 * avgR^lambda for this cluster

   real(dp), pointer :: vertex(:,:)  ! { coordinate# , vertex# } 3xn_(Vertex) array of points
   integer, pointer  :: gPt(:,:)     ! { gcoordinate#, vertex# } (4xn_V) the vertex points in gTab coordinates
   integer, pointer  :: label(:)     ! { vertex# } basis point of the parent lattice for each vertex
   integer, pointer  :: s(:)         ! { vertex# } chebychev index for a vertex ("s-vector", when taken for all vertices)
                                     ! See the first couple of pages of PRB 49, 8627 (1993)
   type(figure), pointer :: reference(:)
endtype figure

type figRep  ! Representative figure (with all its sym-equivalent brothers)
   type(figure), pointer :: Rep(:) ! Store all the brothers together
   integer     :: nV        
   real(dp)    :: damping ! although "damping" is already included in each Rep(:) it is convenient to store
                          ! it separately
   integer(li) :: count ! Number of representatives (just size of Rep)
                        ! BEWARE: for the pizza figures it is NOT the size of represenations but the size of equivalent figures (take a look
                        ! at generate_multilattice_pizza
   integer(li) :: norm !  Normalization (Ncells x Nrep figs)
   type(figRep), pointer :: reference(:)
endtype figRep

type t_pizza
   integer               :: nCl
   type(figure), pointer :: cluster(:)
end type t_pizza

! This type stores variables for the CS algorithm.
type CS_variables

   real(dp) musigma_start,musigma_finish, delta_musigma   ! Stores starting, ending values of mu/sigma and also how to change these values from iteration to iteration
   integer nSubsets, nInpStruct   ! How many subsets of structures do we want; How many input structures do we want.
   integer musigma_change_mode
   real(dp) J_cutoff     ! Cutoff for keeping Js
   character(6) penaltyfxn

   ! For when you are adding input structures
   character(1) loop   ! Should we loop over values of mu/sigma2?
   character(1) mode   ! What mode are we running in
   character(1) pltone   ! Use p<1 norm?
   logical which_structures_to_add
   integer, pointer :: InpStruct(:),avail_Structs(:)
   real(dp), pointer :: enum_PI_matrix(:,:)
   real(dp), pointer :: input_PI_matrix(:,:)

   
endtype CS_variables


! This structure store various information about each input structure, including the results
! of the CE prediction and ground state search
type crystal

   real(dp) :: LV(3,3)           ! Lattice vectors

   real(dp), pointer :: pos(:,:)  ! position of the basis atoms
   integer, pointer  :: spin(:)   ! occupation variables (= -1, +1, etc...)
   integer, pointer  :: pBVlab(:) ! parent basis vector label
   integer           :: nAtoms    ! total number of atoms
   integer, pointer  :: nTyp(:,:) ! The total number of each (spin) type of atom (index: 1 ~ 1st atomic species, 2, 3...)

   logical          :: atomsMissing
   integer          :: missingTyp(2), missingN

   integer, pointer :: xnTyp(:,:)! The total number of each (spin) type of atom. Only the x-relevant atoms are counted
   integer, pointer :: ATyp(:)   ! The (spin) type of each atom (=rank) in the structure = 1, 2, 3...  
                                 ! (1-to-1 mapping to the spin(:))
   integer, pointer :: fixedOcc(:)  ! { fixed# }: the indices of those atoms that have fixed occupation in this CE
   character(1) :: aBasCoord ! Coordinate system used to describe atomic basis vectors
   character(80) :: title        ! Name of the structure
   real(dp) :: pLV(3,3) ! parent lattice vectors
   real(dp), pointer :: dBasCoord(:,:)
   real(dp), pointer :: intPt (:,:) ! Interior points of the parent lattice
   real(dp) :: a0                ! Lattice constant (not really necessary but...)
   integer  :: weight            ! weight of the structure during the fit
   real(dp) :: energy            ! Energy of the structure (DFT)
   real(dp) :: ce_energy         ! Energy predicted by CE
   real(dp) :: ce_energy_upper   ! Upper bound on energy predicted by CE
   real(dp) :: ce_energy_lower   ! Lower bound on energy predicted by CE
   real(dp) :: ce_refEnergy      ! Total Reference energy as provided in all reference structures
   real(dp) :: ce_upperEnergy      ! Total Reference energy as provided in all reference structures
   real(dp) :: ce_lowerEnergy      ! Total Reference energy as provided in all reference structures
   real(dp) :: fitEnergy         ! the energy that will be used for fitting
   real(dp) :: formationenergy   ! Enthalpy of formation of the structure
   real(dp) :: ce_formationenergy! Enthalpy of formation of the structure as predicted by CE
   real(dp) :: quasibinform      ! Quasibinary formation enthalpy, only used in ternaries
   integer ::  gstleft           ! index in dftgstlist(:) of groundstates to left
   real(dp) :: distanceToInputGSL! distance to ground-state hyperplanes
   real(dp), pointer :: PI(:)
   real(dp) :: eps ! epsilon for finite precision checking
   real(dp) :: conc ! binary concentration, in ternary case: quasibinary conc.
   real(dp) :: con(3) ! ternary concentrations
   real(dp), pointer :: x(:,:)  ! { atomic type, CE# }: concentration of each component type (will
                                ! eventually supercede conc and con(3)
                                ! Note: CE# is for coupled CEs. For a "normal" CE, the second index is always ==1
   real(dp), pointer :: delta(:,:) ! Garbulsky Ceder Constraints: { constraint number, observable }

   logical :: isHelpStructure

   type(crystal), pointer :: reference(:)
   logical                :: automaticReference

   integer, pointer :: HNFlist(:,:,:) ! { entry, entry, HNF# }: for structure comparison
   integer, pointer :: pLabel(:,:)    ! { labeling#, entry }  : for structure comparison

endtype crystal

type surface_variables
  logical  :: isSurf             ! switch: is this a surface CE?

  real(dp) :: unitNormal(3)      ! surface normal vector (user input, then normalised)
  real(dp) :: distT, distB       ! distance of "top" and "bottom" surface from origin (calculated)

!   character(1), pointer :: BasType(:)      ! type of the basis point (bulk "B"/surface "S"/adsorbate "A")
!   integer, pointer      :: BasLayer(:)     ! in which layer is this basis atom?
!   integer, pointer      :: BasLayerAtom(:) ! which atom is it in this layer?
!   integer, pointer      :: BasId(:)        ! lookup table for the site identities of the slab, e.g. atom 1 ~ atom 10 due
!                                            ! to the symmetric slab (and those site always have the same occupancy)
!                                            ! (this does imply that those sites have the same onsite correlations,
!                                            ! see CE_variables%BasEq. However, it does not imply that the sites
!                                            ! with the same onsite correlations, see CE_variables%BasEq, are
!                                            ! identical in the surface slab sense)

!  integer, pointer :: BulkList(:)  ! { bulk atom index, i.e. (1) is the first bulk atom }: which basis atoms are bulk like?
!  integer, pointer :: BulkSpin(:)  ! { dset#, i.e. (1) is the first dset member }: spin variables of the bulk atoms
!  character(1), pointer :: BulkEnumRank(:) ! { dset# }: rank corresponding to spin as a character, enum counting
!                                           ! (i.e. starting at 0 for the first species)

  integer  :: LVpar2normal        ! which surface LV is parallel to the surface normal
!  logical  :: isSymmetric         ! is the slab symmetric?
!  real(dp) :: distSymPlane        ! distance of symmetry plane
                                 
!  real(dp) :: ShiftVec(3)         ! surface shift vector
                                 
!  integer  :: nSlabLayers         ! number of slab layers
!  integer  :: nBulkLayers         ! number of bulk layers
!  integer  :: nSurfLayers         ! number of surface layers


endtype surface_variables

type adsorbate_variables

  logical :: isAds
  logical, pointer  :: isAdsD(:)   ! { dset# }: which of the dset members are adsorbate places?
  
  integer           :: freeSiteRank     ! the rank of the free adsorbate site
  real(dp), pointer :: singleEnergy(:)  ! { (coupling) rank# }: energies of the single adsorbate atoms
!  real(dp), pointer :: cleanEnergy(:)   ! { rank# }: energies of the clean 100% surfaces

end type adsorbate_variables


type t_lat2ref
  real(dp)           :: Mlv(3,3)           ! transformation matrix for the lattice vectors
  real(dp)           :: Mds(3,3), Ods(3)   ! transformation matrix and offset for the dset members
  integer            :: nDel
  integer, pointer   :: dsetDel(:)
end type t_lat2ref

type t_MCg2gref
!  integer, pointer :: grefijk(:,:,:,:)   ! { gref coord (1=i, 2=j, 3=k), gcoord i, gcoord j, gcoord k }: returns the
!                                         ! 1, 2, and 3 coordinate of (i,j,k)_original, i.e., (i,j,k)_ref
!  integer, pointer :: grefd(:,:,:,:)     ! { gcoord i, gcoord j, gcoord k, gcoord iD }: returns the reference
!                                         ! dlabel of (i,j,k,iD)_original
   integer, pointer :: gref(:,:,:,:,:)
end type t_MCg2gref


type t_individual ! stores fitting parameters, cv score, etc, for one individual
   logical, pointer :: genome(:)  ! T = fitting parameter is activated, F = not activated
   integer, pointer :: onGenes(:) ! only the (indices of the) genes that are switched ON
   integer          :: nOnGenes   ! number of entries in the onGenes list
   real(dp), pointer:: J(:,:)     ! { cluster#, subset# } effective cluster interactions
   logical, pointer :: subsetFitSuccessful(:) ! { subset# }, not really used
   real(dp) :: cvs
   integer  :: age
end type t_individual


type t_fitResults
  integer                     :: bestIndividualIndex , bestPopulationIndex
  type(t_individual)          :: bestIndividual

  real(dp), pointer           :: Jav(:)     ! { cluster# } effective cluster interactions (averaged)
  real(dp), pointer           :: Jff(:)     ! { cluster# } effective cluster interactions (full fit, i.e. all structures, no subsets)

  real(dp), pointer           :: fittingError(:)     ! { structure# }
  real(dp)                    :: fittingErrorMax, fittingErrorAv, fittingErrorRMS
  real(dp), pointer           :: predictionError(:,:)  ! { structure#, subset# }
  real(dp)                    :: CVS
end type t_fitResults


type CE_variables             ! This data type will contain all lat.in data
                              ! and will simplify shoving data back and
                              ! forth information between routines
                              ! a lot more general data should be put 
                              ! in there later for convenience and clarity 
   character(30) :: projectname ! Title of CE-project
   character(10) :: title     ! Title of lattice
   character(1) :: CEtype     ! Bulk or surface?
   logical :: scalar          ! Do you need an energy expansion (scalar) only, or something more complicated (multidimensional)?
   integer :: nObservables    ! number of configuration dependant observables, that is, _including_ energy
   logical      :: isBulk 
   character(1) :: CoordSys   ! parent basis vectors / surface normals in Cart. or Direct coords?
   integer :: CErank          ! Binary, ternary, ... ?
   real(dp) :: pLV(3,3)       ! parent lattice vectors

   type(surface_variables)   :: surface   ! surface variables
   type(adsorbate_variables) :: adsorbate ! adsorbate variables

   integer           :: nBas          ! number of parent basis vectors
   integer           :: nBasIneq      ! number of inequivalent basis vectors / onsitecorrs
   integer,  pointer :: digit(:), label(:,:), equivalencies(:), fixedOcc(:)
   integer,  pointer :: BasEq(:)      ! lookup table showing which basis corresponds to which onsite
   real(dp), pointer :: pBas(:,:)     ! parent basis vectors (3xnBas)
   logical,  pointer :: xRelevant(:)  ! specifies if a dset member is relevant
                                      ! for calculating the concentration

!   integer           :: nBasEnum      ! number of parent basis vectors that will be enumerated
!   logical,  pointer :: BasEnum(:)    ! enumerate this basis point?
!   real(dp), pointer :: pBasEnum(:,:) ! the parent basis vectors for which BasEnum=.true.

   real(dp) :: eps            ! accuracy for handling of coordinates
   real(dp) :: xeps= 0.00000001_dp   ! accuracy for handling of concentrations
   real(dp) :: cutoff(10) = 0      ! cutoffs for figure definition

   logical  :: GSS_override      ! go on even if struct_enum.out does not seem to be consistent with CE data
   logical  :: GSS_concOnly      ! perform full GSS search or only such meeting the required concentrations?
   integer  :: GSS_concN         ! element number for the next variable:
   integer,pointer :: GSS_concRange(:,:)  ! min/max/divisor concentration of element GSS_concN in GSS

   logical :: GSS_full        ! perform full GSS search or only "true" multi-naries?
   integer :: GSSspan(2)      ! Starting and ending cell sizes for the ground state search (GSS)
   logical :: GSS_MPIcat      ! after the GSS: do you wan to concatenate the different MPI gss.out.* files to a
                              ! single gss.out?
   logical :: GSS_MPIdelete   ! after the concatenation (see above): delete original gss.out.* files?
   
   integer :: totalstructures ! how many structures are defined in struct.in
   integer :: leftoutStr      ! how many structures at the end of structures.in
                              ! are not used for fitting?
   integer, pointer :: pureList(:) ! the pure element structures

   integer, pointer :: GSlist(:)   ! the structures that are ground states
   logical :: quasibin        ! is the ternary CE quasibinary?
   logical :: cstrain, locstrain
   logical :: concdep         ! concentration dependent ECI?
   integer :: polynome        ! polynomial order of ECI x-dependency -1 (!)
   integer, pointer :: dftgstlist(:) ! indices of the groundstate structures

   real(dp), pointer :: delta(:,:,:)             ! { constraint# , update# , observable# }: Garbulsky-Ceder deltas.
                                                 !   constraint# : 1,2,3 for the 3 Deltas
                                                 !   update#     : we can supply a list of Deltas in fitpar.in
                                                 !   observable# : 1st is energy
   integer, pointer  :: gen4DeltaSwitch(:)       ! switch to delta*new after x Generations

   real(dp) :: pureDelta1     ! GC-constraint 1 for the pure elements (i.e. the total deviation)
   real(dp) :: gstweight      ! weight of groundstates when performing fit
   real(dp) :: ktfit          ! fit parameter
   real(dp) :: cone, ctwo, lambda1, lambda2      ! fit parameters
   
   real(dp) :: qp_eps         ! epsilon for quadprog routine
   real(dp) :: inf            ! infinity for quadprog routine
   integer  :: maxcounts      ! maxcounts for quadprog routine


   type(figure), pointer       :: fig(:) => null()
   type(figRep), pointer       :: eqvClusters(:) => null()
   type(figRep), pointer       :: eqvClustersList(:,:) => null()
   real(dp), pointer           :: J(:), Jff(:), ebars(:)

   type(CE_variables), pointer :: reference(:)
   integer                     :: nreferences
   character(80)               :: lat_file, J_file, cluster_file, str_file
   type(t_lat2ref), pointer    :: l2r(:)

   integer                     :: ncouplings   ! how many CEs will be coupled?
   integer, pointer            :: CEbas(:)     ! { dset# }: to what CE do the dvectors belong to?

endtype CE_variables

type MD_variables !Metadynamics variables
   ! Some general constants first

   real(dp) :: Temperature
   real(dp) :: gaussian_height, gaussian_variance
   real(dp) :: bias_potential_division,bias_potential_division_x,bias_potential_division_y
   integer :: num_of_flips_between_gaussian_adds
   real(dp) :: current_bias_value,proposed_bias_value
   real(dp) :: deltaT, deltaI, minheight
   real(dp),pointer :: Gaussians(:,:)
   integer :: num_of_discrete_cv_values,num_of_discrete_cv_values_x,num_of_discrete_cv_values_y
   integer :: nGaussians
   integer   :: nCV   ! Number of collective variables to be used.
   real(dp), pointer :: bias_potential(:,:)
   real(dp),pointer:: CV(:) ! Values of CV over the course of the simulation

endtype MD_variables

type t_subset
  integer            :: nFit, nPred
  integer,   pointer :: fitList(:)       ! structures in subset (fit)
  integer,   pointer :: predList(:)      ! structures that are predicted
  integer            :: nGS              ! number of ground states
  integer,   pointer :: GSList(:)        ! ground states
  integer,   pointer :: pureList(:)      ! { element } the indices of the pure elements
  real(dp),  pointer :: GSLDist(:)       ! distance to GS convex hull
  integer,   pointer :: GSofPlane(:,:)   ! { GS idx, structure } the GS that span the hyperplane                          
                                         ! above which the structure is floating. 
  real(dp),  pointer :: GSMixture(:,:)   ! { mixture idx, structure }
  integer,   pointer :: StrLowE(:)       ! { structure } the idx of the str that is lowest in energy below the current struc
  real(dp),  pointer :: StrLowEDist(:)   ! { structure } distance to structure lowest in energy (fitenergy)
  real(dp),  pointer :: StrDeltaDist(:)  ! { structure } distance to structure lowest in energy (total energy), used for
                                         ! GC deltas
  real(dp),  pointer :: Delta(:,:,:)     ! { delta idx, structure, observable }
  
end type t_subset

type t_GA_parents
  integer :: father, fpop   ! index of the father and the father's population
  integer :: mother, mpop   ! index of the mother and the mother's population
end type t_GA_parents

type t_GA_continue
  logical       :: yes
  character(80) :: filename
  logical       :: best
  logical       :: dead
  logical       :: all
  logical       :: indv
  integer       :: iIndv
end type t_GA_continue


type GA_variables             ! This data type will contain all GApar.in data
                              ! Not sure if it makes sense NOT to put this into
                              ! CE_variables, but I like it better right now.

   ! General variables ................................................................................
   logical :: simplefit       ! perform simple fit or apply GA?
   type(t_GA_continue) :: continue        ! continue a GA run?

   ! Clusters..........................................................................................
   integer :: nAlwaysOn        ! total number of clusters, number of clusters that are always on
                                          ! (i.e. the Constant and the Onsites, e.g. monolattice 1+1=2)
   integer :: nTotClusters                ! total number of clusters available
   integer :: nMaxClusters                ! maximum number of clusters to activate
                                          ! (i.e. figures + s-vectors, i.e. fitting variables!)
   integer :: nPairClusters               ! number of pair clusters
   logical :: selectSmallestPairsOnly = .true.    
                                          ! Decide whether to treat the selection of the pairs differently than the >2-body-clusters:
                                          ! - If set to T then always the smallest pair clusters only will be selected and not an arbitrary selection
                                          ! of pairs clusters. This is the way the old UNCLE treated the pairs.
                                          ! E.g., genome = T T TTTTFFFFFFFFF TFTFFFTTFFTF with 4 activated pair clusters
                                          !                    <-- pairs --> 
                                          ! - If set to F then the pairs clusters will be treated in the same way as >2-body-clusters, and an
                                          ! arbitrary selection of those clusters can be turned on.
                                          ! E.g., genome = T T TFTFFFFTFFTFF TFTFFFTTFFTF with 4 activated pairs clusters
                                          !                    <-- pairs -->


   ! Populations & Individuals of those populations....................................................
   integer :: nPopulations                ! number of populations
   integer :: nIndividuals                ! size of each population
   integer :: maximumAgeForIndividual     ! the individual will be killed after a certain number of generations
   integer :: nChildren                   ! nr. of children generated in each generation
   integer :: nReplacements               ! nr. of children to replace parents
   real(dp):: mutationProbability
   real(dp):: childWithOtherPopulation    ! probability that one or the other parent is not from the current population
   real(dp):: diversityLimit   ! how diverse must the population be? If the population is less diverse, some countermeasures
                               ! are taken, e.g.,...
   real(dp):: diversityEX, diversityMB    ! ... exchange with other populations (EX), or a mutation burst (MB), 
                                                   ! or a resetup (RS) of the whole population. Those variables give
                                                   ! the respective probabilities for that
   integer :: maxGenerations  ! nr. of generations in each GA run

   ! Cross Validation...................................................................................
   type(t_subset), pointer :: subset(:,:)   ! prediction/fit subsets
   logical                 :: CrossValidationLeaveOneOut  ! enable "leave-one-out CV"?
   integer                 :: CrossValidationK            ! k-fold cross validation...
   integer                 :: CrossValidationRep          ! ... with repetitions

   ! Damping............................................................................................
   real(dp) :: dampingC(0:9), dampingLambda(0:9)  ! the damping parameters for each cluster-type: 0 ... 9-body clusters

   ! Other..............................................................................................
   integer :: rndseed         ! seed for random number generator
   integer :: fixedseed=99    ! seed from input instead from system clock
                              ! 99=clock, all others: fixed at this value
                              ! Important: this is for OUR RNG ONLY, for MC
                              ! we use the fortran intrinsic random_number
   logical :: includeend = .false.  ! include endpoints of conc-range in predictions? (deprecated, must be false)
endtype GA_variables


type CE_results             ! This data type will contain all final and temporary
                            ! results of the CE
   real(dp) :: CVS                     ! resulting CVS
   real(dp) :: EVS                     ! resulting EVS
   real(dp), pointer :: ecival(:)  ! the ECI resulting from fit
   real(dp) :: maxerr, averror, rms    ! some errors
   real(dp) :: maxerrpred, averrorpred, rmspred ! more errors for predictions
   real(dp), pointer :: jconst(:)        ! 
   real(dp), pointer :: jcorr(:,:)         ! Just nicer access to final eci
   real(dp), pointer :: mbeci(:,:)       !
   real(dp), pointer :: jconst_first(:)  !
   real(dp), pointer :: jcorr_first(:,:)   !
   real(dp), pointer :: mbeci_first(:,:) !
   integer, pointer :: bestfitset(:)     ! array containing the oprimized figureset
endtype CE_results

type t_list
  integer          :: nList
  integer, pointer :: list(:)
end type t_list

type MC_variables ! variables for Monte Carlo
   ! Some general constants first
   integer :: maxBoxSizeForCalcTotalE             = 50
   integer ::time2printSpeedInfo                 = 30
   integer :: timeAvgStep2printProgressBar        = 600
   integer :: SIsteps                             = 100000

   ! the MC cell
   integer(si), pointer :: site(:,:,:,:)  ! { i, j, k, d }: (i,j,k)=unit cell, d=dset member; 
                                          ! contains the spins of the sites

   ! geometry of cell
   integer  :: boxSize(3)                      ! Size of the simulations cell (times unit cell in 3 directions)
   integer  :: simOffset(3), simSize(3)        ! Restricted simulation cell: here, within "boxsize", the actual simulation takes place
                                               ! This actually only affects the site selection of the MC algorithm. Everything else,
                                               ! like the setup of the total simulation cell and the dependendy marking, takes place
                                               ! in the total cell.
   real(dp) :: sLV(3,3)                        ! super lattice vectors
   integer  :: L(3,3), D(3,3), D4(4)
   real(dp) :: invA(3,3)

   ! Pizza
   type(figRep), pointer :: pizza(:,:) => null() ! This is the list of symmetrically-equivalent figures that touch each site
   type(figRep), pointer :: cv_pizza(:,:) => null() ! This is the list of symmetrically-equivalent figures that touch each site
   
   ! References
   type(MC_variables), pointer :: reference(:) ! { reference # }
   type(t_MCg2gref)  , pointer :: g2gref(:)    ! { reference # }
   real(dp)                    :: scaleEnergy

   ! MC variables 
   character(30) :: cboxSize                  ! character equivalent of boxSize (for printing)
   real(dp), pointer :: x(:,:)                ! { species#, coupling# } Concentration of each atom type in the simulation
   real(dp), pointer :: dx(:,:)               ! { species#, dvector# }  dvector dependent concentrations

   integer           :: nMuTemperatures, nMuX
   real(dp), pointer :: muTemperatures(:)         ! { temperature# }: the temperatures where a mu is supplied
   real(dp), pointer :: muX(:,:)                  ! { rank#, conc# }: for each rank, the concentrations where a mu is supplied
   real(dp), pointer :: mu(:,:,:)                 ! { atom spin (!), temperature#, conc# } chemical potential for each atom type
                                       !E!!!!!! NOTE: have to generalise the concentration and mu for coupled CEs
   real(dp), pointer :: Dmu(:,:)                ! { atom spin (!), atom spin (!) } Difference in chemical potential (lookup-table)
   real(dp), pointer :: currMu(:)                 ! current chemical potential
   integer(si), pointer :: possibleNewSpins(:,:)  ! { index, current spin }: list (1st index)  of possible new spin
                                                  ! occupations if current spin is given (2nd index)
                                                  ! e.g.,  possibleNewSpins(:,-1) = (/ 0,1 /) in a ternary case

   real(dp) :: eps              ! epsilon for finite precision equality checking

   integer :: maxgFigLen        ! maximal length of an interaction, measured in gspace (used for parallel MC
                                ! dependency markings)

   integer              :: nTemperatureSchedule
   integer(si), pointer :: TemperatureMode(:) ! 0 = ramp for T,  i.e., use temperature decrease factor
                                              ! 1 = Delta for T, i.e., go down linearly
   real(dp), pointer    :: changeT(:)
   real(dp), pointer    :: changeScheduleT(:)
   real(dp), pointer    :: startT(:), endT(:) ! Start and end temperature for MC-Run

   real(dp) :: E_conv           ! Energy criterion for E(T) convergence
   integer(si) :: nDepLists     ! number of dependency lists
   integer     :: depBoxLength(3)  ! coarse graining size of the dependency marking

   logical     :: onlyTestSpeed ! do we only want to test the speed (and then abort)?
   real(dp)    :: SIstepsFactor

   integer(li)  :: maxnSteps ! maximal number of MC steps
   real(dp), pointer :: norm_vec(:,:)  ! (3xn) normal vectors for MCcell evaluation
   real(dp), pointer :: basis_vec(:,:) ! (3xn) basis vectors for MCcell evaluation
   integer  :: planes ! how many planes should be evaluated from MCcell?
   integer  :: rank ! Number of atom types (binary, ternary, etc.)
   integer  :: ncouplings 

   type(t_list), pointer :: variableDVectorsOfCoupling(:)   ! { coupling# } contains a list with
                                             ! the indices of the dvectors that can have variable Occupation.
                                             ! = all dvectors \ CE%fixedOcc
   type(t_list), pointer :: occupationEquivalenciesOfDVector(:) ! { dvector# } contains a list with
                                             ! the indices of dvectors that have the same occupation as dvector#
   
   integer(li)           :: nAt            ! number of atoms in the sim. cell, prod(d)
   integer(li), pointer  :: nAtVariable(:) ! { coupling# } number of atoms in the sim. cell that can actually change their spin
   integer, pointer, dimension(:)  :: rk, sp ! Rank space and spin space mapping (-k/2..k/2, 1..k)
!   integer :: nBasIneq
!   integer,  pointer :: BasEq(:)      ! lookup table showing which basis corresponds to which onsite
   integer(li), pointer :: onsiteNorm(:) ! Number of occurrences of each type of atom (for onsite normalization)

   integer(li) :: avsteps   ! over how many MC-steps will be averaged?
   integer :: nTsteps_write_MCcell ! after how many temperature steps shall a MCcell file be written?
   character(30) :: fname

   ! continue parameters
   logical                :: continue                 ! Continue a previous MC-run?
   character(1)           :: continue_fileOrStruc     ! continue with a given structure ('S') or from a file ('F')
   character(30)          :: continue_FileName        ! MCcell-file to read in 
   type(crystal), pointer :: continue_structure(:)    ! structure that serves as a starting point for the MC simulation
                                                      ! (must be pointer(:) in order to use check_input_structure routine)
   integer                :: continue_TemperatureStep ! stepnumber at which calculation will be continued
   real(dp)               :: continue_Temperature     ! Temperature at which continued MC-cell is read in
   real(dp)               :: continue_Energy
   logical                :: continue_EnergyConverged
   integer                :: continue_iTemperatureSchedule
   integer                :: continue_scale(3)        ! Scale factor of the MC continue run (e.g. 2 2 2 => MC cell will be
                                                      ! twice the sice of continue cell in each direction


   logical :: GC        ! GrandCanonical Ensemble or Canonical?     ! should be redundant, see mode
   logical :: pGC       ! partial GrandCanonical                    ! should be redundant, see mode
   integer(si) :: mode  ! 1 = gc, 2 = c, 3 = pgc
   integer(si), pointer :: pGCspecies(:)                            ! should be redundant, see GCspins
   integer(si), pointer :: GCspins(:)      ! spins for which gcMC runs, for the other spins a cMC runs

endtype MC_variables


type structure_prediction_t
character(1)              :: coordsys       ! coordinate system used for real space output
logical                   :: compareStrucs  ! should the input structures be compared with the predictions?
real(dp), pointer         :: xi(:,:), xf(:,:)  ! concentration ranges: xi = start, xf = stop
logical                   :: automatic      ! automatic search?
integer                   :: nSPX           ! number of structures per concentration (for automatic search)
integer                   :: nS             ! number of structures (for manual search)
integer, pointer          :: struclist(:)   ! the list of enum-numbers of the structures
real(dp), pointer         :: energylist(:)  ! the energy (DHf) of the structures
type(crystal), pointer    :: strucs(:)      ! the structures itself (will be built later than struclist and
                                            ! energylist, using the information in those lists.)
integer, pointer          :: GSList(:)      ! the list of groundstates
end type structure_prediction_t


! TK @ GH: can we remove this again???
!
!!
!!CONTAINS
!!! Don't call this routine, that is, don't deallocate things that 
!!! have been copied *from*. The copied things just contain pointers to the original elements. So if
!!  ! you deallocate the elements, those copied structures will become corrupted.
!!subroutine deallocate_figure_array(fig)
!!type(figure), pointer :: fig(:)
!!integer i
!!do i = 1, size(fig)
!!   print *,"deallocating fig #",i
!!   call deallocate_figure_elements(fig(i))
!!enddo
!!deallocate(fig)
!!end subroutine deallocate_figure_array
!!subroutine deallocate_figure_elements(fig)
!!type(figure) fig
!!if(associated(fig%vertex))deallocate(fig%vertex)
!!if(associated(fig%label))deallocate(fig%label)
!!if(associated(fig%gPt))deallocate(fig%gPt)
!!if(associated(fig%S))deallocate(fig%S)
!!end subroutine deallocate_figure_elements


END MODULE
