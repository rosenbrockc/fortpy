!!<summary>Provides functionality to track multiple order parameters in a 
!!thermodynamic monte carlo simultion. Provides access to a single public
!!type cluster_tracker that tracks the correlations for one or more clusters
!!across mulitple iterations of the simulation.</summary>
module monte_orderparam
  use ce_types
  use io_utilities, only: read_clusters, newunit
  use setup_cluster_expansion, only: generate_eqvClusters
  use structure_utilities, only: get_gTab_correlations
  use monte_support
  implicit none
  private

  !!<summary>Implements tracking of correlations for one or more clusters across multiple
  !!iterations of a thermodynamic monte carlo simulation. Clusters are read from a file
  !!named clustertrack.in using @CREF[io_utilities.read_clusters]. The correlations are
  !!averaged over the same number of flips as the parent simulation. Correlations will be
  !!written to a file clustertrack.out as the simulation progresses.</summary>
  !!<usage>Allocate a pointer to the cluster tracker and allocate it before the iterations begin.
  !! - call @CREF[this.enabled] to check whether the tracker can run.
  !! - call @CREF[this.initialize] before calling track. Only needs be called once.
  !! - deallocate(tracker) to clean up the memory and file system before the program exits.
  !! - At each iteration, call @CREF[this.track].</usage>
  type, public :: cluster_tracker
     private
     !!<member name="clusters_file">The name of the file from which the clusters to track will be extracted.
     !!Defaults to "clustertrack.in".</member>
     !!<member name="correlations">The vector of correlations calculated at the current temperature step. The indices
     !!correspond to the relevant cluster in @CREF[this.clusters].</member>
     !!<member name="lastcorrs">The last set of correlations that was calculated at a flip.</member>
     !!<member name="basecorrs">The pre-simulation total correlations for the entire MC cell.</member>
     !!<member name="eqvClusters">The vector of symmetrically equivalent clusters that will be tracked.</member>
     !!<member name="pizza">The pizza of correlations ordered by vertex for computation speed.</member>
     !!<member>The number of clusters being tracked.</member>
     !!<member name="fileunit">I/O unit for the file to write the output to.</member>
     !!<member name="stepcount">The number of flips that actually took place where a correlation value
     !!was calculated and recorded by the tracker.</member>
     character(30), public :: clusters_file = "clustertrack.in"
     real(dp), allocatable :: correlations(:,:)
     real(dp), allocatable :: lastcorrs(:)
     real(dp), pointer :: basecorrs(:) => null()
     type(figRep), pointer :: eqvclusters(:) => null()
     type(figRep), pointer :: pizza(:,:) => null()
     integer :: nclusters
     integer :: fileunit
     integer :: stepcount

     !!<group name="INTERNAL STATUS TRACKING LOGICALS">
     !!<member name="is_enabled">A value indicating whether this simulation has enabled cluster tracking by providing
     !!a clustertrack.in file in the working directory. Use @CREF[this.enabled] to check status.</member>
     !!<member name="file_checked">Indicates whether the instance has already checked the file system for existence
     !!of the clustertrack.in file.</member>
     !!<member name="initialized">Indicates whether the tracker has been initialized and can track correlations.
     !!Public in case someone needs to check or re-initialize.</member>
     !!</group>
     logical, public :: initialized = .FALSE.
     logical :: file_checked = .FALSE.
     logical :: is_enabled = .FALSE.

     contains
       procedure, public :: initialize => tracker_initialize
       procedure, public :: track
       procedure, public :: enabled
       procedure, public :: flip_before
       procedure, public :: flip_after
       procedure, public :: end_track
       final :: tracker_finalize
  end type cluster_tracker

contains
  !!<summary>Initializes the tracker object by allocating any allocatable vars and extracting
  !!the clusters from clustertrack.in, prepares the track to call @CREF[this.track].</summary>
  !!<parameter name="CE" target="stuff">Instance required to compute equivalent clusters. Since the clusters being tracked are a
  !!subset of those in the monte carlo simulation, the lattice definitions stored in CE (and extracted
  !!from the lat.in file) apply well to the tracked clusters.</parameter>
  !!<parameter name="GA">Instance required to call @CREF[io_utilities.read_clusters].
  !!Would be better to modify read_clusters() to not require the GA instance, but for now
  !!this is required for backward compatibility. [QUIRK]</parameter>
  !!<parameter name="eps">When the clusters are imported and sorted by vertices, eps determines how close
  !!two points must be together to be considered equal.</parameter>
  !!<parameter name="MC">Details of the configuration for the overall monte carlo simulation.</parameter>
  subroutine tracker_initialize(this, CE, GA, eps, MC)
    class(cluster_tracker), intent(inout) :: this
    type(CE_variables), intent(in) :: CE
    type(GA_variables), intent(inout) :: GA
    real(dp), intent(in) :: eps
    type(MC_variables), intent(inout) :: MC

    !!<local name="clusters">The vector of clusters to track that were extracted from the @CREF[this.clusters_file].</local>
    !!<local name="maxFigLen">The maximum extent of any of the tracked clusters in g-space. Used to calculate the pizza.</local>
    !!<local name="origFigLen">Because the figure length for the pizza calculation is buried inside the MC object, we need
    !!to overwrite it and then restore it to enable code re-use. This is the original value of the MC figLen.</local>
    integer :: maxFigLen
    type(figure), pointer :: clusters(:)
    integer :: origFigLen
    
    !We only want to initialize the tracker if it is enabled
    if (this%enabled() .and. .not. this%initialized) then
       call read_clusters(clusters, this%clusters_file, GA, eps)
       this%nclusters = size(clusters, 1)

       !The correlation calculation routines require a list of equivalent clusters.
       call generate_eqvClusters(clusters, this%eqvClusters, CE, 0)
       !Before we can generate a pizza, we need to know the extent of the any of the clusters in
       !g-space. This also updates the equivalent clusters with some g-space info that we need for 
       !for the pizza calculation.
       call get_maxgFigLen(this%eqvClusters, MC%L, MC%D, MC%invA, CE%pBas, MC%boxsize, eps, maxFigLen)

       !Now, we can generate a pizza for speeding up the computation. This method needs the maxFigLen
       !for our clusters, not the entire cluster expansion. Since the value for the entire simulation
       !is buried inside MC object, we need to overwrite that value, calculate the pizza and then set

       !it back to its original value before this subroutine exits.
       origFigLen = MC%maxgFigLen
       MC%maxgFigLen = maxFigLen
       call generate_multilattice_pizza(this%eqvClusters, MC, CE, MC%D, this%pizza, 0)
       MC%maxgFigLen = origFigLen

       !Start the output file and add the header string
       open(newunit(this%fileunit), file="clustertrack.out", status="replace")
       write(this%fileunit, *) "# Temperature   Correlations (in the order presented in clustertrack.in)"

       !Allocate the correlations matrix that keeps track of the changes at each flip.
       !We need a +1 on the avsteps because entry 1 is always the total correlations for
       !the entire cell that is calculated at the start of the simulation.
       allocate(this%correlations(MC%avsteps, this%nclusters))
       allocate(this%lastcorrs(this%nclusters))

       !Calculate the correlations for the current temperature step. We only need to do this once.
       call get_gTab_correlations(CE, MC%site, MC%L, MC%D, MC%invA, this%eqvClusters, this%basecorrs, .true.)

       !Set the initialization flag so we don't do this again
       this%lastcorrs = this%basecorrs
       this%stepcount = 1
       this%initialized = .TRUE.
    else
       print *, "WARNING: tried to initialize a tracker without a clustertrack.in file, or the tracker was already initialized."
    end if
  end subroutine tracker_initialize

  !!<summary>Returns a value indicating whether tracking is enabled for this simulation.</summary>
  logical function enabled(this)
    class(cluster_tracker), intent(inout) :: this
    integer :: fileio, ioerr

    !If we have already checked for the existence of the file, just return the enabled value
    if (this%file_checked) then
       enabled = this%is_enabled
    else       
       open(newunit(fileio), file=this%clusters_file, status='old', iostat=ioerr)
       !If there wasn't an error, the file must exist, so the tracker has been enabled
       if (ioerr .eq. 0) then
          this%is_enabled = .TRUE.
          close(fileio)
          enabled = .TRUE.
       end if
       !else, the default value is already false: the tracker is disabled.

       !Register that we have run this check once already so that we don't slow
       !the program down with repeated checks
       this%file_checked = .TRUE.
    end if
  end function enabled

  !!<summary>Calculates the change in correlation *before* the swap. Does not affect
  !!the MC simulation variables. Adds the change to the current total correlation from the
  !!@CREF[this.track] method.</summary>
  !!<parameter name="g1">The g-vector of the first site involved in the spin swap.</parameter>
  !!<parameter name="g2">The g-vector of the second site involved in the spin swap.</parameter>
  !!<parameter name="CE">Cluster expansion details required by the correlation calculator. Should be the
  !!same CE variables used in all correlation calculations for the simulation (not just tracking).</parameter>
  !!<parameter name="MC">Details of the configuration for the overall monte carlo simulation.</parameter>
  subroutine flip_before(this, g1, g2, CE, MC)    
    class(cluster_tracker), intent(inout) :: this
    integer, intent(in) :: g1(4), g2(4)
    type(CE_variables), intent(in) :: CE
    type(MC_variables), intent(inout) :: MC

    !!<local name="pizzaPI1">The pizza correlations for the state of the site.</local>
    real(dp), pointer :: pizzaPI1(:), pizzaPI2(:)

    !Calculate the difference for each of the correlations involved in the flip
    call single_corr(this, g1, CE, MC, pizzaPI1)
    call single_corr(this, g2, CE, MC, pizzaPI2)

    !Update the value of the total correlation for the cell and these clusters.
    this%correlations(this%stepcount, :) = this%correlations(this%stepcount,:) - pizzaPI1 - pizzaPI2
  end subroutine flip_before

  !!<summary>Calculates the change in correlation *after* the swap. Does not affect
  !!the MC simulation variables. Adds the change to the current total correlation from the
  !!@CREF[this.track] method.</summary>
  !!<parameter name="g1">The g-vector of the first site involved in the spin swap.</parameter>
  !!<parameter name="g2">The g-vector of the second site involved in the spin swap.</parameter>
  !!<parameter name="CE">Cluster expansion details required by the correlation calculator. Should be the
  !!same CE variables used in all correlation calculations for the simulation (not just tracking).</parameter>
  !!<parameter name="MC">Details of the configuration for the overall monte carlo simulation.</parameter>
  subroutine flip_after(this, g1, g2, CE, MC)    
    class(cluster_tracker), intent(inout) :: this
    integer, intent(in) :: g1(4), g2(4)
    type(CE_variables), intent(in) :: CE
    type(MC_variables), intent(inout) :: MC

    !!<local name="pizzaPI1, pizzaPI2">The pizza correlations for the state of the site.</local>
    real(dp), pointer :: pizzaPI1(:), pizzaPI2(:)

    !Calculate the difference for each of the correlations involved in the flip
    call single_corr(this, g1, CE, MC, pizzaPI1)
    call single_corr(this, g2, CE, MC, pizzaPI2)
    
    !Update the value of the total correlation for the cell and these clusters.
    this%correlations(this%stepcount, :) = this%correlations(this%stepcount,:) + pizzaPI1 + pizzaPI2

    !Store the value of these correlations that we just calculated to use for the next flip.
    this%correlations(this%stepcount,:) = this%correlations(this%stepcount,:) + this%lastcorrs
    this%lastcorrs = this%correlations(this%stepcount, :)

    !Increment the number of flips that we have calculated correlations for.
    this%stepcount = this%stepcount + 1
  end subroutine flip_after
  
  !!<summary>Calculates the correlation at the specified site.</summary>
  !!<parameter name="g">The g-vector of the site involved in the spin swap.</parameter>
  !!<parameter name="CE">Cluster expansion details required by the correlation calculator. Should be the
  !!same CE variables used in all correlation calculations for the simulation (not just tracking).</parameter>
  !!<parameter name="MC">Details of the configuration for the overall monte carlo simulation.</parameter>
  !!<parameter name="pizzaPI">The pizza correlations that are calculated for the state of the site.</parameter>
  subroutine single_corr(this, g, CE, MC, pizzaPI)   
    class(cluster_tracker), intent(inout) :: this
    integer, intent(in) :: g(4)
    type(CE_variables), intent(in) :: CE
    type(MC_variables), intent(inout) :: MC
    real(dp), pointer :: pizzaPI(:)
    
    !!<local name="L_">When a specific g-vector is specified for correlation calculations, the L
    !!and A matrices are not required in the calculation. These are dummmy variables.</local>
    integer :: L_(3,3)
    real(dp):: invA_(3,3)

    !We need to calculate the correlations. This may be before or after the swap.
    call get_gTab_correlations(CE, MC%site, L_, MC%D, invA_, this%pizza(g(4), :), pizzaPI, mustGetG=.false., atG3=g(1:3))
  end subroutine single_corr

  subroutine new(a, b)
    type(CE_variables), intent(inout) :: Re
  end subroutine

  !!<summary>Calculates the correlations for the clusters being tracked at the current simulation
  !!step. Stores the values in the correlations array and appends the values to the output file if specified.</summary>
  subroutine track(this)
    class(cluster_tracker), intent(inout) :: this

    !Re-initialize the value of the correlations
    this%correlations = 0

    !We can only track if the tracker is enabled and initialized.
    if (this%is_enabled .and. this%initialized) then
       this%stepcount = 2
    else
       print *, "WARNING: attempt to track correlations without initializing the tracker object. Call tracker.initialize() first."
    end if
  end subroutine track

  !!<summary>Saves the averaged correlation values over the entire flip simulation to file.</summary>
  !!<parameter name="temperature">The temperature of the current monte carlo step.</parameter>
  subroutine end_track(this, temperature)
    class(cluster_tracker), intent(inout) :: this   
    real(dp), intent(in)           :: temperature

    !----------------------------------------------------------------------
    !LOCALS
    !----------------------------------------------------------------------
    !real(dp) :: averaged(this%nclusters)
    integer :: i, writeargs
    !----------------------------------------------------------------------

    !Record theses correlations, the next time track() is called, they will be overwritten.
    !We subtract 1 from stepcount because it is incremented at the end of each flip
    !do i = 1, this%nclusters
    !   averaged(i) = SUM(this%correlations(:,i)) / (this%stepcount - 1)
    !end do    

    !Write the simulation averaged correlations for this step to file.
    writeargs = 1 + this%nclusters
    write(this%fileunit, '(<writeargs>f12.7)') temperature, this%lastcorrs !averaged   
    !print *, "Tracked Correlations:", this%lastcorrs
    flush(this%fileunit)

    !Reset all the values to zero for the next temperature swap
    this%correlations = 0
    this%stepcount = 1
  end subroutine

  !!<summary>Releases any allocated pointers and arrays.</summary>
  subroutine tracker_finalize(this)
    type(cluster_tracker), intent(inout) :: this  
    
    if (allocated(this%correlations)) then
       deallocate(this%correlations)
       deallocate(this%lastcorrs)
       deallocate(this%basecorrs)
    end if
    
    if (this%fileunit /= 0) then
       !Close the file that we had open for recording the correlations.
       close(this%fileunit)
    end if
  end subroutine
end module monte_orderparam
