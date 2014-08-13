module monte_support
  use structure_utilities
  use ce_types
  use vector_matrix_utilities, only : norm
  implicit none

contains
!********************************************************************************
! see also: generate_pizza_figures
!
! This routine is supposed to supersede generate_pizza_figures. Its ansatz is to have the
! star of all equivalent figures attached to all possible adjacent sites of the origin site
! and then to check whether a one of those figures connects to the origin site. This also
! works well with multilattices.
!
! In the first part of the routine, a different structure for the temporary pizza is used. 
! I (tk) regard the data structure of pizza(:,:) as sub-optimal. However, changing pizza(:,:)
! everywhere in the code is too much for the time being. Thus, having constructed the
! temporary pizza figures with a new data structure (which also has some flaws, as I cann see
! now), I resort to the legacy data structure of pizza(:,:) for constructing the pizza figures
! out of the temporary pizza figures (the temporary ones  might not be unique).
!
! --> IN
!            eqvFig   :   the star of equivalent figures
!                MC   :   MC data, essentially used for the onsite norms and
!                         for determining the number of the adjacent sites that might
!                         connect to the origin site: maxgFigLen
!                CE   :   CE data
!                 D   :   Smith normal form of the MC cell
! <-- OUT
!              pizza  :   the pizza figures. Note that the whole algorithm proceeds in 
!                         g-space, so no more real space coordinates are necessary.
!
subroutine generate_multilattice_pizza(eqvClusters, MC, CE, D, pizza, MPIr)
type(figRep), intent(in)            :: eqvClusters(0:)
type(MC_variables), intent(in)      :: MC
type(CE_variables), intent(in)      :: CE
integer, intent(in), dimension(3,3) :: D ! SNF
type(figRep), pointer               :: pizza(:,:)  ! { d-vector label, cluster# } intent(OUT), legacy data structure
integer, intent(in)                 :: MPIr

type(t_pizza), pointer         :: tpizza(:)   ! { d-vector label } temporary storage
type(figure)                   :: rep         ! current figure representative
type(figure)                   :: trep        ! temporary figure representative
integer, pointer               :: origFigi(:,:) ! original figure index
integer, pointer               :: tFigi(:,:)

integer iLV1, iLV2, iLV3, nCl, iCl, iRep, nRep, iV, nV, iD, jD, nD, cPF, ioD, ioF, i, tnCl, tnoF, loc(1), jF
integer maxnCl
integer g4(4), g3(3)

logical addToPizza
logical, pointer :: unique(:,:)
integer, pointer :: uniquetpizze(:)

if (MPIr==0) write(*,'(A)') "Generating MC clusters." !. Multilattice extension."

nCl = size(eqvClusters)
nD  = CE%nBas


! allocate the temporary pizza storage
allocate(tpizza(CE%nBas)) ! there's a pizza for every dvector
maxnCl = (2*MC%maxgFigLen+1)**3*nCl  ! maximal number of clusters that can attach from somewhere else to my site

do iD=1,nD; allocate(tpizza(iD)%cluster(maxnCl)); tpizza(iD)%nCl=0; enddo
allocate(origFigi(iD,maxnCl))


! Failsafe: check if all representatives of the eqvClusters start in g-coordinate 0/0/0
do iCl=0,nCl-1 ! Loop over all figures
  nRep = size(eqvClusters(iCl)%Rep)  
  do iRep=1,nRep ! Loop over equivalent representations of that figure
    nV  = eqvClusters(iCl)%Rep(iRep)%nV
    do iV=1,nV
      g3 = eqvClusters(iCl)%Rep(iRep)%gPt(1:3,iV)
      if (all(g3==0)) exit ! g3 is in the 0 cell
      if (iV==nV) stop "ERROR: bad programming in generate_pizza_figures (code: 308jd)"
    enddo
  enddo
enddo

do iD=1,nD

   !write(*,'(/,"d-vector ", I5,/)') iD

   do iCl=0,nCl-1 ! Loop over all clusters (ECIs!)
  
     !write(*,'(I5,1x)',advance='no') iCl

     nRep = size(eqvClusters(iCl)%Rep)  
     do iRep=1,nRep ! Loop over equivalent representations of that figure
       if (.not. any(eqvClusters(iCl)%Rep(iRep)%label==iD)) cycle ! if the figures does not attach to the
                                        ! dvector that we're looking at, we can skip the rest
       rep = eqvClusters(iCl)%Rep(iRep)
       nV  = eqvClusters(iCl)%Rep(iRep)%nV
!       if (associated(trep%gPt)) deallocate(trep%gPt);
       if (associated(trep%gPt)) trep%gPt => null()
       allocate(trep%gPt(4,nV))
   
       ! In the following, we rely on the presence of the g
       if (.not. associated(rep%gPt)) stop "ERROR: Bad programming in generate_pizza_figures"

       ! Brute force loops. Something better anybody?
       do iLV1=-MC%maxgFigLen,MC%maxgFigLen
         do iLV2=-MC%maxgFigLen,MC%maxgFigLen
           do iLV3=-MC%maxgFigLen,MC%maxgFigLen
             ! take the g-coordinates of the current representative and shift them to a new origin
             ! (leave the 4th coordinate as is)
             do iV=1,nV
               trep%gPt(:,iV) = rep%gPt(:,iV) + (/ iLV1,iLV2,iLV3, 0 /)  
               call mod_g3(trep%gPt(:3,iV),D,trep%gpt(:3,iV)) ! bring coordinates back into 1st supercell
             enddo
! Don't sort the vertices, since if you do you also have to take care of the s vector!!
!             call sort_gPt_vertices(trep%gPt)
             
             ! the temporary representative trep is now shifted and (w.r.t. the original figure
             ! representative as given in figures.out) in general rotated, since it stems from the
             ! "star" of figures that constitute the eqvClusters.
   
             ! Check now whether this representative touches the original site. If yes, build a pizza with it.
             addToPizza = .false.
             do iV=1,nV
               if (all(trep%gpt(1:3,iV)==0 .and. trep%gpt(4,iV)==iD)) then ! the temp represenentative touches the
                                                                          ! origin box and also the respective dvector
                 addToPizza = .true.
                 exit ! exit the loop
               endif
             enddo
             if (addToPizza) then
               tpizza(iD)%nCl = tpizza(iD)%nCl + 1
               cPF            = tpizza(iD)%nCl        ! counter for pizza figures
               allocate(tpizza(iD)%cluster(cPF)%gPt(4,nV))
!               allocate(tpizza(iD)%cluster(cPF)%s( size(rep%s,1), size(rep%s,2) )) ! rep, sic!
               allocate(tpizza(iD)%cluster(cPF)%s( nV ))
               tpizza(iD)%cluster(cPF)%gPt   = trep%gPt
               tpizza(iD)%cluster(cPF)%s    =  rep%s ! rep, sic!
!               tpizza(iD)%cluster(cPF)%norm  = eqvClusters(iCl)%count *D(1,1)*D(2,2)*D(3,3)
               origFigi(iD,cPF) = iCl  ! store original figure index
             endif

           enddo ! shift representative by iLV3
         enddo ! shift representative by iLV2
       enddo ! shift representative by iLV11
   
     enddo ! representatives
   enddo ! figures

enddo ! iD

! here we have all the pizza figures, stored in tpizza(iD) for each dvector (iD). However,
! lots of them will be EQUAL (not only equivalent, but really equal). Have to sort those
! out in the following.

allocate(unique(nD,maxnCl))   ;   unique=.true.
do iD=1,nD
  tnCl = tpizza(iD)%nCl
  do iCl=1,tnCl ! loop over figures
    if (.not. unique(iD,iCl)) cycle
    do ioF=1,tnCl ! loop over the _o_ther figures

      if ( .not. ( iCl==ioF ) ) then
      if ( size(tpizza(iD)%cluster(iCl)%gPt,2) == size(tpizza(iD)%cluster(ioF)%gPt,2) ) then  ! does the number of vertices coincide?
      if ( all(tpizza(iD)%cluster(iCl)%gPt == tpizza(iD)%cluster(ioF)%gPt) ) then ! if yes: are all the g-coordinates equal?
      if ( all(tpizza(iD)%cluster(iCl)%s  == tpizza(iD)%cluster(ioF)%s) ) then  ! if yes: are all the s vectors equal?
      if ( origFigi(iD,iCl) == origFigi(iD,ioF) ) then ! if yes: do they hail from the same original cluster? (usually, they do, 
                                                       ! but if you happen to have two *identical* clusters in your J.out file    
                                                       ! with different ECI values, you do not want to remove them)               
                                                       ! (How can that happen anyway? Manually!)                                  
        unique(iD,ioF) = .false.  ! if yes: ==> this pizza is not unique.
      endif
      endif
      endif
      endif
      endif
    enddo

  enddo ! iCl
enddo ! iD


! tk: I don't really like the old pizza(:,:) structure. That's why I setup a different structure above.
! However, at the moment I'm not going to change everything pizza-related in the code. So the following
! is a legacy code, converting the above structure to the old pizza(:,:) structure.
! (I know, it looks rather overcomplicated...)
allocate(pizza(nD,0:nCl-1))
do iD=1,nD

  tnCl = tpizza(iD)%nCl   ! (temporary) number of pizza figures that attach to iD

  ! set the pizza for cluster #0
  allocate( pizza(iD,0)%Rep(1)         )
  allocate( pizza(iD,0)%Rep(1)%s(0)   )
  allocate( pizza(iD,0)%Rep(1)%label(0))
  allocate( pizza(iD,0)%Rep(1)%gPt(4,0))
  pizza(iD,0)%norm = 1

  ! set pizza for all other clusters
  do iCl=1,nCl-1

    ! determine some general feature of figure iCl:
    nV = eqvClusters(iCl)%Rep(1)%nV            ! number of vertices of that figure

    ! number of unique figures that attach to iD and come from figure (ECI) iCl: this is the
    ! number of representatives of this pizza
    ! original code:
    !    where (.not. unique(iD,1:tnCl) )
    !      origFigi(iD,:) = -1 ! we do not want to look at figures that are not unique. Set the corresponding
    !                          ! figure index to -1 to mark it as useless.
    !    end where
    ! new code: (if you have tons of figures, where breaks...)
    do jF=1,tnCl; if (.not. unique(iD,jF)) origFigi(iD,jF)=-1; enddo

    nRep = count(origFigi(iD,1:tnCl) == iCl) ! count those (unique ones) that come from figure iCl (the -1 marked
                                             ! are obiously not counted since there is not figure -1)

    ! if THIS figure (iCl) does not attach to iD (can happen!), then nRep==0. Fortunately,
    ! this doesn't embarass allocate.
    allocate(pizza(iD,iCl)%Rep(nRep))
    pizza(iD,iCl)%count = eqvClusters(iCl)%count ! eqvClusters, sic! BEWARE: this means, for the pizza, %count is NOT the number of the pizza's representatives!
    pizza(iD,iCl)%norm  = eqvClusters(iCl)%count *D(1,1)*D(2,2)*D(3,3)
    !print *, "pizza(iD,iCl)%norm",icl,pizza(id,icl)%count,D(1,1)*D(2,2)*D(3,3),pizza(id,icl)%norm

    ! now fill all the representatives with data
    do iRep=1, nRep

      pizza(iD,iCl)%Rep(iRep)%vertex => null()  ! do not want to allocate vertex. Those are real space
                                               ! coordinates that we don't want to use from here on
      allocate(pizza(iD,iCl)%Rep(iRep)%s(nV))
      allocate(pizza(iD,iCl)%Rep(iRep)%label(nV))   ! don't know why we need this, its already included in gPt
      allocate(pizza(iD,iCl)%Rep(iRep)%gPt(4,nV))

      ! number of vertices
      pizza(iD,iCl)%Rep(iRep)%nV = nV

      ! now: get the next tpizza(iD)%cluster(???) item that belongs to figure iCl. I use maxloc to 
      ! get the first occurence of origFigi==iCl inside origFigi. Since I don't want to find this
      ! position next time I'm here, I mark it by -2. (could also use -1, but for debug reasons like
      ! to keep things separate from the marking with -1 above)
      ! * find
      loc = maxloc(origFigi(iD,1:tnCl), origFigi(iD,1:tnCl)==iCl)     ! I love fortran :-)
      ! * mark
      origFigi(iD,loc(1)) = -2
      ! * set pizza accordingly
      pizza(iD,iCl)%Rep(iRep)%gPt = tpizza(iD)%cluster( loc(1) )%gPt
      pizza(iD,iCl)%Rep(iRep)%s  = tpizza(iD)%cluster( loc(1) )%s
      pizza(iD,iCl)%Rep(iRep)%label(:) = pizza(iD,iCl)%Rep(iRep)%gPt(4,:)  ! set the label

    enddo
  enddo ! iCl
enddo ! iD


! deallocate all tpizza things                                                                                                    
!--don't do that--!>do iD=1,size(tpizza)                                                                                          
!--don't do that--!>   do iCl=1,size(tpizza(iD)%cluster)                                                                          
!--don't do that--!>      deallocate(tpizza(iD)%cluster(iCl)%gPt)                                                                 
!--don't do that--!>      deallocate(tpizza(iD)%cluster(iCl)%s)                                                                   
!--don't do that--!>   enddo                                                                                                      
!--don't do that--!>   deallocate(tpizza(iD)%cluster)                                                                             
!--don't do that--!>enddo                                                                                                         
deallocate(tpizza)

!do iD=1,nD
!  do iF=1,nF ! loop over figures
!    if (.not. unique(iD,iF)) cycle
!    print *, "================================================================================"
!    print *, "iD, Figure #",iD,iF
!    do iV=1,size(tpizza(iD)%figure(iF)%gPt,2)
!      print *, tpizza(iD)%figure(iF)%gPt(:,iV)
!    enddo
!  enddo
!enddo

!do iD=1,nD
!  do iF=1,newpizza(iD)%nF
!    print *, "================================================================================"
!    print *, "iD, Figure #",iD,iF
!    do iV=1,size(newpizza(iD)%figure(iF)%gPt,2)
!      print *, newpizza(iD)%figure(iF)%gPt(:,iV)
!    enddo
!  enddo
!enddo


if (any(pizza(:,:)%norm<=0)) then
  print *, "norm=",pizza(:,:)%norm
  stop "ERROR: bad programming, integer overflow!"
endif



end subroutine generate_multilattice_pizza

!................................................................................

! Helper subroutine to set the site variable at some g-point (g4).
! (It looks stupid to have that routine, but take a look at the form of g4 for references, see
! above in the main routine!)
pure subroutine set_spin_at_g(site,g4,spin)
integer(si), pointer       :: site(:,:,:,:)  ! must be a pointer, intent(inout) does not suffice :-(
integer    , intent(in)    :: g4(:)
integer(si), intent(in)    :: spin

site( g4(1),g4(2),g4(3),g4(4) ) = spin
end subroutine set_spin_at_g


!*******************************************************************************
! Subroutine to find the maximum extent of any figure in g-space
! --> IN
!          eqvFig     : equivalent figures. Will be looped over to find the maximum extent
!      L, D, invA     : g-space information
!            pBas     : parent basis vectors
!         boxsize     : the size of the MC simulation cell
!             eps     : finite precision checking for the conversion to g-space
!
! <-- OUT
!      maxgFigLen     : maximum extent of any figure in g-space
!
subroutine get_maxgFigLen(eqvFig,L,D,invA,pBas,boxsize,eps,  maxgFigLen  )
integer, intent(out) :: maxgFigLen           ! return value
integer, pointer     :: maxgFigLenForRep(:)  ! helper variable

type(figRep), intent(inout)          :: eqvFig(:)
integer, dimension(3,3), intent(in)  :: L, D
real(dp), dimension(3,3), intent(in) :: invA
real(dp), pointer                    :: pBas(:,:) ! parent basis vectors (3xnBas)
integer, intent(in)                  :: boxSize(3)
real(dp), intent(in)                 :: eps

real(dp) :: vertex(3)
integer iCl, nCl, iRep, nRep, iV, nV, iPbas
integer :: g4(4), g3(3)
integer :: tmp_gPt(3)


nCl = size(eqvFig)

!________________________________________________________________________________
! Find the maximum extent of a figure in g-space
!
! 1) Translate real space figures into g-space
! 2) Find the maximum
! 3) Clean up
!
! Remark: both 1) and 2) could be combined into 1 single step, but I want
!         it to be more modular, and it is certainly no time issue here
!
! 1) Translate real space figures into g-space
!
do iCl = 1, nCl
  nV   = eqvFig(iCl)%Rep(1)%nV
  nRep = size(eqvFig(iCl)%Rep)
  do iRep=1, nRep
    allocate(eqvFig(iCl)%Rep(iRep)%gPt(4,nV))
    do iV=1, nV
      iPbas = eqvFig(iCl)%Rep(iRep)%label(iV)
      vertex = eqvFig(iCl)%Rep(iRep)%vertex(:,iV)
      call find_g(vertex,invA,L,D,pBas,iPbas,eps,g3)
      g4 = (/g3, iPbas/)
      eqvFig(iCl)%rep(iRep)%gPt(:,iV) = g4
    enddo ! iV
  enddo ! iRep
enddo ! iCl

!
! 2) find the maximum extent of all figures in g-space
!
maxgFigLen = 0
do iCl = 1, nCl
    nRep = size(eqvFig(iCl)%Rep)
    allocate(maxgFigLenForRep(nRep))
    maxgFigLenForRep=0
    do iRep=1, nRep
      do iV= 1, eqvFig(iCl)%Rep(1)%nV
        tmp_gPt=eqvFig(iCl)%Rep(iRep)%gPt(1:3,iV)+boxSize/2 ! unfortunately, that's a cast to real(dp), multilattices??
                                                       ! TK: that should be ok for multilattices, box size is for the
                                                       ! unit cells only, is it not?
        tmp_gPt=mod(tmp_gPt,boxSize)-boxSize/2
!DEB        write(*,'(I5,I5,I5,A,3I5)') iCl,iRep, iV, "...", tmp_gPt
        maxgFigLenForRep(iRep) = max(maxgFigLenForRep(iRep),ceiling(norm(real(tmp_gPt,dp))))
      enddo ! iV
    enddo ! iRep
!    maxgFigLen = max(maxgFigLen,minval(maxgFigLenForRep))
    maxgFigLen = max(maxgFigLen,maxval(maxgFigLenForRep))
    deallocate(maxgFigLenForRep)
enddo

! No longer clean up; keep it.
!>!
!>! 3) Clean up
!>!
!>do iCl = 1, size(eqvFig)
!>  do iRep=1, size(eqvFig(iCl)%Rep)
!>    deallocate(eqvFig(iCl)%Rep(iRep)%gPt)
!>  enddo
!>enddo

!
! Uncomment the following line for debug purposes in order to fix the figure length
! (-> wrong results!)
!
!maxgFigLen=20


end subroutine get_maxgFigLen

end module monte_support
