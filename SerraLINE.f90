!						22/07/2019
! SerraLINE: MAIN
! Calculates bending angles. 

  program SerraLINE

  use IO_mod
  use parms
  use functions_mod

  implicit none


! PARAMETERS
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! circle_ str -> Logical to indicate if structured being analysed is closed or opened
  logical :: circle_str

! nbp : number of bp
! box : If there's a box in topology file
! str : parameter needed for reading trajectory. Tell if structure to be analysed is
!       double or single stranded, but for SerraLINE, this parameter is always single stranded
! frames : number of frames
! ierror : for checking errors
! ndim : used for obtaining mdim
! mdim : number of parameters calculated for opened structure and number of parameters per length
!        for closed structures (= nbp)
! max_width  : index of maximum width calculated (c_ is for circular structures)
! max_height : index of maximum height calculated (c_ is for circular structures)
  integer :: nbp, box, str, frames, ierror, &
           & ndim, mdim, & !mdim is the extended dimension
           & i, l

! seq : sequence
  character(1), allocatable :: seq(:)

! top : topology file
! traj : trajectory file
  character(360) :: top, traj

! coords : coordinates
! tangents: tangent vectors
! bends : bending angles
! avstd_ : are average (,1) and standard deviations (,2)
! c_ : indicate parameter for closed strucutres.
  real(dp), allocatable :: coords(:,:,:), &
          & tangents(:,:,:), bends(:,:), & 
          & avstd_bends(:,:), &             
          & c_bends(:,:,:), c_avstd_bends(:,:,:)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------


  !INPUT SECTION
  !-----------------------------------------------------------------------------------
  !Get the directories of trajectory and topology, and 
  !read if its a closed structure or not. 
  call SerraLINE_inputs(traj,top,circle_str,str,nbp)

  !READING SECTION
  !-----------------------------------------------------------------------------------
  !If we have a topology file
  if (str /= 0) then 
    write(6,*) "Reading topology file"
    call topology_amber(nbp,box,seq,top,str)

  !If not, the whole sequence will be " " (Nothing)
  else
    allocate(seq(nbp), stat=ierror)
    if (ierror /= 0) stop "Error in allocating sequence"
    seq = " " 
  end if

  !We do not care about ring atoms, we suppose that the topology file only 
  !contains the atoms needed to draw the curvature

  !Reads all the coordinates from the amber trajectory file
  !It assumes that the number of bp is the same number of atoms
  ! e.g. there's not more atoms in the trajectory/topology
  write(6,*) "Reading trajectory file"
  call coordinates_amber_WrLINE(coords,nbp,frames,traj,box,str)

  !-----------------------------------------------------------------------------------

  !DIMENSION SECTION
  !-----------------------------------------------------------------------------------
 
  !Get dimensions for opened or closed structures 

  if (circle_str) then
    ! closed structure
    mdim = nbp
  else         
    ! opened structure 
    ndim = nbp-1
    mdim = ndim*(ndim-1)/2
  end if

  !-----------------------------------------------------------------------------------

  !TANGENT VECTORS SECTION
  !-----------------------------------------------------------------------------------
 
  !Get tangent vectors to the curve (rj-ri)

  !If the structure is not a circle then the last bp does not has tangent vector
  write(6,*) "Calculating tangent vectors"
  if (circle_str) then
    allocate(tangents(3,nbp,frames), stat=ierror)
  else
    allocate(tangents(3,nbp-1,frames), stat=ierror)
  end if
  if(ierror/=0) stop "Error in allocating tangent vectors"
 
  call get_tangent_vectors(coords,nbp,frames,circle_str,tangents)
 
  !Coordinates are not longer needed
  deallocate(coords, stat=ierror)
  if(ierror/=0) stop "Error in deallocating coordinates"

  !-----------------------------------------------------------------------------------
 
  !BEDINGS SECTION 
  !-----------------------------------------------------------------------------------

  write(6,*) "Calculating bending angles"

  if (circle_str) then
    ! closed structure

    allocate(c_bends(frames,mdim,mdim-1), c_avstd_bends(2,mdim,mdim-1), stat=ierror)
    if(ierror/=0) stop "Error in allocating bendings"

    !Calculate bending angles for closed structure
    call get_c_bendings(mdim,frames,tangents,c_bends)    

  else         
    ! opened structure 

    allocate(bends(frames,mdim), avstd_bends(2,mdim), stat=ierror)
    if(ierror/=0) stop "Error in allocating bendings"

    !Calculate bending angles for opened structure
    call get_bendings(ndim,frames,tangents,bends)    

  end if

  deallocate(tangents, stat=ierror) !No longer needed
  if (ierror/=0) stop "Error in deallocating tangent vectors"

 
  !Averages, standard deviations SECTION 
  !-----------------------------------------------------------------------------------

  !Now calculate averages and standard deviations using function
  !average_std()
  write(6,*) "Calculating averages and standard deviations"

  if (circle_str) then
    !Closed structure
    
    !bends
    do l=1,mdim-1
      do i=1,mdim
        c_avstd_bends(:,i,l) = average_std(c_bends(1:frames,i,l),frames)
      end do
    end do

    deallocate(c_bends, stat=ierror) !No longer needed
    if (ierror/=0) stop "Error in deallocating bendings"

  else    

    !Opened structure

    !bends
    do i=1,mdim
      avstd_bends(:,i) = average_std(bends(1:frames,i),frames)
    end do

    deallocate(bends, stat=ierror) !No longer needed
    if (ierror/=0) stop "Error in deallocating bendings"

  end if   !circle_str
  !-----------------------------------------------------------------------------------

  !WRITING SECTION
  !-----------------------------------------------------------------------------------
  if (circle_str) then

    !Write bendings
    call write_c_av_bends(nbp,frames,seq,c_avstd_bends)

  else
  
    !Write bendings
    call write_av_bends(nbp,ndim,mdim,frames,seq,avstd_bends)
      
  end if   !circle_str
 
  !CLEAR MEMORY
  !-----------------------------------------------------------------------------------
  if (circle_str) then

    deallocate(c_avstd_bends, stat=ierror)
    if (ierror /=0) stop "Error in deallocating parameters"

  else

    deallocate(avstd_bends, stat=ierror)
    if (ierror /=0) stop "Error in deallocating parameters"

  end if

  !THE END
  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------

  end program SerraLINE
