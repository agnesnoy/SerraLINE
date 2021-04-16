!						16/04/2021
! SerraLINE: MAIN
! Calculates bending angles and can also fit a best plane to the structure or
! to some specific atoms. If a fitting is performed, then the program also 
! outpus height, width and aspect ratio (width/height).

!    -----------------------------------------------------------------------
!    Written by Victor Velasco
!
!    SerraLINE is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SerraLINE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!    -----------------------------------------------------------------------

  program SerraLINE

  use IO_mod
  use parms
  use functions_mod

  implicit none


! PARAMETERS
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! circle_ str -> Logical to indicate if structured being analysed is closed or opened
! fitplane -> If true, then fitting wil be performed and if not, real bendings will be
!             calculated
! print_proj -> print projections trajectory if true
! xyz -> if print_proj and xyz are true, then prints the trajectory in xyz format.
!        And if xyz false, it prints the trayectory in crd
  logical :: circle_str, fitplane, print_proj, xyz

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
! max_height : index of maximum height calculatedo (c_ is for circular structures)
! t_length : is the length in which tangent vectors will be defined
  integer :: nbp, box, str, frames, n_fitting, t_length, ierror, &
           & ndim, mdim, &   !mdim is the extended dimension for sizes
           & nldim, mldim, & !dimensions for tangents and bends of linear structure
           & i, l

! best : best plane fitted for k snapshot
! bp_fitting: Indicate which bp to select for fitting the best plane
  integer, allocatable :: best(:), bp_fitting(:), temporal(:)

! seq : sequence
  character(1), allocatable :: seq(:)

! top : topology file
! traj : trajectory file
  character(1200) :: top, traj

! coords : coordinates
! G_n : normal vectors of fitted planes. One for each snapshot k
! tangents: tangent vectors
! bends : bending angles
! width : width between bp i and j in the projected structure
! height : height '' '' '' ....
! aratio: aspect ratio (width/height)
! ** there width is assigned as the first dimension in G_n which is perpendicular
!    to the fitted plane
! avstd_ : are average (,1) and standard deviations (,2)
! c_ : indicate parameter for closed strucutres.
! width, height and aspect_ratio are calculated whne projection method is chosen. 
! They have the  same dimension for eather closed or opened structures. 
! That's why they don't have c_*
! These three quantities have one value per snapshot which are the maximum values (maximum
! width, maximum height and maximum aspect ratio)
  real(dp), allocatable :: coords(:,:,:), G_n(:,:,:), &
          & tangents(:,:,:), bends(:,:), & 
          & avstd_bends(:,:), &             
          & c_bends(:,:,:), c_avstd_bends(:,:,:), & 
          & width(:), height(:), aspect_ratio(:), &
          & avstd_width(:), avstd_height(:), avstd_aspect_ratio(:), &
          & avg_dist_frame(:), max_dist_frame(:), &
          & avg_dist_rel_frame(:), max_dist_rel_frame(:)

! avg_dist : average of average distances between base-pairs and fitted plane
! max_dist : average maximum distance between base-pairs and fitted plane
! both of these have 1-average and 2-standard deviation
! *_frame  : (look up), means the distance per frame.
! *_rel    : means relative to the longest distance in projected structure (Height)
  real(dp) :: avg_dist(2), max_dist(2), &
            & avg_dist_rel(2), max_dist_rel(2)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------


  !INPUT SECTION
  !-----------------------------------------------------------------------------------
  !Get the directories of trajectory and topology, and 
  !read if its a closed structure or not. Also, read if structure to be analysed is
  !closed or linear, and if there are specific atoms to fit.
  !There is an option where no fitting can be performed
  call SerraLINE_inputs(traj,top,circle_str,str,nbp,bp_fitting, &
                      & fitplane,print_proj,xyz,t_length)
  !str=1 !Reading will be performed as if it were a single stranded structure

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
    seq = "X"
  end if

  !We do not care about ring atoms, we suppose that the topology file only 
  !contains the atoms needed to draw the curvature

  !Reads all the coordinates from the amber trajectory file
  !It assumes that the number of bp is the same number of atoms
  ! e.g. there's not more atoms in the trajectory/topology
  write(6,*) "Reading trajectory file"
  call coordinates_amber_WrLINE(coords,nbp,frames,traj,box,str)

  !Small warning for tangent lengths!!!!!!1
  if ( .not. circle_str .and. t_length > nbp-2 ) stop "Tangent length must be higher than N-2 for opened structures"
  if ( circle_str .and. t_length > nbp-1 ) stop "Tangent length must be higher than N-1 for closed structures"

  !If the fitting will be performed, let's clean bp_fitting (in case of bp greater than nbp)
  if (fitplane) then

    if (allocated(bp_fitting) ) then      !If it is not allocated is fine...

      if ( any(bp_fitting(:) > nbp) ) then

        !It is not fine... We need to ignore bad bp
        i = 0        !re-using i and l...
        do l = 1, size(bp_fitting)
          if ( bp_fitting(l) <= nbp) i = i+1
        end do
        if ( i <= 1) stop "Number of bp for fitting plane, must be greater than 1"

        allocate(temporal(i), stat = ierror)
        if (ierror /= 0 ) stop "Error in fixing fitting bps index..."

        !Let's assign values
        i = 0
        do l = 1, size(bp_fitting)
          if ( bp_fitting(l) <= nbp) then 
            i = i+1
            temporal(i) = bp_fitting(l)
          end if
        end do

        deallocate(bp_fitting, stat = ierror)
        if (ierror /= 0 ) stop "Error in fixing fitting bps index..."
        allocate(bp_fitting(i), stat = ierror)
        if (ierror /= 0 ) stop "Error in allocating fitting bps indeces..."
        bp_fitting = temporal
        deallocate(temporal, stat = ierror)
        if (ierror /= 0 ) stop "Error in deallocating temporal array for fixing bp indices for fitting"
        n_fitting = i
      else
        n_fitting = size(bp_fitting) !In case is all fine...
      end if                      !any(bp_fitting(:) > nbp)

    else
      n_fitting = 0 !In case there a not specific fitting bps
    end if                        !allocated(bp_fitting

  end if                          !fitplane
  !-----------------------------------------------------------------------------------

  !DIMENSION SECTION
  !-----------------------------------------------------------------------------------
 
  !Get dimensions for opened or closed structures 

  if (circle_str) then
    ! closed structure
    mdim = nbp
  else         
    ! opened structure 
    ndim = nbp-1          !These first two will be used in sizes since they lose one distance
    mdim = ndim*(ndim-1)/2
    nldim = nbp-t_length  !But in tangents and bends, they lose t_length bps
    mldim = nldim*(nldim-1)/2

  end if

  !-----------------------------------------------------------------------------------

  !PLANE SECTION
  !-----------------------------------------------------------------------------------

  !In case of fitting, let's fit
  if (fitplane) then 

    write(6,*) "Fitting planes"
 
    !A global plane is calculated, then coordinates are projected onto it
    !One plane for snap shot. best(k) is the best fitted plane (index)
    !bp_fitting indicates if base-pairs that the user indicated to be fitted.
    allocate(G_n(3,3,frames), best(frames), stat=ierror)
    if(ierror/=0) stop "Error in allocating tangent vectors"


    !Allocate the distances
    allocate(avg_dist_frame(frames), max_dist_frame(frames), stat=ierror)
    if(ierror/=0) stop "Error in allocating distances"

    call project_coordinates_G_plane(coords, nbp, frames, n_fitting, bp_fitting, G_n, best, &
                                     avg_dist_frame, max_dist_frame)

    !AND PRINT THE PROJECTED TRAJECTORIES
    if (print_proj) then
      call write_projected_coords(coords, G_n, best, seq, nbp, frames, xyz )
    end if

    if (allocated(bp_fitting)) then
      deallocate(bp_fitting, stat=ierror)
      if(ierror/=0) stop "Error in deallocating bp fitting index"
    end if

  end if
  !-----------------------------------------------------------------------------------

  !WIDTH, HEIGHT AND ASPECT RATIO SECTION
  !-----------------------------------------------------------------------------------
 
  !In case of fitting, let's get width and height
  if (fitplane) then 

 
    write(6,*) "Calculating widths, heights and aspect ratios"
 
    if (circle_str) then
      ! closed structure

      allocate(width(frames), avstd_width(2), stat=ierror)
      if(ierror/=0) stop "Error in allocating widths"

      allocate(height(frames), avstd_height(2), stat=ierror)
      if(ierror/=0) stop "Error in allocating heights"

      allocate(aspect_ratio(frames), avstd_aspect_ratio(2), stat=ierror)
      if(ierror/=0) stop "Error in allocating aspect ratios"

      !Calculate width, height and aspect ratios for closed structure
      call get_c_width_height_aratio(mdim,frames,G_n,best,coords,width,height,&
                                   & aspect_ratio)    

    else         
      ! opened structure 

      allocate(width(frames), avstd_width(2), stat=ierror)
      if(ierror/=0) stop "Error in allocating widths"

      allocate(height(frames), avstd_height(2), stat=ierror)
      if(ierror/=0) stop "Error in allocating heights"

      allocate(aspect_ratio(frames), avstd_aspect_ratio(2), stat=ierror)
      if(ierror/=0) stop "Error in allocating aspect ratios"

      !Calculate width and height for opened structure
      call get_width_height_aratio(ndim,mdim,frames,G_n,best,coords,width,height, &
                                   & aspect_ratio)    

    end if !circle_str

  end if   !fitplane

  !-----------------------------------------------------------------------------------
 
  !TANGENT VECTORS SECTION
  !-----------------------------------------------------------------------------------
 
  !Get tangent vectors to the curve (rj-ri)
 
  !Note that this process is the same even if a plane was not fitted

  !If the structure is not a circle then the last t_length bps do not have tangent vectors

  write(6,*) "Calculating tangent vectors"
  if (circle_str) then
    allocate(tangents(3,nbp,frames), stat=ierror)
  else
    allocate(tangents(3,nldim,frames), stat=ierror)
  end if
  if(ierror/=0) stop "Error in allocating tangent vectors"
 
  call get_tangent_vectors(coords,nbp,frames,circle_str,tangents,t_length)
 
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
    if (fitplane) then                                   !If a fitting was performed
      call get_c_bendings(mdim,frames,G_n,best,tangents,c_bends)    
    else
      call get_c_bendings_nofit(mdim,frames,tangents,c_bends)    
    end if

  else         
    ! opened structure 

    allocate(bends(frames,mldim), avstd_bends(2,mldim), stat=ierror)
    if(ierror/=0) stop "Error in allocating bendings"

    !Calculate bending angles for opened structure
    if (fitplane) then                                   !If a fitting was performed
      call get_bendings(nldim,frames,G_n,best,tangents,bends) 
    else
      call get_bendings_nofit(nldim,frames,tangents,bends)    
    end if

  end if

  deallocate(tangents, stat=ierror) !No longer needed
  if (ierror/=0) stop "Error in deallocating tangent vectors"

  if (fitplane) then
    deallocate(G_n, best, stat=ierror) !No longer needed
    if (ierror/=0) stop "Error in deallocating planes"
  end if

  !Relative sizes SECTION 
  !-----------------------------------------------------------------------------------
  !Before calculating averages, I need to get the relative sizes in (%) for each
  !frame. This is done for the projection method
  if (fitplane) then

    !Allocate arrays first
    allocate(avg_dist_rel_frame(frames), max_dist_rel_frame(frames), stat=ierror)
    if(ierror/=0) stop "Error in allocating relative distances"

    !avg_dist first
    do i = 1, frames
       avg_dist_rel_frame(i) = 100.0_dp*(avg_dist_frame(i)/height(i))
    end do !i
    !Now max_dist
    do i = 1, frames
       max_dist_rel_frame(i) = 100.0_dp*(max_dist_frame(i)/height(i))
    end do !i 

    !I separated in loops to make it faster (actually, it won't make a big difference...)
  end if
  !-----------------------------------------------------------------------------------

 
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

    !Parameters that depend on fitting
    if (fitplane) then

      !widths
      avstd_width(:) = average_std(width,frames)

      !heights
      avstd_height(:) = average_std(height,frames)

      !aspect ratios
      avstd_aspect_ratio = average_std(aspect_ratio, frames  )

      deallocate(width, height, aspect_ratio, stat=ierror) !No longer needed
      if (ierror/=0) stop "Error in deallocating widths, heights and aspect ratios"

      !avg_dist, max_dist and relative sizes
      avg_dist = average_std(avg_dist_frame, frames)
      max_dist = average_std(max_dist_frame, frames)
      avg_dist_rel = average_std(avg_dist_rel_frame, frames)
      max_dist_rel = average_std(max_dist_rel_frame, frames)

      deallocate( avg_dist_frame, max_dist_frame, &
                & avg_dist_rel_frame, max_dist_rel_frame, &
                & stat=ierror) !No longer needed as well

    end if  !fitplane

  else      !circle_str

    !Opened structure

    !bends
    do i=1,mldim
      avstd_bends(:,i) = average_std(bends(1:frames,i),frames)
    end do

    deallocate(bends, stat=ierror) !No longer needed
    if (ierror/=0) stop "Error in deallocating bendings"


    !Parameters that depend on fitting
    if (fitplane) then

      !widths
      avstd_width(:) = average_std(width,frames)

      !heights
      avstd_height(:) = average_std(height,frames)

      !aspect ratios
      avstd_aspect_ratio = average_std(aspect_ratio, frames  )

      deallocate(width, height, aspect_ratio, stat=ierror) !No longer needed
      if (ierror/=0) stop "Error in deallocating widths, heights and aspect ratios"

      !avg_dist, max_dist and relative sizes
      avg_dist = average_std(avg_dist_frame, frames)
      max_dist = average_std(max_dist_frame, frames)
      avg_dist_rel = average_std(avg_dist_rel_frame, frames)
      max_dist_rel = average_std(max_dist_rel_frame, frames)

      deallocate( avg_dist_frame, max_dist_frame, &
                & avg_dist_rel_frame, max_dist_rel_frame, &
                & stat=ierror) !No longer needed as well

    end if !fitplane

  end if   !circle_str
  !-----------------------------------------------------------------------------------

  !WRITING SECTION
  !-----------------------------------------------------------------------------------
  if (circle_str) then

    if (fitplane) then

      !Write parameters
      call write_c_av_parms(mdim,frames,seq,c_avstd_bends,avstd_width,avstd_height, &
                          & avstd_aspect_ratio,t_length,avg_dist,max_dist, &
                          & avg_dist_rel, max_dist_rel)
      
    else

      !Write bendings
      call write_c_av_bends(nbp,frames,seq,c_avstd_bends,t_length)

    end if   !fitplane

  else
  
    if (fitplane) then

      !Write parameters
      call write_av_parms(nbp,nldim,mldim,frames,seq,avstd_bends,avstd_width,avstd_height, &
                        & avstd_aspect_ratio,t_length,avg_dist,max_dist, &
                        & avg_dist_rel, max_dist_rel)

    else

      !Write bendings
      call write_av_bends(nbp,nldim,mldim,frames,seq,avstd_bends,t_length)
      
    end if !fitplane

  end if   !circle_str
 
  !CLEAR MEMORY
  !-----------------------------------------------------------------------------------
  if (circle_str) then

    deallocate(c_avstd_bends, stat=ierror)
    if (ierror /=0) stop "Error in deallocating parameters"

    if (fitplane) then
      deallocate(avstd_width, avstd_height, avstd_aspect_ratio, stat=ierror)
      if (ierror /=0) stop "Error in deallocating parameters"
    end if

  else

    deallocate(avstd_bends, stat=ierror)
    if (ierror /=0) stop "Error in deallocating parameters"

    if (fitplane) then
      deallocate(avstd_width, avstd_height, avstd_aspect_ratio, stat=ierror)
      if (ierror /=0) stop "Error in deallocating parameters"
    end if

  end if

  !THE END
  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------

  end program SerraLINE
