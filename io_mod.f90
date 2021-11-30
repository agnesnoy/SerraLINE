!------------------------------------------------------------- 20/04/2021
! This module contains functions and subroutines needed for reading
! and writing.
!
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

  module IO_mod

  use parms
  implicit none

  private arrdims
  public topology_amber, coordinates_amber_WrLINE, SerraLINE_inputs

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  contains

!-----------------------------------------------------------------------
!  Reads topology file from amber topology file
!    Gets number bp (nbp), checks if there's a "box" and the sequence on
!    both strands

  subroutine topology_amber(nbp,box,seq,top,str)
  implicit none
  integer ::  nbp, box, str                !if box == 1 then there is a box
  integer :: i, j, ierror, n_atoms, n_res, iaux(10), nlines, l, f_lines, &
           & pur, pyr  !purines, pyrimidines and type of structure
  character(360) :: top, caux                 !caux helps in reading
  character(1), allocatable, intent(out) :: seq(:)
  character(4), allocatable :: res_names(:)
  character(1), allocatable :: aseq(:) !auxiliar sequence

  !Get dimension of file top. Here, "l" is a dummy variable
  call arrdims(top,f_lines,l)

! READING TOPOLOGY----------------------------------------------------------------
  !Open topology file
  open(unit=10, file=trim(top), status="old",action="read",iostat=ierror)
    if (ierror/=0) stop 'Error in opening topology file'


    n_res = -100  !Just in case
    do j=1,f_lines
      read(10,"(A)") caux

      !Escape statement: if three quarters of file have been readed, then it
      !  it would automatically exit the loop. Hopefully, this will never happend,
      !  hence, it only is an emergency escape before collapsing.
      if (j > 3*f_lines / 4) exit

      !FLAG POINTERS
      if ( trim(caux) == "%FLAG POINTERS" ) then

        do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        read(10, 1001) iaux                                 !Has dimension 10 
        n_atoms=iaux(1)                                     !n_atoms
        read(10, 1001) iaux  
        n_res=iaux(2)                                       !n_res
        read(10, 1001) iaux 
        box=iaux(8)                                         !box

        !If there are no atoms stop.
        if (n_atoms <= 0) stop "No atoms in topology"
        !In case of no residues the program will try to
        !identify them (it could go wrong)

        !Print info
        write(6,"(1A20,1I10)") "Number of atoms = ", n_atoms
        write(6,"(1A20,1I10)") "Number of residues = ", n_res
        write(6,"(1A20,1I10)") "BOX = ", box

      end if
      !End of POINTERS section 

      !FLAG RESIDUE_LABEL
      if ( trim(caux) == "%FLAG RESIDUE_LABEL") then

         do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        !calculate lines in RESIDUE_LABEL section
        if ( mod(n_res,20) == 0 ) then
          nlines = n_res / 20                        !20 residues per line
        else
          nlines = n_res / 20 + 1                 
        end if
        
        !allocate res_names (is a temporary variable)
        allocate(res_names(n_res), stat=ierror)
        if (ierror /= 0 ) stop "Error in allocating res_names"

        !Read residues names
        l=1
        if (nlines > 1) then                         !in case (nlines < 20) which is unlikely
          do i=1,nlines-1
            read(10, 1002) res_names(l:19+l)
            l=l+20
          end do
        end if
        read(10, 1002) res_names(l:n_res)
      end if
      !End of RESIDUE_LABEL SECTION
      
    end do
  close(10, iostat=ierror)
  if (ierror /= 0) stop "Error in closing topology file"
! END OF READING -------------------------------------------------------


  !Check if data obtained from topology file is enough to go on
  if ( n_res < 1) then
    stop "Couldn't read number of residues"
  end if 

! COUNT NBP -------------------------------------------------------
  !The structure will be treated as a double stranded structure or a single strandedi one.
  !In case of dsDNA, if the number of residues in both strands is not the same, then
  !this will be treated as an error. 

  nbp=0
  do i=1,n_res
    !Identify Adenine
    if (any( A_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Guanine
    if (any( G_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Cytosine
    if (any( C_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Thymine
    if (any( T_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Uracil
    if (any( U_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
  end do   !close i
 
       
! GET SEQUENCE -------------------------------------------------------

  !Allocate sequence 
  allocate(seq(nbp), stat=ierror) 
  if(ierror/=0) stop 'Error in allocating sequence'

  !Get the correct number of nbp and check if structure is incomplete in case of 
  !a double stranded structure
  if ( str == 2 ) then
    if ( mod(nbp,2) /= 0 ) stop "Incomplete double stranded structure"
    nbp = nbp/2                 !We counted bases from both strands
  end if  
  !if str=1 then its a single stranded structure. If its 2 then is double stranded

  !Similar loop as before but now identifying residues
  l=0 !will help us count bps
  pur = 0 ; pyr = 0 !purines and pyrimidines counters
  do i=1,n_res
    !Identify Adenine
    if (any( A_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "A" 
      pur = pur +1
    end if
    !Identify Guanine
    if (any( G_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "G" 
      pur = pur +1
    end if
    !Identify Cytosine
    if (any( C_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "C" 
      pyr = pyr +1
    end if
    !Identify Thymine
    if (any( T_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "T" 
      pyr = pyr +1
    end if
    !Identify Uracil
    if (any( U_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "U" 
      pyr = pyr +1
    end if
  end do   !close i
     
  !Check we read all correctly the sequence and indices
  if (l /= nbp*str) stop 'Error in extracting sequence and res pointers'

  !Get sequence of first strand only

  !allcocate auxiliar sequence
  !     note: I was lasy and created a new array because I didn't want to
  !           modify the file, sorry....
  allocate(aseq(nbp*str), stat=ierror)
  if(ierror/=0) stop 'Error in allocating auxiliar sequence'

  aseq = seq
  deallocate(seq, stat=ierror)
  if(ierror/=0) stop 'Error in deallocating sequence'

  !Allocate sequence
  allocate(seq(nbp), stat=ierror)
  if(ierror/=0) stop 'Error in allocating sequence'

  !And then the actual sequence
  do i = 1, nbp
    seq(i) = aseq(i)
  end do

  write(6,*) "Base pairs read = ", nbp

 1001 format (10I8)   !Format used in FLAG pointers and residue pointers
 1002 format (20A4)   !Format used in atom and residue names

  end subroutine topology_amber
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
 !Reads the atom coordinates from WrLINE coordinate file (amber style)

 subroutine coordinates_amber_WrLINE(coords,N,frames,traj,box,str)
 real(dp), allocatable :: coords(:,:,:), b(:)
 integer :: nlin, ncol, frames, N, box, str, ierror, i, k, l, rl, &
          & rn, an
 character(360) :: traj, aux
 logical :: amber, col3, xyz, good  !Logicals to determine types of trajectory

 !It is a mess, there are 3 possible formats!

 !At the beginning, everything is false!
 amber = .false.
 col3 = .false.
 xyz = .false.
 good = .false. !Good will be true if we could read the trajectory
                ! hence, determining the format
 
 an = len(trim(traj)) !Length of string

 !Let's determine the format judging by the file type (only analysing the ending
 if (trim(traj(an:an)) .eq. "x" .or. trim(traj(an-2:an)) .eq. "crd") amber = .true.
 if (trim(traj(an-2:an)) .eq. "col") col3 = .true.
 if (trim(traj(an-2:an)) .eq. "xyz") xyz = .true.

 !First one, which is Amber trajectory style which ends with crd or x

 if (amber) then

   if (str == 0) then
     box = 0
     write(6,*) "Assuming there is no box in amber trajectory file"
   end if

   !get dimensions of trajectory file
   call arrdims(traj,nlin,ncol)

   ncol=3*N/10 !Let's recycle ncol
   if (3*N-10*ncol > 0) then
     ncol = ncol+1   !the extra line
   end if

   !Get number of frames in file

   if (box /= 0) then           !Just in case let's put it here
     frames = nlin/(ncol+1)
   else
     frames = nlin/ncol
   end if

   write(6,"(1A20,1I10)") "Number of frames =", frames
   !Lets recycle nlin, it will be the number of lines that can be readed completly by 3
   nlin = N/10  

   !rl is the number of remaining lines
   rl=3*mod(N,10)
   if ( mod(rl,10) == 0 ) then
     rl = rl/10
   else
     rl = 1+ rl/10
   end if

   !rn is number of atoms remaining
   rn = N - 10*nlin    !Its 10 atoms per nlin (3 lines)
 
   allocate(coords(3,N,frames), b(rn*3), stat=ierror)
   if (ierror /= 0) stop "Error in allocating coordinates array"
 
   open(unit=10, file=trim(traj), status="old",action="read",iostat=ierror)
     if (ierror/=0) stop 'Error in opening trajectory file'

   !Begin reading
   read(10,*) aux !skips first line

   do k=1,frames

   l = 0

     do i=1,nlin
       read(10, 1003) coords(1:3,l+1,k), coords(1:3,l+2,k), &
                    & coords(1:3,l+3,k), coords(1,l+4,k)
 
       read(10, 1003) coords(2:3,l+4,k), coords(1:3,l+5,k), &
                    & coords(1:3,l+6,k), coords(1:2,l+7,k)

       read(10, 1003) coords(3,l+7,k), coords(1:3,l+8,k), &
                    & coords(1:3,l+9,k), coords(1:3,l+10,k)
       l = l+10
     end do !close i

     ! If there are more lines to read
     if (rl > 0) then
       !If there's one more
       if (rl >= 1) then
         if (rl < 2) then
           read(10, 1003) b(1:rn*3)
         else
           read(10, 1003) b(1:10)
         end if
       end if 
    
       !If there's more than one
       if (rl >= 2) then
         if (rl < 3) then
           read(10, 1003) b(11:rn*3)
         else
           read(10, 1003) b(11:20)
         end if
       end if 
    
       !If there's 3 left
       if (rl == 3) then
         read(10, 1003) b(21:rn*3)
       end if 
      
       !Now, collect the remaining coordinates
       do i=1,rn
         coords(1,l+i,k) = b((i-1)*3+1)
         coords(2,l+i,k) = b((i-1)*3+2)
         coords(3,l+i,k) = b((i-1)*3+3)
       end do

     end if !close rl>0 
    
     !Just in case that there is a box, skip this line
     if (box /= 0) then
       read(10,*) aux
     end if
 
   end do !close k

   close(10, iostat=ierror) !close trajectory
   if (ierror /= 0) stop "Error in deallocating coordinates array"

   good = .true. !We could read the trajectory!!!!!!!!

 end if

 !Format 3col, this is the fastest format!
 if (col3) then
   !get dimensions of trajectory file
   call arrdims(traj,nlin,ncol)

   !Get number of frames in file

   frames = nlin/N 

   !Frames and bp should be a multiple of lines in 3col
   if (frames*N /= nlin) stop "Coordinates in .3col does not match number of bp"

   write(6,"(1A20,1I10)") "Number of frames =", frames

   allocate(coords(3,N,frames), stat=ierror)
   if (ierror /= 0) stop "Error in allocating coordinates array"


   !Time to read! 
   open(unit=10, file=trim(traj), status="old",action="read",iostat=ierror)
     if (ierror/=0) stop 'Error in opening trajectory file'

     do k=1,frames
       do i=1,N
         read(10,*) coords(:,i,k)
       end do
     end do
   
   close(10, iostat=ierror) !close trajectory
   if (ierror /= 0) stop "Error in deallocating coordinates array"

 
   good = .true. !We could read the trajectory!!!!!!!!

 end if !Close 3col reading  !!


 !xyz format! Not as fast as 3col but better than crd
 if (xyz) then

   !get dimensions of trajectory file
   call arrdims(traj,nlin,ncol)

   !Get number of frames in file

   frames = nlin/(N+2) !This format has +2 lines for each frame 

   !Frames and N+2 should be a multiple of lines in xyz
   if (frames*(N+2) /= nlin) stop "Coordinates in .xyz does not match number of bp"

   write(6,"(1A20,1I10)") "Number of frames =", frames

   allocate(coords(3,N,frames), stat=ierror)
   if (ierror /= 0) stop "Error in allocating coordinates array"


   !Time to read! 
   open(unit=10, file=trim(traj), status="old",action="read",iostat=ierror)
     if (ierror/=0) stop 'Error in opening trajectory file'

     do k=1,frames
         read(10,*) !Ignoring first 2 lines of each frame
         read(10,*) 
       do i=1,N
         read(10,*) aux, coords(:,i,k)
       end do
     end do
   
   close(10, iostat=ierror) !close trajectory
   if (ierror /= 0) stop "Error in deallocating coordinates array"

 
   good = .true. !We could read the trajectory!!!!!!!!
 end if

 !Warnings!!!!!!1
 if (.not. good) then
   write(6,*)  "Couldn't identify format in trajectory file"
   write(6,*)  "Please be sure that file ends with crd or x if it has an amber style"
   write(6,*) "trajectory format"
   write(6,*) "Or if 3col for *.3col format (3 columns"
   write(6,*) "Of xyz which ends with *.xyz"
 end if

 1003 format (10F8.3) !Format used in amber trajectory file 

 end subroutine coordinates_amber_WrLINE
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This subroutine writes coordinates from last snapshot. Both
 !Projected and centred (not projected but centred at the origin, where
 !the centroid of the fitted atoms is zero

 !F_AMB_T is amber trajectory format. Variable is stored in parms.f90
 subroutine write_crds( centred, projected, N )

 implicit none
 integer, intent(in) :: N
 real(dp), intent(in) :: centred(3,N), projected(3,N)
 character(360) :: c_crd, p_crd
 integer :: ierror

 c_crd = "centred.crd"
 p_crd = "projected.crd"

 !Let's write centred coordinates
 open(unit=10, file=trim(c_crd), status="replace", action="write", iostat=ierror)
 if (ierror /= 0) stop "Error in openning centred.crd"
 
 !Header
 write(10,*) "Created by SerraLINE"

 !This is risky, but let's do it simple
 write(10,trim(F_AMB_T)) centred(:,:)

 close(10, iostat=ierror)
 if (ierror /= 0) stop "Error in closing centred.crd"

 !Let's write projected
 open(unit=11, file=trim(p_crd), status="replace", action="write", iostat=ierror)
 if (ierror /= 0) stop "Error in openning projected.crd"

 !Header
 write(11,*) "Created by SerraLINE"

 !Again, a bit risky
 write(11,trim(F_AMB_T)) projected(:,:)

 close(11, iostat=ierror)
 if (ierror /= 0) stop "Error in closing projected.crd"

 end subroutine write_crds
!----------------------------------------------------------------------- 


!-----------------------------------------------------------------------
 !This subroutine is  writes the output file with the bending angles, 
 !widths and heights, This is for opened structures
 subroutine write_av_parms(nbp,ndim,mdim,frames,seq,avstd_bends, &
          & avstd_width, avstd_height, avstd_aratio,t_length, &
          & avg_dist, max_dist, avg_dist_rel, max_dist_rel)
 
 implicit none

 integer, intent(in) :: nbp, ndim, mdim, frames, t_length
 character(1), intent(in) :: seq(nbp)
 real(dp), intent(in) :: avstd_bends(2,mdim), avstd_width(2), &
                       & avstd_height(2), avstd_aratio(2), &
                       & avg_dist(2), max_dist(2), &
                       & avg_dist_rel(2), max_dist_rel(2)

 character(360) :: output_file
 integer :: i, j, l, ierror

 output_file = "SerraLINE.out"

 open(unit=10, file=trim(output_file), status="replace", action="write", iostat=ierror)
 if (ierror/=0) stop 'Error in opening output bending file'

  write(10,*) "PARAMETERS"
  write(10,*) ""
  write(10,*) "OPENED STRUCTURE"
  write(10,*) "METHOD: PROJECTION"
  write(10,*) "BASE PAIRS ", nbp!ndim
  write(10,*) "SEQUENCE"
  write(10,*) "", seq(1:nbp)
  write(10,F_PARM_5) "SNAPSHOTS ANALYSED", frames
  write(10,F_PARM_5) "TANGENT LENGTH:", t_length
  write(10,*) "First column averages, second column standard deviations"
!  write(10,F_PARM_4) "WIDTH: ", avstd_width(:)
!  write(10,F_PARM_4) "HEIGHT: ", avstd_height(:)
!  write(10,F_PARM_4) "ASPECT RATIO:", avstd_aratio(:)
  write(10,F_PARM_6) "WIDTH (Angstroms): ", avstd_width(:)
  write(10,F_PARM_6) "HEIGHT (Angstroms): ", avstd_height(:)
  write(10,F_PARM_6) "ASPECT RATIO:", avstd_aratio(:)
  write(10,F_PARM_6) "AVERAGE OF DISTANCES TO PLANE (Angstroms):", avg_dist(:)
  write(10,F_PARM_6) "AVERAGE OF MAXIMUM DISTANCES TO PLANE (Angstroms):", max_dist(:)
  write(10,F_PARM_6) "RELATIVE AVERAGE OF DISTANCES TO PLANE (%):", avg_dist_rel(:)
  write(10,F_PARM_6) "RELATIVE AVERAGE OF MAXIMUM DISTANCES TO PLANE (%):", max_dist_rel(:)
  write(10,*) ""
  l=0
  do j=1,ndim-1
    write(10,trim(F_PARM_3)) j+1, "mer"
    !write(10,trim(F_PARM_3)) j, "bp step"
    write(10,trim(F_PARM_1)) "base-step", "Bending angle"
    write(10,*) "--------------------------------------"

    do i=1,ndim-j
      l=l+1
      write(10,trim(F_PARM_2)) i, "-", i+j, " "," ",seq(i),seq(i+j), " "," ", &
               & avstd_bends(1,l), avstd_bends(2,l)
    end do
    write(10,*) ""
  end do

 close(10) 

 end subroutine write_av_parms

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This subroutine is similar to write_av_parms subroutine but this one is for
 ! closed structures, which means that per every length l, there are 
 ! nbp parameters. (There are n lengths also)
 subroutine write_c_av_parms(nbp,frames,seq,avstd_bends,avstd_width, &
                           & avstd_height,avstd_aratio,t_length, &
                           & avg_dist, max_dist, &
                           & avg_dist_rel, max_dist_rel)
 implicit none

 integer, intent(in) :: nbp, frames, t_length
 character(1), intent(in) :: seq(nbp)
 real(dp), intent(in) :: avstd_bends(2,nbp,nbp-1), avstd_width(2), &
                       & avstd_height(2), avstd_aratio(2), &
                       & avg_dist(2), max_dist(2), &
                       & avg_dist_rel(2), max_dist_rel(2)

 character(1) :: c_seq(2*nbp)
 character(360) :: output_file
 integer :: i, l, s, ierror

 output_file = "SerraLINE.out"
 
 c_seq(1:nbp) = seq(1:nbp)

 c_seq(nbp+1:2*nbp) = seq(1:nbp)

 open(unit=10, file=trim(output_file), status="replace", action="write", iostat=ierror)
 if (ierror/=0) stop 'Error in opening output bending file'

  write(10,*) "PARAMETERS"
  write(10,*) ""
  write(10,*) "CLOSED STRUCTURE"
  write(10,*) "METHOD: PROJECTION"
  write(10,*) "BASE PAIRS ", nbp
  write(10,*) "SEQUENCE"
  write(10,*) "", seq(1:nbp)
  write(10,F_PARM_5) "SNAPSHOTS ANALYSED", frames
  write(10,F_PARM_5) "TANGENT LENGTH:", t_length
  write(10,*) "First column averages, second column standard deviations"
!  write(10,F_PARM_4) "WIDTH: ", avstd_width(:)
!  write(10,F_PARM_4) "HEIGHT: ", avstd_height(:)
!  write(10,F_PARM_4) "ASPECT RATIO:", avstd_aratio(:)
  write(10,F_PARM_6) "WIDTH (Angstroms): ", avstd_width(:)
  write(10,F_PARM_6) "HEIGHT (Angstroms): ", avstd_height(:)
  write(10,F_PARM_6) "ASPECT RATIO:", avstd_aratio(:)
  write(10,F_PARM_6) "AVERAGE OF DISTANCES TO PLANE (Angstroms):", avg_dist(:)
  write(10,F_PARM_6) "AVERAGE OF MAXIMUM DISTANCES TO PLANE (Angstroms):", max_dist(:)
  write(10,F_PARM_6) "RELATIVE AVERAGE OF DISTANCES TO PLANE (%):", avg_dist_rel(:)
  write(10,F_PARM_6) "RELATIVE AVERAGE OF MAXIMUM DISTANCES TO PLANE (%):", max_dist_rel(:)
  write(10,*) ""
  do l=1,nbp-1
    write(10,trim(F_PARM_3)) l+1, "mer"
    !write(10,trim(F_PARM_3)) l, "bp step"
    write(10,trim(F_PARM_1)) "base-step", "Bending angle"
    write(10,*) "--------------------------------------"

    do i=1,nbp

      s = i + l
      if (s > nbp) then
        s = s - nbp
      end if

      write(10,trim(F_PARM_2)) i, "-", s, " "," ",c_seq(i),c_seq(s), " "," ", &
               & avstd_bends(1,i,l), avstd_bends(2,i,l)
    end do
    write(10,*) ""
  end do

 close(10) 

 end subroutine write_c_av_parms

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 ! This subroutine  writes the output file with the bending angles
 ! and is for opened structures
 subroutine write_av_bends(nbp,ndim,mdim,frames,seq,avstd_bends,t_length)
 implicit none

 integer, intent(in) :: nbp, ndim, mdim, frames, t_length
 character(1), intent(in) :: seq(nbp)
 real(dp), intent(in) :: avstd_bends(2,mdim)

 character(360) :: output_file
 integer :: i, j, l, ierror

 output_file = "SerraLINE.out"

 open(unit=10, file=trim(output_file), status="replace", action="write", iostat=ierror)
 if (ierror/=0) stop 'Error in opening output bending file'

  write(10,*) "PARAMETERS"
  write(10,*) ""
  write(10,*) "OPENED STRUCTURE"
  write(10,*) "METHOD: WITHOUT PROJECTION"
  write(10,*) "BASE PAIRS ", nbp!ndim
  write(10,*) "SEQUENCE"
  write(10,*) "", seq(1:nbp)
  write(10,F_PARM_5) "SNAPSHOTS ANALYSED", frames
  write(10,F_PARM_5) "TANGENT LENGTH:", t_length
  write(10,*) "First column averages, second column standard deviations"
  write(10,*) ""
  l=0
  do j=1,ndim-1
    write(10,trim(F_AVB_3)) j+1, "mer"
    !write(10,trim(F_AVB_3)) j, "bp step"
    write(10,trim(F_AVB_1)) "base-step", "Bending angle"
    write(10,*) "--------------------------------------"

    do i=1,ndim-j
      l=l+1
      write(10,trim(F_AVB_2)) i, "-", i+j, " "," ",seq(i),seq(i+j), " "," ", &
               & avstd_bends(1,l), avstd_bends(2,l)
    end do
    write(10,*) ""
  end do

 close(10)

 end subroutine write_av_bends

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 ! This subroutine is similar to write_av_bends subroutine but this one is for
 ! closed structures, which means that per every length l, there are 
 ! nbp parameters. (There are n lengths also)
 subroutine write_c_av_bends(nbp,frames,seq,avstd_bends,t_length)
 implicit none

 integer, intent(in) :: nbp, frames, t_length
 character(1), intent(in) :: seq(nbp)
 real(dp), intent(in) :: avstd_bends(2,nbp,nbp-1)

 character(1) :: c_seq(2*nbp)
 character(360) :: output_file
 integer :: i, l, s, ierror

 output_file = "SerraLINE.out"

 c_seq(1:nbp) = seq(1:nbp)

 c_seq(nbp+1:2*nbp) = seq(1:nbp)

 open(unit=10, file=trim(output_file), status="replace", action="write", iostat=ierror)
 if (ierror/=0) stop 'Error in opening output bending file'

  write(10,*) "PARAMETERS"
  write(10,*) ""
  write(10,*) "CLOSED STRUCTURE"
  write(10,*) "METHOD: WITHOUT PROJECTION"
  write(10,*) "BASE PAIRS ", nbp
  write(10,*) "SEQUENCE"
  write(10,*) "", seq(1:nbp)
  write(10,F_PARM_5) "SNAPSHOTS ANALYSED", frames
  write(10,F_PARM_5) "TANGENT LENGTH:", t_length
  write(10,*) "First column averages, second column standard deviations"
  write(10,*) ""
  do l=1,nbp-1
    write(10,trim(F_AVB_3)) l+1, "mer"
    !write(10,trim(F_AVB_3)) l, "bp step"
    write(10,trim(F_AVB_1)) "base-step", "Bending angle"
    write(10,*) "--------------------------------------"

    do i=1,nbp

      s = i + l
      if (s > nbp) then
        s = s - nbp
      end if

      write(10,trim(F_AVB_2)) i, "-", s, " "," ",c_seq(i),c_seq(s), " "," ", &
               & avstd_bends(1,i,l), avstd_bends(2,i,l)
    end do
    write(10,*) ""
  end do

 close(10)

 end subroutine write_c_av_bends

!-----------------------------------------------------------------------
! Obtains the lines and columns of the file "filen"
 subroutine arrdims(filen,nlin,ncol)
 character(260), intent(in) :: filen
 integer, intent(out) :: nlin,ncol
 character(260) :: fwc
 integer :: ierr

   fwc='.functionsmod_wc-ou.tmp'   ! auxilliary file, where we store number
                                   ! of lines and columns separated by newline

   call system('head -n 1 '//filen//' | wc -w  > '//trim(fwc),ierr)
   if (ierr == 0) then
     call system('cat '//filen//' | wc -l  >> '//trim(fwc),ierr)
     if (ierr == 0) then
       open(unit=10,file=trim(fwc),status='old',action='read')
       read(10,*) ncol
       read(10,*) nlin        ! We obtain the number of lines and number of
       close(10)              ! columns of the data file
     else                
       stop                
     end if
   else
     stop
   end if

   if (max(nlin,ncol) == 0) then
     write(6,'(a)') 'Error in arrdims. Did not get valid dimensions'
     write(6,'(a,i10,a,i10)') '  nlin=',nlin,'   ncol=',ncol
     stop
   end if
 end subroutine arrdims

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !Subroutine used by SerraLINE.
  ! Reads directory of trajectory (traj) and topology (top) files.
  ! Also, it returns a logical variable circle_str which indicates
  ! If the structure being analysed is a closed structure.
  ! and if the topology contains single or double stranded DNA, even
  ! if the WrLINE only analysis double-stranded DNA.
  ! bp_fitting indicates which bp to fit the plane
  ! print_print indicates if we are going to print the projected trajectory
  ! xyz indicates if in xyz format or crd (if falsE)
  ! t_length specifies the length between base-pairs to calculate the
  ! tangent vectors.
  ! All this information is indicated in the input file "s_line.in
  ! If topology is not provided, SerraLINE can still run but without
  ! Detail in the sequence

  subroutine SerraLINE_inputs(traj,top,circle_str,str,nbp,bp_fitting,fitplane, &
                            & print_proj,xyz,t_length)
  implicit none

  !top, traj : directories
  !circle_str : true if is closed structure
  !fitplane : if true, then a plane will be fitted
  character(1200), intent(out) :: top, traj
  logical, intent(out) :: circle_str, fitplane, print_proj, xyz
  integer, intent(out) :: str, t_length, nbp
  logical :: not_good
  integer, allocatable :: bp_fitting(:)
  character(1200) :: strinng, aux
  integer :: stra, i, auxi, ierror

  !Closed structure?
  call nextpar
  read(5,*,iostat=ierror) stra
  if (ierror/=0) stop "Invalid input. Tell me if the fragment is closed or opened"
  if ( (stra /= 0) .and. (stra /= 1) ) stop "Please type 0 for opened or 1 for closed structure"

  !Single or double-stranded topology?
  call nextpar
  read(5,*,iostat=ierror) str
  if (ierror/=0) stop "Invalid input. Tell me if there is a topology"
  if ( (str < 0) .or. (str > 2) ) stop "Please type 1 for single stranded or 2 for double " // &
                                      & "or 0 if topology is not being provided"

  !Number of base-pairs
  call nextpar
  read(5,*) aux
  if (str == 0) then
    read(aux,*,iostat=ierror) nbp
    if (ierror /= 0) stop "Please type a valid input for number of base-pairs"
    if (nbp <= 0) stop "Number of base-pairs should be greater than 0"
  end if

  !Mask
  call nextpar
  read(5,"(A)",iostat=ierror) strinng
  if (ierror /=0) stop "Error in reading fitting options"

  not_good = .true.
  fitplane = .true.
  do i=1,len(strinng)
    if ( strinng(i:i) == "'"  ) then
      call read_fitting_bps(strinng, bp_fitting)   !If yes, read string (will be complicated)
      not_good = .false.
      exit
    else if (strinng(i:i) .eq. "0" ) then          !If not
      not_good = .false.
      exit
    else if (strinng(i:i) .eq. "1" ) then          !If not
      not_good = .false.
      fitplane = .false.                           !Fitting will not be performed
      exit
    end if
  end do

  if (not_good) stop "Could not determined bps for fitting"
 
  !Tangent length
  call nextpar
  read(5,*,iostat=ierror) t_length
  if (ierror /=0) stop "Error in reading tangent length input"
  if ( t_length <=0 ) stop "Please type tangent length > 0"

  !Topology
  call nextpar
  read(5,"(A)",iostat=ierror) top
  if (ierror /=0) stop "Error in reading topology input"
  top = adjustl(top)
 
  !Trajectory
  call nextpar
  read(5,"(A)",iostat=ierror) traj
  if (ierror /=0) stop "Error in reading trajectory input"
  traj = adjustl(traj)

  !Printing options
  call nextpar
  read(5,*,iostat=ierror) auxi, strinng
  if (ierror /=0) stop "Error in reading printing trajectory options"

  !Let's check printing options
  if (fitplane) then           !If a fitting is performed

     strinng = adjustl(strinng) !This will help

     !If printing projection
     if (auxi == 1) then

        if ( trim(strinng) .eq. "xyz") then
           print_proj = .true.
           xyz = .true.
        else if ( trim(strinng) .eq. "crd") then
           print_proj = .true.
           xyz = .false.
        else
           write(6,"(A)") "Warning, could not determine format printing option"
           write(6,"(A)") "The projected trajectory will not be written"
           print_proj = .false.
           xyz = .false. !Whatever
        end if

     !If not printing
     else if (auxi == 0) then

        print_proj = .false.
        xyz = .false. !Whatever, this logical will be ignored

     !Any other value of auxi will do nothing
     else

        write(6,"(A)") "Warning, could not determine printing option"
        write(6,"(A)") "The projected trajectory will not be written"
        print_proj = .false.
        xyz = .false. 

     end if

  !If we are not fitting we are not printing
  else

     print_proj = .false.
     xyz = .false.

  end if  

  !Finally, just write the info...
  !!--------------

  !Closed or opened....
  if (stra == 1) then
    circle_str = .true.
    write(6,"(A)") "CLOSED STRUCTURE"
  else if (stra == 0) then
    circle_str = .false.
    write(6,"(A)") "OPENED STRUCTURE"
  else
    stop "Tell me if its a closed or open structure (1 or 0)"
  end if           

  !Single or double...
  if (str == 1) then
    write(6,"(A)") "SINGLE STRANDED TOPOLOGY"
  !No point of more if's since we already have the warning above
  else if (str == 2) then
    write(6,"(A)") "DOUBLE STRANDED TOPOLOGY"
  else
    write(6,"(A)") "TOPOLOGY NOT PROVIDED"
  end if

  !Fitting?
  if (fitplane) then
    if (allocated(bp_fitting)) then
      write(6,"(A)") "Fitting will be applied considering specifc base-pairs"
    else 
      write(6,"(A)") "Fitting will be applied considering all base-pairs"
    end if
  else
    write(6,"(A)") "No fitting will be performed, width and height will not be obtained"    
  end if

  !Topology and trajectory directories 
  write(6,"(2A)") "Topology file =", trim(top) 
  write(6,"(2A)") "Trajectory file =", trim(traj) 

  !Printing projection?
  if (fitplane) then

    if (print_proj) then
      write(6,"(2A)") "Projected trajectory will be printed"
    else
      write(6,"(2A)") "Projected trajectory will not be printed"
    end if

  end if

  !Tangents length
  write(6,"(A,I5)") "Tangent length =", t_length

  !And it is done!

  end subroutine SerraLINE_inputs 
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  !Reads input file "extract.in" and it gets length to filter and
  !path to "serraline.out"
  subroutine extract_inputs(in_file, my_bp )
  implicit none

  character(1200), intent(out) :: in_file
  integer, intent(out) :: my_bp 

  call nextpar

  read(5,"(A)") in_file

  in_file = adjustl(in_file)

  !Gets bp of interest
  call nextpar

  read(5,*) my_bp

  end subroutine extract_inputs
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
  !Reads parameters in serraline.out and extracts them for a particular
  !length of interest (my_bp)
  subroutine read_parms( in_file, my_bp, N, nbp, parms, t_length, planefit, F_FORM )

  implicit none
  character(1200), intent(in) :: in_file
  character(1200), intent(out) :: F_FORM
  integer, intent(in) :: my_bp
  integer, intent(out) :: N, nbp, t_length
  real(dp), allocatable, intent(out) :: parms(:,:)
  logical, intent(out) :: planefit
  integer ::  i, j, l, n_parms, bp, ierror
  logical :: c_str
  character(360) :: aux
  character(50) :: F_FORM_1, F_FORM_2  !Formats
 
  !We don't care about any information, just about the parameters that we
  !want

  open(unit=10, file=trim(in_file), status="old", action="read", iostat=ierror)
  if (ierror/=0) then
    write(6,*) "Error in opening file", trim(in_file)
    stop
  end if

  read(10,"(A)") aux !PARAMETERS
  read(10,"(A)") aux !""
  !READ CLOSED OR OPENED STRUCTURE
  read(10,"(A)") aux
  if (aux(2:7) .eq. "CLOSED") then
    c_str = .true.
  else
    c_str = .false.
  end if
  !AND METHOD
  read(10,"(A)") aux
  if (aux(10:19) .eq. "PROJECTION") then
    planefit = .true.
  else
    planefit = .false.
  end if
  read(10,*) aux, aux, nbp !N bp
  read(10,"(A)") aux !SEQUENCE
  read(10,"(A)") aux !actual sequence
  read(10,"(A)") aux !snaps
  read(10,F_PARM_5) aux, t_length ! TANGENT LENGTH
  read(10,"(A)") aux !First column...
  if (planefit) then
    read(10,"(A)") aux ! WIDTH
    read(10,"(A)") aux ! HEIGHT
    read(10,"(A)") aux ! ASPECT RATIO
    read(10,"(A)") aux ! AVG DIST 
    read(10,"(A)") aux ! MAX DIST
    read(10,"(A)") aux ! AVG DIST R
    read(10,"(A)") aux ! MAX DIST R
  end if

  read(10,"(A)") aux !""

  !Let's check if we have a valid my_bp
  if (my_bp < 1) stop "Distance of interest must be positive and greater than 0"

  if (c_str) then    !If it is a circle
    if (my_bp > nbp-1) then
      write(6,*) "Distance of interest must not be greater than ",nbp-1
      stop
    end if
  else               !If linear
    if (my_bp > nbp-t_length) then
      write(6,*) "Distance of interest must not be greater than ",nbp-t_length
      stop
    end if
  end if

  !We have what we need for allocating our parms array
  !Let's check how many steps we have (depends on opened or closed structure)
  if (c_str) then
    N = nbp
  else
    N = nbp-my_bp - t_length !+ 1
  end if
  !And if we'll read sizes or just bendings
  !If a plane was fitted, then we have bendings, widths, heigths and aspect ratios
  !if not only bendings
  !Let's also determine the format that well use
  if (planefit) then
    n_parms = 2
    F_FORM_1 = F_PARM_2
    F_FORM_2 = F_PARM_3
  else
    n_parms = 2
    F_FORM_1 = F_AVB_2
    F_FORM_2 = F_AVB_3
  end if
  !NOTE that it is not arranged in the same way that in serraLINE. This is
  !because we only want values for one length
  allocate(parms(N,n_parms), stat=ierror)
  if (ierror/=0) stop "Error in allocating parameters"

  !Now, extract what we want
  do l=1,my_bp-1

    read(10,*) aux !bp
    read(10,*) aux !base-step...
    read(10,*) aux !-----

    !Check how many parameters we're going to ignore
    if (c_str) then
      j = nbp
    else
      j = nbp-l-t_length
    end if

    do i=1,j
      read(10,"(A)") aux !we don't care about this
    end do ! i

    read(10,"(A1)") aux !""

  end do ! l

  !Now, read what we actually want
  read(10,trim(F_FORM_2)) bp, aux !bp step or mer
  !Just to be sure
  !if (bp /= l .or. bp /= my_bp) then !In case of bp step
  if (bp-1 /= l .or. bp-1 /= my_bp) then !in case of bp+1 mer
    stop "Reading wrong data"
  end if

  read(10,*) aux !base-step...
  read(10,*) aux !-----
  do i=1,N
    !Let's ignore the first 4 terms (an recycle j)
    read(10,trim(F_FORM_1)) j, aux(1:1), j, aux(1:1), aux(1:1), &
        &  aux(1:1), aux(1:1), aux(1:1), aux(1:1), parms(i,:) 
  end do
  close(10)

  !Prepare output format

!  F_FORM(1:5) = "(I10,"
  F_FORM(1:7) = "(F10.1,"
  F_FORM(8:) = F_FORM_1(15:)

  !That;s it, we have what we need

  end subroutine read_parms
!-----------------------------------------------------------------------
 

!-----------------------------------------------------------------------
  !Reads chr, if first character is different from 'c' (comment line), then 
  !it exits and then next line should be readed
  subroutine nextpar
  implicit none
  character(len=1) :: chr
    do 
      read(5,'(A1)') chr
      if (chr(1:1) /= 'c') exit
    end do

  end subroutine nextpar
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !This subroutine basically allocates bp_fitting
  subroutine read_fitting_bps(strinng, bp_fitting)
  implicit none
  character(1200) :: strinng
  character(len=1200), dimension(100) :: more_c
  integer, allocatable, intent(out) :: bp_fitting(:)
  integer, allocatable :: temporal(:)
  integer :: i, j, k, l, s, p, a ,b, m, first, middle, last, ierror
  logical :: first_comma, escape

  escape = .false.
  k = 0
  more_c(1)(:) = strinng

  !First, find how many coordinates
  do i=1,len(strinng)

    if ( more_c(1)(i:i) == "'")  then   !Begins!
      p = i+1

      do l=1,100               !l line
        if (l > 1) then        !First line we already have it (strinng)
          read(5,"(A)") strinng!more_c(l)
          more_c(l)(:) = strinng
        end if

 !       first_comma = .true. !It is the first comma

        first = p
        do j= p, len(strinng)

          !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
          if (more_c(l)(j:j) == "," ) then !if comma
!            if (first_comma) k=k+1
            k = k+1
 !           first_comma = .false.
            first = j+1
 
          !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
          else if (more_c(l)(j:j) == ":") then
            middle = j
            do s=j+1,len(strinng)
              if (more_c(l)(s:s) == "," .or. more_c(l)(s:s) == "'" .or. more_c(l)(s:s) == "&" ) then
                last = s-1
              end if
            end do        !close s

            !Count how many numbers from a bp to b bp (a:b)
            read(more_c(l)(first:middle-1),*) a
            read(more_c(l)(middle+1:last),*) b
!            if (a > b) stop "Write fitting bp from lowest to highest value"
            k = k + abs(b - a)
            first = j+1
!            first_comma = .false.

          !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          else if (more_c(l)(j:j) == "&") then
            k = k+1
            exit                                !Just exit loop

          !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
          else if (more_c(l)(j:j) == "'") then
            k = k+1
            escape = .true. 
            exit

          else if (j==len(strinng) ) then
            stop "Write ' or & in bp fittings"

          end if
          !------------------------------------------------------------------------------------------ 

        end do           !close j
        
        p = 1

        if (escape) exit !escape l
      end do             !close l


      if (escape) exit   !escape i

    end if

  end do                 !close i

  if (k <= 2) stop "Invalid number of bp to be fitted"

  allocate(temporal(k), stat=ierror)
  if (ierror /= 0) stop "Error in allocating indeces for bp fittings"

  !Now, lest get the values
  escape = .false.
  k = 0
  temporal = 0   !To don't get
  !First, find how many coordinates
  do i=1,len(strinng)

    if ( more_c(1)(i:i) == "'")  then   !Begins!
      p = i+1

      do l=1,100               !l line

        first_comma = .true. !It is the first comma

        first = p
        do j= p, len(strinng)

          !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
          if (more_c(l)(j:j) == "," ) then !if comma
            k = k+1
            read(more_c(l)(first:j-1),*) temporal(k)
            first = j+1

          !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
          else if (more_c(l)(j:j) == ":") then
            middle = j
            do s=j+1,len(strinng)
              if (more_c(l)(s:s) == "," .or. more_c(l)(s:s) == "'" .or. more_c(l)(s:s) == "&" ) then
                last = s-1
              end if
            end do        !close s

            !Count how many numbers from a bp to b bp (a:b)
            read(more_c(l)(first:middle-1),*) a
            read(more_c(l)(middle+1:last),*) b
            do m = a, b-1
              k = k + 1
              temporal(k) = m
            end do
            first = j+1

          !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          else if (more_c(l)(j:j) == "&") then
            k = k + 1
            read(more_c(l)(first:j-1),*) temporal(k)
            exit                                !Just exit loop

          !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
          else if (more_c(l)(j:j) == "'") then
            k = k + 1
            read(more_c(l)(first:j-1),*) temporal(k)
            escape = .true. 
            exit

          else if (j==len(strinng) ) then
            stop "Write ' or & in bp fittings"

          end if
          !------------------------------------------------------------------------------------------ 
        end do           !close j
        
        p = 1

        if (escape) exit !escape l
      end do             !close l


      if (escape) exit   !escape i

    end if

  end do                 !close i


  !Let's clean up our indices (remove repetead, zero and negative indices)
  k = 0
  do i = 1, size(temporal)
    if (i > 1) then
      if (temporal(i) > 0 .and. all( temporal(1:i-1) /= temporal(i) ) ) k = k+1
    else
      if (temporal(i) > 0 ) k = k+1
    end if
  end do

  if (k <= 2) stop "Invalid number bp to be fitted"
 
  allocate(bp_fitting(k), stat=ierror)
  if (ierror /= 0) stop "Error in allocating bp fitting indices"

  !Now, let's fill them
  k = 0
  do i = 1, size(temporal)
    if (i > 1) then
      if (temporal(i) > 0 .and. all( temporal(1:i-1) /= temporal(i) ) ) then
        k = k+1
        bp_fitting(k) = temporal(i)
      end if
    else
      if (temporal(i) > 0 ) then
         k = k+1
        bp_fitting(k) = temporal(i)
      end if
    end if
  end do

  end subroutine read_fitting_bps
!----------------------------------------------------------------------

  end module IO_mod
