!------------------------------------------------------------- 22/07/2019
! This module contains functions and subroutines needed for reading
! and writing.
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
 ! This subroutine  writes the output file with the bending angles
 ! and is for opened structures
 subroutine write_av_bends(nbp,ndim,mdim,frames,seq,avstd_bends)
 implicit none

 integer, intent(in) :: nbp, ndim, mdim, frames
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
!  write(10,*) "METHOD: WITHOUT PROJECTION"
  write(10,*) "BASE PAIRS ", ndim
  write(10,*) "SEQUENCE"
  write(10,*) "", seq(1:nbp)
  write(10,*) "SNAPSHOTS ANALYSED", frames
  write(10,*) ""
  write(10,*) "First column averages, second column standard deviations"
  l=0
  do j=1,ndim-1
    write(10,trim(F_AVB_3)) j, "bp step"
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
 subroutine write_c_av_bends(nbp,frames,seq,avstd_bends)
 implicit none

 integer, intent(in) :: nbp, frames
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
!  write(10,*) "METHOD: WITHOUT PROJECTION"
  write(10,*) "BASE PAIRS ", nbp
  write(10,*) "SEQUENCE"
  write(10,*) "", seq(1:nbp)
  write(10,*) "SNAPSHOTS ANALYSED", frames
  write(10,*) ""
  write(10,*) "First column averages, second column standard deviations"
  do l=1,nbp-1
    write(10,trim(F_AVB_3)) l, "bp step"
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
  ! Subroutine used by SerraLINE.
  ! Reads directory of trajectory (traj) and topology (top) files.
  ! Also, it returns a logical variable circle_str which indicates
  ! If the structure being analysed is a closed structure.
  ! and if the topology contains single or double stranded DNA, even
  ! if the WrLINE only analysis double-stranded DNA.
  ! All this information is indicated in the input file "s_line.in
  subroutine SerraLINE_inputs(traj,top,circle_str,str,nbp)
  implicit none

  !top, traj : directories
  !circle_str : true if is closed structure
  character(360), intent(out) :: top, traj
  logical, intent(out) :: circle_str
  integer, intent(out) :: str, nbp
  character(360) :: aux
  integer :: stra,ierror

  !Closed structure?
  call nextpar
  read(5,*,iostat=ierror) stra
  if (ierror/=0) stop "Invalid input. Tell me if the fragment is closed or opened"
  if ( (stra < 0) .or. (stra > 1) ) stop "Please type 0 for opened and 1 for closed"
 
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

  !Topology
  call nextpar
  read(5,"(A)") top
  top = adjustl(top)
 
  !Trajectory
  call nextpar
  read(5,"(A)",iostat=ierror) traj
  if (ierror /= 0) stop "Please type a the path to the trajectory file"
  traj = adjustl(traj)

  !Write info...

  if (stra == 1) then
    circle_str = .true.
    write(6,"(A)") "CLOSED STRUCTURE"
  else if (stra == 0) then
    circle_str = .false.
    write(6,"(A)") "OPENED STRUCTURE"
  else
    stop "Tell me if its a closed or open structure (1 or 0)"
  end if           

  if (str == 1) then
    write(6,"(A)") "SINGLE STRANDED TOPOLOGY"
  !No point of more if's since we already have the warning above
  else if (str == 2) then
    write(6,"(A)") "DOUBLE STRANDED TOPOLOGY"
  else 
    write(6,"(A)") "TOPOLOGY NOT PROVIDED"
  end if

  if (str /=0) write(6,"(2A)") "Topology file =", trim(top) 
  write(6,"(2A)") "Trajectory file =", trim(traj) 

  end subroutine SerraLINE_inputs 
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  !Reads input file "extract.in" and it gets length to filter and
  !path to "serraline.out"
  subroutine extract_inputs(in_file, my_bp )
  implicit none

  character(360), intent(out) :: in_file
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
  subroutine read_parms( in_file, my_bp, N, nbp, parms, F_FORM ) 
  implicit none
  character(360), intent(in) :: in_file
  character(360), intent(out) :: F_FORM
  integer, intent(in) :: my_bp
  integer, intent(out) :: N, nbp
  real(dp), allocatable, intent(out) :: parms(:,:)
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
  read(10,*) aux, aux, nbp !N bp
  read(10,"(A)") aux !SEQUENCE
  read(10,"(A)") aux !actual sequence
  read(10,"(A)") aux !snaps
  read(10,"(A)") aux !""

  read(10,"(A)") aux !First column...

  !Let's check if we have a valid my_bp
  if (my_bp < 1) stop "Distance of interest must be positive and greater than 0"
  if (my_bp > nbp-1) then
    write(6,*) "Distance of interest must not be greater than ",nbp-1
  end if

  !We have what we need for allocating our parms array
  !Let's check how many steps we have (depends on opened or closed structure)
  if (c_str) then
    N = nbp
  else
    N = nbp-my_bp
  end if
  !And if we'll read sizes or just bendings
  !Let's also determine the format that well use
  n_parms = 2
  F_FORM_1 = F_AVB_2
  F_FORM_2 = F_AVB_3
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
      j = nbp-l
    end if

    do i=1,j
      read(10,"(A)") aux !we don't care about this
    end do ! i

    read(10,"(A1)") aux !""

  end do ! l

  !Now, read what we actually want
  read(10,trim(F_FORM_2)) bp, aux !bp
  !Just to be sure
  if (bp /= l .or. bp /= my_bp) then
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

  end module IO_mod
