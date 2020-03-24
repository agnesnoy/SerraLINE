!                                                              22/07/2019
! This program extracts information from serraline.out for a particular length,
  program Extract
  
  use IO_mod
  use parms

  implicit none

  integer :: i, j, N, nbp, my_bp, ierror
  character(360) :: in_file, out_file, F_FORM !Format
  real(dp), allocatable :: parms(:,:), mid(:), & !Point in the middle of two bp
                         & aux_r1(:)
  real(dp) :: aux_r2

 !INPUTS SECTION
 !-----------------------------------------------
 !Read extract.in and determine path to serraline.out and length to filter

 call extract_inputs( in_file, my_bp)

 !READING SECTION
 !-----------------------------------------------
 !Read serraline.out
 !The next subroutine figures out if is closed structure
 call read_parms( in_file, my_bp, N, nbp, parms, F_FORM )

 !WRITING SECTION
 !-----------------------------------------------
 !Now let's just write the results with no particular format

  !Let's first prepare the output file
  ! 1 digit
  if (my_bp < 10) then
    write(out_file,"(A12,I1,A4)") "subfragment_", my_bp, ".out"
  ! 2 digits
  else if (my_bp > 9 .and. my_bp < 100) then
    write(out_file,"(A12,I2,A4)") "subfragment_", my_bp, ".out"
  ! 3 digits
  else if (my_bp > 99 .and. my_bp < 1000) then
    write(out_file,"(A12,I3,A4)") "subfragment_", my_bp, ".out"
  ! 4 difits
  else if (my_bp > 999 .and. my_bp < 10000) then
    write(out_file,"(A12,I4,A4)") "subfragment_", my_bp, ".out"
  ! too big
  else
    stop "BP distance too big, modify Extract.f90"
  end if

  !Now, write
  open(unit=11, file=trim(out_file), status="replace", action="write", iostat=ierror)
  if (ierror/=0) then
    write(6,*) "Error in openining output file: ", trim(out_file)
    stop
  end if

  allocate(aux_r1( size(parms,2) ), mid(N), stat=ierror)
  if (ierror/=0) stop "Error in allocating mid points"

  !Let's rearranged our parms to get the angles from smaller to heigher positions
  do i=1,N

    !Calculate mid point
    mid(i) = real( 2*i + my_bp, dp) / 2.0_dp
    mid(i) = mid(i) + 0.5_dp
    if (  mid(i) - real(nbp,dp)  > eps ) mid(i) = mid(i) - real(nbp,dp)

  end do

  !Sort data
  do j = N-1, 1, -1
    do i = 1, j
      if (mid(i) > mid(i+1)) then
        aux_r2 = mid(i)
        aux_r1 = parms(i,:)
        mid(i) = mid(i+1)
        parms(i,:) = parms(i+1,:)
        mid(i+1) = aux_r2
        parms(i+1,:) = aux_r1
      end if
    end do
  end do

  !Now do the real writing
  do i=1,N
    write(11,trim(F_FORM)) mid(i), parms(i,:)
  end do
  close(11)
 
 !CLEAR DATA
 !-----------------------------------------------
 
  deallocate(parms, aux_r1, mid, stat=ierror)
  if (ierror/=0) stop "Error in deallocating parms"
  
 !FORMATS
 !-----------------------------------------------
!  1006 format(A15,2F10.3)
!  1007 format(I10,A10)
!  1008 format(I10,2F10.3)
!-----------------------------------------------

  end program Extract
