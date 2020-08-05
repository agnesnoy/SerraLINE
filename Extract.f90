!                                                              04/07/2020
! This program extracts information from serraline.out for a particular length,
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

  program Extract
  
  use IO_mod
  use parms

  implicit none

  integer :: i, j, N, nbp, my_bp, ierror
  character(1200) :: in_file, out_file, F_FORM !Format
  real(dp), allocatable :: parms(:,:), mid(:), & !Point in the middle of two bp
                         & aux_r1(:)
  real(dp) :: aux_r2                             
  logical :: planefit

 !INPUTS SECTION
 !-----------------------------------------------
 !Read extract.in and determine path to serraline.out and length to filter

 call extract_inputs( in_file, my_bp)

 !READING SECTION
 !-----------------------------------------------
 !Read serraline.out
 !The next subroutine figures out if is closed structure and if a plane was fitted
 call read_parms( in_file, my_bp, N, nbp, parms, planefit, F_FORM )

 !we didn't use planefit, but could be useful in the future (or not?)

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

  !Now write!
  do i=1,N
    write(11,trim(F_FORM)) mid(i), parms(i,:)
  end do
  close(11)
 
 !CLEAR DATA
 !-----------------------------------------------
 
  deallocate(parms, stat=ierror)
  if (ierror/=0) stop "Error in deallocating parms"
  
 !FORMATS
 !-----------------------------------------------
!  1006 format(A15,2F10.3)
!  1007 format(I10,A10)
!  1008 format(I10,2F10.3)
!-----------------------------------------------

  end program Extract
