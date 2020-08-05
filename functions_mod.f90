! Module containing all the mathematical functions needed for      03/04/2019
! SerraLINE
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

 module functions_mod

 use parms
 use IO_mod

 implicit none

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 contains

!-----------------------------------------------------------------------
 !Fits a global plane onto the coordinates, then project the coordinates
 !onto it
 !coords(3,N,frames) is the coordinates of N atoms in a number of frames
 !G_n is normal vector to the fitted plane and c is an auxiliar vector.
 !This subroutine uses the Jacobi diagonalization method to solve 
 !the SVD matrices (just one part), this is done to fit the best plane.
 subroutine project_coordinates_G_plane(coords, N, frames, n_fitting, bp_fitting, G_n, best)
 real(dp), intent(inout) :: coords(:,:,:)
 integer,  intent(in) :: N, frames, n_fitting
 integer, allocatable, intent(in) :: bp_fitting(:)
 integer, intent(out) :: best(:)
 real(dp), intent(out) :: G_n(:,:,:)
 real(dp) :: b(3), c(3), A(3,N), U(3,3), S(3,3), &
           & d(3)
 real(dp), allocatable :: A_fit(:,:)
 integer :: i, j, k, ierror

 !A are the coordinates centered at the origin.
 !A_p are the projected coordinates to some plane (U perpendicular)
 !d will be 3 distances, that will help to decide which projection is the
 !best fit

  if ( allocated(bp_fitting) ) then

    allocate(A_fit(3,n_fitting), stat = ierror)
    if (ierror/=0) stop "Error in allocating A_fit"

    ! Fit plane to all coordinates--------------------------------------------
    do k=1,frames
      !1.- Get centroid
      c = 0.0_dp
      do i = 1,n_fitting
        c = coords(:,bp_fitting(i),k) + c
      end do
      c = c/real(n_fitting,dp)

      !2.- Center coordinates
      do i = 1,N
        A(:,i) = coords(:,i,k) - c(:)
      end do

      do i = 1,n_fitting
        A_fit(:,i) = A(:,bp_fitting(i))
      end do

      !3.-Solve SVD
      ! A is not symmetric, but A*At and At*A are...
      !Only U is neededl...
      call diagonalization_Jacobi( matmul(A_fit,transpose(A_fit)),S,U,3)

      d(1) = S(1,1)
      d(2) = S(2,2)
      d(3) = S(3,3)

      !Find eigenvalue
      best(k) = 0
      do i=1,3
        if (d(i) == minval(d)) best(k) = i
      end do
      if (best(k) <= 0 .or. best(k) >3) stop "Error finding best plane"


     !Project coordinates to the plane j
      do j=1,3
        do i=1,N
          b = A(:,i) - projv_u( A(:,i), U(:,best(k)) )
 
          coords(:,i,k) = b                 !New coordinates
        end do
      end do

      !4.- We know the best plane [best(k)], now let's store
      ! the planes
      G_n(:,:,k) = U(:,:)

     end do !close k

  else
    ! Fit plane to all coordinates--------------------------------------------
    do k=1,frames
      !1.- Get centroid
      c(1) = sum( coords(1,:,k) )
      c(2) = sum( coords(2,:,k) )
      c(3) = sum( coords(3,:,k) )
      c = c/real(N,dp)

      !2.- Center coordinates
      do i = 1,N
        A(:,i) = coords(:,i,k) - c(:)
      end do

      !3.-Solve SVD
      ! A is not symmetric, but A*At and At*A are...
      !Only U is neededl...
      call diagonalization_Jacobi( matmul(A,transpose(A)),S,U,3)

      d(1) = S(1,1)
      d(2) = S(2,2)
      d(3) = S(3,3)

      !Find eigenvalue
      best(k) = 0
      do i=1,3
        if (d(i) == minval(d)) best(k) = i
      end do
      if (best(k) <= 0 .or. best(k) >3) stop "Error finding best plane"


     !Project coordinates to the plane j
      do j=1,3
        do i=1,N
          b = A(:,i) - projv_u( A(:,i), U(:,best(k)) )
 
          coords(:,i,k) = b                 !New coordinates
        end do
      end do

      !4.- We know the best plane [best(k)], now let's store
      ! the planes
      G_n(:,:,k) = U(:,:)
   
    end do !close k

  end if !End if (allocated)-----------------------------------------------

  !Now, let's jump to module IO_mod to write crd files before erasing
  !what is not needed
!  call write_crds( A, coords(:,:,frames), N ) Uncomment this line to
                                             ! print centered coordinate

 end subroutine project_coordinates_G_plane
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Only works with real symmetric matrix "M" (ndim x ndim)
! This one returns the Diagonalization matrix S with eigenvalues in the diagonal
! and the eigenvectors V. Eigenvalue S(i,i) corresponds to eigenvector V(:,i)
! It uses the Jacobi diagonalization method

 subroutine diagonalization_Jacobi(M,S,V,ndim)
 integer, intent(in) :: ndim
 real(dp), intent(in) :: M(ndim,ndim)
 real(dp), intent(out) :: S(ndim,ndim),V(ndim,ndim)
 integer :: i,j,l,k, r, cont
 real(dp) :: larg, th, n, cond, co, si,  &
             upS(ndim,ndim), Image(ndim,ndim) ,&
             G(ndim,ndim)

 n = real(ndim,dp)
 S = M
 cont = 50
 Image = 0.0_dp
 do k=1,ndim
   Image(k,k) = 1.0_dp
 end do

 V = Image
 cond = 4 !Just because
 r = 0
 do
   if (cond <= eps .or. r==cont) then
     S = upS
     exit
   end if
   if (r/=0) then
     S=upS
   end if
   r=r+1
   larg=-1.0_dp
   do k=1,ndim-1
     do l=k+1,ndim
       if ( larg < abs( S(k,l) ) ) then
         larg = abs( S(k,l) )
         i = k
         j = l
       end if
     end do
   end do
   if ( S(i,i) == S(j,j) ) then
     th= pi/4.0_dp
   else
     th=0.5_dp*datan( 2.0_dp*S(i,j)/( S(j,j) - S(i,i) ) )
   end if
   co = dcos(th)
   si = dsin(th)
   G = Image
   G(i,i) = co
   G(j,j) = co
   G(i,j) = si
   G(j,i) = -si
   V = MATMUL(V,G)
   upS = MATMUL(transpose(G),S)
   upS = MATMUL(upS,G)
   cond = 0.0_dp
   do k=1,ndim-1
     do l=k+1,ndim
       cond = cond + abs( upS(k,l) )
     end do
   end do
   cond=cond*2.0_dp
 end do

 end subroutine diagonalization_Jacobi
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This subroutine is used in SerraLINE
 subroutine get_tangent_vectors(coords,N,frames,circle_str,tangents)
 
 implicit none
 
 logical, intent(in) :: circle_str
 integer, intent(in) :: N,frames
 real(dp), intent(in) :: coords(:,:,:)
 real(dp), intent(out) :: tangents(:,:,:)
 integer :: i,k

 do k=1,frames

  do i=1,N-1
    tangents(:,i,k) = coords(:,i+1,k) - coords(:,i,k)
    tangents(:,i,k) = normalize_vector(tangents(:,i,k),3)
  end do

  !If its a closed structure then the last bp has tangent vector
  if (circle_str) then
    tangents(:,N,k) = coords(:,1,k) - coords(:,N,k)
    tangents(:,N,k) = normalize_vector(tangents(:,N,k),3)
  end if

 end do

 end subroutine get_tangent_vectors
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Calculates widths, heights and aspect ratios (absolute values) of each projected 
 !trajectory (frame). The vector from bp i to bp j is projected to
 !the two perpendicular components of plane G_n.
 !Coordinates (coord) are already projected to plane G_n(:,best,k)
 !ndim -> number of maximum parameters at length 1
 !mdim -> number of parameters in frame k
 subroutine get_width_height_aratio(ndim,mdim,frames,G_n,best,coords,width,height,aratio)

 implicit none
 
 integer, intent(in) :: ndim, mdim, frames, best(:)
 real(dp), intent(in) :: coords(:,:,:), G_n(:,:,:)
 real(dp), intent(out) :: width(:), height(:), aratio(:)
 real(dp) :: r(3),xy(3,2), G_k(3,2), & !x=1,y=2, G_k rotated plane
           & d, d_op, d_r(3), & !square distances: absolute of normal, largest, and vector of largest
           & awi(mdim), ahe(mdim)   !auxiliars for calculating sizes
 integer :: i, j, k, l, s


 do k=1,frames       !frame k

    !First, let's rotate plane so the height axis aligns with
    !the largest square value distance
    d_op = -1.0_dp !Just the get the first distance as the greatest
    d_r = 0.0_dp

    s = 0
    do l=1,ndim-1     !length l

      do i=1,ndim-l   !from bp i to j

        s = s + 1
        j = i + l

        !r is the vector from bp i to j, l is the length in bp
        r = coords(:,j,k) - coords(:,i,k)

        d =  absv2(r)

        if (d > d_op) then
          d_op = d
          d_r = r
        end if

      end do !close i
    end do   !close l
 
    !We have the largest distance locations d_loc.
    !No need for rotation, the height axis is parallel to r and
    !width is perpendicualr

    !Direction of largest (normalize d_r)
    r = normalize_vector(d_r,3)

    G_k(:,1) = r
    G_k(:,2) = cross_product3( r, G_n(:,best(k),k) )

    !Now, get width, height and aspect ratios between all bps
    s=0
    do l=1,ndim-1     !length l
      do i=1,ndim-l   !from bp i to j
        s=s+1         !parameter s
        j = i+l

        !r is the vector from bp i to j
        r = coords(:,j,k) - coords(:,i,k)

        !calculate components of r perpendicular to plane G_n(best)
        xy(:,1) = projv_u( r, G_k(:,1) ) 
        xy(:,2) = projv_u( r, G_k(:,2) ) 

        ahe(s) = absv( xy(:,1) )        !First
        awi(s)  = absv( xy(:,2) )        !Second projection

      end do !close i
    end do   !close l

    !Save maximums width and height
    width(k)  = maxval(awi)
    height(k) = maxval(ahe)

    !Now the global
    aratio(k) = width(k)/height(k)

  end do     !close k

 end subroutine get_width_height_aratio
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Subroutine like the one above but for closed structures.
 !Calculates width, height and aspect ratio (absolute values) of each projected 
 !trajectory (frame). The vector from bp i to bp j is projected to
 !the two perpendicular components of plane G_n.
 !Coordinates (coord) are already projected to plane G_n(:,best,k)
 !Now, each plane is rotated so a component aligns with the highest distance
 subroutine get_c_width_height_aratio(ndim,frames,G_n,best,coords,width,height,aratio)

 implicit none
 
 integer, intent(in) :: ndim, frames, best(:)
 real(dp), intent(in) :: coords(:,:,:), G_n(:,:,:)
 real(dp), intent(out) :: width(:), height(:), aratio(:)
 real(dp) :: r(3),xy(3,2), G_k(3,2), & !x=1,y=2, G_k rotated plane
           & d, d_op, d_r(3), & !square distances: absolute of normal, largest, and vector of largest
           & awi(ndim,ndim-1), ahe(ndim,ndim-1)
 integer :: i, j, k, s, l

  do k=1,frames       !frame k

    !First, let's rotate plane so the height axis aligns with
    !the largest square value distance
    d_op = -1.0_dp !Just the get the first distance as the greatest
    d_r = 0.0_dp

    do l=1,ndim-1     !length l

      do i=1,ndim   !from bp i to j

        s = i + l

        !Check if bp i will be > j
        if (s > ndim) then
          j = s - ndim
        else
          j = s
        end if

        !r is the vector from bp i to j, l is the length in bp
        r = coords(:,j,k) - coords(:,i,k)


        d =  absv2(r)

        if (d > d_op) then
          d_op = d
          d_r = r
        end if
      end do !close i
    end do   !close l
 
    !We have the largest distance locations d_loc.
    !No need for rotation, the height axis is parallel to r and
    !width is perpendicualr

    !Direction of largest (normalize d_r)
    r = normalize_vector(d_r,3)

    G_k(:,1) = r
    G_k(:,2) = cross_product3( r, G_n(:,best(k),k) )
    
    !Now, get real width and heights
    do l=1,ndim-1     !length l

      do i=1,ndim   !from bp i to j

        s = i + l

        !Check if bp i will be > j
        if (s > ndim) then
          j = s - ndim
        else
          j = s
        end if

        !r is the vector from bp i to j
        r = coords(:,j,k) - coords(:,i,k)

        !calculate components of r perpendicular to plane G_n(best)
        xy(:,1) = projv_u( r, G_k(:,1) ) 
        xy(:,2) = projv_u( r, G_k(:,2) ) 

        ahe(i,l) = absv( xy(:,1) )       !First
        awi(i,l) = absv( xy(:,2) )       !Second projection

      end do !close i
    end do   !close l

    !Get sizes (maximums
    width(k)  = maxval(awi)
    height(k) = maxval(ahe)

    !Get global aspect ratio. max width divided by max height in frame k
    aratio(k) = width(k)/height(k)

  end do     !close k

 end subroutine get_c_width_height_aratio
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
 !This subroutine calculates the bending angles between tangent vectors
 ! The output is in degrees
 ! This subroutine is fo a opened structure, that means that it will have
 ! n-1 possible lengths and for each length it have n(n-1)/2 parameters
 ! Bending angles can be greater than 180 degrees, since these angles
 ! are in the same plane, orientations are calculated using the normal
 ! vector to the plane
 subroutine get_bendings(ndim,frames,G_n,best,tangents,bends)
 
 implicit none
 
 integer, intent(in) :: ndim, frames, best(:)
 real(dp), intent(in) :: G_n(:,:,:), tangents(:,:,:)
 real(dp), intent(out) :: bends(:,:)
 real(dp) :: a, b, c(3), d
 integer :: i, j, k, l

  do k=1,frames
    l = 0
    j = 1  !First lengths
    do i=1,ndim-j
      l=l+1
      !Define sign
      c = cross_product3( tangents(:,i,k), tangents(:,i+j,k) )
      d = dot_product( tangents(:,i,k),tangents(:,i+j,k) ) !auxiliar
      !The first two ifs are unlikely, but sometimes happen and
      !cause error...
      if ( d > 1.0_dp  ) then
        b = 0.0_dp
      else if ( d < -1.0_dp  ) then
        b = pi
      else
        a = dot_product(tangents(:,i,k),tangents(:,i+j,k))
        b = dacos(a)
      end if

      !Change sign if needed
      if ( dot_product( c, G_n(:,best(k),k) ) < 0.0_dp ) b = -b

      bends(k,l) = rad_to_deg*b !Convert to degrees

    end do   !i

    do j=2,ndim-1
      do i=1,ndim-j
        l = l + 1
        bends(k,l) = sum( bends(k,i:i+j) )
      end do !i
    end do   !j

  end do     !k

  bends = abs(bends)

 end subroutine get_bendings
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This subroutine calculates the bending angles between tangent vectors
 ! The output is in degrees
 ! This subroutine is for a opened structure, that means that it will have
 ! n-1 possible lengths and for each length it have n(n-1)/2 parameters
 ! This subroutine is for not projected trajectories, this means that
 ! bendings won't be greater than 180 degrees.
 subroutine get_bendings_nofit(ndim,frames,tangents,bends)

 implicit none

 integer, intent(in) :: ndim,frames
 real(dp), intent(in) :: tangents(:,:,:)
 real(dp), intent(out) :: bends(:,:)
 real(dp) :: a, b
 integer :: i, j, k, l

  do k=1,frames
    l=0
    do j=1,ndim-1
      do i=1,ndim-j
        l=l+1
        a = dot_product(tangents(:,i,k),tangents(:,i+j,k))
        b = dacos(a)
        bends(k,l) = rad_to_deg*b !Convert to degrees
      end do
    end do
  end do

 end subroutine get_bendings_nofit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
 !This subroutine calculates the bending angles between tangent vectors
 ! The output is in degrees
 ! get_c_bends_bendings works for closed structures, each snapshot has
 ! n possible lengths and n parameters for each length. This is because
 ! the bending between i and j=i+l, goes from i=1 to i=n, for every
 ! length l.
 ! UPDATE: calculates the bending angle from 0 to 360. It calculates the
 ! bending angle of atoms projected to a plane. Bending have directions

 subroutine get_c_bendings(ndim,frames,G_n,best,tangents,bends)
 
 implicit none
 
 integer, intent(in) :: ndim,frames, best(:)
 real(dp), intent(in) :: G_n(:,:,:), tangents(:,:,:)
 real(dp), intent(out) :: bends(:,:,:)
 real(dp) :: a, b, d, c(3)
 integer :: i, k, l, j

 !l = length
 !i = initial bp

  do k=1,frames

    !-----------------------------
    !length l = 1 (first bendings)
    l = 1
    do i=1,ndim-l
      j=i+l
      !Define sign
      c = cross_product3( tangents(:,i,k), tangents(:,j,k) )
      d = dot_product( tangents(:,i,k),tangents(:,j,k) ) !auxiliar
      !The first two ifs are unlikely, but sometimes happen and
      !cause error...
      if ( d > 1.0_dp  ) then
        b = 0.0_dp
      else if ( d < -1.0_dp  ) then
        b = pi
      else
        a = dot_product(tangents(:,i,k),tangents(:,j,k))
        b = dacos(a)
      end if

      !Change sign if needed
      if ( dot_product( c, G_n(:,best(k),k) ) < 0.0_dp ) b = -b

      bends(k,i,l) = rad_to_deg*b
    end do !i

    !i+l is outside boundaries, this is fixed with another loop
    do i=ndim-l+1,ndim
      j=i+l-ndim
      !Define sign
      c = cross_product3( tangents(:,i,k), tangents(:,j,k) )
      d = dot_product( tangents(:,i,k),tangents(:,j,k) ) !auxiliar
      !The first two ifs are unlikely, but sometimes happen and
      !cause error...
      if ( d > 1.0_dp  ) then
        b = 0.0_dp
      else if ( d < -1.0_dp  ) then
        b = pi
      else
        a = dot_product(tangents(:,i,k),tangents(:,j,k))
        b = dacos(a)
      end if

      !Change sign if needed
      if ( dot_product( c, G_n(:,best(k),k) ) < 0.0_dp ) b = -b

      bends(k,i,l) = rad_to_deg*b
    end do !i

    !-----------------------------
    !For lengths greater than 1 bp, an auxiliar variable T_bend,
    !which is the total bend, will be used. This variable will help to
    !determina the correct bending. The way in which the bending here 
    !is calculated, is by performing the previous operations again, 
    !but it could be obtained through the sum of bendings
    do l=2,ndim-1

      do i=1,ndim-l
        j=i+l

        bends(k,i,l) = sum( bends(k,i:i+l-1,1) )

      end do !i

      do i=ndim-l+1,ndim
        j = i+l-ndim

        bends(k,i,l) = sum( bends(k,i:ndim,1) ) + sum( bends(k,1:j-1,1) )

      end do !i

    end do !l

  end do !k

  bends = abs(bends)

 end subroutine get_c_bendings
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
 !This subroutine calculates the bending angles between tangent vectors
 ! The output is in degrees
 ! get_c_bends_bendings works for closed structures, each snapshot has
 ! n possible lengths and n parameters for each length. This is because
 ! the bending between i and i+l, goes from i=1 to i=n, for every
 ! length l.
 ! This subroutine is for not projected trajectories, this means that
 ! bendings won't be greater than 180 degrees.
 subroutine get_c_bendings_nofit(ndim,frames,tangents,bends)

 implicit none

 integer, intent(in) :: ndim,frames
 real(dp), intent(in) :: tangents(:,:,:)
 real(dp), intent(out) :: bends(:,:,:)
 real(dp) :: a, b
 integer :: i, k, l

  do k=1,frames
    do l=1,ndim-1
      do i=1,ndim-l
        a = dot_product(tangents(:,i,k),tangents(:,i+l,k))
        b = dacos(a)
        bends(k,i,l) = rad_to_deg*b
      end do
      !i+l is outside boundaries, this is fixed with another loop
      do i=ndim-l+1,ndim
        a = dot_product(tangents(:,i,k),tangents(:,i+l-ndim,k))
        b = dacos(a)
        bends(k,i,l) = rad_to_deg*b
      end do
    end do
  end do

 end subroutine get_c_bendings_nofit
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
! This subroutines writes the projected structure where the
! normal vector to the plane is align with the zth axis

! Even though functions_mod performs all the mathematical functions,
! this subroutine is a mixed of io_mod and functions_mod...

  subroutine write_projected_coords(coords,G_n,best,seq,nbp,frames,xyz)
  real(dp), intent(in) :: coords(:,:,:), G_n(:,:,:)
  integer,  intent(in) :: best(:), nbp, frames
  logical, intent(in) :: xyz
  character(1), intent(in) :: seq(:)

  !Auxiliar variables
  real(dp) :: b(3), d(3), Rot(3,3), angle
  character(160) :: caux
  integer :: i, k, ierror

  !Decide name of file
  if (xyz) then
    open(unit=10, file="projected.xyz", status="replace",action="write", iostat=ierror)
    if (ierror/=0) stop "Error opening output projection trajectory"
  else
    open(unit=11, file="projected.crd", status="replace",action="write", iostat=ierror)
    if (ierror/=0) stop "Error opening output projection trajectory"

    write(11,"(A)") "Generated by SerraLINE" !Small message
  end if

  !Let's begin!

  do k = 1, frames

   !Let's see the angle between normal and z axis (e3 = z axis)
    b = G_n(:,best(k),k) !This is the normal vector
    angle = dacos( dot_product( b, e3 ) )

    !Now the perpendicular vector to those two
    d = cross_product3( b, e3)
    d = normalize_vector(d,3)

    !Create rotation matrix (Rot)
    call general_rotation_matrix(Rot,d,angle)

    !Now let's rotate and write it
    if (xyz) then                                   !XYZ FORMAT

       write(10,"(I10)", iostat=ierror) nbp
       write(10,*, iostat=ierror) "   "

       do i = 1, nbp
          write(caux,*) i
          caux = adjustl(caux)
          caux = seq(i)//trim(caux)
          write(10, trim(F_XYZ_T) ,iostat=ierror) trim(caux), matmul( Rot, coords(:,i,k) )
          if (ierror/=0) stop "Error printing projection trajectory"
       end do

    else                                             !CRD FORMAT

       write(11, trim(F_AMB_T), iostat=ierror) matmul( Rot, coords(:,:,k) )
       if (ierror/=0) stop "Error printing projection trajectory"

    end if !xyz

  end do !k

  !Close file
  if (xyz) then
    close(10, iostat=ierror)
    if (ierror/=0) stop "Error closing output projection trajectory"
  else
    close(11, iostat=ierror)
    if (ierror/=0) stop "Error closeing output projection trajectory"
  end if



  end subroutine write_projected_coords
!----------------------------------------------------------------------



!-----------------------------------------------------------------------
 !This subroutine fits a plane to a set of s points
 !It returns a normal vector G_n to the plane 
 subroutine get_G_plane( points, s, G_n  )
 implicit none
 real(dp), intent(in) :: points(:,:)
 integer, intent(in) :: s
 real(dp), intent(out) :: G_n(3)
 real(dp) :: At(3,s), n(3), B(s), C(3,3), I_C(3,3), D(3)
! integer :: i

 At(1,:) = points(1,:) !x components
 At(2,:) = points(2,:) !y components
 At(3,:) = 1.0_dp      !1...
 B = points(3,:)       !z components

 C = matmul(At,transpose(At))

 I_C = matrix_inversion_3x3(C)

 D = matmul(At,B)

 n = matmul(I_C,D)

 G_n = normalize_vector(n,3)

 end subroutine get_G_plane
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Calculates the inverse of 
 function matrix_inversion_3x3(A) result(I_A)
 implicit none
 real(dp) :: A(3,3), I_A(3,3), det_A

 !First calculate determinant det_A
 
 det_A = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
     & + A(1,3)*A(2,1)*A(3,2) - A(3,1)*A(2,2)*A(1,3) &
     & - A(3,2)*A(2,3)*A(1,1) - A(3,3)*A(2,1)*A(1,2)

 if (abs(det_A) < eps) stop "0 determinant of matrix A"
 
 I_A(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
 I_A(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
 I_A(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
 I_A(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
 I_A(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
 I_A(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
 I_A(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
 I_A(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
 I_A(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

 I_A = I_A/det_A
 
 end function matrix_inversion_3x3
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This function normalizes vector v with dimension n
 function normalize_vector(v,n)
 implicit none
 integer, intent(in) :: n
 real(dp), intent(in) :: v(n)
 real(dp) :: normalize_vector(n)
 integer :: i
 real(dp) :: v_sum

 v_sum=0._dp

 do i=1,n
  v_sum = v_sum + v(i)*v(i)
 end do

 v_sum=sqrt(v_sum)

 do i=1,n
  normalize_vector(i)=v(i)/v_sum
 end do

 end function 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Calculates cross product =  U X V
 function cross_product3(U,V)
 real(dp) :: U(3), V(3)
 real(dp) :: cross_product3(3)

 cross_product3(1) = U(2)*V(3)-U(3)*V(2)
 cross_product3(2) = U(3)*V(1)-U(1)*V(3)
 cross_product3(3) = U(1)*V(2)-U(2)*V(1)

 end function cross_product3
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This function calculates the projection of vector u on v.
 !Note that this subroutine assumes that u and v are unit vectors, this
 !knowing this the calculation will be speed up.
 function projv_u(u,v)
 implicit none
 real(dp) :: u(3), v(3), projv_u(3), a

 a = dot_product(u,v)

 projv_u = a*v/absv(v) !No division needed since absolute of v is 1

 end function projv_u

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This function calculates the projection of vector u on v.
 !Note that this subroutine assumes that u and v are unit vectors, this
 !knowing this the calculation will be speed up.
 function nprojv_u(u,v)
 implicit none
 real(dp) :: u(3), v(3), nprojv_u(3), a

 a = dot_product(u,v)

 nprojv_u = a*v !No division needed since absolute of v is 1

 end function nprojv_u

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !Average and standard deviations are stored in the same arrays in SerraLINE
 function average_std(X,n)
 implicit none
 real(dp) :: average_std(2)
 integer :: n
 real(dp) :: X(n)
 real(dp) :: a !auxiliar variable
 integer :: i

 !Get Average
 average_std(1) = sum(X)/real(n,dp)

 !Now Standard deviation
 average_std(2)=0.0_dp

 do i=1,n
   a = X(i) - average_std(1)
   average_std(2) = average_std(2) + a*a
 end do
 average_std(2) = average_std(2)/real(n,dp)
 average_std(2) = sqrt(average_std(2))

 end function 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculates the square value of vector a and stores it in b
 function absv2(a) result (b)
 implicit none
 real(dp), intent(in) :: a(:)
 real(dp) :: b
 integer :: i
  b=0.
  do i=1,size(a)
   b=b+a(i)**2
  end do
 end function absv2

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculates the absolute value of vector a and stores it in b
 function absv(a) result (b)
 implicit none
 real(dp), intent(in) :: a(:)
 real(dp) :: b
 integer :: i
  b=0.
  do i=1,size(a)
   b=b+a(i)**2
  end do
  b=sqrt(b)
 end function absv

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Constructs the general rotation matrix 'R' ('a')
! describes a rotation of magnitude 'a' about an
! unit vector 'u'
! (In R3)
 subroutine general_rotation_matrix(R,u,a)
 real(dp), intent(in) :: u(3), a
 real(dp), intent(out) :: R(3,3)

 R(1,1)=dcos(a)+(1-dcos(a))*u(1)**2
 R(1,2)=(1-dcos(a))*u(1)*u(2)-u(3)*dsin(a)
 R(1,3)=(1-dcos(a))*u(1)*u(3)+u(2)*dsin(a)
 R(2,1)=(1-dcos(a))*u(1)*u(2)+u(3)*dsin(a)
 R(2,2)=dcos(a)+(1-dcos(a))*u(2)**2
 R(2,3)=(1-dcos(a))*u(2)*u(3)-u(1)*dsin(a)
 R(3,1)=(1-dcos(a))*u(1)*u(3)-u(2)*dsin(a)
 R(3,2)=(1-dcos(a))*u(2)*u(3)+u(1)*dsin(a)
 R(3,3)=dcos(a)+(1-dcos(a))*u(3)**2

 end subroutine general_rotation_matrix
!-----------------------------------------------------------------------

 end module functions_mod
