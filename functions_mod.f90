! Module containing all the mathematical functions needed for      22/07/2019
! SerraLINE
 module functions_mod

 use parms
 use IO_mod

 implicit none

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 contains

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
 !This subroutine calculates the bending angles between tangent vectors
 ! The output is in degrees
 ! This subroutine is for a opened structure, that means that it will have
 ! n-1 possible lengths and for each length it have n(n-1)/2 parameters
 subroutine get_bendings(ndim,frames,tangents,bends)

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

 end subroutine get_bendings
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This subroutine calculates the bending angles between tangent vectors
 ! The output is in degrees
 ! get_c_bends_bendings works for closed structures, each snapshot has
 ! n possible lengths and n parameters for each length. This is because
 ! the bending between i and i+l, goes from i=1 to i=n, for every
 ! length l.
 subroutine get_c_bendings(ndim,frames,tangents,bends)

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

 end subroutine get_c_bendings
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
