!------------------------------------------------------------- 04/09/2021
! This module contains parameters needed for SerraLINE
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

  module parms
  implicit none

  integer, parameter :: sp= kind(0.0), dp = kind(0.d0), &
                        nlinR=9, nlinY=6    !purine and pyrimidine ring atoms

  real(dp), parameter :: pi = datan(1.0_dp)*2.0_dp*2.0_dp, &
                       & rad_to_deg = 180.0_dp/pi, eps = 10.0E-16_dp !enough small

  !BASES
!-----------------------------------------------------------------------

  character(len=4), dimension(6), parameter :: &

  !ADENINE
  & A_l = ["A   ", "A5  ", "A3  ", "DA  ", "DA5 ", "DA3 "], &

  !GUANINE
  & G_l = ["G   ", "G5  ", "G3  ", "DG  ", "DG5 ", "DG3 "], &

  !CYTOSINE
  & C_l = ["C   ", "C5  ", "C3  ", "DC  ", "DC5 ", "DC3 "], &

  !THYMINE
  & T_l = ["T   ", "T5  ", "T3  ", "DT  ", "DT5 ", "DT3 "], &

  !URACIL
  & U_l = ["U   ", "U5  ", "U3  ", "DU  ", "DU5 ", "DU3 "]
!-----------------------------------------------------------------------


  !RING ATOMS
!-----------------------------------------------------------------------
  character(len=4), parameter :: &

  N1 = "N1  ", C2 = "C2  ", N3 = "N3  ", C4 = "C4  ", C5 = "C5  ", &
  C6 = "C6  ", N7 = "N7  ", C8 = "C8  ", N9 = "N9  "


  !INPUT/OUTPUT FORMATS
!-----------------------------------------------------------------------
  character(50), parameter :: &

  !AMBER TRAJECTORY FORMAT
  & F_AMB_T = "(10F8.3)", &

  !OUTPUT PARAMETERS (including bending)
  & F_PARM_1 = "(A15,A20)", F_PARM_2 = "(I4,A1,I4,6A1,2F10.3)", &
  & F_PARM_3 = "(I10,A10)", F_PARM_4 = "(A15,2F10.3)", &
  & F_PARM_5 = "(A20,1I10)", &

  !BENDINGS (no widths, heights, ...)
  & F_AVB_1 = "(A15,A20)", F_AVB_2 = "(I4,A1,I4,6A1,2F10.3)", &
  & F_AVB_3 = "(I10,A10)",&

  !OUTPUT XYZ
  & F_XYZ_T = "(A5,3F9.3)"

  !UNITARY VECTORS
!----------------------------------------------------------------------- 

  real(dp), dimension(3), parameter :: &

  !X
  e1 = [1.0_dp, 0.0_dp, 0.0_dp], &  

  !Y
  e2 = [0.0_dp, 1.0_dp, 0.0_dp], &  

  !Z
  e3 = [0.0_dp, 0.0_dp, 1.0_dp]

! Everything public !!!!!!!
  private
  public sp, dp, nlinR, nlinY, &
       & pi, rad_to_deg, eps, &
       & A_l, G_l, C_l, T_l, U_l, N1, C2, N3, C4, C5, C6, N7, C8, N9, &
       & F_AMB_T, F_PARM_1, F_PARM_2, F_PARM_3, F_PARM_4, F_PARM_5, &
       & F_AVB_1, F_AVB_2, F_AVB_3, F_XYZ_T, &
       & e1, e2, e3

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!   contains

!-----------------------------------------------------------------------

  end module parms
