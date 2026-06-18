!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: physical constant in cgs
! yr        -> year in second
! u         -> standard atomic weight
! mu_gas    -> mean molecular weight
! mh        -> hydrogen atomic mass 
! gamma_gas -> adiabatic index
! grav      -> gravitational constant
! kB        -> Bolztamn constant
!------------------------------------------- 
module phy_cst
   use precision
   implicit none

   real(wp), parameter :: pi        = 4._wp*atan(1._wp)
   real(wp), parameter :: yr        = 31556926._wp
   real(wp), parameter :: u         = 1660538921e-33_wp
   real(wp), parameter :: mu_gas    = 23e-1_wp
   real(wp), parameter :: mh        = 100749e-5_wp*1660538921e-33_wp
   real(wp), parameter :: gamma_gas = 5._wp/3._wp
   real(wp), parameter :: grav      = 6.67e-8_wp
   real(wp), parameter :: kB        = 1.38e-16_wp

end module phy_cst