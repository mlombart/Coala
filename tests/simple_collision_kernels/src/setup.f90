!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: Setup parameters
! kernel           -> choose among simple kernels
! K0               -> normalisation coefficient for simple kernels
! nbins            -> number of grain size bins
! kpol             -> order of polynomials
! Q                -> number of Gauss points for Gauss-quadrature method
! eps              -> minimum value for mass density distribution
! coeff_CFL        -> coefficient for SSPRK 3 time solver
! dthydro          -> hydro timestep 
! ndthydro         -> number of dthydro
! massmax          -> maximal mass of grains to consider
! minmass          -> minimal mass of grains to consider
! optimised_solver -> to use optimised version of solver in Coala
! isave            -> to save data for plots
!-------------------------------------------
module setup
   use precision
   implicit none


   !> @brief choose collision kernel
   !! among kconst, kadd, kmul
   character(len=*),  parameter :: kernel = "kadd"
   real(wp),          parameter :: K0 = 1._wp

   integer,  parameter :: nbins = 20
   integer,  parameter :: kpol = 0
   integer,  parameter :: Q = 5
   real(wp), parameter :: eps = 1e-30_wp
   real(wp), parameter :: coeff_CFL = 3e-1_wp

   !default setup simple kernels (value defined in main program coala.f90)
   real(wp) :: dthydro
   integer  :: ndthydro

   !limit mass grid 
   real(wp) :: massmax = 1e6_wp
   real(wp) :: massmin = 1e-3_wp

   !choose optimised version for flux caculations and then solver
   logical :: optimised_solver = .true.


   !save data
   logical :: isave = .true.



end module setup