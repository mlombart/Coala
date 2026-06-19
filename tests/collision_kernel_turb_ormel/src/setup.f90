!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: Setup parameters
! Coala parameters
! kernel           -> only k_cross_section
! K0               -> normalisation coefficient for physical kernel
! nbins            -> number of grain size bins
! kpol             -> order of polynomials
! Q                -> number of Gauss points for Gauss-quadrature method
! eps              -> minimum value for mass density distribution
! coeff_CFL        -> coefficient for SSPRK 3 time solver
! optimised_solver -> to use optimised version of solver in Coala
! isave            -> to save data for plots
!
! Physical parameters
! smax             -> maximum grainsize
! scut             -> initial largest grainsize 
! smin             -> minimum grainsize
! coeff_pl         -> coefficient pfor initial power-law size distribution
! rhograin         -> intrinsic grainsize density
! dtg              -> initial dust-to-gas ratio
! rho_gas          -> gas density
! temp             -> gas temperature
! alpha_turb       -> level of turbulence for grain-grain collision
! dynamical_time   -> time scale  
! dthydro          -> hydro timestep 
! n_tdyn           -> number of dynamical time
! ndthydro         -> number of dthydro
!-------------------------------------------
module setup
   use precision
   use phy_cst
   use ISO_FORTRAN_ENV
   implicit none


   !> @brief choose collision kernel
   !! k_cross_section -> ballistic collision kernel with 2D array dv from hydro or subgrid model
   character(len=*),  parameter :: kernel = "k_cross_section"

   integer,  parameter :: nbins = 20
   integer,  parameter :: kpol = 0
   integer,  parameter :: Q = 15
   real(wp), parameter :: eps = 1e-30_wp
   real(wp), parameter :: coeff_CFL = 3e-1_wp

   !choose optimised version for flux caculations and then solver in Coala
   logical :: optimised_solver = .true.

   !Save data
   logical :: isave = .false.


   !limit grainsizes in cm
   real(wp), parameter :: smax = 1._wp
   real(wp), parameter :: scut = 250e-7_wp
   real(wp), parameter :: smin = 5e-7_wp

   !power law for initial distribution in size
   real(wp), parameter :: coeff_pl = -3.5_wp !mrn

   !physical quantities
   real(wp), parameter :: rhograin    = 2.3_wp
   real(wp), parameter :: dtg         = 1e-2_wp
   real(wp), parameter :: rho_gas     = 1e-15_wp
   real(wp), parameter :: temp        = 10._wp
   real(wp), parameter :: alpha_turb  = 1.5_wp

   !time variables
   character(len=*),  parameter :: dynamical_time = "tff"
   real(wp),          parameter :: dthydro = 1e2_wp * yr 
   integer,           parameter :: n_tdyn = 5 ! in dynamical time unit
   integer                      :: ndthydro


   


end module setup