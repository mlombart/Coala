


# COALA code

The COALA code is a high-order scheme based on the discontinuous Galerkin method to solve the Smoluchowski coagulation equation in its conservative form (see details in [Lombart & Laibe, 2021, MNRAS](https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.4298L/abstract)).

This new scheme allows to treat accurately dust coagulation in 3D hydrodynamics simulations.

Current version: 

- Only coagulation
- DG scheme with integrals calculated with Gauss-Legendre quadrature and polynomials approximation up to order 10
- any collision kernels can be used

Development in progress:
- Fragmentation process
- C and Python version of COALA

## Citing the code
<!-- [![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg)](https://arxiv.org/abs/2011.12298) -->

[![nasa ADS](https://img.shields.io/badge/nasa%20ADS-2021MNRAS.501.4298L-blue)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.4298L/abstract) [![arXiv](https://img.shields.io/badge/arXiv-2011.12298-red)](https://arxiv.org/abs/2011.12298)


# Documentation
Proper documentation is coming soon !
COALA is detailed in section [Code description](#code_coala).



# Code desciption <a name="code_coala"></a>
## Table of contents
1. [Compilation](#compilation)
2. [Source files of COALA](#src)
3. [COALA tests](#coala_tests)
    1. [Tests simple kernels](#tests_simple_kernels)
        1. [Setup file](#setup_tests_simple_kernels)
        2. [Run tests](#run_tests_simple_kernels)
        3. [Plots](#plots_simple_kernels)
    2. [Test Brownian collision kernel](#test_Br_kernel)
        1. [Setup file](#setup_test_Br_kernel)
        2. [Run test](#run_test_Br_kernel)
        3. [Plot](#plot_Br_kernel)
    3. [Test collision kernel from Ormel's model](#test_ormel_kernel)
        1. [Setup file](#setup_test_ormel_kernel)
        2. [Run test](#run_test_ormel_kernel)
        3. [Plot](#plot_ormel_kernel)

4. [Use COALA in your code](#use_coala)
    1. [To define in your setup file](#setup_file)
    2. [Generate ```massgrid``` and ```massbins```](#compute_massgrid)
    3. [Pre-calculation for time solver](#precomputing)
        1. [Scheme order 0](#order_0)
        2. [Scheme higher orders](#higher_orders)

    4. [Run COALA](#run_coala)
    	1. [Time solver order 0](#scheme_k0)
    	2. [Time solver for higher order](#scheme_high_k)

    5. [Test interface COALA and hydro codes](#coala_hydro)
    	1. [Setup file](#setup_test_coala_hydro)
        2. [Initial dust density distribution](#init_test_coala_hydro) 
        3. [Run test](#run_test_coala_hydro)
        4. [Plot](#plot_coala_hydro)


*** 

## Compilation <a name="compilation"></a>

Requirements: 
- Cmake version >= 3.10
- gfortran
- Python installation

The tests of the COALA code use CMake to compile. According to your operating system, it can be required to adapt the path to the fortran compiler in the files 
```tests/simple_collision_kernels/src/CMakelists.txt```
```tests/physical_kernel/src/CMakelists.txt``` 

such as
```cmake
set (CMAKE_Fortran_COMPILER /path/to/gfortran)
```

## Source files of COALA <a name="src"></a>
The sources of the code are located in ``` src/``` 


```\precision.f90```
: define the precision (double precision or other) for the code

```\GQ_legendre_nodes_weights.f90```
: nodes and weights for the Gauss-Legendre quadrature method

```\progress_bar.f90```
: provide a progress bar for the calculations

```\tension_module.f90```
: file for interpolation with splines under tension

```\init_massgrid.f90```
: generate mass grid in log scale from maximum and minimum values of dust masses or dust size

```\polynomials_legendre.f90```
: functions for Legendre polynomials approximation

```\limiter.f90```
: to calculate the gamma coefficient for scaling limiter (positivity preserved)

```\reconstruction_g.f90```
: functions to reconstruct the polynomials approximating the solution

```\collision_kernels.f90```
: file containing the kernel functions (any kernels can by added through a function)

```/coag_term/coagflux_functions_GQ.f90```
: coagulation flux evaluated with GQ method

```/coag_term/coagintflux_functions_GQ.f90```
: term integral of the coagulation flux evaluated with GQ method

```\functions_tabflux_tabintflux.f90```
: functions to pre-calculate the arrays needed for ```flux``` and ```intflux``` for DG scheme

```\functions_flux_intflux.f90```
: functions to compute ```flux``` and ```intflux``` for DG scheme

```\solver_DG.f90```
: functions for the time solver for DG scheme

```\compute_solver_coag.f90```
: functions to compute coagulation solver for 1 hydro time-step

```\initial_condition_L2proj_GQ.f90```
: file to interpolate dust mass density with splines under tension and to generate the polynomials coefficients by L2 projection 

```\interface_coag.f90```
: functions to interface COALA with other codes



## COALA tests <a name="coala_tests"></a>

The code COALA can be tested with two different tests. The first one is the comparaison with the analytical solutions for the Smoluchowski coagulation equation (constant, additive and multiplicative kernel). The second test is with a physical kernels with the Brownian motion and the Ormel's turbulence model as source of differential velocity between grains.


### Tests simple kernels <a name="tests_simple_kernels"></a>
The files are located in ```test/simple_collision_kernels```. The directory ```src``` is composed of 4 fortran files to run COALA


```setup.f90```
: file defining all the parameters to run COALA

```iterate_coag_k0.f90```
: function to iterate COALA with DG scheme order 0 for all hydro timesteps

```iterate_coag.f90```
: function to iterate COALA with DG scheme higher orders for all hydro timesteps

```coala.f90```
: main program to run COALA for tests with simple kernels

#### Setup file <a name="setup_tests_simple_kernels"></a>
The input parameters are defined in the setup file ```setup.f90``` and in the main file ```coala.f90```, here an example showing the parameters
```fortran
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

```

In ```setup.f90```, only the parameters ```kernel,K0,nbins,kpol,Q,eps,coeff_CFL,massmax,massmin``` are defined. The values of ```dthydro,ndthydro``` are defined according to the choice of kernel at the beginning of the main file ```coala.f90```:

```fortran
!set dthydro and ndthydro according to kernel
select case (kernel)
case ("kconst")
    dthydro = 1e2_wp
    ndthydro = 200

case ("kadd")
    dthydro = 1e-2_wp
    ndthydro = 300

case default
    print*,"Missing dthydro and ndthydro for kmul"
    stop
end select
```

#### Initial mass distribution for simple kernels <a name="init_tests_simple_kernels"></a>
The initial mass distribution for the simple kernels is $g(x,0) = x exp(-x)$ with $x$ the dimensional mass. Two functions to generate the initial mass distribution g are available in ```src/initial_condition_L2proj_GQ.f90```:
```fortran
//initial mass distribution for DG scheme order 0
L2proj_GQ_k0(eps,nbins,massgrid,massbins,Q,vecnodes,vecweights,gij)

//initial mass distribution for DG scheme higher orders
L2proj_GQ(eps,nbins,kpol,massgrid,massbins,Q,vecnodes,vecweights,gij)

```

#### Run test  <a name="run_tests_simple_kernels"></a>
To run the test, you first need to run Cmake:
```shell
cd test/simple_collision_kernels
cmake -S src/ -B build
```

Then you compile and run the code:
```shell
cd build/
make -j
./coala
```

#### Plots <a name="plots_simple_kernels"></a>
The python file for the plots are located in ```scripts_plots/```. For each test, you have a python file to generate the log-log scale plots of the mass density ditribution and the mass conservation. The variable `nbins` and `kernel` need to be specify in python file and then run the file such as:
```shell
python plots.py
```
By default plots are displayed but saving options are present in the file.


### Test Brownian collision kernel <a name="test_Br_kernel"></a>
The physical kernel is called the ballistic kernel and writes $K(x,y) = \sigma(x,y) \Delta v(x,y)$, with $\sigma(x,y)$ the cross-section and $\Delta v$ the differential velocity between grains. To include the ballistic kernel in the DG scheme, only the cross-section is integrated to calculate the flux and integral of the flux. The differential velocity is used in the DG scheme as a 2D histogram, then just implemented by multiplication to calculate flux and integral of the flux, in order to mimic a coupling between COALA and any hydro solver.

The files are located in ```test/collision_kernel_brownian```. The directory ```src``` is composed of 5 fortran files to run COALA

```setup.f90```
: file defining all the parameters to run COALA

```iterate_coag_k0.f90```
: function to iterate COALA with DG scheme order 0 for all hydro timestep

```iterate_coag.f90```
: function to iterate COALA with DG scheme higher orders for all hydro timestep

```dust_dv.f90```
: function giving the differential velocity between grains, here from the Brownian motion

```coala.f90```
: main program to run COALA for tests with the physical kernel


#### Setup file <a name="setup_test_Br_kernel"></a>
The input parameters are defined in the setup file ```setup.f90``` 
```fortran
!-------------------------------------------
! MODULE: Setup parameters
! kernel           -> collision kernel from Brownian motion (analytic or approximated)
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
   !! k_brownian      -> Brownian motion ballistic kernel (analytic expression)
   !! k_cross_section -> ballistic collision kernel with 2D array dv approximated from Brownian motion
   character(len=*),  parameter :: kernel = "k_cross_section"
   real(wp),          parameter :: K0 = 1._wp

   integer,  parameter :: nbins = 20
   integer,  parameter :: kpol = 0
   integer,  parameter :: Q = 5
   real(wp), parameter :: eps = 1e-30_wp
   real(wp), parameter :: coeff_CFL = 3e-1_wp

   !default setup brownian kernel
   real(wp) :: dthydro = 1e-1_wp
   integer  :: ndthydro = 500

   !limit mass grid
   real(wp) :: massmax = 1e6_wp
   real(wp) :: massmin = 1e-3_wp

   !choose optimised version for flux caculations and then solver
   logical :: optimised_solver = .true.

   !save data
   logical :: isave = .true.

```

#### Initial mass distribution <a name="init_test_Br_kernel"></a>
The initial mass distribution for the physical kernel is $g(x,0) = x exp(-x)$ with $x$ the dimensional mass, similar to the tests with simple kernels.

#### Run test  <a name="run_test_Br_kernel"></a>
To run the test, you first need to run Cmake:
```shell
cd test_coala/physical_kernel
cmake -S src/ -B build
```

Then you compile and run the code:
```shell
cd build/
make -j
./coala
```

#### Plot <a name="plot_Br_kernel"></a>
The python file for the plots are located in ```scripts_plot/```. Three plots are displayed, $g(x,t)$, $x.g(x,t)$ and the mass conservation.


### Test collision kernel from Ormel's turbulence model <a name="test_ormel_kernel"></a>
Similarly, only the geometric cross-section is integrated in Coala solver. The grain-grain differential velocity is given by Ormel's model and used as input parameter to calculate the collision kernel in the solver. 


The files are located in ```test/collision_kernel_turb_ormel```. The directory ```src``` is composed of 7 fortran files to run COALA

```setup.f90```
: file defining all the parameters to run COALA

```phy_cst.f90```
: file defining physical quantities

```dust_dv.f90```
: function giving the differential velocity between grains from Ormel's model

```initial_condition_L2proj_power_law.f90```
: function to compute the initial condition with power-law mass distribution

```iterate_coag_k0.f90```
: function to iterate COALA with DG scheme order 0 for all hydro timestep

```iterate_coag.f90```
: function to iterate COALA with DG scheme higher orders for all hydro timestep

```coala.f90```
: main program to run COALA for tests with the physical kernel


#### Setup file <a name="setup_test_ormel_kernel"></a>
The input parameters are defined in the setup file ```setup.f90``` 
```fortran
!-------------------------------------------
! MODULE: Setup parameters
! Coala parameters
! kernel           -> only k_cross_section
! K0               -> normalisation coefficient for simple kernels
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
   logical :: isave = .true.


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

```

#### Initial mass distribution <a name="init_test_ormel_kernel"></a>
The initial mass distribution is defined as a power-law in order to use MRN distribution.

#### Run test  <a name="run_test_ormel_kernel"></a>
To run the test, you first need to run Cmake:
```shell
cd test/collision_kernel_turb_ormel
cmake -S src/ -B build
```

Then you compile and run the code:
```shell
cd build/
make -j
./coala
```

#### Plot <a name="plot_ormel_kernel"></a>
The python file for the plots are located in ```scripts_plot/```. Three plots are displayed mass density distribution, mass fraction distribution and the mass conservation.



## Use COALA in your code <a name="use_coala"></a>
The COALA files in ``` src/``` has to be compiled in the following order:

```\precision.f90```
```\GQ_legendre_nodes_weights.f90```
```\tension_module.f90```
```\init_massgrid.f90```
```\polynomials_legendre.f90```
```\limiter.f90```
```\reconstruction_g.f90```
```\collision_kernels.f90```
```/coag_term/coagflux_functions_GQ.f90```
```/coag_term/coagintflux_functions_GQ.f90```
```\functions_tabflux_tabintflux.f90```
```\functions_flux_intflux.f90```
```\solver_DG.f90```
```\compute_solver_coag.f90```
```\initial_condition_L2proj_GQ.f90```
```\interface_coala_coag.f90```

In the following, steps are detailed to interface COALA to hydro code. An example is given in the file ```tests/interface_coala_hydro/src/coala_hydro```.


### To define in hydro code setup file <a name="setup_file"></a>

The first step is to define in hydro code setup file the input parameters for COALA
```fortran
character(len=30) :: kernel
integer           :: nbins
integer           :: kpol
integer           :: Q
double precision  :: rhodust
double precision  :: dv(nbins,nbins)
double precision  :: eps_rhodust
```

Example for values: 
```fortran
kernel = "k_cross_section"
nbins  = 20
kpol   = 0 ! up to 10
Q      = 5 ! 1<Q<20
```


### Generate ```massgrid``` and ```massbins``` <a name="compute_massgrid"></a>

It is assumed that the variables ```massgrid``` (the mass grid -> boundaries of mass bins) and ```massbins``` (arithmetic mean of the mass grid) can be provided by the hydro solver. 

If it not the case, the function ```compute_sizegrid_massgrid(nbins,smax,smin,pi,rhograin,sizegrid,sizemeanlog,massgrid,massbins,massmeanlog)``` in ```init_massgrid.f90``` can be used to compute ```massgrid, massbins``` by giving the following variables

```smax```
: maximum grain size considered in the simulation

```smin```
: minimum grain size considered in the simulation

```rhograin```
: intrinsic density of grains


### Precomputing for COALA solver <a name="precomputing"></a>

COALA is designed to be as fast as possible thanks to the precomputing part (functions depending only on the massgrid) for the time solver. The precomputing part is done only once for one simulation. This can be done in your dust_init file for instance.

In COALA, the file ```src/coala_GQ_legendre_nodes_weights.f90``` provide the nodes and weights for the Gauss-Legendre quadrature rule up to ```Q=20``` to evaluate integrals. It is required to generate the two variables ```vecnodes,vecweights``` such as
```fortran
double precision :: venodes(Q),vecweights(Q)

!for GQ
vecnodes = 0._wp
vecweights = 0._wp
call GQLeg_nodes(Q,vecnodes)
call GQLeg_weights(Q,vecweights)
```

#### COALA order 0 <a name="order_0"></a>
The precomputing part is to generate ```tabflux_coag``` needed to compute the coagulation flux for the DG scheme order 0. The function is located in ```function_tabflux_tabintflux.f90```.
```fortran
!precompute array for DG scheme
double precision :: tabflux_coag_k0(nbins,nbins,nbins)
double precision :: mat_coeffs_leg(kpol+1,kpol+1)

call compute_mat_coeffs(kpol,mat_coeffs_leg)

!coeff to convert cross-section in mass into size^2
K0 = pi*(4._wp/3._wp*pi*rhograin)**(-2._wp/3._wp)

call compute_coagtabflux_GQ_k0(kernel,K0,Q,vecnodes,vecweights,ndust,kpol,massgrid,mat_coeffs_leg,tabflux_coag_k0)
```

#### COALA higher-orders <a name="higher_orders"></a>

The precomputing part is to generate ```tabflux_coag, tabintflux_coag``` needed to compute coagulation flux and integral of flux for the DG scheme higher order (kpol>0). The functions are located in ```functions_tabflux_tabintflux.f90```.
```fortran
!precompute arrays for DG scheme
double precision :: tabflux_coag(nbins,nbins,nbins,kpol+1,kpol+1)
double precision :: tabintflux_coag(nbins,kpol+1,nbins,nbins,kpol+1,kpol+1)
double precision :: mat_coeffs_leg(kpol+1,kpol+1)

call compute_mat_coeffs(kpol,mat_coeffs_leg)


!coeff to convert cross-section in mass into size^2
K0 = pi*(4._wp/3._wp*pi*rhograin)**(-2._wp/3._wp)


call compute_coagtabflux_GQ(kernel,K0,Q,vecnodes,vecweights,ndust,kpol,massgrid,mat_coeffs_leg,tabflux_coag)
call compute_coagtabintflux_GQ(kernel,K0,Q,vecnodes,vecweights,ndust,kpol,massgrid,mat_coeffs_leg,tabintflux_coag)
   
```

### Run COALA <a name="run_coala"></a>

The functions ```coala_coag_k0``` and ```coala_coag``` located in ```src/interface_coala_coag.f90``` give the evolution of the dust mass density after a given hydrodynamique time-step.

Input requirement from hydro code: 
```rhodust```       => the dust mass density

```dv```            => the differential velocities between grains (dim (nbins,nbins))

```dthydro```       => the hydro time-step

```eps_rhodust```   => the lowest value considered for rhodust.


#### Coala solver order 0 <a name="coala_k0"></a>

In your file to make evolve the dust mass density after one hydro step, you just need to call the function ```coala_coag_k0``` such as
```fortran
call coala_coag_k0(ndust,massgrid,tabflux_coag_k0,rhodust,eps_rhodust,dv,dthydro,new_rhodust)
```
where ```new_rhodust``` is the evolved ```rhodust``` after coagulation process.

The dust mass density by mass unit needed in the Smoluchowski coagulation equation ```gij``` (the unknown variable) is calculated from rhodust in the function.

#### Coala solver higher orders <a name="coala_k>0"></a>

In your file to make evolve the dust mass density after one hydro step, you just need to call the function ```coala_coag``` such as

```fortran
call coala_coag(ndust,kpol,massgrid,massbins,mat_coeffs_leg,Q,vecnodes,vecweights,tabflux_coag,tabintflux_coag,rhodust,eps_rhodust,dv,dthydro,new_rhodust)
```
where ```new_rhodust``` is the evolved ```rhodust``` after coagulation process.

The continuous dust mass density by mass unit needed in the Smoluchowski coagulation equation (the unknown variable) is calculated in the function from rhodust by interpolation with splines in tension. Then the interpolation is projected by L2-norm on Legendre polynomials basis with Gauss-Legendre quadrature method to obtain the polynomials component ```gij```.

#### Test interface COALA and hydro codes <a name="coala_hydro"></a>

##### Setup file <a name="setup_test_coala_hydro"></a>
The physical parameters and the one needed for COALA solver are defined in the main file ```coala_hydro.f90```. All parameters defined above are used in this file.

In this test, you can use among three sources of grain-grain differential velocity: Ormel's model, Brownian motion or both.

```fortran
!choice of grain-grain differential velocity
!"dv_ormel"          -> dv from Ormel's turbulence model
!"dv_brownian"       -> dv from Brownian motion in physical unit
!"dv_ormel+dv_brownian" -> both sources
source_dv = "dv_ormel+dv_brownian"
```

#### Initial dust density distribution <a name="init_test_coala_hydro"></a>
The initial mass distribution is defined as a power-law in order to use MRN distribution.

#### Run test  <a name="run_test_coala_hydro"></a>
To run the test, you first need to run Cmake:
```shell
cd test/interface_coala_hydro
cmake -S src/ -B build
```

Then you compile and run the code:
```shell
cd build/
make -j
./coala_hydro
```

#### Plot <a name="plot_coala_hydro"></a>
The python file for the plots are located in ```scripts_plot/```. Three plots are displayed mass density distribution, mass fraction distribution and the mass conservation.



# License
COALA is licensed under the `CeCILL Free Software License Agreement v2.1`.

Copyright 2026 Maxime Lombart

