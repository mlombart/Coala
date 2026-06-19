cmake_minimum_required(VERSION 3.10)


project(coala)
enable_language (Fortran)
option(Coala_DEBUG "Enable Coala debug features" ON)

if(Coala_DEBUG)
   set(CMAKE_Fortran_FLAGS "-g -ffree-line-length-none -pedantic -fcheck=all -fbacktrace -std=f2008 -Werror ")
else()
   set(CMAKE_Fortran_FLAGS "-O3 -ffree-line-length-none -std=f2008 ")
endif()

message("flag -> ${CMAKE_Fortran_FLAGS}")


set(DIR_COALA ../../../src)
include_directories(${DIR_COALA})

set(SOURCES
   ${DIR_COALA}/precision.f90
   ${DIR_COALA}/progress_bar.f90
   ${DIR_COALA}/tension_module.f90
   ${DIR_COALA}/init_massgrid.f90
   ${DIR_COALA}/polynomials_legendre.f90
   ${DIR_COALA}/limiter.f90
   ${DIR_COALA}/reconstruction_g.f90
   ${DIR_COALA}/GQ_legendre_nodes_weights.f90
   ${DIR_COALA}/coag_term/coagflux_functions_GQ.f90
   ${DIR_COALA}/coag_term/coagintflux_functions_GQ.f90
   ${DIR_COALA}/collision_kernels.f90
   ${DIR_COALA}/functions_tabflux_tabintflux.f90
   ${DIR_COALA}/functions_flux_intflux.f90
   ${DIR_COALA}/solver_DG.f90
   ${DIR_COALA}/initial_condition_L2proj_GQ.f90
   ${DIR_COALA}/compute_solver_coag.f90
   ${DIR_COALA}/interface_coala_coag.f90
   ../../utils.f90
   phy_cst.f90
   dust_dv.f90
   coala_hydro.f90
)
add_executable(coala_hydro ${SOURCES})




