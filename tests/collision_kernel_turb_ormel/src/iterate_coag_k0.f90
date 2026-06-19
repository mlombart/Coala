!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: run coala coag for k=0
!-------------------------------------------

module iterate_coag_k0_DG
   use precision
   use phy_cst
   use setup, only: rhograin,kernel,isave,scut,coeff_pl,optimised_solver,coeff_CFL
   use GQLeg_nodes_weights
   use initial_condition_L2proj_power_law
   use functions_tabflux_tabintflux
   use functions_flux_intflux
   use compute_solver_coag
   implicit none
   
   contains
   
!> @brief Iterate coagulation solver to reach the time ndthydro x dthydro  
!!
!! DG scheme k=0, piecewise constant approximation
!!
!! @param[in]    path_data     path to save data
!! @param[in]    path_log      path to save logfile
!! @param[in]    chnbins       number of dust bins in string format
!! @param[in]    nbins         number of dust bins
!! @param[in]    kpol          degree of polynomials for approximation
!! @param[in]    sizegrid      grid of sizes, borders value of bins in size
!! @param[in]    massgrid      grid of masses, borders value of mass bins
!! @param[in]    massbins      arithmetic mean value of massgrid for each mass bins
!! @param[in]    Q             number of points for Gauss-Legendre quadrature
!! @param[in]    vecnodes      nodes of the Legendre polynomials
!! @param[in]    vecweights    weights coefficients for the Gauss-Legendre polynomials
!! @param[in]    ndthydro      number of hydro timestep
!! @param[in]    dthydro       hydro timestep
!! @param[in]    dv            2D array of the differential velocity between grains in bins
subroutine iterate_coag_k0(path_data,path_log,chnbins,nbins,kpol,sizegrid,massgrid,massbins,Q,vecnodes,vecweights,ndthydro,dthydro,dv)
   implicit none
   character(len=*),  intent(in) :: path_data
   character(len=100),intent(in) :: path_log
   character(len=3),  intent(in) :: chnbins
   integer,           intent(in) :: ndthydro,nbins,kpol,Q
   real(wp),          intent(in) :: dthydro
   real(wp),          intent(in) :: sizegrid(nbins+1),massgrid(nbins+1),massbins(nbins)
   real(wp),          intent(in) :: vecnodes(Q),vecweights(Q)
   real(wp),          intent(in) :: dv(nbins,nbins)

   integer  :: j,i,iprogress,tot_nsub,tot_ndt,ndt,nsub,bin_scut
   real(wp) :: gij(nbins),gijnew(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp) :: start,finish
   real(wp) :: time,M1_t0,M1_tend
   real(wp) :: K0
   real(wp) :: tabflux_coag(nbins,nbins,nbins)     
   
   character(len=100) :: path_gtend_massmeanlog,path_gt0_massmeanlog,&
                         path_gij_t0,path_gij_tend,path_time,path_gij

   path_gt0_massmeanlog   = trim(path_data)//"gt0_massmeanlog.txt"
   path_gtend_massmeanlog = trim(path_data)//"gtend_massmeanlog.txt"
   path_gij_t0            = trim(path_data)//"gij_t0.txt"
   path_gij_tend          = trim(path_data)//"gij_tend.txt"
   path_gij               = trim(path_data)//"gij.txt"
   path_time              = trim(path_data)//"time.txt"


   print*, ""
   print*, ''//achar(27)//'[96m>>>Precomputing arrays<<<'//achar(27)//'[0m'

   !generate Legendre polynomial coefficient
   call compute_mat_coeffs(kpol,mat_coeffs_leg)

   !precompute array for DG scheme k=0
   print*,"Computing arrays ..."

   !coeff to cross-section in mass to size^2
   K0 = pi*(4._wp/3._wp*pi*rhograin)**(-2._wp/3._wp)

   call cpu_time(start)
   call compute_coagtabflux_GQ_k0(kernel,K0,Q,vecnodes,vecweights,nbins,kpol,massgrid,mat_coeffs_leg,tabflux_coag)
   call cpu_time(finish)
   print '("Array generated in ",f10.3,"s")',finish-start


   !log simu add tabflux generated time
   if (isave) then
      open(unit=1,file=trim(path_log),action='write',position='append')
      write(1,'("Array generated in ",f10.3,"s")') finish-start
      close(unit=1)
   endif


   print*, ""
   print*, ''//achar(27)//'[96m>>>Time solver<<<'//achar(27)//'[0m'

   print*,"Computing initial power law initial distribution"
   ! power law initial distribution
   do j=1,nbins
      if ((sizegrid(j) < scut) .and. (scut <= sizegrid(j+1))) then
         bin_scut = j
      endif
   enddo
   !generate gij component on Legendre polynomials basis
   call cpu_time(start)
   call L2proj_powerlaw_GQ_k0(nbins,massgrid,massbins,Q,vecnodes,vecweights,coeff_pl,bin_scut,gij)
   call cpu_time(finish)


   print '("gij init in",f6.3,"s.")',finish-start
   print*,"gij init ="
   do j=1,nbins
      print*,gij(j)
   enddo

   call cpu_time(finish)
   print '("gij generated in ",f10.3,"s")',finish-start

   



   !write gij and gij massmeanlog start
   if (isave) then
      open(unit=1,file=trim(path_gij_t0))
      write(1,'('//chnbins//'(e30.16E3,2x))') gij
      close(unit=1)

      open(unit=1,file=trim(path_gt0_massmeanlog))
      write(1,'('//chnbins//'(e30.16E3,2x))') gij
      close(unit=1)
   endif


   !write gij and time for all time
   time=0._wp
   if (isave) then

      open(unit=1,file=trim(path_gij))
      write(1,'('//chnbins//'(e30.16E3,2x))') gij

      open(unit=2,file=trim(path_time))
      write(2,'(e30.16E3)') time

   endif

   !compute initial total dust mass density
   M1_t0 = 0._wp
   do j=1,nbins
      M1_t0 = M1_t0 + (massgrid(j+1)-massgrid(j))*gij(j)
   enddo
   print*,"Total dust mass density t0 =",M1_t0

   iprogress=1

   tot_nsub = 0
   tot_ndt  = 0

   !time solver for coagulation 
   print*,"Running coala ..."
   call cpu_time(start)
   do i=1,ndthydro

      if (optimised_solver) then
         call compute_coag_k0_kdv_opt(eps,coeff_CFL,nbins,massgrid,massbins,&
                                       tabflux_coag,dv,dthydro,gij,gijnew,nsub,ndt)
      else
         call compute_coag_k0_kdv(eps,coeff_CFL,nbins,massgrid,tabflux_coag,dv,dthydro,gij,gijnew,nsub,ndt)
      endif


      gij = gijnew
      tot_nsub = tot_nsub + nsub
      tot_ndt  = tot_ndt + ndt

      time = time + dthydro

      if (isave) then
         !write time, gij
         write(2,'(e30.16E3)') time
         write(1,'('//chnbins//'(e30.16E3,2x))') gij
      endif

      call display_progress_bar(ndthydro,iprogress)
      iprogress = iprogress+1

   enddo


   close(1)
   close(2)

   call cpu_time(finish)

   print*,""
   print*,""
   print '("Number of sub-timestep < dthydro for coagulation = ",i10)',tot_nsub
   print '("Number of timestep at dthydro = ",i10)',tot_ndt
   print '("Total number timesteps = ",i10)',tot_ndt+tot_nsub
   print '("Coala run in ",f10.3,"s")',finish-start
   print '("Time per timestep in ",f10.5,"ms")',(finish-start)/(tot_ndt+tot_nsub) * 1000.



   !log simu add Time solver time
   if (isave) then
      open(unit=1,file=trim(path_log),action='write',position='append')
      write(1,'("Number of sub-cyclings dt ",i10)') tot_nsub
      close(unit=1)

      open(unit=1,file=trim(path_log),action='write',position='append')
      write(1,'("Number of dt ",i10)') tot_ndt
      close(unit=1)

      open(unit=1,file=trim(path_log),action='write',position='append')
      write(1,'("Total number of time-step ",i10)') tot_ndt+tot_nsub
      close(unit=1)


      open(unit=1,file=trim(path_log),action='write',position='append')
      write(1,'("Coala run in ",f10.3,"s")') finish-start
      close(unit=1)

      open(unit=1,file=trim(path_log),action='write',position='append')
      write(1,'("Time per timestep in ",f10.5,"ms")') (finish-start)/(tot_ndt+tot_nsub) * 1000.
      close(unit=1)

   endif


   print*,""
   print*,"gij end"
   do j=1,nbins
      print*,gij(j)
   enddo


   M1_tend = 0._wp
   do j=1,nbins
      M1_tend = M1_tend + (massgrid(j+1)-massgrid(j))*gij(j)
   enddo
   print*,"Total dust mass density t0 =",M1_t0
   print*,"Total dust mass density tend =",M1_tend
   print*,"Abs error total dust mass density =",abs(M1_tend - M1_t0)/M1_t0



   !write gij end
   if (isave) then
      open(unit=1,file=trim(path_gij_tend))
      write(1,'('//chnbins//'(e30.16E3,2x))') gij
      close(unit=1)
   endif

   !write gmassmeanlog end
   if (isave) then
      open(unit=1,file=trim(path_gtend_massmeanlog))
      write(1,'('//chnbins//'(e30.16E3,2x))') gij
      close(unit=1)
   endif

   !log simu add mass
   if (isave) then
      open(unit=1,file=trim(path_log),action='write',position='append')
      write(1,'("Initial total mass ",e30.16E3,"s")') M1_t0
      write(1,'("Final total mass ",e30.16E3,"s")') M1_tend
      write(1,'("Abs error total dust mass density ",e30.16E3,"s")') abs(M1_tend - M1_t0)/M1_t0
      close(unit=1)
   endif

      

end subroutine iterate_coag_k0

end module iterate_coag_k0_DG
