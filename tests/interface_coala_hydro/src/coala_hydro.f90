!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


program coala_hydro
   use precision
   use phy_cst
   use utils
   use init_massgrid
   use GQLeg_nodes_weights
   use dust_dv
   use functions_tabflux_tabintflux
   use functions_flux_intflux
   use interface_coala_coag

   implicit none

   integer           :: ndust,bin_scut,n_tff
   real(wp)          :: smax,scut,smin,coeff_pl,coeff_norm,hj,eps_rhodust
   real(wp)          :: rhograin,dtg,rho_gas,temp,alpha_turb,tff,dthydro
   character(len=30) :: kernel,source_dv

   real(wp), allocatable :: sizegrid(:),sizemeanlog(:)
   real(wp), allocatable :: dv(:,:),dv_ormel(:,:),dv_br(:,:)
   real(wp), allocatable :: rhodust(:),new_rhodust(:)

   real(wp) :: K0
   integer  :: kpol,Q
   logical  :: isave
   real(wp), allocatable :: massgrid(:),massbins(:),massmeanlog(:)
   real(wp), allocatable :: mat_coeffs_leg(:,:),vecnodes(:),vecweights(:)
   real(wp), allocatable :: tabflux_coag_k0(:,:,:)
   real(wp), allocatable :: tabflux_coag(:,:,:,:,:),tabintflux_coag(:,:,:,:,:,:)

   integer  :: j
   real(wp) :: start,finish,start_tot,finish_tot

   !for saving data
   character(len=100) :: path_data,path_massgrid,path_sizegrid,path_sizemeanlog,path_rhodust_t0,path_rhodust_tend,path_time
   character(len=10)  :: chndust,ch_dim_massgrid

   

   

   !-------------------------------------------
   ! physical parameters (cgs) defined in hydro code
   !-------------------------------------------
   ! ndust        -> number of dust size bins
   ! smax         -> maximum grainsize
   ! scut         -> initial largest grainsize 
   ! smin         -> minimum grainsize
   ! coeff_pl     -> coefficient for initial power-law size distribution
   ! rhograin     -> intrinsic grainsize density
   ! dtg          -> initial dust-to-gas ratio
   ! rho_gas      -> gas density
   ! temp         -> gas temperature
   ! alpha_turb   -> level of turbulence for grain-grain collision
   ! tff          -> free fall time

   ndust = 20

   smax        = 1._wp  
   scut        = 250e-7_wp
   smin        = 5e-7_wp
   coeff_pl    = -3.5_wp
   eps_rhodust = 1e-30_wp

   rhograin    = 2.3_wp
   dtg         = 1e-2_wp
   rho_gas     = 1e-15_wp
   temp        = 10._wp
   alpha_turb  = 1.5_wp

   tff     = sqrt(3._wp*pi/(32._wp*grav*rho_gas))
   n_tff   = 10
   dthydro = n_tff*tff


   !choice of grain-grain differential velocity
   !"dv_ormel"          -> dv from Ormel's turbulence model
   !"dv_brownian"       -> dv from Brownian motion in physical unit
   !"dv_ormel+dv_brownian" -> both sources
   source_dv = "dv_ormel+dv_brownian"

   !-------------------
   ! COALA parameters 
   !-------------------
   ! kernel           -> only "k_cross_section" for ballistic collision kernel with 2D array dv from hydro or subgrid model
   ! K0               -> normalisation coefficient for simple kernels
   ! kpol             -> order of polynomials ()
   ! Q                -> number of Gauss points for Gauss-quadrature method
   ! isave            -> to save data for plots

   kernel = "k_cross_section"
   kpol  = 0
   Q     = 15
   isave = .true.


   call cpu_time(start_tot)

   print*, ''//achar(27)//'[96m>>>Compute grid in size<<<'//achar(27)//'[0m'

   allocate(sizegrid(ndust+1),sizemeanlog(ndust))
   allocate(massgrid(ndust+1),massbins(ndust),massmeanlog(ndust))
   
   call cpu_time(start)
   call compute_sizegrid_massgrid(ndust,smax,smin,pi,rhograin,sizegrid,sizemeanlog,massgrid,massbins,massmeanlog)
   call cpu_time(finish)

   print '("Init sizegrid in ",f10.3,"s.")',finish-start


   !compute dv 2D array for k_cross_section
   print*,""
   print*, ''//achar(27)//'[96m>>>Compute dv 2D array<<<'//achar(27)//'[0m'
   
   allocate(dv(ndust,ndust))

   print*, "source dv -> ",source_dv

   call cpu_time(start)
   select case (source_dv)
   case ("dv_ormel")
      call compute_dv_ormel(ndust,rho_gas,temp,rhograin,sizemeanlog,tff,alpha_turb,dv)
      
   case ("dv_brownian")
      call compute_dv_brownian(ndust,temp,massmeanlog,dv)

   case ("dv_ormel+dv_brownian")
      allocate(dv_ormel(ndust,ndust))
      allocate(dv_br(ndust,ndust))

      call compute_dv_ormel(ndust,rho_gas,temp,rhograin,sizemeanlog,tff,alpha_turb,dv_ormel)
      call compute_dv_brownian(ndust,temp,massmeanlog,dv_br)
      dv = sqrt(dv_ormel**2 + dv_br**2)

      deallocate(dv_ormel)
      deallocate(dv_br)

   case default
      stop "Need to choose available kernel"

   end select

   call cpu_time(finish)
   print '(" dv array generated in ",f10.3,"s.")',finish-start

   
   !Precomputing part has be made once before time solver
   !Computing depend only on massgrid
   print*, ""
   print*, ''//achar(27)//'[96m>>>Precomputing arrays for COALA<<<'//achar(27)//'[0m'

   !generate Legendre polynomial coefficient
   allocate(mat_coeffs_leg(kpol+1,kpol+1))
   call compute_mat_coeffs(kpol,mat_coeffs_leg)



   !precompute array for DG scheme 
   print*,"Computing arrays ..."

   !for Gauss quadrature
   allocate(vecnodes(Q),vecweights(Q))
   call GQLeg_nodes(Q,vecnodes)
   call GQLeg_weights(Q,vecweights)

   !coeff to convert cross-section in mass into size^2
   K0 = pi*(4._wp/3._wp*pi*rhograin)**(-2._wp/3._wp)

   call cpu_time(start)

   if (kpol==0) then
      allocate(tabflux_coag_k0(ndust,ndust,ndust))
      call compute_coagtabflux_GQ_k0(kernel,K0,Q,vecnodes,vecweights,ndust,kpol,massgrid,mat_coeffs_leg,tabflux_coag_k0)

   else

      stop "Need to fix bug for kpol > 0"
      ! allocate(tabflux_coag(ndust,ndust,ndust,kpol+1,kpol+1))
      ! allocate(tabintflux_coag(ndust,kpol+1,ndust,ndust,kpol+1,kpol+1))
      ! call compute_coagtabflux_GQ(kernel,K0,Q,vecnodes,vecweights,ndust,kpol,massgrid,mat_coeffs_leg,tabflux_coag)
      ! call compute_coagtabintflux_GQ(kernel,K0,Q,vecnodes,vecweights,ndust,kpol,massgrid,mat_coeffs_leg,tabintflux_coag)
   
   endif

   call cpu_time(finish)
   print '("Array generated in ",f10.3,"s")',finish-start


   print*, ""
   print*, ''//achar(27)//'[96m>>>Rhodust initial power law distribution<<<'//achar(27)//'[0m'
   
   call cpu_time(start)
   ! power law initial distribution up to scut
   do j=1,ndust
      if ((sizegrid(j) < scut) .and. (scut <= sizegrid(j+1))) then
         bin_scut = j
      endif
   enddo

   allocate(rhodust(ndust))
   rhodust = eps_rhodust

   coeff_norm = dtg*rho_gas * ((4._wp+coeff_pl)/3._wp)/(massgrid(bin_scut+1)**((4._wp+coeff_pl)/3._wp) - massgrid(1)**((4._wp+coeff_pl)/3._wp))
   do j=1,bin_scut
      hj = massgrid(j+1) - massgrid(j)
      rhodust(j) = coeff_norm * hj * massbins(j)**((1._wp + coeff_pl)/3._wp)
   enddo

   call cpu_time(finish)
   print '("initial rhodust generated in ",f10.3,"s")',finish-start
   print*,"rhodust init"
   do j=1,ndust
      print*,"j=",j,", rhodust(j)=",rhodust(j)
   enddo


   print*, ""
   print*, ''//achar(27)//'[96m>>>COALA time solver<<<'//achar(27)//'[0m'   

   print*, "Running COALA for dthydro = ",n_tff," tff"
   !Run coala for 1 hydro timestep
   allocate(new_rhodust(ndust))
   new_rhodust = 0._wp

   call cpu_time(start)
   select case (kpol)
   case (0)

      call coala_coag_k0(ndust,massgrid,tabflux_coag_k0,rhodust,eps_rhodust,dv,dthydro,new_rhodust)
      
      
   case(1:10)
      stop "Need to fix bug for kpol > 0"
      ! call coala_coag(ndust,kpol,massgrid,massbins,mat_coeffs_leg,Q,vecnodes,vecweights,tabflux_coag,tabintflux_coag,&
      !                   rhodust,eps_rhodust,dv,dthydro,new_rhodust)

   case default
      print*,"Need kpol <= 10"
      stop
   end select
   call cpu_time(finish)

   print'("Run coala coag in ",f10.3,"s")',finish-start


   print*,"rhodust end"   
   do j=1,ndust
      print*,"j=",j,", rhodust(j)=",new_rhodust(j)
   enddo

   call cpu_time(finish_tot)
   print*, ""
   print*, ''//achar(27)//'[96m>>>Total computing time<<<'//achar(27)//'[0m'
   print '("Elapsed Time: ",f10.3,"s")',finish_tot-start_tot


   !save data
   if (isave) then

      call init_path_files(ndust,kpol,source_dv,path_data,chndust)

      path_massgrid     = trim(path_data)//"massgrid.txt"
      path_sizegrid     = trim(path_data)//"sizegrid.txt"
      path_sizemeanlog  = trim(path_data)//"sizemeanlog.txt"
      path_rhodust_t0   = trim(path_data)//"rhodust_t0.txt"
      path_rhodust_tend = trim(path_data)//"rhodust_tend.txt"
      path_time         = trim(path_data)//"time.txt"

      if (log10(real(ndust+1,wp))<1) then
         write( ch_dim_massgrid,'(i1)' ) ndust+1
      else if (1<=log10(real(ndust+1,wp)) .and. log10(real(ndust+1,wp))<2) then
         write( ch_dim_massgrid,'(i2)' ) ndust+1
      else if (2<=log10(real(ndust+1,wp)) .and. log10(real(ndust+1,wp))<3) then
         write( ch_dim_massgrid,'(i3)' ) ndust+1
      else if (3<=log10(real(ndust+1,wp)) .and. log10(real(ndust+1,wp))<4) then
         write( ch_dim_massgrid,'(i4)' ) ndust+1
      else
         print*,"Need to set proper definition to convert ndust to string in coala.f90"
         stop
      endif


      !write data

      open(unit=1,file=trim(path_massgrid))
      write(1,'('//ch_dim_massgrid//'(e30.16E3,2x))') massgrid
      close(unit=1)
   
      open(unit=1,file=trim(path_sizegrid))
      write(1,'('//ch_dim_massgrid//'(e30.16E3,2x))') sizegrid
      close(unit=1)

      open(unit=1,file=trim(path_sizemeanlog))
      write(1,'('//chndust//'(e30.16E3,2x))') sizemeanlog
      close(unit=1)

      open(unit=1,file=trim(path_rhodust_t0))
      write(1,'('//chndust//'(e30.16E3,2x))') rhodust
      close(unit=1)

      open(unit=1,file=trim(path_rhodust_tend))
      write(1,'('//chndust//'(e30.16E3,2x))') new_rhodust
      close(unit=1)

      open(unit=1,file=trim(path_time))
      write(1,'(e30.16E3)') tff
      write(1,'(e30.16E3)') dthydro
      close(unit=1)

   endif

   

end program coala_hydro
