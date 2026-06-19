!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


program coala
   use precision
   use phy_cst
   use utils
   use setup
   use init_massgrid
   use reconstruction_g
   use GQLeg_nodes_weights
   use iterate_coag_k0_DG
   use iterate_coag_DG
   use dust_dv

   implicit none

   real(wp) :: sizegrid(nbins+1),sizemeanlog(nbins)
   real(wp) :: massgrid(nbins+1),massbins(nbins),massmeanlog(nbins)
   real(wp) :: vecnodes(Q),vecweights(Q)
   real(wp) :: dv(nbins,nbins),t_dyn

   character(len=100) :: path_data,path_massgrid,path_massbins,path_massmeanlog,path_sizegrid,path_sizemeanlog,path_log
   character(len=10)  :: chnbins,ch_dim_massgrid

   real(wp) :: start,finish,start_tot,finish_tot


   !set ndthydro according to dynamical time
   select case (dynamical_time)
   case ("tff")
      !free fall time
      t_dyn = sqrt(3._wp*pi/(32._wp*grav*rho_gas))
      ndthydro = floor(n_tdyn*t_dyn/dthydro)

      if (ndthydro == 0) then
         print*,"t_dyn [yr] =",t_dyn/yr," n_tdyn=",n_tdyn,",dthydro [yr] =",dthydro/yr
         stop "Error in dthydro value"
      endif
   case default
      print*,"Need to add other dynamical time in cgs"
      stop

   end select


   if (log10(real(nbins+1,wp))<1) then
      write( ch_dim_massgrid,'(i1)' ) nbins+1
   else if (1<=log10(real(nbins+1,wp)) .and. log10(real(nbins+1,wp))<2) then
      write( ch_dim_massgrid,'(i2)' ) nbins+1
   else if (2<=log10(real(nbins+1,wp)) .and. log10(real(nbins+1,wp))<3) then
      write( ch_dim_massgrid,'(i3)' ) nbins+1
   else if (3<=log10(real(nbins+1,wp)) .and. log10(real(nbins+1,wp))<4) then
      write( ch_dim_massgrid,'(i4)' ) nbins+1
   else
      print*,"Need to set proper definition to convert nbins to string in coala.f90"
      stop
   endif



   call init_path_files(nbins,kpol,kernel,path_data,chnbins)

   path_massgrid    = trim(path_data)//"massgrid.txt"
   path_massbins    = trim(path_data)//"massbins.txt"
   path_massmeanlog = trim(path_data)//"massmeanlog.txt"
   path_sizegrid    = trim(path_data)//"sizegrid.txt"
   path_sizemeanlog = trim(path_data)//"sizemeanlog.txt"
   path_log         = trim(path_data)//"log.txt"
   
   
   print*, ''//achar(27)//'[96m>>>Setup Ormel collision kernel<<<'//achar(27)//'[0m'


   !check correct name for kernel
   select case (kernel)
   case ("k_cross_section")
      print*, "collision kernel -> ",kernel

   case default
      print*,"Need to choose correct kernel name in setup.f90"
      stop
      
   end select

   if (wp == real32) then
      print*, 'Valuetype = FP32'
   else if (wp == real64) then
      print*, 'Valuetype = FP64'
   else if (wp == real128) then
      print*, 'Valuetype = FP128'
   else
      print*, 'Error in type of wp'
      stop
   endif

   print*, "nbins = ",nbins
   print*, "kpol=",kpol
   print*, "smin=",smin
   print*, "smax=",smax
   print*, "dthydro = ",dthydro/yr," [yr]"
   print*, "ndthydro =",ndthydro
   print*, "tdyn ->"//trim(dynamical_time)//" = ",t_dyn/yr," [yr]"
   print*, "coeff CFL=",coeff_CFL
   print*, "Q=",Q
   print*, "eps = ",eps
   print*, "Gij and gxmeanlog saving: ",isave

   
   !log simu
   if (isave) then
      open(unit=1,file=trim(path_log))

      write(1,'("Kernel = ",a50)') 'cross section with approximated dv for Ormel model'


      if (wp == real32) then
         write(1,'("Precision DG scheme = ",a2)') 'sp'
      else if (wp == real64) then
         write(1,'("Precision DG scheme = ",a2)') 'dp'
      else if (wp == real128) then
         write(1,'("Precision DG scheme = ",a2)') 'qp'
      else
         print*, 'Error in type of wp'
         stop
      endif

      write(1,'("nbins = ",i3)') nbins
      write(1,'("kpol = ",i1)') kpol
      write(1,'("smin = ",es10.3E2)') smin
      write(1,'("smax = ",es10.3E2)') smax
      write(1,'("dthydro = ",es10.3E2," [yr]")') dthydro/yr
      write(1,'("ndthydro = ",i10)') ndthydro
      write(1,'("tdyn -> ",a5," = ",es10.3E2,"[yr]")') dynamical_time,t_dyn/yr
      write(1,'("coeff CFL = ",f10.3)') coeff_CFL
      write(1,'("nb Gauss points = ",i2)') Q
      write(1,'("eps  = ",es10.3E3)') eps
      close(unit=1)
   endif

   !for GQ
   vecnodes = 0._wp
   vecweights = 0._wp
   call GQLeg_nodes(Q,vecnodes)
   call GQLeg_weights(Q,vecweights)


   
   print*, ""
   print*, ''//achar(27)//'[96m>>>Compute grid in mass<<<'//achar(27)//'[0m'
   
   call cpu_time(start_tot)
   call cpu_time(start)
   call compute_sizegrid_massgrid(nbins,smax,smin,pi,rhograin,sizegrid,sizemeanlog,massgrid,massbins,massmeanlog)
   call cpu_time(finish)

   print '("Init massgrid and massbins in ",f10.3,"s.")',finish-start

   ! print*,"sizegrid=",sizegrid
   ! print*,"massgrid=",massgrid
   ! stop

   

   !write massgrid and sizegrid
   if (isave) then
      open(unit=1,file=trim(path_massgrid))
      write(1,'('//ch_dim_massgrid//'(e30.16E3,2x))') massgrid
      close(unit=1)

      open(unit=1,file=trim(path_massbins))
      write(1,'('//chnbins//'(e30.16E3,2x))') massbins
      close(unit=1)

      open(unit=1,file=trim(path_massmeanlog))
      write(1,'('//chnbins//'(e30.16E3,2x))') massmeanlog
      close(unit=1)

      open(unit=1,file=trim(path_sizegrid))
      write(1,'('//ch_dim_massgrid//'(e30.16E3,2x))') sizegrid
      close(unit=1)

      open(unit=1,file=trim(path_sizemeanlog))
      write(1,'('//chnbins//'(e30.16E3,2x))') sizemeanlog
      close(unit=1)

   endif

   !compute dv 2D array for k_cross_section
   print*, ''//achar(27)//'[96m>>>Compute dv 2D array<<<'//achar(27)//'[0m'
   call cpu_time(start)
   call dv_ormel(nbins,sizemeanlog,t_dyn,alpha_turb,dv)
   call cpu_time(finish)
   print '(" dv array generated in ",f10.3,"s.")',finish-start


   !Run coala with kpol <= 10
   select case (kpol)
   case (0)

      call cpu_time(start)
      call iterate_coag_k0(path_data,path_log,chnbins,nbins,kpol,sizegrid,massgrid,massbins,Q,vecnodes,vecweights,ndthydro,dthydro,dv)
      call cpu_time(finish)
      
   case(1:10)

      call cpu_time(start)
      call iterate_coag(path_data,path_log,chnbins,nbins,kpol,sizegrid,massgrid,massbins,massmeanlog,Q,vecnodes,vecweights,ndthydro,dthydro,dv)
      call cpu_time(finish)

            
   case default
      print*,"Need kpol <= 10"
      stop
   end select

   print'("Run coala coag in ",f10.3,"s")',finish-start

   !log simu add run coala time
   if (isave) then
      open(unit=1,file=trim(path_log),action='write',position='append')
      write(1,'("Run coala coag in ",f10.3,"s")') finish-start
      close(unit=1)
   endif
   
   call cpu_time(finish_tot)
   print*, ""
   print*, ''//achar(27)//'[96m>>>Total computing time<<<'//achar(27)//'[0m'
   print '("Elapsed Time: ",f10.3,"s")',finish_tot-start_tot

   !log simu add total time
   if (isave) then
      open(unit=1,file=trim(path_log),action='write',position='append')
      write(1,'("Elapsed Time: ",f10.3,"s")') finish_tot-start_tot
      close(unit=1)
   endif
   

end program coala
