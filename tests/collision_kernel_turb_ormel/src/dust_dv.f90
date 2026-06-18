!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: grain-grain differential velocity
!         from Ormel's model, see Lebreuilly et al., 2023
!-------------------------------------------

module dust_dv
   use precision
   use setup, only:rho_gas,temp,rhograin
   implicit none

contains




!> @brief Compute 2D array grain-grain differential velocity from Ormel's model for turbulence
!!
!!
!! @param[in]    nbins         number of dust bins
!! @param[in]    sizemeanlog   geometric mean value of sizegrid for each size bins
!! @param[in]    t_dyn         dynamical time for dust localised spatially
!! @param[in]    alpha_turb    coeff for the level of turbulence
!! @param[out]   dv            2D array of the differential velocity between grains in bins
subroutine dv_ormel(nbins,sizemeanlog,t_dyn,alpha_turb,dv)
   use phy_cst
   implicit none 
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: sizemeanlog(nbins),t_dyn,alpha_turb
   real(wp), intent(out) :: dv(nbins,nbins)

   real(wp) :: nh,cs,Re,t_eta
   real(wp) :: ts_i,ts_j,ts_1,St_1,St_2,x_St,beta_St
   real(wp) :: res
   integer  :: i,j

   nh = rho_gas/(mu_gas*mh)
   cs = sqrt(gamma_gas*kB*temp/(mu_gas*mh))        !sound speed 
   
   Re = 62e6_wp*sqrt(nh/1e5_wp)*sqrt(temp/10._wp)  !Reynolds number
   t_eta = t_dyn/sqrt(Re)

   res = 0._wp
   dv  = 0._wp
   do i=1,nbins
      do j=1,nbins

         !stopping times
         ts_i = sqrt(pi*gamma_gas/8._wp) * rhograin*sizemeanlog(i)/(rho_gas*cs)
         ts_j = sqrt(pi*gamma_gas/8._wp) * rhograin*sizemeanlog(j)/(rho_gas*cs)
         ts_1 = ts_i

         !stokes numbers
         St_1    = ts_i/t_dyn
         St_2    = ts_j/t_dyn
         

         !to symmetrize dv
         if (j > i) then
            ts_1    = ts_j
            St_1    = ts_j/t_dyn
            St_2    = ts_i/t_dyn
         endif

         x_St    = St_2/St_1
         
         beta_St = 3.2_wp - (1._wp + x_St) + 2._wp/(1._wp + x_St) * (1._wp/2.6_wp + x_St**3/(1.6_wp + x_St))
         

         if (ts_1 < t_eta) then

            if (abs(St_1 - St_2) < epsilon(St_1)) then
               res = 0._wp
            else
               res = alpha_turb * cs**2 * (St_1 - St_2)/(St_1 + St_2) * (St_1**2/(St_1 + 1._wp/sqrt(Re)) + St_2**2/(St_2 + 1._wp/sqrt(Re)))
            endif

            ! res = alpha_turb * cs**2 * (St_1 - St_2)/(St_1 + St_2) * (St_1**2/(St_1 + 1._wp/sqrt(Re)) + St_2**2/(St_2 + 1._wp/sqrt(Re)))

         else if ( (t_eta <= ts_1) .and. (ts_1 < t_dyn) ) then
               
            res = alpha_turb * cs**2 * beta_St * St_1

         else
            res = alpha_turb * cs**2 * (1._wp/(St_1 + 1._wp) + 1._wp/(St_2 + 1._wp))

         endif

         dv(i,j) = sqrt(res)

      enddo
   enddo

end subroutine dv_ormel



end module dust_dv