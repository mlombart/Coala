!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: L2 projection on Legendre polynomials basis
!         for initial condition with GQ method
!-------------------------------------------
module initial_condition_L2proj_power_law
   use precision
   use setup, only: dtg,rho_gas,eps
   use phy_cst
   use polynomials_legendre
   implicit none

contains


subroutine L2proj_powerlaw_GQ_k0(nbins,massgrid,massbins,Q,vecnodes,vecweights,coeff_pl,bin_scut,gij)
   implicit none
   integer,  intent(in)  :: nbins,Q,bin_scut
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins)
   real(wp), intent(in)  :: coeff_pl
   real(wp), intent(in)  :: vecnodes(Q),vecweights(Q)
   real(wp), intent(out) :: gij(nbins)

   real(wp) :: c,xj,hj,term_sum,xjalpha,node,coeff_norm
   integer  :: alpha,j


   coeff_norm = dtg*rho_gas * ((4._wp+coeff_pl)/3._wp)/(massgrid(bin_scut+1)**((4._wp+coeff_pl)/3._wp) - massgrid(1)**((4._wp+coeff_pl)/3._wp))

   gij = eps

   do j=1,bin_scut
   
      xj = massbins(j)
      hj = massgrid(j+1)-massgrid(j)

      call coeffnorm(0,c)
      term_sum = 0._wp

      do alpha=1,Q
         node = vecnodes(alpha)
         xjalpha = xj + hj*node/2._wp

         term_sum = term_sum + vecweights(alpha)*xjalpha**((1._wp + coeff_pl)/3._wp)
      enddo

      gij(j)=coeff_norm*term_sum/c

   enddo


end subroutine L2proj_powerlaw_GQ_k0




subroutine L2proj_powerlaw_GQ(nbins,kpol,massgrid,massbins,mat_coeffs_leg,Q,vecnodes,vecweights,coeff_pl,bin_scut,gij)
   implicit none
   integer,  intent(in)  :: nbins,kpol,Q,bin_scut
   real(wp), intent(in)  :: massgrid(nbins+1),massbins(nbins),mat_coeffs_leg(kpol+1,kpol+1)
   real(wp), intent(in)  :: coeff_pl
   real(wp), intent(in)  :: vecnodes(Q),vecweights(Q)
   real(wp), intent(out) :: gij(nbins,kpol+1) 

   real(wp) :: c,xj,hj,term_sum,xjalpha,node,coeff_norm
   integer :: alpha,k,j
   real(wp), allocatable :: ak(:)

   coeff_norm = dtg*rho_gas * ((4._wp+coeff_pl)/3._wp)/(massgrid(bin_scut+1)**((4._wp+coeff_pl)/3._wp) - massgrid(1)**((4._wp+coeff_pl)/3._wp))

   gij(:,1) = eps
   gij(:,2:) = 0._wp

   do j=1,bin_scut
      xj = massbins(j)
      hj = massgrid(j+1)-massgrid(j)

      do k=0,kpol
         allocate(ak(k+1))
         ak = mat_coeffs_leg(k+1,:k+1)

         call coeffnorm(0,c)
         term_sum = 0._wp

         do alpha=1,Q
            node = vecnodes(alpha)
            xjalpha = xj + hj*node/2._wp

            term_sum = term_sum + vecweights(alpha)*xjalpha**((1._wp + coeff_pl)/3._wp) * phi_pol(k,ak,vecnodes(alpha))
         enddo

         gij(j,k+1) = coeff_norm*term_sum/c

         deallocate(ak)

      enddo

   enddo

end subroutine L2proj_powerlaw_GQ



end module initial_condition_L2proj_power_law
