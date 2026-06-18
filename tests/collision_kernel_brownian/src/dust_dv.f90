!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: functions for grain-grain differential velocity
!-------------------------------------------

module dust_dv
   use precision
   implicit none

contains

!> @brief Compute 2D array grain-grain differential velocity 
!!        for dimensionless Brownian collision kernel
!!
!! function depending on mass and function evaluated at geometric mean of massgrid
!!
!! @param[in]    nbins         number of dust bins
!! @param[in]    massmeanlog   geometric mean value of massgrid for each mass bins
!! @param[out]   dv            2D array of the differential velocity between grains in bins
subroutine dv_brownian(nbins,massmeanlog,dv)
   implicit none 
   integer,  intent(in)  :: nbins
   real(wp), intent(in)  :: massmeanlog(nbins)
   real(wp), intent(out) :: dv(nbins,nbins)

   integer :: i,j

   do i=1,nbins
      do j=1,nbins
         dv(i,j) = sqrt(1._wp/massmeanlog(i) + 1._wp/massmeanlog(j))
      enddo
   enddo

end subroutine dv_brownian



end module dust_dv