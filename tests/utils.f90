!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************
module utils
    use precision
    contains

    !> @brief Compute path for data
    !!
    !!
    !! @param[in]    nbins        number of dust bins
    !! @param[in]    kernel       integer to choose collision kernel
    !! @param[out]   path_data    path to save data
    !! @param[out]   chnbins      string for nbins
    subroutine init_path_files(nbins,kpol,kernel,path_data,chnbins)
        implicit none
        integer,            intent(in)  :: nbins,kpol
        character(len=*),   intent(in)  :: kernel
        character(len=100), intent(out) :: path_data
        character(len=10),  intent(out) :: chnbins

        character(len=3)  :: chkpol
        integer :: i


        write( chkpol,'(i1)' ) kpol

        if (log10(real(nbins,wp))<1) then
            write( chnbins,'(i1)' ) nbins
        else if (1<=log10(real(nbins,wp)) .and. log10(real(nbins,wp))<2) then
            write( chnbins,'(i2)' ) nbins
        else if (2<=log10(real(nbins,wp)) .and. log10(real(nbins,wp))<3) then
            write( chnbins,'(i3)' ) nbins
        else if (3<=log10(real(nbins,wp)) .and. log10(real(nbins,wp))<4) then
            write( chnbins,'(i4)' ) nbins
        else
            print*,"Need to set proper definition to convert nbins to string in setup.f90"
            stop
        endif


        path_data = "../data/"//trim(kernel)//"/nbins="//trim(chnbins)//"/kpol="//trim(chkpol)//"/"
        call execute_command_line ("[ ! -f "//path_data//" ] && mkdir -p "//path_data, exitstat=i)

    end subroutine init_path_files

end module utils