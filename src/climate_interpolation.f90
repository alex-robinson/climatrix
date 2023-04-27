module climate_interpolation

    ! Methods that can be used for determining the value of a given
    ! climatic variable like SMB from a reference field, given 
    ! current elevation and the elevation of the reference field, ie.,
    !
    ! var(i,j) = f(z_srf,var_ref,z_srf_ref)
    !

    use climatrix_defs 

    implicit none

contains

    subroutine climinterp_elevation_analog(var,z_srf,mask,var_ref,z_srf_ref,mask_ref,dx,dist_max,dz)

        implicit none

        real(wp), intent(OUT) :: var(:,:) 
        real(wp), intent(IN)  :: z_srf(:,:)
        real(wp), intent(IN)  :: mask(:,:) 
        real(wp), intent(IN)  :: var_ref(:,:) 
        real(wp), intent(IN)  :: z_srf_ref(:,:)
        real(wp), intent(IN)  :: mask_ref(:,:)
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dz
        real(wp), intent(IN)  :: dist_max
        

        ! Local variables
        integer  :: i, j, i1, j1, nx, ny 
        real(wp) :: z_now, z_min, z_max
        real(wp) :: dist, eps
        real(wp), allocatable :: wt(:,:)

        nx = size(var,1)
        ny = size(var,2)

        allocate(wt(nx,ny))

        eps = dx*1e-3

        var = MV 

        ! Loop over all points and determine var value from reference fields
        do j = 1, ny 
        do i = 1, nx 

            z_now = z_srf(i,j) 
            z_min = z_now - dz 
            z_max = z_now + dz 

            ! Get weighting of all points within z_now ± dz
            wt = 0.0 

            do j1 = 1, ny 
            do i1 = 1, nx 

                if (z_srf_ref(i1,j1) .ge. z_min .and. &
                    z_srf_ref(i1,j1) .le. z_max .and. &
                    mask_ref(i1,j1)  .eq. mask(i,j) ) then 
                    ! Point is within elevation band and of the same surface type (ice, ice-free),
                    ! calculate distance. 

                    dist = sqrt( (dx*(i1-i))**2 + (dx*(j1-j))**2 + eps**2)

                    if (dist .le. dist_max) then
                        ! Point is within distance of interest, calculate weight

                        ! Calculate weight as distance weighting in a radius, as a Modified Shepard's weighting
                        ! https://en.wikipedia.org/wiki/Inverse_distance_weighting
                        wt(i1,j1) = ((dist_max - dist)/(dist_max*dist))**2

                    end if 

                end if

            end do 
            end do 

            ! Check that some weights exist
            if (maxval(wt) .gt. 0.0) then 
                
                var(i,j) = sum(var_ref*wt,mask=wt.gt.0.0) / sum(wt,mask=wt.gt.0.0) 

            end if  

        end do
        end do 

        return

    end subroutine climinterp_elevation_analog

end module climate_interpolation

