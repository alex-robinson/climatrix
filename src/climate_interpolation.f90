module climate_interpolation

    ! Methods that can be used for determining the value of a given
    ! climatic variable like SMB from a reference field, given 
    ! current elevation and the elevation of the reference field, ie.,
    !
    ! var(i,j) = f(z_srf,var_ref,z_srf_ref)
    !

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use climatrix_defs 

    implicit none

contains

    subroutine climinterp_elevation_analog(var,z_srf,mask,var_ref,z_srf_ref,mask_ref,dx,dist_max,dz,mask_interp)

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
        logical,  intent(IN), optional :: mask_interp(:,:) 

        ! Local variables
        integer  :: i, j, i1, j1, nx, ny 
        integer  :: nij_dx  
        real(wp) :: z_now, z_min, z_max
        real(wp) :: dist, eps
        real(wp), allocatable :: wt(:,:)
        logical,  allocatable :: mask_interpolate(:,:)
        
        nx = size(var,1)
        ny = size(var,2)

        allocate(wt(nx,ny))
        allocate(mask_interpolate(nx,ny))

        eps = dx*1e-3
        nij_dx = ceiling(dist_max / dx)

        if (present(mask_interp)) then
            mask_interpolate = mask_interp
        else
            mask_interpolate = .TRUE. 
        end if 

        where (mask_interpolate) var = MV

        ! Loop over all points and determine var value from reference fields
        do j = 1, ny 
        do i = 1, nx 

            if (mask_interpolate(i,j)) then 
                ! Interpolation desired at this point, proceed

                z_now = z_srf(i,j) 
                z_min = z_now - dz 
                z_max = z_now + dz 

                ! Get weighting of all points within z_now ± dz
                wt = 0.0 

                do j1 = j-nij_dx, j+nij_dx
                do i1 = i-nij_dx, i+nij_dx

                    if (i1 .ge. 1 .and. i1 .le. nx .and. j1 .ge. 1 .and. j1 .le. ny) then
                        ! Current neighbor is within domain boundaries

                        if (z_srf_ref(i1,j1) .ge. z_min .and. z_srf_ref(i1,j1) .le. z_max .and. &
                            mask_ref(i1,j1)  .eq. mask(i,j) .and. var_ref(i1,j1) .ne. MV ) then 
                            ! Point is within elevation band, of the same surface type (ice, ice-free),
                            ! and not a missing value, so calculate distance. 

                            dist = sqrt( (dx*(i1-i))**2 + (dx*(j1-j))**2 + &
                                            (z_now-z_srf_ref(i1,j1))**2 + eps**2 )

                            if (dist .le. dist_max) then
                                ! Point is within distance of interest, calculate weight

                                ! Calculate weight as distance weighting in a radius, as a Modified Shepard's weighting
                                ! https://en.wikipedia.org/wiki/Inverse_distance_weighting
                                wt(i1,j1) = ((dist_max - dist)/(dist_max*dist))**2

                            end if 

                        end if

                    end if 

                end do 
                end do 

                ! Check that some weights exist
                if (maxval(wt) .gt. 0.0) then 
                    
                    var(i,j) = sum(var_ref*wt,mask=wt.gt.0.0) / sum(wt,mask=wt.gt.0.0) 

                end if  

            end if 

        end do
        end do 

        return

    end subroutine climinterp_elevation_analog

    subroutine climinterp_gen_lookup(var,z_srf,var_ref,z_srf_ref,mask_interp_ref)

        implicit none

        real(wp), intent(OUT) :: var(:)         ! [nz] Variable estimated at target elevation levels
        real(wp), intent(IN)  :: z_srf(:)       ! [nz] Target elevation levels
        real(wp), intent(IN)  :: var_ref(:,:) 
        real(wp), intent(IN)  :: z_srf_ref(:,:)
        logical,  intent(IN)  :: mask_interp_ref(:,:)   ! Allowed source points to use in calculation of var

        ! Local variables
        integer  :: i, j, k, nx, ny, nz
        real(wp) :: zmin, zmax
        real(wp) :: wt_tot
        real(wp), allocatable :: wt(:,:) 

        ! Get total bin centers we will estimate
        nx = size(var_ref,1)
        ny = size(var_ref,2)
        nz = size(z_srf,1)

        ! Loop over bins and calculate var value

        do k = 1, nz 

            ! Get boundaries of bin
            if (k .eq. 1) then 
                zmin = z_srf(k) - 0.5*(z_srf(k+1)-z_srf(k))
            else
                zmin = 0.5*(z_srf(k-1)+z_srf(k))
            end if

            if (k .eq. nz) then 
                zmax = z_srf(k) + 0.5*(z_srf(k)-z_srf(k-1))
            else
                zmax = 0.5*(z_srf(k)+z_srf(k+1))
            end if

            ! Get the mean value of all points in region of interest within vertical bin
            
            wt = 0.0

            do j = 1, ny 
            do i = 1, nx 

                if (mask_interp_ref(i,j)) then 
                    if (z_srf_ref(i,j) .gt. zmin .and. z_srf_ref(i,j) .le. zmax) then

                        wt = 1.0

                    end if
                end if

            end do 
            end do

            wt_tot = sum(wt)

            if (wt_tot .gt. 0.0) then 
                var(k) = sum(var_ref*wt) / wt_tot
            else
                var(k) = MV
            end if

        end do

        return

    end subroutine climinterp_gen_lookup

    subroutine climinterp_interp_lookup(var,z_srf,var_ref,z_srf_ref)

        implicit none

        real(wp), intent(OUT) :: var
        real(wp), intent(IN)  :: z_srf
        real(wp), intent(IN)  :: var_ref(:)
        real(wp), intent(IN)  :: z_srf_ref(:)

        ! Local variables
        integer  :: k, nz 
        integer  :: k0, k1 
        real(wp) :: wt 

        nz = size(var_ref,1)

        if (count(var_ref .eq. MV) .eq. nz) then
             ! No interpolation values available, retain missing value

             var = MV 
        
        else
            ! Some interpolation values are available, proceed

            ! Check extreme boundaries
            do k = 1 , nz 
                if (var_ref(k) .ne. MV) then
                    k0 = k 
                    exit
                end if
            end do 

            do k = nz, 1, -1
                if (var_ref(k) .ne. MV) then
                    k1 = k 
                    exit
                end if
            end do 
            
            if (z_srf .lt. z_srf_ref(k0)) then
                ! Current elevation is too low, set minimum elevation variable value

                var = var_ref(k0)
            
            else if (z_srf .gt. z_srf_ref(k1)) then 
                ! Current elevation is too high, set maximum elevation variable value

                var = var_ref(k1)

            else if (k0 .eq. k1) then
                ! Only one value avaiable, us it

                var = var_ref(k1)
            
            else
                ! Current elevation is within the lookup table, proceed to linear interpolation

                ! First find upper index
                do k = 1, nz
                    if (var_ref(k) .ne. MV .and. z_srf_ref(k) .gt. z_srf) then
                        k1 = k 
                        exit
                    end if
                end do

                ! Next find lower index
                do k = k1, 1, -1
                    if (var_ref(k) .ne. MV .and. z_srf_ref(k) .le. z_srf) then
                        k0 = k 
                        exit
                    end if
                end do

                if (k0 .ge. k1) then
                    write(error_unit,*) "climinterp_interp_lookup:: Error with bounding indices."
                    write(error_unit,*) "k0, k1: ", k0, k1 
                    stop 
                end if

                ! Perform linear interpolation to current elevation

                wt = (z_srf-z_srf_ref(k0)) / (z_srf_ref(k1)-z_srf_ref(k0))

                var = (1.0-wt)*z_srf_ref(k0) + wt*z_srf_ref(k1)

            end if 

        end if

        return

    end subroutine climinterp_interp_lookup

end module climate_interpolation

