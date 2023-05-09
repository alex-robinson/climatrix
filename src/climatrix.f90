module climatrix

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use climatrix_defs
    use climate_interpolation
    use ncio 

    implicit none

    type climatrix_param_class

        integer :: ng           ! Number of geometries
        integer :: nc           ! Number of climates
        integer :: nx 
        integer :: ny 
        
        real(wp) :: dx          ! [m] Grid resolution
        real(wp) :: dist_max    ! [m] Maximum radius for neighborhood
        real(wp) :: dz          ! [m] Vertical distance for bins
    end type

    type climatrix_field
        character(len=56) :: name
        real(wp), allocatable :: z_srf(:,:,:,:)     ! [nx,ny,ng,nc]
        real(wp), allocatable :: mask(:,:,:,:)      ! [nx,ny,ng,nc]
        real(wp), allocatable :: var(:,:,:,:)       ! [nx,ny,ng,nc]

        real(wp), allocatable :: reg_mask(:,:)      ! [nx,ny]
        real(wp), allocatable :: reg_z_srf(:)       ! [nz]
        real(wp), allocatable :: reg_var(:,:)       ! [nr,nz]
    end type

    type climatrix_class

        type(climatrix_param_class) :: p 

        ! Axis variables [geometry axis, climate axis] that define interpolation matrix
        real(wp), allocatable :: x_geom(:)          ! [ng]
        real(wp), allocatable :: x_clim(:)          ! [nc]

        ! All fields are for a given month, or can represent annual fields too
        !type(climatrix_field) :: t2m
        !type(climatrix_field) :: pr
        type(climatrix_field) :: smb

    end type

    public 

contains

    subroutine climatrix_tables_gen(cax,name)

        implicit none

        type(climatrix_class), intent(INOUT) :: cax
        character(len=*), intent(IN) :: name 

        ! Local variables
        type(climatrix_field) :: fld

        ! ==================================================================
        ! Determine which field variable is being interpolated,
        ! extract field information from climatrix object
        
        select case(trim(name))

            case("smb")
                fld = cax%smb 

            case DEFAULT
                write(error_unit,*) "climatrix_interp:: field name not recognized."
                write(error_unit,*) "name = ", trim(name)
                stop

        end select

        ! ==================================================================
        ! Populate look up tables

        ! To do 




        ! ==================================================================
        ! Finally populate the right field object at output

        select case(trim(name))

            case("smb")
                cax%smb = fld

            case DEFAULT
                write(error_unit,*) "climatrix_interp:: field name not recognized."
                write(error_unit,*) "name = ", trim(name)
                stop

        end select

        return

    end subroutine climatrix_tables_gen

    subroutine climatrix_interp(var,z_srf,mask,x_geom,x_clim,name,cax,x_geom_subset,x_clim_subset)

        implicit none

        real(wp), intent(OUT) :: var(:,:)           ! [nx,ny] Calculated field
        real(wp), intent(IN)  :: z_srf(:,:)         ! [nx,ny] Current surface topography
        real(wp), intent(IN)  :: mask(:,:)          ! [nx,ny] Current mask (0: ice-free, 1: ice)
        real(wp), intent(IN)  :: x_geom             ! Current x_geom value
        real(wp), intent(IN)  :: x_clim             ! Current x_clim value
        character(len=*), intent(IN) :: name
        type(climatrix_class), intent(IN) :: cax 
        real(wp), intent(IN), optional :: x_geom_subset(:)
        real(wp), intent(IN), optional :: x_clim_subset(:)
        
        ! Local variables
        integer  :: i, j, l, m, ng, nc 
        integer, allocatable :: ii(:)
        integer, allocatable :: jj(:)
        integer :: i_iter(4), j_iter(4)
        integer  :: i0, i1, j0, j1 
        real(wp) :: wt_geom, wt_clim
        real(wp) :: var_lo, var_hi 
        type(climatrix_field) :: fld
        
        real(wp), allocatable :: var_geom_lo(:,:)
        real(wp), allocatable :: var_geom_hi(:,:)
        real(wp), allocatable :: var_geom_lo_now(:,:)
        real(wp), allocatable :: var_geom_hi_now(:,:)

        if (present(x_geom_subset)) then 

            ! Get indices of cax%x_geom that match the subset of interest
            call vec_in_vec(ii,cax%x_geom,x_geom_subset)

        else 
            allocate(ii(cax%p%ng))
            do i = 1, cax%p%ng 
                ii(i) = i 
            end do
        end if 

        if (present(x_clim_subset)) then 

            ! Get indices of cax%x_clim that match the subset of interest
            call vec_in_vec(jj,cax%x_clim,x_clim_subset)

        else 
            allocate(jj(cax%p%nc))
            do j = 1, cax%p%nc
                jj(j) = j
            end do
        end if 
        
        ng = size(ii)
        nc = size(jj) 

        ! Determine which field variable is being interpolated,
        ! extract field information from climatrix object
        select case(trim(name))

            case("smb")
                fld = cax%smb 

            case DEFAULT
                write(error_unit,*) "climatrix_interp:: field name not recognized."
                write(error_unit,*) "name = ", trim(name)
                stop

        end select

        ! Determine relevant axis indices bounding (x_geom,x_clim)

        if (x_geom .lt. minval(cax%x_geom(ii))) then 
            i0 = ii(1)
            i1 = ii(1)
        else if (x_geom .gt. maxval(cax%x_geom(ii))) then 
            i0 = ii(ng)
            i1 = ii(ng)
        else 
            do i = 1, ng 
                if (cax%x_geom(ii(i)) .ge. x_geom) then 
                    i0 = ii(i-1)
                    i1 = ii(i)
                    exit 
                end if
            end do 
        end if

        if (x_clim .lt. minval(cax%x_clim(jj))) then 
            j0 = jj(1)
            j1 = jj(1)
        else if (x_clim .gt. maxval(cax%x_clim(jj))) then 
            j0 = jj(nc)
            j1 = jj(nc)
        else 
            do j = 1, nc
                if (cax%x_clim(jj(j)) .ge. x_clim) then 
                    j0 = jj(j-1)
                    j1 = jj(j) 
                    exit 
                end if
            end do 
        end if

        ! Get interpolation weights
        wt_geom = (x_geom - cax%x_geom(i0)) / (cax%x_geom(i1) - cax%x_geom(i0))
        wt_clim = (x_clim - cax%x_clim(j0)) / (cax%x_clim(j1) - cax%x_clim(j0))
        
        write(output_unit,*)
        write(output_unit,*) "Current interpolation characteristics:"
        write(output_unit,*) "x_geom: ", wt_geom, x_geom, cax%x_geom(i0), cax%x_geom(i1)
        write(output_unit,*) "x_clim: ", wt_clim, x_clim, cax%x_clim(j0), cax%x_clim(j1)

        ! === Get new variable estimate =====

if (.FALSE.) then
    ! ajr: testing individual snapshot interpolation
        i = i1
        j = j1
        call climinterp_elevation_analog(var,z_srf,mask, &
                            fld%var(:,:,i,j),fld%z_srf(:,:,i,j),fld%mask(:,:,i,j), &
                            cax%p%dx,cax%p%dist_max,cax%p%dz)
else
        
        ! Perform full interpolation over climate axis and geometry axis 

        allocate(var_geom_lo(cax%p%nx,cax%p%ny))
        allocate(var_geom_hi(cax%p%nx,cax%p%ny))
        allocate(var_geom_lo_now(cax%p%nx,cax%p%ny))
        allocate(var_geom_hi_now(cax%p%nx,cax%p%ny))
        
        ! geom_lo, geom_hi 
        var_geom_lo = (1.0-wt_clim)*fld%var(:,:,i0,j0) + wt_clim*fld%var(:,:,i0,j1)
        var_geom_hi = (1.0-wt_clim)*fld%var(:,:,i1,j0) + wt_clim*fld%var(:,:,i1,j1)

        ! Now two arrays that are at the right climatic temperature, with different geometries

        ! geom_lo => geom_now
        i = i0 
        j = j0 

        write(output_unit,*) "Interpolating: ", cax%x_geom(i), cax%x_clim(j)

        call climinterp_elevation_analog(var_geom_lo_now,z_srf,mask, &
                        var_geom_lo,fld%z_srf(:,:,i,j),fld%mask(:,:,i,j), &
                        cax%p%dx,cax%p%dist_max,cax%p%dz)

        ! geom_hi => geom_hi_now
        i = i1
        j = j0 

        write(output_unit,*) "Interpolating: ", cax%x_geom(i), cax%x_clim(j)

        call climinterp_elevation_analog(var_geom_hi_now,z_srf,mask, &
                        var_geom_hi,fld%z_srf(:,:,i,j),fld%mask(:,:,i,j), &
                        cax%p%dx,cax%p%dist_max,cax%p%dz)

        ! Merge geom arrays
        do j = 1, cax%p%ny 
        do i = 1, cax%p%nx 

            if (var_geom_lo_now(i,j) .ne. MV .and. var_geom_hi_now(i,j) .ne. MV) then
                ! All values exist, get weighted average of geometry estimates

                var(i,j) = (1.0-wt_geom)*var_geom_lo_now(i,j) + (wt_geom)*var_geom_hi_now(i,j) 

            else if (var_geom_lo_now(i,j) .ne. MV) then 
                ! Only low value exists, use it

                var(i,j) = var_geom_lo_now(i,j)

            else if (var_geom_hi_now(i,j) .ne. MV) then 
                ! Only high value exists, use it

                var(i,j) = var_geom_hi_now(i,j)

            else 
                ! All neighbors are missing, keep missing for now.

                var(i,j) = MV 

            end if 

        end do 
        end do 

        var(42:52,50:60) = MV 

        ! Finally interpolate field with itself to ensure there are
        ! no missing values

        var_geom_lo_now = var 

        call climinterp_elevation_analog(var,z_srf,mask, &
                        var_geom_lo_now,z_srf,mask, &
                        cax%p%dx,cax%p%dist_max,cax%p%dz, &
                        mask_interp=var_geom_lo_now.eq.MV)

end if

        return

    end subroutine climatrix_interp

    subroutine climatrix_init(cax,filename,group)

        implicit none

        type(climatrix_class),  intent(INOUT) :: cax 
        character(len=*), intent(IN) :: filename            ! Parameter filename
        character(len=*), intent(IN) :: group               ! Parameter group name

        call climatrix_par_load(cax%p,filename,group)

        ! === Initialize matrix information ===

        if (allocated(cax%x_geom))  deallocate(cax%x_geom)
        if (allocated(cax%x_clim))  deallocate(cax%x_clim)

        allocate(cax%x_geom(cax%p%ng))
        allocate(cax%x_clim(cax%p%nc))

        cax%x_geom = MV
        cax%x_clim = MV

        ! === Initialize climatrix field variables === 

        cax%smb%name = "smb"
        call climatrix_field_alloc(cax%smb,cax%p%nx,cax%p%ny,cax%p%ng,cax%p%nc)

        return

    end subroutine climatrix_init

    subroutine climatrix_load_field_slice(cax,filename,group,x_geom)
        ! Load a slice of a field variable

        use nml 

        implicit none

        type(climatrix_class),  intent(INOUT) :: cax
        character(len=*), intent(IN) :: filename            ! Parameter filename
        character(len=*), intent(IN) :: group               ! Parameter group name
        real(wp),         intent(IN) :: x_geom              ! Current slice x_geom value
        
        ! Local variables 
        integer  :: i, j, nc
        real(wp) :: xg
        real(wp) :: xc(100)
        character(len=512) :: path 
        character(len=56)  :: fldrs(100)
        character(len=256) :: file 
        character(len=56)  :: zs_name 
        character(len=56)  :: mask_name 
        character(len=56)  :: var_name 
        character(len=512) :: file_path 
        real(wp)           :: mask_val 

        real(wp), allocatable :: mask_in(:,:) 

        allocate(mask_in(cax%p%nx,cax%p%ny))

        xc = MV

        call nml_read(filename,group,"xg",xg)
        call nml_read(filename,group,"xc",xc)
        call nml_read(filename,group,"path",path)
        call nml_read(filename,group,"fldrs",fldrs)
        call nml_read(filename,group,"file",file)
        call nml_read(filename,group,"zs_name",zs_name)
        call nml_read(filename,group,"mask_name",mask_name)
        call nml_read(filename,group,"var_name",var_name)
        call nml_read(filename,group,"mask_val",mask_val)
        
        ! Determine how many x_clim values there are
        nc = count(xc .ne. MV)

        ! === Print summary === 

        write(output_unit,*)
        write(output_unit,*) "slice parameters: x_geom = ", x_geom 
        write(output_unit,*)

        call nml_print("xg",xg)
        call nml_print("xc",xc(1:nc))
        call nml_print("path",path)
        call nml_print("fldrs",fldrs(1:nc))
        call nml_print("file",file)
        call nml_print("zs_name",zs_name)
        call nml_print("mask_name",mask_name)
        call nml_print("var_name",var_name)
        call nml_print("mask_val",mask_val)
        
        write(output_unit,*)
        
        if (xg .ne. x_geom) then 
            write(error_unit,*) "x_geom value does not match that of the slice parameters."
            write(error_unit,*) "x_geom = ", x_geom 
            write(error_unit,*) "xg     = ", xg 
            stop 
        end if

        ! === Load data for this slice [x_geom=xg,x_clim=xc] ===

        ! Determine x_geom where data should be added

        i = minloc(abs(cax%x_geom-x_geom),1)

        if (abs(cax%x_geom(i)-x_geom) .gt. TOL) then
            write(error_unit,*) "x_geom value does not match nearest cax%x_geom value."
            write(error_unit,*) "x_geom = ", x_geom 
            write(error_unit,*) "cax%x_geom(i), i = ", cax%x_geom(i), i 
            stop 
        end if

        do j = 1, nc 

            file_path = trim(path)//"/"//trim(fldrs(j))//"/"//trim(file)

            call nc_read(trim(file_path),zs_name,  cax%smb%z_srf(:,:,i,j),start=[1,1,1],count=[cax%p%nx,cax%p%ny,1])
            call nc_read(trim(file_path),var_name, cax%smb%var(:,:,i,j),  start=[1,1,1],count=[cax%p%nx,cax%p%ny,1])
            call nc_read(trim(file_path),mask_name,mask_in,start=[1,1,1],count=[cax%p%nx,cax%p%ny,1])
            
            ! Get ice mask
            cax%smb%mask(:,:,i,j) = 0.0
            where( abs(cax%smb%mask(:,:,i,j)-mask_in) .lt. TOL) cax%smb%mask(:,:,i,j) = 1.0

            write(output_unit,*) "Loaded slice, i, j: ", i, j 
            write(output_unit,*) "    z_srf: ", minval(cax%smb%z_srf(:,:,i,j)), maxval(cax%smb%z_srf(:,:,i,j))
            write(output_unit,*) "     mask: ", minval(cax%smb%mask(:,:,i,j)),  maxval(cax%smb%mask(:,:,i,j))
            write(output_unit,*) "      var: ", minval(cax%smb%var(:,:,i,j)),   maxval(cax%smb%var(:,:,i,j))
            
        end do 

        return

    end subroutine climatrix_load_field_slice

    subroutine climatrix_par_load(par,filename,group)

        use nml 

        implicit none

        type(climatrix_param_class), intent(INOUT) :: par 
        character(len=*), intent(IN) :: filename            ! Parameter filename
        character(len=*), intent(IN) :: group               ! Parameter group name

        call nml_read(filename,group,"ng",     par%ng)
        call nml_read(filename,group,"nc",     par%nc)
        call nml_read(filename,group,"nx",     par%nx)
        call nml_read(filename,group,"ny",     par%ny)
        
        call nml_read(filename,group,"dx",      par%dx)
        call nml_read(filename,group,"dist_max",par%dist_max)
        call nml_read(filename,group,"dz",      par%dz)
        
        write(output_unit,*) "climatrix parameters: ", trim(filename), " : ", trim(group)
        
        call nml_print("ng",     par%ng)
        call nml_print("nc",     par%nc)
        call nml_print("nx",     par%nx)
        call nml_print("ny",     par%ny)
        
        call nml_print("dx",        par%dx)
        call nml_print("dist_max",  par%dist_max)
        call nml_print("dz",        par%dz)
        
        return
        
    end subroutine climatrix_par_load

    subroutine climatrix_field_alloc(fld,nx,ny,ng,nc)

        implicit none

        type(climatrix_field), intent(INOUT) :: fld
        integer, intent(IN) :: nx
        integer, intent(IN) :: ny
        integer, intent(IN) :: ng
        integer, intent(IN) :: nc
        
        ! First ensure field is deallocated
        call climatrix_field_dealloc(fld)

        ! Allocate field variables
        allocate(fld%z_srf(nx,ny,ng,nc))
        allocate(fld%mask(nx,ny,ng,nc))
        allocate(fld%var(nx,ny,ng,nc))
        
        fld%z_srf = MV 
        fld%mask  = MV 
        fld%var   = MV 

        return

    end subroutine climatrix_field_alloc

    subroutine climatrix_field_dealloc(fld)

        implicit none

        type(climatrix_field), intent(INOUT) :: fld

        if (allocated(fld%z_srf)) deallocate(fld%z_srf)
        if (allocated(fld%mask))  deallocate(fld%mask)
        if (allocated(fld%var))   deallocate(fld%var)
        
        return

    end subroutine climatrix_field_dealloc



    subroutine vec_in_vec(ii,x,x_subset)
        
        implicit none

        integer, allocatable, intent(OUT) :: ii(:)
        real(wp), intent(IN) :: x(:)
        real(wp), intent(IN) :: x_subset(:)
        
        ! Local variables
        integer :: i, j, n 
        integer, allocatable :: ii_tmp(:)

        allocate(ii_tmp(size(x)))

        n = 0 
        do i = 1, size(x)

            do j = 1, size(x_subset)
                if (abs(x(i)-x_subset(j)) .lt. TOL) then 
                    ! Index found in subset add it to list 
                    n = n+1
                    ii_tmp(n) = i 
                    exit
                end if
            end do
        end do 

        if (n .eq. 0) then 
            write(error_unit,*) "vec_in_vec: no values found in subset."
            write(error_unit,*) "x: ", x 
            write(error_unit,*) "x_subset: ", x_subset 
            stop 
        else
            if (allocated(ii)) deallocate(ii)
            allocate(ii(n))
            ii = ii_tmp(1:n)
        end if 

        return

    end subroutine vec_in_vec

    subroutine which(ind,x,stat)
        ! Analagous to R::which function
        ! Returns indices that match condition x==.TRUE.

        implicit none 

        integer, allocatable, intent(OUT) :: ind(:)
        logical :: x(:)
        integer, optional :: stat 

        ! Local variables
        integer, allocatable :: tmp(:)
        integer :: n, i  

        n = count(x)
        allocate(tmp(n))
        tmp = 0 

        n = 0
        do i = 1, size(x) 
            if (x(i)) then 
                n = n+1
                tmp(n) = i 
            end if
        end do 

        if (present(stat)) stat = n 

        if (allocated(ind)) deallocate(ind)

        if (n .eq. 0) then 
            allocate(ind(1))
            ind = -1 
        else
            allocate(ind(n))
            ind = tmp(1:n)
        end if 
        
        return 

    end subroutine which

end module climatrix