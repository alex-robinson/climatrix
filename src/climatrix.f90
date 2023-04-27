module climatrix

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use climatrix_defs
    use ncio 

    implicit none

    type climatrix_param_class

        integer :: ng           ! Number of geometries
        integer :: nc           ! Number of climates
        integer :: nx 
        integer :: ny 
        
        real(wp) :: dx          ! [m] Grid resolution
        real(wp) :: rad_max     ! [m] Maximum radius for neighborhood

    end type

    type climatrix_field
        character(len=56) :: name
        real(wp), allocatable :: z_srf(:,:,:,:)     ! [ng,nc,nx,ny]
        real(wp), allocatable :: mask(:,:,:,:)      ! [ng,nc,nx,ny]
        real(wp), allocatable :: var(:,:,:,:)       ! [ng,nc,nx,ny]
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

    subroutine climatrix_interp(var,z_srf,x_geom,x_clim,name,cax)

        implicit none

        real(wp), intent(OUT) :: var(:,:)           ! [nx,ny] Calculated field
        real(wp), intent(IN)  :: z_srf(:,:)         ! [nx,ny] Current surface topography
        real(wp), intent(IN)  :: x_geom             ! Current x_geom value
        real(wp), intent(IN)  :: x_clim             ! Current x_clim value
        character(len=*), intent(IN) :: name
        type(climatrix_class), intent(IN) :: cax 

        ! Local variables
        integer  :: i, j
        integer  :: i0, i1, j0, j1 
        real(wp) :: wt_geom, wt_clim
        type(climatrix_field) :: fld
        
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

        if (x_geom .lt. minval(cax%x_geom)) then 
            i0 = 1
            i1 = 1 
        else if (x_geom .gt. maxval(cax%x_geom)) then 
            i0 = cax%p%ng 
            i1 = cax%p%ng 
        else 
            do i = 1, cax%p%ng 
                if (cax%x_geom(i) .ge. x_geom) then 
                    i0 = i-1
                    i1 = i 
                    exit 
                end if
            end do 
        end if

        if (x_clim .lt. minval(cax%x_clim)) then 
            j0 = 1
            j1 = 1 
        else if (x_clim .gt. maxval(cax%x_clim)) then 
            j0 = cax%p%ng 
            j1 = cax%p%ng 
        else 
            do j = 1, cax%p%ng 
                if (cax%x_clim(j) .ge. x_clim) then 
                    j0 = j-1
                    j1 = j 
                    exit 
                end if
            end do 
        end if

        write(output_unit,*) "x_geom: ", x_geom, cax%x_geom(i0), cax%x_geom(i1)
        write(output_unit,*) "x_clim: ", x_clim, cax%x_clim(i0), cax%x_clim(i1)


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
        call climatrix_field_alloc(cax%smb,cax%p%ng,cax%p%nc,cax%p%nx,cax%p%ny)

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

            call nc_read(trim(file_path),zs_name,  cax%smb%z_srf(i,j,:,:),start=[1,1,1],count=[cax%p%nx,cax%p%ny,1])
            call nc_read(trim(file_path),var_name, cax%smb%var(i,j,:,:),  start=[1,1,1],count=[cax%p%nx,cax%p%ny,1])
            call nc_read(trim(file_path),mask_name,mask_in,start=[1,1,1],count=[cax%p%nx,cax%p%ny,1])
            
            ! Get ice mask
            cax%smb%mask(i,j,:,:) = 0.0
            where( abs(cax%smb%mask(i,j,:,:)-mask_in) .lt. TOL) cax%smb%mask(i,j,:,:) = 1.0

            write(output_unit,*) "Loaded slice, i, j: ", i, j 
            write(output_unit,*) "    z_srf: ", minval(cax%smb%z_srf(i,j,:,:)), maxval(cax%smb%z_srf(i,j,:,:))
            write(output_unit,*) "    mask: ", minval(cax%smb%mask(i,j,:,:)), maxval(cax%smb%mask(i,j,:,:))
            write(output_unit,*) "      var: ", minval(cax%smb%var(i,j,:,:)),   maxval(cax%smb%var(i,j,:,:))
            
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
        
        call nml_read(filename,group,"dx",     par%dx)
        call nml_read(filename,group,"rad_max",par%rad_max)
        
        write(output_unit,*) "climatrix parameters: ", trim(filename), " : ", trim(group)
        
        call nml_print("ng",     par%ng)
        call nml_print("nc",     par%nc)
        call nml_print("nx",     par%nx)
        call nml_print("ny",     par%ny)
        
        call nml_print("dx",     par%dx)
        call nml_print("rad_max",par%rad_max)
        
        return
        
    end subroutine climatrix_par_load

    subroutine climatrix_field_alloc(fld,ng,nc,nx,ny)

        implicit none

        type(climatrix_field), intent(INOUT) :: fld
        integer, intent(IN) :: ng
        integer, intent(IN) :: nc
        integer, intent(IN) :: nx
        integer, intent(IN) :: ny
        
        ! First ensure field is deallocated
        call climatrix_field_dealloc(fld)

        ! Allocate field variables
        allocate(fld%z_srf(ng,nc,nx,ny))
        allocate(fld%mask(ng,nc,nx,ny))
        allocate(fld%var(ng,nc,nx,ny))
        
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

end module climatrix