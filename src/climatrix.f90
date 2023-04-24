module climatrix

    use climatrix_defs

    implicit none

    type climatrix_par_class

        integer :: ng           ! Number of geometries
        integer :: nc           ! Number of climates
        integer :: nx 
        integer :: ny 
        
        real(wp) :: dx          ! [m] Grid resolution
        real(wp) :: rad_max     ! [m] Maximum radius for neighborhood

    end type

    type climatrix_field
        real(wp), allocatable :: mask(:,:,:,:)      ! [ng,nc,nx,ny]
        real(wp), allocatable :: z_srf(:,:,:,:)     ! [ng,nc,nx,ny]
        real(wp), allocatable :: var(:,:,:,:)       ! [ng,nc,nx,ny]
    end type

    type climatrix_class

        type(climatrix_par_class) :: p 

        ! Axis variables [geometry axis, climate axis] that define interpolation matrix
        real(wp), allocatable :: x_geom(:)          ! [ng]
        real(wp), allocatable :: x_clim(:)          ! [nc]

        ! All fields are for a given month, or can represent annual fields too
        !type(climatrix_field) :: t2m
        !type(climatrix_field) :: pr
        type(climatrix_field) :: smb

    end type


contains

    subroutine climatrix_init(cax,filename,group)

        implicit none

        type(climatrix_class),  intent(INOUT) :: cax 
        character(len=*), intent(IN) :: filename            ! Parameter filename
        character(len=*), intent(IN) :: group               ! Parameter group name





        return

    end subroutine climatrix_init

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
        call allocate(fld%mask(ng,nc,nx,ny))
        call allocate(fld%z_srf(ng,nc,nx,ny))
        call allocate(fld%var(ng,nc,nx,ny))
        
        fld%mask  = 0.0 
        fld%z_srf = 0.0 
        fld%var   = 0.0 

        return

    end subroutine climatrix_field_alloc

    subroutine climatrix_field_dealloc(fld)

        implicit none

        type(climatrix_field), intent(INOUT) :: fld

        if (allocated(fld%mask))  deallocate(fld%mask)
        if (allocated(fld%z_srf)) deallocate(fld%z_srf)
        if (allocated(fld%var))   deallocate(fld%var)
        
        return

    end subroutine climatrix_field_dealloc

end module climatrix