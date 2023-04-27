program test_climatrix

    use climatrix

    implicit none

    type(climatrix_class) :: cax
    character(len=256) :: par_path
    character(len=256) :: par_slice_path

    real(wp) :: x_geom
    real(wp) :: x_clim 
    real(wp), allocatable :: z_srf(:,:) 
    real(wp), allocatable :: mask(:,:) 
    real(wp), allocatable :: smb(:,:) 

    ! ===========================================
    write(*,*) " "
    write(*,*) " Program: test_climatrix"
    write(*,*) 

    par_path       = "input/climatrix_greenland.nml"
    par_slice_path = "input/rembo-tipmip01.nml"

    ! Initialize cax object
    call climatrix_init(cax,par_path,"greenland")

    ! Define matrix axes
    cax%x_geom  = [0.0, 20.0, 50.0, 80.0, 100.0]
    cax%x_clim  = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

    ! Load matrix data

    call climatrix_load_field_slice(cax,par_slice_path,"sims100",x_geom=100.0_wp)
    call climatrix_load_field_slice(cax,par_slice_path,"sims080",x_geom=80.0_wp)
    call climatrix_load_field_slice(cax,par_slice_path,"sims050",x_geom=50.0_wp)
    call climatrix_load_field_slice(cax,par_slice_path,"sims020",x_geom=20.0_wp)
    call climatrix_load_field_slice(cax,par_slice_path,"sims000",x_geom=0.0_wp)

    write(*,*)
    write(*,*) "All input data loaded."
    write(*,*) "    z_srf: ", minval(cax%smb%z_srf), maxval(cax%smb%z_srf)
    write(*,*) "     mask: ", minval(cax%smb%mask),  maxval(cax%smb%mask)
    write(*,*) "      smb: ", minval(cax%smb%var),   maxval(cax%smb%var)
    
    ! Perform interpolation to a given location (x_geom,x_clim)

    allocate(z_srf(cax%p%nx,cax%p%ny))
    allocate(mask(cax%p%nx,cax%p%ny))
    allocate(smb(cax%p%nx,cax%p%ny))

    x_geom = 70.0 
    x_clim =  2.0 

    ! Load variable info

    !call climatrix_interp(var,z_srf,H_ice,x_geom,x_clim,name,cax)


    write(*,*)
    write(*,*) " test_climatrix complete."
    write(*,*)

contains



end program test_climatrix
