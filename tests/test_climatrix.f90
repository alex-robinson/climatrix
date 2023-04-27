program test_climatrix

    use ncio 
    use climatrix

    implicit none

    type(climatrix_class) :: cax
    character(len=256) :: par_path
    character(len=256) :: par_slice_path
    character(len=256) :: file_test 
    
    integer  :: i, j, inow, jnow
    real(wp) :: x_geom
    real(wp) :: x_clim 
    real(wp), allocatable :: smb(:,:) 

    
    ! ===========================================
    write(*,*) " "
    write(*,*) " Program: test_climatrix"
    write(*,*) 

    par_path       = "input/climatrix_greenland.nml"
    par_slice_path = "input/rembo-tipmip01.nml"
    
    file_test      = "output/test_climinterp_analog.nc" 

    ! Initialize cax object
    call climatrix_init(cax,par_path,"greenland")

    ! Define matrix axes
    cax%x_geom  = [0.0, 20.0, 50.0, 70.0, 80.0, 100.0]
    cax%x_clim  = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

    ! Load matrix data

    call climatrix_load_field_slice(cax,par_slice_path,"sims100",x_geom=100.0_wp)
    call climatrix_load_field_slice(cax,par_slice_path,"sims080",x_geom=80.0_wp)
    call climatrix_load_field_slice(cax,par_slice_path,"sims070",x_geom=70.0_wp)
    call climatrix_load_field_slice(cax,par_slice_path,"sims050",x_geom=50.0_wp)
    call climatrix_load_field_slice(cax,par_slice_path,"sims020",x_geom=20.0_wp)
    call climatrix_load_field_slice(cax,par_slice_path,"sims000",x_geom=0.0_wp)

    write(*,*)
    write(*,*) "All input data loaded."
    write(*,*) "    z_srf: ", minval(cax%smb%z_srf), maxval(cax%smb%z_srf)
    write(*,*) "     mask: ", minval(cax%smb%mask),  maxval(cax%smb%mask)
    write(*,*) "      smb: ", minval(cax%smb%var),   maxval(cax%smb%var)
    
    ! Write climatrix data to file
    call climatrix_write_init(cax,file_test,time_init=0.0_wp,units="years")

    ! Perform interpolation to a given location (x_geom,x_clim)

    allocate(smb(cax%p%nx,cax%p%ny))

    x_geom = 70.0 
    x_clim =  2.0 
    inow   = 4
    jnow   = 3 

    ! Load variable info

    call climatrix_interp(smb,cax%smb%z_srf(:,:,inow,jnow),cax%smb%mask(:,:,inow,jnow), &
                            x_geom,x_clim,"smb",cax, &
                            x_geom_subset=[0.0_wp,20.0_wp,50.0_wp,80.0_wp,100.0_wp])


    write(*,*)
    write(*,*) " test_climatrix complete."
    write(*,*)

contains

    subroutine climatrix_write_init(cax,filename,time_init,units)

        implicit none 

        type(climatrix_class), intent(IN) :: cax 
        character(len=*),      intent(IN) :: filename
        !real(wp),              intent(IN) :: xc(:)
        !real(wp),              intent(IN) :: yc(:)
        character(len=*),      intent(IN) :: units 
        real(wp),              intent(IN) :: time_init
        
        ! Local variables 
        character(len=16) :: xnm 
        character(len=16) :: ynm 
        character(len=16) :: grid_mapping_name

        xnm = "xc"
        ynm = "yc" 
        
        ! Create the empty netcdf file
        call nc_create(filename)

        ! Add grid axis variables to netcdf file
        ! call nc_write_dim(filename,xnm,x=xc*1e-3,units="km")
        ! call nc_write_dim(filename,ynm,x=yc*1e-3,units="km")
        call nc_write_dim(filename,xnm,x=1.0_wp,dx=1.0_wp,nx=cax%p%nx,units="")
        call nc_write_dim(filename,ynm,x=1.0_wp,dx=1.0_wp,nx=cax%p%ny,units="")
        call nc_write_dim(filename,"x_geom",x=cax%x_geom,units="%")
        call nc_write_dim(filename,"x_clim",x=cax%x_clim,units="K")
        
        call nc_write_dim(filename,"time",x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)

        ! Static information
        call nc_write(filename,"smb_var",  cax%smb%var,   dim1="xc",dim2="yc",dim3="x_geom",dim4="x_clim",long_name="Surface mass balance",units="m/yr")
        call nc_write(filename,"smb_mask", cax%smb%mask,  dim1="xc",dim2="yc",dim3="x_geom",dim4="x_clim",long_name="Surface mask",units="")
        call nc_write(filename,"smb_z_srf",cax%smb%z_srf, dim1="xc",dim2="yc",dim3="x_geom",dim4="x_clim",long_name="Surface elevation",units="m")
        
        return

    end subroutine climatrix_write_init

end program test_climatrix
