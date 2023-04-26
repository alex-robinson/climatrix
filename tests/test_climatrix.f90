program test_climatrix

    use climatrix

    implicit none

    type(climatrix_class) :: cax
    character(len=256) :: par_path
    character(len=256) :: par_slice_path

    ! ===========================================
    write(*,*) " "
    write(*,*) " Program: test_climatrix"
    write(*,*) 

    par_path       = "input/climatrix_greenland.nml"
    par_slice_path = "input/rembo-tipmip01.nml"

    call climatrix_init(cax,par_path,"greenland")

    call climatrix_load_field_slice(cax,par_slice_path,"sims100",x_geom=100.0_wp)


    write(*,*)
    write(*,*) " test_climatrix complete."
    write(*,*)

contains



end program test_climatrix
