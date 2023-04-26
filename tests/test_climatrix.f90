program test_climatrix

    use climatrix

    implicit none

    type(climatrix_class) :: cax

    ! ===========================================
    write(*,*) " "
    write(*,*) " Program: test_climatrix"
    write(*,*) 

    call climatrix_init(cax,"input/climatrix_greenland.nml","greenland")

    



    write(*,*)
    write(*,*) " test_climatrix complete."
    write(*,*)

contains



end program test_climatrix
