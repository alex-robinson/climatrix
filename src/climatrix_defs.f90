module climatrix_defs

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp

    private
    public :: dp 
    public :: sp 
    public :: wp 
    
end module climatrix_defs