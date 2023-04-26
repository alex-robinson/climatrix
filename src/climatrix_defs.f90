module climatrix_defs

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp

    real(wp), parameter :: MV = -9999.0 
    
    private
    public :: dp 
    public :: sp 
    public :: wp 
    public :: MV 
    
end module climatrix_defs