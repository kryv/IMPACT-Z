!----------------------------------------------------------------
! (c) Copyright, 2004 by the Regents of the University of California.
! PhysConstclass: physical constants class in CONSTANTS module of 
!                 DATA STRUCTURE layer. (PC version)
! Version: 2.0
! Author: Ji Qiang, LBNL, 1/9/04
! Description: This class defines the physical constant parameters used
!              in the simulation.
! Comments:
!----------------------------------------------------------------
      module PhysConstclass
!        use mpistub
        implicit none
      
        !physical parameters and constants ---------------------
        double precision :: Pi = 2.0*asin(1.0)
        double precision :: Clight = 299792458.0 !speed of light in vaccum
        double precision :: Epsilon0 = 8.854187817e-12 !permittivity of vacuum
        double precision :: Rad2deg !conversion factor from radian to degree
        double precision :: Scxl !length scale
        double precision :: Scfreq !time scale
      contains
        subroutine construct_PhysConst(freq)
        double precision, intent(in) :: freq

        Scxl = Clight/(2.0*Pi*freq)
        Scfreq = freq
        Rad2deg = 180.0/Pi

        end subroutine construct_PhysConst
 
      end module PhysConstclass
