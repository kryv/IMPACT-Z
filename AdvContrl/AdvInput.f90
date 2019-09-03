!!QZ20120210 read "ooutputflag" the last entry for output emittance %
!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! Inputclass: Input class in I/O module of CONTROL layer. 
! Version: 2.0
! Author: Ji Qiang, LBNL, 7/05/04
! Description: This class defines functions to input the global
!              beam and computational parameters and the lattice input
!              parameters in the accelerator.
! Comments: J. Q. modified the source code so that the user can put
!           coments line starting with "!" for each number line in the
!           input file "test.in".
!----------------------------------------------------------------
module AdvInputclass
    integer :: ioinpfile
    
    interface in_Input
        module procedure in1_Input, in2_Input
    end interface
    
    contains
! Input all parameters except beam line element parameters.
    subroutine in1_Input(fname,odim,onp,onx,ony,onz,oflagbc,oflagdist, &
                         orstartflg,oflagmap,distparam,nparam,obcurr,&
                         obkenergy,obmass,obcharge,obfreq,oxrad,oyrad,&
                         operdlen,onblem,onpcol,onprow,oflagerr,oflagdiag,&
                         oflagsbstp,ophsini,onchrg,onptlist,ocurrlist,&
                         oqmcclist,ooutputflag)

        implicit none
        include 'mpif.h'
        character, intent(in) :: fname*256
        integer, intent(out) :: odim,onp,onx,ony,onz,oflagbc,oflagdist
        integer, intent(out) :: orstartflg,oflagmap,onblem,onpcol,onprow 
        integer, intent(out) :: oflagerr,oflagdiag,oflagsbstp,onchrg
        integer, intent(in) :: nparam
        double precision, dimension(nparam), intent(out) :: distparam
        double precision, intent(out) :: obcurr,obkenergy,obmass
        double precision, intent(out) :: obcharge,obfreq,operdlen,&
                                         oxrad,oyrad,ophsini
        double precision, dimension(:), intent(inout) :: ocurrlist,oqmcclist
        integer, dimension(:), intent(inout) :: onptlist
        double precision, intent(out) :: ooutputflag
        integer :: myrank,ierr,i,j,inloc,ios,itot
        
        double precision :: xjunk
        integer :: njunk1,njunk2,njunk3
        character*1 comst

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        inloc = 0
        ios = 0
        ioinpfile = 0
        
        if(myrank.eq.0) then
            open(unit=13,file=trim(fname),status='old',iostat=ioinpfile)
        endif

        call MPI_BCAST(ioinpfile,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (ioinpfile.ne.0) then
            if(myrank.eq.0) print*, 'AdvInput.f90 open(): Input file (',trim(fname),') does not found.'
            return
        endif
        
        if(myrank.eq.0) then
            j = 0
            do while (inloc.ne.11)
                ios = 5010 !Error Code 5010: LIBERROR_READ_VALUE
                do while (ios.eq.5010)
                    read(13,*,iostat=ios) xjunk     
                    j = j + 1
                enddo
                backspace(13)
                
                select case(inloc)
                    case (0)
                        read(13,*) onpcol,onprow
                        inloc = 1
                    case (1)
                        read(13,*) odim, onp, oflagmap, oflagerr, oflagdiag
                        inloc = 2
                    case (2)
                        read(13,*) onx, ony, onz, oflagbc, oxrad,oyrad,operdlen
                        inloc = 3
                    case (3)
                        read(13,*) oflagdist, orstartflg, oflagsbstp, onchrg
                        inloc = 4
                    case (4)                    
                        read(13,*) onptlist(1:onchrg)
                        inloc = 5
                    case (5)
                        read(13,*) ocurrlist(1:onchrg)
                        inloc = 6
                    case (6)
                        read(13,*) oqmcclist(1:onchrg)
                        inloc = 7
                    case (7)
                        read(13,*) distparam(1:7)
                        inloc = 8
                    case (8)
                        read(13,*) distparam(8:14)
                        inloc = 9                        
                    case (9)
                        read(13,*) distparam(15:21)
                        inloc = 10                          
                    case (10)
                        read(13,*) obcurr,obkenergy,obmass,obcharge,obfreq,ophsini,ooutputflag
                        inloc = 11
                end select
            enddo
                
            itot = 0
            ios = 0
            do while(ios.ge.0)
                ios = 5010
                do while (ios.eq.5010)
                    read(13,*,iostat=ios) xjunk
                enddo
                backspace(13)
                read(13,*,iostat=ios) xjunk,njunk1,njunk2,njunk3
                itot = itot + 1
                if (njunk3.eq.-99) then
                    ios = -99
                    itot = itot + 1
                endif
            enddo
            
            onblem=itot-1
            !write(6,*)'onblem = ',onblem
            rewind(13)
            
            do i = 1, j
                read(13,*)comst
            enddo
        endif  
        
        call MPI_BCAST(onpcol,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onprow,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(odim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagmap,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagerr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagdiag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onchrg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onptlist(1),onchrg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ocurrlist(1),onchrg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oqmcclist(1),onchrg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ony,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagbc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oxrad,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oyrad,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(operdlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagdist,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(orstartflg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagsbstp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onblem,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(distparam(1),nparam,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(obcurr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(obkenergy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(obmass,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(obcharge,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(obfreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ophsini,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
    end subroutine in1_Input

    
    ! Input beam line element parameters.
    !        subroutine in2_Input(onblem,operd,oblength,obnseg,obmpstp,&
    !                             obtype,value1,value2,value3,value4,value5)
    subroutine in2_Input(onblem,oblength,obnseg,obmpstp,obtype,varr)
        implicit none
        include 'mpif.h'
        !integer,intent(in) :: onblem,operd
        integer,intent(in) :: onblem
        integer,intent(out) :: obnseg(onblem)
        integer,intent(out) :: obmpstp(onblem)
        integer,intent(out) :: obtype(onblem)
        double precision,intent(out) :: oblength(onblem)
        double precision,dimension(28,onblem),intent(out) :: varr
        integer :: i,irf
        integer :: myrank,ierr
        
        double precision :: xjunk
        integer :: ios

        ios = 0
        varr = 0.0
        
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
        
        if(myrank.eq.0) then
            i = 1
            do while(ios.ge.0)
                ios = 5010
                do while (ios.eq.5010)
                    read(13,*,iostat=ios) xjunk
                enddo
                backspace(13)
                read(13,*,iostat=ios) oblength(i),obnseg(i),obmpstp(i),obtype(i), varr(2:,i)
                if (obtype(i).eq.-99) ios = -99
                i = i+1
            enddo
            close(13)
        endif
        
        call MPI_BCAST(oblength,onblem,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(varr,onblem*28,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(obnseg,onblem,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(obmpstp,onblem,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(obtype,onblem,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        
    end subroutine in2_Input

end module AdvInputclass

