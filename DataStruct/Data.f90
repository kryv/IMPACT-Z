!----------------------------------------------------------------
! (c) Copyright, 2003 by the Regents of the University of California.
! Dataclass: Field data class in DATA STRUCTURE layer (PC version). 
! Version: 2.0
! Author: Ji Qiang, LBNL, 6/3/03
! Description: This class stores the rf cavity data Ez, Ez', Ez'' on the 
!              axis; Fourier coefficients of Ez on the axis; Ez(r,z),
!              Er(r,z), Htheta(r,z) on the r-z grid plane; and Ex(x,y,z),
!              Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z), Bz(x,y,z) on
!              uniform x, y, z grid.
! Comments: 
!----------------------------------------------------------------
      module Dataclass
!        use mpistub
        save
!-----------------------------------------------------------------------
! using the x-y-z field data (Ex,Ey,Ez,Bx,By,Bz) directly.
        !number of grid points along x, y, and z direction.
        integer :: NxIntvRfg = 1
        integer :: NyIntvRfg = 1
        integer :: NzIntvRfg = 1
        integer :: Nmultigrid = 1
        !range in x, y, and zdirections.
        double precision :: XmaxRfg,XminRfg,YmaxRfg,YminRfg,ZmaxRfg,ZminRfg
        ! discrete Ex(x,y,z), Ey(x,y,z), Ez(x,y,z) and Bx(x,y,z), By(x,y,z), and 
        ! Bz(x,y,z) rf data. Here, the grid is uniform in x, y and z.
        double precision,allocatable,dimension(:,:,:) :: &
                        Exgrid,Eygrid,Ezgrid,Bxgrid,Bygrid,Bzgrid,&
                        Exgrid2,Eygrid2,Ezgrid2,Bxgrid2,Bygrid2,Bzgrid2,&
                        Exgrid3,Eygrid3,Ezgrid3,Bxgrid3,Bygrid3,Bzgrid3
!-----------------------------------------------------------------------
! using the r-z field data (Er,Ez,Htheta) directly.
        !number of grid points along r direction.
        integer :: NrIntvRf = 1
        !number of grid points along z direction.
        integer :: NzIntvRf = 1
        !range in r and z directions.
        double precision :: RmaxRf,RminRf,ZmaxRf,ZminRf
        ! discrete Ez(r,z), Er(r,z) and Htheta(r,z) rf data. Here, the grid
        ! is uniform in r and z.
        double precision,allocatable,dimension(:,:) :: &
               ezdata,erdata,btdata
!-----------------------------------------------------------------------
! using only on-axis field data and its derivities.
        !initial number of grid points on the axis.
        integer, parameter :: Ndataini = 5000
        ! discrete Ez(0,z), Ez'(0,z), Ez''(0,z) rf data.
        double precision,dimension(Ndataini) :: zdat,edat,epdat,eppdat
!----------------------------------------------------------------------
! using the Fourier coefficients
        !number of Fourier expansion coefficients.
        integer, parameter :: NcoefF = 101
        double precision,dimension(NcoefF) :: Fcoef
        !Fcoef(1): constant
        !Fcoef(2k): cosine term
        !Fcoef(2k+1): sine term
!----------------------------------------------------------------------
! using the RFQ data
!   i : cell id
!   Rfqpara (1,i) : scale
!           (2,i) : synchronous phase
!           (3,i) : radius
!           (4,i) : modulation
!           (5,i) : r0
!           (6,i) : A10
!           (7,i) : A0
!           (8,i) : A12I4
!           (9,i): A1
!           (10,i): A30I0
!           (11,i): A21I2
!           (12,i): A32I4
!           (13,i): A23I6
!           (14,i): standard cell count
        integer :: Rfqid, Rfqntot, Rfqtseg, tmprfqncl
        integer,allocatable,dimension(:) :: Rfqnseg,Rfqpstp,Rfqltyp,Rfqsegs,Rfqppos
        integer,allocatable,dimension(:,:) :: Rfqprof
        double precision,allocatable,dimension(:) :: Rfqzlen, Rfqzpos
        double precision,allocatable,dimension(:,:) :: Rfqpara
!----------------------------------------------------------------------
        ! practical number of grid data on the axis or Fourier coefficients.
        integer :: Ndata, iodatafile
        character(256) :: dataclass_dir = ''
      contains
        !Initialize the data storage arrays.
        subroutine init_Data()
        implicit none
        include 'mpif.h' 
        integer :: i,j

        NzIntvRf = 1
        NrIntvRf = 1
        allocate(ezdata(NzIntvRf+1,NrIntvRf+1))
        allocate(erdata(NzIntvRf+1,NrIntvRf+1))
        allocate(btdata(NzIntvRf+1,NrIntvRf+1))

        do j = 1, NrIntvRf+1
          do i = 1, NzIntvRf+1
            ezdata(i,j) = 0.0
            erdata(i,j) = 0.0
            btdata(i,j) = 0.0
          enddo
        enddo

        do i = 1, Ndataini
          zdat(i) = 0.0
          edat(i) = 0.0
          epdat(i) = 0.0
          eppdat(i) = 0.0
        enddo

        Ndata = 1
        RminRf = 0.0
        RmaxRf = 1.0
        ZminRf = 0.0
        ZminRf = 1.0

!initialization of Ex,Ey,Ez,Bx,By,Bz
        NxIntvRfg = 1
        NyIntvRfg = 1
        NzIntvRfg = 1
        
        allocate(Exgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Eygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Ezgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bxgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bzgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))

        Exgrid = 0.0
        Eygrid = 0.0
        Ezgrid = 0.0
        Bxgrid = 0.0
        Bygrid = 0.0
        Bzgrid = 0.0
        
        XminRfg = 0.0
        XmaxRfg = 1.0
        YminRfg = 0.0
        YmaxRfg = 1.0
        ZminRfg = 0.0
        ZmaxRfg = 1.0

        end subroutine init_Data

        subroutine destruct_Data()
        implicit none
        include 'mpif.h' 

        if (allocated(ezdata)) deallocate(ezdata)
        if (allocated(erdata)) deallocate(erdata)
        if (allocated(btdata)) deallocate(btdata)
        if (allocated(Exgrid)) deallocate(Exgrid)
        if (allocated(Eygrid)) deallocate(Eygrid)
        if (allocated(Ezgrid)) deallocate(Ezgrid)
        if (allocated(Bxgrid)) deallocate(Bxgrid)
        if (allocated(Bygrid)) deallocate(Bygrid)
        if (allocated(Bzgrid)) deallocate(Bzgrid)

        if (allocated(Exgrid2)) deallocate(Exgrid3)
        if (allocated(Eygrid2)) deallocate(Eygrid3)
        if (allocated(Ezgrid2)) deallocate(Ezgrid3)
        if (allocated(Bxgrid2)) deallocate(Bxgrid3)
        if (allocated(Bygrid2)) deallocate(Bygrid3)
        if (allocated(Bzgrid2)) deallocate(Bzgrid3)

        if (allocated(Exgrid3)) deallocate(Exgrid3)
        if (allocated(Eygrid3)) deallocate(Eygrid3)
        if (allocated(Ezgrid3)) deallocate(Ezgrid3)
        if (allocated(Bxgrid3)) deallocate(Bxgrid3)
        if (allocated(Bygrid3)) deallocate(Bygrid3)
        if (allocated(Bzgrid3)) deallocate(Bzgrid3)
        
        end subroutine destruct_Data

        !read in the discrete field data Ez(0,z), Ez'(0,z), Ez''(0,z) 
        !distribution along axis zdat from files "rfdatax or rfdataxx 
        !or rfdataxxx".
        subroutine read1_Data(ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1
        character(256) :: fileid

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          write(fileid,*), ifile
          open(14, file=trim(dataclass_dir)//'rfdata'//trim(adjustl(fileid)), status='old', iostat=iodatafile)
        endif
        
        call MPI_BCAST(iodatafile,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (iodatafile.ne.0) then
          if(myrank.eq.0) print*, 'Data.f90 open(): RF data file (',&
                           trim(dataclass_dir)//'rfdata'//trim(adjustl(fileid)),&
                           ') does not found.'
          close(14)
          return
        endif
        
        if(myrank.eq.0) then
          n = 0
50        continue
            read(14,*,end=77)tmp1,tmp2,tmp3,tmp4
            n = n + 1
            zdat(n) = tmp1
            edat(n) = tmp2
            epdat(n) = tmp3
            eppdat(n) = tmp4
            !write(15,100)zdat(n),edat(n),epdat(n),eppdat(n)
            !divided by 100 is due to the unit in rf data is cm.
          goto 50
77        continue
          close(14)
!          close(15)
          Ndata = n
          zdat1 = zdat(1)
          do i = 1, Ndata
            zdat(i) = zdat(i) - zdat1
          enddo
        endif
100     format(6(1x,1pe15.8))

        call MPI_BCAST(Ndata,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(zdat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(edat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(epdat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(eppdat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndata
        end subroutine read1_Data

        ! read in discrete Ez(r,z), Er(r,z) and Btheta(r,z) rf data from
        ! files "1Tx.T7 or 1Txx.T7 or 1Txxx.T7. Here, the grid
        ! is uniform in r and z.
        subroutine read2_Data(ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,nn,Ndatalc,j,tmpint
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0
        character(256) :: fileid

        mu0 = 4*2*asin(1.0)*1.0e-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          write(fileid,*), ifile
          open(14, file=trim(dataclass_dir)//'1T'//trim(adjustl(fileid))//'.T7', status='old',iostat=iodatafile)
        endif
        
        call MPI_BCAST(iodatafile,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (iodatafile.ne.0) then
          if(myrank.eq.0) print*, 'Data.f90 open(): Field file (',&
                           trim(dataclass_dir)//'1T'//trim(adjustl(fileid))//'.T7',&
                           ') does not found.'
          close(14)
          return
        endif
        
        if(myrank.eq.0) then
          ! the input range units are cm
          read(14,*,end=33)tmp1,tmp2,tmpint
          ZminRf = tmp1/100.0
          ZmaxRf = tmp2/100.0
          NzIntvRf = tmpint
          if(tmpint.ne.NzIntvRf) then
            print*,"input data wrong in Z: ",NzIntvRf,tmpint
            stop
          endif
          ! the input range units are cm
          read(14,*,end=33)tmp1
          read(14,*,end=33)tmp1,tmp2,tmpint
          RminRf = tmp1/100.0
          RmaxRf = tmp2/100.0
          NrIntvRf = tmpint
          if(tmpint.ne.NrIntvRf) then
            print*,"input data wrong in R: ",NrIntvRf,tmpint
            stop
          endif
          !print*,"Nz: ",NzIntvRf,ZminRf,ZmaxRf
          !print*,"Nr: ",NrIntvRf,RminRf,RmaxRf
        endif
33      continue
        call MPI_BCAST(NrIntvRf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(NzIntvRf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        deallocate(ezdata)
        deallocate(erdata)
        deallocate(btdata)
        allocate(ezdata(NzIntvRf+1,NrIntvRf+1))
        allocate(erdata(NzIntvRf+1,NrIntvRf+1))
        allocate(btdata(NzIntvRf+1,NrIntvRf+1))

        if(myrank.eq.0) then
          n = 0
50        continue
            if(mod(n,2).eq.0) then
              read(14,*,end=77)tmp1,tmp2,tmp3
              nn = n/2+1
              j  = (nn-1)/(NzIntvRf+1) + 1
              i = mod((nn-1),NzIntvRf+1) + 1
              ezdata(i,j) = tmp1
              erdata(i,j) = tmp2
              n = n + 1
              write(15,100)float(i-1),ezdata(i,j)
            else
              read(14,*,end=77)tmp1
              nn = (n+1)/2
              j  = (nn-1)/(NzIntvRf+1) + 1
              i = mod((nn-1),NzIntvRf+1) + 1
              btdata(i,j) = tmp1
              n = n + 1
            endif
          goto 50
77        continue
          close(14)
          close(15)
          Ndatalc = n/2
          !print*,"Ndata in 0: ",Ndatalc
        endif
100     format(2(1x,1pe15.8))

        call MPI_BCAST(Ndatalc,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZminRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZmaxRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(RminRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(RmaxRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ezdata,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(erdata,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(btdata,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndatalc
        end subroutine read2_Data

        !readin the Fourier coefficients of the RF field from files 
        !"rfdatax or rfdataxx or rfdataxxx".
        subroutine read3_Data(ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1
        character(256) :: fileid

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          write(fileid,*), ifile
          open(14, file=trim(dataclass_dir)//'rfdata'//trim(adjustl(fileid)), status='old', iostat=iodatafile)
        endif
        
        call MPI_BCAST(iodatafile,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (iodatafile.ne.0) then
          if(myrank.eq.0) print*, 'Data.f90 open(): RF data file (',&
                           trim(dataclass_dir)//'rfdata'//trim(adjustl(fileid)),&
                           ') does not found.'
          close(14)
          return
        endif
        
        if(myrank.eq.0) then
          n = 0
50        continue
            read(14,*,end=77)tmp1
            n = n + 1
            Fcoef(n) = tmp1
          goto 50
77        continue
          close(14)
          Ndata = n
        endif

        call MPI_BCAST(Ndata,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Fcoef,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        !enforce the total # of coefficients to be an odd number.
        if(mod(Ndata,2).eq.0) then
          Ndata = Ndata + 1
          Fcoef(Ndata) = 0.0
        endif

        !print*,"Ndata: ",Ndata
        end subroutine read3_Data

        ! read in discrete Ex(x,y,z), Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z),
        ! Bz(x,y,z) rf data from
        ! files "1Tx.T7 or 1Txx.T7 or 1Txxx.T7. Here, the grid
        ! is uniform in x, y and z.
        subroutine read4_Data(ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,j,k,n,m,tmpint,ios, tmpsize
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0,tmp5,tmp6, nan, zero
        double precision, dimension(6,6) :: tmpls
        double precision, allocatable, dimension(:,:) :: tmpgrid, tmpgrid1, tmpgrid2, tmpgrid3
        character(256) :: fileid

        mu0 = 4*2*asin(1.0)*1.0e-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          write(fileid,*), ifile
          open(14, file=trim(dataclass_dir)//'1T'//trim(adjustl(fileid))//'.T7', status='old' ,iostat=iodatafile)
        endif
        
        call MPI_BCAST(iodatafile,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (iodatafile.ne.0) then
          if(myrank.eq.0) print*, 'Data.f90 open(): Field file (',&
                           trim(dataclass_dir)//'1T'//trim(adjustl(fileid))//'.T7',&
                           ') does not found.'
          close(14)
          return
        endif
        
        if(myrank.eq.0) then
          ! the input range units are m
          read(14,*,iostat=iodatafile)tmp1,tmp2,tmpint
          XminRfg = tmp1
          XmaxRfg = tmp2
          NxIntvRfg = tmpint
          read(14,*,iostat=iodatafile)tmp1,tmp2,tmpint
          YminRfg = tmp1
          YmaxRfg = tmp2
          NyIntvRfg = tmpint
          read(14,*,iostat=iodatafile)tmp1,tmp2,tmpint
          ZminRfg = tmp1
          ZmaxRfg = tmp2
          NzIntvRfg = tmpint

        endif
        
        call MPI_BCAST(iodatafile,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (iodatafile.ne.0) then
            if(myrank.eq.0) print*, 'Data.f90 read(): Field file (',&
                           trim(dataclass_dir)//'1T'//trim(adjustl(fileid))//'.T7',&
                           ') has incorrect format. (1)'
            close(14)
            return
        endif
        
        call MPI_BCAST(NxIntvRfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(NyIntvRfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(NzIntvRfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(XminRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(XmaxRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(YminRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(YmaxRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZminRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZmaxRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        tmpsize = (NxIntvRfg+1)*(NyIntvRfg+1)*(NzIntvRfg+1)
        
        if(myrank.eq.0) then
            allocate(tmpgrid(6, (tmpsize+2)*3))
            allocate(tmpgrid1(6, tmpsize))
            allocate(tmpgrid2(6, tmpsize))
            allocate(tmpgrid3(6, tmpsize))
            tmpgrid = 0.0
            tmpgrid1 = 0.0
            tmpgrid2 = 0.0
            tmpgrid3 = 0.0
            ios = 0

            if (tmpsize .lt. 6) then
                zero = 0.0
                nan = zero/zero
                tmpgrid = nan
                n = 1
                do while(ios.eq.0)
                    tmpls = nan
                    read(14,*,iostat=ios) tmpls
                    tmpgrid(:,n:n+5) = tmpls
                    n = n+6
                enddo
                
                do i = 1,(tmpsize+2)*3
                    if (tmpgrid(1,i).eq.tmpgrid(1,i)) m = i
                enddo
                
                if (tmpsize .eq. m) then
                    Nmultigrid = 1
                    tmpgrid1 = tmpgrid(1:6,1:tmpsize)
                else if (tmpsize*2 .eq. m) then
                    Nmultigrid = 2
                    tmpgrid1 = tmpgrid(1:6,1:tmpsize*2:2)
                    tmpgrid2 = tmpgrid(1:6,2:tmpsize*2:2)
                else if (tmpsize*3 .eq. m) then
                    Nmultigrid = 3
                    tmpgrid1 = tmpgrid(1:6,1:tmpsize*3:3)
                    tmpgrid2 = tmpgrid(1:6,2:tmpsize*3:3)
                    tmpgrid3 = tmpgrid(1:6,3:tmpsize*3:3)
                else
                    print*, 'Data.f90 read(): Field file (',&
                            trim(dataclass_dir)//'1T'//trim(adjustl(fileid))//'.T7',&
                            ') has incorrect format. (2)'
                    close(14)
                    return
                endif
                
            else
                n = 1
                do while(ios.eq.0) 
                    read(14,*,iostat=ios) tmpls
                    tmpgrid(:,n:n+5) = tmpls
                    n = n+6
                enddo
                
                if ((tmpsize/6) .eq. ((n-1)/6-1)) then
                    Nmultigrid = 1
                    tmpgrid1 = tmpgrid(1:6,1:tmpsize)
                else if ((tmpsize/3) .eq. ((n-1)/6-1)) then
                    Nmultigrid = 2
                    tmpgrid1 = tmpgrid(1:6,1:tmpsize*2:2)
                    tmpgrid2 = tmpgrid(1:6,2:tmpsize*2:2)
                else if((tmpsize/2) .eq. ((n-1)/6-1)) then
                    Nmultigrid = 3
                    tmpgrid1 = tmpgrid(1:6,1:tmpsize*3:3)
                    tmpgrid2 = tmpgrid(1:6,2:tmpsize*3:3)
                    tmpgrid3 = tmpgrid(1:6,3:tmpsize*3:3)
                else
                    print*, 'Data.f90 read(): Field file (',&
                            trim(dataclass_dir)//'1T'//trim(adjustl(fileid))//'.T7',&
                            ') has incorrect format. (2)'
                    close(14)
                    return
                endif
            endif
            
            close(14)
            deallocate(tmpgrid)
        endif 

        call MPI_BCAST(Nmultigrid,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        if (allocated(Exgrid)) deallocate(Exgrid)
        if (allocated(Eygrid)) deallocate(Eygrid)
        if (allocated(Ezgrid)) deallocate(Ezgrid)
        if (allocated(Bxgrid)) deallocate(Bxgrid)
        if (allocated(Bygrid)) deallocate(Bygrid)
        if (allocated(Bzgrid)) deallocate(Bzgrid)
        allocate(Exgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Eygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Ezgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bxgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bzgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        Exgrid = 0.0
        Eygrid = 0.0
        Ezgrid = 0.0
        Bxgrid = 0.0
        Bygrid = 0.0
        Bzgrid = 0.0
        
        if (Nmultigrid.eq.2 .or.Nmultigrid.eq.3) then
            if (allocated(Exgrid2)) deallocate(Exgrid2)
            if (allocated(Eygrid2)) deallocate(Eygrid2)
            if (allocated(Ezgrid2)) deallocate(Ezgrid2)
            if (allocated(Bxgrid2)) deallocate(Bxgrid2)
            if (allocated(Bygrid2)) deallocate(Bygrid2)
            if (allocated(Bzgrid2)) deallocate(Bzgrid2)
            allocate(Exgrid2(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Eygrid2(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Ezgrid2(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Bxgrid2(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Bygrid2(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Bzgrid2(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            Exgrid2 = 0.0
            Eygrid2 = 0.0
            Ezgrid2 = 0.0
            Bxgrid2 = 0.0
            Bygrid2 = 0.0
            Bzgrid2 = 0.0
        endif
        
        if (Nmultigrid.eq.3) then
            if (allocated(Exgrid3)) deallocate(Exgrid3)
            if (allocated(Eygrid3)) deallocate(Eygrid3)
            if (allocated(Ezgrid3)) deallocate(Ezgrid3)
            if (allocated(Bxgrid3)) deallocate(Bxgrid3)
            if (allocated(Bygrid3)) deallocate(Bygrid3)
            if (allocated(Bzgrid3)) deallocate(Bzgrid3)
            allocate(Exgrid3(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Eygrid3(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Ezgrid3(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Bxgrid3(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Bygrid3(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            allocate(Bzgrid3(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
            Exgrid3 = 0.0
            Eygrid3 = 0.0
            Ezgrid3 = 0.0
            Bxgrid3 = 0.0
            Bygrid3 = 0.0
            Bzgrid3 = 0.0
        endif
        
        if(myrank.eq.0) then
            do n=1,tmpsize
                k = (n-1)/((NxIntvRfg+1)*(NyIntvRfg+1))+1
                j = (n-1-(k-1)*(NxIntvRfg+1)*(NyIntvRfg+1))/(NxIntvRfg+1) + 1
                i = n - (k-1)*(NxIntvRfg+1)*(NyIntvRfg+1) - (j-1)*(NxIntvRfg+1)
                Exgrid(i,j,k) = tmpgrid1(1,n)
                Eygrid(i,j,k) = tmpgrid1(2,n)
                Ezgrid(i,j,k) = tmpgrid1(3,n)
                Bxgrid(i,j,k) = tmpgrid1(4,n)
                Bygrid(i,j,k) = tmpgrid1(5,n)
                Bzgrid(i,j,k) = tmpgrid1(6,n)

                if (Nmultigrid.eq.2 .or.Nmultigrid.eq.3) then
                    Exgrid2(i,j,k) = tmpgrid2(1,n)
                    Eygrid2(i,j,k) = tmpgrid2(2,n)
                    Ezgrid2(i,j,k) = tmpgrid2(3,n)
                    Bxgrid2(i,j,k) = tmpgrid2(4,n)
                    Bygrid2(i,j,k) = tmpgrid2(5,n)
                    Bzgrid2(i,j,k) = tmpgrid2(6,n)
                endif
                
                if (Nmultigrid.eq.3) then
                    Exgrid3(i,j,k) = tmpgrid3(1,n)
                    Eygrid3(i,j,k) = tmpgrid3(2,n)
                    Ezgrid3(i,j,k) = tmpgrid3(3,n)
                    Bxgrid3(i,j,k) = tmpgrid3(4,n)
                    Bygrid3(i,j,k) = tmpgrid3(5,n)
                    Bzgrid3(i,j,k) = tmpgrid3(6,n)
                endif
            enddo
            
            deallocate(tmpgrid1, tmpgrid2, tmpgrid3)
        endif

        call MPI_BCAST(tmpsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Exgrid,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Eygrid,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Ezgrid,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Bxgrid,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Bygrid,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Bzgrid,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        if (Nmultigrid.eq.2 .or.Nmultigrid.eq.3) then
            call MPI_BCAST(Exgrid2,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Eygrid2,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Ezgrid2,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Bxgrid2,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Bygrid2,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Bzgrid2,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif

        if (Nmultigrid.eq.3) then
            call MPI_BCAST(Exgrid3,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Eygrid3,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Ezgrid3,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Bxgrid3,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Bygrid3,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Bzgrid3,tmpsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        end subroutine read4_Data
        
        subroutine read_RFQline(ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        integer :: myrank, ierr, ios, ietot, imtot, i, j, k, tmp1, tmp2, tmp3
        character(256) :: fileid
        double precision :: ztmp, zpos
        double precision, dimension(13) :: tmpdata
        integer, allocatable, dimension(:) :: tmpsegs

        if (allocated(Rfqzlen)) deallocate(Rfqzlen)
        if (allocated(Rfqzpos)) deallocate(Rfqzpos)
        if (allocated(Rfqnseg)) deallocate(Rfqnseg)
        if (allocated(Rfqpstp)) deallocate(Rfqpstp)
        if (allocated(Rfqltyp)) deallocate(Rfqltyp)
        if (allocated(Rfqpara)) deallocate(Rfqpara)
        if (allocated(Rfqsegs)) deallocate(Rfqsegs)
        if (allocated(Rfqprof)) deallocate(Rfqprof)
        if (allocated(Rfqppos)) deallocate(Rfqppos)
        
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          write(fileid,*), ifile
          open(14, file=trim(dataclass_dir)//'rfqline'//trim(adjustl(fileid)), status='old' ,iostat=iodatafile)
        endif

        call MPI_BCAST(iodatafile,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (iodatafile.ne.0) then
          if(myrank.eq.0) print*, 'Data.f90 open(): RFQ file (',&
                           trim(dataclass_dir)//'rfqline'//trim(adjustl(fileid)),&
                           ') does not found.'
          close(14)
          return
        endif

        if (myrank.eq.0) then
            ietot = 0
            imtot = 0
            ios = 0
            do while(ios.ge.0)
                ios = 5010
                do while (ios.eq.5010)
                    read(14,*,iostat=ios) ztmp
                enddo
                backspace(14)
                ztmp = 0.0
                tmp3 = 0
                read(14,*,iostat=ios) ztmp, tmp1, tmp2, tmp3
                if (tmp3.eq.-2) then
                    imtot = imtot + 1
                else if (tmp3.ge.0) then
                    ietot = ietot + 1
                endif
            enddo
            rewind(14)
        endif

        call MPI_BCAST(ietot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        
        allocate(Rfqzlen(ietot))
        allocate(Rfqzpos(ietot))
        allocate(Rfqnseg(ietot))
        allocate(Rfqpstp(ietot))
        allocate(Rfqltyp(ietot))
        allocate(Rfqpara(14,ietot))
        allocate(tmpsegs(ietot+1))

        allocate(Rfqprof(2,imtot))
        allocate(Rfqppos(ietot))
        
        Rfqntot = ietot
        
        Rfqzlen = 0.0
        Rfqnseg = 0
        Rfqpstp = 0
        Rfqltyp = 0
        Rfqpara = 0.0
        tmpsegs = 0
        Rfqtseg = 0
        
        Rfqprof = 0
        Rfqppos = 0
        
        if (myrank.eq.0) then
            tmpdata = 0.0
            i = 1
            j = 1
            k = 1
            ios = 0
            do while(ios.ge.0)
                ios = 5010
                do while (ios.eq.5010)
                    read(14,*,iostat=ios) ztmp
                enddo
                backspace(14)
                ztmp = 0.0
                tmp1 = 0
                tmp3 = 0
                read(14,*,iostat=ios) ztmp, tmp1, tmp2, tmp3, tmpdata
                if (tmp3.ge.0) then
                    Rfqzlen(i) = ztmp
                    Rfqnseg(i) = tmp1
                    Rfqpstp(i) = tmp2
                    Rfqltyp(i) = tmp3
                    Rfqpara(1:13,i) = tmpdata
                    Rfqtseg = Rfqtseg + tmp1
                    tmpsegs(i+1) = Rfqtseg
                    if (any(tmp3.eq.(/0,1/))) then
                        Rfqpara(14,i) = k
                        k = k + 1
                    endif
                    i = i + 1
                else if (tmp3.eq.-2) then
                    Rfqprof(1,j) = tmp2
                    Rfqprof(2,j) = int(tmpdata(2)+0.1)
                    Rfqppos(i) = j
                    j = j +1
                endif
            enddo
            close(14)
        endif
        
        call MPI_BCAST(Rfqtseg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Rfqzlen,ietot,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Rfqnseg,ietot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Rfqpstp,ietot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Rfqltyp,ietot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Rfqpara,14*ietot,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(tmpsegs,ietot+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Rfqprof,2*imtot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Rfqppos,ietot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        
        allocate(Rfqsegs(Rfqtseg))
        zpos = 0.0
        do i = 1, ietot
            Rfqzpos(i) = zpos
            Rfqsegs(tmpsegs(i)+1:tmpsegs(i+1)) = i
            zpos = zpos + Rfqzlen(i)
        enddo
        
        deallocate(tmpsegs)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        Rfqid = ifile

        end subroutine read_RFQline

      end module Dataclass
