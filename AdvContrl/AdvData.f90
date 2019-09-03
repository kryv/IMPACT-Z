module AdvDataclass
    use Dataclass

    implicit none

    integer :: file_num = 1
    integer, parameter :: nmaxfile = 100
    integer, dimension(nmaxfile) :: file_list
    integer :: current_file = 0

    type data_grid
        integer :: Ndata
        double precision,allocatable,dimension(:) :: z,e,ep,epp
        double precision,allocatable,dimension(:) :: Fc
        integer :: Nx, Ny, Nz, Nmg
        double precision :: Xmax,Xmin,Ymax,Ymin,Zmax,Zmin
        double precision,allocatable,dimension(:,:,:) :: Ex,Ey,Ez,Bx,By,Bz,&
                              Ex2,Ey2,Ez2,Bx2,By2,Bz2,Ex3,Ey3,Ez3,Bx3,By3,Bz3

        integer :: qntot, qid, qtseg
        integer,allocatable,dimension(:) :: qnseg,qpstp,qltyp,qsegs,qppos
        integer,allocatable,dimension(:,:) :: qprof
        double precision,allocatable,dimension(:) :: qzlen
        double precision,allocatable,dimension(:,:) :: qpara
    end type

    type(data_grid) :: dgrid(nmaxfile)

contains
    subroutine init_AdvData()
        !Inititalize AdvData module
        implicit none
        file_num = 1
        file_list = 0
    end subroutine init_AdvData

    subroutine destruct_AdvData()
        !Destruct storaged grid data in AdvData
        implicit none
        integer :: i
        file_num = 1
        file_list = 0
        do i = 1, nmaxfile
            if (allocated(dgrid(i)%z)) deallocate(dgrid(i)%z)
            if (allocated(dgrid(i)%e)) deallocate(dgrid(i)%e)
            if (allocated(dgrid(i)%ep)) deallocate(dgrid(i)%ep)
            if (allocated(dgrid(i)%epp)) deallocate(dgrid(i)%epp)
            if (allocated(dgrid(i)%Fc)) deallocate(dgrid(i)%Fc)
            if (allocated(dgrid(i)%Ex)) deallocate(dgrid(i)%Ex)
            if (allocated(dgrid(i)%Ey)) deallocate(dgrid(i)%Ey)
            if (allocated(dgrid(i)%Ez)) deallocate(dgrid(i)%Ez)
            if (allocated(dgrid(i)%Bx)) deallocate(dgrid(i)%Bx)
            if (allocated(dgrid(i)%By)) deallocate(dgrid(i)%By)
            if (allocated(dgrid(i)%Bz)) deallocate(dgrid(i)%Bz)
            if (allocated(dgrid(i)%Ex2)) deallocate(dgrid(i)%Ex2)
            if (allocated(dgrid(i)%Ey2)) deallocate(dgrid(i)%Ey2)
            if (allocated(dgrid(i)%Ez2)) deallocate(dgrid(i)%Ez2)
            if (allocated(dgrid(i)%Bx2)) deallocate(dgrid(i)%Bx2)
            if (allocated(dgrid(i)%By2)) deallocate(dgrid(i)%By2)
            if (allocated(dgrid(i)%Bz2)) deallocate(dgrid(i)%Bz2)
            if (allocated(dgrid(i)%Ex3)) deallocate(dgrid(i)%Ex3)
            if (allocated(dgrid(i)%Ey3)) deallocate(dgrid(i)%Ey3)
            if (allocated(dgrid(i)%Ez3)) deallocate(dgrid(i)%Ez3)
            if (allocated(dgrid(i)%Bx3)) deallocate(dgrid(i)%Bx3)
            if (allocated(dgrid(i)%By3)) deallocate(dgrid(i)%By3)
            if (allocated(dgrid(i)%Bz3)) deallocate(dgrid(i)%Bz3)
            if (allocated(dgrid(i)%qnseg)) deallocate(dgrid(i)%qnseg)
            if (allocated(dgrid(i)%qpstp)) deallocate(dgrid(i)%qpstp)
            if (allocated(dgrid(i)%qltyp)) deallocate(dgrid(i)%qltyp)
            if (allocated(dgrid(i)%qsegs)) deallocate(dgrid(i)%qsegs)
            if (allocated(dgrid(i)%qzlen)) deallocate(dgrid(i)%qzlen)
            if (allocated(dgrid(i)%qpara)) deallocate(dgrid(i)%qpara)
            if (allocated(dgrid(i)%qprof)) deallocate(dgrid(i)%qprof)
            if (allocated(dgrid(i)%qppos)) deallocate(dgrid(i)%qppos)

        enddo
    end subroutine destruct_AdvData

    subroutine store_Data(flagmap, bitype, ifile)
        !Read and store grid data to memory
        implicit none
        include 'mpif.h'
        integer, intent(in) :: flagmap, bitype, ifile
        integer :: myrank, ierr

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if (ALL(file_list.ne.ifile)) then
            iodatafile = 0
            if (flagmap.eq.1) then
                if (bitype.eq.106 .or. bitype.eq.110 .or. bitype.eq.114) then
                    if(myrank.eq.0) print*, 'Element type: ', bitype, ' does not support in linear map integrator.'
                    return
                endif

                call read1_Data(ifile)
                dgrid(file_num)%Ndata = Ndata
                dgrid(file_num)%z = zdat
                dgrid(file_num)%e = edat
                dgrid(file_num)%ep = epdat
                dgrid(file_num)%epp = eppdat
            else if (bitype.eq.106) then
                call read_RFQline(ifile)
                dgrid(file_num)%qntot = Rfqntot
                dgrid(file_num)%qid = Rfqid
                dgrid(file_num)%qtseg =Rfqtseg
                dgrid(file_num)%qzlen = Rfqzlen
                dgrid(file_num)%qnseg = Rfqnseg
                dgrid(file_num)%qpstp = Rfqpstp
                dgrid(file_num)%qltyp = Rfqltyp
                dgrid(file_num)%qsegs = Rfqsegs
                dgrid(file_num)%qpara = Rfqpara
                dgrid(file_num)%qprof = Rfqprof
                dgrid(file_num)%qppos = Rfqppos
            else if (bitype.eq.110 .or. bitype.eq.114) then
                call read4_Data(ifile)
                dgrid(file_num)%Nx = NxIntvRfg
                dgrid(file_num)%Ny = NyIntvRfg
                dgrid(file_num)%Nz = NzIntvRfg
                dgrid(file_num)%Nmg = Nmultigrid
                dgrid(file_num)%Xmin = XminRfg
                dgrid(file_num)%Xmax = XmaxRfg
                dgrid(file_num)%Ymin = YminRfg
                dgrid(file_num)%Ymax = YmaxRfg
                dgrid(file_num)%Zmin = ZminRfg
                dgrid(file_num)%Zmax = ZmaxRfg
                dgrid(file_num)%Ex = Exgrid
                dgrid(file_num)%Ey = Eygrid
                dgrid(file_num)%Ez = Ezgrid
                dgrid(file_num)%Bx = Bxgrid
                dgrid(file_num)%By = Bygrid
                dgrid(file_num)%Bz = Bzgrid

                if (Nmultigrid.eq.2 .or. Nmultigrid.eq.3) then
                    dgrid(file_num)%Ex2 = Exgrid2
                    dgrid(file_num)%Ey2 = Eygrid2
                    dgrid(file_num)%Ez2 = Ezgrid2
                    dgrid(file_num)%Bx2 = Bxgrid2
                    dgrid(file_num)%By2 = Bygrid2
                    dgrid(file_num)%Bz2 = Bzgrid2
                endif

                if (Nmultigrid.eq.3) then
                    dgrid(file_num)%Ex3 = Exgrid3
                    dgrid(file_num)%Ey3 = Eygrid3
                    dgrid(file_num)%Ez3 = Ezgrid3
                    dgrid(file_num)%Bx3 = Bxgrid3
                    dgrid(file_num)%By3 = Bygrid3
                    dgrid(file_num)%Bz3 = Bzgrid3
                endif
            else
                call read3_Data(ifile)
                dgrid(file_num)%Ndata = Ndata
                dgrid(file_num)%Fc = Fcoef
            endif
            if (iodatafile.ne.0) return
            file_list(file_num) = ifile
            current_file = ifile
            file_num = file_num + 1
        endif
    end subroutine store_Data

    subroutine use_Data(flagmap, bitype, ifile)
        !Copy grid data from memory storage to calculation area
        implicit none
        integer, intent(in) :: flagmap, bitype, ifile
        integer :: i,fn

        fn = 0

        do i = 1, nmaxfile
            if (file_list(i).eq.ifile) fn = i
        enddo

        if (fn.eq.0) then
            call store_Data(flagmap, bitype, ifile)
        else if (current_file.ne.ifile) then
            if (flagmap.eq.1) then
                Ndata = dgrid(fn)%Ndata
                zdat = dgrid(fn)%z
                edat = dgrid(fn)%e
                epdat = dgrid(fn)%ep
                eppdat = dgrid(fn)%epp
            else if (bitype.eq.106) then
                Rfqntot = dgrid(fn)%qntot
                Rfqid = dgrid(fn)%qid
                Rfqtseg = dgrid(fn)%qtseg
                Rfqzlen = dgrid(fn)%qzlen
                Rfqnseg = dgrid(fn)%qnseg
                Rfqpstp = dgrid(fn)%qpstp
                Rfqltyp = dgrid(fn)%qltyp
                Rfqsegs = dgrid(fn)%qsegs
                Rfqpara = dgrid(fn)%qpara
                Rfqprof = dgrid(fn)%qprof
                Rfqppos = dgrid(fn)%qppos
            else if (bitype.eq.110 .or. bitype.eq.114) then
                NxIntvRfg = dgrid(fn)%Nx
                NyIntvRfg = dgrid(fn)%Ny
                NzIntvRfg = dgrid(fn)%Nz
                Nmultigrid = dgrid(fn)%Nmg
                XminRfg = dgrid(fn)%Xmin
                XmaxRfg = dgrid(fn)%Xmax
                YminRfg = dgrid(fn)%Ymin
                YmaxRfg = dgrid(fn)%Ymax
                ZminRfg = dgrid(fn)%Zmin
                ZmaxRfg = dgrid(fn)%Zmax
                Exgrid = dgrid(fn)%Ex
                Eygrid = dgrid(fn)%Ey
                Ezgrid = dgrid(fn)%Ez
                Bxgrid = dgrid(fn)%Bx
                Bygrid = dgrid(fn)%By
                Bzgrid = dgrid(fn)%Bz

                if (Nmultigrid.eq.2 .or. Nmultigrid.eq.3) then
                    Exgrid2 = dgrid(fn)%Ex2
                    Eygrid2 = dgrid(fn)%Ey2
                    Ezgrid2 = dgrid(fn)%Ez2
                    Bxgrid2 = dgrid(fn)%Bx2
                    Bygrid2 = dgrid(fn)%By2
                    Bzgrid2 = dgrid(fn)%Bz2
                endif

                if (Nmultigrid.eq.3) then
                    Exgrid3 = dgrid(fn)%Ex3
                    Eygrid3 = dgrid(fn)%Ey3
                    Ezgrid3 = dgrid(fn)%Ez3
                    Bxgrid3 = dgrid(fn)%Bx3
                    Bygrid3 = dgrid(fn)%By3
                    Bzgrid3 = dgrid(fn)%Bz3
                endif
            else
                Ndata = dgrid(fn)%Ndata
                Fcoef = dgrid(fn)%Fc
            endif
            current_file = ifile
        endif

    end subroutine use_Data

    subroutine store_RFQline
        implicit none


    end subroutine store_RFQline

end module AdvDataclass