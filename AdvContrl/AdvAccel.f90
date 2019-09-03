module AdvAccelclass
    use Pgrid2dclass
    use CompDomclass
    use FieldQuantclass
    use BeamLineElemclass
    use Ptclmgerclass
    use BeamBunchclass
    use Timerclass
    use AdvInputclass
    use AdvOutputclass
    use Dataclass
    use PhysConstclass
    use NumConstclass
    use Distributionclass
    use Besselclass
    use AdvDataclass
    use MTrndclass
    use Messaging

    implicit none

    !1d logical processor array.
    type (Pgrid2d), private :: Grid2d

    !beam particle object and array.
    type (BeamBunch), private :: Bpts

    !beam charge density and field potential arrays.
    type (FieldQuant), private :: Potential

    !geometry object.
    type (CompDom), private :: Ageom

    !beam line element array.
    type (BPM),target,private,dimension(Nbpmmax) :: beamln0
    type (DriftTube),target,private,dimension(Ndriftmax) :: beamln1
    type (Quadrupole),target,private,dimension(Nquadmax) :: beamln2
    type (DTL),target,private,dimension(Ndtlmax) :: beamln3
    type (CCDTL),target,private,dimension(Nccdtlmax) :: beamln4
    type (CCL),target,private,dimension(Ncclmax) :: beamln5
    type (SC),target,private,dimension(Nscmax) :: beamln6
    type (ConstFoc),target,private,dimension(Ncfmax) :: beamln7
    type (SolRF),target,private,dimension(Nslrfmax) :: beamln8
    type (Sol),target,private,dimension(Nslmax) :: beamln9
    type (Dipole),target,private,dimension(Ndipolemax) :: beamln10
    type (EMfld),target,private,dimension(Ncclmax) :: beamln11
    type (Multipole),target,private,dimension(Nquadmax) :: beamln12
    type (RFQ),target,private,dimension(Ncclmax) :: beamln13
    type (Sol2),target,dimension(Nslmax) :: beamln14
    type (BeamLineElem),private,dimension(Nblemtmax)::Blnelem

    ! number of particles of charge state.
    integer, dimension(100) :: Nptlist0

    ! current list of charge state, charge/mass list of charge state.
    double precision, dimension(100) :: Currlist0, Qmcclist0

    interface configure
        module procedure get_Beamln, set_Beamln1, set_Beamln2
    end interface

contains
    subroutine load_Input(filename)
        implicit none
        include 'mpif.h'
        character, intent(in) :: filename*256
        double precision :: t0

        ! initialize Timer.
        call construct_Timer(0.0d0)

        call starttime_Timer(t0)

        Nptlistin = 0
        Currlistin = 0.0
        Qmcclistin = 0.0

        ! get all global input parameters.
        call in_Input(filename, Dm, Npin, Nx, Ny, Nz, Flagbc, Flagdist, Rstartflg,&
                      Flagmap, Distparam, 21, Bcurr, Bkenergy, Bmass, Bchargein,&
                      Bfreq, Xwallrad, Ywallrad, Perdlen, Nblem, Npcol, Nprow, Flagerr,&
                      Flagdiag, Flagsubstep, Phsini, Nchrgin, Nptlistin, Currlistin,&
                      Qmcclistin,Outputflag)
        if (ioinpfile.ne.0) return

        call construct_Beamline(1)

        !construct Constants class.
        call construct_PhysConst(Bfreq)

        t_init = t_init + elapsedtime_Timer(t0)

    end subroutine load_Input

    subroutine construct_Beamline(flaginput)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: flaginput
        integer :: i, ierr
        double precision :: zlength

        if (flaginput.eq.1) then
        ! construct beam line elements.

            if (allocated(Blnlen)) deallocate(Blnlen)
            if (allocated(Blnseg)) deallocate(Blnseg)
            if (allocated(Blnstp)) deallocate(Blnstp)
            if (allocated(Blntyp)) deallocate(Blntyp)
            if (allocated(Blnparams)) deallocate(Blnparams)

            allocate(Blnlen(Nblem),Blnseg(Nblem),Blnstp(Nblem),Blntyp(Nblem))
            allocate(Blnparams(28,Nblem))

            call in_Input(Nblem,Blnlen,Blnseg,Blnstp,Blntyp,Blnparams)
        endif

        if (allocated(Ielems)) deallocate(Ielems)

        allocate(Ielems(Nblem))

        zlength = 0.0

        iccl = 0
        iccdtl = 0
        idtl = 0
        isc = 0
        idr = 0
        iqr = 0
        ibpm = 0
        icf = 0
        islrf = 0
        isl = 0
        isl2 = 0
        idipole = 0
        iemfld = 0
        imultpole = 0
        irfq = 0

        Idiag = 0
        Idiagbpm = 0
        Idstout = 0

        do i = 1, Nblem
            call configure_Beamln(i,1,Blntyp(i),Blnseg(i),Blnstp(i),Blnlen(i))
        enddo

        call count_diags(1,Nblem)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    end subroutine construct_Beamline

    subroutine count_diags(ini,fin)
        implicit none
        integer, intent(in) :: ini, fin
        integer :: i,n, etyp, eseg, nfile, nlen
        double precision :: zlength, rfile

        if (ini.le.0) then
            Lcini = 1
        else
            Lcini = ini
        endif

        if (fin.le.0) then
            Lcfin = Nblem
        else
            Lcfin = fin
        endif

        nlen = Lcfin-Lcini+2

        Idiag = 0
        Idiagbpm = 0
        Idstout = 0

        zlength = 0.0

        if (ini .eq. 1 .and. fin .eq. Nblem) then
            if (allocated(Blnpos)) deallocate(Blnpos)
            if (allocated(Blnidg)) deallocate(Blnidg)
            allocate(Blnpos(nlen), Blnidg(nlen))
            Blnpos = zlength
            Blnidg = 0
        endif

        nfile = 0

        n = 1
        do i = Lcini, Lcfin
            etyp = Blntyp(i)
            eseg = Blnseg(i)
            if (etyp.eq.-23) then
                Idiagbpm = Idiagbpm + 1
                Idiag = Idiag + 1
            else if (any(etyp.eq.(/-11, -13, -14, -21, -25 /))) then
                Idiag = Idiag + 1
            else if (etyp.eq.-2) then
                Idstout = Idstout + 1
            else if (any(etyp.eq.(/0, 1, 2, 3, 5, 13, 101, 102, 103, 104, 105, 110, 114/))) then
                Idiag = Idiag + eseg
            else if (etyp.eq.4) then
                if (Blnparams(4,i).le.100.0) then
                    Idiag = Idiag + 1
                else
                    Idiag = Idiag + 1
                endif
            else if (etyp.eq.106) then
                if (Flagmap.eq.1) then
                    Idiag = Idiag + eseg
                else
                    call getparam_BeamLineElem(Blnelem(i),5, rfile)
                    nfile = int(rfile + 0.1)
                    call use_Data(Flagmap, etyp, nfile)
                    Idiag = Idiag + Rfqtseg
                    Idstout = Idstout + size(Rfqprof)/2
                endif
            else
            endif

            if (ini .eq. 1 .and. fin .eq. Nblem) then
                zlength = zlength + Blnlen(i)
                Blnpos(n+1) = zlength
                Blnidg(n+1) = Idiag
                n = n + 1
            endif
        enddo

    end subroutine count_diags

    subroutine configure_Beamln(i,flagconst,bitype,bnseg,bmpstp,blength)
        implicit none
        integer, intent(in) :: i, flagconst
        integer, optional, intent(in) :: bnseg,bmpstp,bitype
        double precision, optional, intent(in) :: blength
        integer :: ide
        double precision :: wk, x1, x2, xx, aa, besi1, besi2
        double precision, dimension(27) :: tmpvarr27
        double precision, dimension(30) :: tmpvarr30

        if(bitype.lt.0) then
            if (flagconst.eq.1) then
                ibpm = ibpm + 1
                Ielems(i) = ibpm
                call construct_BPM(beamln0(ibpm),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_BPM(beamln0(ide),Blnparams(1:8,i))
            Blnelem(i) = assign_BeamLineElem(beamln0(ide))
        else if(bitype.eq.0) then
            if (flagconst.eq.1) then
                idr = idr + 1
                Ielems(i) = idr
                call construct_DriftTube(beamln1(idr),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_DriftTube(beamln1(ide),Blnparams(1:2,i))
            Blnelem(i) = assign_BeamLineElem(beamln1(ide))
        else if(bitype.eq.1) then
            if (flagconst.eq.1) then
                iqr = iqr + 1
                Ielems(i) = iqr
                call construct_Quadrupole(beamln2(iqr),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_Quadrupole(beamln2(ide),Blnparams(1:9,i))
            Blnelem(i) = assign_BeamLineElem(beamln2(ide))
        else if(bitype.eq.2) then
            if (flagconst.eq.1) then
                icf = icf + 1
                Ielems(i) = icf
                call construct_ConstFoc(beamln7(icf),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_ConstFoc(beamln7(ide),Blnparams(1:5,i))
            Blnelem(i) = assign_BeamLineElem(beamln7(ide))
        else if(bitype.eq.3) then
            if (flagconst.eq.1) then
                isl = isl + 1
                Ielems(i) = isl
                call construct_Sol(beamln9(isl),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_Sol(beamln9(ide),Blnparams(1:9,i))
            Blnelem(i) = assign_BeamLineElem(beamln9(ide))
        else if(bitype.eq.4) then
            if (flagconst.eq.1) then
                idipole = idipole + 1
                Ielems(i) = idipole
                call construct_Dipole(beamln10(idipole),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_Dipole(beamln10(ide),Blnparams(1:11,i))
            Blnelem(i) = assign_BeamLineElem(beamln10(ide))
        else if(bitype.eq.5) then
            if (flagconst.eq.1) then
                imultpole = imultpole + 1
                Ielems(i) = imultpole
                call construct_Multipole(beamln12(imultpole),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_Multipole(beamln12(ide),Blnparams(1:14,i))
            Blnelem(i) = assign_BeamLineElem(beamln12(ide))
        else if(bitype.eq.13) then
            if (flagconst.eq.1) then
                isl2 = isl2 + 1
                Ielems(i) = isl2
                call construct_Sol2(beamln14(isl2),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_Sol2(beamln14(ide),Blnparams(1:11,i))
            Blnelem(i) = assign_BeamLineElem(beamln14(ide))
        else if(bitype.eq.101) then
            if (flagconst.eq.1) then
                idtl = idtl + 1
                Ielems(i) = idtl
                call construct_DTL(beamln3(idtl),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_DTL(beamln3(ide),Blnparams(1:25,i))
            Blnelem(i) = assign_BeamLineElem(beamln3(ide))
        else if(bitype.eq.102) then
            if (flagconst.eq.1) then
                iccdtl = iccdtl + 1
                Ielems(i) = iccdtl
                call construct_CCDTL(beamln4(iccdtl),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_CCDTL(beamln4(ide),Blnparams(1:11,i))
            Blnelem(i) = assign_BeamLineElem(beamln4(ide))
        else if(bitype.eq.103) then
            if (flagconst.eq.1) then
                iccl = iccl + 1
                Ielems(i) = iccl
                call construct_CCL(beamln5(iccl),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            tmpvarr27 = 0.0
            tmpvarr27(1:25) = Blnparams(1:25,i)
            call setparam_CCL(beamln5(ide),tmpvarr27)
            Blnelem(i) = assign_BeamLineElem(beamln5(ide))
        else if(bitype.eq.104) then
            if (flagconst.eq.1) then
                isc = isc + 1
                Ielems(i) = isc
                call construct_SC(beamln6(isc),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_SC(beamln6(ide),Blnparams(1:11,i))
            Blnelem(i) = assign_BeamLineElem(beamln6(ide))
        else if(bitype.eq.105) then
            if (flagconst.eq.1) then
                islrf = islrf + 1
                Ielems(i) = islrf
                call construct_SolRF(beamln8(islrf),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            call setparam_SolRF(beamln8(ide),Blnparams(1:12,i))
            Blnelem(i) = assign_BeamLineElem(beamln8(ide))
        else if(bitype.eq.106) then
            irfq = irfq + 1
            Ielems(i) = irfq
            wk = 2*PI/blength
            x1 = wk*Blnparams(6,i)
            x2 = wk*Blnparams(7,i)*Blnparams(6,i)
            besi1 = bessi0(x1)
            besi2 = bessi0(x2)
            xx = (besi1+besi2)/(Blnparams(7,i)**2*besi1+besi2)
            aa = (Blnparams(7,i)**2-1)/(Blnparams(7,i)**2*besi1+besi2)
            call construct_RFQ(beamln13(irfq),bnseg,bmpstp,bitype,blength,aa,xx)
            call setparam_RFQ(beamln13(irfq),Blnparams(1:12,i))
            Blnelem(i) = assign_BeamLineElem(beamln13(irfq))
        else if(bitype.eq.110 .or. bitype.eq.114) then
            if (flagconst.eq.1) then
                iemfld = iemfld + 1
                Ielems(i) = iemfld
                call construct_EMfld(beamln11(iemfld),bnseg,bmpstp,bitype,blength)
            endif
            ide = Ielems(i)
            tmpvarr30 = 0.0
            tmpvarr30(1:28) = Blnparams(1:28,i)
            call setparam_EMfld(beamln11(ide),tmpvarr30)
            Blnelem(i) = assign_BeamLineElem(beamln11(ide))
        else
        endif
    end subroutine configure_Beamln

    subroutine get_Beamln(i, length, seg, step, itype, varr)
        implicit none
        integer, intent(in) :: i
        integer, intent(out) :: seg, step, itype
        double precision, intent(out) :: length
        double precision, dimension(27), intent(out) :: varr
        integer :: ierr

        call check_Index(i,ierr)

        if (ierr.eq.0) then
            call getparam_BeamLineElem(Blnelem(i), length, seg, step, itype)
            varr = Blnparams(2:,i)
        endif
    end subroutine get_Beamln

    subroutine set_Beamln1(i,varr)
        implicit none
        integer, intent(in) :: i
        double precision, dimension(:), intent(in) :: varr
        integer :: seg, step, itype, ierr
        double precision :: length

        call check_Index(i,ierr)

        if (ierr.eq.0) then
            if (size(varr).gt.28) then
                print*, "Input array size is too big."
            else
                call getparam_BeamLineElem(Blnelem(i), length, seg, step, itype)
                Blnparams(1,i) = -1.0
                Blnparams(2:,i) = varr
                call configure_Beamln(i,0,itype)
            endif
        endif
    end subroutine set_Beamln1

    subroutine set_Beamln2(i,pid,v)
        implicit none
        integer, intent(in) :: i, pid
        double precision, intent(in) :: v
        integer :: loc, seg, step, itype, ierr
        double precision :: length

        call check_Index(i,ierr)

        if (ierr.eq.0) then
            if (pid .le. 3) then
                call getparam_BeamLineElem(Blnelem(i), length, seg, step, itype)
                select case (pid)
                    case(1)
                        length = v
                        Blnlen(i) = v
                    case(2)
                        seg = nint(v)
                        Blnseg = nint(v)
                    case(3)
                        step = nint(v)
                        Blnstp = nint(v)
                end select
                call setparam_BeamLineElem(Blnelem(i), length, seg, step, itype)
                call count_diags(1,Nblem)
            else
                loc = pid - 3
                Blnparams(1,i) = -1.0
                Blnparams(loc,i) = v
                call getparam_BeamLineElem(Blnelem(i), length, seg, step, itype)
                call configure_Beamln(i, 0, itype, seg, step, length)
            endif
        endif
    end subroutine set_Beamln2

    subroutine search_Index(stype)
        implicit none
        integer, intent(in) :: stype
        integer, allocatable, dimension(:) :: tmp
        double precision :: length
        integer :: i, cnt, seg, step, itype, ierr

        if (allocated(searchres)) deallocate(searchres)
        cnt = 0
        call check_Index(Nblem,ierr)
        if (ierr.eq.0) then
            allocate(tmp(Nblem))
            do i=1,Nblem
                call getparam_BeamLineElem(Blnelem(i), length, seg, step, itype)
                if (itype.eq.stype) then
                    cnt = cnt + 1
                    tmp(cnt) = i
                endif
            enddo
            allocate(searchres(cnt))
            searchres = tmp(1:cnt)
            deallocate(tmp)
        else
            allocate(searchres(0))
            searchres = -1
        endif

    end subroutine search_Index

    subroutine check_Index(i, ierr)
        implicit none
        integer, intent(in) :: i
        integer, intent(out) :: ierr

        ierr = 0
        if ((Nblem.le.0) .or. (.not.allocated(Ielems))) then
            print*, 'Beam line have not been constructed.'
            ierr = 1
        else if  (Nblem.ne.size(Ielems)) then
            print*, 'Beam line have not been constructed.'
            ierr = 1
        else if ((i.gt.Nblem).or.(i.le.0)) then
            print*, 'Inputted index ',i,' is out of range (', Nblem,').'
            ierr = 1
        endif

    end subroutine check_Index

    subroutine export_Setting(fname)
        implicit none
        include 'mpif.h'
        integer :: myrank, ierr
        character, intent(in) :: fname*256
        double precision :: rtlength
        integer :: i, rtseg, rtstep, rtitype, parasize
        double precision, dimension(28) :: rtvarr
        character(256) :: ffm, num

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if (myrank.eq.0) then

            open(12, file=trim(fname))
            ffm = '(2(1x,i0))'
            write(12,ffm), Npcol, Nprow
            ffm = '(5(1x,i0))'
            write(12,ffm), Dm, Npin, Flagmap, Flagerr, Flagdiag
            ffm = '(4(1x,i0),3(1x,e21.15))'
            write(12,ffm), Nx, Ny, Nz, Flagbc, Xwallrad, Ywallrad, Perdlen
            ffm = '(4(1x,i0))'
            write(12,ffm), Flagdist, Rstartflg, Flagsubstep, Nchrgin
            write(num,*) Nchrgin
            ffm = '('//trim(adjustl(num))//'(1x,i0))'
            write(12,ffm), Nptlistin(1:Nchrgin)
            ffm = '('//trim(adjustl(num))//'(1x,e21.15))'
            write(12,ffm), Currlistin(1:Nchrgin)
            write(12,ffm), Qmcclistin(1:Nchrgin)
            ffm = '(7(1x,e21.15))'
            write(12,ffm), Distparam(1:7)
            write(12,ffm), Distparam(8:14)
            write(12,ffm), Distparam(15:21)
            ffm = '(6(1x,e21.15),1x,e12.6)'
            write(12,ffm), Bcurr, Bkenergy, Bmass, Bchargein, Bfreq, Phsini, Outputflag

            do i=1,Nblem
                call get_Beamln(i, rtlength, rtseg, rtstep, rtitype, rtvarr)
                parasize = 1
                if (rtitype.lt.0) then
                    parasize = 7
                else if (rtitype.eq.0) then
                    parasize = 1
                else if (rtitype.eq.1) then
                    parasize = 8
                else if (rtitype.eq.2) then
                    parasize = 4
                else if (rtitype.eq.3) then
                    parasize = 8
                else if (rtitype.eq.4) then
                    parasize = 9
                else if (rtitype.eq.5) then
                    parasize = 13
                else if (rtitype.eq.101) then
                    parasize = 24
                else if (rtitype.eq.103) then
                    parasize = 10
                else if (rtitype.eq.104) then
                    parasize = 10
                else if (rtitype.eq.105) then
                    parasize = 11
                else if (rtitype.eq.106) then
                    parasize = 11
                else if (rtitype.eq.110) then
                    parasize = 26
                else
                endif

                write(num,*) parasize
                ffm = '(e21.15,3(1x,i0),'//trim(adjustl(num))//'(1x,e21.15),a2)'
                write(12,ffm), rtlength, rtseg, rtstep, rtitype, rtvarr(1:parasize), ' /'

            enddo
            write(12,*), ''

            close(12)

        endif

    end subroutine export_Setting

    subroutine load_AllData()
        implicit none
        integer :: i,bnseg,bmpstp,bitype,nfile
        double precision :: blength,rfile

        nfile = 0
        do i = 1, Nblem
            call getparam_BeamLineElem(Blnelem(i),blength,bnseg,bmpstp,bitype)
            if(bitype.gt.100) then
                call getparam_BeamLineElem(Blnelem(i),5,rfile)
                nfile = int(rfile + 0.1)
                call store_Data(Flagmap, bitype, nfile)
                if (iodatafile.ne.0) return
            endif
        enddo

    end subroutine load_AllData

    subroutine init_Particle(usrbcurr, usrbmass, usrbkenergy, usrphsini)
        implicit none
        include 'mpif.h'
        double precision, intent(in), optional :: usrbcurr, usrbmass, usrbkenergy, usrphsini
        integer :: myid, myidx, myidy, ierr, tmpnp
        double precision :: tmpbcurr, tmpbmass, tmpbkenergy, tmpphsini

        if (idistnp .ge. 0) then
            Bcharge = Bchargein
            Nptlist0 =  Nptlist
            Currlist0 =  Currlist
            Qmcclist0 =  Qmcclist
            Np = idistnp
            tmpnp = idistnp
        else
            Bcharge = Bchargein
            Nchrg = Nchrgin
            Nptlist =  Nptlistin
            Currlist =  Currlistin
            Qmcclist =  Qmcclistin
            Nptlist0 =  Nptlistin
            Currlist0 =  Currlistin
            Qmcclist0 =  Qmcclistin
            Np = Npin
            tmpnp = Npin
        endif

        if (present(usrbcurr) .and. idistnp .ge. 0) then
            tmpbcurr = usrbcurr
        else
            tmpbcurr = Bcurr
        endif

        if (present(usrbmass) .and. idistnp .ge. 0) then
            tmpbmass = usrbmass
        else
            tmpbmass = Bmass
        endif

        if (present(usrbkenergy) .and. idistnp .ge. 0) then
            tmpbkenergy = usrbkenergy
        else
            tmpbkenergy = Bkenergy
        endif

        if (present(usrphsini) .and. idistnp .ge. 0) then
            tmpphsini = usrphsini
        else
            tmpphsini = Phsini
        endif

        if (associated(Bpts%Pts1)) call destruct_BeamBunch(Bpts)
        if (associated(Potential%FieldQ)) call destruct_FieldQuant(Potential)
        if (associated(Ageom%LcTabrg).and.associated(Ageom%LcTabnm)) then
            call destruct_CompDom(Ageom)
        endif

        ! construct 2D logical processor Cartesian coordinate
        call construct_Pgrid2d(Grid2d, MPI_COMM_WORLD, Nprow, Npcol)

        !construct Constants class.
        call construct_PhysConst(Bfreq)

        call construct_CompDom(Ageom,Distparam,21,Flagdist,&
           Nx,Ny,Nz,Grid2d,Nprow,Npcol,Flagbc,Xwallrad,Ywallrad,Perdlen)

        ! construct BeamBunch class.
        !call construct_BeamBunch(Bpts,Bcurr,Bkenergy,Bmass,Bchargein,Npin,Phsini)
        call construct_BeamBunch(Bpts,tmpbcurr,tmpbkenergy,tmpbmass,Bchargein,tmpnp,tmpphsini)

        call getpost_Pgrid2d(Grid2d, myid, myidy, myidx)
        ! reset random number seed for init_Particle()
        call mtseed(isrd+21+myid)

        ! sample initial particle distribution.
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (idistnp .ge. 0) then
            call set_idistdata(Bpts)
        else
            if (Flagdist .eq. 19) then
                Nchrg = 100
                ! Maximum number of the charge states
                ! This will be overwrite in sample_Dist
            endif

            call sample_Dist(Bpts,Distparam,21,Flagdist,Ageom,Grid2d,Flagbc,&
                             Nchrg,Nptlist0,Qmcclist0,Currlist0)

            if (Flagdist .eq. 19) then
                Qmcclist = Qmcclist0
                Nptlist = Nptlist0
            endif

            if (iodistfile .ne. 0) return
        endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !get local particle number and mesh number on each processor.
        call getnpt_BeamBunch(Bpts,Nplocal)

        ! construct FieldQuant class objects.
        call construct_FieldQuant(Potential,Nx,Ny,Nz,Ageom,Grid2d)

        call phase_init_Store(Bpts)
        if (ninedata_id.eq.0) ninedata_id = -1

    end subroutine init_Particle

    subroutine set_idistdata(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer :: i,j,avgpts,myid,totnp,ierr,npsum,restnpt
        integer, allocatable, dimension(:) :: shead, scount

        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,totnp,ierr)

        npsum = idistnp

        call MPI_BCAST(npsum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        this%Npt = npsum
        avgpts = npsum/totnp
        restnpt = mod(this%Npt, totnp)
        allocate(shead(totnp), scount(totnp))
        shead(1) = 0
        do i = 1, totnp
            if ((i-1).lt.restnpt) then
                scount(i) = 9*(avgpts + 1)
            else
                scount(i) = 9*avgpts
            endif

            if (i.gt.1) shead(i) = shead(i-1) + scount(i-1)
        end do
        if (myid .lt. restnpt) avgpts = avgpts + 1

        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        call MPI_SCATTERV(idistd(1,1), scount, shead, MPI_DOUBLE_PRECISION,&
                          this%Pts1(1,1), 9*avgpts, MPI_DOUBLE_PRECISION, &
                          0, MPI_COMM_WORLD, ierr)

        this%Nptlocal = avgpts

    end subroutine set_idistdata

    !Run beam dynamics simulation through accelerator.
    subroutine run_Accel(ini,fin)
        implicit none
        include 'mpif.h'
        integer, intent(in), optional :: ini, fin
        integer :: i,j,bnseg,bmpstp,bitype,myid,myidx,myidy,comm2d,commcol,commrow,&
            npx,npy,totnp,ierr,nylcr,nzlcr,nbal,ibal,nstep,nfile,&
            nmod,k,ii,jj,ipt,ncl,nprof,celltype,stdcnt,nsubstep,integerSamplePeriod,&
            Flagbctmp,nev,iev,ichg,nptmp2,nptottmp,nplistlc,nchrgtot,ndiags,cid,flg
        integer, dimension(3) :: lcgrid
        integer, allocatable, dimension(:) :: modth,pydisp,lctabnmx,lctabnmy
        integer, allocatable, dimension(:,:,:) :: temptab
        double precision :: z0,z,tau1,tau2,blength,t0,t1,hy,hz,ymin,zmin,&
            zend,zedge,tmp1,tmp2,tmp3,tmp4,rfile,fncl,realSamplePeriod,&
            tg,tv,gam,piperad,piperad2,piperad0,beta0,gamma0,gambetz,rwkinq,rwkinf,&
            avgw,ri,r0,thetap,rr0,rra,cosang0,sinang0,drr,E0,El0,El,Thx,Thy,En,Th,&
            xradmin,xradmax,yradmin,yradmax,dpi,ang0,hd0,hd1,dstr1,dstr2,&
            angF,tanphiF,tanphiFb,angB,tanphiB,tanphiBb,hF,hB,qm0,qmi,deltax,&
            deltat,dstr1i,hdi,gambet,beta00,gamma00,gambet00,dk1, hgap,tmppt,&
            psi1,psi2,gamn,xrad,yrad,ma,wk,zz,qmrel,dlength,xerr,yerr,&
            anglerrx,anglerry,anglerrz,b0,gam0,gambetz0,finv,bnddeg,&
            bndfrn,bndrad,crvtr,zlfrg,zlbnd,zlength,refcen0,refcen,drefcen,&
            refang,refdz,zadv,fsync,flagv,strp0,strgam0,strqm,tmptheta,gtf
        double precision, dimension(3) :: msize,al0,ga0,epson0
        double precision, dimension(6) :: lcrange,srange,ptrange,ptref,ptarry,ptarry2
        double precision, dimension(8) :: drange
        double precision, dimension(11) :: dparam
        double precision, dimension(100) :: qmcclc
        double precision, allocatable, dimension(:,:) :: lctabrgx,lctabrgy,bessnorm,gml,tmpptcs
        double precision, allocatable, dimension(:,:,:) :: chgdens,besscoef

        type(BunchStats), allocatable, dimension(:) :: bstats

        !-------------------------------------------------------------------
        ! prepare initial parameters, allocate temporary array.
        if (present(ini) .and. ini.ge.1) then
            Lcini = ini
        else
            Lcini = 1
        endif

        if (present(fin) .and. fin.ge.1) then
            Lcfin = fin
        else
            Lcfin = Nblem
        endif

        iodatafile = 0
        iostrpfile = 0

        call getpost_Pgrid2d(Grid2d,myid,myidy,myidx)
        call getcomm_Pgrid2d(Grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(Grid2d,totnp,npy,npx)

        ! reset random number seed for run_Accel()
        call mtseed(isrd+21240+10*myid)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        allocate(lctabnmx(0:npx-1))
        allocate(lctabnmy(0:npy-1))
        allocate(lctabrgx(2,0:npx-1))
        allocate(lctabrgy(2,0:npy-1))
        allocate(temptab(2,0:npx-1,0:npy-1))
        allocate(tmpptcs(1,1))

        nbal = 5
        ibal = 0
        nstep = 0
        z = Zposini

        call count_diags(Lcini,Lcfin)
        nchrgtot = Nchrg
        qmcclc(1:Nchrg) = Qmcclist0

        if (Flagdiag.eq.1 .or. Flagdiag.eq.2) then
            ndiags = Idiag+1
        elseif (Flagdiag.eq.3 .or. Flagdiag.eq.4) then
            ndiags = Idiagbpm
        elseif (Flagdiag.eq.5 .or. Flagdiag.eq.6) then
            ndiags = Lcfin-Lcini+2
        endif

        allocate(bstats(ndiags))
        ciStats = 1

        call init_SpData(Idstout)

        if (any(Flagdiag.eq.(/1, 2, 5, 6/))) then
            call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                 Flagdiag, Outputflag, bstats(ciStats))
        endif

        allocate(chgdens(1,1,1))

        ninedata_id = -1

        !-------------------------------------------------------------------
        ! prepare for round pipe, open longitudinal
        if((Flagbc.eq.3).or.(Flagbc.eq.7)) then
            allocate(pydisp(0:npy-1))
            call getlctabnm_CompDom(Ageom,temptab)
            pydisp(0) = 0
            do i = 1, npy-1
                pydisp(i) = pydisp(i-1) + temptab(2,0,i-1)
            enddo
            call getlcmnum_CompDom(Ageom,lcgrid)
            Nxlocal = lcgrid(1)
            if(Npcol.gt.1) then
                Nylocal = lcgrid(2) + 2
            else
                Nylocal = lcgrid(2)
            endif
            if(Nprow.gt.1) then
                Nzlocal = lcgrid(3) + 2
            else
                Nzlocal = lcgrid(3)
            endif

            if(Nprow.gt.1) then
                nzlcr = Nzlocal-2
            else
                nzlcr = Nzlocal
            endif
            if(Npcol.gt.1) then
                nylcr = Nylocal-2
            else
                nylcr = Nylocal
            endif
            if(myidy.eq.(Npcol-1)) nylcr = nylcr - 1
            allocate(modth(nylcr))
                do i = 1, nylcr
                    modth(i) = (pydisp(myidy)+i-2+1)/2
                enddo
            if(myidy.eq.0) then
                modth(1) = 0
                modth(2) = (ny-1)/2
                nmod = modth(nylcr) - modth(1) + 2
            else
                nmod = modth(nylcr) - modth(1) + 1
            endif
            allocate(besscoef(lcgrid(1),lcgrid(1),nmod))
            allocate(bessnorm(lcgrid(1),nmod))
            allocate(gml(lcgrid(1),nmod))

            call getmsize_CompDom(Ageom,msize)
            call Bessprep_Bessel(msize(1),lcgrid(1),nylcr,nmod,modth,&
                        besscoef,bessnorm,gml)
        endif

        ! start looping through 'Nblem' beam line elements.
        do i = Lcini, Lcfin
            zedge = z

            call getparam_BeamLineElem(Blnelem(i),blength,bnseg,bmpstp,bitype)
            call getradius_BeamLineElem(Blnelem(i),piperad,piperad2)
            if(Flagerr.eq.0 .or. Flagmap.eq.1 .or. bitype.eq.3) then
                xradmin = -piperad
                xradmax = piperad
                yradmin = -piperad2
                yradmax = piperad2
            else
                call geterr_BeamLineElem(Blnelem(i),xerr,yerr,anglerrx,anglerry,anglerrz)
                xradmin = -piperad + xerr
                xradmax = piperad + xerr
                yradmin = -piperad2 + yerr
                yradmax = piperad2 + yerr
            endif

            nfile = 0
            if (bitype.eq.106 .and. Flagmap.eq.2) then
                call getparam_BeamLineElem(Blnelem(i),5,rfile)
                nfile = int(rfile + 0.1)
                call use_Data(Flagmap, bitype, nfile)
                if (blength .ne. 0.0) then
                    Rfqzlen = (blength/sum(Rfqzlen))*Rfqzlen
                endif
                bnseg = Rfqtseg
                tmprfqncl = 0
            else if (bitype.gt.100) then
                call starttime_Timer(t1)
                call getparam_BeamLineElem(Blnelem(i),5,rfile)
                nfile = int(rfile + 0.1)
                call use_Data(Flagmap, bitype, nfile)
                t_load3dfile = t_load3dfile + elapsedtime_Timer(t1)
            endif

            nsubstep = bmpstp
            tau1 = 0.0
            if (bitype.ge.0) tau1 = 0.5*blength/bnseg
            tau2 = blength/bnseg

            if (bitype.eq.114 .and. Flagmap.eq.2) then
                call getparam_BeamLineElem(Blnelem(i), 3, bnddeg)
                call getparam_BeamLineElem(Blnelem(i), 4, bndfrn)
                piperad0 = 0.5*piperad
                piperad2 = 0.5*piperad2
                bndrad = Pi*bnddeg/180
                crvtr = (blength-2.0*bndfrn)/bndrad
                zlfrg = bndfrn*cos(0.5*bndrad)
                zlbnd = 2.0*crvtr*sin(0.5*bndrad)
                zlength = 2.0*zlfrg + zlbnd
                zadv = 0.0
                refcen0 = 0.0
                tau2 = zlength/bnseg
                zedge = zedge - 0.5*((ZmaxRfg - ZminRfg)-zlength)
                gam = -Bpts%refptcl(6)
                call yrotation_BPM(Bpts%Pts1,Nplocal,0.5*bnddeg,gam)
            endif

            ! print out beam information using BPM
            if(bitype.eq.-1) then
                call shift_BPM(Bpts%Pts1,bitype,Nplocal,Np)
            else if(bitype.eq.-2) then
                call getparam_BeamLineElem(Blnelem(i),drange)
                realSamplePeriod = drange(2)
                integerSamplePeriod = realSamplePeriod
                !call phase_Output(bmpstp,Bpts,integerSamplePeriod)
                call phase_Store(bmpstp,Bpts,integerSamplePeriod)
            !else if(bitype.eq.-3) then
            !    call getparam_BeamLineElem(Blnelem(i),drange)
            !    call accdens1d_Output(nstep,8,Bpts,Np,drange(2),-drange(3),&
            !                          drange(3),-drange(5),drange(5))
            !else if(bitype.eq.-4) then
            !    call getparam_BeamLineElem(Blnelem(i),drange)
            !    call dens1d_Output(nstep,8,Bpts,Np,drange(2),-drange(3),&
            !                       drange(3),-drange(5),drange(5))
            !else if(bitype.eq.-5) then
            !    call getparam_BeamLineElem(Blnelem(i),drange)
            !    call dens2d_Output(nstep,8,Bpts,Np,-drange(3),drange(3),&
            !                       -drange(4),drange(4),-drange(5),drange(5),&
            !                       -drange(6),drange(6),-drange(7),drange(7),&
            !                       -drange(8),drange(8))
            !else if(bitype.eq.-6) then
            !    call dens3d_Output(nstep,8,Bpts,Np,-drange(3),drange(3),&
            !                       -drange(5),drange(5),-drange(7),drange(7))
            !else if(bitype.eq.-7) then
            !    call outpoint_Output(myid+31,Bpts,z,i,j,npx,npy,Ageom)
            else if(bitype.eq.-10) then
                !mismatch the beam at given location.
                !here, drange(3:8) stores the mismatch factor.
                call getparam_BeamLineElem(Blnelem(i),drange)
                call scale_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4),drange(5),&
                               drange(6),drange(7),drange(8))
            else if(bitype.eq.-11) then !MSU stripper model
                deallocate(tmpptcs)
                allocate(tmpptcs(8,Nplocal))
                gam = -Bpts%refptcl(6)
                gambetz = sqrt(gam**2-1.0)
                nfile = bmpstp
                beta0 = sqrt(1.0-(1.0/gam)**2)
                do ii = 1, Nplocal
                  tmpptcs(1,ii) = Bpts%Pts1(9,ii)
                  tmpptcs(2,ii) = Bpts%Pts1(1,ii)*Scxl*100
                  tmpptcs(3,ii) = Bpts%Pts1(3,ii)*Scxl*100
                  tmpptcs(4,ii) = Bpts%Pts1(5,ii)*90/asin(1.0)
                  tmpptcs(5,ii) = Bpts%Pts1(2,ii)/gambetz * 1000
                  tmpptcs(6,ii) = Bpts%Pts1(4,ii)/gambetz * 1000
                  tmpptcs(7,ii) = -Bpts%Pts1(6,ii)/(gam-1)*100
                  tmpptcs(8,ii) = Bpts%Pts1(7,ii)*Bpts%mass
                enddo

                call MStripF_BPM(beta0,tmpptcs,nfile,Nplocal,Nplocal,rwkinf,&
                                 rwkinq,Nchrg,Nptlist,Currlist,Qmcclist)

                Nptlist0 =  Nptlist
                Currlist0 =  Currlist
                Qmcclist0 =  Qmcclist/Bpts%mass

                do jj = 1, Nchrg
                    ii = 1
                    do while (ii .le. nchrgtot)
                        if (qmcclc(ii) .ne. Qmcclist0(jj)) then
                            nchrgtot = nchrgtot + 1
                            qmcclc(nchrgtot) = Qmcclist0(jj)
                            exit
                        endif
                        ii = ii + 1
                    enddo
                enddo

                gam = (1.0+rwkinq/Bmass)
                gambetz = sqrt(gam**2-1.0)
                avgw = 0.0
                do ii = 1, Nplocal
                    Bpts%Pts1(6,ii) = -(gam-1.0d0)*tmpptcs(7,ii)/100
                    Bpts%Pts1(2,ii) = tmpptcs(5,ii)/1000*gambetz
                    Bpts%Pts1(4,ii) = tmpptcs(6,ii)/1000*gambetz
                    Bpts%Pts1(7,ii) = tmpptcs(8,ii)/Bpts%mass
                    avgw = avgw + tmpptcs(7,ii)
                enddo
                avgw = avgw/Nplocal
                call chgupdate_BeamBunch(Bpts,Nchrg,Nptlist0,Qmcclist0)

                !update reference particle information
                Bcharge = Bmass*Qmcclist0(1)
                Bpts%Charge = Bcharge
                Bpts%refptcl(6) = -(1.0+rwkinq/Bmass)

                if (any(Flagdiag.eq.(/1, 2, 5, 6/))) then
                    call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                         Flagdiag, Outputflag, bstats(ciStats))
                endif

            else if(bitype.eq.-12) then !ANL stripper model
                nfile = bmpstp
                read(nfile,*)Nchrg,E0,El0
                read(nfile,*)(Nptlist(ichg),ichg=1,Nchrg)
                read(nfile,*)(Currlist(ichg),ichg=1,Nchrg)
                read(nfile,*)(Qmcclist(ichg),ichg=1,Nchrg)

                Bpts%refptcl(6) = -(1.0+(E0-El0)*1.0d6/Bmass)
                gam = -Bpts%refptcl(6)
                gambetz = sqrt(gam**2-1.0)

                nptottmp = 0
                do ichg = 1, Nchrg
                    nptottmp = nptottmp + Nptlist(ichg)
                enddo
                nptmp2 = 0
                do ichg = 1, Nchrg
                    if(ichg.ne.Nchrg) then
                        nplistlc = Nplocal*Nptlist(ichg)*1.d0/nptottmp
                        nptmp2 = nptmp2 + nplistlc
                    else
                        nplistlc = Nplocal - nptmp2
                    endif
                    do iev = 1, nplistlc
                        ! ..... Generate the energy loss and angle projections
                        !call gen_strip(2,0.d0,El,Thx,Thy)
                        ! ..... Calculate the ion exit energy
                        En  = E0 -El0 -El
                        Bpts%Pts1(2,iev) = Thx/1000*gambetz
                        Bpts%Pts1(4,iev) = Thy/1000*gambetz
                        Bpts%Pts1(6,iev) = gam-(1.0+En*1.0d6/Bmass)
                        Bpts%Pts1(7,iev) = Qmcclist(ichg)
                        Bpts%Pts1(8,iev) = Currlist(ichg)/Scfreq/Nptlist(i)*Qmcclist(i)/abs(Qmcclist(i))
                    enddo
                enddo

                Nptlist0 =  Nptlist
                Currlist0 =  Currlist
                Qmcclist0 =  Qmcclist

                do jj = 1, Nchrg
                    ii = 1
                    do while (ii .le. nchrgtot)
                        if (qmcclc(ii) .ne. Qmcclist0(jj)) then
                            nchrgtot = nchrgtot + 1
                            qmcclc(nchrgtot) = Qmcclist0(jj)
                            exit
                        endif
                        ii = ii + 1
                    enddo
                enddo

                !call chgupdate_BeamBunch(Bpts,Nchrg,Nptlist0,Qmcclist0)
                Bcharge = Bmass*Qmcclist0(1)
                Bpts%Charge = Bcharge
                nchrgtot = nchrgtot + Nchrg

            else if(bitype.eq.-13) then !collimator slit
                call getparam_BeamLineElem(Blnelem(i),drange)
                xradmin = drange(3)
                xradmax = drange(4)
                yradmin = drange(5)
                yradmax = drange(6)
                call lostcountXY_BeamBunch(Bpts,Nplocal,Np,xradmin,&
                      xradmax,yradmin,yradmax)
                call chgupdate_BeamBunch(Bpts,Nchrg,Nptlist0,Qmcclist0)

                if (any(Flagdiag.eq.(/1, 2/))) then
                    call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                         Flagdiag, Outputflag, bstats(ciStats))
                endif
            else if(bitype.eq.-14) then !Free Beam Parameter Selector
                call getparam_BeamLineElem(Blnelem(i),drange)
                cid = nint(drange(2))
                    ! {1:x, 2:xp, 3:y, 4:yp, 5:z, 6:zp, 7:q/m, 8:q, 9:id}
                    ! if cid = 0, skip selection
                tmp1 = drange(3)
                    ! lower limit  [m, rad, m, rad, deg, MeV, q/m, q, 1]
                tmp2 = drange(4)
                    ! higher limit [m, rad, m, rad, deg, MeV, q/m, q, 1]
                flg = nint(drange(6))
                    ! reference energy switch for cid=6.
                    ! 0: use relative energy (default)
                    ! 1: use absolute energy
                if (xradmin.ne.xradmax .and. cid.ge.1 .and. cid.le.9) then
                    call lostcountFP_BeamBunch(Bpts,Nplocal,Np,tmp1,tmp2,cid,flg)
                    call chgupdate_BeamBunch(Bpts,Nchrg,Nptlist0,Qmcclist0)
                endif

                if (any(Flagdiag.eq.(/1, 2/))) then
                    call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                         Flagdiag, Outputflag, bstats(ciStats))
                endif
            else if(bitype.eq.-21)then
                !shift the beam centroid in the 6D phase space.
                !This element can be used to model steering magnets etc.
                !here, drange(3:8) stores the amount of shift.
                !drange(3) : shift in x (m)
                !drange(4) : shift in Px (rad)
                !drange(5) : shift in y (m)
                !drange(6) : shift in Py (rad)
                !drange(7) : shift in z (deg)
                !drange(8) : shift in Pz (MeV)
                call getparam_BeamLineElem(Blnelem(i),drange)
                call kick_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4),drange(5),&
                              drange(6),drange(7),drange(8),-Bpts%refptcl(6),&
                              Bpts%Mass,Bpts%Charge)

                if (any(Flagdiag.eq.(/1, 2/))) then
                    call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                         Flagdiag, Outputflag, bstats(ciStats))
                endif

            else if(bitype.eq.-22)then
                !rotate the particle coordiante by drange(3) in x-y plane
                ! by drange(4) degrees in px-py plane. This function is
                ! most used in order to tanfer the beam through a vertic bend.
                call getparam_BeamLineElem(Blnelem(i),drange)
                call xypxpyrot_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4))
            else if(bitype.eq.-23)then
                !output flag
                if (any(Flagdiag.eq.(/1, 2, 3, 4/))) then
                    call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                         Flagdiag, Outputflag, bstats(ciStats))
                endif

            else if(bitype.eq.-25)then
                !based on kick_BPM, shift same amount of value for all particles
                !to simulate misalignment of element.
                !here, drange(3:8) stores the amount of shift.
                !drange(3) : shift in x (m)
                !drange(4) : shift in Px (rad)
                !drange(5) : shift in y (m)
                !drange(6) : shift in Py (rad)
                !drange(7) : shift in z (deg)
                !drange(8) : shift in Pz (MeV)
                call getparam_BeamLineElem(Blnelem(i),drange)
                call qzmis_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4),drange(5),&
                               drange(6),drange(7),drange(8),-Bpts%refptcl(6),&
                               Bpts%Mass,Bpts%Charge)
                !if(Flagdiag.eq.1) then
                !    call diagnostic1_Output(z,Bpts,Nchrg,Nptlist0)
                !elseif(Flagdiag.eq.2) then
                !    call diagnostic2_Output(Bpts,z,Nchrg,Nptlist0,Outputflag)
                !endif
                if (any(Flagdiag.eq.(/1, 2/))) then
                    call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                         Flagdiag, Outputflag, bstats(ciStats))
                endif
            else if(bitype.eq.-26) then
                call shift_BPM(Bpts%Pts1,bitype,Nplocal,Np)
            else if(bitype.eq.-27)then
                !frequency jump
                !drange(3): new frequency (Hz)
                !drange(4): phase shift (deg)
                call getparam_BeamLineElem(Blnelem(i),drange)
                tmp1 = drange(3)
                tmp2 = drange(4)

                do ii = 1, Nplocal
                  ! frequency change
                  Bpts%Pts1(1,ii) = Bpts%Pts1(1,ii)*tmp1/Bfreq
                  Bpts%Pts1(3,ii) = Bpts%Pts1(3,ii)*tmp1/Bfreq
                  Bpts%Pts1(5,ii) = Bpts%Pts1(5,ii)*tmp1/Bfreq

                  !phase shift
                  Bpts%Pts1(5,ii) = Bpts%Pts1(5,ii) + tmp2/Rad2deg
                  ! scaling
                  Bpts%Pts1(5,ii) = mod(Bpts%Pts1(5,ii)+Pi,2.0*Pi)-Pi
                enddo

                call construct_PhysConst(tmp1)

            else if(bitype.eq.-28) then
                !space charge switch
                call getparam_BeamLineElem(Blnelem(i),drange)
                Bpts%Current = drange(3)
            else if((bitype.eq.-40).and.(Flagmap.ne.1)) then
                !hard-edge solenoid wrapper
                !drange(3) : field strength
                !drange(4) : 0:in, 1:out
                call getparam_BeamLineElem(Blnelem(i),drange)
                gamma0 = -Bpts%refptcl(6)
                call hedge_Sol2(Bpts%Pts1,Nplocal,gamma0,drange)
            else if(bitype.eq.-99) then
                exit
            endif

            if (bitype.eq.103) then
                call getparam_BeamLineElem(Blnelem(i),12,fsync)
                if (nint(fsync).eq.1) then
                    call getparam_BeamLineElem(Blnelem(i),1,flagv)
                    call getparam_BeamLineElem(Blnelem(i),26,strp0)
                    call getparam_BeamLineElem(Blnelem(i),27,strgam0)
                    call getparam_BeamLineElem(Blnelem(i),28,strqm)
                    if (flagv.lt.0.0 .or. strp0.ne.Bpts%refptcl(5) .or. &
                        strgam0.ne.-Bpts%refptcl(6) .or. strqm.ne.dble(Bpts%Charge/Bpts%Mass)) then
                        call phasescan_BeamBunch(Blnelem(i),Flagmap,Bpts%Charge,Bpts%Mass,Bpts%refptcl)
                        Blnparams(25,i) = tmptheta
                    endif
                else
                    call setparam_BeamLineElem(Blnelem(i), 25, dble(0.0))
                endif
                call getparam_BeamLineElem(Blnelem(i), 25, tmptheta)
                Blnparams(25,i) = tmptheta
            elseif (bitype.eq.110) then
                call getparam_BeamLineElem(Blnelem(i),15,fsync)
                if (nint(fsync).eq.1) then
                    call getparam_BeamLineElem(Blnelem(i),1,flagv)
                    call getparam_BeamLineElem(Blnelem(i),29,strp0)
                    call getparam_BeamLineElem(Blnelem(i),30,strgam0)
                    call getparam_BeamLineElem(Blnelem(i),31,strqm)
                    if (flagv.lt.0.0 .or. strp0.ne.Bpts%refptcl(5) .or. &
                        strgam0.ne.-Bpts%refptcl(6) .or. strqm.ne.dble(Bpts%Charge/Bpts%Mass)) then
                        call phasescan_BeamBunch(Blnelem(i),Flagmap,Bpts%Charge,Bpts%Mass,Bpts%refptcl)
                    endif
                else
                    call setparam_BeamLineElem(Blnelem(i), 28, dble(0.0))
                endif
                call getparam_BeamLineElem(Blnelem(i), 28, tmptheta)
                Blnparams(28,i) = tmptheta
            endif

            ! loop through 'bnseg' numerical segments in each beam element
            ! using 2 step symplectic integeration (ie. leap frog).
            call setparam_BeamLineElem(Blnelem(i),1,zedge)
            if(bitype.ne.4) then  !//no bend
                do j = 1, bnseg
                    if (bitype.eq.106 .and. Flagmap.eq.2) then
                        if (Flagmap.eq.1) then
                            wk = 2*Pi/blength
                            call getparam_BeamLineElem(Blnelem(i),1,zedge)
                            call getparam_BeamLineElem(Blnelem(i),7,ma)
                            zz = z - zedge
                            xrad = piperad*((ma+1)/2-(ma-1.0)/2*cos(wk*zz))
                            yrad = piperad*((ma+1)/2+(ma-1.0)/2*cos(wk*zz))
                            piperad = xrad
                            piperad2 = yrad
                        else if (Flagmap.eq.2) then
                            ncl = Rfqsegs(j)
                            if (ncl.ne.tmprfqncl) then
                                fncl = float(ncl)
                                call setparam_BeamLineElem(Blnelem(i), 13, fncl)
                                tau2 = Rfqzlen(ncl)/Rfqnseg(ncl)
                                nsubstep = Rfqpstp(ncl)
                                stdcnt = int(Rfqpara(14,ncl) + 0.1)
                                piperad0 = Rfqpara(3,ncl)
                                ma = Rfqpara(4,ncl)
                                celltype = Rfqltyp(ncl)
                                if (any(celltype.eq.(/0, 1/))) then
                                    if (stdcnt.ge.2) then
                                        wk = 2.0*Pi/(Rfqzlen(ncl)+Rfqzlen(ncl-1))
                                    else
                                        wk = Pi/Rfqzlen(ncl)
                                    endif
                                else if (celltype.eq.5) then
                                    wk = 0.0
                                else
                                    wk = 0.5*Pi/Rfqzlen(ncl)
                                endif
                                tmprfqncl = ncl
                            endif

                            if (any(celltype.eq.(/0,1/))) then
                                if (mod(stdcnt,2).eq.1) then
                                    zz = z - zedge + Rfqzlen(ncl)
                                else
                                    zz = z - zedge
                                endif
                            else
                                zz = z - zedge
                            endif

                            piperad  = piperad0*0.5*((ma+1.0)-(ma-1.0)*cos(wk*zz))
                            piperad2 = piperad0*0.5*((ma+1.0)+(ma-1.0)*cos(wk*zz))
                            nprof = Rfqppos(ncl)

                            if (nprof.ne.0) then
                                call phase_Store(Rfqprof(1,nprof),Bpts,Rfqprof(2,nprof))
                            endif
                        endif
                    else if (bitype.eq.114 .and. Flagmap.eq.2) then
                        zadv = zadv + 0.5*tau2
                        call getbendref_EMfld(zadv,bndrad,zlfrg,zlbnd,refang,refcen,refdz)
                        piperad = piperad0/cos(Pi*refang/180.0)
                        drefcen = refcen0 - refcen
                        refcen0 = refcen
                        Bpts%Pts1(1,:) = Bpts%Pts1(1,:) + drefcen/Scxl
                        call setparam_BeamLineElem(Blnelem(i), 15, refcen)
                    endif

                    if(Flagerr.ge.1 .and. bitype.ne.0) then
                        call geomerrL2_BeamBunch(Bpts,Blnelem(i),z-zedge,Flagerr,Flagmap)
                    end if

                    !Hard edge Solenoid option (entrance)
                    if (bitype.eq.13 .and. j.eq.1) then
                        call getparam_BeamLineElem(Blnelem(i),3,tmp2)
                        if (nint(tmp2).eq.1) then
                            call getparam_BeamLineElem(Blnelem(i),2,tmp1)
                            drange(3) = tmp1
                            drange(4) = 0
                            gamma0 = -Bpts%refptcl(6)
                            call hedge_Sol2(Bpts%Pts1,Nplocal,gamma0,drange)
                        endif
                    endif

                    ! use linear map or nonlinear Lorentz integrator to advance particles.
                    ! spatial drift.
                    !linear map integrator
                    if(Flagmap.eq.1) then
                        call map1_BeamBunch(Bpts,Blnelem(i),z,tau1)
                    else
                        if(bitype.ne.3) then
                            call map1_BeamBunch(Bpts,z,tau2)
                        else
                            z = z + 0.5*tau2
                            gamma0 = -Bpts%refptcl(6)
                            beta0 = sqrt(1.0-1.0/gamma0/gamma0)
                            Bpts%refptcl(5) = Bpts%refptcl(5)+0.5*tau2/(Scxl*beta0)
                            call getparam_BeamLineElem(Blnelem(i),2,gtf)
                            if(Bpts%Current.ge.1.0e-30)  then
                                ! half advance for space-charge kick
                                call linearmod_Sol(Bpts%Nptlocal,Bpts%Pts1,0.5*tau2,&
                                                   gamma0,0.5*blength*100.0,gtf*10000.0)
                            else
                                call linearmod_Sol(Bpts%Nptlocal,Bpts%Pts1,tau2,&
                                                   gamma0,blength*100.0,gtf*10000.0)
                            endif
                        endif
                    endif

                    ! goto the space charge calculation
                    if(Bpts%Current.ge.1.0e-30)  then !space-charge flag
                        call conv1st_BeamBunch(Bpts,tau2,Nplocal,Np,ptrange,&
                                               Flagbc,Perdlen,piperad,piperad2)
                        call chgupdate_BeamBunch(Bpts,Nchrg,Nptlist0,Qmcclist0)
                        !fix the global range for sub-cycle of space charge potential.
                        if(Flagsubstep.eq.1) then
                            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                                ptrange(1) = 0.0
                                ptrange(2) = piperad
                                ptrange(3) = 0.0
                                ptrange(4) = 4*asin(1.0)
                            else
                                ptrange(1) = -piperad
                                ptrange(2) =  piperad
                                ptrange(3) = -piperad2
                                ptrange(4) =  piperad2
                            endif
                            ptrange(5) = -Perdlen/2
                            ptrange(6) =  Perdlen/2
                        endif

                        ! get new boundary from the range of beam particles.
                        if(Flagbc.ne.4) then
                            call update_CompDom(Ageom,ptrange,Grid2d,Flagbc)
                        endif

                        call getlcmnum_CompDom(Ageom,lcgrid)
                        Nxlocal = lcgrid(1)
                        if(npy.gt.1) then
                            Nylocal = lcgrid(2) + 2
                        else
                            Nylocal = lcgrid(2)
                        endif
                        if(npx.gt.1) then
                            Nzlocal = lcgrid(3) + 2
                        else
                            Nzlocal = lcgrid(3)
                        endif

                        call getlcrange_CompDom(Ageom,lcrange)
                        if(totnp.gt.1) then
                            ! pass particles to local space domain on new processor.
                            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,Grid2d,Pdim,&
                                                      Nplcmax,lcrange)
                            else
                                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,Grid2d,Pdim,&
                                                    Nplcmax,lcrange)
                            endif
                        endif

                        ! assign new 'Nplocal' local particles on each processor.
                        call setnpt_BeamBunch(Bpts,Nplocal)
                        if((mod(j-1,nsubstep).eq.0).or.(Flagsubstep.ne.1)) then
                            call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,Grid2d,npx,npy)
                            deallocate(chgdens)
                            allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
                            ! deposit particles onto grid to obtain charge density.
                            call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,Ageom,&
                                                  Grid2d,chgdens,Flagbc,Perdlen)


                            ! start load balance. (at the location of new space charge calculation)
                            if((mod(ibal,nbal).eq.0).and.(totnp.gt.1) ) then
                                call MPI_BARRIER(comm2d,ierr)
                                !if(myid.eq.0) then
                                !    print*," load balance! "
                                !endif

                                call getlctabnm_CompDom(Ageom,temptab)
                                lctabnmx(0:npx-1) = temptab(1,0:npx-1,0)
                                lctabnmy(0:npy-1) = temptab(2,0,0:npy-1)
                                call getmsize_CompDom(Ageom,msize)
                                hy = msize(2)
                                hz = msize(3)
                                call getrange_CompDom(Ageom,srange)
                                ymin = srange(3)
                                zmin = srange(5)
                                if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                                    call balance_CompDom(chgdens,lctabnmx,&
                                    lctabrgx,npx,npy,commrow,commcol, &
                                    lcgrid(1),lcgrid(2),lcgrid(3),Ny,Nz,hz,zmin)
                                    call setlctab_CompDom(Ageom,lctabnmx,lctabrgx,&
                                                          npx,npy,myidx,myidy)
                                else
                                    call balance_CompDom(chgdens,lctabnmx,&
                                    lctabnmy,lctabrgx,lctabrgy,npx,npy,commrow,commcol, &
                                    lcgrid(1),lcgrid(2),lcgrid(3),Ny,Nz,hy,hz,ymin,zmin)
                                    call setlctab_CompDom(Ageom,lctabnmx,lctabnmy,lctabrgx,&
                                                          lctabrgy,npx,npy,myidx,myidy)
                                endif

                                call getlcmnum_CompDom(Ageom,lcgrid)
                                Nxlocal = lcgrid(1)
                                if(npy.gt.1) then
                                    Nylocal = lcgrid(2) + 2
                                else
                                    Nylocal = lcgrid(2)
                                endif

                                if(npx.gt.1) then
                                    Nzlocal = lcgrid(3) + 2
                                else
                                    Nzlocal = lcgrid(3)
                                endif

                                call getlcrange_CompDom(Ageom,lcrange)
                                deallocate(chgdens)
                                allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
                                call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,Grid2d,npx,npy)

                                ! pass particles to local space domain on new processor.
                                if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                                    call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,Grid2d,Pdim,&
                                                          Nplcmax,lcrange)
                                else
                                    call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,Grid2d,Pdim,&
                                                         Nplcmax,lcrange)
                                endif
                                ! assign new 'Nplocal' local particles on each processor.
                                call setnpt_BeamBunch(Bpts,Nplocal)

                                ! deposit particles onto grid to obtain charge density.
                                call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,&
                                                      Ageom,Grid2d,chgdens,Flagbc,Perdlen)
                            endif

                            ibal = ibal + 1

                            if(npx.gt.1) then
                                nzlcr = Nzlocal-2
                            else
                                nzlcr = Nzlocal
                            endif
                            if(npy.gt.1) then
                                nylcr = Nylocal-2
                            else
                                nylcr = Nylocal
                            endif

                            ! solve 3D Poisson's equation
                            if(Flagbc.eq.1) then
                                ! solve Poisson's equation using 3D isolated boundary condition.
                                call update3O_FieldQuant(Potential,chgdens,Ageom,&
                                                         Grid2d,Nxlocal,Nylocal,Nzlocal,&
                                                         npx,npy,nylcr,nzlcr)
                            else if(Flagbc.eq.2) then
                                ! solve Poisson's equation using 2D isolated 1D periodic
                                ! boundary condition.
                                if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
                                call update2O1P_FieldQuant(Potential,chgdens,Ageom,&
                                                           Grid2d,Nxlocal,Nylocal,Nzlocal,&
                                                           npx,npy,nylcr,nzlcr)
                            else if(Flagbc.eq.3) then
                                if(myidy.eq.(npy-1)) nylcr = nylcr - 1
                                call update3_FieldQuant(Potential,chgdens,Ageom,&
                                                        Grid2d,Nxlocal,Nylocal,Nzlocal,&
                                                        npx,npy,nylcr,nzlcr,&
                                                        besscoef,bessnorm,gml,modth,nmod)
                            else if(Flagbc.eq.4) then
                                if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
                                if(myidy.eq.(npy-1)) nylcr = nylcr - 1
                                call update4_FieldQuant(Potential,chgdens,Ageom,&
                                                        Grid2d,Nxlocal,Nylocal,Nzlocal,&
                                                        npx,npy,nylcr,nzlcr,Perdlen)
                            else if(Flagbc.eq.5) then
                                ! solve Poisson's equation using 2D rectangular pipe,
                                ! 1D open boundary condition
                                if(myidy.eq.(npy-1)) nylcr = nylcr-1
                                call update5_FieldQuant(Potential,chgdens,Ageom,&
                                                        Grid2d,Nxlocal,Nylocal,Nzlocal,&
                                                        npx,npy,nylcr,nzlcr)
                            else if(Flagbc.eq.6) then
                                if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
                                if(myidy.eq.(npy-1)) nylcr = nylcr - 1
                                call update6_FieldQuant(Potential,chgdens,Ageom,&
                                                        Grid2d,Nxlocal,Nylocal,Nzlocal,&
                                                        npx,npy,nylcr,nzlcr)
                            else
                                print*,"no such boundary condition type!!!"
                                stop
                            endif
                        endif

                        call cvbkforth1st_BeamBunch(Bpts)
                        if(totnp.gt.1) then
                            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,Grid2d,Pdim,&
                                                      Nplcmax,lcrange)
                            else
                                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,Grid2d,Pdim,&
                                                     Nplcmax,lcrange)
                            endif
                        endif
                        call setnpt_BeamBunch(Bpts,Nplocal)
                    endif

                    ! use linear map or nonlinear Lorentz integrator to advance particles.
                    if(Flagmap.eq.1) then
                        ! kick particles in velocity space.
                        call map2_BeamBunch(Bpts,tau2,Nxlocal,Nylocal,Nzlocal,&
                                            Potential%FieldQ,Ageom,Grid2d,Flagbc,Perdlen)
                        call map1_BeamBunch(Bpts,Blnelem(i),z,tau1)
                    else
                        if (bitype.eq.106) then
                            call lostcountRFQ_BeamBunch(Bpts,Nplocal,Np,&
                                                    piperad,piperad2,piperad0)
                        else if (bitype.eq.114) then
                            xradmin = -piperad
                            xradmax =  piperad
                            yradmin = -piperad2
                            yradmax =  piperad2
                            call lostcountXY_BeamBunch(Bpts,Nplocal,Np,&
                                        xradmin, xradmax, yradmin, yradmax)
                        else
                            call lostcount_BeamBunch(Bpts,Nplocal,Np,piperad,piperad2)
                        endif

                        call map2_BeamBunch(Bpts,Blnelem(i),z,tau2,Nxlocal,Nylocal,&
                                            Nzlocal,Potential%FieldQ,Ageom,Grid2d,Flagbc,Flagerr)
                        if(bitype.ne.3) then
                            call map1_BeamBunch(Bpts,z,tau2)
                        else
                            z = z + 0.5*tau2
                            gamma0 = -Bpts%refptcl(6)
                            beta0 = sqrt(1.0-1.0/gamma0/gamma0)
                            Bpts%refptcl(5) = Bpts%refptcl(5) +0.5*tau2/(Scxl*beta0)
                            call getparam_BeamLineElem(Blnelem(i),2,gtf)
                            if(Bpts%Current.ge.1.0e-30)  then
                                ! half advance for space-charge kick
                                call linearmod_Sol(Bpts%Nptlocal,Bpts%Pts1,0.5*tau2,&
                                                   gamma0,0.5*blength*100.0,gtf*10000.0)
                            endif
                        endif
                    endif

                    !Hard edge Solenoid option (exit)
                    if (bitype.eq.13 .and. j.eq.bnseg) then
                        call getparam_BeamLineElem(Blnelem(i),3,tmp2)
                        if (nint(tmp2).eq.1) then
                            call getparam_BeamLineElem(Blnelem(i),2,tmp1)
                            drange(3) = tmp1
                            drange(4) = 1
                            gamma0 = -Bpts%refptcl(6)
                            call hedge_Sol2(Bpts%Pts1,Nplocal,gamma0,drange)
                        endif
                    endif

                    if(Flagerr.ge.1 .and. bitype.ne.0) then
                        call geomerrT2_BeamBunch(Bpts,Blnelem(i),z-zedge,Flagerr,Flagmap)
                    end if

                    if (bitype.eq.114 .and. Flagmap.eq.2) then
                        zadv = zadv + 0.5*tau2
                        call getbendref_EMfld(zadv,bndrad,zlfrg,zlbnd,refang,refcen,refdz)
                        drefcen = refcen0 - refcen
                        refcen0 = refcen
                        Bpts%Pts1(1,:) = Bpts%Pts1(1,:) + drefcen/Scxl
                        if (j.eq.bnseg) then
                            z = z + blength-zlength
                            print*, blength-zlength
                        else
                            z = z + refdz
                        endif
                        gam = -Bpts%refptcl(6)
                        call yrotation_BPM(Bpts%Pts1,Nplocal,refang,gam)
                    endif

                    if (any(Flagdiag.eq.(/1, 2/))) then
                        call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                             Flagdiag, Outputflag, bstats(ciStats))
                    endif

                    if (bitype.eq.114 .and. Flagmap.eq.2 .and. j.ne.bnseg) then
                        z = z - refdz
                        call yrotation_BPM(Bpts%Pts1,Nplocal,-refang,gam)
                    endif

                    nstep = nstep + 1
                end do

            else !//bend magnet using Time domain
                call getparam_BeamLineElem(Blnelem(i),dparam)
                if(dparam(4).gt.100.0) then !file id > 100 to use 2nd order matrix
                    if(Flagerr.ge.1) then
                        call geomerrL2_BeamBunch(Bpts,Blnelem(i),z-zedge,Flagerr,Flagmap)
                    end if
                    z = z + blength
                    dpi = 2*asin(1.0)
                    gambet00 = dparam(3)
                    gamma00 = sqrt(1+gambet00*gambet00)
                    beta00 = sqrt(1.0-1.0/gamma00/gamma00)
                    gamma0 = -Bpts%refptcl(6)
                    beta0 = sqrt(1.0-1.0/gamma0/gamma0)
                    Bpts%refptcl(5) = Bpts%refptcl(5)+blength/(Scxl*beta0)
                    gambet = beta0*gamma0
                    hgap = 2*dparam(5)
                    ang0 = dparam(2)*dpi/180.0
                    hd1=0.0
                    ! 20130607  this entry is used by "reference energy of beta-gamma"
                    ! therefore, hd1 is hardwire into 0, as it should be in most of the normal cases
                    angF = dparam(6)*dpi/180.0 !e1
                    angB = dparam(7)*dpi/180.0 !e2
                    hF = dparam(8) !1/r1
                    hB = dparam(9) !/1/r2
                    dstr1 = dparam(10) !fringe field K of entrance
                    dstr2 = dstr1 !fringe field K of exit. here, assume Kb = Kf
                    dk1 = dparam(11)

                    hd0 = ang0/blength !k0
                    tanphiF = tan(angF)
                    psi1 = hd0*hgap*dstr1*(1.0+sin(angF)*sin(angF))/cos(angF)
                    tanphiFb = tan(angF-psi1)
                    tanphiB = tan(angB)
                    psi2 = hd0*hgap*dstr2*(1.0+sin(angB)*sin(angB))/cos(angB)
                    tanphiBb = tan(angB-psi2)
                    qm0 = Bpts%Charge/Bpts%Mass
                    cosang0 = cos(ang0)
                    sinang0 = sin(abs(ang0))

                    r0  = 1.0/hd0

                    do ipt = 1, Nplocal
                        ptarry(1) = Bpts%Pts1(1,ipt)*Scxl
                        gamn = gamma0 - Bpts%Pts1(6,ipt)
                        gambetz = sqrt(gamn**2-1.0-Bpts%Pts1(2,ipt)**2-Bpts%Pts1(4,ipt)**2)
                        ptarry(2) = Bpts%Pts1(2,ipt)/gambetz
                        ptarry(3) = Bpts%Pts1(3,ipt)*Scxl
                        ptarry(4) = Bpts%Pts1(4,ipt)/gambetz
                        ptarry(5) = -Bpts%Pts1(5,ipt)*beta00*Scxl
                        ptarry(6) = (-Bpts%Pts1(6,ipt)+(gamma0-gamma00))/beta00/gambet00 - &
                                    (Bpts%Pts1(7,ipt)-qm0)/qm0
                        qmrel = (Bpts%Pts1(7,ipt)-qm0)/qm0

                        !transformed input array for the "sector": ptarray -> ptarray2
                        if(dparam(4) .eq. 400.0) then
                            call Fpol_Dipole(hd0,hF,tanphiF,tanphiFb,hd1,&
                                             psi1,ptarry,ptarry2,angF)
                        else
                            ptarry2 = ptarry
                        endif

                        ! ptarray2 -> ptarray
                        if(dparam(4) .eq. 1500.0 ) then
                            call ESector_Dipole(blength,beta00,hd0,dk1,ptarry2,ptarry,0)
                        else if(dparam(4) .eq. 1501.0 ) then
                            call ESector_Dipole(blength,beta00,hd0,dk1,ptarry2,ptarry,1)
                        else
                            call Sector_Dipole(blength,beta00,hd0,dk1,ptarry2,ptarry,qmrel)
                        endif

                        !transformed array for skipping "exit edge" : ptarray -> ptarray2
                        if(dparam(4) .eq. 600.0) then
                            call Bpol_Dipole(hd0,hB,tanphiB,tanphiBb,hd1,&
                                             psi2,ptarry,ptarry2,angB)
                        else
                            ptarry2 = ptarry
                        endif

                        Bpts%Pts1(1,ipt) = ptarry2(1)/Scxl
                        Bpts%Pts1(3,ipt) = ptarry2(3)/Scxl
                        Bpts%Pts1(5,ipt) = -ptarry2(5)/(Scxl*beta00)
                        gambetz = sqrt((gamn**2-1)/(1+ptarry2(2)**2+ptarry2(4)**2))
                        Bpts%Pts1(2,ipt) = ptarry2(2)*gambetz
                        Bpts%Pts1(4,ipt) = ptarry2(4)*gambetz
                    enddo

                    if(Flagerr.ge.1) then
                        call geomerrT2_BeamBunch(Bpts,Blnelem(i),z-zedge,Flagerr,Flagmap)
                    end if

                    if (any(Flagdiag.eq.(/1, 2/))) then
                        call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                             Flagdiag, Outputflag, bstats(ciStats))
                    endif

                else
                    ! go to T frame
                    call convZT_BeamBunch(Bpts)
                    ! loop through bnseg steps
                    gam = sqrt(1.0+Bpts%refptcl(6)**2)
                    ! normalized t (omega t)
                    tau2 = 2*pi*Scfreq*blength/(Clight*Bpts%refptcl(6)/gam)/bnseg
                    tg = Bpts%refptcl(5)
                    tv = 0.0
                    !call diagnosticT_Output(tv,Bpts)
                    Flagbctmp = 1
                    do j = 1, bnseg
                        call drifthalfT_BeamBunch(Bpts,tv,tau2)
                        ptref = Bpts%refptcl
                        ! go to the local "ptref" coordinates
                        call rotto_BeamBunch(Bpts,ptref,ptrange)
                        if(Bpts%Current.ge.1.0e-30) then
                            call update_CompDom(Ageom,ptrange,Grid2d,Flagbctmp)
                            call getlcmnum_CompDom(Ageom,lcgrid)
                            Nxlocal = lcgrid(1)
                            if(npy.gt.1) then
                                Nylocal = lcgrid(2) + 2
                            else
                                Nylocal = lcgrid(2)
                            endif
                            if(npx.gt.1) then
                                Nzlocal = lcgrid(3) + 2
                            else
                                Nzlocal = lcgrid(3)
                            endif
                            call getlcrange_CompDom(Ageom,lcrange)
                            call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,Grid2d,Pdim,&
                                                 Nplcmax,lcrange)
                            ! assign new 'Nplocal' local particles on each processor.
                            call setnpt_BeamBunch(Bpts,Nplocal)
                            call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,Grid2d,npx,npy)
                            deallocate(chgdens)
                            allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
                            ! deposit particles onto grid to obtain charge density.
                            call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,Ageom,&
                                                  Grid2d,chgdens,Flagbc,Perdlen)
                            ibal = ibal + 1

                            if(npx.gt.1) then
                                nzlcr = Nzlocal-2
                            else
                                nzlcr = Nzlocal
                            endif

                            if(npy.gt.1) then
                                nylcr = Nylocal-2
                            else
                                nylcr = Nylocal
                            endif

                            call update3O_FieldQuant(Potential,chgdens,Ageom,&
                                                     Grid2d,Nxlocal,Nylocal,Nzlocal,&
                                                     npx,npy,nylcr,nzlcr)
                        endif
                        !//kick particles in the rotated local coordinates.
                        call kickT_BeamBunch(Bpts,Blnelem(i),tv,tau2,Nxlocal,Nylocal,&
                                             Nzlocal,Potential%FieldQ,Ageom,Grid2d,Flagbc,Flagerr)
                        call rotback_BeamBunch(Bpts,ptref)
                        call drifthalfT_BeamBunch(Bpts,tv,tau2)
                        tv = tv + tau2
                        !call diagnosticT_Output(tv,Bpts)
                    enddo

                    tg = tg+tv
                    call convTZ_BeamBunch(Bpts,tg) !//go back to Z frame
                    z = z + blength
                    if (any(Flagdiag.eq.(/1, 2/))) then
                        call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                             Flagdiag, Outputflag, bstats(ciStats))
                    endif
                endif
            endif  !//end bend magnets

            if (any(Flagdiag.eq.(/5, 6/))) then
                call diag_std_Output(z, Bpts, Nchrg, Qmcclist0, Nptlist0,&
                                     Flagdiag, Outputflag, bstats(ciStats))
            endif

        enddo

        call phase_init_Store(Bpts)
        if (ninedata_id.eq.0) ninedata_id = -1

        call gatherStats_Output(ndiags, nchrgtot, qmcclc, bstats, Flagdiag)

        ! final output.
        call MPI_BARRIER(comm2d,ierr)
        t_integ = t_integ + elapsedtime_Timer(t0)
        !call showtime_Timer()

        deallocate(lctabnmx,lctabnmy)
        deallocate(lctabrgx,lctabrgy)
        deallocate(temptab)
        deallocate(tmpptcs)
        deallocate(chgdens)
        if(Flagbc.eq.3) then
            deallocate(besscoef)
            deallocate(bessnorm)
            deallocate(gml)
            deallocate(pydisp)
            deallocate(modth)
        endif

    end subroutine run_Accel

    subroutine destruct_Accel(time)
        implicit none
        include 'mpif.h'
        double precision :: time

        !call destruct_Data()
        if (associated(Bpts%Pts1)) call destruct_BeamBunch(Bpts)
        if (associated(Potential%FieldQ)) call destruct_FieldQuant(Potential)
        if (associated(Ageom%LcTabrg).and.associated(Ageom%LcTabnm)) then
            call destruct_CompDom(Ageom)
        endif

        if (allocated(Blnparams)) deallocate(Blnparams)
        if (allocated(Ielems)) deallocate(Ielems)

    end subroutine destruct_Accel


    subroutine mpi_Setup(time)
        implicit none
        include 'mpif.h'
        integer :: ierr
        double precision :: time
        call MPI_INIT(ierr)
        time = MPI_WTIME()
        call init_Data()
    end subroutine mpi_Setup


    subroutine mpi_Teardown(time)
        implicit none
        include 'mpif.h'
        integer :: myrank, ierr
        double precision :: time, endtime, mtime
        call destruct_Data()

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
        endtime = MPI_WTIME()
        time = endtime - time
        call MPI_REDUCE(time,mtime,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        !if (myrank.eq.0) print*,"time: ",mtime
        call MPI_Finalize(ierr)

    end subroutine mpi_Teardown

end module AdvAccelclass





