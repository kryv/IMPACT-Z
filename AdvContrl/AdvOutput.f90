!!20120224	change the output format of tag "1" to that of tag "2" for 24,25 & 26
!!QZ20120210  changed into "real, intent(in) :: outputflag",
!!	if "10<outputflag<100", then use "outputflag" emittance, otherwise use "99.9%"
!!QZ06022011 no changes as to 03/22/2010, but checked currently using 99.9%
!!--QZ only calculate 99% values instead of 99.9% for ~2k particles 12/15/2009
!!--format output 32 with "132 format" line 10/21/2009
!!--QZ only calculate 99.9% values instead of 99.99% 10/21/2009
!!--QZ only calculate 99.99% values instead of 99.5% 02/05/08
!----------------------------------------------------------------
! (c) Copyright, 2003 by the Regents of the University of California.
! Outputclass: Output class in I/O module of CONTROL layer.
! Version: 2.0
! Author: Ji Qiang, LBNL, 7/24/03
! Description: This class defines functions to print out the charged
!              particle beam information in the accelerator.
! Comments: We have added 3 more attributes to the particle:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
      module AdvOutputclass
        use BeamBunchclass
        use Pgrid2dclass
        use Fldmgerclass
        use PhysConstclass
        use Messaging

        type BunchStorage
          double precision :: Current,Mass,Charge, Scxl
          integer :: Npt,Nptlocal
          double precision, allocatable, dimension(:,:) :: Pts1
          double precision, dimension(6) :: refptcl
        end type BunchStorage

        type BunchStats
            integer :: lpmin, lpmax, nchrg
            integer, allocatable, dimension(:) :: ptot, npchg

            double precision,dimension(5) :: ref
            ! ref(1) : z        (2) : phi       (3) : gamma
            !    (4) : energy   (5) : beta
            double precision, allocatable, dimension(:) :: qmcc
            double precision, allocatable, dimension(:,:) :: xyz
            ! cs = 1 : whole beam stats
            ! xyz( 1, cs) : xcen    ( 2, cs) : ycen    ( 3, cs) : zcen
            !    ( 4, cs) : xrms    ( 5, cs) : yrms    ( 6, cs) : zrms
            !    ( 7, cs) : xpcen   ( 8, cs) : ypcen   ( 9, cs) : zpcen
            !    (10, cs) : xprms   (11, cs) : yprms   (12, cs) : zprms
            !    (13, cs) : xtwsa   (14, cs) : ytwsa   (15, cs) : ztwsa
            !    (16, cs) : xtwsb   (17, cs) : ytwsb   (18, cs) : ztwsb
            !    (19, cs) : xepsn   (20, cs) : yepsn   (21, cs) : zepsn
            !    (22, cs) : xmax    (23, cs) : ymax    (24, cs) : phimax
            !    (25, cs) : xpmax   (26, cs) : ypmax   (27, cs) : demax
            !    (28, cs) : refrmax
            !    (29, cs) : x3rd    (30, cs) : y3rd    (31, cs) : z3rd
            !    (32, cs) : xp3rd   (33, cs) : yp3rd   (34, cs) : zp3rd
            !    (35, cs) : x4th    (36, cs) : y4th    (37, cs) : z4th
            !    (38, cs) : xp4th   (39, cs) : yp4th   (40, cs) : zp4th
            !    (41, cs) : xepsnp  (42, cs) : yepsnf  (43, cs) : zepsnp
            !    (44, cs) : xepsnf  (45, cs) : yepsnf  (46, cs) : zepsnf
            !    (47, 1 ) : rrms    (48, cs) : rnpr   (49, cs) : rmax
        end type BunchStats

        integer :: ciStats
        integer, private :: ciHist, ciSp

        type(BunchStorage), private :: inibunch
        type(BunchStorage), allocatable, private, dimension(:) :: hbunch

      contains

        subroutine Output_Hist(fname, outflag, qid)
            implicit none
            include 'mpif.h'
            integer, intent(in) :: outflag, qid
            character(256), intent(in) :: fname
            integer :: i, fnum, hlen, myrank, ierr
            character(256) :: f18, f24, f28, f29, f32, fout

            call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

            if (myrank.eq.0) then
                f18 = '(1x,e13.6,1x,e15.8,6(1x,e13.6))'
                f24 = '(10(1x,e13.6))'
                f28 = '(1x,e13.6,3I10)'
                f29 = '(7(1x,e13.6))'
                f32 = '(1x,e13.6,10(1x,I7))'
                hlen = size(hrefz)

                open(18, file=trim(fname)//'.18')
                open(27, file=trim(fname)//'.27')
                open(28, file=trim(fname)//'.28')
                open(32, file=trim(fname)//'.32')
                do i = 1, hlen
                    write(18,f18), hrefz(i), hrefphi(i), hrefgma(i), hrefeng(i),&
                                   hrefbeta(i), hrefr(qid,i)
                    write(27,f24), hrefz(i), hxmax(qid,i), hxpmax(qid,i), hymax(qid,i),&
                                   hypmax(qid,i), hphimax(qid,i), hdemax(qid,i)
                    write(28,f28), hrefz(i), hlpmin(i), hlpmax(i), hptot(qid,i)

                    if (hnchg(i) .eq. 1) then
                        write(32,f32), hrefz(i), hnpchg(1:hnchg(i),i)
                    else
                        write(32,f32), hrefz(i), hnpchg(2:hnchg(i)+1,i)
                    endif
                enddo

                close(18)
                close(27)
                close(28)
                close(32)

                if (outflag.eq.1 .or. outflag.eq.3) then

                    open(24, file=trim(fname)//'.24')
                    open(25, file=trim(fname)//'.25')
                    open(26, file=trim(fname)//'.26')
                    open(29, file=trim(fname)//'.29')
                    open(30, file=trim(fname)//'.30')
                    do i = 1, hlen
                        write(24,f24), hrefz(i), hxcen(qid,i), hxrms(qid,i), hxpcen(qid,i), hxprms(qid,i),&
                                       hxtwsa(qid,i), hxtwsb(qid,i), hxepsn(qid,i), hxepsn(qid,i), hxepsn(qid,i)
                        write(25,f24), hrefz(i), hycen(qid,i), hyrms(qid,i), hypcen(qid,i), hyprms(qid,i),&
                                       hytwsa(qid,i), hytwsb(qid,i), hyepsn(qid,i), hyepsn(qid,i), hyepsn(qid,i)
                        write(26,f24), hrefz(i), hzcen(qid,i), hzrms(qid,i), hzpcen(qid,i), hzprms(qid,i),&
                                       hztwsa(qid,i), hztwsb(qid,i), hzepsn(qid,i), hzepsn(qid,i), hzepsn(qid,i)
                        write(29,f24), hrefz(i), hx3rd(qid,i), hxp3rd(qid,i), hy3rd(qid,i), hyp3rd(qid,i),&
                                       hz3rd(qid,i), hzp3rd(qid,i)
                        write(30,f24), hrefz(i), hx4th(qid,i), hxp4th(qid,i), hy4th(qid,i), hyp4th(qid,i),&
                                       hz4th(qid,i), hzp4th(qid,i)
                    enddo
                    close(24)
                    close(25)
                    close(26)
                    close(29)
                    close(30)

                else if (outflag.eq.2 .or. outflag.eq.4) then

                    open(24, file=trim(fname)//'.24')
                    open(25, file=trim(fname)//'.25')
                    open(26, file=trim(fname)//'.26')
                    open(29, file=trim(fname)//'.29')
                    do i = 1, hlen
                        write(24,f24), hrefz(i), hxcen(qid,i), hxrms(qid,i), hxpcen(qid,i), hxprms(qid,i),&
                                       hxtwsa(qid,i), hxtwsb(qid,i), hxepsn(qid,i), hxepsnp(qid,i), hxepsnf(qid,i)
                        write(25,f24), hrefz(i), hycen(qid,i), hyrms(qid,i), hypcen(qid,i), hyprms(qid,i),&
                                       hytwsa(qid,i), hytwsb(qid,i), hyepsn(qid,i), hyepsnp(qid,i), hyepsnf(qid,i)
                        write(26,f24), hrefz(i), hzcen(qid,i), hzrms(qid,i), hzpcen(qid,i), hzprms(qid,i),&
                                       hztwsa(qid,i), hztwsb(qid,i), hzepsn(qid,i), hzepsnp(qid,i), hzepsnf(qid,i)
                        write(29,f29), hrefz(i), hrrms(qid,i), hrnpr(qid,i), hrmax(qid,i)

                    enddo
                    close(24)
                    close(25)
                    close(26)
                    close(29)

                endif
            endif
        end subroutine Output_Hist

        subroutine init_SpData(plen)
            implicit none
            integer :: iflag = 0
            integer, intent(in) :: plen

            ciSp = 1

            if (iflag.eq.1) then
                deallocate(hbunch, hbunchid, hbunchsmp)
            endif
            iflag = 1
            allocate(hbunch(plen))
            allocate(hbunchid(plen))
            allocate(hbunchsmp(plen))
        end subroutine init_SpData

        subroutine get_Dist(rid)
            implicit none
            integer, intent(in) :: rid

            if (rid.eq.0) then
                call gather_Particles(inibunch)
            else
                call gather_Particles(hbunch(rid))
            endif

            ninedata_id = rid

        end subroutine get_Dist

        subroutine gather_Particles(this)
            implicit none
            include 'mpif.h'
            type(BunchStorage), intent(in) :: this
            integer :: i, nprocs, nptot, lpos, myrank, ierr
            integer, allocatable, dimension(:) :: nptlst, shead
            double precision :: gam

            call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
            call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)


            allocate(nptlst(nprocs))
            allocate(shead(nprocs))

            call MPI_GATHER(this%Nptlocal, 1, MPI_INTEGER, nptlst, 1,&
                            MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            lpos = 0
            do i=1,nprocs
                shead(i) = lpos
                lpos = lpos + 9*nptlst(i)
            enddo

            nptot = sum(nptlst)

            if (myrank.eq.0) then
                if (allocated(ninedata)) deallocate(ninedata)
                allocate(ninedata(9,nptot))

                gam = -this%refptcl(6)
                ninedata_bg = sqrt(gam**2-1.0)
                ninedata_mass = this%Mass
                ninedata_scxl = this%Scxl
            endif

            call MPI_GATHERV(this%Pts1(1,1), 9*this%Nptlocal, &
                             MPI_DOUBLE_PRECISION, ninedata(1,1), 9*nptlst, shead,&
                             MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        end subroutine gather_Particles

        subroutine Output_Dist(fname)
            implicit none
            include 'mpif.h'
            character(256), intent(in) :: fname
            character(256) :: fileid, frmt
            integer :: myrank, ierr, i, j, rid, dlen, pnum, smpstp

            call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

            dlen = size(hbunch)
            do i=1,dlen
                rid = hbunchid(i)
                smpstp = hbunchsmp(i)
                if (smpstp.eq.0) smpstp = 1
                write(fileid,*), rid
                call get_Dist(i)
                if (myrank.eq.0) then
                    frmt = '(9(1x,e14.7))'
                    open(rid, file=trim(fname)//'.'//trim(adjustl(fileid)))
                    pnum = size(ninedata)/9
                    do j=1,pnum,smpstp
                        write(rid,frmt), ninedata(:,j)
                    enddo
                endif

                close(rid)
            enddo

        end subroutine Output_Dist

        subroutine diag_std_Output(z, this, nchrg, qmcclist, nptlist, fdiag, fop, stats)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), intent(in) :: this
        integer, intent(in) :: nchrg, fdiag
        double precision, dimension(nchrg), intent(in) :: qmcclist
        integer, dimension(nchrg), intent(in) :: nptlist
        double precision, intent(in) :: fop
        type (BunchStats), intent(out) :: stats

        integer :: i, j, k, my_rank, ierr, innp, nptot, ci, nlen,&
                   blen, npctmin, npctmax
        integer, dimension(nchrg) :: lnpc, bnpc

        integer, allocatable, dimension(:) :: cidx, ptot

        double precision :: qmc, xl, xt, weight, tden1, twa, twb, twg,&
                            r2lc, r2avg, rrlc, rravg, gam, gambet2, gambet, &
                            bet, energy, mevu2rad, deg2m
        double precision, dimension(6) :: txyz, txyz0, txyz1
        double precision, dimension(6*4+4) :: tlxyz
        double precision, dimension(nchrg) :: den1, lwgh, bwgh
        double precision, dimension(6*4+4, nchrg) :: lxyz, bxyz
        double precision, dimension(3, nchrg) :: lcpl
        double precision, dimension(7, nchrg) :: lmax, bmax

        double precision, allocatable, dimension(:) :: &
                                xpx, ypy, zpz, epsx, epsy, epsz, &
                                xtmp, ytmp, ptmp, pptmp, epstmp
        double precision, allocatable, dimension(:,:) :: &
                                dxyz, x0n, y0n, z0n, px0n, py0n, pz0n,&
                                npcepsx, npcepsy, npcepsz, npcr

        ciStats = ciStats + 1

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg

        innp = this%Nptlocal
        nptot = this%Npt

        lnpc = 0
        lwgh = 0.0
        lxyz = 0.0
        lcpl = 0.0
        lmax = 0.0

        allocate(cidx(innp))
        do i = 1, innp
            do j = 1, nchrg
                if (this%Pts1(7,i) .eq. qmcclist(j)) then
                    cidx(i) = j
                    exit
                else
                    cidx(i) = 1
                endif
            enddo
        enddo

        do i = 1, innp

            ci = cidx(i)

            lnpc(ci) = lnpc(ci) + 1
            if (this%Pts1(8, i) .le. 0.0) then
                weight = 1.0
            else
                weight = this%Pts1(8, i)/this%Pts1(7, i)/this%Mass
            endif

            ! first order stats
            txyz0 = this%Pts1(1:6, i)
            tlxyz(1:6) = txyz0
            tlxyz(7:9) = txyz0(1:6:2)*txyz0(2:6:2)
            do j = 1, 6
                if (lmax(j, ci) .lt. abs(txyz0(j))) then
                    lmax(j, ci) = abs(txyz0(j))
                endif
            enddo
            ! second order stats
            txyz = txyz0*txyz0
            tlxyz(10:15) = txyz
            if (lmax(7,ci) .lt. (txyz(1)+txyz(3))) then
                lmax(7,ci) = (txyz(1)+txyz(3))
            endif
            ! third order stats
            txyz = txyz*txyz0
            txyz1 = txyz
            txyz1(5) = abs(txyz1(5))
            tlxyz(16:21) = txyz1
            ! forth order stats
            txyz = txyz*txyz0
            tlxyz(22:27) = txyz
            tlxyz(28) = 1.0
            lxyz(:,ci) = lxyz(:,ci) + tlxyz*weight
        enddo

        if (mod(fdiag, 2) .eq. 0) then
            call MPI_ALLREDUCE(lxyz, bxyz, 28*nchrg, MPI_DOUBLE_PRECISION,&
                            MPI_SUM, MPI_COMM_WORLD, ierr)
        else
            call MPI_REDUCE(lxyz, bxyz, 28*nchrg, MPI_DOUBLE_PRECISION,&
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        endif

        call MPI_ALLREDUCE(lnpc, bnpc, nchrg, MPI_INTEGER, MPI_SUM,&
                           MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(lmax, bmax, 7*nchrg, MPI_DOUBLE_PRECISION,&
                        MPI_MAX, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        den1 = 1.0/bxyz(28,:)
        tden1 = 1.0/sum(bxyz(28,:))

        if (nchrg .eq. 1) then
            ci = 1
            nlen = 1
            allocate(dxyz(27,1))
        else
            ci = 2
            nlen = nchrg + 1

            allocate(dxyz(27, nchrg+1))
            do i = 1, 27
                dxyz(i,1) = sum(bxyz(i,:))*tden1
            end do
        endif

        do i = 1, 27
            dxyz(i,ci:) = bxyz(i,:)*den1
        end do

        allocate(x0n(4,nlen))
        allocate(px0n(4,nlen))
        allocate(y0n(4,nlen))
        allocate(py0n(4,nlen))
        allocate(z0n(4,nlen))
        allocate(pz0n(4,nlen))

        call calc_moments(nlen, dxyz(1,:), dxyz(10,:), dxyz(16,:), dxyz(22,:), x0n)
        call calc_moments(nlen, dxyz(2,:), dxyz(11,:), dxyz(17,:), dxyz(23,:), px0n)
        call calc_moments(nlen, dxyz(3,:), dxyz(12,:), dxyz(18,:), dxyz(24,:), y0n)
        call calc_moments(nlen, dxyz(4,:), dxyz(13,:), dxyz(19,:), dxyz(25,:), py0n)
        call calc_moments(nlen, dxyz(5,:), dxyz(14,:), dxyz(20,:), dxyz(26,:), z0n)
        call calc_moments(nlen, dxyz(6,:), dxyz(15,:), dxyz(21,:), dxyz(27,:), pz0n)

        allocate(xpx(nlen))
        allocate(ypy(nlen))
        allocate(zpz(nlen))
        allocate(epsx(nlen))
        allocate(epsy(nlen))
        allocate(epsz(nlen))

        xpx = dxyz(7,:) - x0n(1,:)*px0n(1,:)
        ypy = dxyz(8,:) - y0n(1,:)*py0n(1,:)
        zpz = dxyz(9,:) - z0n(1,:)*pz0n(1,:)

        do i = 1, nlen
            epsx(i) = sqrt(max(x0n(2,i)*px0n(2,i) - xpx(i)*xpx(i), 0.0))
            epsy(i) = sqrt(max(y0n(2,i)*py0n(2,i) - ypy(i)*ypy(i), 0.0))
            epsz(i) = sqrt(max(z0n(2,i)*pz0n(2,i) - zpz(i)*zpz(i), 0.0))
        enddo

        allocate(ptot(nlen))
        if (nlen .eq. 1) then
            ptot(1) = bnpc(1)
        else
            ptot(1) = sum(bnpc(:))
            ptot(2:) = bnpc(:)
        endif

        gam = -this%refptcl(6)
        gambet2 = gam**2 - 1.0
        gambet = sqrt(gambet2)
        bet = gambet/gam

        if (mod(fdiag, 2) .eq. 0) then
            allocate(xtmp(innp))
            allocate(ytmp(innp))
            allocate(ptmp(innp))
            allocate(pptmp(innp))
            allocate(epstmp(innp))
            if (flagmcbin .eq. 0) then
                blen = 1
            else
                blen = nlen
            endif
            allocate(npcepsx(2, blen))
            allocate(npcepsy(2, blen))
            allocate(npcepsz(2, blen))
            allocate(npcr(2, blen))

            xtmp = this%Pts1(1,:) - x0n(1,1)
            pptmp = (this%Pts1(2,:) - px0n(1,1))
            twa = -xpx(1)/epsx(1)
            twb = x0n(2,1)*gambet/epsx(1)
            twg = (1.0+twa*twa)/twb
            twa = twa/gambet
            twb = twb/gambet2
            epstmp = twg*xtmp*xtmp + 2.0*twa*xtmp*pptmp + twb*pptmp*pptmp
            call calc_npc(innp, blen, ptot(1:blen), cidx, epstmp, fop, npcepsx)

            ytmp = this%Pts1(3,:) - y0n(1,1)
            pptmp = (this%Pts1(4,:) - py0n(1,1))
            twa = -ypy(1)/epsy(1)
            twb = y0n(2,1)*gambet/epsy(1)
            twg = (1.0+twa*twa)/twb
            twa = twa/gambet
            twb = twb/gambet2
            epstmp = twg*ytmp*ytmp + 2.0*twa*ytmp*pptmp + twb*pptmp*pptmp
            call calc_npc(innp, blen, ptot(1:blen), cidx, epstmp, fop, npcepsy)

            ptmp = (this%Pts1(5,:) - z0n(1,1))
            pptmp = (this%Pts1(6,:) - pz0n(1,1))
            twa = -zpz(1)/epsz(1)
            twb = z0n(2,1)/epsz(1)*gambet**3
            twg = (1.0+twa*twa)/twb
            twa = twa*bet/(gam*gambet2)
            twb = twb/(gam*gam*gambet2*gambet2)
            twg = twg*bet*bet
            epstmp = twg*ptmp*ptmp + 2.0*twa*ptmp*pptmp + twb*pptmp*pptmp
            call calc_npc(innp, blen, ptot(1:blen), cidx, epstmp, fop, npcepsz)

            epstmp = xtmp*xtmp + ytmp*ytmp
            r2lc = sum(epstmp)
            call MPI_REDUCE(r2lc,r2avg,1,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,0,MPI_COMM_WORLD,ierr)
            epstmp = sqrt(epstmp)
            rrlc = sum(epstmp)
            call MPI_REDUCE(rrlc,rravg,1,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call calc_npc(innp, blen, ptot(1:blen), cidx, epstmp, fop, npcr)
        endif

        if (my_rank .eq. 0) then
            allocate(stats%ptot(nlen))
            stats%ptot = ptot

            stats%nchrg = nchrg
            stats%lpmin = npctmin
            stats%lpmax = npctmax

            allocate(stats%qmcc(nchrg))
            stats%qmcc = qmcclist

            allocate(stats%npchg(nlen))
            if (nlen.eq.1) then
                stats%npchg = nptlist
            else
                stats%npchg(1) = sum(nptlist)
                stats%npchg(2:nlen) = nptlist
            endif

            energy = (gam - 1.0)*qmc
            mevu2rad = 1.0/energy*gam/(gam + 1.0)
            deg2m = 1.0/360.0*gambet/gam*Clight/Scfreq

            stats%ref(1) = z
            stats%ref(2) = this%refptcl(5)
            stats%ref(3) = gam
            stats%ref(4) = energy
            stats%ref(5) = bet

            allocate(stats%xyz(49, nlen))
            stats%xyz = 0.0

            stats%xyz(1,:) = x0n(1,:)*xl
            stats%xyz(2,:) = y0n(1,:)*xl
            stats%xyz(3,:) = z0n(1,:)*xt
            stats%xyz(4,:) = sqrt(x0n(2,:))*xl
            stats%xyz(5,:) = sqrt(y0n(2,:))*xl
            stats%xyz(6,:) = sqrt(z0n(2,:))*xt
            stats%xyz(7,:) = px0n(1,:)/gambet
            stats%xyz(8,:) = py0n(1,:)/gambet
            stats%xyz(9,:) = pz0n(1,:)*qmc
            stats%xyz(10,:) = sqrt(px0n(2,:))/gambet
            stats%xyz(11,:) = sqrt(py0n(2,:))/gambet
            stats%xyz(12,:) = sqrt(pz0n(2,:))*qmc
            stats%xyz(13,:) = -xpx/epsx
            stats%xyz(14,:) = -ypy/epsy
            stats%xyz(15,:) = zpz/epsz
            stats%xyz(16,:) = x0n(2,:)/epsx*xl*gambet
            stats%xyz(17,:) = y0n(2,:)/epsy*xl*gambet
            stats%xyz(18,:) = z0n(2,:)/epsz*xt/qmc*deg2m/mevu2rad
            stats%xyz(19,:) = epsx*xl
            stats%xyz(20,:) = epsy*xl
            stats%xyz(21,:) = epsz*qmc*xt*mevu2rad*deg2m*gambet

            stats%xyz(22,1) = maxval(bmax(1,:))*xl
            stats%xyz(23,1) = maxval(bmax(3,:))*xl
            stats%xyz(24,1) = maxval(bmax(5,:))*xt
            stats%xyz(25,1) = maxval(bmax(2,:))/gambet
            stats%xyz(26,1) = maxval(bmax(4,:))/gambet
            stats%xyz(27,1) = maxval(bmax(6,:))*qmc
            stats%xyz(28,1) = sqrt(maxval(bmax(7,:)))*xl
            if (ci .eq. 2) then
                stats%xyz(22,ci:) = bmax(1,:)*xl
                stats%xyz(23,ci:) = bmax(3,:)*xl
                stats%xyz(24,ci:) = bmax(5,:)*xt
                stats%xyz(25,ci:) = bmax(2,:)/gambet
                stats%xyz(26,ci:) = bmax(4,:)/gambet
                stats%xyz(27,ci:) = bmax(6,:)*qmc
                stats%xyz(28,ci:) = sqrt(bmax(7,:))*xl
            endif

            stats%xyz(29,:) = x0n(3,:)*xl
            stats%xyz(30,:) = y0n(3,:)*xl
            stats%xyz(31,:) = z0n(3,:)*xt
            stats%xyz(32,:) = px0n(3,:)/gambet
            stats%xyz(33,:) = py0n(3,:)/gambet
            stats%xyz(34,:) = pz0n(3,:)*qmc

            stats%xyz(35,:) = x0n(4,:)*xl
            stats%xyz(36,:) = y0n(4,:)*xl
            stats%xyz(37,:) = z0n(4,:)*xt
            stats%xyz(38,:) = px0n(4,:)/gambet
            stats%xyz(39,:) = py0n(4,:)/gambet
            stats%xyz(40,:) = pz0n(4,:)*qmc

            if (mod(fdiag, 2) .eq. 0) then
                if (flagmcbin .eq. 0 ) then
                    stats%xyz(41,:) = npcepsx(1,1)*gambet*xl
                    stats%xyz(42,:) = npcepsy(1,1)*gambet*xl
                    stats%xyz(43,:) = npcepsz(1,1)*gambet**2*gam**2*qmc*xt*mevu2rad*deg2m

                    stats%xyz(44,:) = npcepsx(2,1)*gambet*xl
                    stats%xyz(45,:) = npcepsy(2,1)*gambet*xl
                    stats%xyz(46,:) = npcepsz(2,1)*gambet**2*gam**2*qmc*xt*mevu2rad*deg2m

                    rravg = rravg/nptot
                    stats%xyz(47,:) = sqrt(r2avg/nptot - rravg*rravg)*xl
                    stats%xyz(48,:) = npcr(1,1)*xl
                    stats%xyz(49,:) = npcr(2,1)*xl
                else
                    stats%xyz(41,:) = npcepsx(1,:)*gambet*xl
                    stats%xyz(42,:) = npcepsy(1,:)*gambet*xl
                    stats%xyz(43,:) = npcepsz(1,:)*gambet**2*gam**2*qmc*xt*mevu2rad*deg2m

                    stats%xyz(44,:) = npcepsx(2,:)*gambet*xl
                    stats%xyz(45,:) = npcepsy(2,:)*gambet*xl
                    stats%xyz(46,:) = npcepsz(2,:)*gambet**2*gam**2*qmc*xt*mevu2rad*deg2m

                    rravg = rravg/nptot
                    stats%xyz(47,:) = sqrt(r2avg/nptot - rravg*rravg)*xl
                    stats%xyz(48,:) = npcr(1,:)*xl
                    stats%xyz(49,:) = npcr(2,:)*xl
                endif
            endif
        endif

        end subroutine diag_std_Output

        subroutine calc_moments(n, x0, sqx, cux, ftx, x0n)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: x0, sqx, cux, ftx
        double precision, dimension(4,n), intent(out) :: x0n
        double precision, dimension(n) :: tmp1, tmp2, tmp3
            tmp1 = x0*x0
            tmp2 = sqx*x0
            tmp3 = x0*tmp1
            x0n(1,:) = x0
            x0n(2,:) = sqx - tmp1
            x0n(3,:) = (abs(cux-3.0*tmp2+2.0*tmp3))**(1.0/3.0)
            x0n(4,:) = sqrt(sqrt(abs(ftx-4.0*cux*x0+6.0*tmp2*x0-3.0*tmp3*x0)))
        end subroutine calc_moments

        subroutine calc_npc(innp, nlen, ptot, cidx, epst, npc, npceps)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp, nlen
        integer, dimension(nlen) :: ptot
        integer, dimension(innp) :: cidx
        double precision, dimension(innp), intent(in) :: epst
        double precision, intent(in) :: npc
        double precision, dimension(2, nlen), intent(out) :: npceps

        integer, parameter :: nbin = 10000
        integer, dimension(nbin, nlen) :: binlc, bin
        integer, dimension(nlen) :: f995
        integer, dimension(innp) :: nilst
        integer :: i, j, k, my_rank, ierr, itmp, ni, ai, bi ,ci
        double precision :: eps, ex995
        double precision, dimension(nlen) :: epsmxlc, epsmx, hxeps

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        if (nlen .eq. 1) then
            epsmxlc(1) = maxval(epst)
            epsmxlc(1) = max(epsmxlc(1), -1.0e10)
        else
            epsmxlc = 0.0
            do i = 1, innp
                ci = cidx(i)
                epsmxlc(ci+1) = max(epsmxlc(ci+1), epst(i))
            enddo
            epsmxlc(1) = maxval(epsmxlc)
        endif

        call MPI_ALLREDUCE(epsmxlc,epsmx,nlen,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        eps = 1.0e-8
        hxeps = (epsmx + eps)/float(nbin)

        binlc = 0

        if (nlen .eq. 1) then
            nilst = int(epst/hxeps(1)) + 1
            do i = 1, innp
                ni = nilst(i)
                binlc(ni,1) = binlc(ni,1) + 1
            enddo
        else
            do i = 1, innp
                ni = int(epst(i)/hxeps(1)) + 1
                binlc(ni,1) = binlc(ni,1) + 1
                ci = cidx(i)
                ni = int(epst(i)/hxeps(ci+1)) + 1
                binlc(ni,ci+1) = binlc(ni,ci+1) + 1
            enddo
        endif

        call MPI_REDUCE(binlc,bin,nbin*nlen,MPI_INTEGER,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (npc.ge.10 .and. npc.lt.100) then
          f995 = npc/100.*ptot !!
        else
          f995 = 0.999*ptot    !!99.9% by default
        endif

        if(my_rank.eq.0) then
            do k = 1, nlen
                if (bin(1,k) .gt. f995(k)) then
                    npceps(1,k) = 0.0
                else
                    bi = bin(1,k)
                    do i = 2, nbin
                        ai = bi
                        bi = bi + bin(i,k)
                        if (bi .gt. f995(k)) then
                            npceps(1,k) = ((f995(k) - ai))/(bi-ai)*hxeps(k) + hxeps(k)*(i-1)
                            exit
                        endif
                    enddo
                endif
            enddo
            npceps(2,:) = epsmx(:)
        endif

        end subroutine calc_npc

        subroutine gatherStats_Output(n, nctot, qmcclc, bsl, fdiag)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: n, nctot, fdiag
        double precision, dimension(100) :: qmcclc
        type (BunchStats), dimension(n), intent(in) :: bsl
        integer :: iflag = 0
        integer :: i, j, k, cj, ck, my_rank, ierr, nlen

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        if (my_rank .eq. 0) then

            if (nctot.eq.1) then
                nlen = 1
            else
                nlen = nctot + 1
            endif

            if (allocated(hchrgpos)) deallocate(hchrgpos)
            allocate(hchrgpos(nlen, 2))
            hchrgpos = 0

            if (allocated(qmlabel)) deallocate(qmlabel)
            allocate(qmlabel(nctot))
            qmlabel = qmcclc(1:nctot)

            if (iflag .eq. 1) then
                deallocate(&
                    hrefz, hrefphi, hrefgma, hrefeng, hrefbeta, hrefr, &
                    hxcen, hycen, hzcen, hxrms, hyrms, hzrms, &
                    hxpcen, hypcen, hzpcen, hxprms, hyprms, hzprms, &
                    hxtwsa, hytwsa, hztwsa, hxtwsb, hytwsb, hztwsb, &
                    hxepsn, hyepsn, hzepsn, hxepsnp, hyepsnp, hzepsnp, &
                    hxepsnf, hyepsnf, hzepsnf, hxmax, hymax, hphimax, &
                    hxpmax, hypmax, hdemax, hlpmin, hlpmax, hptot, &
                    hrrms, hrnpr, hrmax, &
                    hx3rd, hy3rd, hxp3rd, hyp3rd, hz3rd, hzp3rd, &
                    hx4th, hy4th, hxp4th, hyp4th, hz4th, hzp4th, &
                    hnchg, hnpchg)
            endif

            iflag = 1

            allocate(hrefz(n))
            allocate(hrefphi(n))
            allocate(hrefgma(n))
            allocate(hrefeng(n))
            allocate(hrefbeta(n))

            allocate(hrefr(nlen, n))
            allocate(hxcen(nlen, n))
            allocate(hycen(nlen, n))
            allocate(hzcen(nlen, n))
            allocate(hxrms(nlen, n))
            allocate(hyrms(nlen, n))
            allocate(hzrms(nlen, n))
            allocate(hxpcen(nlen, n))
            allocate(hypcen(nlen, n))
            allocate(hzpcen(nlen, n))
            allocate(hxprms(nlen, n))
            allocate(hyprms(nlen, n))
            allocate(hzprms(nlen, n))
            allocate(hxtwsa(nlen, n))
            allocate(hytwsa(nlen, n))
            allocate(hztwsa(nlen, n))
            allocate(hxtwsb(nlen, n))
            allocate(hytwsb(nlen, n))
            allocate(hztwsb(nlen, n))
            allocate(hxepsn(nlen, n))
            allocate(hyepsn(nlen, n))
            allocate(hzepsn(nlen, n))
            allocate(hx3rd(nlen, n))
            allocate(hy3rd(nlen, n))
            allocate(hz3rd(nlen, n))
            allocate(hxp3rd(nlen, n))
            allocate(hyp3rd(nlen, n))
            allocate(hzp3rd(nlen, n))
            allocate(hx4th(nlen, n))
            allocate(hy4th(nlen, n))
            allocate(hz4th(nlen, n))
            allocate(hxp4th(nlen, n))
            allocate(hyp4th(nlen, n))
            allocate(hzp4th(nlen, n))
            allocate(hxmax(nlen, n))
            allocate(hymax(nlen, n))
            allocate(hphimax(nlen, n))
            allocate(hxpmax(nlen, n))
            allocate(hypmax(nlen, n))
            allocate(hdemax(nlen, n))
            allocate(hxepsnp(nlen, n))
            allocate(hyepsnp(nlen, n))
            allocate(hzepsnp(nlen, n))
            allocate(hxepsnf(nlen, n))
            allocate(hyepsnf(nlen, n))
            allocate(hzepsnf(nlen, n))
            allocate(hrrms(nlen, n))
            allocate(hrnpr(nlen, n))
            allocate(hrmax(nlen, n))

            allocate(hptot(nlen, n))
            allocate(hnpchg(nlen, n))

            allocate(hlpmin(n))
            allocate(hlpmax(n))
            allocate(hnchg(n))

            hrefz = 0d0
            hrefphi = 0d0
            hrefgma = 0d0
            hrefeng = 0d0
            hrefbeta = 0d0
            hrefr = 0d0
            hxcen = 0d0
            hycen = 0d0
            hzcen = 0d0
            hxrms = 0d0
            hyrms = 0d0
            hzrms = 0d0
            hxpcen = 0d0
            hypcen = 0d0
            hzpcen = 0d0
            hxprms = 0d0
            hyprms = 0d0
            hzprms = 0d0
            hxtwsa = 0d0
            hytwsa = 0d0
            hztwsa = 0d0
            hxtwsb = 0d0
            hytwsb = 0d0
            hztwsb = 0d0
            hxepsn = 0d0
            hyepsn = 0d0
            hzepsn = 0d0
            hx3rd = 0d0
            hy3rd = 0d0
            hxp3rd = 0d0
            hyp3rd = 0d0
            hz3rd = 0d0
            hzp3rd = 0d0
            hx4th = 0d0
            hy4th = 0d0
            hxp4th = 0d0
            hyp4th = 0d0
            hz4th = 0d0
            hzp4th = 0d0
            hxmax = 0d0
            hymax = 0d0
            hphimax = 0d0
            hxpmax = 0d0
            hypmax = 0d0
            hdemax = 0d0
            hptot = 0
            hnpchg = 0

            hxepsnp = 0d0
            hyepsnp = 0d0
            hzepsnp = 0d0
            hxepsnf = 0d0
            hyepsnf = 0d0
            hzepsnf = 0d0
            hrrms = 0d0
            hrnpr = 0d0
            hrmax = 0d0

            hlpmin = 0
            hlpmax = 0
            hnchg = 0

            do i = 1, n
                hrefz(i) = bsl(i)%ref(1)
                hrefphi(i) = bsl(i)%ref(2)
                hrefgma(i) = bsl(i)%ref(3)
                hrefeng(i) = bsl(i)%ref(4)
                hrefbeta(i) = bsl(i)%ref(5)
                hlpmin(i) = bsl(i)%lpmin
                hlpmax(i) = bsl(i)%lpmax
                hnchg(i) = bsl(i)%nchrg

                hptot(:, i) = bsl(i)%ptot
                hnpchg(:, i) = bsl(i)%npchg

                do j = 0, bsl(i)%nchrg
                    if (j .eq. 0) then
                        cj = 1
                    else
                        do k = 1, nctot
                            if (bsl(i)%qmcc(j) .eq. qmcclc(k)) then
                                cj = k + 1
                                exit
                            endif
                        enddo
                    endif

                    if (hchrgpos(cj, 1) .eq. 0) then
                        hchrgpos(cj, 1) = i
                    else
                        hchrgpos(cj, 2) = i
                    endif

                    hxcen(cj, i) = bsl(i)%xyz(1, j+1)
                    hycen(cj, i) = bsl(i)%xyz(2, j+1)
                    hzcen(cj, i) = bsl(i)%xyz(3, j+1)
                    hxrms(cj, i) = bsl(i)%xyz(4, j+1)
                    hyrms(cj, i) = bsl(i)%xyz(5, j+1)
                    hzrms(cj, i) = bsl(i)%xyz(6, j+1)
                    hxpcen(cj, i) = bsl(i)%xyz(7, j+1)
                    hypcen(cj, i) = bsl(i)%xyz(8, j+1)
                    hzpcen(cj, i) = bsl(i)%xyz(9, j+1)
                    hxprms(cj, i) = bsl(i)%xyz(10, j+1)
                    hyprms(cj, i) = bsl(i)%xyz(11, j+1)
                    hzprms(cj, i) = bsl(i)%xyz(12, j+1)
                    hxtwsa(cj, i) = bsl(i)%xyz(13, j+1)
                    hytwsa(cj, i) = bsl(i)%xyz(14, j+1)
                    hztwsa(cj, i) = bsl(i)%xyz(15, j+1)
                    hxtwsb(cj, i) = bsl(i)%xyz(16, j+1)
                    hytwsb(cj, i) = bsl(i)%xyz(17, j+1)
                    hztwsb(cj, i) = bsl(i)%xyz(18, j+1)
                    hxepsn(cj, i) = bsl(i)%xyz(19, j+1)
                    hyepsn(cj, i) = bsl(i)%xyz(20, j+1)
                    hzepsn(cj, i) = bsl(i)%xyz(21, j+1)
                    hxmax(cj, i) = bsl(i)%xyz(22, j+1)
                    hymax(cj, i) = bsl(i)%xyz(23, j+1)
                    hphimax(cj, i) = bsl(i)%xyz(24, j+1)
                    hxpmax(cj, i) = bsl(i)%xyz(25, j+1)
                    hypmax(cj, i) = bsl(i)%xyz(26, j+1)
                    hdemax(cj, i) = bsl(i)%xyz(27, j+1)
                    hrefr(cj, i) = bsl(i)%xyz(28, j+1)
                    hx3rd(cj, i) = bsl(i)%xyz(29, j+1)
                    hy3rd(cj, i) = bsl(i)%xyz(30, j+1)
                    hz3rd(cj, i) = bsl(i)%xyz(31, j+1)
                    hxp3rd(cj, i) = bsl(i)%xyz(32, j+1)
                    hyp3rd(cj, i) = bsl(i)%xyz(33, j+1)
                    hzp3rd(cj, i) = bsl(i)%xyz(34, j+1)
                    hx4th(cj, i) = bsl(i)%xyz(35, j+1)
                    hy4th(cj, i) = bsl(i)%xyz(36, j+1)
                    hz4th(cj, i) = bsl(i)%xyz(37, j+1)
                    hxp4th(cj, i) = bsl(i)%xyz(38, j+1)
                    hyp4th(cj, i) = bsl(i)%xyz(39, j+1)
                    hzp4th(cj, i) = bsl(i)%xyz(40, j+1)

                    if (mod(fdiag,2) .eq. 0) then
                        hxepsnp(cj, i) = bsl(i)%xyz(41, j+1)
                        hyepsnp(cj, i) = bsl(i)%xyz(42, j+1)
                        hzepsnp(cj, i) = bsl(i)%xyz(43, j+1)
                        hxepsnf(cj, i) = bsl(i)%xyz(44, j+1)
                        hyepsnf(cj, i) = bsl(i)%xyz(45, j+1)
                        hzepsnf(cj, i) = bsl(i)%xyz(46, j+1)
                        hrrms(cj, i) = bsl(i)%xyz(47, j+1)
                        hrnpr(cj, i) = bsl(i)%xyz(48, j+1)
                        hrmax(cj, i) = bsl(i)%xyz(49, j+1)
                    endif

                    if (bsl(i)%nchrg.eq.1) exit
                enddo

            enddo
        endif
        end subroutine gatherStats_Output

        ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance.
        subroutine diagnosticT_Output(z,this)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), intent(in) :: this
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet,gambet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(27) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax
        integer :: npctmin,npctmax

        call starttime_Timer(t0)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt
        !print*,"pts1: ",this%Pts1(1,1),this%Pts1(2,1),my_rank,innp

        den1 = 1.0/float(nptot)
        den2 = den1*den1
        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        z0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0

        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          x0lc3 = x0lc3 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)
          x0lc4 = x0lc4 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)*&
                  this%Pts1(1,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*this%Pts1(2,i)
          px0lc = px0lc + this%Pts1(2,i)
          sqsum2local = sqsum2local + this%Pts1(2,i)*this%Pts1(2,i)
          px0lc3 = px0lc3 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)
          px0lc4 = px0lc4 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)*&
                   this%Pts1(2,i)
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          y0lc3 = y0lc3 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)
          y0lc4 = y0lc4 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)*&
                  this%Pts1(3,i)
          ypylocal = ypylocal + this%Pts1(3,i)*this%Pts1(4,i)
          py0lc = py0lc + this%Pts1(4,i)
          sqsum4local = sqsum4local + this%Pts1(4,i)*this%Pts1(4,i)
          py0lc3 = py0lc3 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)
          py0lc4 = py0lc4 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)*&
                   this%Pts1(4,i)
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)
          z0lc3 = z0lc3 + abs(this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i))
          z0lc4 = z0lc4 + this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i)*&
                          this%Pts1(5,i)

          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc3 = pz0lc3 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc4 = pz0lc4 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)*&
                            this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4

        call MPI_REDUCE(tmplc,tmpgl,27,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        gambet = sqrt(gam**2-1.0)
        bet = sqrt(gam**2-1.0)/gam
        energy = (gam-1.)*qmc

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
          sqsum5 = sqz - z0*z0
          sqpz = tmpgl(12)*den1
          sqsum6 = sqpz - pz0*pz0
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
          zpz = tmpgl(15)*den1 - z0*pz0
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
                 3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(sqsum1)
          pxrms = sqrt(sqsum2)
          yrms = sqrt(sqsum3)
          pyrms = sqrt(sqsum4)
          zrms = sqrt(sqsum5)
          pzrms = sqrt(sqsum6)
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          write(48,100)z,this%refptcl(1)*xl,this%refptcl(2),this%refptcl(3)*xl,&
                      this%refptcl(4),this%refptcl(5)*xl,this%refptcl(6)
          write(34,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl
          write(35,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl
          write(36,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz/epz,epz*xl

          write(37,100)z,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                       glmax(4),glmax(5)*xl,glmax(6)
          write(38,101)z,npctmin,npctmax,nptot
          write(39,100)z,x03*xl,px03,y03*xl,py03,z03*xl,&
                       pz03
          write(40,100)z,x04*xl,px04,y04*xl,py04,z04*xl,&
                       pz04

          call flush(48)
          call flush(34)
          call flush(35)
          call flush(36)
          call flush(37)
          call flush(38)
          call flush(39)
          call flush(40)
        endif

99      format(6(1x,e13.6))
100      format(7(1x,e13.6))
101     format(1x,e13.6,3I10)

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnosticT_Output


        subroutine outpoint_Output(nfile,this,z,inb,jstp,nprocrow,nproccol,&
                                   geom)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: z
        type(CompDom), intent(in) :: geom
        integer, intent(in) :: inb,jstp,nprocrow,nproccol
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,Localnum)

        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))
        call getlctabrg_CompDom(geom,Localrange)

        call getmsize_CompDom(geom,msize)
        call getrange_CompDom(geom,range)

        open(nfile,status="unknown",form="unformatted")

        write(nfile)z
        write(nfile)inb,jstp
        write(nfile)msize(1:3)
        write(nfile)range(1:6)
        write(nfile)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        write(nfile)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        deallocate(Localnum)
        deallocate(Localrange)

        write(nfile)this%Nptlocal
        write(nfile)this%refptcl(1:6)
        write(nfile)this%Pts1(1:9,1:this%Nptlocal)

        close(nfile)

        end subroutine outpoint_Output

        subroutine phase_Store(nfile,this,smpstp)
            implicit none
            integer, intent(in) :: nfile, smpstp
            type (BeamBunch), intent(in) :: this

            hbunch(ciSp)%Current = this%Current
            hbunch(ciSp)%Mass = this%Mass
            hbunch(ciSp)%Charge = this%Charge
            hbunch(ciSp)%Npt = this%Npt
            hbunch(ciSp)%Nptlocal = this%Nptlocal
            hbunch(ciSp)%Pts1 = this%Pts1
            hbunch(ciSp)%refptcl = this%refptcl
            hbunch(ciSp)%Scxl = Scxl

            hbunchid(ciSp) = nfile
            hbunchsmp(ciSp) = smpstp

            ciSp = ciSp + 1
        end subroutine phase_Store

        subroutine phase_init_Store(this)
            implicit none
            type (BeamBunch), intent(in) :: this

            inibunch%Current = this%Current
            inibunch%Mass = this%Mass
            inibunch%Charge = this%Charge
            inibunch%Npt = this%Npt
            inibunch%Nptlocal = this%Nptlocal
            inibunch%Pts1 = this%Pts1
            inibunch%refptcl = this%refptcl
            inibunch%Scxl = Scxl

        end subroutine phase_init_Store

        subroutine phase_Output(nfile,this,samplePeriod)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: samplePeriod
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf

        if (samplePeriod .eq. 0) then
           samplePeriod = 1
        endif

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        call MPI_ALLREDUCE(this%Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)

        allocate(nptlist(0:np-1))
        nptlist = 0
        allocate(recvbuf(9,mnpt))
        sixnpt = 9*this%Nptlocal

        call MPI_GATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        nptlist = 9*nptlist

        if(my_rank.eq.0) then
          open(nfile,status='unknown')
          do i = 1, this%Nptlocal,samplePeriod
            write(nfile,100)this%Pts1(1,i),this%Pts1(2,i),this%Pts1(3,i),&
                            this%Pts1(4,i),this%Pts1(5,i),this%Pts1(6,i),&
                            this%Pts1(7,i),this%Pts1(8,i),this%Pts1(9,i)
          enddo
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                          i,1,MPI_COMM_WORLD,status,ierr)

            do j = 1, nptlist(i)/9,samplePeriod
              write(nfile,100)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),&
                              recvbuf(4,j),recvbuf(5,j),recvbuf(6,j),&
                              recvbuf(7,j),recvbuf(8,j),recvbuf(9,j)
            enddo
          enddo
          close(nfile)
        else
          call MPI_SEND(this%Pts1(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
                        MPI_COMM_WORLD,ierr)
        endif

100     format(9(1x,e14.7))

        deallocate(nptlist)
        deallocate(recvbuf)

        end subroutine phase_Output

        !//output the 3D particle number density
        subroutine dens3d_Output(nstep,nfile,this,totnptcls,xmni,xmxi,ymni,&
          ymxi,zmni,zmxi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nstep
        type (BeamBunch), intent(in) :: this
        integer, intent(in) :: totnptcls
        double precision, intent(in) :: xmni,xmxi,ymni,ymxi,zmni,zmxi
        integer, parameter :: nxcell = 100
        integer, parameter :: nycell = 100
        integer, parameter :: nzcell = 100
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,k,ii,jj,kk,icc
        integer :: ix,iy,iz,ic,nsend,req
        integer, dimension(nxcell,nycell,nzcell) :: count
        integer, dimension(4,nxcell*nycell*nzcell) :: sendbuf,recvbuf
        integer, allocatable,dimension(:) :: nptlist
        double precision :: hxx,hyy,hzz,invvol,tmp1,tmp2,tmp3,tmp4,eps
        double precision, dimension(6) :: range
        double precision, dimension(3) :: lcrange1,lcrange2,temp1
        double precision :: xl,xt
        double precision  :: xmnin,xmxin,ymnin,ymxin,zmnin,zmxin,xmin,xmax,&
        ymin,ymax,zmin,zmax

        xl = Scxl
        xt = 90.0/asin(1.0)

        xmnin = xmni/xl
        xmxin = xmxi/xl
        ymnin = ymni/xl
        ymxin = ymxi/xl
        zmnin = zmni/xt
        zmxin = zmxi/xt

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        eps = 1.0e-4
        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange2(i) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i-1) = temp1(i) + eps*temp1(i)
        enddo
        call MPI_ALLREDUCE(lcrange2,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + eps*temp1(i)
        enddo

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        zmin = min(zmnin,range(5))
        zmax = max(zmxin,range(6))

        hxx = (xmax-xmin)/nxcell
        hyy = (ymax-ymin)/nycell
        hzz = (zmax-zmin)/nzcell
        count = 0
        do i = 1, this%Nptlocal
          ix = int((this%Pts1(1,i)-range(1))/hxx) + 1
          iy = int((this%Pts1(3,i)-range(3))/hyy) + 1
          iz = int((this%Pts1(5,i)-range(5))/hzz) + 1
          count(ix,iy,iz) = count(ix,iy,iz) + 1
        enddo

        invvol = 1.0/(hxx*hyy*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do k = 1, nzcell
            do j =1, nycell
              do i =1, nxcell
                ic = (k-1)*nycell*nxcell+(j-1)*nxcell + i
                sendbuf(4,ic) = count(i,j,k)
              enddo
            enddo
          enddo
        else
          do k = 1, nzcell
            do j =1, nycell
              do i =1, nxcell
!                ic = (k-1)*nycell*nxcell+(j-1)*nxcell + i
!                sendbuf(4,ic) = 0
                if(count(i,j,k).gt.0) then
                  icc = icc + 1
                  sendbuf(1,icc) = i
                  sendbuf(2,icc) = j
                  sendbuf(3,icc) = k
                  sendbuf(4,icc) = count(i,j,k)
                endif
              enddo
            enddo
          enddo
        endif
        nsend = 4*icc
        allocate(nptlist(0:np-1))
        nptlist = 0

        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/4
              ii = recvbuf(1,j)
              jj = recvbuf(2,j)
              kk = recvbuf(3,j)
              icc = (kk-1)*nycell*nxcell+(jj-1)*nxcell+ii
              sendbuf(4,icc) = sendbuf(4,icc) + recvbuf(4,j)
            enddo
          enddo

          open(nfile,status='unknown')
          tmp4 = invvol/float(totnptcls)
          do k = 1, nzcell
            do j = 1, nycell
              do i = 1, nxcell
                tmp1 = (range(1) + i*hxx - 0.5*hxx)*xl
                tmp2 = (range(3) + j*hyy - 0.5*hyy)*xl
                tmp3 = (range(5) + k*hzz - 0.5*hzz)*xt
                ic = (k-1)*nycell*nxcell+(j-1)*nxcell+i
                write(nfile,100)tmp1,tmp2,tmp3,sendbuf(4,ic)*tmp4
              enddo
            enddo
          enddo
          close(nfile)
        else
          call MPI_ISEND(sendbuf(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

100     format(4(1x,e14.7))

        deallocate(nptlist)

        end subroutine dens3d_Output

        !//output 2D particle number density.
        subroutine dens2d_Output(nstep,nfile,this,totnptcls,xmnin,xmxin,&
        pxmnin,pxmxin,ymnin,ymxin,pymnin,pymxin,zmnin,zmxin,pzmnin,pzmxin)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nstep
        type (BeamBunch), intent(in) :: this
        integer, intent(in) :: totnptcls
        double precision, intent(in) :: xmnin,xmxin,pxmnin,pxmxin,ymnin,&
        ymxin,pymnin,pymxin,zmnin,zmxin,pzmnin,pzmxin
        integer, parameter :: nxcell = 128
        integer, parameter :: npxcell = 128
        integer, parameter :: nycell = 128
        integer, parameter :: npycell = 128
        integer, parameter :: nzcell = 128
        integer, parameter :: npzcell = 128
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,k,ii,req,jj,kk,icc
        integer :: ix,ipx,iy,ipy,iz,ipz,ic,nsend,nfile1,nfile2
        integer :: nfile3,nfile4,nfile5
        integer, dimension(nxcell,npxcell) :: count1
        integer, dimension(nxcell,npxcell) :: count2
        integer, dimension(nxcell,npxcell) :: count3
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(nxcell,nzcell) :: count5
        integer, dimension(nycell,nzcell) :: count6
        integer, dimension(3,nxcell*npxcell) :: sendbuf1,recvbuf1
        integer, dimension(3,nycell*npycell) :: sendbuf2,recvbuf2
        integer, dimension(3,nzcell*npzcell) :: sendbuf3,recvbuf3
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, dimension(3,nxcell*nzcell) :: sendbuf5,recvbuf5
        integer, dimension(3,nycell*nzcell) :: sendbuf6,recvbuf6
        integer, allocatable, dimension(:) :: nptlist
        double precision, dimension(6) :: range,prange
        double precision, dimension(6) :: lcrange1,lcrange2,temp1
        double precision :: hxx,hyy,hzz,hpxx,hpyy,hpzz,invvol
        double precision ::eps,tmp1,tmp2,tmp3,xmax,xmin,pxmax,pxmin,&
        ymax,ymin,pymax,pymin,zmin,zmax,pzmin,pzmax
        character*3 name
        integer :: ioerr
        integer gethostname

!        ioerr = gethostname(name,3)
        eps = 1.0e-4
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange1(i+3) = this%Pts1(2*i,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
            lcrange2(i+3) = this%Pts1(2*i,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange1(i+3) = 5000.0
            lcrange2(i) = -5000.0
            lcrange2(i+3) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
            ii = i + 3
            if(lcrange1(ii).gt.this%Pts1(2*i,j)) then
              lcrange1(ii) = this%Pts1(2*i,j)
            endif
            if(lcrange2(ii).lt.this%Pts1(2*i,j)) then
              lcrange2(ii) = this%Pts1(2*i,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,6,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)

        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo
        do i = 1, 3
          prange(2*i-1) = temp1(i+3) + temp1(i+3)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,6,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo
        do i = 1, 3
          prange(2*i) = temp1(i+3) + temp1(i+3)*eps
        enddo

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        zmin = min(zmnin,range(5))
        zmax = max(zmxin,range(6))
        pxmin = min(pxmnin,prange(1))
        pxmax = max(pxmxin,prange(2))
        pymin = min(pymnin,prange(3))
        pymax = max(pymxin,prange(4))
        pzmin = min(pzmnin,prange(5))
        pzmax = max(pzmxin,prange(6))

        hxx = (xmax-xmin)/nxcell
        hpxx = (pxmax-pxmin)/npxcell
        hyy = (ymax-ymin)/nycell
        hpyy = (pymax-pymin)/npycell
        hzz = (zmax-zmin)/nzcell
        hpzz = (pzmax-pzmin)/npzcell

        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        do i = 1, this%Nptlocal
          ix = int((this%Pts1(1,i)-range(1))/hxx) + 1
          ipx = int((this%Pts1(2,i)-prange(1))/hpxx) + 1
          count1(ix,ipx) = count1(ix,ipx) + 1
          iy = int((this%Pts1(3,i)-range(3))/hyy) + 1
          ipy = int((this%Pts1(4,i)-prange(3))/hpyy) + 1
          count2(iy,ipy) = count2(iy,ipy) + 1
          iz = int((this%Pts1(5,i)-range(5))/hzz) + 1
          ipz = int((this%Pts1(6,i)-prange(5))/hpzz) + 1
          count3(iz,ipz) = count3(iz,ipz) + 1
          count4(ix,iy) = count4(ix,iy) + 1
          count5(ix,iz) = count5(ix,iz) + 1
          count6(iy,iz) = count6(iy,iz) + 1
        enddo

        invvol = 1.0/(hxx*hpxx)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npxcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf1(3,ic) = count1(i,j)
            enddo
          enddo
        else
          do j = 1, npxcell
            do i = 1, nxcell
              if(count1(i,j).gt.0) then
                icc = icc + 1
                sendbuf1(1,icc) = i
                sendbuf1(2,icc) = j
                sendbuf1(3,icc) = count1(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        allocate(nptlist(0:np-1))
        nptlist = 0

        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf1(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf1(1,j)
              jj = recvbuf1(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf1(3,icc) = sendbuf1(3,icc) + recvbuf1(3,j)
            enddo
          enddo

!          open(unit=nfile,file="/n/"//name//&
!               "/scratch/jiqiang/XPx.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile,file='XPx.data',status='unknown',&
               position='append')
          tmp3 = invvol/float(totnptcls)
          write(nfile,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,prange(1)+0.5*hpxx,prange(1)+npxcell*hpxx-&
           0.5*hpxx
          do j = 1, npxcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = prange(1) + j*hpxx - 0.5*hpxx
!              write(nfile)sendbuf1(3,ic)
!              write(nfile)tmp1,tmp2,sendbuf1(3,ic)*tmp3
              write(nfile,100)tmp1,tmp2,sendbuf1(3,ic)*tmp3
            enddo
          enddo
          close(nfile)
        else
          call MPI_ISEND(sendbuf1(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hyy*hpyy)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npycell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              sendbuf2(3,ic) = count2(i,j)
            enddo
          enddo
        else
          do j = 1, npycell
            do i = 1, nycell
              if(count2(i,j).gt.0) then
                icc = icc + 1
                sendbuf2(1,icc) = i
                sendbuf2(2,icc) = j
                sendbuf2(3,icc) = count2(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf2(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf2(1,j)
              jj = recvbuf2(2,j)
              icc = (jj-1)*nycell+ii
              sendbuf2(3,icc) = sendbuf2(3,icc) + recvbuf2(3,j)
            enddo
          enddo

          nfile1 = nfile + 1
!          open(unit=nfile1,file="/n/"//name//&
!               "/scratch/jiqiang/YPy.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile1,file='YPy.data',status='unknown',&
               position='append')
          tmp3 = invvol/float(totnptcls)
          write(nfile1,*)nstep,range(3)+0.5*hyy,range(3)+nycell*hyy &
           -0.5*hyy,prange(3)+0.5*hpyy,prange(3)+npycell*hpyy-&
           0.5*hpyy
          do j = 1, npycell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              tmp1 = range(3) + i*hyy - 0.5*hyy
              tmp2 = prange(3) + j*hpyy - 0.5*hpyy
!              write(nfile1)sendbuf2(3,ic)
!              write(nfile1)tmp1,tmp2,sendbuf2(3,ic)*tmp3
              write(nfile1,100)tmp1,tmp2,sendbuf2(3,ic)*tmp3
            enddo
          enddo
          close(nfile1)
        else
          call MPI_ISEND(sendbuf2(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hzz*hpzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npzcell
            do i = 1, nzcell
              ic = (j-1)*nzcell + i
              sendbuf3(3,ic) = count3(i,j)
            enddo
          enddo
        else
          do j = 1, npzcell
            do i = 1, nzcell
              if(count3(i,j).gt.0) then
                icc = icc + 1
                sendbuf3(1,icc) = i
                sendbuf3(2,icc) = j
                sendbuf3(3,icc) = count3(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf3(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf3(1,j)
              jj = recvbuf3(2,j)
              icc = (jj-1)*nzcell + ii
              sendbuf3(3,icc) = sendbuf3(3,icc) + recvbuf3(3,j)
            enddo
          enddo

          nfile2 = nfile + 2
!          open(unit=nfile2,file="/n/"//name//&
!               "/scratch/jiqiang/ZPz.data",status="unknown",&
!               position="append",form="unformatted")
         open(unit=nfile2,file='ZPz.data',status='unknown',&
              position='append')
          tmp3 = invvol/float(totnptcls)
          write(nfile2,*)nstep,range(5)+0.5*hzz,range(5)+nzcell*hzz &
           -0.5*hzz,prange(5)+0.5*hpzz,prange(5)+npzcell*hpzz-&
           0.5*hpzz
          do j = 1, npzcell
            do i = 1, nzcell
              ic = (j-1)*nzcell + i
              tmp1 = range(5) + i*hzz - 0.5*hzz
              tmp2 = prange(5) + j*hpzz - 0.5*hpzz
!              write(nfile2)sendbuf3(3,ic)
!              write(nfile2)tmp1,tmp2,sendbuf3(3,ic)*tmp3
              write(nfile2,100)tmp1,tmp2,sendbuf3(3,ic)*tmp3
            enddo
          enddo
          close(nfile2)
        else
          call MPI_ISEND(sendbuf3(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hxx*hyy)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
!          open(unit=nfile3,file="/n/"//name//&
!                "/scratch/jiqiang/XY.data",status="unknown",&
!                position="append",form="unformatted")
          open(unit=nfile3,file='XY.data',status='unknown',&
                position='append')
          tmp3 = invvol/float(totnptcls)
          write(nfile3,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,range(3)+0.5*hyy,range(3)+nycell*hyy-&
           0.5*hyy
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = range(3) + j*hyy - 0.5*hyy
!              write(nfile3)sendbuf4(3,ic)
!              write(nfile3)tmp1,tmp2,sendbuf4(3,ic)*tmp3
              write(nfile3,100)tmp1,tmp2,sendbuf4(3,ic)*tmp3
            enddo
          enddo
          close(nfile3)
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hxx*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nzcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf5(3,ic) = count5(i,j)
            enddo
          enddo
        else
          do j = 1, nzcell
            do i = 1, nxcell
              if(count5(i,j).gt.0) then
                icc = icc + 1
                sendbuf5(1,icc) = i
                sendbuf5(2,icc) = j
                sendbuf5(3,icc) = count5(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf5(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf5(1,j)
              jj = recvbuf5(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf5(3,icc) = sendbuf5(3,icc) + recvbuf5(3,j)
            enddo
          enddo

          nfile4 = nfile + 4
!          open(unit=nfile4,file="/n/"//name//&
!               "/scratch/jiqiang/XZ.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile4,file='XZ.data',status='unknown',&
               position='append')
          tmp3 = invvol/float(totnptcls)
          write(nfile4,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,range(5)+0.5*hzz,range(5)+nzcell*hzz-&
           0.5*hzz
          do j = 1, nzcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = range(5) + j*hzz - 0.5*hzz
!              write(nfile4)sendbuf5(3,ic)
!              write(nfile4)tmp1,tmp2,sendbuf5(3,ic)*tmp3
              write(nfile4,100)tmp1,tmp2,sendbuf5(3,ic)*tmp3
            enddo
          enddo
          close(nfile4)
        else
          call MPI_ISEND(sendbuf5(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hyy*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nzcell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              sendbuf6(3,ic) = count6(i,j)
            enddo
          enddo
        else
          do j = 1, nzcell
            do i = 1, nycell
              if(count6(i,j).gt.0) then
                icc = icc + 1
                sendbuf6(1,icc) = i
                sendbuf6(2,icc) = j
                sendbuf6(3,icc) = count6(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf6(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf6(1,j)
              jj = recvbuf6(2,j)
              icc = (jj-1)*nycell + ii
              sendbuf6(3,icc) = sendbuf6(3,icc) + recvbuf6(3,j)
            enddo
          enddo

          nfile5 = nfile + 5
!          open(unit=nfile5,file="/n/"//name//&
!               "/scratch/jiqiang/YZ.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile5,file='YZ.data',status='unknown',&
               position='append')
          tmp3 = invvol/float(totnptcls)
          write(nfile5,*)nstep,range(3)+0.5*hyy,range(3)+nycell*hyy &
           -0.5*hyy,range(5)+0.5*hzz,range(5)+nzcell*hzz-&
           0.5*hzz
          do j = 1, nzcell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              tmp1 = range(3) + i*hyy - 0.5*hyy
              tmp2 = range(5) + j*hzz - 0.5*hzz
!              write(nfile5)sendbuf6(3,ic)
!              write(nfile5)tmp1,tmp2,sendbuf6(3,ic)*tmp3
              write(nfile5,100)tmp1,tmp2,sendbuf6(3,ic)*tmp3
            enddo
          enddo
          close(nfile5)
        else
          call MPI_ISEND(sendbuf6(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        deallocate(nptlist)

100     format(3(1x,e14.7))

        end subroutine dens2d_Output

        !//output the 1D accumulated particle density
        subroutine accdens1d_Output(nstep,nfile,this,nptot,rmxi,xmni,xmxi,&
                                 ymni,ymxi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nptot,nstep
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: rmxi,xmni,xmxi,ymni,ymxi
        integer, parameter :: nrcell = 128
        integer, parameter :: nxcell = 128
        integer, parameter :: nycell = 128
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, allocatable, dimension(:) :: nptlist
        integer, dimension(nxcell) :: densx,acdensx
        integer, dimension(nycell) :: densy,acdensy
        integer :: np,my_rank,ierr,i,ir,ic,icc,sumx,sumy,ix,iy,nsend,&
                   nfile3,ii,jj,req,j,testsum
        integer status(MPI_STATUS_SIZE)
        integer, dimension(nrcell) :: count,glcount
        double precision :: hrr,rmax,rad,xl,hxx,hyy,xmin,xmax,ymin,ymax
        double precision, dimension(6) :: range
        double precision, dimension(3) :: lcrange1,lcrange2,temp1
        double precision :: rmxin,xmnin,xmxin,ymnin,ymxin,eps

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        xl = Scxl

        rmxin = rmxi/xl
        xmnin = xmni/xl
        xmxin = xmxi/xl
        ymnin = ymni/xl
        ymxin = ymxi/xl

        eps = 0.001
        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange2(i) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        rmax = max(rmxin,sqrt(xmax*xmax+ymax*ymax))
        call MPI_BARRIER(MPI_COMM_WORLD,ierr);
        !print*,"xmin: ",xmin,xmax,ymin,ymax,rmax

        hxx = (xmax-xmin)/nxcell
        hyy = (ymax-ymin)/nycell
        hrr = rmax/nrcell
        count = 0
        count4 = 0
        do i = 1, this%Nptlocal
          rad = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+this%Pts1(3,i)*&
                     this%Pts1(3,i))
          ir = int(rad/hrr) + 1
          count(ir) = count(ir) + 1
          ix = int((this%Pts1(1,i)-xmin)/hxx) + 1
          iy = int((this%Pts1(3,i)-ymin)/hyy) + 1
          count4(ix,iy) = count4(ix,iy) + 1
        enddo

        call MPI_REDUCE(count,glcount,nrcell,MPI_INTEGER,MPI_SUM,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nrcell
            glcount(i) = glcount(i) + glcount(i-1)
          enddo
          !print*,"glcount: ",glcount(nrcell)

!          open(unit=nfile,file='RadDens.data',status='unknown',&
!               position='append',form='unformatted')
          open(unit=nfile,file='RadDens.data',status='unknown',&
               position='append')

          write(nfile,*)nstep
          do i = 1, nrcell
            !write(nfile,100)float(glcount(i))/glcount(nrcell)
            write(nfile,100)(i*hrr-0.5*hrr)*xl,float(glcount(i))/nptot
          enddo
          close(nfile)
        endif

        allocate(nptlist(0:np-1))
        nptlist = 0
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
          open(unit=nfile3,file='Yprof.data',status='unknown',&
                position='append')
          !testsum = 0
          do j = 1, nycell
            sumy = 0
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sumy = sumy + sendbuf4(3,ic)
            enddo
            !write(nfile3,100)float(sumy)/nptot
            !testsum = testsum + sumy
            densy(j) = sumy
          enddo
          acdensy = 0
          do j = nycell/2,1,-1
            acdensy(j) = densy(j) + densy(nycell-j+1) + acdensy(j+1)
          enddo
          do j = 1,nycell/2
            acdensy(j+nycell/2) = acdensy(nycell/2-j+1)
          enddo
          write(nfile3,*)nstep
          do j = 1, nycell
            write(nfile3,100)((j*hyy-0.5*hyy)+ymin)*xl,float(acdensy(j))/nptot
          enddo

          close(nfile3)
          !print*,"testsumy: ",testsum

          open(unit=nfile3,file='Xprof.data',status='unknown',&
                position='append')
          !testsum = 0
          do i = 1, nxcell
            sumx = 0
            do j = 1, nycell
              ic = (j-1)*nxcell + i
              sumx = sumx + sendbuf4(3,ic)
            enddo
            !write(nfile3,100)float(sumx)/nptot
            !testsum = testsum + sumx
            densx(i) = sumx
          enddo
          acdensx = 0
          do j = nxcell/2,1,-1
            acdensx(j) = densx(j) + densx(nxcell-j+1) + acdensx(j+1)
          enddo
          do j = 1,nxcell/2
            acdensx(j+nxcell/2) = acdensx(nxcell/2-j+1)
          enddo
          write(nfile3,*)nstep
          do j = 1, nxcell
            write(nfile3,100)((i*hxx-0.5*hxx)+xmin)*xl,float(acdensx(j))/nptot
          enddo
          close(nfile3)
          !print*,"testsumx: ",testsum
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

100     format(2(1x,e14.7))

        deallocate(nptlist)

        end subroutine accdens1d_Output

        !//output 1d particle density
        subroutine dens1d_Output(nstep,nfile,this,nptot,rmxi,xmni,xmxi,&
                                  ymni,ymxi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nptot,nstep
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: rmxi,xmni,xmxi,ymni,ymxi
        integer, parameter :: nrcell = 64
        integer, parameter :: nxcell = 64
        integer, parameter :: nycell = 64
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, allocatable, dimension(:) :: nptlist
        integer :: np,my_rank,ierr,i,ir,ic,icc,sumx,sumy,ix,iy,nsend,&
                   nfile3,ii,jj,req,j,testsum
        integer status(MPI_STATUS_SIZE)
        integer, dimension(nrcell) :: count,glcount
        double precision :: hrr,rmax,rad,xl,hxx,hyy,xmin,xmax,ymin,ymax,eps
        double precision, dimension(6) :: range
        double precision, dimension(3) :: lcrange1,lcrange2,temp1
        double precision :: rmxin,xmnin,xmxin,ymnin,ymxin

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        xl = Scxl
        rmxin = rmxi/xl
        xmnin = xmni/xl
        xmxin = xmxi/xl
        ymnin = ymni/xl
        ymxin = ymxi/xl

        eps = 0.0001
        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange2(i) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        rmax = max(rmxin,sqrt(xmax*xmax+ymax*ymax))

        hxx = (xmax-xmin)/nxcell
        hyy = (ymax-ymin)/nycell
        hrr = rmax/nrcell
        count = 0
        count4 = 0
        do i = 1, this%Nptlocal
          rad = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+this%Pts1(3,i)*&
                     this%Pts1(3,i))
          ir = int(rad/hrr) + 1
          count(ir) = count(ir) + 1
          ix = int((this%Pts1(1,i)-xmin)/hxx) + 1
          iy = int((this%Pts1(3,i)-ymin)/hyy) + 1
          count4(ix,iy) = count4(ix,iy) + 1
        enddo

        call MPI_REDUCE(count,glcount,nrcell,MPI_INTEGER,MPI_SUM,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
!          do i = 2, nrcell
!            glcount(i) = glcount(i) + glcount(i-1)
!          enddo
!          print*,"glcount: ",glcount(nrcell)

!          open(unit=nfile,file='RadDens2.data',status='unknown',&
!               position='append',form='unformatted')
          open(unit=nfile,file='RadDens2.data',status='unknown',&
               position='append')
          write(nfile,*)nstep
          do i = 1, nrcell
            !write(nfile,100)float(glcount(i))/glcount(nrcell)
            write(nfile,100)(i*hrr-0.5*hrr)*xl,float(glcount(i))/nptot
          enddo
          close(nfile)
        endif

        allocate(nptlist(0:np-1))
        nptlist = 0
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
          open(unit=nfile3,file='Yprof2.data',status='unknown',&
                position='append')
          testsum = 0
          write(nfile3,*)nstep
          do j = 1, nycell
            sumy = 0
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sumy = sumy + sendbuf4(3,ic)
            enddo
            write(nfile3,100)((j*hyy-0.5*hyy)+ymin)*xl,float(sumy)/nptot
            testsum = testsum + sumy
          enddo
          close(nfile3)
!          print*,"testsumy: ",testsum

          open(unit=nfile3,file='Xprof2.data',status='unknown',&
                position='append')
          testsum = 0
          write(nfile3,*)nstep
          do i = 1, nxcell
            sumx = 0
            do j = 1, nycell
              ic = (j-1)*nxcell + i
              sumx = sumx + sendbuf4(3,ic)
            enddo
            write(nfile3,100)((i*hxx-0.5*hxx)+xmin)*xl,float(sumx)/nptot
            testsum = testsum + sumx
          enddo
          close(nfile3)
!          print*,"testsumx: ",testsum
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

100     format(2(1x,e14.7))

        deallocate(nptlist)

        end subroutine dens1d_Output

      end module AdvOutputclass
