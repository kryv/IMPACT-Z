!!GS20150112 add marker "-23" to enable fort.* output at each marker position
!!QZ20130607 energy dependent in dipole
!!20130607  the 6th entry (after bending angle) in the "4" line is "beta-gamma" of the ideal reference particle
!!QZ20120603 add output for dipoe "4", slit "-13", stripper "-11", corretor "-21", element shift "-25", element rotation "-22"
!!QZ20120530 if dipole ID=500, no fringe fields calculation, 400 & 600 only for entrance & exit edge, respectively 
!!QZ20120210 introduced "tmpmultip(12)" instead of sharing with "tmpdipole"
!!QZ06022011 add tmpmultip(11) = val10(i) and tmpmultip(12) = val11(i) 
!!QZ06022011 and tmpmultip(13) = val12(i) and tmpmultip(14) = val13(i) for Multipoles
!!based on 03/22/2010
!----------------------------------------------------------------
! (c) Copyright, 2002 by the Regents of the University of California.
! AccSimulatorclass: Linear accelerator simulator class in CONTROL layer.
! Version: 2.0
! Author: Ji Qiang, LBNL, 1/05/02
! Description: This class defines functions to set up the initial beam 
!              particle distribution, field information, computational
!              domain, beam line element lattice and run the dynamics
!              simulation through the system.
! Comments:
!----------------------------------------------------------------
      module AccSimulatorclass
        use Pgrid2dclass
        use CompDomclass
        use FieldQuantclass
        use BeamLineElemclass
        use Ptclmgerclass
        use BeamBunchclass
        use Timerclass
        use Inputclass
        use Outputclass
        use Dataclass
        use PhysConstclass
        use NumConstclass
        use Distributionclass
        use Besselclass
        use MTrndclass
        implicit none
        !# of phase dim., num. total and local particles, int. dist. 
        !and restart switch, error study switch, substep for space-charge
        !switch
        integer, private :: Dim, Np, Nplocal,Flagdist,Rstartflg,Flagerr,&
                            Flagsubstep 

        !# of num. total x, total and local y mesh pts., type of BC, 
        !# of beam elems, type of integrator.
        integer, private :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal,Flagbc,&
                            Nblem,Flagmap,Flagdiag

        !# of processors in column and row direction.
        integer, private :: npcol, nprow

        !beam current, kin. energy, part. mass, and charge.
        double precision, private :: Bcurr,Bkenergy,Bmass,Bcharge,Bfreq,&
                                     Perdlen

        !conts. in init. dist.
        double precision, private, dimension(21) :: distparam

        !1d logical processor array.
        type (Pgrid2d), private :: grid2d

        !beam particle object and array.
        type (BeamBunch), private :: Bpts

        !beam charge density and field potential arrays.
        type (FieldQuant), private :: Potential

        !geometry object.
        type (CompDom), private :: Ageom

        !beam line element array.
        type (BPM),target,dimension(Nbpmmax) :: beamln0
        type (DriftTube),target,dimension(Ndriftmax) :: beamln1
        type (Quadrupole),target,dimension(Nquadmax) :: beamln2
        type (DTL),target,dimension(Ndtlmax) :: beamln3
        type (CCDTL),target,dimension(Nccdtlmax) :: beamln4
        type (CCL),target,dimension(Ncclmax) :: beamln5
        type (SC),target,dimension(Nscmax) :: beamln6
        type (ConstFoc),target,dimension(Ncfmax) :: beamln7
        type (SolRF),target,dimension(Nslrfmax) :: beamln8
        type (Sol),target,dimension(Nslmax) :: beamln9
        type (Dipole),target,dimension(Ndipolemax) :: beamln10
        type (EMfld),target,dimension(Ncclmax) :: beamln11
        type (Multipole),target,dimension(Nquadmax) :: beamln12
        type (RFQ),target,dimension(Ncclmax) :: beamln13
        type (BeamLineElem),private,dimension(Nblemtmax)::Blnelem
        !beam line element period.
        interface construct_AccSimulator
          module procedure init_AccSimulator
        end interface

        !//total # of charge state
        integer :: nchrg
        !//current list of charge state.
        !//charge/mass list of charge state.
        double precision, dimension(100) :: currlist,qmcclist
        !//number of particles of charge state.
        integer, dimension(100) :: nptlist
        integer, allocatable, dimension(:) :: nptlist0
        double precision, allocatable, dimension(:) :: currlist0,qmcclist0
!!QZ output emittance spec
        real :: outputflag
        
        integer :: isrd = 1010          ! Seed for random number generators
        
      contains
        !set up objects and parameters.
        subroutine init_AccSimulator(time)
        implicit none
        include 'mpif.h'
        integer :: i,test1,test2,j
        integer :: myid,myidx,myidy,ierr,inb,jstp
        integer, allocatable, dimension(:) :: bnseg,bmpstp,bitype
        double precision, allocatable, dimension(:) :: blength,val1,&
        val2, val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24
        double precision :: time
        double precision :: t0
        double precision :: z,xrad,yrad,phsini
        double precision, dimension(2) :: tmpdr 
        double precision, dimension(5) :: tmpcf 
        double precision, dimension(8) :: tmpbpm 
        double precision, dimension(9) :: tmpquad
        double precision, dimension(10) :: tmpdipole
		double precision, dimension(14) :: tmpmultip  
        double precision, dimension(11) :: tmprf
        double precision, dimension(12) :: tmpslrf
        double precision, dimension(14) :: tmp13
        double precision, dimension(25) :: tmpdtl
        integer :: iqr,idr,ibpm,iccl,iccdtl,idtl,isc,icf,islrf,isl,idipole,&
                   iemfld,myrank,imultpole,irfq
        real*8 :: wk,x1,x2,xx,aa,besi1,besi2

        !start up MPI.
        call init_Input(time)

        ! initialize Timer.
        call construct_Timer(0.0d0)

        call starttime_Timer(t0)

        !Flagmap = 0

        nptlist = 0
        currlist = 0.0
        qmcclist = 0.0
!-------------------------------------------------------------------
! get all global input parameters.
        call in_Input(Dim,Np,Nx,Ny,Nz,Flagbc,Flagdist,Rstartflg,&
              Flagmap,distparam,21,Bcurr,Bkenergy,Bmass,Bcharge,&
		Bfreq,xrad,yrad,Perdlen,Nblem,npcol,nprow,Flagerr,Flagdiag,&
		Flagsubstep,phsini,nchrg,nptlist,currlist,qmcclist,outputflag)
!!QZ        Flagsubstep,phsini,nchrg,nptlist,currlist,qmcclist)
 
        allocate(nptlist0(nchrg))
        allocate(currlist0(nchrg))
        allocate(qmcclist0(nchrg))
        do i = 1, nchrg
          nptlist0(i) =  nptlist(i)
          currlist0(i) =  currlist(i)
          qmcclist0(i) =  qmcclist(i)
        enddo
!        print*,"npt0: ",nptlist0
!        print*,"qmcc0: ",qmcclist0
!-------------------------------------------------------------------
! construct 2D logical processor Cartesian coordinate
        call construct_Pgrid2d(grid2d,MPI_COMM_WORLD,nprow,npcol)
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        if(myid.eq.0) then
          print*,"Start simulation:"
        endif

        !construct Constants class.
        call construct_PhysConst(Bfreq)

!-------------------------------------------------------------------
! construct computational domain CompDom class and get local geometry 
! information on each processor.
        !if(Rstartflg.eq.1) then
        !  call ingeom_Output(1500,z,inb,jstp,nprow,npcol,Ageom,Nx,Ny,Nz,&
        !                    myidx,myidy)
        !  if(myid.eq.0) print*,"rstart at: ",z,inb,jstp
        !  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !else
          !xrad = 0.1363243029*0.2
          call construct_CompDom(Ageom,distparam,21,Flagdist,&
               Nx,Ny,Nz,grid2d,nprow,npcol,Flagbc,xrad,yrad,Perdlen)
        !endif

!-------------------------------------------------------------------
! initialize Data class.
        call init_Data()

!-------------------------------------------------------------------
! construct BeamBunch class.
        call construct_BeamBunch(Bpts,Bcurr,Bkenergy,Bmass,Bcharge,&
                            Np,phsini)

!-------------------------------------------------------------------
! sample initial particle distribution.
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
        call mtseed(isrd+21+myid)        
        if(Rstartflg.eq.1) then
          !call phasein_Output(1500,Bpts)
          !call phasein2_Output(myrank+31,Bpts)
          call inpoint_Output(myid+31,Bpts,z,inb,jstp,nprow,npcol,&
               Ageom,Nx,Ny,Nz,myidx,myidy,Np)
          if(myid.eq.0) print*,"rstart at: ",z,inb,jstp
        else
          call sample_Dist(Bpts,distparam,21,Flagdist,Ageom,grid2d,Flagbc,&
                           nchrg,nptlist0,qmcclist0,currlist0)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !if(myid.eq.0) print*,"pass generating initial distribution..."
!        print*,"qmcclist1: ",qmcclist0

        !get local particle number and mesh number on each processor.
        call getnpt_BeamBunch(Bpts,Nplocal)

!-------------------------------------------------------------------
! construct FieldQuant class objects.
        call construct_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d)

!-------------------------------------------------------------------
! construct beam line elements.
        allocate(blength(Nblem),bnseg(Nblem),bmpstp(Nblem))
        allocate(bitype(Nblem))
        allocate(val1(Nblem),val2(Nblem),val3(Nblem),val4(Nblem))
        allocate(val5(Nblem),val6(Nblem),val7(Nblem),val8(Nblem))
        allocate(val9(Nblem),val10(Nblem),val11(Nblem),val12(Nblem))
        allocate(val13(Nblem),val14(Nblem),val15(Nblem),val16(Nblem))
        allocate(val17(Nblem),val18(Nblem),val19(Nblem),val20(Nblem))
        allocate(val21(Nblem),val22(Nblem),val23(Nblem),val24(Nblem))

        call in_Input(Nblem,blength,bnseg,bmpstp,bitype,val1,val2,val3,&
        val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24)

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
        idipole = 0
        iemfld = 0
        imultpole = 0
        irfq = 0
        do i = 1, Nblem
          if(bitype(i).lt.0) then
            ibpm = ibpm + 1
            call construct_BPM(beamln0(ibpm),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpbpm(1) = 0.0
            tmpbpm(2) = val1(i)
            tmpbpm(3) = val2(i)
            tmpbpm(4) = val3(i)
            tmpbpm(5) = val4(i)
            tmpbpm(6) = val5(i)
            tmpbpm(7) = val6(i)
            tmpbpm(8) = val7(i)
            call setparam_BPM(beamln0(ibpm),tmpbpm)
            Blnelem(i) = assign_BeamLineElem(beamln0(ibpm))
          else if(bitype(i).eq.0) then
            idr = idr + 1
            call construct_DriftTube(beamln1(idr),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpdr(1) = 0.0
            tmpdr(2) = val1(i)
            call setparam_DriftTube(beamln1(idr),tmpdr)
            Blnelem(i) = assign_BeamLineElem(beamln1(idr))
          else if(bitype(i).eq.1) then
            iqr = iqr + 1
            call construct_Quadrupole(beamln2(iqr),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpquad(1) = 0.0
            tmpquad(2) = val1(i)
            tmpquad(3) = val2(i)
            tmpquad(4) = val3(i)
            tmpquad(5) = val4(i)
            tmpquad(6) = val5(i)
            tmpquad(7) = val6(i)
            tmpquad(8) = val7(i)
            tmpquad(9) = val8(i)
            call setparam_Quadrupole(beamln2(iqr),tmpquad)
            Blnelem(i) = assign_BeamLineElem(beamln2(iqr))
          else if(bitype(i).eq.2) then
            icf = icf + 1
            call construct_ConstFoc(beamln7(icf),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpcf(1) = 0.0
            tmpcf(2) = val1(i)
            tmpcf(3) = val2(i)
            tmpcf(4) = val3(i)
            tmpcf(5) = val4(i)
            call setparam_ConstFoc(beamln7(icf),tmpcf)
            Blnelem(i) = assign_BeamLineElem(beamln7(icf))
          else if(bitype(i).eq.3) then
            isl = isl + 1
            call construct_Sol(beamln9(isl),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpquad(1) = 0.0
            tmpquad(2) = val1(i)
            tmpquad(3) = val2(i)
            tmpquad(4) = val3(i)
            tmpquad(5) = val4(i)
            tmpquad(6) = val5(i)
            tmpquad(7) = val6(i)
            tmpquad(8) = val7(i)
            tmpquad(9) = val8(i)
            call setparam_Sol(beamln9(isl),tmpquad)
            Blnelem(i) = assign_BeamLineElem(beamln9(isl))
          else if(bitype(i).eq.4) then
            idipole = idipole + 1
            call construct_Dipole(beamln10(idipole),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpdipole(1) = 0.0
            tmpdipole(2) = val1(i)
            tmpdipole(3) = val2(i)
            tmpdipole(4) = val3(i)
            tmpdipole(5) = val4(i)
            tmpdipole(6) = val5(i)
            tmpdipole(7) = val6(i)
            tmpdipole(8) = val7(i)
            tmpdipole(9) = val8(i)
            tmpdipole(10) = val9(i)
            call setparam_Dipole(beamln10(idipole),tmpdipole)
            Blnelem(i) = assign_BeamLineElem(beamln10(idipole))
          else if(bitype(i).eq.5) then
            imultpole = imultpole + 1
            call construct_Multipole(beamln12(imultpole),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpmultip(1) = 0.0
            tmpmultip(2) = val1(i)
            tmpmultip(3) = val2(i)
            tmpmultip(4) = val3(i)
            tmpmultip(5) = val4(i)
            tmpmultip(6) = val5(i)
            tmpmultip(7) = val6(i)
            tmpmultip(8) = val7(i)
            tmpmultip(9) = val8(i)		!!Radius
            tmpmultip(10) = val9(i)		!!dX
			tmpmultip(11) = val10(i)	!!dY
			tmpmultip(12) = val11(i)	!!Rx
			tmpmultip(13) = val12(i)	!!Ry
			tmpmultip(14) = val13(i)	!!Rz
!!QZ06022011 add tmpmultip(11) = val10(i) and tmpmultip(12) = val11(i)
!!QZ06022011 add tmpmultip(13) = val12(i) and tmpmultip(14) = val13(i)
            call setparam_Multipole(beamln12(imultpole),tmpmultip)
            Blnelem(i) = assign_BeamLineElem(beamln12(imultpole))
          else if(bitype(i).eq.101) then
            idtl = idtl + 1
            call construct_DTL(beamln3(idtl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpdtl(1) = 0.0
            tmpdtl(2) = val1(i) 
            tmpdtl(3) = val2(i) 
            tmpdtl(4) = val3(i) 
            tmpdtl(5) = val4(i) 
            tmpdtl(6) = val5(i) 
            tmpdtl(7) = val6(i) 
            tmpdtl(8) = val7(i) 
            tmpdtl(9) = val8(i) 
            tmpdtl(10) = val9(i) 
            tmpdtl(11) = val10(i) 
            tmpdtl(12) = val11(i) 
            tmpdtl(13) = val12(i) 
            tmpdtl(14) = val13(i) 
            tmpdtl(15) = val14(i) 
            tmpdtl(16) = val15(i) 
            tmpdtl(17) = val16(i) 
            tmpdtl(18) = val17(i) 
            tmpdtl(19) = val18(i) 
            tmpdtl(20) = val19(i) 
            tmpdtl(21) = val20(i) 
            tmpdtl(22) = val21(i) 
            tmpdtl(23) = val22(i) 
            tmpdtl(24) = val23(i) 
            tmpdtl(25) = val24(i) 
            call setparam_DTL(beamln3(idtl),tmpdtl)
            Blnelem(i) = assign_BeamLineElem(beamln3(idtl))
          else if(bitype(i).eq.102) then
            iccdtl = iccdtl + 1
            call construct_CCDTL(beamln4(iccdtl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmprf(1) = 0.0
            tmprf(2) = val1(i) 
            tmprf(3) = val2(i) 
            tmprf(4) = val3(i) 
            tmprf(5) = val4(i) 
            tmprf(6) = val5(i) 
            tmprf(7) = val6(i)
            tmprf(8) = val7(i) 
            tmprf(9) = val8(i) 
            tmprf(10) = val9(i) 
            tmprf(11) = val10(i) 
            call setparam_CCDTL(beamln4(iccdtl),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln4(iccdtl))
          else if(bitype(i).eq.103) then
            iccl = iccl + 1
            call construct_CCL(beamln5(iccl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmprf(1) = 0.0
            tmprf(2) = val1(i) 
            tmprf(3) = val2(i) 
            tmprf(4) = val3(i) 
            tmprf(5) = val4(i) 
            tmprf(6) = val5(i) 
            tmprf(7) = val6(i) 
            tmprf(8) = val7(i) 
            tmprf(9) = val8(i) 
            tmprf(10) = val9(i) 
            tmprf(11) = val10(i) 
            call setparam_CCL(beamln5(iccl),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln5(iccl))
          else if(bitype(i).eq.104) then
            isc = isc + 1
            call construct_SC(beamln6(isc),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmprf(1) = 0.0
            tmprf(2) = val1(i) 
            tmprf(3) = val2(i) 
            tmprf(4) = val3(i) 
            tmprf(5) = val4(i) 
            tmprf(6) = val5(i) 
            tmprf(7) = val6(i)
            tmprf(8) = val7(i) 
            tmprf(9) = val8(i) 
            tmprf(10) = val9(i) 
            tmprf(11) = val10(i) 
            call setparam_SC(beamln6(isc),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln6(isc))
          else if(bitype(i).eq.105) then
            islrf = islrf + 1
            call construct_SolRF(beamln8(islrf),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpslrf(1) = 0.0
            tmpslrf(2) = val1(i) 
            tmpslrf(3) = val2(i) 
            tmpslrf(4) = val3(i) 
            tmpslrf(5) = val4(i) 
            tmpslrf(6) = val5(i) 
            tmpslrf(7) = val6(i) 
            tmpslrf(8) = val7(i) 
            tmpslrf(9) = val8(i) 
            tmpslrf(10) = val9(i) 
            tmpslrf(11) = val10(i) 
            tmpslrf(12) = val11(i) 
            call setparam_SolRF(beamln8(islrf),tmpslrf)
            Blnelem(i) = assign_BeamLineElem(beamln8(islrf))
          else if(bitype(i).eq.106) then
            irfq = irfq + 1
            wk = 2*PI/blength(i)
            x1 = wk*val5(i)
            x2 = wk*val6(i)*val5(i)
            besi1 = bessi0(x1)
            besi2 = bessi0(x2)
            xx = (besi1+besi2)/(val6(i)**2*besi1+besi2)
            aa = (val6(i)**2-1)/(val6(i)**2*besi1+besi2)
            call construct_RFQ(beamln13(irfq),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i),aa,xx)
            tmpslrf(1) = 0.0
            tmpslrf(2) = val1(i) 
            tmpslrf(3) = val2(i) 
            tmpslrf(4) = val3(i) 
            tmpslrf(5) = val4(i) 
            tmpslrf(6) = val5(i) 
            tmpslrf(7) = val6(i) 
            tmpslrf(8) = val7(i) 
            tmpslrf(9) = val8(i) 
            tmpslrf(10) = val9(i) 
            tmpslrf(11) = val10(i) 
            tmpslrf(12) = val11(i) 
            call setparam_RFQ(beamln13(irfq),tmpslrf)
            Blnelem(i) = assign_BeamLineElem(beamln13(irfq))
          else if(bitype(i).eq.110) then
            iemfld = iemfld + 1
            call construct_EMfld(beamln11(iemfld),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmp13(1) = 0.0
            tmp13(2) = val1(i)
            tmp13(3) = val2(i)
            tmp13(4) = val3(i)
            tmp13(5) = val4(i)
            tmp13(6) = val5(i)
            tmp13(7) = val6(i)
            tmp13(8) = val7(i)
            tmp13(9) = val8(i)
            tmp13(10) = val9(i)
            tmp13(11) = val10(i)
            tmp13(12) = val11(i)
            tmp13(13) = val12(i)
            tmp13(14) = val13(i)
            call setparam_EMfld(beamln11(iemfld),tmp13)
            Blnelem(i) = assign_BeamLineElem(beamln11(iemfld))
            !print*,"tmp13_12: ",tmp13(12),tmp13(13)
          else
          endif 
        enddo
!-------------------------------------------------------------------
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !if(myid.eq.0) print*,"pass setting up lattice..."

        deallocate(blength,bnseg,bmpstp,bitype)
        deallocate(val1,val2,val3,val4,val5,val6,val7,val8,val9)
        deallocate(val10,val11,val12,val13,val14,val15,val16)
        deallocate(val17,val18,val19,val20,val21,val22,val23)
        deallocate(val24)

        t_init = t_init + elapsedtime_Timer(t0)

        end subroutine init_AccSimulator

        !Run beam dynamics simulation through accelerator.
        subroutine run_AccSimulator()
        implicit none
        include 'mpif.h'
        integer :: i,j,bnseg,bmpstp,bitype
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,npx,npy,&
                   totnp,ierr,nylcr,nzlcr
        integer :: nbal,ibal,nstep,ifile
        integer :: iend,ibalend,nstepend,tmpfile,nfile
        integer, dimension(3) :: lcgrid
        integer, allocatable, dimension(:) :: lctabnmx,lctabnmy
        integer, allocatable, dimension(:,:,:) :: temptab
        double precision :: z0,z,tau1,tau2,blength,t0,t1
        double precision, allocatable, dimension(:,:) :: lctabrgx, lctabrgy
        double precision, dimension(6) :: lcrange, range, ptrange,ptref
        double precision, dimension(3) :: msize
        double precision :: hy,hz,ymin,zmin,zend,piperad,zedge
        double precision :: tmp1,tmp2,tmp3,tmp4,rfile
        double precision, allocatable, dimension(:,:,:) :: chgdens
        double precision, allocatable, dimension(:,:,:) :: besscoef
        double precision, allocatable, dimension(:,:) :: bessnorm,gml
        integer, allocatable, dimension(:) :: modth,pydisp
        integer :: nmod,k,ii,jj
        !double precision :: sumtest, sumtest2, sumtest3
        double precision, dimension(8) :: drange
        double precision, dimension(3) :: al0,ga0,epson0
        double precision :: realSamplePeriod,tg,tv,gam,piperad2
        integer :: nsubstep,integerSamplePeriod,Flagbctmp
        double precision :: beta0,gamma0,gambetz
        double precision, allocatable, dimension(:,:) :: tmpptcs
        double precision :: rwkinq,rwkinf,avgw

        double precision :: dpi,ang0,hd0,hd1,dstr1,dstr2,angF,tanphiF,tanphiFb,&
           angB,tanphiB,tanphiBb,hF,hB,qm0,qmi,deltax,deltat,dstr1i,&
           hdi,gambet,beta00,gamma00,gambet00
        double precision, dimension(6) :: ptarry,ptarry2
        double precision, dimension(10) :: dparam
        integer :: ipt
        double precision :: ri,r0,thetap,rr0,rra,cosang0,sinang0,drr
        !ANL stripper model parameters
! ... Incident beam energy and average energy loss
        real*8      E0             ! Incident beam energy
        real*8      El0            ! Average energy loss.
!        parameter  (E0 = 13.3d0, El0 = 1.275d-1)
!       parameter  (E0 =  85.d0, El0 = 3.292d0)
! ... For the stripper generator
        real*8      El             ! Energy lost by the ion at exit
        real*8      Thx            ! Angle projection on (X,Z)
        real*8      Thy            ! Angle projection on (Y,Z)
! ... Variables calculated from generated ones
        real*8      En             ! Energy of the ion at exit
        real*8      Th                 ! Resultant transverse angle
! ... Generating loop variables
        integer     nev                ! number of events to generate
        integer     iev                ! event index
        integer :: ichg,nptmp2,nptottmp,nplistlc
        double precision :: xradmin,xradmax,yradmin,yradmax
        !for bending magnet
        double precision :: hgap,tmppt,psi1,psi2,gamn
        real*8 :: xrad,yrad,ma,wk,zz
        double precision :: qmrel,dlength
        double precision  :: xerr,yerr,anglerrx,anglerry,anglerrz

!-------------------------------------------------------------------
! prepare initial parameters, allocate temporary array.
        iend = 0
        ibalend = 0
        nstepend = 0
        zend = 0.0

        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid2d,totnp,npy,npx)
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
        ibal = ibalend
        nstep = nstepend
        z = zend

        if(Rstartflg.eq.1) then
          call restart_AccSimulator(iend,nstepend,ibalend,zend)
          ibal = ibalend
          nstep = nstepend
          z = zend
        else
          if(Flagdiag.eq.1) then
            call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
          elseif(Flagdiag.eq.2) then
            call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
          endif
        endif

        allocate(chgdens(1,1,1))
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
          if(npcol.gt.1) then
            Nylocal = lcgrid(2) + 2
          else
            Nylocal = lcgrid(2)
          endif
          if(nprow.gt.1) then
            Nzlocal = lcgrid(3) + 2
          else
            Nzlocal = lcgrid(3)
          endif

          if(nprow.gt.1) then
            nzlcr = Nzlocal-2
          else
            nzlcr = Nzlocal
          endif
          if(npcol.gt.1) then
            nylcr = Nylocal-2
          else
            nylcr = Nylocal
          endif
          if(myidy.eq.(npcol-1)) nylcr = nylcr - 1
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

!        print*,"qmcclist: ",nchrg,qmcclist0(1),qmcclist0(2),nptlist0(1),nptlist0(2)
!-------------------------------------------------------------------
! start looping through 'Nblem' beam line elements.
        tmpfile = 0
        do i = iend+1, Nblem
          !if(myid.eq.0) then
          !  print*,"enter elment: ",i
          !endif

          call getparam_BeamLineElem(Blnelem(i),blength,bnseg,bmpstp,&
                                     bitype)
          call getradius_BeamLineElem(Blnelem(i),piperad,piperad2)
          if(Flagerr.ne.1 .or. Flagmap.eq.1 .or. bitype.eq.3) then
            xradmin = -piperad
            xradmax = piperad
            yradmin = -piperad2
            yradmax = piperad2
          else
            call geterr_BeamLineElem(Blnelem(i),xerr,yerr,anglerrx,&
                                     anglerry,anglerrz)
            xradmin = -piperad + xerr
            xradmax = piperad + xerr
            yradmin = -piperad2 + yerr
            yradmax = piperad2 + yerr
          endif
          !if(myid.eq.0) then
          !  print*,"element: ",i,piperad,piperad2
          !endif  
          nsubstep = bmpstp
          nfile = 0
          tau1 = 0.0
          if(bitype.ge.0) tau1 = 0.5*blength/bnseg
          tau2 = 2.0*tau1

          !for RFQ, the aperature size varies within a period
          if(bitype.eq.106) then
            wk = 2*PI/blength
            call getparam_BeamLineElem(Blnelem(i),1,zedge)
            call getparam_BeamLineElem(Blnelem(i),7,ma)
            zz = z - zedge
            xrad = piperad*((ma+1)/2-(ma-1.0)/2*cos(wk*zz))
            yrad = piperad*((ma+1)/2+(ma-1.0)/2*cos(wk*zz))
            piperad = xrad
            piperad2 = yrad
          endif
!-------------------------------------------------------------------
! read in the on axis E field for rf cavities.
!          if(myid.eq.0) print*,"bitype: ",bitype
          if(bitype.gt.100) then
            call getparam_BeamLineElem(Blnelem(i),5,rfile)
            nfile = int(rfile + 0.1)
            ifile = nfile
!            if(myid.eq.0) print*,"ifile: ",nfile
            if(ifile.ne.tmpfile)then
              call starttime_Timer(t1)
              !for linear map integrator
              if(Flagmap.eq.1) then
                !read in Ez, Ez', Ez'' on the axis
                call read1_Data(ifile)
              else
                !read in Er, Ez, H_theta on r-z grid 
                !call read2_Data(ifile)
                !read in Fourier coefficients of Ez on the axis
                !call read3_Data(ifile)
                if(bitype.eq.110) then
                  call read4_Data(ifile)
                else
                  call read3_Data(ifile)
                endif
                t_load3dfile = t_load3dfile + elapsedtime_Timer(t1)
              endif
              tmpfile=ifile
            endif
          endif

!-------------------------------------------------------------------
! print out beam information using BPM 
          if(bitype.eq.-1) then
            call shift_BPM(Bpts%Pts1,bitype,Nplocal,Np)
          endif
!QZ		07/23/07
!QZ		type=-26: shift both the transverse centroid positions and angles to 0
          if(bitype.eq.-26) then
            call shift_BPM(Bpts%Pts1,bitype,Nplocal,Np)
          endif
!QZ
          if(bitype.eq.-2) then
            !call phase_Output(bmpstp,Bpts)
            !call phaseleda_Output(bmpstp,Bpts)
            call getparam_BeamLineElem(Blnelem(i),drange)
            realSamplePeriod = drange(2)
            integerSamplePeriod = realSamplePeriod
            call phase_Output(bmpstp,Bpts,integerSamplePeriod)
          else if(bitype.eq.-3) then
            call getparam_BeamLineElem(Blnelem(i),drange)
            call accdens1d_Output(nstep,8,Bpts,Np,drange(2),-drange(3),&
            drange(3),-drange(5),drange(5))
          else if(bitype.eq.-4) then
            call getparam_BeamLineElem(Blnelem(i),drange)
            call dens1d_Output(nstep,8,Bpts,Np,drange(2),-drange(3),&
            drange(3),-drange(5),drange(5))
          else if(bitype.eq.-5) then
            call getparam_BeamLineElem(Blnelem(i),drange)
            call dens2d_Output(nstep,8,Bpts,Np,-drange(3),drange(3),&
            -drange(4),drange(4),-drange(5),drange(5),-drange(6),drange(6),&
            -drange(7),drange(7),-drange(8),drange(8))
          else if(bitype.eq.-6) then
            call dens3d_Output(nstep,8,Bpts,Np,-drange(3),drange(3),&
              -drange(5),drange(5),-drange(7),drange(7))
          else if(bitype.eq.-7) then
            !output all particles in 6d phase space in file "xxnstep"
            !call phaseout_Output(nstep,Bpts)
            !output all geomtry information in file "xxnstep".
            !call outgeom_Output(nstep,z,i,j,npx,npy,Ageom)
            call outpoint_Output(myid+31,Bpts,z,i,j,npx,npy,Ageom)
          else if(bitype.eq.-10) then
            !mismatch the beam at given location.
            !here, drange(3:8) stores the mismatch factor.
            call getparam_BeamLineElem(Blnelem(i),drange)
            call scale_BPM(Bpts%Pts1,Nplocal,drange(3),&
            drange(4),drange(5),drange(6),drange(7),drange(8))
          else if(bitype.eq.-11) then !MSU stripper model
            deallocate(tmpptcs)
            allocate(tmpptcs(8,Nplocal))
            gam = -Bpts%refptcl(6)
            gambetz = sqrt(gam**2-1.0)
            nfile = bmpstp
            beta0 = sqrt(1.0-(1.0/gam)**2)
            !beta0 = 0.1581611
            do ii = 1, Nplocal
              tmpptcs(1,ii) = Bpts%Pts1(9,ii)
              tmpptcs(2,ii) = Bpts%Pts1(1,ii)*Scxl*100
              tmpptcs(3,ii) = Bpts%Pts1(3,ii)*Scxl*100
              tmpptcs(4,ii) = Bpts%Pts1(5,ii)*90/asin(1.0)
              !gambetz = sqrt((gam-Bpts%Pts1(6,ii))**2-Bpts%Pts1(2,ii)**2-&
              !               Bpts%Pts1(4,ii)**2-1.0)
              tmpptcs(5,ii) = Bpts%Pts1(2,ii)/gambetz * 1000
              tmpptcs(6,ii) = Bpts%Pts1(4,ii)/gambetz * 1000
              tmpptcs(7,ii) = -Bpts%Pts1(6,ii)/(gam-1)*100
              tmpptcs(8,ii) = Bpts%Pts1(7,ii)*Bpts%mass
            enddo
            !test only
            !open(11,file='incord.dat')
!            open(12,file='outcord.dat')
            !beta0 = 0.1581611
            !print*,"before reading particles: ",nfile,Nplocal,bmpstp
            !do ii = 1, Nplocal
            !  read(11,*)(tmpptcs(jj,ii),jj=1,8)
            !enddo
            !print*,"after reading particles: " 
            call MStripF_BPM(beta0,tmpptcs,nfile,Nplocal,Nplocal,rwkinf,&
                 rwkinq,nchrg,nptlist,currlist,qmcclist)
            !print*,"after stripper ",rwkinf,rwkinq,nchrg 
!            do ii = 1, Nplocal
!              write(12,13)(tmpptcs(jj,ii),jj=1,8)
!            enddo
13          format(1x,f6.0,1x,7e15.7)
            deallocate(nptlist0)
            deallocate(currlist0)
            deallocate(qmcclist0)
            allocate(nptlist0(nchrg))
            allocate(currlist0(nchrg))
            allocate(qmcclist0(nchrg))
            do ii = 1, nchrg
              nptlist0(ii) =  nptlist(ii)
              currlist0(ii) =  currlist(ii)
              qmcclist0(ii) =  qmcclist(ii)/Bpts%mass
            enddo

            gam = (1.0+rwkinq/Bmass)
            gambetz = sqrt(gam**2-1.0)
            avgw = 0.0
            do ii = 1, Nplocal
               Bpts%Pts1(6,ii) = -(gam-1.0d0)*tmpptcs(7,ii)/100
               !gambetz = sqrt((gam-Bpts%Pts1(6,ii))**2-1.0) /  &
               !sqrt(1.0+(tmpptcs(5,ii)/1000)**2+(tmpptcs(6,ii)/1000)**2)
               Bpts%Pts1(2,ii) = tmpptcs(5,ii)/1000*gambetz
               Bpts%Pts1(4,ii) = tmpptcs(6,ii)/1000*gambetz
               Bpts%Pts1(7,ii) = tmpptcs(8,ii)/Bpts%mass
               avgw = avgw + tmpptcs(7,ii)
            enddo
            avgw = avgw/Nplocal
            call chgupdate_BeamBunch(Bpts,nchrg,nptlist0,qmcclist0)

            !update reference particle information
            !print*,"Bcharge0: ",Bcharge,Bpts%refptcl(6)
            Bcharge = Bmass*qmcclist0(1)
            Bpts%Charge = Bcharge
            Bpts%refptcl(6) = -(1.0+rwkinq/Bmass)
            !print*,"Bcharge: ",Bcharge,Bpts%refptcl(6)
!!add 2013050603 ------
            if(Flagdiag.eq.1) then
            	call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
             elseif(Flagdiag.eq.2) then
            	call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
            endif
!!add 2013050603 ------
          else if(bitype.eq.-12) then !ANL stripper model
            nfile = bmpstp
            read(nfile,*)nchrg,E0,El0
            read(nfile,*)(nptlist(ichg),ichg=1,nchrg)
            read(nfile,*)(currlist(ichg),ichg=1,nchrg)
            read(nfile,*)(qmcclist(ichg),ichg=1,nchrg)

            !isrd=214748364 + 50*myid
            !call rnset (isrd,0,0)
            Bpts%refptcl(6) = -(1.0+(E0-El0)*1.0d6/Bmass)
            gam = -Bpts%refptcl(6)
            gambetz = sqrt(gam**2-1.0)

            nptottmp = 0
            do ichg = 1, nchrg
              nptottmp = nptottmp + nptlist(ichg) 
            enddo 
            nptmp2 = 0
            do ichg = 1, nchrg
              if(ichg.ne.nchrg) then
                nplistlc = Nplocal*nptlist(ichg)*1.d0/nptottmp
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
                Bpts%Pts1(7,iev) = qmcclist(ichg)
                Bpts%Pts1(8,iev) = currlist(ichg)/Scfreq/nptlist(i)*qmcclist(i)/abs(qmcclist(i))
              enddo
            enddo

            deallocate(nptlist0)
            deallocate(currlist0)
            deallocate(qmcclist0)
            allocate(nptlist0(nchrg))
            allocate(currlist0(nchrg))
            allocate(qmcclist0(nchrg))
            do ii = 1, nchrg
              nptlist0(ii) =  nptlist(ii)
              currlist0(ii) =  currlist(ii)
              qmcclist0(ii) =  qmcclist(ii)
            enddo

            !call chgupdate_BeamBunch(Bpts,nchrg,nptlist0,qmcclist0)
            Bcharge = Bmass*qmcclist0(1)
            Bpts%Charge = Bcharge

!            nev = 1000
!            !test only
!            ! ... Open output file
!            open(7,file='tst_strip.out',status='unknown')
!            do iev = 1, nev
!              ! ..... Generate the energy loss and angle projections
!              call gen_strip(2,0.d0,El,Thx,Thy)
!              ! ..... Calculate the ion exit energy
!              En  = E0 -El0 -El
! 
!              ! ----- Polar angle (Z)
!              Th = dsqrt(Thx**2 + Thy**2)
! 
!              ! ..... No units Change is needed : energy in MeV/u, angle in mrad
!              ! ..... Print generated variables
!              !        write(6,*) ' En, Th, Thx, Thy = ', En, Th, Thx, Thy
!              write(7,*) En, Th, Thx, Thy
!            enddo
!            close(7)
          else if(bitype.eq.-13) then !collimator slit
            call getparam_BeamLineElem(Blnelem(i),drange)
            xradmin = drange(3) 
            xradmax = drange(4) 
            yradmin = drange(5) 
            yradmax = drange(6) 
            call lostcountXY_BeamBunch(Bpts,Nplocal,Np,xradmin,&
                  xradmax,yradmin,yradmax)
            call chgupdate_BeamBunch(Bpts,nchrg,nptlist0,qmcclist0)
!!add 2013050603 ------
            if(Flagdiag.eq.1) then
            	call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
            elseif(Flagdiag.eq.2) then
            	call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
            endif
!!add 2013050603 ------
          else if(bitype.eq.-21)then
            !shift the beam centroid in the 6D phase space.
            !This element can be used to model steering magnets etc.
            !here, drange(3:8) stores the amount of shift.
            !drange(3); shift in x (m)
            !drange(4); shift in Px (rad)
            !drange(5); shift in y (m)
            !drange(6); shift in Py (rad)
            !drange(7); shift in z (deg)
            !drange(8); shift in Pz (MeV)
            call getparam_BeamLineElem(Blnelem(i),drange)
            call kick_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4),drange(5),&
            drange(6),drange(7),drange(8),-Bpts%refptcl(6),Bpts%Mass,&
            Bpts%Charge)
!!add 2013050603 ------
            if(Flagdiag.eq.1) then
            	call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
            elseif(Flagdiag.eq.2) then

            	call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
            endif
!!add 2013050603 ------
          else if(bitype.eq.-22)then
            !rotate the particle coordiante by drange(3) in x-y plane
            ! by drange(4) degrees in px-py plane. This function is 
            ! most used in order to tanfer the beam through a vertic bend.
            call getparam_BeamLineElem(Blnelem(i),drange)
            call xypxpyrot_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4))
!!add 2013050603 ------
            if(Flagdiag.eq.1) then
            	call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
            elseif(Flagdiag.eq.2) then
            	call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
            endif
!!add 2013050603 ------
          else if(bitype.eq.-23)then
            !!add 20150112 ------
            ! add saving data at -23 marker
            if(Flagdiag.eq.1 .or. Flagdiag.eq.3) then
            	call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
            elseif(Flagdiag.eq.2 .or. Flagdiag.eq.4) then
            	call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
            endif
          else if(bitype.eq.-25)then
!QZ		05/18/2007
!QZ		based on kick_BPM, shift same amount of value for all particles 
!QZ     to simulate misalignment of element.
            !here, drange(3:8) stores the amount of shift.
            !drange(3); shift in x (m)
            !drange(4); shift in Px (rad)
            !drange(5); shift in y (m)
            !drange(6); shift in Py (rad)
            !drange(7); shift in z (deg)
            !drange(8); shift in Pz (MeV)
            call getparam_BeamLineElem(Blnelem(i),drange)
            call qzmis_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4),drange(5),&
            drange(6),drange(7),drange(8),-Bpts%refptcl(6),Bpts%Mass,&
            Bpts%Charge)
!!add 2013050603 ------
            if(Flagdiag.eq.1) then
            	call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
            elseif(Flagdiag.eq.2) then
            	call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
            endif
!!add 2013050603 ------
	   endif
          if(bitype.eq.-99) then
            exit
          endif

!          print*,"nptlist: ",nptlist0,avgw
!          print*,"currlist: ",currlist0
!-------------------------------------------------------------------
! loop through 'bnseg' numerical segments in each beam element
! using 2 step symplectic integeration (ie. leap frog).
          zedge = z
          call setparam_BeamLineElem(Blnelem(i),1,zedge)
          !if(myid.eq.0) print*,"zedge: ",zedge,bitype
          if(bitype.ne.4) then  !//no bend

          do j = 1, bnseg
            if((Flagerr.eq.1).and.(Flagmap.eq.1)) then
              call geomerrL_BeamBunch(Bpts,Blnelem(i)) 
            end if
            if((Flagerr.eq.1).and.(bitype.eq.3)) then
              call geomerrL_BeamBunch(Bpts,Blnelem(i)) 
            end if
            !if(nstep.eq.1) then  !distort the distribution.
            !  call gammaepson_Dist(Bpts,al0,ga0,epson0)
            !  call Distort_Dist(Bpts,al0,ga0,epson0)
            !endif

!-------------------------------------------------------------------
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
              endif
!              call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
            endif

!-------------------------------------------------------------------
! escape the space charge calculation for 0 current case
            if(Bcurr.lt.1.0e-30)  then !no space-charge
!              call lostcount_BeamBunch(Bpts,Nplocal,Np,piperad,piperad2)
!              call chgupdate_BeamBunch(Bpts,nchrg,nptlist0,qmcclist0)
!              print*,"nchrg: ",nchrg,nptlist0,qmcclist0
              !print*,"diag1: ",z,Nplocal,Np
!              call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
              goto 200
            else !calculate space charge forces
              call conv1st_BeamBunch(Bpts,tau2,Nplocal,Np,ptrange,&
                                   Flagbc,Perdlen,piperad,piperad2)
              call chgupdate_BeamBunch(Bpts,nchrg,nptlist0,qmcclist0)
              !fix the global range for sub-cycle of space charge potential.
              if(Flagsubstep.eq.1) then
                if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                  ptrange(1) = 0.0
                  ptrange(2) = piperad
                  ptrange(3) = 0.0
                  ptrange(4) = 4*asin(1.0)
                else
                  ptrange(1) = -piperad
                  ptrange(2) = piperad
                  ptrange(3) = -piperad2
                  ptrange(4) = piperad2
                endif
                ptrange(5) = -Perdlen/2
                ptrange(6) = Perdlen/2
              endif

            ! get new boundary from the range of beam particles.
            if(Flagbc.eq.4) then
            else
              call update_CompDom(Ageom,ptrange,grid2d,Flagbc)
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
!                print*,"before ptcl manager: ...",myid
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                      Nplcmax,lcrange)
!                print*,"after ptcl manager: ...",myid
!              else if(Flagbc.eq.2) then
!                call ptsmv2perd_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
!                                Nplcmax,lcrange,Perdlen)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                Nplcmax,lcrange)
              endif
            endif
            ! assign new 'Nplocal' local particles on each processor.
            call setnpt_BeamBunch(Bpts,Nplocal)

            if((mod(j-1,nsubstep).eq.0).or.(Flagsubstep.ne.1)) then

            call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,&
                                npy)
            deallocate(chgdens)
            allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
            ! deposit particles onto grid to obtain charge density.
            call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,Ageom,&
                                  grid2d,chgdens,Flagbc,Perdlen)

!-------------------------------------------------------------------
! start load balance. (at the location of new space charge calculation)
            if((mod(ibal,nbal).eq.0).and.(totnp.gt.1) ) then
              call MPI_BARRIER(comm2d,ierr)
              if(myid.eq.0) then 
                print*," load balance! "
              endif

              call getlctabnm_CompDom(Ageom,temptab)
              lctabnmx(0:npx-1) = temptab(1,0:npx-1,0)
              lctabnmy(0:npy-1) = temptab(2,0,0:npy-1)
              call getmsize_CompDom(Ageom,msize)
              hy = msize(2)
              hz = msize(3) 
              call getrange_CompDom(Ageom,range)
              ymin = range(3)
              zmin = range(5)
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
              call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,npy)

              ! pass particles to local space domain on new processor.
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                      Nplcmax,lcrange)
!              else if(Flagbc.eq.2) then
!                call ptsmv2perd_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
!                                Nplcmax,lcrange,Perdlen)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                 Nplcmax,lcrange)
              endif
              ! assign new 'Nplocal' local particles on each processor.
              call setnpt_BeamBunch(Bpts,Nplocal)

              ! deposit particles onto grid to obtain charge density.
              call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,&
                                    Ageom,grid2d,chgdens,Flagbc,Perdlen)
            endif
            ibal = ibal + 1
!end load balance.
!-------------------------------------------------------------------

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
           
!-------------------------------------------------------------------------
! solve 3D Poisson's equation
            if(Flagbc.eq.1) then
              ! solve Poisson's equation using 3D isolated boundary condition.
              call update3O_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.2) then
              ! solve Poisson's equation using 2D isolated 1D periodic 
              ! boundary condition.
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              call update2O1P_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.3) then
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update3_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,&
              besscoef,bessnorm,gml,modth,nmod)
            else if(Flagbc.eq.4) then
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update4_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,Perdlen)
            else if(Flagbc.eq.5) then
              ! solve Poisson's equation using 2D rectangular pipe, 1D open
              ! boundary condition
              if(myidy.eq.(npy-1)) nylcr = nylcr-1
              call update5_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.6) then
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update6_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else
              print*,"no such boundary condition type!!!"
              stop
            endif

            endif

              call cvbkforth1st_BeamBunch(Bpts)
              if(totnp.gt.1) then
                if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                  call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                   Nplcmax,lcrange)
!                else if(Flagbc.eq.2) then
!                  call ptsmv2perd_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
!                                  Nplcmax,lcrange,Perdlen)
                else
                  call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                   Nplcmax,lcrange)
                endif
              endif
              call setnpt_BeamBunch(Bpts,Nplocal)

            endif
200         continue
! end space charge field calcualtion.
!-------------------------------------------------------------------

!-------------------------------------------------------------------
! use linear map or nonlinear Lorentz integrator to advance particles.
            if(Flagmap.eq.1) then
            ! kick particles in velocity space.
              call map2_BeamBunch(Bpts,tau2,Nxlocal,Nylocal,Nzlocal,&
                   Potential%FieldQ,Ageom,grid2d,Flagbc,Perdlen)
              call map1_BeamBunch(Bpts,Blnelem(i),z,tau1)
            else
              call lostcount_BeamBunch(Bpts,Nplocal,Np,piperad,piperad2)
!!changed into round pipe instead of rectanglar one, 08/21/08, Q.Z.
!!              call lostcount_BeamBunch(Bpts,Nplocal,Np,xradmin,xradmax,&
!!                   yradmin,yradmax)
              !call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
              call map2_BeamBunch(Bpts,Blnelem(i),z,tau2,Nxlocal,Nylocal,&
                           Nzlocal,Potential%FieldQ,Ageom,grid2d,Flagbc,Flagerr)
              if(bitype.ne.3) then
                call map1_BeamBunch(Bpts,z,tau2)
              else
                z = z + 0.5*tau2
                gamma0 = -Bpts%refptcl(6)
                beta0 = sqrt(1.0-1.0/gamma0/gamma0)
                Bpts%refptcl(5) = Bpts%refptcl(5) +0.5*tau2/(Scxl*beta0)
              endif
              !call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
            endif

            if((Flagerr.eq.1).and.(Flagmap.eq.1)) then
              call geomerrT_BeamBunch(Bpts,Blnelem(i)) 
            end if
            if((Flagerr.eq.1).and.(bitype.eq.3)) then
              call geomerrT_BeamBunch(Bpts,Blnelem(i)) 
            end if

            if(Flagdiag.eq.1) then
!              print*,"diag2: ",z,Nplocal
              call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
            elseif(Flagdiag.eq.2) then
              call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
            endif

            nstep = nstep + 1
            !if(myid.eq.0) then 
            !  print*,"j, nstep, z",j,nstep,z
            !endif
!            if(mod(nstep,5).eq.0) then
!              !output all particles in 6d phase space in file "xxnstep"
!              call phaseout_Output(nstep/5,Bpts)
!              !output all geomtry information in file "xxnstep".
!              !call outgeom_Output(nstep,z,i,j,npx,npy,Ageom)
!              !call outpoint_Output(myid+31,Bpts,z,i,j,npx,npy,Ageom)
!            endif
            !output six 2-D phase projections.
            !call phase2d_Output(Bpts,Np)
            !output accumulated density distribution. (r,x,y)
            !call phaserd_Output(20,Bpts,Np)
            !call phaseleda_Output(10,Bpts)
          end do

          else !//bend magnet using Time domain
!            print*,"inside bend0..",bitype
            call getparam_BeamLineElem(Blnelem(i),dparam)
            if(dparam(4).gt.100.0) then !file id > 100 to use 2nd order matrix
              if(Flagerr.eq.1) then
                call geomerrL_BeamBunch(Bpts,Blnelem(i))
              end if
              z = z + blength
              dpi = 2*asin(1.0)
!!20130607  introduce input momentum of "beta-gamma" in the input line              
              gambet00 = dparam(3)
              gamma00 = sqrt(1+gambet00*gambet00)
              beta00 = sqrt(1.0-1.0/gamma00/gamma00)
              gamma0 = -Bpts%refptcl(6)
              beta0 = sqrt(1.0-1.0/gamma0/gamma0)
              Bpts%refptcl(5) = Bpts%refptcl(5)+blength/(Scxl*beta0)
              gambet = beta0*gamma0
              hgap = 2*dparam(5)
              ang0 = dparam(2)*dpi/180.0
!!20130607              hd1 = dparam(3) !k1
              hd1=0.0
!!20130607  this entry is used by "reference energy of beta-gamma"
!!          therefore, hd1 is hardwire into 0, as it should be in most of the normal cases
              angF = dparam(6)*dpi/180.0 !e1
              angB = dparam(7)*dpi/180.0 !e2
              hF = dparam(8) !1/r1
              hB = dparam(9) !/1/r2
              dstr1 = dparam(10) !fringe field K of entrance 
              dstr2 = dstr1 !fringe field K of exit. here, assume Kb = Kf
!              print*,"inside the bend:",dparam,blength

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

!              print*,"qm0: ",qm0,hd0,angF,angB,beta0,gambet,dstr1,&
!                             gamma0,scxl,blength,Nplocal
              r0  = 1.0/hd0

              !print*,"before: ",qm0,Bpts%Pts1(7,1)
              do ipt = 1, Nplocal
                !print*,"gambetz: ",ipt,gambetz
                ptarry(1) = Bpts%Pts1(1,ipt)*Scxl
                gamn = gamma0 - Bpts%Pts1(6,ipt) 
                gambetz = sqrt(gamn**2-1.0-Bpts%Pts1(2,ipt)**2-&
                               Bpts%Pts1(4,ipt)**2)
                ptarry(2) = Bpts%Pts1(2,ipt)/gambetz
                ptarry(3) = Bpts%Pts1(3,ipt)*Scxl
                ptarry(4) = Bpts%Pts1(4,ipt)/gambetz
!!20130607                ptarry(5) = -Bpts%Pts1(5,ipt)*beta0*Scxl
!!20130607                ptarry(6) = -Bpts%Pts1(6,ipt)/beta0/gambet - &
!!20130607                            (Bpts%Pts1(7,ipt)-qm0)/qm0
                ptarry(5) = -Bpts%Pts1(5,ipt)*beta00*Scxl
                ptarry(6) = (-Bpts%Pts1(6,ipt)+(gamma0-gamma00))/beta00/gambet00 - &
                            (Bpts%Pts1(7,ipt)-qm0)/qm0

                qmrel = (Bpts%Pts1(7,ipt)-qm0)/qm0
                !hdi = hd0
                !dstr1i = dstr1
                !print*,"sumpt1: ",sum(ptarry),qm0
                if(dparam(4) .eq. 500.0 .or. dparam(4) .eq. 600.0) then
			ptarry2 = ptarry	!transformed input array for the "sector"
			goto 501 
		  endif
!!skip the entrance fringe field for "500" - the middle piece(s) for a break-down dipole
!!skip the entrance fringe field for "600" - the last piece for a break-down dipole
                call Fpol_Dipole(hd0,hF,tanphiF,tanphiFb,hd1,&
                                 psi1,ptarry,ptarry2,angF)
!!subroutine Fpol_Dipole(h0,h1,tanphi,tanphib,k1,psi,ptarry1,ptarry2,ang)
                !print*,"sumpt2: ",sum(ptarry2),qm0
!                ptarry2 = ptarry
501             continue
!!20130607                call Sector_Dipole(blength,beta0,hd0,hd1,ptarry2,&
!!20130607                                   ptarry,qmrel)
                call Sector_Dipole(blength,beta00,hd0,hd1,ptarry2,&
                                   ptarry,qmrel)
!!subroutine Sector_Dipole(len,beta,h0,k1,ptarry1,ptarry2,qmrel)
                !print*,"sumpt3: ",sum(ptarry),qm0
!                ptarry2 = ptarry
                if(dparam(4) .eq. 500.0 .or. dparam(4) .eq. 400.0) then
		 	ptarry2 = ptarry	!transformed array for skipping "exit edge"
			goto 502 
		  endif
!!skip the entrance fringe field for "500" - the middle piece(s) for a break-down dipole
!!skip the entrance fringe field for "400" - the first piece for a break-down dipole
                call Bpol_Dipole(hd0,hB,tanphiB,tanphiBb,hd1,&
                                 psi2,ptarry,ptarry2,angB)
!!subroutine Bpol_Dipole(h0,h1,tanphi,tanphib,k1,psi,ptarry1,ptarry2,ang)
                !print*,"sumpt4: ",sum(ptarry2),qm0
502             continue                
                Bpts%Pts1(1,ipt) = ptarry2(1)/Scxl 
                Bpts%Pts1(3,ipt) = ptarry2(3)/Scxl 
!!20130607                Bpts%Pts1(5,ipt) = -ptarry2(5)/(Scxl*beta0) 
                Bpts%Pts1(5,ipt) = -ptarry2(5)/(Scxl*beta00)
                gambetz = sqrt((gamn**2-1)/(1+ptarry2(2)**2+ptarry2(4)**2))
                Bpts%Pts1(2,ipt) = ptarry2(2)*gambetz 
                Bpts%Pts1(4,ipt) = ptarry2(4)*gambetz 
              enddo
!              print*,"after: ",qmrel
              if(Flagerr.eq.1) then
                call geomerrT_BeamBunch(Bpts,Blnelem(i))
              end if

!              do ipt = 1, Nplocal
!                ptarry(1) = Bpts%Pts1(1,ipt)*Scxl
!                ptarry(2) = Bpts%Pts1(2,ipt)/gambet*(qm0/Bpts%Pts1(7,ipt))
!                ptarry(3) = Bpts%Pts1(3,ipt)*Scxl
!                ptarry(4) = Bpts%Pts1(4,ipt)/gambet*(qm0/Bpts%Pts1(7,ipt))
!                ptarry(5) = -Bpts%Pts1(5,ipt)*Scxl
!                tmppt = gamma0 - Bpts%Pts1(6,ipt)
!                ptarry(6) = (tmppt*(qm0/Bpts%Pts1(7,ipt))-gamma0)/gambet
!                hdi = hd0
!                dstr1i = dstr1
!                print*,"ptout0: ",ipt,ptarry(1:6)
!                call Fpol_Dipole(hdi,hF,tanphiF,tanphiFb,dstr1i,&
!                                 ptarry,ptarry2)
!                print*,"ptout1: ",ipt,ptarry2(1:6)
!                call Sector_Dipole(blength,beta0,hdi,dstr1i,ptarry2,&
!                                   ptarry)
!                print*,"ptout2: ",ipt,ptarry(1:6)
!                call Bpol_Dipole(hdi,hB,tanphiB,tanphiBb,dstr1i,&
!                                 ptarry,ptarry2)
!                print*,"ptout3: ",ipt,ptarry2(1:6)
!                Bpts%Pts1(1,ipt) = ptarry2(1)/Scxl 
!                Bpts%Pts1(2,ipt) = ptarry2(2)*gambet*(Bpts%Pts1(7,ipt)/qm0)
!                Bpts%Pts1(3,ipt) = ptarry2(3)/Scxl 
!                Bpts%Pts1(4,ipt) = ptarry2(4)*gambet*(Bpts%Pts1(7,ipt)/qm0)
!                Bpts%Pts1(5,ipt) = -ptarry2(5)/Scxl 
!                Bpts%Pts1(6,ipt) = gamma0 - (ptarry2(6)*gambet+gamma0)*&
!                                  (Bpts%Pts1(7,ipt)/qm0)
!              enddo

!              do ipt = 1, Nplocal
!                ptarry(1) = Bpts%Pts1(1,ipt)*Scxl
!                ptarry(2) = Bpts%Pts1(2,ipt)/gambet
!                ptarry(3) = Bpts%Pts1(3,ipt)*Scxl
!                ptarry(4) = Bpts%Pts1(4,ipt)/gambet
!                ptarry(5) = -Bpts%Pts1(5,ipt)*Scxl
!                ptarry(6) = -Bpts%Pts1(6,ipt)/gambet
!                qmi = Bpts%Pts1(7,ipt)
!                hdi = qmi/qm0*hd0
!                ri = 1.0/hdi
!                drr = abs(ri)-abs(r0)
!                rra = -drr*cosang0+sqrt(ri**2+drr**2*cosang0**2-drr**2)
!                if(ang0.ge.0) then
!                  rr0 = rra
!                else
!                  rr0 = -rra
!                endif
!                dstr1i = qmi/qm0*dstr1
!                !deltax = (1.0/hdi - 1.0/hd0)/Scxl
!                deltax = (rr0-r0)/Scxl
!                thetap = asin(abs(rr0/ri)*sinang0)
!                deltat = (abs(ri*thetap)-abs(r0*ang0))/(Scxl*beta0)
!                !deltat = deltax*ang0/beta0 
!                call Fpol_Dipole(hdi,hF,tanphiF,tanphiFb,dstr1i,&
!                                 ptarry,ptarry2)
!                !print*,"ptout1: ",ipt,ptarry2(1:6)
!                call Sector_Dipole(blength,beta0,hdi,dstr1i,ptarry2,&
!                                   ptarry)
!                !print*,"ptout2: ",ipt,ptarry(1:6)
!                call Bpol_Dipole(hdi,hB,tanphiB,tanphiBb,dstr1i,&
!                                 ptarry,ptarry2)
!                !print*,"ptout3: ",ipt,ptarry2(1:6)
!                Bpts%Pts1(1,ipt) = ptarry2(1)/Scxl + deltax
!                Bpts%Pts1(2,ipt) = ptarry2(2)*gambet
!                Bpts%Pts1(3,ipt) = ptarry2(3)/Scxl 
!                Bpts%Pts1(4,ipt) = ptarry2(4)*gambet
!                Bpts%Pts1(5,ipt) = -ptarry2(5)/Scxl + deltat
!                Bpts%Pts1(6,ipt) = -ptarry2(6)*gambet
!                !print*,"deltax, deltat: ",qm0,qmi,ipt,deltax,deltat
!              enddo

            else

            !//go to T frame
            call convZT_BeamBunch(Bpts)
            !//loop through bnseg steps
            gam = sqrt(1.0+Bpts%refptcl(6)**2)
            !//normalized t (omega t)
            tau2 = 2*pi*Scfreq*blength/(Clight*Bpts%refptcl(6)/gam)/bnseg
            tg = Bpts%refptcl(5)
            tv = 0.0
            call diagnosticT_Output(tv,Bpts)
            Flagbctmp = 1
            do j = 1, bnseg
              call drifthalfT_BeamBunch(Bpts,tv,tau2)
              ptref = Bpts%refptcl
              !//go to the local "ptref" coordinates
              call rotto_BeamBunch(Bpts,ptref,ptrange)
              if(Bcurr.lt.1.0e-10) goto 400
              call update_CompDom(Ageom,ptrange,grid2d,Flagbctmp)
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
              call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                Nplcmax,lcrange)
              ! assign new 'Nplocal' local particles on each processor.
              call setnpt_BeamBunch(Bpts,Nplocal)
              call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,&
                                npy)
              deallocate(chgdens)
              allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
              ! deposit particles onto grid to obtain charge density.
              call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,Ageom,&
                                  grid2d,chgdens,Flagbc,Perdlen)
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
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
400           continue
              !//kick particles in the rotated local coordinates.
              call kickT_BeamBunch(Bpts,Blnelem(i),tv,tau2,Nxlocal,Nylocal,&
                      Nzlocal,Potential%FieldQ,Ageom,grid2d,Flagbc,Flagerr)
              call rotback_BeamBunch(Bpts,ptref)
              call drifthalfT_BeamBunch(Bpts,tv,tau2)
              tv = tv + tau2
              call diagnosticT_Output(tv,Bpts)
            enddo
            tg = tg+tv
            call convTZ_BeamBunch(Bpts,tg) !//go back to Z frame
            z = z + blength
     
            endif
          endif  !//end bend magnets
        enddo

! final output.
        call MPI_BARRIER(comm2d,ierr)
        !output six 2-D phase projections.
        !call phase2dold_Output(30,Bpts,Np)
        !output all particles in 6d phase space.
        !call phase_Output(30,Bpts)
!        call phaserd2_Output(20,Bpts,Np,0.0,0.0,0.0,0.0,0.0)
        t_integ = t_integ + elapsedtime_Timer(t0)
        call showtime_Timer()

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
        

        end subroutine run_AccSimulator

        subroutine restart_AccSimulator(iend,nstepend,ibalend,zend)
        implicit none
        include 'mpif.h'
        integer, intent(inout) :: iend,nstepend,ibalend
        double precision, intent(inout) :: zend
        integer :: i,j,bnseg,bmpstp,bitype
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,npx,npy,&
                   totnp,ierr,nylcr,nzlcr
        integer :: nbal,ibal,nstep,ifile
        integer :: tmpfile,nfile
        integer, dimension(3) :: lcgrid
        integer, allocatable, dimension(:) :: lctabnmx,lctabnmy
        integer, allocatable, dimension(:,:,:) :: temptab
        double precision :: z0,z,tau1,tau2,blength,t0
        double precision, allocatable, dimension(:,:) :: lctabrgx, lctabrgy
        double precision, dimension(6) :: lcrange, range, ptrange
        double precision, dimension(3) :: msize
        double precision :: hy,hz,ymin,zmin,piperad,zedge
        double precision :: tmp1,tmp2,tmp3,tmp4,rfile
        double precision, allocatable, dimension(:,:,:) :: chgdens
        double precision, allocatable, dimension(:,:,:) :: besscoef
        double precision, allocatable, dimension(:,:) :: bessnorm,gml
        integer, allocatable, dimension(:) :: modth,pydisp
        integer :: nmod,k
        !double precision :: sumtest, sumtest2, sumtest3
        integer :: rffile,jend
        double precision :: piperad2
  
!-------------------------------------------------------------------
! prepare initial parameters, allocate temporary array.
        rffile = 0
        iend = 0
        jend = 0
        ibalend = 0
        nstepend = 0
        zend = 0.0

        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid2d,totnp,npy,npx)
        if(myid.eq.0) then
          !print*,&
          !"please input rffile, iend, jend, ibalend, nstepend, zend:"
          !read(*,*)rffile,iend,jend,ibalend,nstepend,zend
          print*,&
          "Input rffile, iend, jend, ibalend, nstepend, zend from restart.in"
          open(2,file="restart.in",status="old")
          read(2,*)rffile,iend,jend,ibalend,nstepend,zend
          close(2)
        endif

        call MPI_BCAST(rffile,1,MPI_INTEGER,0,comm2d,ierr)
        call MPI_BCAST(iend,1,MPI_INTEGER,0,comm2d,ierr)
        call MPI_BCAST(jend,1,MPI_INTEGER,0,comm2d,ierr)
        call MPI_BCAST(ibalend,1,MPI_INTEGER,0,comm2d,ierr)
        call MPI_BCAST(nstepend,1,MPI_INTEGER,0,comm2d,ierr)
        call MPI_BCAST(zend,1,MPI_DOUBLE_PRECISION,0,comm2d,ierr)
        call MPI_BARRIER(comm2d,ierr)

        call starttime_Timer(t0)

        allocate(lctabnmx(0:npx-1))
        allocate(lctabnmy(0:npy-1))
        allocate(lctabrgx(2,0:npx-1))
        allocate(lctabrgy(2,0:npy-1))
        allocate(temptab(2,0:npx-1,0:npy-1))

        nbal = 5
        ibal = ibalend
        nstep = nstepend
        z = zend

        if(Flagdiag.eq.1) then
          call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
        else
          call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
        endif

        allocate(chgdens(1,1,1))
!-------------------------------------------------------------------
! prepare for round pipe, open longitudinal
        if(Flagbc.eq.3) then
          allocate(pydisp(0:npy-1))
          call getlctabnm_CompDom(Ageom,temptab)
          pydisp(0) = 0
          do i = 1, npy-1
            pydisp(i) = pydisp(i-1) + temptab(2,0,i-1)
          enddo
          call getlcmnum_CompDom(Ageom,lcgrid)
          Nxlocal = lcgrid(1)
          if(npcol.gt.1) then
            Nylocal = lcgrid(2) + 2
          else
            Nylocal = lcgrid(2)
          endif
          if(nprow.gt.1) then
            Nzlocal = lcgrid(3) + 2
          else
            Nzlocal = lcgrid(3)
          endif

          if(nprow.gt.1) then
            nzlcr = Nzlocal-2
          else
            nzlcr = Nzlocal
          endif
          if(npcol.gt.1) then
            nylcr = Nylocal-2
          else
            nylcr = Nylocal
          endif
          if(myidy.eq.(npcol-1)) nylcr = nylcr - 1
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

!-------------------------------------------------------------------
! start looping through 'Nblem' beam line elements.
          if(myid.eq.0) then
            print*,"enter elment: ",iend
          endif

          call getparam_BeamLineElem(Blnelem(iend),blength,bnseg,bmpstp,&
                                     bitype)
          call getradius_BeamLineElem(Blnelem(iend),piperad,piperad2)
          nfile = 0
          tau1 = 0.0
          if(bitype.ge.0) tau1 = 0.5*blength/bnseg
          tau2 = 2.0*tau1

!-------------------------------------------------------------------
! read in the on axis E field for rf cavities.
          if(bitype.gt.100) then
            call getparam_BeamLineElem(Blnelem(iend),5,rfile)
            nfile = int(rfile + 0.1)
            tmpfile = 0
            ifile = nfile
            if(ifile.ne.tmpfile)then
              !for linear map integrator
              if(Flagmap.eq.1) then
                !read in Ez, Ez', Ez'' on the axis
                call read1_Data(ifile)
              else
                !read in Er, Ez, H_theta on r-z grid 
                !call read2_Data(ifile)
                !read in Fourier coefficients of Ez on the axis
                call read3_Data(ifile)
              endif
              tmpfile=ifile
            endif
          endif

!-------------------------------------------------------------------
! loop through 'bnseg' numerical segments in each beam element
! using 2 step symplectic integeration (ie. leap frog).
          zedge = z
          call setparam_BeamLineElem(Blnelem(iend),1,zedge)
          if(myid.eq.0) print*,"zedge: ",zedge
          do j = jend+1, bnseg
            if((Flagerr.eq.1).and.(Flagmap.eq.1)) then
              call geomerrL_BeamBunch(Bpts,Blnelem(iend)) 
            end if

!-------------------------------------------------------------------
! use linear map or nonlinear Lorentz integrator to advance particles.
            ! spatial drift.
            !linear map integrator
            if(Flagmap.eq.1) then
              call map1_BeamBunch(Bpts,Blnelem(iend),z,tau1)
            else
              call map1_BeamBunch(Bpts,z,tau2)
            endif

!-------------------------------------------------------------------
! escape the space charge calculation for 0 current case
            if(Bcurr.lt.1.0e-30) goto 200

            call conv1st_BeamBunch(Bpts,tau2,Nplocal,Np,ptrange,&
                                   Flagbc,Perdlen,piperad,piperad2)

            !fix the global range for sub-cycle of space charge potential.
            if(Flagsubstep.eq.1) then
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                ptrange(1) = 0.0
                ptrange(2) = piperad
                ptrange(3) = 0.0
                ptrange(4) = 4*asin(1.0)
              else
                ptrange(1) = -piperad
                ptrange(2) = piperad
                ptrange(3) = -piperad2
                ptrange(4) = piperad2
              endif
              ptrange(5) = -Perdlen/2
              ptrange(6) = Perdlen/2
            endif

            ! get new boundary from the range of beam particles.
            if(Flagbc.eq.4) then
            else
              call update_CompDom(Ageom,ptrange,grid2d,Flagbc)
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
            deallocate(chgdens)
            allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
            call getlcrange_CompDom(Ageom,lcrange)
            call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,&
                                npy)

            if(totnp.gt.1) then
              ! pass particles to local space domain on new processor.
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                      Nplcmax,lcrange)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                Nplcmax,lcrange)
              endif
            endif
            ! assign new 'Nplocal' local particles on each processor.
            call setnpt_BeamBunch(Bpts,Nplocal)

            ! deposit particles onto grid to obtain charge density.
            call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,Ageom,&
                                  grid2d,chgdens,Flagbc,Perdlen)

!-------------------------------------------------------------------
! start load balance.
            if((mod(ibal,nbal).eq.0).and.(totnp.gt.1)) then
              call MPI_BARRIER(comm2d,ierr)
              if(myid.eq.0) then 
                print*," load balance! "
              endif

              call getlctabnm_CompDom(Ageom,temptab)
              lctabnmx(0:npx-1) = temptab(1,0:npx-1,0)
              lctabnmy(0:npy-1) = temptab(2,0,0:npy-1)
              call getmsize_CompDom(Ageom,msize)
              hy = msize(2)
              hz = msize(3) 
              call getrange_CompDom(Ageom,range)
              ymin = range(3)
              zmin = range(5)
              if(Flagbc.eq.3) then
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
              call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,npy)

              ! pass particles to local space domain on new processor.
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                      Nplcmax,lcrange)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                 Nplcmax,lcrange)
              endif
              ! assign new 'Nplocal' local particles on each processor.
              call setnpt_BeamBunch(Bpts,Nplocal)

              ! deposit particles onto grid to obtain charge density.
              call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,&
                                    Ageom,grid2d,chgdens,Flagbc,Perdlen)
            endif
            ibal = ibal + 1
!end load balance.
!-------------------------------------------------------------------

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
           
!-------------------------------------------------------------------------
! solve 3D Poisson's equation
            if(Flagbc.eq.1) then
              ! solve Poisson's equation using 3D isolated boundary condition.
              call update3O_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.2) then
              ! solve Poisson's equation using 2D isolated 1D periodic 
              ! boundary condition.
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              call update2O1P_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.3) then
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update3_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,&
              besscoef,bessnorm,gml,modth,nmod)
            else if(Flagbc.eq.4) then
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update4_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,Perdlen)
            else if(Flagbc.eq.5) then
              ! solve Poisson's equation using 2D rectangular pipe, 1D open
              ! boundary condition
              if(myidy.eq.(npy-1)) nylcr = nylcr-1
              call update5_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.6) then
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update6_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else
              print*,"no such boundary condition type!!!"
              stop
            endif

            call cvbkforth1st_BeamBunch(Bpts)
            if(totnp.gt.1) then
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                 Nplcmax,lcrange)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                 Nplcmax,lcrange)
              endif
            endif
            call setnpt_BeamBunch(Bpts,Nplocal)

200         continue
! end space charge field calcualtion.
!-------------------------------------------------------------------

!-------------------------------------------------------------------
! use linear map or nonlinear Lorentz integrator to advance particles.
            if(Flagmap.eq.1) then
            ! kick particles in velocity space.
              call map2_BeamBunch(Bpts,tau2,Nxlocal,Nylocal,Nzlocal,&
                   Potential%FieldQ,Ageom,grid2d,Flagbc,Perdlen)
              call map1_BeamBunch(Bpts,Blnelem(iend),z,tau1)
            else
              call map2_BeamBunch(Bpts,Blnelem(iend),z,tau2,Nxlocal,Nylocal,&
                           Nzlocal,Potential%FieldQ,Ageom,grid2d,Flagbc,Flagerr)
              call map1_BeamBunch(Bpts,z,tau2)
            endif

            if(Flagerr.eq.1) then
              call geomerrT_BeamBunch(Bpts,Blnelem(iend)) 
            end if

            if(Flagdiag.eq.1) then
              call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
            else
              call diagnostic2_Output(Bpts,z,nchrg,nptlist0,outputflag)
            endif

            nstep = nstep + 1
            !if(myid.eq.0) then 
            !  print*,"j, nstep, z",j,nstep,z
            !endif

          end do

          zend = z

        call MPI_BARRIER(comm2d,ierr)

        deallocate(lctabnmx,lctabnmy)
        deallocate(lctabrgx,lctabrgy)
        deallocate(temptab)
        deallocate(chgdens)
        if(Flagbc.eq.3) then
          deallocate(besscoef)
          deallocate(bessnorm)
          deallocate(gml)
          deallocate(pydisp)
          deallocate(modth)
        endif

        end subroutine restart_AccSimulator

        subroutine destruct_AccSimulator(time)
        implicit none
        include 'mpif.h'
        double precision :: time
 
        call destruct_Data()
        call destruct_BeamBunch(Bpts)
        call destruct_FieldQuant(Potential)
        call destruct_CompDom(Ageom)

        deallocate(nptlist0)
        deallocate(currlist0)
        deallocate(qmcclist0)
 
        call end_Output(time)

        end subroutine destruct_AccSimulator

      end module AccSimulatorclass
