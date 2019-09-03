!----------------------------------------------------------------
! (c) Copyright, 2003 by the Regents of the University of California.
! BeamBunchclass: Charged beam bunch class in Beam module of APPLICATION
!                 layer.
! Version: 2.0
! Author: Ji Qiang, Robert Ryne, LBNL, 7/24/03
! Description: This class defines the charged particle beam bunch
!              information in the accelerator.
! Comments: 1) We have added the 3 attributes to the particle array:
!           x,px,y,py,t,pt,charge/mass,charge weight,id. We have moved
!           the charge*curr/freq/Ntot into the charge density calculation,
!           which is represented by the "charge weight" of a particle.
!           2) The map integrator does not work for multiple species, only
!              the Lorenze integrator works for the multiple species.
!----------------------------------------------------------------
      module BeamBunchclass
        use CompDomclass
        use Pgrid2dclass
        use BeamLineElemclass
        use Timerclass
        use Fldmgerclass
        use PhysConstclass
        use Dataclass

        integer, private :: ix_bb, iy_bb, iz_bb, kadd_bb, jadd_bb, flagbc_bb
        double precision, private, allocatable, dimension(:,:,:) :: egx_bb, egy_bb, egz_bb
        double precision, private, dimension(3) :: msize_bb, msizei_bb, sizem_bb

        type BeamBunch
!          private
          !beam freq, current, part. mass and charge.
          double precision :: Current,Mass,Charge
          !# of total global macroparticles and local particles
          integer :: Npt,Nptlocal
          !particles type one.
          double precision, pointer, dimension(:,:) :: Pts1
          !reference particle
          double precision, dimension(6) :: refptcl
        end type BeamBunch
        interface map1_BeamBunch
          module procedure drift1_BeamBunch,drift2_BeamBunch
        end interface
        interface map2_BeamBunch
          module procedure kick1_BeamBunch,kick2_BeamBunch
        end interface
      contains
        subroutine construct_BeamBunch(this,incurr,inkin,inmass,incharge,innp,&
                                       phasini)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: incurr,inkin,inmass,&
                                        incharge,phasini
        integer, intent(in) :: innp
        integer :: myid, myidx, myidy,comm2d,commrow,commcol,ierr
        integer :: nptot,nprocrow,nproccol

        this%Current = incurr
        this%Mass = inmass
        this%Charge = incharge
        this%Npt = innp

        this%refptcl = phasini
        this%refptcl(6) = -(inkin/this%Mass + 1.0)

        end subroutine construct_BeamBunch

        subroutine geomerrL2_BeamBunch(this,beamln,zpos,ferr,fmap)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), intent(in) :: beamln
        double precision, intent(in) :: zpos
        integer, intent(in) :: ferr,fmap
        integer :: i,bnseg,bmapstp,btype
        double precision :: err1,err2,err3,err4,err5,blength
        double precision :: xl,dx,dy,pitch,yaw,roll
        double precision :: gam0,gam,gambetz,cosz,sinz
        double precision, dimension(4) :: ps

        call getparam_BeamLineElem(beamln,blength,bnseg,bmapstp,btype)
        call geterr_BeamLineElem(beamln, err1, err2, err3, err4, err5)

        xl = Scxl

        cosz = cos(roll)
        sinz = sin(roll)

        gam0 = -this%refptcl(6)
        if (fmap.eq.1) then
            gambetz = sqrt(gam0*gam0-1.0)
        endif

        if (ferr.eq.1) then
            pitch = err3
            yaw = err4
            roll = err5
            dx = (err1 - (0.5*blength - zpos)*yaw)/xl
            dy = (err2 - (0.5*blength - zpos)*pitch)/xl
        else if (ferr.eq.2) then
            pitch = (err2 - err4)/blength
            yaw = (err1 - err3)/blength
            roll = err5
            dx = -(err1 - (err1 - err3)*zpos/blength)/xl
            dy = -(err2 - (err2 - err4)*zpos/blength)/xl
        endif

        do i = 1, this%Nptlocal
            ps(1) = this%Pts1(1,i)
            ps(2) = this%Pts1(2,i)
            ps(3) = this%Pts1(3,i)
            ps(4) = this%Pts1(4,i)

            if (fmap.eq.2) then
                gam = gam0 - this%Pts1(6,i)
                gambetz = sqrt(gam*gam - 1.0 - ps(2)*ps(2) - ps(4)*ps(4))
            endif

            if (roll.eq.0.0) then
                this%Pts1(1,i) = ps(1) - dx
                this%Pts1(2,i) = ps(2) - yaw*gambetz
                this%Pts1(3,i) = ps(3) - dy
                this%Pts1(4,i) = ps(4) - pitch*gambetz
            else
                ps(1) = ps(1) - dx
                ps(2) = ps(2) - yaw*gambetz
                ps(3) = ps(3) - dy
                ps(4) = ps(4) - pitch*gambetz

                this%Pts1(1,i) = ps(1)*cosz - ps(3)*sinz
                this%Pts1(2,i) = ps(2)*cosz - ps(4)*sinz
                this%Pts1(3,i) = ps(3)*cosz + ps(1)*sinz
                this%Pts1(4,i) = ps(4)*cosz + ps(2)*sinz
            endif

        enddo

        end subroutine geomerrL2_BeamBunch

        subroutine geomerrT2_BeamBunch(this,beamln,zpos,ferr,fmap)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), intent(in) :: beamln
        double precision, intent(in) :: zpos
        integer, intent(in) :: ferr,fmap
        integer :: i,bnseg,bmapstp,btype
        double precision :: err1,err2,err3,err4,err5,blength
        double precision :: xl,dx,dy,pitch,yaw,roll,pdx,pdy, ppx, ppy, diff
        double precision :: gam0,beta0,gam,gambetz,cosz,sinz,tau,gambetz0,gambetz1
        double precision, dimension(4) :: ps

        call getparam_BeamLineElem(beamln,blength,bnseg,bmapstp,btype)
        call geterr_BeamLineElem(beamln, err1, err2, err3, err4, err5)

        xl = Scxl

        cosz = cos(-roll)
        sinz = sin(-roll)

        gam0 = -this%refptcl(6)
        beta0 = sqrt(1.0-1.0/gam0/gam0)

        if (fmap.eq.1) then
            gambetz = sqrt(gam0*gam0-1.0)
            pdx = 0.0
            pdy = 0.0
            tau = 0.0
        else if (fmap.eq.2) then
            tau = blength/bnseg
        endif

        if (ferr.eq.1) then
            pitch = err3
            yaw = err4
            roll = err5
            dx = (err1 - (0.5*blength - zpos + tau)*yaw)/xl
            dy = (err2 - (0.5*blength - zpos + tau)*pitch)/xl
        else if (ferr.eq.2) then
            pitch = (err2 - err4)/blength
            yaw = (err1 - err3)/blength
            roll = err5
            dx = -(err1 - (err1 - err3)*(zpos-tau)/blength)/xl
            dy = -(err2 - (err2 - err4)*(zpos-tau)/blength)/xl
        endif

        do i = 1, this%Nptlocal
            ps(1) = this%Pts1(1,i)
            ps(2) = this%Pts1(2,i)
            ps(3) = this%Pts1(3,i)
            ps(4) = this%Pts1(4,i)

            if (fmap.eq.2) then
                gam = gam0 - this%Pts1(6,i)
                gambetz1 = sqrt(gam*gam - 1.0 - ps(2)*ps(2) - ps(4)*ps(4))
                ppx = ps(2) + yaw*gambetz1
                ppy = ps(4) + pitch*gambetz1
                gambetz = sqrt(gam*gam - 1.0 - ppx*ppx - ppy*ppy)
                diff = 1.0/gambetz - 1.0/gambetz1
                pdx = (yaw + ps(2)*diff)*tau/xl
                pdy = (pitch + ps(4)*diff)*tau/xl
                this%Pts1(5,i) = this%Pts1(5,i) + tau*(gam0-this%Pts1(6,i))*diff/xl
            endif

            if (roll.eq.0.0) then
                this%Pts1(1,i) = ps(1) + dx + pdx
                this%Pts1(2,i) = ps(2) + yaw*gambetz
                this%Pts1(3,i) = ps(3) + dy + pdy
                this%Pts1(4,i) = ps(4) + pitch*gambetz
            else
                this%Pts1(1,i) = ps(1)*cosz - ps(3)*sinz + dx + pdx
                this%Pts1(2,i) = ps(2)*cosz - ps(4)*sinz + yaw*gambetz
                this%Pts1(3,i) = ps(3)*cosz + ps(1)*sinz + dy + pdy
                this%Pts1(4,i) = ps(4)*cosz + ps(2)*sinz + pitch*gambetz
            endif
        enddo

        end subroutine geomerrT2_BeamBunch

        !shift and rotate the beam at the leading edge of the element
        subroutine geomerrL_BeamBunch(this,beamln)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), intent(in) :: beamln
        double precision  :: xerr,yerr,anglerrx,anglerry,anglerrz
        double precision :: dx,dy,anglex,xl,angley,anglez
        double precision, dimension(6) :: temp,tmp
        integer :: i

        xl = Scxl

        call geterr_BeamLineElem(beamln,xerr,yerr,anglerrx,anglerry,anglerrz)
        dx = xerr/xl
        dy = yerr/xl
        anglex = anglerrx
        angley = anglerry
        anglez = anglerrz

!        print*,"errorL: ",xerr,yerr,anglerrx

        do i = 1, this%Nptlocal
          temp(1) = this%Pts1(1,i) - dx
          temp(2) = this%Pts1(2,i)
          temp(3) = this%Pts1(3,i) - dy
          temp(4) = this%Pts1(4,i)
          tmp(1) = temp(1)*cos(anglez) + temp(3)*sin(anglez)
          tmp(2) = temp(2)*cos(anglez) + temp(4)*sin(anglez)
          tmp(3) = -temp(1)*sin(anglez) + temp(3)*cos(anglez)
          tmp(4) = -temp(2)*sin(anglez) + temp(4)*cos(anglez)
          !5 and 6 corresponds to relative phase and energy. you can not
          !simply transform them and combine with space and momentum.
          tmp(5) = this%Pts1(5,i)
          tmp(6) = this%Pts1(6,i)
          temp(1) = tmp(1)*cos(angley)+tmp(5)*sin(angley)
          temp(2) = tmp(2)*cos(angley)+tmp(6)*sin(angley)
          temp(3) = tmp(3)
          temp(4) = tmp(4)
          temp(5) = -tmp(1)*sin(angley)+tmp(5)*cos(angley)
          temp(6) = -tmp(2)*sin(angley)+tmp(6)*cos(angley)
          this%Pts1(1,i) = temp(1)
          this%Pts1(2,i) = temp(2)
          this%Pts1(3,i) = temp(3)*cos(anglex)+temp(5)*sin(anglex)
          this%Pts1(4,i) = temp(4)*cos(anglex)+temp(6)*sin(anglex)
          this%Pts1(5,i) = -temp(3)*sin(anglex)+temp(5)*cos(anglex)
          this%Pts1(6,i) = -temp(4)*sin(anglex)+temp(6)*cos(anglex)
        enddo

        end subroutine geomerrL_BeamBunch

        !shift and rotate the beam at the tail edge of the element
        subroutine geomerrT_BeamBunch(this,beamln)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), intent(in) :: beamln
        double precision  :: xerr,yerr,anglerrx,anglerry,anglerrz
        double precision :: dx,dy,anglex,xl,angley,anglez
        double precision, dimension(6) :: temp,tmp
        integer :: i

        xl = Scxl

        call geterr_BeamLineElem(beamln,xerr,yerr,anglerrx,anglerry,anglerrz)
        dx = xerr/xl
        dy = yerr/xl
        !rotate back
        anglex = -anglerrx
        angley = -anglerry
        anglez = -anglerrz

!        print*,"errorT: ",xerr,yerr,anglerrx

        do i = 1, this%Nptlocal
          !5 and 6 corresponds to relative phase and energy. you can not
          !simply transform them and combine with space and momentum.
          temp(1) = this%Pts1(1,i)
          temp(2) = this%Pts1(2,i)
          temp(3) = this%Pts1(3,i)*cos(anglex)+this%Pts1(5,i)*sin(anglex)
          temp(4) = this%Pts1(4,i)*cos(anglex)+this%Pts1(6,i)*sin(anglex)
          temp(5) = -this%Pts1(3,i)*sin(anglex)+this%Pts1(5,i)*cos(anglex)
          temp(6) = -this%Pts1(4,i)*sin(anglex)+this%Pts1(6,i)*cos(anglex)
          tmp(1) = temp(1)*cos(angley)+temp(5)*sin(angley)
          tmp(2) = temp(2)*cos(angley)+temp(6)*sin(angley)
          tmp(3) = temp(3)
          tmp(4) = temp(4)
          tmp(5) = -temp(1)*sin(angley)+temp(5)*cos(angley)
          tmp(6) = -temp(2)*sin(angley)+temp(6)*cos(angley)
          temp(1) = tmp(1)*cos(anglez) + tmp(3)*sin(anglez)
          temp(2) = tmp(2)*cos(anglez) + tmp(4)*sin(anglez)
          temp(3) = -tmp(1)*sin(anglez) + tmp(3)*cos(anglez)
          temp(4) = -tmp(2)*sin(anglez) + tmp(4)*cos(anglez)
          temp(5) = tmp(5)
          temp(6) = tmp(6)
          this%Pts1(1,i) = temp(1) + dx
          this%Pts1(2,i) = temp(2)
          this%Pts1(3,i) = temp(3) + dy
          this%Pts1(4,i) = temp(4)
          this%Pts1(5,i) = temp(5)
          this%Pts1(6,i) = temp(6)
        enddo

        end subroutine geomerrT_BeamBunch

        ! Drift beam half step using linear map for external field.
        subroutine drift1_BeamBunch(this,beamln,z,tau)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), intent(inout) :: beamln
        double precision, intent(inout) :: z
        double precision, intent (in) :: tau
        double precision, dimension(6) :: temp
        double precision, dimension(6,6) :: xm
        double precision :: t0
        integer :: i,j,k

        call starttime_Timer(t0)

        call maplinear_BeamLineElem(beamln,z,tau,xm,this%refptcl,this%Charge,&
                       this%Mass)

        do i = 1, this%Nptlocal
          do j = 1, 6
            temp(j) = 0.0
            do k = 1, 6
              temp(j) = temp(j) + this%Pts1(k,i)*xm(j,k)
            enddo
          enddo
          do j = 1,6
            this%Pts1(j,i) = temp(j)
          enddo
        enddo

        z=z+tau
!        do i = 1, 6
!          write(11,*)xm(i,1),xm(i,2),xm(i,3),xm(i,4),xm(i,5),xm(i,6)
!        enddo

        t_map1 = t_map1 + elapsedtime_Timer(t0)

        end subroutine drift1_BeamBunch

! drift half step using nonlinear Lorentz integrator
        subroutine drift2_BeamBunch(this,z,tau)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: z
        double precision, intent (in) :: tau
        double precision :: xl,pzi,xp,yp,zp
        double precision :: t0,beta0,gamma0
        double precision, dimension(9) :: blparam
        integer :: mapstp,itype,bsg,ierr,my_rank
        integer :: i

        call starttime_Timer(t0)

        xl = Scxl
        gamma0 = -this%refptcl(6)
        beta0 = sqrt(1.0-1.0/gamma0/gamma0)

        do i = 1, this%Nptlocal
          xp = this%Pts1(2,i)
          yp = this%Pts1(4,i)
          zp = gamma0 - this%Pts1(6,i)

          pzi = 1.0/xl/sqrt(zp*zp - 1.0 - xp*xp - yp*yp)

          this%Pts1(1,i) = this%Pts1(1,i)+0.5*tau*xp*pzi
          this%Pts1(3,i) = this%Pts1(3,i)+0.5*tau*yp*pzi
          this%Pts1(5,i) = this%Pts1(5,i)+0.5*tau*(zp*pzi-1.0/(beta0*xl))
        enddo

        this%refptcl(5) = this%refptcl(5) + 0.5*tau/(xl*beta0)

        z=z+0.5*tau

        t_map1 = t_map1 + elapsedtime_Timer(t0)

        end subroutine drift2_BeamBunch

        subroutine phasescan_BeamBunch(beamln,fmap,chg,mass,refptcl)
        implicit none
        include 'mpif.h'
        type (BeamLineElem), intent(inout) :: beamln
        integer, intent(in) :: fmap
        double precision, intent(in) :: chg,mass
        double precision, dimension(6), intent(in) :: refptcl
        double precision :: refp0,gam0,p0v1,p0v2,p0v3,ekout,p0off

        refp0 = refptcl(5)
        gam0 = -refptcl(6)

        p0v1 = 0.0
        call findsync_BeamBunch(beamln,fmap,chg,mass,refp0,gam0,p0v1)

        p0v2 = p0v1 + 180.0
        call findsync_BeamBunch(beamln,fmap,chg,mass,refp0,gam0,p0v2)

        p0v3 = (p0v1 + p0v2)/2.0
        call engain_BeamBunch(beamln,fmap,chg,mass,refp0,gam0,p0v3,ekout)

        if (ekout .gt. 0.0) then
            p0off = dmod(p0v1+90.0, dble(360.0))
        else
            p0off = dmod(p0v2+90.0, dble(360.0))
        endif

        if (associated(beamln%pccl)) then
            call setparam_BeamLineElem(beamln,25,p0off)
            call setparam_BeamLineElem(beamln,26,refp0)
            call setparam_BeamLineElem(beamln,27,gam0)
            call setparam_BeamLineElem(beamln,28,dble(chg/mass))
        else if (associated(beamln%pemfld)) then
            call setparam_BeamLineElem(beamln,28,p0off)
            call setparam_BeamLineElem(beamln,29,refp0)
            call setparam_BeamLineElem(beamln,30,gam0)
            call setparam_BeamLineElem(beamln,31,dble(chg/mass))
        endif

        end subroutine phasescan_BeamBunch

        subroutine findsync_BeamBunch(beamln,fmap,chg,mass,refp0,gam0,p0io)
        implicit none
        include 'mpif.h'
        type (BeamLineElem), intent(inout) :: beamln
        integer, intent(in) :: fmap
        double precision, intent(in) :: chg,mass,refp0,gam0
        double precision, intent(inout) :: p0io
        integer :: i, maxitr
        double precision :: eps,h,x0,x1,x2,y0,y1,y2,h0,h1,d0,d1,a,b,c,rad,den,dx


        ! root finding by muller method
        maxitr = 1000
        eps = 1d-10
        h = 0.1

        x2 = p0io
        x1 = p0io + h
        x0 = p0io - h

        call engain_BeamBunch(beamln,fmap,chg,mass,refp0,gam0,x0,y0)
        call engain_BeamBunch(beamln,fmap,chg,mass,refp0,gam0,x1,y1)

        do i = 1,maxitr
            h0 = x1 - x0
            h1 = x2 - x1
            call engain_BeamBunch(beamln,fmap,chg,mass,refp0,gam0,x2,y2)

            d0 = (y1 - y0)/h0
            d1 = (y2 - y1)/h1

            a = (d1 - d0)/(h0 + h1)
            b = a*h1 + d1
            c = y2

            rad = abs(sqrt(cmplx(b*b - 4.0*a*c)))
            if (abs(b + rad) > abs(b - rad)) then
                den = b + rad
            else
                den = b - rad
            endif

            dx = -2.0*c/den
            p0io = x2 + dx
            if (abs(dx) < eps) exit
            x0 = x1
            y0 = y1
            x1 = x2
            y1 = y2
            x2 = p0io
        enddo

        end subroutine findsync_BeamBunch

        subroutine engain_BeamBunch(beamln,fmap,chg,mass,refp0,gam0in,p0v,ekout)
        implicit none
        include 'mpif.h'
        type (BeamLineElem), intent(inout) :: beamln
        integer, intent(in) :: fmap
        double precision, intent(in) :: chg,mass,refp0,p0v,gam0in
        double precision, intent(out) :: ekout

        integer :: i,bnseg,bmapstp,btype,FlagDisc,FlagCart
        double precision :: blength,xl,z,ez0,gam0,gam1,beta0,qmcc,&
                            ek0,ek1,p0tmp,tau,tg
        double precision, dimension(6) :: refptcl0
        double precision, dimension(9,0) :: drays
        double precision, dimension(6,6) :: dummy

        qmcc = chg/mass
        xl = Scxl
        call getparam_BeamLineElem(beamln,blength,bnseg,bmapstp,btype)
        z = 0.0
        call setparam_BeamLineElem(beamln,1,z)
        tau = blength/bnseg

        ek0 = (gam0in-1.0)*mass/1.0e6

        call getparam_BeamLineElem(beamln,4,p0tmp)

        if (associated(beamln%pccl)) then
            call setparam_BeamLineElem(beamln,25,p0v-p0tmp)
        else if (associated(beamln%pemfld)) then
            call setparam_BeamLineElem(beamln,28,p0v-p0tmp)
        endif

        if (fmap.eq.1) then
            refptcl0 = 0.0
            refptcl0(5) = refp0
            refptcl0(6) = -gam0in
            do i = 1, bnseg
                call maplinear_BeamLineElem(beamln,z,tau,dummy,refptcl0,chg,mass)
                z = z + tau
            enddo
            gam1 = -refptcl0(6)

        else if (fmap.eq.2) then
            tg = refp0
            gam0 = gam0in
            do i = 1,bnseg
                beta0 = sqrt(1.0 - 1.0/gam0/gam0)
                tg = tg + 0.5*tau/(xl*beta0)
                z = z + 0.5*tau

                call scatter2tt_BeamBunch(0,0,drays,tg,gam0,chg,mass,tau,z,beamln)

                beta0 = sqrt(1.0 - 1.0/gam0/gam0)
                tg = tg + 0.5*tau/(xl*beta0)
                z = z + 0.5*tau
            enddo

            gam1 = gam0
        endif

        ek1 = (gam1-1.0)*mass/1.0e6

        ekout = ek1 - ek0

        end subroutine engain_BeamBunch

!!	Subroutine lostcount_BeamBunch copied directly from BeamBunch.f90 for 020606
!!	to include round aperture instead of rectanglar one, and no longitudinal losses
!!	08/21/2008 Q.Z.
        !counter the particles get lost outside the xrad and yrad.
        !we have not put the lost through rf bucket yet.
        subroutine lostcount_BeamBunch(this,nplc,nptot,xrad,yrad)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: xrad,yrad
        integer, intent(out) :: nplc,nptot
        integer :: i,ilost,i0,ierr,ntmp5
        double precision :: tmpx,tmpy,xl,rad,rad2,gamma0,tmp5, radtest


        xl = Scxl
        rad = (xrad+yrad)/xl/2
        rad2 = rad*rad

        gamma0 = -this%refptcl(6)

        ilost = 0
        do i0 = 1, this%Nptlocal
            i = i0 - ilost

            tmpx = this%Pts1(1,i0)
            this%Pts1(1,i) = tmpx
            this%Pts1(2,i) = this%Pts1(2,i0)

            tmpy = this%Pts1(3,i0)
            this%Pts1(3,i) = tmpy
            this%Pts1(4,i) = this%Pts1(4,i0)

            this%Pts1(5,i) = this%Pts1(5,i0)
            this%Pts1(6,i) = this%Pts1(6,i0)
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)

            ntmp5 = this%Pts1(5,i0)/Pi
            tmp5 = this%Pts1(5,i0) - ntmp5*Pi
            if(mod(ntmp5,2).eq.0) then
                this%Pts1(5,i) = tmp5
            else
                if(tmp5.gt.0.0) then
                    this%Pts1(5,i) = tmp5 - Pi
                else
                    this%Pts1(5,i) = tmp5 + Pi
                endif
            endif

            radtest = tmpx * tmpx+tmpy*tmpy
            if(radtest.ge.rad2) then
                ilost = ilost + 1
            endif
        enddo

        this%Nptlocal = this%Nptlocal - ilost
        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Npt = nptot

        end subroutine lostcount_BeamBunch

        !counter the particles get lost outside the xrad and yrad.
        !we have not put the lost through rf bucket yet.
        subroutine lostcountXY_BeamBunch(this,nplc,nptot,x1,x2,y1,y2)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: x1,x2,y1,y2
        integer, intent(out) :: nplc,nptot
        integer :: i
        double precision :: tmpx,tmpy,xl
        integer :: ilost,i0,ierr
        double precision :: xapmin,xapmax,yapmin,yapmax

        xl = Scxl
        xapmin = x1/xl
        xapmax = x2/xl
        yapmin = y1/xl
        yapmax = y2/xl

        ilost = 0

        do i0 = 1, this%Nptlocal
            i = i0 - ilost

            this%Pts1(1,i) = this%Pts1(1,i0)
            this%Pts1(2,i) = this%Pts1(2,i0)
            this%Pts1(3,i) = this%Pts1(3,i0)
            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = this%Pts1(5,i0)
            this%Pts1(6,i) = this%Pts1(6,i0)
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)
            tmpx = this%Pts1(1,i0)
            tmpy = this%Pts1(3,i0)

            if(tmpx.le.xapmin) then
                ilost = ilost + 1
            else if(tmpx.ge.xapmax) then
                ilost = ilost + 1
            else if(tmpy.le.yapmin) then
                ilost = ilost + 1
            else if(tmpy.ge.yapmax) then
                ilost = ilost + 1
            endif
        enddo

        this%Nptlocal = this%Nptlocal - ilost
        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Npt = nptot

        end subroutine lostcountXY_BeamBunch

        !Beam Parameter Selector
        subroutine lostcountFP_BeamBunch(this,nplc,nptot,e1,e2,cid,flg)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: e1,e2
        integer, intent(in) :: cid,flg
        integer, intent(out) :: nplc,nptot
        integer :: i
        double precision :: tmpp,xl,xt,refeng,qmc,gam,gambet
        integer :: ilost,i0,ierr
        double precision :: fpmin,fpmax

        xl = Scxl
        xt = Rad2deg
        qmc = this%Mass/1.0e6
        gam = -this%refptcl(6)
        gambet = sqrt(gam*gam - 1.0)

        select case (cid)
            case(1, 3)
                fpmin = e1/xl
                fpmax = e2/xl
            case(2, 4)
                fpmin = e1*gambet
                fpmax = e2*gambet
            case(5)
                fpmin = e1/xt
                fpmax = e2/xt
            case(6)
                if (flg.eq.0) then
                    refeng = 0.0
                else
                    refeng = (gam - 1.0)*qmc
                endif
                fpmin = -(e1 - refeng)/qmc
                fpmax = -(e2 - refeng)/qmc
            case default
                fpmin = e1
                fpmax = e2
        end select

        ilost = 0

        do i0 = 1, this%Nptlocal
            i = i0 - ilost

            this%Pts1(1,i) = this%Pts1(1,i0)
            this%Pts1(2,i) = this%Pts1(2,i0)
            this%Pts1(3,i) = this%Pts1(3,i0)
            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = this%Pts1(5,i0)
            this%Pts1(6,i) = this%Pts1(6,i0)
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)
            tmpp = this%Pts1(cid,i0)

            if(tmpp.le.fpmin) then
                ilost = ilost + 1
            else if(tmpp.ge.fpmax) then
                ilost = ilost + 1
            endif
        enddo

        this%Nptlocal = this%Nptlocal - ilost
        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Npt = nptot

        end subroutine lostcountFP_BeamBunch

        subroutine lostcountRFQ_BeamBunch(this,nplc,nptot,xrad,yrad,rad0)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(out) :: nplc,nptot
        double precision, intent(in) :: xrad, yrad,rad0
        integer :: i, i0, ilost, ntmp5, ierr
        double precision :: xl, xmod, ymod, gamma0, tmp5, tmpx, tmpy, tmpz, thrshld

        xl = Scxl
        xmod = xrad - rad0
        ymod = yrad - rad0
        gamma0 = -this%refptcl(6)

        ilost = 0
        do i0 = 1, this%Nptlocal
            i = i0 - ilost

            this%Pts1(1,i) = this%Pts1(1,i0)
            this%Pts1(2,i) = this%Pts1(2,i0)
            this%Pts1(3,i) = this%Pts1(3,i0)
            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = this%Pts1(5,i0)
            this%Pts1(6,i) = this%Pts1(6,i0)
            this%Pts1(7,i) = this%Pts1(7,i0)
            this%Pts1(8,i) = this%Pts1(8,i0)
            this%Pts1(9,i) = this%Pts1(9,i0)
            tmpx = this%Pts1(1,i0)*xl
            tmpy = this%Pts1(3,i0)*xl

            ntmp5 = this%Pts1(5,i0)/Pi
            tmp5 = this%Pts1(5,i0) - ntmp5*Pi
            if(mod(ntmp5,2).eq.0) then
                this%Pts1(5,i) = tmp5
            else
                if(tmp5.gt.0.0) then
                    this%Pts1(5,i) = tmp5 - Pi
                else
                    this%Pts1(5,i) = tmp5 + Pi
                endif
            endif

            if (abs(tmpx) .le. abs(tmpy)) then
                tmpx = max(abs(tmpx) - xmod, 0.0)
            else
                tmpy = max(abs(tmpy) - ymod, 0.0)
            endif

            thrshld = abs(tmpx*tmpx - tmpy*tmpy) - rad0*rad0
            thrshld = max(abs(tmpx) - 2.0*rad0, thrshld)
            thrshld = max(abs(tmpy) - 2.0*rad0, thrshld)
            thrshld = max(this%Pts1(6,i) - gamma0, thrshld)

            if (thrshld .gt. 0.0) then
                ilost = ilost + 1
            endif

        enddo

        this%Nptlocal = this%Nptlocal - ilost
        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Npt = nptot

        end subroutine lostcountRFQ_BeamBunch

        !//update the total current fraction of each charge state
        !//update total # of ptcl for each charge state
        subroutine chgupdate_BeamBunch(this,nchge,idchgold,qmcclist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nchge
        integer, dimension(:), intent(inout) :: idchgold
        double precision, dimension(:), intent(in) :: qmcclist
        integer :: i, j, ierr
        integer, dimension(nchge) :: idchglc,idchg
        double precision :: dd,eps

        eps = 1.0e-20
        idchglc = 0 !//local # of ptcl of each charge state
        do i = 1, this%Nptlocal
          do j = 1, nchge
            dd = abs((this%Pts1(7,i)-qmcclist(j))/qmcclist(j))
            if(dd.lt.eps) then
              idchglc(j) = idchglc(j) + 1
            endif
          enddo
        enddo

        !//get total # of ptcl for each charge state
        call MPI_ALLREDUCE(idchglc,idchg,nchge,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        !//update the total current fraction of each charge state
        do i = 1, this%Nptlocal
          do j = 1, nchge
            dd = abs((this%Pts1(7,i)-qmcclist(j))/qmcclist(j))
            if(dd.lt.eps) then
              this%Pts1(8,i) = this%Pts1(8,i)*idchg(j)/idchgold(j)
            endif
          enddo
        enddo

        !//update total # of ptcl for each charge state
        idchgold = idchg

        end subroutine chgupdate_BeamBunch

        !from z to t beam frame 1st order transformation.
        subroutine conv1st_BeamBunch(this,tau,nplc,nptot,ptrange,&
                                     Flagbc,perd,xrad,yrad)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: tau,xrad,yrad,perd
        integer, intent(in) :: Flagbc
        double precision, dimension(6), intent(out) :: ptrange
        integer, intent(out) :: nplc,nptot
        integer :: i
        double precision :: xk,xl
        double precision :: pi,gam,bet,bbyk,rcpgammai,betai,rad
        double precision :: twopi,radtest,tmp0,tmpx,tmpy,tmp5,halfperd
        integer :: ilost,i0,ierr,ntmp5

        pi = 2.0*asin(1.0)
        twopi = 2.0*pi
        xl = Scxl
        xk = 1/xl

        gam = -this%refptcl(6)
        bet = sqrt(gam**2-1.0)/gam
        bbyk = bet/xk
        rad = (xrad+yrad)/2

        do i = 1, 3
          ptrange(2*i-1) = 1.0e20
          ptrange(2*i) = -1.0e20
        enddo

        halfperd = 0.5*perd

        if(Flagbc.eq.1) then ! open 3D
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
            rcpgammai = 1.0/(-this%Pts1(6,i0)+gam)
            betai = sqrt(1.0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
                        this%Pts1(4,i0)**2) )
            tmpx = this%Pts1(1,i0)*xl
            this%Pts1(1,i) = (this%Pts1(1,i0)-this%Pts1(2,i0) &
                             *this%Pts1(5,i0)*rcpgammai)*xl

            tmp0 = min(tmpx,this%Pts1(1,i))
            if(ptrange(1)>tmp0) then
              ptrange(1) = tmp0
            endif
            tmp0 = max(tmpx,this%Pts1(1,i))
            if(ptrange(2)<tmp0) then
              ptrange(2) = tmp0
            endif

            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
            this%Pts1(3,i) = (this%Pts1(3,i0)-this%Pts1(4,i0) &
                             *this%Pts1(5,i0)*rcpgammai)*xl

            tmp0 = min(tmpy,this%Pts1(3,i))
            if(ptrange(3)>tmp0) then
              ptrange(3) = tmp0
            endif
            tmp0 = max(tmpy,this%Pts1(3,i))
            if(ptrange(4)<tmp0) then
              ptrange(4) = tmp0
            endif

            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gam*betai*this%Pts1(5,i0)*xl

            if(ptrange(5)>this%Pts1(5,i)) then
              ptrange(5) = this%Pts1(5,i)
            endif
            if(ptrange(6)<this%Pts1(5,i)) then
              ptrange(6) = this%Pts1(5,i)
            endif

            this%Pts1(6,i) = this%Pts1(6,i0)
            radtest = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+ &
                      this%Pts1(3,i)*this%Pts1(3,i))
            if(radtest.ge.rad) then
              ilost = ilost + 1
            endif
          enddo
        else if(Flagbc.eq.2) then ! open 2D, 1D z periodic
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
            rcpgammai = 1.0/(-this%Pts1(6,i0)+gam)
            betai = sqrt(1.0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
                       this%Pts1(4,i0)**2) )

            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gam*betai*this%Pts1(5,i0)*xl

            ntmp5 = this%Pts1(5,i)/halfperd
            tmp5 = this%Pts1(5,i) - ntmp5*halfperd
            this%Pts1(5,i) = tmp5 - mod(ntmp5,2)*halfperd
            tmp5 = -this%Pts1(5,i)/(gam*betai*xl)

            this%Pts1(6,i) = this%Pts1(6,i0)

            tmpx = this%Pts1(1,i0)*xl
! for perd bunch
!            this%Pts1(1,i) = (this%Pts1(1,i0)-this%Pts1(2,i0) &
!                              *this%Pts1(5,i0)*rcpgammai)*xl
            this%Pts1(1,i) = (this%Pts1(1,i0)-this%Pts1(2,i0) &
                              *tmp5*rcpgammai)*xl
!            this%Pts1(1,i) = tmpx

            tmp0 = min(tmpx,this%Pts1(1,i))
            !tmp0 = this%Pts1(1,i)
            if(ptrange(1)>tmp0) then
               ptrange(1) = tmp0
            endif
            tmp0 = max(tmpx,this%Pts1(1,i))
            !tmp0 = this%Pts1(1,i)
            if(ptrange(2)<tmp0) then
               ptrange(2) = tmp0
            endif

            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
! for perd bunch
!            this%Pts1(3,i) = (this%Pts1(3,i0)-this%Pts1(4,i0) &
!                              *this%Pts1(5,i0)*rcpgammai)*xl
            this%Pts1(3,i) = (this%Pts1(3,i0)-this%Pts1(4,i0) &
                              *tmp5*rcpgammai)*xl
!            this%Pts1(3,i) = tmpy

            tmp0 = min(tmpy,this%Pts1(3,i))
            !tmp0 = this%Pts1(3,i)
            if(ptrange(3)>tmp0) then
              ptrange(3) = tmp0
            endif
            tmp0 = max(tmpy,this%Pts1(3,i))
            !tmp0 = this%Pts1(3,i)
            if(ptrange(4)<tmp0) then
               ptrange(4) = tmp0
            endif

! for perd bunch
!            this%Pts1(4,i) = this%Pts1(4,i0)
!            this%Pts1(5,i) = -gam*betai*this%Pts1(5,i0)*xl
!            this%Pts1(6,i) = this%Pts1(6,i0)

            radtest = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+ &
                      this%Pts1(3,i)*this%Pts1(3,i))
            if(radtest.ge.rad) then
              ilost = ilost + 1
            endif
          enddo
          ptrange(5) = -halfperd
          ptrange(6) = halfperd
        else if(Flagbc.eq.3) then !round pipe, 1D z open
          ptrange(1) = 0.0
          ptrange(2) = xrad
          ptrange(3) = 0.0
          ptrange(4) = 4*asin(1.0)
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
            rcpgammai = 1.0/(-this%Pts1(6,i0)+gam)
            betai = sqrt(1.0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
                        this%Pts1(4,i0)**2) )
            tmpx = this%Pts1(1,i0)*xl
            this%Pts1(1,i) = (this%Pts1(1,i0)-this%Pts1(2,i0) &
                             *this%Pts1(5,i0)*rcpgammai)*xl

            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
            this%Pts1(3,i) = (this%Pts1(3,i0)-this%Pts1(4,i0) &
                             *this%Pts1(5,i0)*rcpgammai)*xl

            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gam*betai*this%Pts1(5,i0)*xl

            if(ptrange(5)>this%Pts1(5,i)) then
              ptrange(5) = this%Pts1(5,i)
            endif
            if(ptrange(6)<this%Pts1(5,i)) then
              ptrange(6) = this%Pts1(5,i)
            endif

            this%Pts1(6,i) = this%Pts1(6,i0)
            radtest = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+ &
                      this%Pts1(3,i)*this%Pts1(3,i))
            tmp5 = sqrt(tmpx*tmpx+tmpy*tmpy)
            if((radtest.ge.rad).or.(tmp5.ge.rad)) then
              ilost = ilost + 1
            endif
          enddo
        else if(Flagbc.eq.4) then
          ptrange(1) = 0.0
          ptrange(2) = xrad
          ptrange(3) = 0.0
          ptrange(4) = 4*asin(1.0)
          ptrange(5) = -halfperd
          ptrange(6) = halfperd
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
            rcpgammai = 1.0/(-this%Pts1(6,i0)+gam)
            betai = sqrt(1.0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
                        this%Pts1(4,i0)**2) )

            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gam*betai*this%Pts1(5,i0)*xl
            ntmp5 = this%Pts1(5,i)/halfperd
            tmp5 = this%Pts1(5,i) - ntmp5*halfperd
            this%Pts1(5,i) = tmp5 - mod(ntmp5,2)*halfperd
            this%Pts1(6,i) = this%Pts1(6,i0)

            tmp5 = -this%Pts1(5,i)/(gam*betai*xl)

            tmpx = this%Pts1(1,i0)*xl
!            this%Pts1(1,i) = (this%Pts1(1,i0)-this%Pts1(2,i0) &
!                             *this%Pts1(5,i0)*rcpgammai)*xl
            this%Pts1(1,i) = (this%Pts1(1,i0)-this%Pts1(2,i0) &
                             *tmp5*rcpgammai)*xl
!            this%Pts1(1,i) = tmpx

            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
!            this%Pts1(3,i) = (this%Pts1(3,i0)-this%Pts1(4,i0) &
!                             *this%Pts1(5,i0)*rcpgammai)*xl
            this%Pts1(3,i) = (this%Pts1(3,i0)-this%Pts1(4,i0) &
                             *tmp5*rcpgammai)*xl
!            this%Pts1(3,i) = tmpy

            radtest = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+ &
                      this%Pts1(3,i)*this%Pts1(3,i))
            tmp5 = sqrt(tmpx*tmpx+tmpy*tmpy)
            if((radtest.ge.rad).or.(tmp5.ge.rad)) then
              ilost = ilost + 1
            endif
          enddo
        else if(Flagbc.eq.5) then !2D rectangular finite, 1D z open
          ptrange(1) = -xrad
          ptrange(2) = xrad
          ptrange(3) = -yrad
          ptrange(4) = yrad
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
            rcpgammai = 1.0/(-this%Pts1(6,i0)+gam)
            betai = sqrt(1.0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
                        this%Pts1(4,i0)**2) )
            tmpx = this%Pts1(1,i0)*xl
            this%Pts1(1,i) = (this%Pts1(1,i0)-this%Pts1(2,i0) &
                             *this%Pts1(5,i0)*rcpgammai)*xl

            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
            this%Pts1(3,i) = (this%Pts1(3,i0)-this%Pts1(4,i0) &
                             *this%Pts1(5,i0)*rcpgammai)*xl

            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gam*betai*this%Pts1(5,i0)*xl

            if(ptrange(5)>this%Pts1(5,i)) then
              ptrange(5) = this%Pts1(5,i)
            endif
            if(ptrange(6)<this%Pts1(5,i)) then
              ptrange(6) = this%Pts1(5,i)
            endif

            this%Pts1(6,i) = this%Pts1(6,i0)

            if(tmpx.le.(-xrad)) then
              ilost = ilost + 1
            else if(tmpx.ge.xrad) then
              ilost = ilost + 1
            else if(tmpy.le.(-yrad)) then
              ilost = ilost + 1
            else if(tmpy.gt.yrad) then
              ilost = ilost + 1
            else if(this%Pts1(1,i).le.(-xrad)) then
              ilost = ilost + 1
            else if(this%Pts1(1,i).ge.xrad) then
              ilost = ilost + 1
            else if(this%Pts1(3,i).le.(-yrad)) then
              ilost = ilost + 1
            else if(this%Pts1(3,i).ge.yrad) then
              ilost = ilost + 1
            else
            endif
          enddo
        else if(Flagbc.eq.6) then
          ptrange(1) = -xrad
          ptrange(2) = xrad
          ptrange(3) = -yrad
          ptrange(4) = yrad
          ptrange(5) = -halfperd
          ptrange(6) = halfperd
          ilost = 0
          do i0 = 1, this%Nptlocal
            i = i0 - ilost
            ! The following steps go from z to t frame.
            ! 2) Linear algorithm to transfer from z to t frame.
            rcpgammai = 1.0/(-this%Pts1(6,i0)+gam)
            betai = sqrt(1.0-rcpgammai*rcpgammai*(1+this%Pts1(2,i0)**2+ &
                        this%Pts1(4,i0)**2) )

            this%Pts1(4,i) = this%Pts1(4,i0)
            this%Pts1(5,i) = -gam*betai*this%Pts1(5,i0)*xl
            ntmp5 = this%Pts1(5,i)/halfperd
            tmp5 = this%Pts1(5,i) - ntmp5*halfperd
            this%Pts1(5,i) = tmp5 - mod(ntmp5,2)*halfperd
            this%Pts1(6,i) = this%Pts1(6,i0)

            tmp5 = -this%Pts1(5,i)/(gam*betai*xl)

            tmpx = this%Pts1(1,i0)*xl
!            this%Pts1(1,i) = (this%Pts1(1,i0)-this%Pts1(2,i0) &
!                             *this%Pts1(5,i0)*rcpgammai)*xl
            this%Pts1(1,i) = (this%Pts1(1,i0)-this%Pts1(2,i0) &
                             *tmp5*rcpgammai)*xl

            this%Pts1(2,i) = this%Pts1(2,i0)
            tmpy = this%Pts1(3,i0)*xl
!            this%Pts1(3,i) = (this%Pts1(3,i0)-this%Pts1(4,i0) &
!                             *this%Pts1(5,i0)*rcpgammai)*xl
            this%Pts1(3,i) = (this%Pts1(3,i0)-this%Pts1(4,i0) &
                             *tmp5*rcpgammai)*xl

            if(tmpx.le.(-xrad)) then
              ilost = ilost + 1
            else if(tmpx.ge.xrad) then
              ilost = ilost + 1
            else if(tmpy.le.(-yrad)) then
              ilost = ilost + 1
            else if(tmpy.gt.yrad) then
              ilost = ilost + 1
            else if(this%Pts1(1,i).le.(-xrad)) then
              ilost = ilost + 1
            else if(this%Pts1(1,i).ge.xrad) then
              ilost = ilost + 1
            else if(this%Pts1(3,i).le.(-yrad)) then
              ilost = ilost + 1
            else if(this%Pts1(3,i).ge.yrad) then
              ilost = ilost + 1
            else
            endif
          enddo
        else
          print*,"no such boundary condition!!!"
          stop
        endif

        this%Nptlocal = this%Nptlocal - ilost
        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        this%Current = this%Current*nptot/this%Npt
        this%Npt = nptot

        end subroutine conv1st_BeamBunch

        subroutine cvbkforth1st_BeamBunch(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer :: i
        double precision :: gamma0,xk,xl
        double precision :: rcpgammai,betai,beta

        xl = Scxl
        xk = 1/xl

        gamma0 = -this%refptcl(6)
        beta = sqrt(gamma0*gamma0 - 1.0)/gamma0

        do i = 1, this%Nptlocal
          rcpgammai = 1.0/(-this%Pts1(6,i)+gamma0)
          betai = sqrt(1.0-rcpgammai*rcpgammai*(1+this%Pts1(2,i)**2+ &
                                            this%Pts1(4,i)**2) )
          this%Pts1(5,i) = this%Pts1(5,i)*xk/(-gamma0*betai)
          this%Pts1(1,i) = this%Pts1(1,i)*xk+this%Pts1(2,i) &
                            *this%Pts1(5,i)*rcpgammai
          this%Pts1(3,i) = this%Pts1(3,i)*xk+this%Pts1(4,i) &
                            *this%Pts1(5,i)*rcpgammai
        enddo

        !do i = 1, this%Nptlocal
        !  rcpgammai = 1.0/(-this%Pts1(6,i)+gamma0)
        !  betai = sqrt(1.0-rcpgammai*rcpgammai*(1+this%Pts1(2,i)**2+ &
        !               this%Pts1(4,i)**2) )
        !  this%Pts1(1,i) = this%Pts1(1,i)*xl
        !  this%Pts1(3,i) = this%Pts1(3,i)*xl
        !  this%Pts1(5,i) = -gamma0*betai*this%Pts1(5,i)*xl
        !enddo

        end subroutine cvbkforth1st_BeamBunch

        subroutine setselffld_BB(innx,inny,innz,Flagbc,temppotent,ptsgeom,grid)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz, Flagbc
        double precision,dimension(innx,inny,innz),intent(inout) :: temppotent
        type (CompDom), intent(in) :: ptsgeom
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: drange
        double precision :: hxi, hyi, hzi
        integer :: totnp,npx,npy,myid,myidx,myidy
        integer :: i,j,k,yadd,zadd,innp

        flagbc_bb = Flagbc

        if (.not. allocated(egx_bb)) then
            allocate(egx_bb(innx,inny,innz))
            allocate(egy_bb(innx,inny,innz))
            allocate(egz_bb(innx,inny,innz))
            ix_bb = innx
            iy_bb = inny
            iz_bb = innz
        endif

        if ((innx.ne.ix_bb).or.(inny.ne.iy_bb).or.(innz.ne.iz_bb)) then
            deallocate(egx_bb, egy_bb, egz_bb)
            allocate(egx_bb(innx,inny,innz))
            allocate(egy_bb(innx,inny,innz))
            allocate(egz_bb(innx,inny,innz))
            ix_bb = innx
            iy_bb = inny
            iz_bb = innz
        endif

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        if(npy.gt.1) then
          yadd = 1
        else
          yadd = 0
        endif
        if(npx.gt.1) then
          zadd = 1
        else
          zadd = 0
        endif

        call getcdpara_BB(ptsgeom,npx,npy,myidx,myidy,msize_bb,drange,kadd_bb,jadd_bb)
        msizei_bb = 1.0/msize_bb
        hxi = msizei_bb(1)
        hyi = msizei_bb(2)
        hzi = msizei_bb(3)

        sizem_bb(1) = drange(1)
        sizem_bb(2) = drange(3)
        sizem_bb(3) = drange(5)

        ! Exchange potential information on guard grid cells.
        if((Flagbc.eq.1).or.(Flagbc.eq.5)) then ! 3D open
            if(totnp.ne.1) then
                call guardexch1_Fldmger(temppotent,innx,inny,innz,grid)
            endif
        else if((Flagbc.eq.2).or.(Flagbc.eq.6)) then ! 2D open, 1D periodic
            call guardexch2_Fldmger(temppotent,innx,inny,innz,grid)
        else if(Flagbc.eq.3) then
            call guardexch3_Fldmger(temppotent,innx,inny,innz,grid)
        else if(Flagbc.eq.4) then
            call guardexch4_Fldmger(temppotent,innx,inny,innz,grid)
        else
            print*,"no such boundary condition!!!!"
            stop
        endif

        ! cache optimization here.
        !Ex
        !egx_bb = 0.0
        do k = zadd+1, innz-zadd
            do j = yadd+1, inny-yadd
                egx_bb(1,j,k) = hxi*(temppotent(1,j,k)-temppotent(2,j,k))
                do i = 2, innx-1
                    egx_bb(i,j,k) = 0.5*hxi*(temppotent(i-1,j,k) - temppotent(i+1,j,k))
                enddo
            egx_bb(innx,j,k) = hxi*(temppotent(innx-1,j,k)- temppotent(innx,j,k))
            enddo
        enddo

        !Ey
        !egy_bb = 0.0
        if((Flagbc.eq.3).or.(Flagbc.eq.4)) then  ! periodic in theta

            if(npx.gt.1) then
                do k = 1, innz
                    do j = 2, inny-1
                        do i = 2, innx
                            egy_bb(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k) -  &
                                temppotent(i,j+1,k))*hxi/(i-1)
                        enddo
                        egy_bb(1,j,k) = 0.0
                    enddo
                enddo
            else
                do k = 1, innz
                    do j = 2, inny-1
                        do i = 2, innx
                            egy_bb(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k) -  &
                                temppotent(i,j+1,k))*hxi/(i-1)
                        enddo
                        egy_bb(1,j,k) = 0.0
                    enddo
                enddo
                do k = 1, innz
                    do i = 2, innx
                        egy_bb(i,1,k) = 0.5*hyi*(temppotent(i,inny-1,k) -  &
                           temppotent(i,2,k))*hxi/(i-1)
                        egy_bb(i,inny,k) = egy_bb(i,1,k)
                    enddo
                    egy_bb(1,1,k) = 0.0
                    egy_bb(1,inny,k) = 0.0
                enddo
            endif

        else  ! open in Y

            if(npx.gt.1) then
                if(myidy.eq.0) then
                    do k = zadd+1, innz-zadd
                        do i = 1, innx
                            egy_bb(i,yadd+1,k) = hyi*(temppotent(i,yadd+1,k)- &
                                temppotent(i,yadd+2,k))
                        enddo
                    enddo
                    do k = zadd+1, innz-zadd
                        do j = yadd+2, inny-yadd
                            do i = 1, innx
                                egy_bb(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                                    temppotent(i,j+1,k))
                            enddo
                        enddo
                    enddo
                else if(myidy.eq.(npx-1)) then
                    do k = zadd+1, innz-zadd
                        do i = 1, innx
                            egy_bb(i,inny-yadd,k) = hyi*(temppotent(i,inny-yadd-1,k)- &
                                temppotent(i,inny-yadd,k))
                        enddo
                    enddo
                    do k = zadd+1, innz-zadd
                        do j = yadd+1, inny-yadd-1
                            do i = 1, innx
                                egy_bb(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                                    temppotent(i,j+1,k))
                            enddo
                        enddo
                    enddo
                else
                    do k = zadd+1, innz-zadd
                        do j = yadd+1, inny-yadd
                            do i = 1, innx
                                egy_bb(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                                    temppotent(i,j+1,k))
                            enddo
                        enddo
                    enddo
                endif
            else
                do k = zadd+1, innz-zadd
                    do i = 1, innx
                        egy_bb(i,1,k) = hyi*(temppotent(i,1,k)-temppotent(i,2,k))
                    enddo
                    do j = 2, inny-1
                        do i = 1, innx
                            egy_bb(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                                temppotent(i,j+1,k))
                        enddo
                    enddo
                    do i = 1, innx
                        egy_bb(i,inny,k) = hyi*(temppotent(i,inny-1,k)- &
                            temppotent(i,inny,k))
                    enddo
                enddo
            endif
        endif

        !Ez
        !egz_bb = 0.0
        if(npy.gt.1) then
            if((Flagbc.eq.1).or.(Flagbc.eq.3).or.(Flagbc.eq.5)) then ! 3D open
                if(myidx.eq.0) then
                    do j = yadd+1, inny-yadd
                        do i = 1, innx
                            egz_bb(i,j,zadd+1) = hzi*(temppotent(i,j,zadd+1)- &
                                temppotent(i,j,zadd+2))
                        enddo
                    enddo
                    do k = zadd+2, innz-zadd
                        do j = yadd+1, inny-yadd
                            do i = 1, innx
                                egz_bb(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                    temppotent(i,j,k+1))
                            enddo
                        enddo
                    enddo
                else if(myidx.eq.(npy-1)) then
                    do k = zadd+1, innz-zadd-1
                        do j = yadd+1, inny-yadd
                            do i = 1, innx
                                egz_bb(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                    temppotent(i,j,k+1))
                            enddo
                        enddo
                    enddo
                    do j = yadd+1, inny-yadd
                        do i = 1, innx
                            egz_bb(i,j,innz-zadd) = hzi*(temppotent(i,j,innz-zadd-1)- &
                                temppotent(i,j,innz-zadd))
                        enddo
                    enddo
                else
                    do k = zadd+1, innz-zadd
                        do j = yadd+1, inny-yadd
                            do i = 1, innx
                                egz_bb(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                    temppotent(i,j,k+1))
                            enddo
                        enddo
                    enddo
                endif

            else if((Flagbc.eq.2).or.(Flagbc.eq.4).or.(Flagbc.eq.6)) then
            ! 2D open, 1D periodic
                do k = zadd+1, innz-zadd
                    do j = yadd+1, inny-yadd
                        do i = 1, innx
                            egz_bb(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                temppotent(i,j,k+1))
                        enddo
                    enddo
                enddo
            else
                print*,"no such boundary condition!!!"
                stop
            endif

        else
            if((Flagbc.eq.1).or.(Flagbc.eq.3).or.(Flagbc.eq.5)) then ! 3D open
                do j = yadd+1, inny-yadd
                    do i = 1, innx
                        egz_bb(i,j,1) = hzi*(temppotent(i,j,1)-temppotent(i,j,2))
                    enddo
                enddo
            else if((Flagbc.eq.2).or.(Flagbc.eq.4).or.(Flagbc.eq.6)) then
            ! 2D open, 1D periodic
                do j = yadd+1, inny-yadd
                    do i = 1, innx
                        egz_bb(i,j,1)=0.5*hzi*(temppotent(i,j,innz-1)-temppotent(i,j,2))
                    enddo
                enddo
            else
            endif

            do k = 2, innz-1
                do j = yadd+1, inny-yadd
                    do i = 1, innx
                        egz_bb(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                            temppotent(i,j,k+1))
                    enddo
                enddo
            enddo

            if((Flagbc.eq.1).or.(Flagbc.eq.3).or.(Flagbc.eq.5)) then
                do j = yadd+1, inny-yadd
                    do i = 1, innx
                        egz_bb(i,j,innz) = hzi*(temppotent(i,j,innz-1)- &
                            temppotent(i,j,innz))
                    enddo
                enddo
            else if((Flagbc.eq.2).or.(Flagbc.eq.4).or.(Flagbc.eq.6)) then
                do j = yadd+1, inny-yadd
                    do i = 1, innx
                        egz_bb(i,j,innz) = 0.5*hzi*(temppotent(i,j,innz-1)- &
                            temppotent(i,j,2))
                    enddo
                enddo
            else
            endif
        endif

        !Send the E field to the neibhoring guard grid to do the CIC
        !interpolation.
        if(totnp.ne.1) then
          call boundint4_Fldmger(egx_bb,egy_bb,egz_bb,innx,inny,innz,grid)
        endif

        end subroutine setselffld_BB

        subroutine getcdpara_BB(ptsgeom,npx,npy,myidx,myidy,msize,drange,kadd,jadd)
        implicit none
        type (CompDom) :: ptsgeom
        integer, intent(in) :: npx, npy, myidx, myidy
        double precision, intent(out), dimension(3) :: msize
        double precision, intent(out), dimension(6) :: drange
        integer, intent(out) :: kadd, jadd
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer :: i

        call getmsize_CompDom(ptsgeom,msize)
        call getrange_CompDom(ptsgeom,drange)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        end subroutine getcdpara_BB

        subroutine getselffld_BB(loc,gam,efld)
        implicit none
        double precision, intent(in), dimension(3) :: loc
        double precision, intent(in) :: gam
        double precision, intent(out), dimension(6) :: efld
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin,ri,&
                            thi,xcnv,ycnv,zcnv,ab,cd,ef,&
                            fac1,fac2,exn,eyn,ezn
        integer :: ix,jx,kx,ix1,jx1,kx1

        hxi = msizei_bb(1)
        hyi = msizei_bb(2)
        hzi = msizei_bb(3)
        hx = msize_bb(1)
        hy = msize_bb(2)
        hz = msize_bb(3)

        xmin = sizem_bb(1)
        ymin = sizem_bb(2)
        zmin = sizem_bb(3)

        xcnv = loc(1)
        ycnv = loc(2)
        zcnv = loc(3)

        ri = sqrt(xcnv*xcnv+ycnv*ycnv)
        if(xcnv.gt.0.0) then
            if(ycnv.gt.0.0) then
                thi = asin(ycnv/ri)
            else
                thi = 2*Pi+asin(ycnv/ri)
            endif
        else
            thi = Pi - asin(ycnv/ri)
        endif
            if ((flagbc_bb.eq.3) .or. (flagbc_bb.eq.4)) then
                ix=(ri-xmin)*hxi + 1
                ab=(ix*ix*hx*hx-(ri-xmin)*(ri-xmin))/(ix*ix*hx*hx-(ix-1)*(ix-1)*hx*hx)
                jx=(thi-ymin)*hyi + 1 + jadd_bb
                cd=((ymin-thi)+(jx-jadd_bb)*hy)*hyi
                kx=(zcnv-zmin)*hzi + 1 + kadd_bb
                ef=((zmin-zcnv)+(kx-kadd_bb)*hz)*hzi
            else
                ix=(xcnv-xmin)*hxi + 1
                ab=((xmin-xcnv)+ix*hx)*hxi
                jx=(ycnv-ymin)*hyi + 1 + jadd_bb
                cd=((ymin-ycnv)+(jx-jadd_bb)*hy)*hyi
                kx=(zcnv-zmin)*hzi + 1 + kadd_bb
                ef=((zmin-zcnv)+(kx-kadd_bb)*hz)*hzi
            endif

            ix1 = ix + 1
            jx1 = jx + 1
            kx1 = kx + 1

            exn = (egx_bb(ix,jx,kx)*ab*cd*ef  &
                +egx_bb(ix,jx1,kx)*ab*(1.-cd)*ef &
                +egx_bb(ix,jx1,kx1)*ab*(1.-cd)*(1.0-ef) &
                +egx_bb(ix,jx,kx1)*ab*cd*(1.0-ef) &
                +egx_bb(ix1,jx,kx1)*(1.0-ab)*cd*(1.0-ef) &
                +egx_bb(ix1,jx1,kx1)*(1.0-ab)*(1.0-cd)*(1.0-ef)&
                +egx_bb(ix1,jx1,kx)*(1.0-ab)*(1.0-cd)*ef &
                +egx_bb(ix1,jx,kx)*(1.0-ab)*cd*ef)*gam

            eyn = (egy_bb(ix,jx,kx)*ab*cd*ef  &
                +egy_bb(ix,jx1,kx)*ab*(1.-cd)*ef &
                +egy_bb(ix,jx1,kx1)*ab*(1.-cd)*(1.0-ef) &
                +egy_bb(ix,jx,kx1)*ab*cd*(1.0-ef) &
                +egy_bb(ix1,jx,kx1)*(1.0-ab)*cd*(1.0-ef) &
                +egy_bb(ix1,jx1,kx1)*(1.0-ab)*(1.0-cd)*(1.0-ef)&
                +egy_bb(ix1,jx1,kx)*(1.0-ab)*(1.0-cd)*ef &
                +egy_bb(ix1,jx,kx)*(1.0-ab)*cd*ef)*gam

            ezn = egz_bb(ix,jx,kx)*ab*cd*ef  &
                +egz_bb(ix,jx1,kx)*ab*(1.-cd)*ef &
                +egz_bb(ix,jx1,kx1)*ab*(1.-cd)*(1.0-ef) &
                +egz_bb(ix,jx,kx1)*ab*cd*(1.0-ef) &
                +egz_bb(ix1,jx,kx1)*(1.0-ab)*cd*(1.0-ef) &
                +egz_bb(ix1,jx1,kx1)*(1.0-ab)*(1.0-cd)*(1.0-ef)&
                +egz_bb(ix1,jx1,kx)*(1.0-ab)*(1.0-cd)*ef &
                +egz_bb(ix1,jx,kx)*(1.0-ab)*cd*ef

            if ((flagbc_bb.eq.3) .or. (flagbc_bb.eq.4)) then
                fac1 = xcnv/ri
                fac2 = ycnv/ri
                efld(1) = exn*fac1-eyn*fac2
                efld(2) = exn*fac2+eyn*fac1
                efld(3) = ezn
            else
                efld(1) = exn
                efld(2) = eyn
                efld(3) = ezn
            endif

            efld(4) = efld(2)
            efld(5) = efld(1)
            efld(6) = 0

        end subroutine getselffld_BB

        ! Here, all indices of potential are local to processor.
        ! Advance the particles in the velocity space using the force
        ! from the external field and the self space charge force
        ! interpolated from the grid to particles. (linear map)
        subroutine kick1_BeamBunch(this,tau,innx,inny,innz,temppotent,&
                              ptsgeom,grid,Flagbc,perdlen)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: innx, inny, innz, Flagbc
        type (CompDom), intent(in) :: ptsgeom
        double precision, dimension(innx,inny,innz), intent(inout) :: temppotent
        type (Pgrid2d), intent(in) :: grid
        double precision, intent(in) :: tau,perdlen
        double precision :: t0,gam,curr,mass,chrg
        integer :: innp

        call starttime_Timer(t0)

        curr = this%Current

        if(curr.ge.1.0e-30) then
            call setselffld_BB(innx,inny,innz,Flagbc,temppotent,ptsgeom,grid)

            gam = -this%refptcl(6)
            mass = this%Mass
            innp = this%Nptlocal
            chrg = this%Charge
            call scatter1tt_BeamBunch(innp,this%Pts1,gam,tau,mass,chrg)
        endif

        t_force = t_force + elapsedtime_Timer(t0)

        end subroutine kick1_BeamBunch

        subroutine scatter1tt_BeamBunch(innp,rays,gam,tau,mass,chrg)
        implicit none
        include 'mpif.h'
        double precision, intent (inout), dimension (9,innp) :: rays
        double precision, intent (in) :: gam,tau,mass,chrg
        integer, intent(in) :: innp
        integer :: n
        double precision :: twopi,fpei,xk,bet,bbyk,gambet,brho,vz0,&
            perv0,xycon,tcon,rcpgammai,betai,tmpsq,tmpscale
        double precision, dimension(3) :: loc
        double precision, dimension(6) :: efld

        tmpscale = 1.0/Scfreq

        twopi = 2.0*Pi
        fpei = Clight*Clight*1.0e-7
        xk = 1/Scxl

        bet = sqrt(gam**2-1.0)/gam
        bbyk = bet/xk
        gambet = gam*bet
        brho = gambet/Clight
        vz0 = bet*Clight
        perv0 = 2.0*fpei/(brho*vz0*vz0*gam*gam)
        xycon = 0.5*perv0*gambet*bbyk*twopi/tmpscale
        tcon = bet*xycon*gam**2

        do n = 1, innp
            rcpgammai = 1.0/(-rays(6,n)+gam)
            betai = sqrt(1.0-rcpgammai*rcpgammai*(1+rays(2,n)**2+rays(4,n)**2))
            loc(1) = rays(1,n)*Scxl
            loc(2) = rays(3,n)*Scxl
            loc(3) = -gam*betai*rays(5,n)*Scxl
            call getselffld_BB(loc,gam,efld)

            rays(2,n) = rays(2,n)+tau*xycon*efld(1)*rays(7,n)
            rays(4,n) = rays(4,n)+tau*xycon*efld(2)*rays(7,n)
            rays(6,n) = rays(6,n)-tau*tcon*efld(3)*rays(7,n)
        enddo

        end subroutine scatter1tt_BeamBunch

        subroutine kick2_BeamBunch(this,beamelem,z,tau,innx,inny,innz, &
                             temppotent,ptsgeom,grid,Flagbc,flagerr)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), intent(in) :: beamelem
        integer, intent(in) :: innx, inny, innz, Flagbc,flagerr
        double precision,dimension(innx,inny,innz),intent(inout) :: temppotent
        type (CompDom), intent(in) :: ptsgeom
        type (Pgrid2d), intent(in) :: grid
        double precision, intent(in) :: z, tau
        double precision :: t0,tg,chge,mass,curr,gam
        integer :: innp,scflag

        call starttime_Timer(t0)

        curr = this%Current
        innp = this%Nptlocal
        tg = this%refptcl(5)
        gam = -this%refptcl(6)
        mass = this%Mass
        chge = this%Charge

        scflag = 0
        if(curr.ge.1.0e-30) then
            scflag = 1
            call setselffld_BB(innx,inny,innz,Flagbc,temppotent,ptsgeom,grid)
        endif

        call scatter2tt_BeamBunch(innp,scflag,this%Pts1,tg,gam,chge,mass,tau,z,beamelem)

        this%refptcl(6) = -gam

        t_force = t_force + elapsedtime_Timer(t0)

        end subroutine kick2_BeamBunch

        pure subroutine getfc1_BB(escale,ww,tg,theta0,hr,loc,extfld6,extfld)
        implicit none
        double precision, intent(in) :: escale, ww, tg, theta0, hr
        double precision, dimension(3), intent(in) :: loc
        double precision, dimension(6,NrIntvRf+1), intent(in) :: extfld6
        double precision, dimension(6), intent(inout) :: extfld

        integer :: ir,ir1
        double precision :: et,bt,rr0,rr0i,rr,efr,er,btheta

        et = escale*cos(ww*(tg+loc(3))+theta0)
        bt = escale*sin(ww*(tg+loc(3))+theta0)
        rr0 = sqrt(loc(1)*loc(1)+loc(2)*loc(2))
        rr0i = 1.0/rr0
        rr = rr0/hr
        ir = rr + 1
        if(ir.eq.NrIntvRf) ir=ir-1
        ir1 = ir+1
        efr = ir - rr
        er = (extfld6(1,ir)*efr+extfld6(1,ir1)*(1.0-efr))*et
        extfld(1) = extfld(1) + er*loc(1)*rr0i
        extfld(2) = extfld(2) + er*loc(2)*rr0i
        extfld(3) = extfld(3) + (extfld6(3,ir)*efr+extfld6(3,ir1)*(1.0-efr))*et
        btheta = (extfld6(5,ir)*efr+extfld6(5,ir1)*(1.0-efr))*bt
        extfld(4) = extfld(4) - btheta*loc(2)*rr0i
        extfld(5) = extfld(5) + btheta*loc(1)*rr0i
        extfld(6) = extfld(6)

        end subroutine getfc1_BB

        pure subroutine getfc2_BB(escale,ww,tg,theta0,hxx,hyy,loc,extfld6xyz,extfld)
        implicit none
        double precision, intent(in) :: escale,ww,tg,theta0,hxx,hyy
        double precision, dimension(3), intent(in) :: loc
        double precision,dimension(6,NxIntvRfg+1,NyIntvRfg+1), intent(in) :: extfld6xyz
        double precision, dimension(6), intent(inout) :: extfld

        double precision :: xx,yy,efx,efy
        double precision, dimension(6) :: scl,ef00,ef10,ef01,ef11
        integer:: ix,ix1,iy,iy1

        if (ww.ne. 0.0) then
            scl(1:3) = escale*cos(ww*(tg+loc(3))+theta0)
            scl(4:6) = escale*sin(ww*(tg+loc(3))+theta0)
        else
            scl = escale
        endif

        xx = (loc(1) - XminRfg)/hxx
        yy = (loc(2) - YminRfg)/hyy

        ix = int(xx) + 1
        ix1 = ix + 1
        efx = -xx + ix
        iy = int(yy) + 1
        iy1 = iy + 1
        efy = -yy + iy

        if (ix.ge.1 .and. iy.ge.1 .and. ix1.le.NxIntvRfg+1 .and. iy1.le.NyIntvRfg+1) then
            ef00 = extfld6xyz(:,ix,iy)*efx*efy
            ef10 = extfld6xyz(:,ix1,iy)*(1.0-efx)*efy
            ef01 = extfld6xyz(:,ix,iy1)*efx*(1.0-efy)
            ef11 = extfld6xyz(:,ix1,iy1)*(1.0-efx)*(1.0-efy)
            extfld = extfld + (ef00 + ef01 + ef10 + ef11)*scl
        else
            extfld = extfld
        endif

        end subroutine getfc2_BB

        pure subroutine updatep_BB(ray,extfldt,ez0qmcct,gam,tgam,tgam2)
        implicit none
        double precision, dimension(9), intent(inout) :: ray
        double precision, dimension(6), intent(in) :: extfldt
        double precision, intent(in) :: ez0qmcct,gam,tgam,tgam2

        integer :: ii
        double precision :: tmp12,tmp13,tmp23,tmp1,tmp2,tmp3,tmp4,tmp5,&
                            a1,a2,a3,b1,b2,b3,p1,p2,p3,&
                            s11,s12,s13,s21,s22,s23,s31,s32,s33,&
                            p1t,p2t,p3t,pz,pz2,det,tmpsq2

        tmp12 = -Clight*extfldt(6)*ray(7)*0.5
        tmp13 = extfldt(1)*ray(7)*0.5
        tmp23 = extfldt(2)*ray(7)*0.5

        tmp1 = tgam*extfldt(1)*ray(7)
        tmp2 = Clight*extfldt(5)*ray(7)
        tmp3 = tgam*extfldt(2)*ray(7)
        tmp4 = Clight*extfldt(4)*ray(7)
        tmp5 = (ez0qmcct-extfldt(3)*ray(7))

        p1 = ray(2)
        p2 = ray(4)
        p3 = ray(6)
        pz2 = (gam-p3)*(gam-p3)-1.0-p1*p1-p2*p2
        pz = sqrt(pz2)

        a1 = tmp12/pz
        a2 = tmp13/pz
        a3 = tmp23/pz
        b1 = p1-a1*p2-a2*p3+tmp1/pz - tmp2
        b2 = a1*p1+p2-a3*p3+tmp3/pz + tmp4
        b3 = -a2*p1-a3*p2+p3 + tmp5
        det = 1.0+a1*a1-a2*a2-a3*a3
        s11 = 1.0-a3*a3
        s12 = -a1+a2*a3
        s13 = a1*a3-a2
        s21 = a1+a2*a3
        s22 = 1.0-a2*a2
        s23 = -a3-a1*a2
        s31 = -a1*a3-a2
        s32 = -a3+a1*a2
        s33 = 1.0+a1*a1
        p1t = (s11*b1+s12*b2+s13*b3)/det
        p2t = (s21*b1+s22*b2+s23*b3)/det
        p3t = (s31*b1+s32*b2+s33*b3)/det
        tmpsq2 = (tgam2-p3t)*(tgam2-p3t)-1.0-p1t*p1t-p2t*p2t
        if(tmpsq2.gt.0.0) then
            if (pz2.ne.tmpsq2) then
                pz = 0.5*pz + 0.5*sqrt(tmpsq2)
                a1 = tmp12/pz
                a2 = tmp13/pz
                a3 = tmp23/pz
                b1 = p1-a1*p2-a2*p3+tmp1/pz - tmp2
                b2 = a1*p1+p2-a3*p3+tmp3/pz + tmp4
                b3 = -a2*p1-a3*p2+p3 + tmp5
                det = 1.0+a1*a1-a2*a2-a3*a3
                s11 = 1.0-a3*a3
                s12 = -a1+a2*a3
                s13 = a1*a3-a2
                s21 = a1+a2*a3
                s22 = 1.0-a2*a2
                s23 = -a3-a1*a2
                s31 = -a1*a3-a2
                s32 = -a3+a1*a2
                s33 = 1.0+a1*a1
                p1t = (s11*b1+s12*b2+s13*b3)/det
                p2t = (s21*b1+s22*b2+s23*b3)/det
                p3t = (s31*b1+s32*b2+s33*b3)/det
            endif
        endif

        ray(2) = p1t
        ray(4) = p2t
        ray(6) = p3t


        end subroutine updatep_BB

        subroutine scatter2tt_BeamBunch(innp,scflag,rays,tg,gam,chge,mass,tau,z,beamelem)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,scflag
        double precision, intent (inout), dimension (9,innp) :: rays
        double precision, intent (in) :: tau,mass,chge,tg,z
        double precision, intent (inout) :: gam
        type (BeamLineElem), intent(in) :: beamelem
        integer :: n, FlagCart, FlagDisc, csid
        double precision, dimension(4) :: pos
        double precision, dimension(3) :: loc
        double precision, dimension(6) :: extfld,selffld,sclfld,selfscale
        double precision, dimension(6,innp) :: textfld
        double precision, dimension(6,NrIntvRf+1) :: extfld6
        double precision :: zedge,t0,tgam,tgam2,qmcc,ez0,hxx,hyy,hr,&
                            beta0,et,bt,ezqmct,gamn,&
                            escl,theta,ww,&
                            escl2,theta2,ww2,&
                            escl3,theta3,ww3,&
                            xref,xoff,yoff,rcpgammai,betai,xl
        double precision,dimension(6,NxIntvRfg+1,NyIntvRfg+1) :: &
                            extfld6xyz,extfld6xyz2,extfld6xyz3

        call starttime_Timer( t0 )

        if (scflag.eq.1) then
            selfscale(1:3) = Clight*Clight*1.0e-7
            beta0 = sqrt(gam**2-1.0)/gam
            selfscale(4) = -beta0*Clight*1.0e-7
            selfscale(5) = -selfscale(4)
            selfscale(6) = 0.0
        else
            selfscale = 0.0
        endif

        qmcc = chge/mass

        csid = 0

        FlagDisc = 0
        FlagCart = 1
        if(associated(beamelem%pemfld)) then
            if (beamelem%pemfld%Itype .eq. 110) then
                csid = 1
                zedge = beamelem%pemfld%Param(1)

                escl = beamelem%pemfld%Param(2)
                ww = beamelem%pemfld%Param(3)/Scfreq
                theta = (beamelem%pemfld%Param(4)+beamelem%pemfld%Param(28))*Pi/180.0

                escl2 = beamelem%pemfld%Param(18)
                ww2 = beamelem%pemfld%Param(19)/Scfreq
                theta2 = theta+(beamelem%pemfld%Param(20))*Pi/180.0

                escl3 = beamelem%pemfld%Param(21)
                ww3 = beamelem%pemfld%Param(22)/Scfreq
                theta3 = theta+(beamelem%pemfld%Param(23))*Pi/180.0

                FlagDisc = 0.1+beamelem%pemfld%Param(13)
                FlagCart = 0.1+beamelem%pemfld%Param(14)

                if(.not. FlagDisc.eq.1) then
                    pos(1) = 0.0
                    pos(2) = 0.0
                    pos(3) = z
                    pos(4) = tg
                    call getfld_BeamLineElem(beamelem,pos,extfld)
                    ez0 = extfld(3)
                else
                    ez0 = 0.0
                endif

                if(FlagCart.eq.1) then
                    call getfld6_EMfld(z-zedge,extfld6)
                    hr = (RmaxRf-RminRf)/NrIntvRf
                    et = escl*cos(ww*tg+theta)
                    ez0 = ez0 + extfld6(3,1)*et
                else if(FlagCart.eq.2) then
                    call getfld6xyz_EMfld3(z-zedge,extfld6xyz, extfld6xyz2, extfld6xyz3)
                    hxx = (XmaxRfg-XminRfg)/NxIntvRfg
                    hyy = (YmaxRfg-YminRfg)/NyIntvRfg
                    loc = 0.0
                    extfld = 0.0
                    call getfc2_BB(escl, ww, tg, theta, hxx, hyy, loc, extfld6xyz, extfld)
                    if (Nmultigrid.eq.2) then
                        call getfc2_BB(escl2, ww2, tg, theta2,&
                                       hxx, hyy, loc, extfld6xyz2, extfld)
                    else if (Nmultigrid.eq.3) then
                        call getfc2_BB(escl2, ww2, tg, theta2,&
                                       hxx, hyy, loc, extfld6xyz2, extfld)
                        call getfc2_BB(escl3, ww3, tg, theta3,&
                                       hxx, hyy, loc, extfld6xyz3, extfld)
                    endif
                    ez0 = ez0 + extfld(3)
                endif
            else if (beamelem%pemfld%Itype .eq. 114) then
                csid = 2
                escl = beamelem%pemfld%Param(2)
                call getfld6xyz_EMfld3(z-zedge,extfld6xyz, extfld6xyz2, extfld6xyz3)
                hxx = (XmaxRfg-XminRfg)/NxIntvRfg
                hyy = (YmaxRfg-YminRfg)/NyIntvRfg
                xoff = beamelem%pemfld%Param(13)
                yoff = beamelem%pemfld%Param(14)
                xref = beamelem%pemfld%Param(15)
                loc(1) = xoff + xref
                loc(2) = yoff
                loc(3) = 0.0
                ww = 0.0
                theta = 0.0
                extfld = 0.0
                call getfc2_BB(escl, ww, tg, theta, hxx, hyy, loc, extfld6xyz, extfld)
                ez0 = extfld(3)
            endif
        else
            csid = 3
            pos(1) = 0.0
            pos(2) = 0.0
            pos(3) = z
            pos(4) = tg
            call getfld_BeamLineElem(beamelem,pos,extfld)
            ez0 = extfld(3)
        endif

        tgam = gam + 0.5*tau*qmcc*ez0
        tgam2 = gam + tau*qmcc*ez0
        ezqmct = ez0*qmcc*tau

        if(associated(beamelem%pccl)) then
            csid = 4
            call getfldpts_CCL(rays,z,tg,textfld,beamelem%pccl,innp)
        else if (associated(beamelem%psl)) then
            csid = 5
        endif

        if (csid.eq.1) then
            escl = escl + beamelem%pemfld%Param(16)
            theta = theta + beamelem%pemfld%Param(17)*Pi/180.0
            escl2 = escl2 + beamelem%pemfld%Param(24)
            theta2 = theta2 + beamelem%pemfld%Param(25)*Pi/180.0
            escl3 = escl3 + beamelem%pemfld%Param(26)
            theta3 = theta3 + beamelem%pemfld%Param(27)*Pi/180.0
        endif

        do n = 1, innp
            select case(csid)
                case(1)
                    if(.not. FlagDisc.eq.1) then
                        pos(1) = rays(1,n)*Scxl
                        pos(2) = rays(3,n)*Scxl
                        pos(3) = z
                        pos(4) = rays(5,n) + tg
                        call getfld_BeamLineElem(beamelem,pos,extfld)
                    else
                        extfld = 0.0
                    endif
                    loc(1) = rays(1,n)*Scxl
                    loc(2) = rays(3,n)*Scxl
                    loc(3) = rays(5,n)
                    if(FlagCart.eq.1) then
                        call getfc1_BB(escl, ww, tg, theta, hr, loc, extfld6, extfld)
                    else if(FlagCart.eq.2) then
                        call getfc2_BB(escl, ww, tg, theta, hxx, hyy, loc, extfld6xyz, extfld)
                        if (Nmultigrid.eq.2) then
                            call getfc2_BB(escl2, ww2, tg, theta2,&
                                           hxx, hyy, loc, extfld6xyz2, extfld)
                        else if (Nmultigrid.eq.3) then
                            call getfc2_BB(escl2, ww2, tg, theta2,&
                                           hxx, hyy, loc, extfld6xyz2, extfld)
                            call getfc2_BB(escl3, ww3, tg, theta3,&
                                           hxx, hyy, loc, extfld6xyz3, extfld)
                        endif
                    endif
                case(2)
                    loc(1) = rays(1,n)*Scxl + xoff + xref
                    loc(2) = rays(3,n)*Scxl + yoff
                    loc(3) = rays(5,n)
                    extfld = 0.0
                    call getfc2_BB(escl, ww, tg, theta, hxx, hyy, loc, extfld6xyz, extfld)
                case(3)
                    pos(1) = rays(1,n)*Scxl
                    pos(2) = rays(3,n)*Scxl
                    pos(3) = z
                    if(associated(beamelem%pquad).and.(beamelem%pquad%Param(3) == -1)) then
                        gamn = gam - rays(6,n)
                        pos(4) = gamn/sqrt(gamn*gamn-1.0)/Clight
                    else
                        pos(4) = rays(5,n) + tg
                    endif
                    call getfld_BeamLineElem(beamelem,pos,extfld)
                case(4)
                    extfld(:) = textfld(:,n)
                case(5)
                    extfld = 0.0
            end select

            if (scflag.eq.1) then
                rcpgammai = 1.0/(-rays(6,n)+gam)
                betai = sqrt(1.0-rcpgammai*rcpgammai*(1+rays(2,n)**2+rays(4,n)**2))
                loc(1) = rays(1,n)*Scxl
                loc(2) = rays(3,n)*Scxl
                loc(3) = -gam*betai*rays(5,n)*Scxl
                call getselffld_BB(loc,gam,selffld)
                sclfld = (selfscale*selffld + extfld)*tau
            else
                sclfld = extfld*tau
            endif

            if (ANY(sclfld .ne. 0.0)) then
                call updatep_BB(rays(:,n), sclfld, ezqmct, gam, tgam, tgam2)
            endif
        enddo

        gam = gam + tau*qmcc*ez0

        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine scatter2tt_BeamBunch

        !All indices here are local to the processor.
        !find charge density on grid from particles.
        subroutine charge_BeamBunch(this,innp,innx,inny,innz,ptsgeom,&
                                    grid,chgdens,Flagbc,perdlen)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: innp,innx,inny,innz,Flagbc
        type (CompDom), intent(in) :: ptsgeom
        type (Pgrid2d), intent(in) :: grid
        double precision, intent(in) :: perdlen
        double precision,dimension(innx,inny,innz),intent(out) :: chgdens
        double precision :: t0
        !double precision, dimension(innx,inny,innz) :: tmpchg
        double precision :: recpnpt,sumtest
        integer :: totnp,nproccol,nprocrow,myid,myidx,myidy
        integer :: comm2d,commcol,commrow,ierr
        integer :: i,j,k
!        integer :: jadd,kadd,jj,jdisp

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,nproccol,nprocrow)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        ! call MPI_BARRIER(comm2d,ierr)
        ! if(myid.eq.0) then
        !   print*,"before deposit"
        ! endif

        !deposition local particles onto grids.
        !call deposit_BeamBunch(innp,innx,inny,innz,this%Pts1, &
        !tmpchg,ptsgeom,nprocrow,nproccol,myidx,myidy)
        if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
          call depositr_BeamBunch(innp,innx,inny,innz,this%Pts1, &
          chgdens,ptsgeom,nprocrow,nproccol,myidx,myidy)
!        else if(Flagbc.eq.2) then
!          call depositperd_BeamBunch(innp,innx,inny,innz,this%Pts1, &
!          chgdens,ptsgeom,nprocrow,nproccol,myidx,myidy,perdlen)
        else
          call deposit_BeamBunch(innp,innx,inny,innz,this%Pts1, &
          chgdens,ptsgeom,nprocrow,nproccol,myidx,myidy)
        endif

        !call MPI_BARRIER(comm2d,ierr)
        !if(myid.eq.0) then
        !  print*,"after deposit, before guard"
        !endif

        !sum up the contribution from neighboring guard grid to
        !the boundary grid in CIC scheme.
        if((Flagbc.eq.1).or.(Flagbc.eq.5)) then  ! 3D open
!          if(totnp.ne.1) then
            !gather contributions from guard cells.
            call guardsum1_Fldmger(chgdens,innx,inny,innz,grid)
!          endif
        else if((Flagbc.eq.2).or.(Flagbc.eq.6)) then  ! 2D open, 1D periodic
          call guardsum2_Fldmger(chgdens,innx,inny,innz,grid)
        else if(Flagbc.eq.3) then
          call guardsum3_Fldmger(chgdens,innx,inny,innz,grid)
        else if(Flagbc.eq.4) then
          call guardsum4_Fldmger(chgdens,innx,inny,innz,grid)
        else
          print*,"no such type of boundary conditions!!!"
          stop
        endif
!        recpnpt = 1.0/this%Npt
        recpnpt = 1.0 !//this%Npt has been included in rays(8,n)
!        sumtest = 0.0
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              !this%ChargeDensity(i,j,k) = tmpchg(i,j,k)*recpnpt
              chgdens(i,j,k) = chgdens(i,j,k) &
                                          *recpnpt
!              sumtest = sumtest + chgdens(i,j,k)
            enddo
          enddo
        enddo

        t_charge = t_charge + elapsedtime_Timer(t0)

        end subroutine charge_BeamBunch

        ! set local # of particles.
        subroutine setnpt_BeamBunch(this,innpt)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innpt
        type (BeamBunch), intent(out) :: this

        this%Nptlocal = innpt

        end subroutine setnpt_BeamBunch

        ! get local # of particles.
        subroutine getnpt_BeamBunch(this,outnpt)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: this
        integer, intent(out) :: outnpt

        outnpt = this%Nptlocal

        end subroutine getnpt_BeamBunch

        ! deposit particles onto grid.
        subroutine deposit_BeamBunch(innp,innx,inny,innz,rays,rho,&
                                ptsgeom,npx,npy,myidx,myidy)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,myidy
        double precision, intent (in), dimension (9, innp) :: rays
        double precision, intent (out), dimension (innx,inny,innz) :: rho
!        logical, intent (in), dimension (innp) :: msk
        integer :: ix,jx,kx,ix1,jx1,kx1
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab,cd,ef
        integer :: ngood, i, j, k, kadd, jadd
        double precision :: t0
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        type (CompDom) :: ptsgeom

        call starttime_Timer( t0 )

        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0/msize(1)
        hyi = 1.0/msize(2)
        hzi = 1.0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        rho=0.
        do i = 1, innp
          ix=(rays(1,i)-xmin)*hxi + 1
          ab=((xmin-rays(1,i))+ix*hx)*hxi
          jx=(rays(3,i)-ymin)*hyi + 1 + jadd
          cd=((ymin-rays(3,i))+(jx-jadd)*hy)*hyi
          kx=(rays(5,i)-zmin)*hzi + 1 + kadd
          ef=((zmin-rays(5,i))+(kx-kadd)*hz)*hzi
          ix1=ix+1
          jx1=jx+1
          kx1=kx+1
          ! (i,j,k):
          rho(ix,jx,kx) = rho(ix,jx,kx) + ab*cd*ef*rays(8,i)
          ! (i,j+1,k):
          rho(ix,jx1,kx) = rho(ix,jx1,kx) + ab*(1.0-cd)*ef*rays(8,i)
          ! (i,j+1,k+1):
          rho(ix,jx1,kx1) = rho(ix,jx1,kx1)+ab*(1.0-cd)*(1.0-ef)*rays(8,i)
          ! (i,j,k+1):
          rho(ix,jx,kx1) = rho(ix,jx,kx1)+ab*cd*(1.0-ef)*rays(8,i)
          ! (i+1,j,k+1):
          rho(ix1,jx,kx1) = rho(ix1,jx,kx1)+(1.0-ab)*cd*(1.0-ef)*rays(8,i)
          ! (i+1,j+1,k+1):
          rho(ix1,jx1,kx1) = rho(ix1,jx1,kx1)+(1.0-ab)*(1.0-cd)*(1.0-ef)*rays(8,i)
          ! (i+1,j+1,k):
          rho(ix1,jx1,kx) = rho(ix1,jx1,kx)+(1.0-ab)*(1.0-cd)*ef*rays(8,i)
          ! (i+1,j,k):
          rho(ix1,jx,kx) = rho(ix1,jx,kx)+(1.0-ab)*cd*ef*rays(8,i)
        enddo

        t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

        end subroutine deposit_BeamBunch

        subroutine depositperd_BeamBunch(innp,innx,inny,innz,rays,rho,&
                                ptsgeom,npx,npy,myidx,myidy,perd)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,myidy
        double precision, intent (in), dimension (9, innp) :: rays
        double precision, intent (out), dimension (innx,inny,innz) :: rho
        double precision, intent(in) :: perd
!        logical, intent (in), dimension (innp) :: msk
        integer :: ix,jx,kx,ix1,jx1,kx1
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab,cd,ef
        integer :: ngood, i, j, k, kadd, jadd
        double precision :: t0
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        type (CompDom) :: ptsgeom
        double precision :: tmpz,halfperd
        integer :: tmpk1

        call starttime_Timer( t0 )

        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0/msize(1)
        hyi = 1.0/msize(2)
        hzi = 1.0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        rho=0.
        halfperd = perd/2
        do i = 1, innp
          ix=(rays(1,i)-xmin)*hxi + 1
          ab=((xmin-rays(1,i))+ix*hx)*hxi
          jx=(rays(3,i)-ymin)*hyi + 1 + jadd
          cd=((ymin-rays(3,i))+(jx-jadd)*hy)*hyi
          tmpk1 = rays(5,i)/halfperd
          tmpz = mod(rays(5,i),halfperd)-mod(tmpk1,2)*halfperd
!          kx=(rays(5,i)-zmin)*hzi + 1 + kadd
!          ef=((zmin-rays(5,i))+(kx-kadd)*hz)*hzi
          kx=(tmpz-zmin)*hzi + 1 + kadd
          ef=((zmin-tmpz)+(kx-kadd)*hz)*hzi
          ix1=ix+1
          jx1=jx+1
          kx1=kx+1
          ! (i,j,k):
          rho(ix,jx,kx) = rho(ix,jx,kx) + ab*cd*ef*rays(8,i)
          ! (i,j+1,k):
          rho(ix,jx1,kx) = rho(ix,jx1,kx) + ab*(1.0-cd)*ef*rays(8,i)
          ! (i,j+1,k+1):
          rho(ix,jx1,kx1) = rho(ix,jx1,kx1)+ab*(1.0-cd)*(1.0-ef)*rays(8,i)
          ! (i,j,k+1):
          rho(ix,jx,kx1) = rho(ix,jx,kx1)+ab*cd*(1.0-ef)*rays(8,i)
          ! (i+1,j,k+1):
          rho(ix1,jx,kx1) = rho(ix1,jx,kx1)+(1.0-ab)*cd*(1.0-ef)*rays(8,i)
          ! (i+1,j+1,k+1):
          rho(ix1,jx1,kx1) = rho(ix1,jx1,kx1)+(1.0-ab)*(1.0-cd)*(1.0-ef)*rays(8,i)
          ! (i+1,j+1,k):
          rho(ix1,jx1,kx) = rho(ix1,jx1,kx)+(1.0-ab)*(1.0-cd)*ef*rays(8,i)
          ! (i+1,j,k):
          rho(ix1,jx,kx) = rho(ix1,jx,kx)+(1.0-ab)*cd*ef*rays(8,i)
        enddo

        t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

        end subroutine depositperd_BeamBunch

        ! deposit particles onto grid.
        subroutine depositr_BeamBunch(innp,innx,inny,innz,rays,rho,&
                                ptsgeom,npx,npy,myidx,myidy)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,myidy
        double precision, intent (in), dimension (9, innp) :: rays
        double precision, intent (out), dimension (innx,inny,innz) :: rho
        type (CompDom) :: ptsgeom
!        logical, intent (in), dimension (innp) :: msk
        integer :: ix,jx,kx,ix1,jx1,kx1
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab,cd,ef
        integer :: ngood, i, j, k, kadd, jadd
        double precision :: t0,ri,thi,pi
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range

        call starttime_Timer( t0 )

        pi = 2*asin(1.0)

        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0/msize(1)
        hyi = 1.0/msize(2)
        hzi = 1.0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        rho=0.
        do i = 1, innp
          ri = sqrt(rays(1,i)*rays(1,i)+rays(3,i)*rays(3,i))
          if(rays(1,i).gt.0.0) then
            if(rays(3,i).gt.0.0) then
              thi = asin(rays(3,i)/ri)
            else
              thi = 2*pi+asin(rays(3,i)/ri)
            endif
          else
            thi = pi - asin(rays(3,i)/ri)
          endif
          ix=(ri-xmin)*hxi + 1
          ab=(ix*ix*hx*hx-(ri-xmin)*(ri-xmin))/ &
             (ix*ix*hx*hx-(ix-1)*(ix-1)*hx*hx)
          jx=(thi-ymin)*hyi + 1 + jadd
          cd=((ymin-thi)+(jx-jadd)*hy)*hyi
          kx=(rays(5,i)-zmin)*hzi + 1 + kadd
          kx=mod(kx,innz)
          ef=((zmin-rays(5,i))+(kx-kadd)*hz)*hzi
          ix1=ix+1
          jx1=jx+1
          kx1=kx+1
          ! (i,j,k):
          rho(ix,jx,kx) = rho(ix,jx,kx) + ab*cd*ef*rays(8,i)
          ! (i,j+1,k):
          rho(ix,jx1,kx) = rho(ix,jx1,kx) + ab*(1.0-cd)*ef*rays(8,i)
          ! (i,j+1,k+1):
          rho(ix,jx1,kx1) = rho(ix,jx1,kx1)+ab*(1.0-cd)*(1.0-ef)*rays(8,i)
          ! (i,j,k+1):
          rho(ix,jx,kx1) = rho(ix,jx,kx1)+ab*cd*(1.0-ef)*rays(8,i)
          ! (i+1,j,k+1):
          rho(ix1,jx,kx1) = rho(ix1,jx,kx1)+(1.0-ab)*cd*(1.0-ef)*rays(8,i)
          ! (i+1,j+1,k+1):
          rho(ix1,jx1,kx1) = rho(ix1,jx1,kx1)+(1.0-ab)*(1.0-cd)*(1.0-ef)*rays(8,i)
          ! (i+1,j+1,k):
          rho(ix1,jx1,kx) = rho(ix1,jx1,kx)+(1.0-ab)*(1.0-cd)*ef*rays(8,i)
          ! (i+1,j,k):
          rho(ix1,jx,kx) = rho(ix1,jx,kx)+(1.0-ab)*cd*ef*rays(8,i)
        enddo

        t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

        end subroutine depositr_BeamBunch

        subroutine destruct_BeamBunch(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(out) :: this

        deallocate(this%Pts1)

        end subroutine destruct_BeamBunch

        !rotate to the particle coordinates to local beam coordinate of "ptref".
        !Here, the local "z" direction has been enlarged by "gamma" for the space-charge
        !calculation.
        !Here, "ptref" coordinate has a rotation "theta" in x-z plane with
        !respect to the orginal coordinate.
        subroutine rotto_BeamBunch(this,ptref,ptrange)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, dimension(6) :: ptref
        double precision, dimension(6), intent(out) :: ptrange
        double precision  :: cs,ss,gamma
        double precision, dimension(6) :: temp
        integer :: i

        gamma = sqrt(1.0+ptref(6)**2+ptref(2)**2)
        cs = ptref(6)/sqrt(ptref(6)**2 + ptref(2)**2)
        ss = ptref(2)/sqrt(ptref(6)**2 + ptref(2)**2)
        do i = 1, 3
          ptrange(2*i-1) = 1.0e20
          ptrange(2*i) = -1.0e20
        enddo

        do i = 1, this%Nptlocal
          temp(1) = this%Pts1(1,i)*Scxl
          temp(2) = this%Pts1(2,i)
          temp(3) = this%Pts1(3,i)*Scxl
          temp(4) = this%Pts1(4,i)
          temp(5) = this%Pts1(5,i)*Scxl
          temp(6) = this%Pts1(6,i)
          this%Pts1(1,i) = temp(1)*cs - temp(5)*ss
          this%Pts1(2,i) = temp(2)*cs - temp(6)*ss
          this%Pts1(3,i) = temp(3)
          this%Pts1(4,i) = temp(4)
          this%Pts1(5,i) = (temp(1)*ss + temp(5)*cs)*gamma
          this%Pts1(6,i) = temp(2)*ss + temp(6)*cs
          if(ptrange(1).gt.this%Pts1(1,i)) then
            ptrange(1) = this%Pts1(1,i)
          endif
          if(ptrange(2).lt.this%Pts1(1,i)) then
            ptrange(2) = this%Pts1(1,i)
          endif
          if(ptrange(3).gt.this%Pts1(3,i)) then
            ptrange(3) = this%Pts1(3,i)
          endif
          if(ptrange(4).lt.this%Pts1(3,i)) then
            ptrange(4) = this%Pts1(3,i)
          endif
          if(ptrange(5).gt.this%Pts1(5,i)) then
            ptrange(5) = this%Pts1(5,i)
          endif
          if(ptrange(6).lt.this%Pts1(5,i)) then
            ptrange(6) = this%Pts1(5,i)
          endif
        enddo

        end subroutine rotto_BeamBunch
        subroutine rotback_BeamBunch(this,ptref)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, dimension(6) :: ptref
        double precision, dimension(6) :: tmp
        double precision  :: cs,ss,gamma
        integer :: i

        gamma = sqrt(1.0+ptref(6)**2+ptref(2)**2)
        cs = ptref(6)/sqrt(ptref(6)**2 + ptref(2)**2)
        ss = ptref(2)/sqrt(ptref(6)**2 + ptref(2)**2)
        do i = 1, this%Nptlocal
          tmp(1) = this%Pts1(1,i)*cs + this%Pts1(5,i)*ss/gamma
          tmp(2) = this%Pts1(2,i)*cs + this%Pts1(6,i)*ss
          tmp(3) = this%Pts1(3,i)
          tmp(4) = this%Pts1(4,i)
          tmp(5) = -this%Pts1(1,i)*ss + this%Pts1(5,i)*cs/gamma
          tmp(6) = -this%Pts1(2,i)*ss + this%Pts1(6,i)*cs
          this%Pts1(1,i) = tmp(1)/Scxl
          this%Pts1(2,i) = tmp(2)
          this%Pts1(3,i) = tmp(3)/Scxl
          this%Pts1(4,i) = tmp(4)
          this%Pts1(5,i) = tmp(5)/Scxl
          this%Pts1(6,i) = tmp(6)
        enddo

        end subroutine rotback_BeamBunch

        !convert from z coordinates (x,px,y,py,phase,pt0-pt) to
        !t frame coordinates (x,px,y,py,z,pz). Here, x normalized by xl, px = gamma beta_x.
        subroutine convZT_BeamBunch(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer :: i
        double precision :: gam,rcpgammai,betai

        gam = -this%refptcl(6)

        ! The following steps go from z to t frame.
        ! 2) Linear algorithm to transfer from z to t frame.
        do i = 1, this%Nptlocal
          rcpgammai = 1.0/(-this%Pts1(6,i)+gam)
          betai = sqrt(1.0-rcpgammai*rcpgammai*(1+this%Pts1(2,i)**2+ &
                      this%Pts1(4,i)**2) )
          this%Pts1(1,i) = this%Pts1(1,i)-this%Pts1(2,i) &
                           *this%Pts1(5,i)*rcpgammai  !// x/L
          this%Pts1(2,i) = this%Pts1(2,i) !//gamma beta_x
          this%Pts1(3,i) = this%Pts1(3,i)-this%Pts1(4,i) &
                           *this%Pts1(5,i)*rcpgammai  !// y/L
          this%Pts1(4,i) = this%Pts1(4,i)  !//gamma beta_y
          this%Pts1(5,i) = -betai*this%Pts1(5,i) !// z/L
          this%Pts1(6,i) = betai/rcpgammai  !//gamma beta_z
        enddo

        this%refptcl = 0.0
        this%refptcl(6) = sqrt(gam**2 - 1.0)
!        do i = 1, this%Nptlocal
!          write(49,110)this%Pts1(1:6,i)
!        enddo
!110     format(6(1x,e14.6))

        end subroutine convZT_BeamBunch

        !Linear algorithm to transfer from t to z frame after bend.
        !The Cartesian coordinate has been rotated following the exit direction
        !of the reference particle. The exit angle of the reference particle
        !should correspond to the bend angle.
        subroutine convTZ_BeamBunch(this,tout)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision :: tout
        integer :: i
        double precision :: gamma0
        double precision :: rcpgammai,betai,cs,ss
        double precision, dimension(6) :: temp

!        do i = 1, this%Nptlocal
!          write(50,110)this%Pts1(1:6,i)
!        enddo
!110     format(6(1x,e14.6))
        !rotate and shift to the local coordinates of the reference particle.
        !However, there is no shift of the momentum
        cs = this%refptcl(6)/sqrt(this%refptcl(6)**2 + this%refptcl(2)**2)
        ss = this%refptcl(2)/sqrt(this%refptcl(6)**2 + this%refptcl(2)**2)
        do i = 1, this%Nptlocal
          temp(1) = this%Pts1(1,i) - this%refptcl(1)
          temp(2) = this%Pts1(2,i)
          temp(3) = this%Pts1(3,i) - this%refptcl(3)
          temp(4) = this%Pts1(4,i)
          temp(5) = this%Pts1(5,i) - this%refptcl(5)
          temp(6) = this%Pts1(6,i)
          this%Pts1(1,i) = temp(1)*cs - temp(5)*ss
          this%Pts1(2,i) = temp(2)*cs - temp(6)*ss
          this%Pts1(3,i) = temp(3)
          this%Pts1(4,i) = temp(4)
          this%Pts1(5,i) = temp(1)*ss + temp(5)*cs
          this%Pts1(6,i) = temp(2)*ss + temp(6)*cs
        enddo
        temp(2) = this%refptcl(2)
        temp(6) = this%refptcl(6)
        this%refptcl(2) = temp(2)*cs - temp(6)*ss
        this%refptcl(6) = temp(2)*ss + temp(6)*cs

!        print*,"refpt: ",this%refptcl
!        do i = 1, this%Nptlocal
!          write(64,110)this%Pts1(1:6,i)
!        enddo
!111     format(6(1x,e14.6))

        !//convert from the local T frame (dx,px,dy,py,dz,pz) to Z frame
        !//(x,px,y,py,phase,pt0-pt).
        gamma0 = sqrt(1+this%refptcl(2)**2+this%refptcl(6)**2)
        do i = 1, this%Nptlocal
          rcpgammai = 1.0/sqrt(1.0+this%Pts1(2,i)**2+this%Pts1(4,i)**2+this%Pts1(6,i)**2)
          betai = this%Pts1(6,i)*rcpgammai
          this%Pts1(6,i) = gamma0 - 1.0/rcpgammai
          this%Pts1(5,i) = this%Pts1(5,i)/(-betai)
          this%Pts1(1,i) = this%Pts1(1,i)+this%Pts1(2,i) &
                            *this%Pts1(5,i)*rcpgammai
          this%Pts1(3,i) = this%Pts1(3,i)+this%Pts1(4,i) &
                            *this%Pts1(5,i)*rcpgammai
        enddo

        this%refptcl = 0.0
        this%refptcl(5) = tout
        this%refptcl(6) = -gamma0

        end subroutine convTZ_BeamBunch

        !//drift half step in positions.
        !//Here, x, y, z are normalized by C * Dt
        !//tau - normalized step size (by Dt).
        !//the time "t" is normalized by the scaling frequency.
        subroutine drifthalfT_BeamBunch(this,t,tau)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: t
        double precision, intent (in) :: tau
        double precision :: xl,pz
        double precision :: t0,recpgam
        integer :: i

        call starttime_Timer(t0)

        do i = 1, this%Nptlocal
          !//get 1.0/gamma of each particle
          recpgam = 1.0/sqrt(1.0+this%Pts1(2,i)**2+this%Pts1(4,i)**2+&
                                 this%Pts1(6,i)**2)
          this%Pts1(1,i) = this%Pts1(1,i)+0.5*tau*this%Pts1(2,i)*recpgam
          this%Pts1(3,i) = this%Pts1(3,i)+0.5*tau*this%Pts1(4,i)*recpgam
          this%Pts1(5,i) = this%Pts1(5,i)+0.5*tau*this%Pts1(6,i)*recpgam
        enddo
        recpgam = 1.0/sqrt(1.0+this%refptcl(2)**2+this%refptcl(4)**2+&
                           this%refptcl(6)**2)
        this%refptcl(1) = this%refptcl(1)+0.5*tau*this%refptcl(2)*recpgam
        this%refptcl(3) = this%refptcl(3)+0.5*tau*this%refptcl(4)*recpgam
        this%refptcl(5) = this%refptcl(5)+0.5*tau*this%refptcl(6)*recpgam

        t_map1 = t_map1 + elapsedtime_Timer(t0)

        end subroutine drifthalfT_BeamBunch

        subroutine kickT_BeamBunch(this,beamelem,z,tau,innx,inny,innz, &
                             temppotent,ptsgeom,grid,Flagbc,flagerr)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), intent(in) :: beamelem
        integer, intent(in) :: innx, inny, innz, Flagbc,flagerr
        type (CompDom), intent(in) :: ptsgeom
        double precision,dimension(innx,inny,innz),intent(inout) :: temppotent
        type (Pgrid2d), intent(in) :: grid
        double precision, intent(in) :: z, tau
        double precision, dimension(innx,inny,innz) :: egx,egy,egz
        double precision, dimension(3) :: msize
        double precision :: hxi, hyi, hzi
        double precision :: t0,tg,chge,mass,curr,gam
        integer :: totnp,nproccol,nprocrow,myid,myidx,myidy
        integer :: i, j, k, yadd, zadd, innp
!        integer :: comm2d,commcol,commrow

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,nproccol,nprocrow)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        if(nproccol.gt.1) then
          yadd = 1
        else
          yadd = 0
        endif
        if(nprocrow.gt.1) then
          zadd = 1
        else
          zadd = 0
        endif

        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0/msize(1)
        hyi = 1.0/msize(2)
        hzi = 1.0/msize(3)

        !call MPI_BARRIER(comm2d,ierr)
        !if(myid.eq.16) then
        !  print*,"before guardexch:"
        !endif

        curr = this%Current
        innp = this%Nptlocal
        tg = z
        mass = this%Mass
        chge = this%Charge

        if(curr.gt.0.0) then
!------------------------------------------------------------------
! current greater than 0

        ! Exchange potential information on guard grid cells.
        if((Flagbc.eq.1).or.(Flagbc.eq.5)) then ! 3D open
          if(totnp.ne.1) then
            call guardexch1_Fldmger(temppotent,innx,inny,innz,grid)
          endif
        else if((Flagbc.eq.2).or.(Flagbc.eq.6)) then ! 2D open, 1D periodic
          call guardexch2_Fldmger(temppotent,innx,inny,innz,grid)
        else if(Flagbc.eq.3) then
          call guardexch3_Fldmger(temppotent,innx,inny,innz,grid)
        else if(Flagbc.eq.4) then
          call guardexch4_Fldmger(temppotent,innx,inny,innz,grid)
        else
          print*,"no such boundary condition!!!!"
          stop
        endif

! cache optimization here.
        !print*,"yadd: ",yadd,zadd,innx,inny,innz,1.0/hxi,1.0/hyi,1.0/hzi
        !Ex
        !egx = 0.0
        do k = zadd+1, innz-zadd
          do j = yadd+1, inny-yadd
            egx(1,j,k) = hxi*(temppotent(1,j,k)-temppotent(2,j,k))
            do i = 2, innx-1
              egx(i,j,k) = 0.5*hxi*(temppotent(i-1,j,k)- &
                           temppotent(i+1,j,k))
            enddo
            egx(innx,j,k) = hxi*(temppotent(innx-1,j,k)- &
                                 temppotent(innx,j,k))
          enddo
        enddo

        !Ey
        !egy = 0.0
        if((Flagbc.eq.3).or.(Flagbc.eq.4)) then  ! periodic in theta

        if(nproccol.gt.1) then
          do k = 1, innz
            do j = 2, inny-1
              do i = 2, innx
                egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k) -  &
                             temppotent(i,j+1,k))*hxi/(i-1)
              enddo
              egy(1,j,k) = 0.0
            enddo
          enddo
        else
          do k = 1, innz
            do j = 2, inny-1
              do i = 2, innx
                egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k) -  &
                             temppotent(i,j+1,k))*hxi/(i-1)
              enddo
                egy(1,j,k) = 0.0
            enddo
          enddo
          do k = 1, innz
            do i = 2, innx
              egy(i,1,k) = 0.5*hyi*(temppotent(i,inny-1,k) -  &
                           temppotent(i,2,k))*hxi/(i-1)
              egy(i,inny,k) = egy(i,1,k)
            enddo
            egy(1,1,k) = 0.0
            egy(1,inny,k) = 0.0
          enddo
        endif

        else  ! open in Y

        if(nproccol.gt.1) then
          if(myidy.eq.0) then
            do k = zadd+1, innz-zadd
              do i = 1, innx
                egy(i,yadd+1,k) = hyi*(temppotent(i,yadd+1,k)- &
                                       temppotent(i,yadd+2,k))
              enddo
            enddo
            do k = zadd+1, innz-zadd
              do j = yadd+2, inny-yadd
                do i = 1, innx
                  egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                               temppotent(i,j+1,k))
                enddo
              enddo
            enddo
          else if(myidy.eq.(nproccol-1)) then
            do k = zadd+1, innz-zadd
              do i = 1, innx
                egy(i,inny-yadd,k) = hyi*(temppotent(i,inny-yadd-1,k)- &
                                          temppotent(i,inny-yadd,k))
              enddo
            enddo
            do k = zadd+1, innz-zadd
              do j = yadd+1, inny-yadd-1
                do i = 1, innx
                  egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                               temppotent(i,j+1,k))
                enddo
              enddo
            enddo
          else
            do k = zadd+1, innz-zadd
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                               temppotent(i,j+1,k))
                enddo
              enddo
            enddo
          endif
        else
          do k = zadd+1, innz-zadd
            do i = 1, innx
              egy(i,1,k) = hyi*(temppotent(i,1,k)-temppotent(i,2,k))
            enddo
            do j = 2, inny-1
              do i = 1, innx
                egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                             temppotent(i,j+1,k))
              enddo
            enddo
            do i = 1, innx
              egy(i,inny,k) = hyi*(temppotent(i,inny-1,k)- &
                                   temppotent(i,inny,k))
            enddo
          enddo
        endif

        endif

        !Ez
        !egz = 0.0
        if(nprocrow.gt.1) then
          if((Flagbc.eq.1).or.(Flagbc.eq.3).or.(Flagbc.eq.5)) then ! 3D open
            if(myidx.eq.0) then
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  egz(i,j,zadd+1) = hzi*(temppotent(i,j,zadd+1)- &
                                         temppotent(i,j,zadd+2))
                enddo
              enddo
              do k = zadd+2, innz-zadd
                do j = yadd+1, inny-yadd
                  do i = 1, innx
                    egz(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                 temppotent(i,j,k+1))
                  enddo
                enddo
              enddo
            else if(myidx.eq.(nprocrow-1)) then
              do k = zadd+1, innz-zadd-1
                do j = yadd+1, inny-yadd
                  do i = 1, innx
                    egz(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                 temppotent(i,j,k+1))
                  enddo
                enddo
              enddo
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  egz(i,j,innz-zadd) = hzi*(temppotent(i,j,innz-zadd-1)- &
                                            temppotent(i,j,innz-zadd))
                enddo
              enddo
            else
              do k = zadd+1, innz-zadd
                do j = yadd+1, inny-yadd
                  do i = 1, innx
                    egz(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                 temppotent(i,j,k+1))
                  enddo
                enddo
              enddo
            endif
          else if((Flagbc.eq.2).or.(Flagbc.eq.4).or.(Flagbc.eq.6)) then
          ! 2D open, 1D periodic
            do k = zadd+1, innz-zadd
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  egz(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                               temppotent(i,j,k+1))
                enddo
              enddo
            enddo
          else
            print*,"no such boundary condition!!!"
            stop
          endif
        else
          if((Flagbc.eq.1).or.(Flagbc.eq.3).or.(Flagbc.eq.5)) then ! 3D open
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,1) = hzi*(temppotent(i,j,1)-temppotent(i,j,2))
              enddo
            enddo
          else if((Flagbc.eq.2).or.(Flagbc.eq.4).or.(Flagbc.eq.6)) then
          ! 2D open, 1D periodic
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,1)=0.5*hzi*(temppotent(i,j,innz-1)-temppotent(i,j,2))
              enddo
            enddo
          else
          endif
          do k = 2, innz-1
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                      temppotent(i,j,k+1))
              enddo
            enddo
          enddo
          if((Flagbc.eq.1).or.(Flagbc.eq.3).or.(Flagbc.eq.5)) then
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,innz) = hzi*(temppotent(i,j,innz-1)- &
                                     temppotent(i,j,innz))
              enddo
            enddo
          else if((Flagbc.eq.2).or.(Flagbc.eq.4).or.(Flagbc.eq.6)) then
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,innz) = 0.5*hzi*(temppotent(i,j,innz-1)- &
                                     temppotent(i,j,2))
              enddo
            enddo
          else
          endif
        endif

        ! find field from potential.
!        egx = 0.5*hxi*(cshift(temppotent,-1,1) -  &
!                       cshift(temppotent,1,1))
!        egy = 0.5*hyi*(cshift(temppotent,-1,2) -  &
!                       cshift(temppotent,1,2))
!        egz = 0.5*hzi*(cshift(temppotent,-1,3) -  &
!                      cshift(temppotent,1,3))

!        if(totnp.eq.1) then
!          egz(:,:,2) = 0.5*hzi*(temppotent(:,:,innz-1)-&
!                                temppotent(:,:,3))
!          egz(:,:,innz-1) = 0.5*hzi*(temppotent(:,:,innz-2)-&
!                                temppotent(:,:,2))
!          egy(:,2,:) = 0.5*hyi*(temppotent(:,inny-1,:)-&
!                                temppotent(:,3,:))
!          egy(:,inny-1,:) = 0.5*hyi*(temppotent(:,inny-2,:)-&
!                                temppotent(:,2,:))
!        endif

        !Send the E field to the neibhoring guard grid to do the CIC
        !interpolation.
        if(totnp.ne.1) then
          call boundint4_Fldmger(egx,egy,egz,innx,inny,innz,grid)
        endif

        call scatterT_BeamBunch(innp,innx,inny,innz,this%Pts1,this%refptcl,&
        egx,egy,egz,ptsgeom,nprocrow,nproccol,myidx,myidy,tg,chge,mass,tau,z,&
        beamelem)

        else
!------------------------------------------------------------------
! current is 0
          call scatterT0_BeamBunch(innp,this%Pts1,this%refptcl,tg,chge,mass,tau,z,&
                                   beamelem)
        endif

        !call MPI_BARRIER(comm2d,ierr)
        !if(myid.eq.0) then
        !  print*,"after kick:"
        !endif

        t_force = t_force + elapsedtime_Timer(t0)

        end subroutine kickT_BeamBunch

        !Here, the update of the momentum is in the rotated coordinates
        !since the boost direction is not along z but in the x-z plane
        !due to the bend magnet.
        subroutine scatterT_BeamBunch(innp,innx,inny,innz,rays,refpt,exg,&
        eyg,ezg,ptsgeom,npx,npy,myidx,myidy,tg,chge,mass,dt,z,&
        beamelem)
        implicit none
        include 'mpif.h'
        double precision, intent (inout), dimension (9,innp) :: rays
        double precision, intent (inout), dimension (6) :: refpt
        double precision, intent (in), dimension (innx,inny,innz) :: exg,eyg,ezg
        double precision, intent (in) :: dt,mass,chge,tg,z
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,myidy
        type (BeamLineElem), intent(in) :: beamelem
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab, cd, ef
        integer :: n,i,ii
        double precision :: t0
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin
        double precision, dimension(3) :: msize
        double precision, dimension(4) :: pos
        double precision, dimension(6) :: range, extfld
        type (CompDom) :: ptsgeom
        integer :: ix,jx,kx,ix1,jx1,kx1,ierr,kadd,jadd
        double precision :: tmpscale,qmcc,&
        a1,a2,a3,a4,s1,s2,s3,umx,umy,umz,upx,upy,upz,tmp,&
        exn,eyn,ezn,ex,ey,ez,bx,by,bz,recpgamma
        double precision :: cs,ss,gam,beta0,coefE0,coefB0

        call starttime_Timer( t0 )

        !tmpscale = Clight*Clight*1.0e-7*curr/bfreq*chge
        tmpscale = Clight*Clight*1.0e-7 !curr*chge/bfreq are included in the charge density
        qmcc = chge/mass
        coefE0 = Scxl
        coefB0 = Scxl*Clight
        gam = sqrt(1.0+refpt(2)**2+refpt(6)**2)
        beta0 = sqrt(gam**2-1.0)/gam
        !cos and sin of the rotation angle. In a straight machine,
        !the rotation angle is 0.
        cs = refpt(6)/sqrt(refpt(2)**2+refpt(6)**2)
        ss = refpt(2)/sqrt(refpt(2)**2+refpt(6)**2)

        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0/msize(1)
        hyi = 1.0/msize(2)
        hzi = 1.0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        !update the momentum of the reference particle using the external fields
        !get the external fields at given x,y,z,t
        pos(1) = refpt(1)*Scxl
        pos(2) = refpt(3)*Scxl
        pos(3) = refpt(5)*Scxl
        pos(4) = tg
        call getfld_BeamLineElem(beamelem,pos,extfld)

        ex = extfld(1)
        ey = extfld(2)
        ez = extfld(3)
        bx = extfld(4)
        by = extfld(5)
        bz = extfld(6)
        umx = refpt(2) + coefE0*qmcc*ex*0.5*dt
        umy = refpt(4) + coefE0*qmcc*ey*0.5*dt
        umz = refpt(6) + coefE0*qmcc*ez*0.5*dt
        recpgamma = 1.0/sqrt(1.0+umx*umx+umy*umy+umz*umz)
        tmp = 0.5*coefB0*qmcc*dt*recpgamma
        a1 = tmp*bx
        a2 = tmp*by
        a3 = tmp*bz
        a4 = 1.0+a1*a1+a2*a2+a3*a3
        s1 = umx + tmp*(umy*bz-umz*by)
        s2 = umy - tmp*(umx*bz-umz*bx)
        s3 = umz + tmp*(umx*by-umy*bx)
        upx = ((1.0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
        upy = ((a1*a2-a3)*s1+(1.0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
        upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0+a3*a3)*s3)/a4
        refpt(2) = upx + coefE0*qmcc*ex*0.5*dt
        refpt(4) = upy + coefE0*qmcc*ey*0.5*dt
        refpt(6) = upz + coefE0*qmcc*ez*0.5*dt

        do n = 1, innp
          ix=(rays(1,n)-xmin)*hxi + 1
          ab=((xmin-rays(1,n))+ix*hx)*hxi
          jx=(rays(3,n)-ymin)*hyi + 1 + jadd
          cd=((ymin-rays(3,n))+(jx-jadd)*hy)*hyi
          kx=(rays(5,n)-zmin)*hzi + 1 + kadd
          ef=((zmin-rays(5,n))+(kx-kadd)*hz)*hzi

          ix1 = ix + 1
          jx1 = jx + 1
          kx1 = kx + 1

          exn = (exg(ix,jx,kx)*ab*cd*ef  &
                  +exg(ix,jx1,kx)*ab*(1.-cd)*ef &
                  +exg(ix,jx1,kx1)*ab*(1.-cd)*(1.0-ef) &
                  +exg(ix,jx,kx1)*ab*cd*(1.0-ef) &
                  +exg(ix1,jx,kx1)*(1.0-ab)*cd*(1.0-ef) &
                  +exg(ix1,jx1,kx1)*(1.0-ab)*(1.0-cd)*(1.0-ef)&
                  +exg(ix1,jx1,kx)*(1.0-ab)*(1.0-cd)*ef &
                  +exg(ix1,jx,kx)*(1.0-ab)*cd*ef)*gam

          eyn = (eyg(ix,jx,kx)*ab*cd*ef  &
                  +eyg(ix,jx1,kx)*ab*(1.-cd)*ef &
                  +eyg(ix,jx1,kx1)*ab*(1.-cd)*(1.0-ef) &
                  +eyg(ix,jx,kx1)*ab*cd*(1.0-ef) &
                  +eyg(ix1,jx,kx1)*(1.0-ab)*cd*(1.0-ef) &
                  +eyg(ix1,jx1,kx1)*(1.0-ab)*(1.0-cd)*(1.0-ef)&
                  +eyg(ix1,jx1,kx)*(1.0-ab)*(1.0-cd)*ef &
                  +eyg(ix1,jx,kx)*(1.0-ab)*cd*ef)*gam

          ezn = ezg(ix,jx,kx)*ab*cd*ef  &
                  +ezg(ix,jx1,kx)*ab*(1.-cd)*ef &
                  +ezg(ix,jx1,kx1)*ab*(1.-cd)*(1.0-ef) &
                  +ezg(ix,jx,kx1)*ab*cd*(1.0-ef) &
                  +ezg(ix1,jx,kx1)*(1.0-ab)*cd*(1.0-ef) &
                  +ezg(ix1,jx1,kx1)*(1.0-ab)*(1.0-cd)*(1.0-ef)&
                  +ezg(ix1,jx1,kx)*(1.0-ab)*(1.0-cd)*ef &
                  +ezg(ix1,jx,kx)*(1.0-ab)*cd*ef

          !get field in Cartesian coordinate from analytical function.
          !we need to rotate the coordinate back to the orginal coordinate.
          pos(1) = cs*rays(1,n) + ss*rays(5,n)/gam
          pos(2) = rays(3,n)
          pos(3) = -ss*rays(1,n) + cs*rays(5,n)/gam
          pos(4) = tg
          call getfld_BeamLineElem(beamelem,pos,extfld)

          !//Here, E and B fields are in real units.
          !//after get the external we need to transform the external field into
          !//the rotated coordinate.
          ex = tmpscale*exn+extfld(1)*cs - extfld(3)*ss
          ey = tmpscale*eyn+extfld(2)
          ez = tmpscale*ezn+extfld(1)*ss + extfld(3)*cs
          bx = -beta0/Clight*tmpscale*eyn+extfld(4)*cs-extfld(6)*ss
          by = beta0/Clight*tmpscale*exn+extfld(5)
          bz = extfld(4)*ss + extfld(6)*cs

          !//advance the momenta of particles using implicit central
          !//difference scheme from Birdall and Longdon's book.
          umx = rays(2,n) + coefE0*rays(7,n)*ex*0.5*dt;
          umy = rays(4,n) + coefE0*rays(7,n)*ey*0.5*dt;
          umz = rays(6,n) + coefE0*rays(7,n)*ez*0.5*dt;
          recpgamma = 1.0/sqrt(1.0+umx*umx+umy*umy+umz*umz);
          tmp = 0.5*coefB0*rays(7,n)*dt*recpgamma;
          a1 = tmp*bx;
          a2 = tmp*by;
          a3 = tmp*bz;
          a4 = 1.0+a1*a1+a2*a2+a3*a3;
          s1 = umx + tmp*(umy*bz-umz*by);
          s2 = umy - tmp*(umx*bz-umz*bx);
          s3 = umz + tmp*(umx*by-umy*bx);
          upx = ((1.0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4;
          upy = ((a1*a2-a3)*s1+(1.0+a2*a2)*s2+(a2*a3+a1)*s3)/a4;
          upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0+a3*a3)*s3)/a4;
          rays(2,n) = upx + coefE0*rays(7,n)*ex*0.5*dt;
          rays(4,n) = upy + coefE0*rays(7,n)*ey*0.5*dt;
          rays(6,n) = upz + coefE0*rays(7,n)*ez*0.5*dt;
        enddo

        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine scatterT_BeamBunch

        !Here, the update of the momentum is in the rotated coordinates
        !since the boost direction is not along z but in the x-z plane
        !due to the bend magnet. (Without Current)
        subroutine scatterT0_BeamBunch(innp,rays,refpt,&
        tg,chge,mass,dt,z,beamelem)
        implicit none
        include 'mpif.h'
        double precision, intent (inout), dimension (9,innp) :: rays
        double precision, intent (inout), dimension (6) :: refpt
        double precision, intent (in) :: dt,mass,chge,tg,z
        integer, intent(in) :: innp
        type (BeamLineElem), intent(in) :: beamelem
        integer :: n,i,ii
        double precision :: t0
        double precision, dimension(4) :: pos
        double precision, dimension(6) :: extfld
        double precision :: qmcc,&
        a1,a2,a3,a4,s1,s2,s3,umx,umy,umz,upx,upy,upz,tmp,&
        exn,eyn,ezn,ex,ey,ez,bx,by,bz,recpgamma
        double precision :: cs,ss,gam,coefE0,coefB0

        call starttime_Timer( t0 )

        qmcc = chge/mass
        coefE0 = Scxl
        coefB0 = Scxl*Clight
        gam = sqrt(1.0+refpt(2)**2+refpt(6)**2)
        !cos and sin of the rotation angle. In a straight machine,
        !the rotation angle is 0.
        cs = refpt(6)/sqrt(refpt(2)**2+refpt(6)**2)
        ss = refpt(2)/sqrt(refpt(2)**2+refpt(6)**2)

        !update the momentum of the reference particle using the external fields
        !get the external fields at given x,y,z,t
        pos(1) = refpt(1)*Scxl
        pos(2) = refpt(3)*Scxl
        pos(3) = refpt(5)*Scxl
        pos(4) = tg
        call getfld_BeamLineElem(beamelem,pos,extfld)

        ex = extfld(1)
        ey = extfld(2)
        ez = extfld(3)
        bx = extfld(4)
        by = extfld(5)
        bz = extfld(6)
!        print*,"ex,..",ex,ey,ez,bx,by,bz,pos(1),pos(2),pos(3)
!        print*,"coef: ",coefE0,coefB0,dt,qmcc
!        print*,"refpt: ",refpt(1),refpt(2),refpt(3),refpt(4),refpt(5),refpt(6)
        umx = refpt(2) + coefE0*qmcc*ex*0.5*dt
        umy = refpt(4) + coefE0*qmcc*ey*0.5*dt
        umz = refpt(6) + coefE0*qmcc*ez*0.5*dt
        recpgamma = 1.0/sqrt(1.0+umx*umx+umy*umy+umz*umz)
        tmp = 0.5*coefB0*qmcc*dt*recpgamma
        a1 = tmp*bx
        a2 = tmp*by
        a3 = tmp*bz
        a4 = 1.0+a1*a1+a2*a2+a3*a3
        s1 = umx + tmp*(umy*bz-umz*by)
        s2 = umy - tmp*(umx*bz-umz*bx)
        s3 = umz + tmp*(umx*by-umy*bx)
        upx = ((1.0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
        upy = ((a1*a2-a3)*s1+(1.0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
        upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0+a3*a3)*s3)/a4
        refpt(2) = upx + coefE0*qmcc*ex*0.5*dt
        refpt(4) = upy + coefE0*qmcc*ey*0.5*dt
        refpt(6) = upz + coefE0*qmcc*ez*0.5*dt

!        print*,"qmcc: ",rays(7,1),rays(7,2),rays(7,3),innp,qmcc
        do n = 1, innp
          !get field in Cartesian coordinate from analytical function.
          !we need to rotate the coordinate back to the orginal coordinate.
          pos(1) = cs*rays(1,n) + ss*rays(5,n)/gam
          pos(2) = rays(3,n)
          pos(3) = -ss*rays(1,n) + cs*rays(5,n)/gam
          pos(4) = tg
          call getfld_BeamLineElem(beamelem,pos,extfld)

          !//Here, E and B fields are in real units.
          !//after get the external we need to transform the external field into
          !//the rotated coordinate.
          ex = extfld(1)*cs - extfld(3)*ss
          ey = extfld(2)
          ez = extfld(1)*ss + extfld(3)*cs
          bx = extfld(4)*cs-extfld(6)*ss
          by = extfld(5)
          bz = extfld(4)*ss+extfld(6)*cs

          !//advance the momenta of particles using implicit central
          !//difference scheme from Birdall and Longdon's book.
          umx = rays(2,n) + coefE0*rays(7,n)*ex*0.5*dt;
          umy = rays(4,n) + coefE0*rays(7,n)*ey*0.5*dt;
          umz = rays(6,n) + coefE0*rays(7,n)*ez*0.5*dt;
          recpgamma = 1.0/sqrt(1.0+umx*umx+umy*umy+umz*umz);
          tmp = 0.5*coefB0*rays(7,n)*dt*recpgamma;
          a1 = tmp*bx;
          a2 = tmp*by;
          a3 = tmp*bz;
          a4 = 1.0+a1*a1+a2*a2+a3*a3;
          s1 = umx + tmp*(umy*bz-umz*by);
          s2 = umy - tmp*(umx*bz-umz*bx);
          s3 = umz + tmp*(umx*by-umy*bx);
          upx = ((1.0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4;
          upy = ((a1*a2-a3)*s1+(1.0+a2*a2)*s2+(a2*a3+a1)*s3)/a4;
          upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0+a3*a3)*s3)/a4;
          rays(2,n) = upx + coefE0*rays(7,n)*ex*0.5*dt;
          rays(4,n) = upy + coefE0*rays(7,n)*ey*0.5*dt;
          rays(6,n) = upz + coefE0*rays(7,n)*ez*0.5*dt;
        enddo

        !print*,"ex2: : ",ex,ey,ez,bx,by,bz
        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine scatterT0_BeamBunch
      end module BeamBunchclass
