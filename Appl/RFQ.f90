!----------------------------------------------------------------
! (c) Copyright, 2005 by the Regents of the University of California.
! RFQclass: Radio-Frequency Quadrupole(RFQ) beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, 7/25/05
! Description: This class defines the linear transfer map and RF field
!              for the RFQ beam line elment.
! Comments:
!----------------------------------------------------------------
      module RFQclass
        use PhysConstclass
        use Dataclass
        use Besselclass
        integer, private, parameter :: Nparam = 13
        type RFQ
          !Itype = 106
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : scale
          !      (3) : RF frequency
          !      (4) : theta0
          !      (5) : file ID
          !      (6) : radius
          !      (7) : modulation
          !      (8) : x misalignment error
          !      (9) : y misalignment error
          !      (10) : rotation error x
          !      (11) : rotation error y
          !      (12) : rotation error z
          !      (13) : cell ID
           double precision :: XX, AA !defined in Wangler's book for two term potential
        end type RFQ
        interface getparam_RFQ
          module procedure getparam1_RFQ,  &
                          getparam2_RFQ,   &
                          getparam3_RFQ
        end interface
        interface setparam_RFQ
          module procedure setparam1_RFQ,  &
                          setparam2_RFQ, setparam3_RFQ
        end interface
      contains
        subroutine construct_RFQ(this,numseg,nmpstp,type,blength,aa,xx)
        implicit none
        type (RFQ), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        real*8 :: aa,xx

        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength
        this%AA = aa
        this%XX = xx

        this%Param = 0.0

        end subroutine construct_RFQ

        subroutine setparam1_RFQ(this,i,value)
        implicit none
        type (RFQ), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_RFQ

        subroutine setparam2_RFQ(this,values)
        implicit none
        type (RFQ), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_RFQ

        subroutine setparam3_RFQ(this,numseg,nmpstp,type,blength)
        implicit none
        type (RFQ), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength

        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_RFQ

        subroutine getparam1_RFQ(this,i,blparam)
        implicit none
        type (RFQ), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_RFQ

        subroutine getparam2_RFQ(this,blparams)
        implicit none
        type (RFQ), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_RFQ

        subroutine getparam3_RFQ(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (RFQ), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_RFQ

        subroutine maplinear_RFQ(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (RFQ), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,escale,ww,theta0,qmcc
        double precision :: ein,e1in,e2in,ef,e1f,e2f
        double precision :: sinthi,costhi,uprimi,qpwi,tfin,sinthf,costhf
        double precision :: uprimf,qpwf
        integer :: mpstp

        xm = 0.0

        zedge = this%Param(1)
        escale = this%Param(2)
        ww = this%Param(3)/Scfreq
        theta0 = this%Param(4)*asin(1.0)/90
        qmcc = Bchg/Bmass
        mpstp = this%Mapstp

        y(1)=0.0
        y(2)=0.0
        y(3)=0.0
        y(4)=0.0
        y(5)=refpt(5)
        y(6)=refpt(6)
        y(7)=1.0
        y(8)=0.0
        y(9)=0.0
        y(10)=1.0
        y(11)=1.0
        y(12)=0.0
        y(13)=0.0
        y(14)=1.0
        y(15)=1.0
        y(16)=0.0
        y(17)=0.0
        y(18)=1.0
        h=tau/mpstp
        t0=t
        gi=-refpt(6)
        betai=sqrt(gi**2-1.)/gi
        ui=gi*betai
        squi=sqrt(ui)
        sqi3=sqrt(ui)**3

        call getaxfldE_RFQ(t,this,ein,e1in,e2in)
        sinthi=sin(ww*refpt(5)+theta0)
        costhi=cos(ww*refpt(5)+theta0)
        uprimi=qmcc*ein/betai*costhi
        qpwi=0.5*qmcc*Scxl/(ui*ww)
        dlti=Scxl*(0.5*uprimi/ui-qpwi*e1in*sinthi)
        thli=1.5*Scxl*(uprimi/ui)

        call rk6i_RFQ(h,mpstp,t0,y,18,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3
        tfin = t + tau
        call getaxfldE_RFQ(tfin,this,ef,e1f,e2f)
        sinthf=sin(ww*refpt(5)+theta0)
        costhf=cos(ww*refpt(5)+theta0)
        uprimf=qmcc*ef/betaf*costhf
        qpwf=0.5*qmcc*Scxl/(uf*ww)
        dltf=Scxl*(0.5*uprimf/uf-qpwf*e1f*sinthf)
        thlf=1.5*Scxl*(uprimf/uf)

        xm(1,1)= y(7)*squi/squf+y(9)*squi/squf*dlti
        xm(2,1)=(y(8)-y(7)*dltf)*squi*squf+(y(10)-y(9)*dltf)*squi*squf*&
                 dlti
        xm(1,2)= y(9)/(squi*squf)
        xm(2,2)=(y(10)-y(9)*dltf)*squf/squi
        xm(3,3)= y(11)*squi/squf+y(13)*squi/squf*dlti
        xm(4,3)= &
        (y(12)-y(11)*dltf)*squi*squf+(y(14)-y(13)*dltf)*squi*squf*dlti
        xm(3,4)= y(13)/(squi*squf)
        xm(4,4)=(y(14)-y(13)*dltf)*squf/squi
        xm(5,5)= y(15)*sqi3/sqf3+y(17)*sqi3/sqf3*thli
        xm(6,5)= &
        (y(16)-y(15)*thlf)*sqi3*sqf3+(y(18)-y(17)*thlf)*sqi3*sqf3*thli
        xm(5,6)= y(17)/(sqi3*sqf3)
        xm(6,6)=(y(18)-y(17)*thlf)*sqf3/sqi3

        end subroutine maplinear_RFQ

        subroutine rk6i_RFQ(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (RFQ), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_RFQ(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_RFQ(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_RFQ(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_RFQ(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_RFQ(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_RFQ(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_RFQ(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_RFQ(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_RFQ

        subroutine intfunc1_RFQ(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (RFQ), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55,s11tmp
        double precision :: zedge,escale,ez1,ezp1,ezpp1,ww,theta0,qmcc,&
                            sinphi,cosphi,rfdsgn
        integer :: my_rank, ierr

        zedge = this%Param(1)
        call getaxfldE_RFQ(t,this,ez1,ezp1,ezpp1)
        ww = this%Param(3)/Scfreq
        theta0 = this%Param(4)*asin(1.0)/90
        qmcc = Bchg/Bmass

        f(1) = 0.0
        f(2) = 0.0
        f(3) = 0.0
        f(4) = 0.0

        ! synchronous particle:
        gamma0=-y(6)
        beta0=sqrt(gamma0**2-1.)/gamma0
        gbet=sqrt((gamma0-1.0)*(gamma0+1.0))
        sinphi=sin(ww*y(5)+theta0)
        cosphi=cos(ww*y(5)+theta0)
        rfdsgn=ez1*cosphi
        f(5)=1.0/(beta0*Scxl)
        f(6)=-qmcc*rfdsgn

        ! matrix elements
        s11tmp=0.5e0*(1.0+0.5*gamma0**2)* &
         (qmcc*rfdsgn/beta0**2/gamma0**2)**2 + &
          qmcc*0.5/beta0**3/gamma0**3*ez1*sinphi*ww/Scxl
        s11=s11tmp*Scxl
        s33=s11tmp*Scxl
        s55=  &
         -1.5e0*qmcc/beta0**2/gamma0*ezp1*cosphi &
         +(beta0**2+0.5)/beta0**3/gamma0*qmcc*ez1*sinphi*ww/Scxl &
         +1.5e0*(1.0-0.5*gamma0**2)* &
         (qmcc*rfdsgn/beta0**2/gamma0**2)**2
        s55=Scxl*s55

        f(7)=y(8)/Scxl
        f(8)=-s11*y(7)
        f(9)=y(10)/Scxl
        f(10)=-s11*y(9)
        f(11)=y(12)/Scxl
        f(12)=-s33*y(11)
        f(13)=y(14)/Scxl
        f(14)=-s33*y(13)
        f(15)=y(16)/Scxl
        f(16)=-s55*y(15)
        f(17)=y(18)/Scxl
        f(18)=-s55*y(17)

        end subroutine intfunc1_RFQ

        !interpolate the field from the RFQ rf cavity onto bunch location.
        subroutine getaxfldE_RFQ(z,this,ez1,ezp1,ezpp1)
        implicit none
        include 'mpif.h'
        type (RFQ), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: ez1,ezp1,ezpp1
        double precision:: zz,hstep,slope,zedge,escale,pi,len,wk
        integer :: klo,khi,k
        integer :: my_rank,ierr

        pi = 2*asin(1.0)
        len = this%Length
        wk = 2*pi/len
        zedge = this%Param(1)
        escale = this%Param(2)
        zz=z-zedge

        ez1 = wk*escale*this%AA/2*sin(wk*zz)
        ezp1 = wk*wk*escale*this%AA/2*cos(wk*zz)
        ezpp1 = -wk*wk*wk*escale*this%AA/2*sin(wk*zz)

        end subroutine getaxfldE_RFQ

        !get external field with displacement and rotation errors.
        subroutine  getflderr_RFQ(pos,extfld,this,dx,dy,anglex,angley,&
                                 anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (RFQ), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl,wk
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz,r2,zmid,ra
        double precision, dimension(3) :: temp,tmp
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        wk = 2*pi/len
        ww = this%Param(3)/Scfreq
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        ra = this%Param(6)

        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = pos(3) - zedge
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

! here, the external field is defined by a term potential in Tom Wangler's book, p. 233.
        r2 = pos(1)**2+pos(2)**2
        tmpex = -this%XX*escale*pos(1)/ra**2 - wk*this%AA*escale*wk*pos(1)*cos(wk*zz)/4
        tmpex = this%XX*escale*pos(2)/ra**2 - wk*this%AA*escale*wk*pos(2)*cos(wk*zz)/4
        tmpez = wk*escale*this%AA/2*(1+wk*wk*r2/4)*sin(wk*zz)

        tt = pos(4)
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)

        !for E field
        tmp(1) = extfld(1)
        tmp(2) = extfld(2)*cos(anglex)-extfld(3)*sin(anglex)
        tmp(3) = extfld(2)*sin(anglex)+extfld(3)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(1) = (temp(1)*cos(anglez) - temp(2)*sin(anglez))*tmpcos
        extfld(2) = (temp(1)*sin(anglez) + temp(2)*cos(anglez))*tmpcos
        extfld(3) = temp(3)*tmpcos

        end subroutine getflderr_RFQ

        !get external field without displacement and rotation errors
        subroutine  getfld_RFQ(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (RFQ), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        integer, save :: nc, stdcnt, celltype
        double precision, save :: zedge, escale, zlen, ww, xl, theta0, ra, m,&
                            r0, P1, P2, P3, P4, P5, P6, P7, P8,&
                            oldA10, oldA0, oldA12I4, oldA21I2, oldA1, oldscl, oldzlen,&
                            A0, A1, A10, A12, A30, A21, A32, A23, T10, T30, U10, U30,&
                            alpha, wk, zz, r, theta, tmper, tmpethe,&
                            swkz, swkz3, cwkz, cwkz3, B1, B2, oldA12, oldA21,&
                            s2t, s4t, s6t, c2t, c4t, c6t, tmpex, tmpey, tmpez,&
                            tmpsin, tmpcos, tt

        extfld = 0.0

        if (nc.ne.tmprfqncl) then
            nc = int(this%Param(13) + 0.1)
            zlen = Rfqzlen(nc)
            celltype = Rfqltyp(nc)
            zedge = this%Param(1) + Rfqzpos(nc)
            escale = this%Param(2)*Rfqpara(1,nc)
            ww = this%Param(3)/Scfreq
            theta0 = (this%Param(4)+Rfqpara(2,nc))*pi/180.0
            ra = Rfqpara(3,nc)
            m = Rfqpara(4,nc)
            r0 = Rfqpara(5,nc)
            P1 = Rfqpara(6,nc) ! A10
            P2 = Rfqpara(7,nc) ! A0
            P3 = Rfqpara(8,nc) ! A12I4
            P4 = Rfqpara(9,nc) ! A1
            P5 = Rfqpara(10,nc) ! A30I0
            P6 = Rfqpara(11,nc) ! A21I2
            P7 = Rfqpara(12,nc) ! A32I4
            P8 = Rfqpara(13,nc) ! A23I6
            stdcnt = int(Rfqpara(14,nc) + 0.1)
            
            if (celltype.eq.1) then
                wk = Pi/zlen
                A12 = P3/bessi4(2.0*wk*r0)
                A30 = P5/bessi0(3.0*wk*r0)
                A21 = P6/bessi2(2.0*wk*r0)
                A32 = P7/bessi4(3.0*wk*r0)
                A23 = P8/bessi6(3.0*wk*r0)

                if (stdcnt.ge.2) then
                    oldscl = this%Param(2)*Rfqpara(1,nc-1)
                    oldA10 = Rfqpara(6,nc-1)
                    oldA0 = Rfqpara(7,nc-1)
                    oldA12I4 = Rfqpara(8,nc-1)
                    oldA1 = Rfqpara(9,nc-1)
                    oldA21I2 = Rfqpara(11,nc-1)
                    oldzlen = Rfqzlen(nc-1)
                    
                    oldA12 = oldA12I4/bessi4(2.0*wk*r0)
                    oldA21 = oldA21I2/bessi2(2.0*wk*r0)

                    A10 = 0.5*(P1 + oldA10)
                    A0  = 0.5*(P2 + oldA0)
                    A12 = 0.5*(A12 + oldA12)
                    A1  = 0.5*(P4 + oldA1)
                    A21 = 0.5*(A21 + oldA21)
                    escale = 0.5*(escale + oldscl)
                    wk = 2.0*Pi/(zlen+oldzlen)
                else
                    A10 = P1
                    A0  = P2
                    A1  = P4
                endif
            else if (celltype.eq.0) then
                wk = Pi/zlen
                
                if (stdcnt.ge.2) then
                    oldscl = this%Param(2)*Rfqpara(1,nc-1)
                    oldA10 = Rfqpara(6,nc-1)
                    oldA0 = Rfqpara(7,nc-1)
                    oldzlen = Rfqzlen(nc-1)
                    
                    A10 = 0.5*(P1 + oldA10)
                    A0  = 0.5*(P2 + oldA0)
                    escale = 0.5*(escale + oldscl)
                    wk = 2.0*Pi/(zlen+oldzlen)
                endif
            else if (celltype.eq.3) then
                wk = 0.5*Pi/zlen
                T10 = m*m*bessi0(wk*ra) + bessi0(m*wk*ra)
                T30 = m*m*bessi0(3.0*wk*ra) + bessi0(3.0*m*wk*ra)
                alpha = bessi0(wk*ra)/bessi0(3.0*wk*ra)

                U10 = (m*m-1.0)/(T10 + alpha*T30/3.0)
                U30 = -alpha/3.0/U10
            else if (celltype.eq.4) then
                wk = Pi/zlen
                A30 = P5/bessi0(3.0*wk*r0)
                wk = 0.5*Pi/zlen
            endif

        endif

        zz = pos(3) - zedge

        r = sqrt(pos(1)*pos(1) + pos(2)*pos(2))
        theta = datan2(pos(1), pos(2))

        if (r < 1.0d-256) then
            r = 1.0d-256
            theta = 0.0
        endif

        if (celltype.eq.0) then ! two term  XY potential for standard acceleration cell
            if (mod(stdcnt,2).eq.1) then
                    zz = zz + zlen
            endif
        
            tmper = 2.0*A0*r/r0/r0*cos(2.0*theta) + A10*bessi1(wk*r)*wk*cos(wk*zz)
            tmpethe = -2.0*A0*r/r0/r0*sin(2.0*theta)
            tmpez = -A10*bessi0(wk*r)*wk*sin(wk*zz)

        else if (celltype.eq.1) then ! eight term  XY potential for standard acceleration cell
            if (mod(stdcnt,2).eq.1) then
                    zz = zz + zlen
            endif
            
            s2t = sin(2.0*theta)
            s4t = sin(4.0*theta)
            s6t = sin(6.0*theta)
            c2t = cos(2.0*theta)
            c4t = cos(4.0*theta)
            c6t = cos(6.0*theta)
            
            tmper = 2.0*A0*r/r0/r0*c2t + 6.0*A1*((r/r0)**5)/r0*c6t &
                 + A10*bessi1(wk*r)*wk*cos(wk*zz) + A12*(bessi3(2.0*wk*r) + bessi5(2.0*wk*r))&
                 *wk*c4t*cos(wk*zz) + A21*(bessi1(2.0*wk*r) + bessi3(2.0*wk*r))&
                 *wk*c2t*cos(2.0*wk*zz) + 1.5*A23*(bessi5(3.0*wk*r) + bessi7(3.0*wk*r))&
                 *wk*c6t*cos(2.0*wk*zz) + 3.0*A30*bessi1(3.0*wk*r)*wk*cos(3.0*wk*zz) &
                 + 1.5*A32*(bessi3(3.0*wk*r)+bessi5(3.0*wk*r))*wk*c4t*cos(3.0*wk*zz)

            tmpethe = -2.0*A0*r/r0/r0*s2t - 6.0*A1*((r/r0)**5)/r0*s6t &
                     + (-4.0*A12*bessi4(2.0*wk*r)*s4t*cos(wk*zz) &
                        -2.0*A21*bessi2(2.0*wk*r)*s2t*cos(2.0*wk*zz) &
                        -6.0*A23*bessi6(2.0*wk*r)*s6t*cos(2.0*wk*zz) &
                        -4.0*A32*bessi4(3.0*wk*r)*s4t*cos(3.0*wk*zz))/r 
            
            tmpez = - A10*bessi0(wk*r)*wk*sin(wk*zz) &
                    - A12*bessi4(2.0*wk*r)*wk*c4t*sin(wk*zz) &
                    - 2.0*A21*bessi2(2.0*wk*r)*wk*c2t*sin(2.0*wk*zz) &
                    - 2.0*A23*bessi6(3.0*wk*r)*wk*c6t*sin(2.0*wk*zz) &
                    - 3.0*A30*bessi0(3.0*wk*r)*wk*sin(3.0*wk*zz) &
                    - 3.0*A32*bessi4(3.0*wk*r)*wk*c4t*sin(3.0*wk*zz)
            
        else if (celltype.eq.2) then ! Radial matching
            wk = 0.5*Pi/zlen

            swkz = sin(wk*zz)
            swkz3 = sin(3.0*wk*zz)
            cwkz = cos(wk*zz)
            cwkz3 = cos(3.0*wk*zz)

            B1 = 0.5*wk*(bessi1(wk*r) + bessi3(wk*r))*swkz &
                -1.5/27.0*wk*(bessi1(3.0*wk*r) + bessi3(3.0*wk*r))*swkz3

            B2 = 0.5*wk*(bessi5(wk*r) + bessi7(wk*r))*swkz &
                -1.5/2187.0*wk*(bessi5(3.0*wk*r) + bessi7(3.0*wk*r))*swkz3

            tmper = P2*B1*cos(2.0*theta) + P4*B2*cos(6.0*theta)
            tmpethe = (-2.0*P2*(bessi2(wk*r)*swkz - 1.0/27.0*bessi2(3.0*wk*r)*swkz3)*sin(2.0*theta)&
                       -6.0*P4*(bessi6(wk*r)*swkz - 1.0/2187.0*bessi6(3.0*wk*r)*swkz3)*sin(6.0*theta))/r
            tmpez = P2*(bessi2(wk*r)*cwkz*wk - 1.0/9.0*bessi2(3.0*wk*r)*cwkz3*wk)*cos(2.0*theta)&
                    + P4*(bessi6(wk*r)*cwkz*wk - 1.0/729.0*bessi6(3.0*wk*r)*cwkz3*wk)*cos(6.0*theta)
                    
        else if (celltype.eq.3) then ! Exit transition cell (inverse)
            tmper = 2.0*r/r0/r0*cos(2.0*theta) + U10*wk*bessi1(wk*r)*cos(wk*zz) &
                  + 3.0*U30*wk*bessi1(3.0*wk*r)*cos(3.0*wk*zz)
            tmpethe = -2.0*r/r0/r0*sin(2.0*theta)
            tmpez = U10*bessi0(wk*r)*wk*sin(wk*zz) - U30*bessi0(3.0*wk*r)*3.0*wk*sin(3.0*wk*zz)

        else if (celltype.eq.4) then ! Exit transition cell (normal)
            tmper = 2.0*r/r0/r0*cos(2.0*theta) - P1*bessi1(wk*r)*wk*cos(wk*zz) &
                   -3.0*A30*bessi1(3.0*wk*r)*wk*cos(3.0*wk*zz)
            tmpethe = -2.0*r/r0/r0*sin(2.0*theta)
            tmpez = P1*bessi0(wk*r)*wk*sin(wk*zz) + A30*bessi0(3.0*wk*r)*3.0*wk*sin(3.0*wk*zz)

        else if (celltype.eq.5) then ! M cell
            wk = 0.0
            tmper = 2.0*P2*r/r0/r0*cos(2.0*theta) + 6.0*P4*((r/r0)**5)/r0*cos(6.0*theta)
            tmpethe =-2.0*P2*r/r0/r0*sin(2.0*theta) - 6.0*P4*((r/r0)**5)/r0*sin(6.0*theta)
            tmpez = 0.0

        else if (celltype.eq.6) then ! Fringe cell
            wk = 0.5*pi/zlen
            tmpez = 0.75*P1*wk*(sin(wk*zz)+sin(3.0*wk*zz))

            zz = zlen - zz

            swkz = sin(wk*zz)
            swkz3 = sin(3.0*wk*zz)

            B1 = 0.5*wk*(bessi1(wk*r) + bessi3(wk*r))*swkz &
                -1.5/27.0*wk*(bessi1(3.0*wk*r) + bessi3(3.0*wk*r))*swkz3

            B2 = 0.5*wk*(bessi5(wk*r) + bessi7(wk*r))*swkz &
                -1.5/2187.0*wk*(bessi5(3.0*wk*r) + bessi7(3.0*wk*r))*swkz3

            tmper = P2*B1*cos(2.0*theta) + P4*B2*cos(6.0*theta)

            tmpethe = (-2.0*P2*(bessi2(wk*r)*swkz - 1.0/27.0*bessi2(3.0*wk*r)*swkz3)*sin(2.0*theta) &
                      -6.0*P4*(bessi6(wk*r)*swkz - 1.0/2187.0*bessi6(3.0*wk*r)*swkz3)*sin(6.0*theta))/r

        else
        endif

        tmpey = tmper*cos(theta) - tmpethe*sin(theta)
        tmpex = tmper*sin(theta) + tmpethe*cos(theta)

        tmpex = 0.5*escale*tmpex
        tmpey = 0.5*escale*tmpey
        tmpez = 0.5*escale*tmpez

        tt = pos(4)
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)

        extfld(1) = tmpex*tmpsin
        extfld(2) = tmpey*tmpsin
        extfld(3) = tmpez*tmpsin
        extfld(4) = 0.0
        extfld(5) = 0.0
        extfld(6) = 0.0

        end subroutine getfld_RFQ

      end module RFQclass
