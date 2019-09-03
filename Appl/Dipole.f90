!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! Dipoleclass: Dipole beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 2.0
! Author: Ji Qiang, LBNL, 11/22/04
! Description: This class defines the linear transfer map and field
!              for the Dipole beam line elment.
! Comments: the misalignment and rotation error do not work
!----------------------------------------------------------------
      module Dipoleclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 11
        type Dipole
          !Itype = 4
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : angle (x field strength)
          !      (3) : reference beta*gamma (y field strength)
          !      (4) : file ID
          !      (5) : radius
          !      (6) : front pole face angle1 (x misalignment error)
          !      (7) : back pole face angle2 (y misalignment error)
          !      (8) : front curvature (rotation error x)
          !      (9) : back curvature (rotation error y)
          !      (10) : fringe quadrupole term (rotation error z)
          !      (11) : main quadrupole term
        end type Dipole
        interface getparam_Dipole
          module procedure getparam1_Dipole,  &
                          getparam2_Dipole,   &
                          getparam3_Dipole
        end interface
        interface setparam_Dipole
          module procedure setparam1_Dipole,  &
                          setparam2_Dipole, setparam3_Dipole
        end interface
      contains
        subroutine construct_Dipole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Dipole), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Dipole
   
        subroutine setparam1_Dipole(this,i,value)
        implicit none
        type (Dipole), intent(out) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Dipole

        subroutine setparam2_Dipole(this,values)
        implicit none
        type (Dipole), intent(out) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_Dipole

        subroutine setparam3_Dipole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Dipole), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Dipole
   
        subroutine getparam1_Dipole(this,i,blparam) 
        implicit none 
        type (Dipole), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Dipole
  
        subroutine getparam2_Dipole(this,blparams)
        implicit none
        type (Dipole), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_Dipole

        subroutine getparam3_Dipole(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Dipole), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Dipole
       
!------------------------------------------------------------------------
!The linear map calculation for the dipole is not correct
        subroutine maplinear_Dipole(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (Dipole), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,escale,tfin
        integer  :: mpstp

        xm = 0.0

        zedge = this%Param(1)
        escale = this%Param(2)
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

        dlti=0.0
        thli=0.0

        call rk6i_Dipole(h,mpstp,t0,y,18,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3

        tfin = t + tau

        dltf = 0.0
        thlf = 0.0

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
        !print*,"xm: ",xm(1,1),xm(2,1),xm(1,2),xm(2,2),xm(3,3),xm(4,3),xm(4,4)

        end subroutine maplinear_Dipole

        subroutine rk6i_Dipole(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (Dipole), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_Dipole(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_Dipole

        subroutine intfunc1_Dipole(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (Dipole), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55
        double precision :: zedge,qmcc,brho,bgrad

        zedge = this%Param(1)
        bgrad = this%Param(2)
        qmcc = Bchg/Bmass

        f(1) = 0.0
        f(2) = 0.0
        f(3) = 0.0
        f(4) = 0.0

        ! synchronous particle:
        gamma0=-y(6)
        beta0=sqrt(gamma0**2-1.)/gamma0
        gbet=sqrt((gamma0-1.0)*(gamma0+1.0))

        f(5)=1.0/(beta0*Scxl)
        f(6)=0.0

        ! matrix elements
        brho=gbet/Clight/qmcc
        s11=bgrad/brho*Scxl
        s33=-bgrad/brho*Scxl
        s55=  0.0

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

        end subroutine intfunc1_Dipole
!------------------------------------------------------------------------

        !get external field with displacement and rotation errors.
        subroutine  getflderr_Dipole(pos,extfld,this,dx,dy,anglex,angley,&
                                     anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        type (Dipole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bx,by,zedge
        double precision, dimension(3) :: temp,tmp

        zedge = this%Param(1)
!        zz = pos(3)-zedge

!        dx = this%Param(6)
!        dy = this%Param(7)
!        anglex = this%Param(8)
!        angley = this%Param(9)
!        anglez = this%Param(10)

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

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = this%Param(2)
        extfld(5) = this%Param(3)
        extfld(6) = 0.0

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        end subroutine getflderr_Dipole
        
        !get external field without displacement and rotation errors.
        subroutine  getfld_Dipole(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Dipole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bx,by,zedge,z1,z2
        integer :: faceid

        faceid = Fcoef(1)+0.01 !switch id of pole face
        !Fcoef(2) is the gamma of the reference particle
        !Fcoef(3)-(6) are the coeficients for the geometry of bend
        extfld = 0.0
        !The pole face is characterized by z = a x + b
        !Fcoef(2) is the gamma of the reference particle
        !Fcoef(3)-(6) are the coeficients for the geometry of bend
        extfld = 0.0
        z1 = Fcoef(3)*pos(1) + Fcoef(4)
        z2 = Fcoef(5)*pos(1) + Fcoef(6)
        if((zz.ge.z1).and.(zz.le.z2)) then
            extfld(5) = this%Param(3)
        endif

        end subroutine getfld_Dipole

        !get external field without displacement and rotation errors.
        subroutine  getfldtmp_Dipole(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Dipole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bx,by,zedge,z1,z2
        integer :: faceid

        faceid = Fcoef(1)+0.01 !switch id of pole face
        !Fcoef(2) is the gamma of the reference particle
        !Fcoef(3)-(6) are the coeficients for the geometry of bend
        extfld = 0.0
        if(faceid.eq.1) then !parallel face
          z1 = Fcoef(3)*pos(3) + Fcoef(4)
          z2 = Fcoef(5)*pos(3) + Fcoef(6)
          if((pos(1).ge.z2).and.(pos(1).le.z1)) then
            extfld(5) = this%Param(3)
          endif 
        else if(faceid.eq.2) then !non-parallel face, bending toward x < 0
          z1 = Fcoef(3)*pos(3) + Fcoef(4)
          z2 = Fcoef(5)*pos(3) + Fcoef(6)
          if((pos(1).ge.z2).and.(pos(1).ge.z1)) then
            extfld(5) = this%Param(3)
          endif 
        else if(faceid.eq.3) then !non-parallel face, bending toward x > 0
          z1 = Fcoef(3)*pos(3) + Fcoef(4)
          z2 = Fcoef(5)*pos(3) + Fcoef(6)
          if((pos(1).le.z2).and.(pos(1).le.z1)) then
            extfld(5) = this%Param(3)
          endif
        else
          print*,"wrong input option!!!!"
          stop
        endif

        end subroutine getfldtmp_Dipole

        !get external field without displacement and rotation errors.
        subroutine  getfldold_Dipole(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Dipole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bx,by,zedge

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = this%Param(2)
        extfld(5) = this%Param(3)
        extfld(6) = 0.0

        end subroutine getfldold_Dipole

        !Transfer matrix follows the Transport: K. L. Brown, SLAC-75.
        subroutine Fpol_Dipole(h0,h1,tanphi,tanphib,k1,psi,ptarry1,ptarry2,ang)
        implicit none
        include 'mpif.h'
        double precision, dimension(6), intent(in) :: ptarry1
        double precision :: h0,h1,tanphi,tanphib,k1,psi,ang
        double precision, dimension(6), intent(out) :: ptarry2
        double precision :: tanphi2,secphi2,secphi,secphi3,secphib2

        tanphi2 = tanphi**2
        secphi2 = 1.0 + tanphi2
        secphi = 1.0/cos(ang)
        secphi3 = secphi2*secphi
        secphib2 = 1.+tanphib**2

        ptarry2 = 0.0d0

        ptarry2(1) = ptarry1(1) - 0.5d0*h0*tanphi2*ptarry1(1)**2 + &
                     0.5d0*h0*secphi2*ptarry1(3)**2
        ptarry2(2) = h0*tanphi*ptarry1(1)+ptarry1(2)+h0*tanphi2* &
                     ptarry1(1)*ptarry1(2)+(0.5d0*h0*h1*secphi3+ &
                     k1*tanphi)*ptarry1(1)**2-(0.5d0*h0*h1*secphi3+&
                     k1*tanphi-h0**2*tanphi*tanphi2-0.5*h0**2*tanphi)*&
                     ptarry1(3)**2 - h0*tanphi2*ptarry1(3)*ptarry1(4)-&
                     h0*tanphi*ptarry1(1)*ptarry1(6)
        ptarry2(3) = ptarry1(3) + h0*tanphi2*ptarry1(1)*ptarry1(3)
        ptarry2(4) = -h0*tanphib*ptarry1(3)+ptarry1(4)-h0*tanphi2*&
                     ptarry1(1)*ptarry1(4)-h0*secphi2*ptarry1(2)*&
                     ptarry1(3)-(h0*h1*secphi3+2*k1*tanphi)*ptarry1(1)*&
                     ptarry1(3)+(h0*tanphi-h0*psi*secphib2)*ptarry1(3)*&
                     ptarry1(6)
        ptarry2(5) = ptarry1(5)
        ptarry2(6) = ptarry1(6)

        end subroutine Fpol_Dipole

        subroutine Bpol_Dipole(h0,h1,tanphi,tanphib,k1,psi,ptarry1,ptarry2,ang)
        implicit none
        include 'mpif.h'
        double precision, dimension(6), intent(in) :: ptarry1
        double precision :: h0,h1,tanphi,tanphib,k1,psi,ang
        double precision, dimension(6), intent(out) :: ptarry2
        double precision :: tanphi2,secphi2,secphi,secphi3,secphib2

        tanphi2 = tanphi**2
        secphi2 = 1.0 + tanphi2
        secphi = 1.0/cos(ang)
        secphi3 = secphi2*secphi
        secphib2 = 1.+tanphib**2

        ptarry2 = 0.0d0

        ptarry2(1) = ptarry1(1) + 0.5d0*h0*tanphi2*ptarry1(1)**2 - &
                     0.5d0*h0*secphi2*ptarry1(3)**2
        ptarry2(2) = h0*tanphi*ptarry1(1)+ptarry1(2)-h0*tanphi2* &
                     ptarry1(1)*ptarry1(2)+(0.5d0*h0*h1*secphi3+ &
                     k1*tanphi-h0**2*tanphi2*tanphi/2)*ptarry1(1)**2-&
                     (0.5d0*h0*h1*secphi3+k1*tanphi+h0**2*tanphi*tanphi2/2)*&
                     ptarry1(3)**2 + h0*tanphi2*ptarry1(3)*ptarry1(4)-&
                     h0*tanphi*ptarry1(1)*ptarry1(6)
        ptarry2(3) = ptarry1(3) - h0*tanphi2*ptarry1(1)*ptarry1(3)
        ptarry2(4) = -h0*tanphib*ptarry1(3)+ptarry1(4)+h0*tanphi2*&
                     ptarry1(1)*ptarry1(4)+h0*secphi2*ptarry1(2)*&
                     ptarry1(3)-(h0*h1*secphi3+2*k1*tanphi-h0**2*secphi2*tanphi)*&
                     ptarry1(1)*ptarry1(3)+&
                     (h0*tanphi-h0*psi*secphib2)*ptarry1(3)*ptarry1(6)
        ptarry2(5) = ptarry1(5)
        ptarry2(6) = ptarry1(6)

        end subroutine Bpol_Dipole

        subroutine Sector_Dipole(len,beta,h0,k1,ptarry1,ptarry2,qmrel)
        implicit none
        include 'mpif.h'
        double precision, dimension(6), intent(in) :: ptarry1
        double precision :: h0,len,beta,k1
        double precision, dimension(6), intent(out) :: ptarry2
        double precision :: kx2,kx,cx,sx,dx,j1,ky2,ky,cy,sy,dy,&
                            gambet,gam2,qmrel

        gambet = beta/sqrt(1.0-beta**2)
        gam2 = 1.0/(1.0-beta**2)
        kx2 = h0**2 + k1
        if(kx2.gt.0.0) then
          kx = sqrt(kx2)
          cx = cos(kx*len)
          dx = (1.0-cx)/kx2
          sx = sin(kx*len)/kx
        else if(kx2.eq.0.0d0) then
          kx = sqrt(kx2)
          cx = cos(kx*len)
          dx = len**2/2
          sx = len
        else
          kx = sqrt(-kx2)
          cx = cosh(kx*len)
          dx = (1.0-cx)/kx2
          sx = sinh(kx)/kx
        endif
        j1 = (len-sx)/kx2
        ky2 = -k1
        if(ky2.gt.0.0) then
          ky = sqrt(ky2)
          cy = cos(ky*len)
          dy = (1.0-cy)/ky2
          sy = sin(ky*len)/ky
        else if(ky2.eq.0.0d0) then
          ky = sqrt(ky2)
          cy = cos(ky*len)
          dy = len**2/2
          sy = len
        else
          ky = sqrt(-ky2)
          cy = cosh(ky*len)
          dy = (1.0-cy)/ky2
          sy = sinh(ky)/ky
        endif

        ptarry2 = 0.0d0
        ptarry2(1) = cx*ptarry1(1)+sx*ptarry1(2)+&
                     h0*dx*ptarry1(6) 
        ptarry2(2) = -kx2*sx*ptarry1(1)+cx*ptarry1(2)+&
                     h0*sx*ptarry1(6) 
        ptarry2(3) = cy*ptarry1(3)+sy*ptarry1(4)
        ptarry2(4) = -ky2*sy*ptarry1(3)+cy*ptarry1(4)
        ! (5) is defined as -v dt = -c beta dt
        ! see p.147 of F. Iselin paper
        ptarry2(5) = -h0*sx*ptarry1(1)-h0*dx*ptarry1(2)+&
                     ptarry1(5) - h0**2*j1*ptarry1(6) + len/gam2*&
                     (ptarry1(6)+qmrel)
        ! Transport defines (5) as path differnce which is v dt
        !        ptarry2(5) = h0*sx*ptarry1(1)+h0*dx*ptarry1(2)+&
        !                     ptarry1(5) + h0**2*j1*ptarry1(6)
        ptarry2(6) = ptarry1(6)
        end subroutine Sector_Dipole

        
        !TRIUMF Beam Physics Note, TRI-DN-05-7, Electrostatic Bend Optics
        subroutine ESector_Dipole(lng,beta,h0,k1,ptarry1,ptarry2,flgsp)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: h0,lng,beta,k1
        integer, intent(in) :: flgsp
        double precision, dimension(6), intent(in) :: ptarry1
        double precision, dimension(6), intent(out) :: ptarry2
        double precision :: kx2, gambet, gam2, zeta, eta, btmp, aa, theta,&
                       coszt, sinzt, sinet, coset, sinetoe

        gambet = beta/sqrt(1.0-beta**2)
        gam2 = 1.0/(1.0-beta**2)
        kx2 = h0**2 + k1
        
        btmp = 0.0
        aa = 1./h0      !rho[m] h0= 1/rho
        theta = lng/aa  !angle[rad]
        
        if (flgsp.eq.0) then
            zeta = sqrt(2.0)
            eta = 0.0
            sinet = 0.0
            sinetoe = 1.0
        else if (flgsp.eq.1) then
            zeta = 1.0
            eta = 1.0
            sinet = sin(eta*theta)
            sinetoe = sinet
        endif
        
        coszt = cos(zeta*theta)
        sinzt = sin(zeta*theta)
        coset = cos(eta*theta)
        
        ptarry2 = 0.0d0
        ptarry2(1) = coszt*ptarry1(1) + aa/zeta*sinzt*ptarry1(2)&
                   + (2.-btmp*btmp)/zeta/zeta*aa*(1.-coszt)*ptarry1(6)
        ptarry2(2) = -zeta/aa*sinzt*ptarry1(1) + coszt*ptarry1(2)&
                    + (2.-btmp*btmp)/zeta*sinzt*ptarry1(6) 
        ptarry2(3) = coset*ptarry1(3)+ aa*sinetoe*ptarry1(4)
        ptarry2(4) = -eta/aa*sinet*ptarry1(3)+coset*ptarry1(4)
        ptarry2(5) = -(2.-btmp*btmp)/zeta*sinzt*ptarry1(1)&
                     -(2.-btmp*btmp)/zeta/zeta*aa*(1.-coszt)*ptarry1(2)&
                     + ptarry1(5) + aa*theta*(1./gam2-(((2.-btmp*btmp)/zeta)**2)&
                     *(1.-sinzt/zeta/theta))*ptarry1(6)
        ptarry2(6) = ptarry1(6)

        end subroutine ESector_Dipole
        
      end module Dipoleclass
