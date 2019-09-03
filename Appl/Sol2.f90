!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! Sol2class: Solenoid beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Authors: Ji Qiang, Robert Ryne, LANL, 7/13/01
! Description: This class defines the linear transfer map and field
!              for the Solenoid beam line elment.
! Comments:
!----------------------------------------------------------------
      module Sol2class
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 11
        type Sol2
          !Itype = 13
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : Bz0
          !      (3) : file ID
          !      (4) : radius
          !      (5) : By
          !      (6) : Bx
          !      (7) : x misalignment error
          !      (8) : y misalignment error
          !      (9) : rotation error x
          !      (10) : rotation error y
          !      (11) : rotation error z
        end type Sol2
        interface getparam_Sol2
          module procedure getparam1_Sol2,  &
                          getparam2_Sol2,   &
                          getparam3_Sol2
        end interface
        interface setparam_Sol2
          module procedure setparam1_Sol2,  &
                          setparam2_Sol2, setparam3_Sol2
        end interface
      contains
        subroutine construct_Sol2(this,numseg,nmpstp,type,blength)
        implicit none
        type (Sol2), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Sol2
   
        subroutine setparam1_Sol2(this,i,value)
        implicit none
        type (Sol2), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Sol2

        subroutine setparam2_Sol2(this,values)
        implicit none
        type (Sol2), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_Sol2

        subroutine setparam3_Sol2(this,numseg,nmpstp,type,blength)
        implicit none
        type (Sol2), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Sol2
   
        subroutine getparam1_Sol2(this,i,blparam) 
        implicit none 
        type (Sol2), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Sol2
  
        subroutine getparam2_Sol2(this,blparams)
        implicit none
        type (Sol2), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_Sol2

        subroutine getparam3_Sol2(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Sol2), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Sol2
       
        subroutine maplinear_Sol2(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (Sol2), intent(in) :: this
        double precision, dimension(22) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,qmcc
        integer :: mpstp

        xm = 0.0

        zedge = this%Param(1)
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
        y(19) = 1.0
        y(20) = 0.0
        y(21) = 0.0
        y(22) = 1.0
        h=tau/mpstp
        t0=t
        gi=-refpt(6)
        betai=sqrt(gi**2-1.)/gi
        ui=gi*betai
        squi=sqrt(ui)
        sqi3=sqrt(ui)**3
        dlti = 0.0
        thli = 0.0

        call rk6i_Sol2(h,mpstp,t0,y,22,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3

        dltf = 0.0
        thlf = 0.0

        xm(1,1)= y(19)*(y(7)*squi/squf+y(9)*squi/squf*dlti)
        xm(2,1)=y(19)*((y(8)-y(7)*dltf)*squi*squf+(y(10)-y(9)*dltf)*squi*squf*&
                 dlti)
        xm(1,2)= y(19)*y(9)/(squi*squf)
        xm(2,2)=y(19)*(y(10)-y(9)*dltf)*squf/squi
        xm(1,3)= y(20)*(y(7)*squi/squf+y(9)*squi/squf*dlti)
        xm(2,3)=y(20)*((y(8)-y(7)*dltf)*squi*squf+(y(10)-y(9)*dltf)*squi*squf*&
                 dlti)
        xm(1,4)= y(20)*y(9)/(squi*squf)
        xm(2,4)=y(20)*(y(10)-y(9)*dltf)*squf/squi
        xm(3,1)= y(21)*(y(11)*squi/squf+y(13)*squi/squf*dlti)
        xm(4,1)= &
        y(21)*((y(12)-y(11)*dltf)*squi*squf+(y(14)-y(13)*dltf)*squi*squf*dlti)
        xm(3,2)= y(21)*y(13)/(squi*squf)
        xm(4,2)=y(21)*(y(14)-y(13)*dltf)*squf/squi
        xm(3,3)= y(22)*(y(11)*squi/squf+y(13)*squi/squf*dlti)
        xm(4,3)= y(22)* &
        ((y(12)-y(11)*dltf)*squi*squf+(y(14)-y(13)*dltf)*squi*squf*dlti)
        xm(3,4)= y(22)*y(13)/(squi*squf)
        xm(4,4)=y(22)*(y(14)-y(13)*dltf)*squf/squi
        xm(5,5)= y(15)*sqi3/sqf3+y(17)*sqi3/sqf3*thli
        xm(6,5)= &
        (y(16)-y(15)*thlf)*sqi3*sqf3+(y(18)-y(17)*thlf)*sqi3*sqf3*thli
        xm(5,6)= y(17)/(sqi3*sqf3)
        xm(6,6)=(y(18)-y(17)*thlf)*sqf3/sqi3

        end subroutine maplinear_Sol2

        subroutine rk6i_Sol2(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (Sol2), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_Sol2(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_Sol2(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_Sol2(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_Sol2(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_Sol2(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_Sol2(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_Sol2(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_Sol2(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_Sol2

        subroutine intfunc1_Sol2(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (Sol2), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55
        double precision :: zedge,qmcc,brho,bgrad,C2alpha,bz
        integer :: my_rank, ierr

        zedge = this%Param(1)
        bz = this%Param(2)
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

        s11=((0.5*bz/brho)**2)*Scxl
        s33=((0.5*bz/brho)**2)*Scxl
        s55= 0.0

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
        C2alpha = -bz/brho/2 !rotation matrix of Solenoid
        f(19)=-C2alpha*y(21)
        f(20)=-C2alpha*y(22)
        f(21)=C2alpha*y(19)
        f(22)=C2alpha*y(20)
 
        end subroutine intfunc1_Sol2

        subroutine getBgradfld_Sol2(z,this,bz,bgrad)
        implicit none
        include 'mpif.h'
        type (Sol2), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: bz,bgrad

        !uniform bz field.
        bz = this%Param(2)
        bgrad = 0.0

        end subroutine getBgradfld_Sol2

        !get external field with displacement and rotation errors.
        subroutine  getflderr_Sol2(pos,extfld,this,dx,dy,anglex,angley,&
                                  anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Sol2), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge
        double precision:: zz,bgrad,bz,bx,by
        double precision, dimension(3) :: temp,tmp

        zedge = this%Param(1)
        bz = this%Param(2)
        by = this%Param(5)
        bx = this%Param(6)
        !move into the tank local coordinate.
!        zz=pos(3)-zedge
        !uniform bz field.

!        dx = this%Param(7)
!        dy = this%Param(8)
!        anglex = this%Param(9)
!        angley = this%Param(10)
!        anglez = this%Param(11)
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

        bgrad = 0.0
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = -bgrad/2*tmp(1)+bx
        extfld(5) = -bgrad/2*tmp(2)+by
        extfld(6) = bz

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        end subroutine getflderr_Sol2
        
        !get external field without displacement and rotation errors and
        !with fringe field of Solenoid. (f(z) = bz + bb*z)
        !here, the length of the solenoid has to the effective hard-edge length + 2*aperature 
        subroutine  getfld_Sol2(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Sol2), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge
        double precision:: zz,bgrad,bx,by,bz,bb,cc,zshift,x0,x02

        zedge = this%Param(1)
!        x0 = this%Param(4)
!        x0 = this%Length
        x0 = 0.
        x02= x0
!        x02 = 2*x0 !range of fringe field (x0 outside the length)
!        x02 = x0 !range of fringe field (x0 outside the length)
        !move into the tank local coordinate.
        zz=pos(3)-zedge
        !uniform bz field.
        bz = this%Param(2)
        by = this%Param(5)
        bx = this%Param(6)
        !bz = this%Param(2)*(1.0+x02/(this%Length-1.0*x02))!ikegami
!        bz = this%Param(2)*(this%Length/(this%Length-(4.0/3.0)*x02))**0.5
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0

        if((zz.ge.0.0).and.(zz.lt.x02)) then
          bb = bz/x02
          bgrad = bb
          extfld(4) = -bgrad/2*pos(1)+bx
          extfld(5) = -bgrad/2*pos(2)+by
          extfld(6) = bb*(zz-this%Length*0.5)
        else if((zz.ge.x02).and.(zz.le.this%Length-x02)) then
          extfld(4) = bx
          extfld(5) = by
          extfld(6) = bz
        else if((zz.gt.this%Length-x02).and.(zz.le.this%Length)) then
          bb = -bz/x02
          bgrad = bb
          zshift = this%Length - x02
          extfld(4) = -bgrad/2*pos(1)+bx
          extfld(5) = -bgrad/2*pos(2)+by
          extfld(6) = bz + bb*(zz-zshift)

        else
          !stop
          extfld(4) = 0.0
          extfld(5) = 0.0
          extfld(6) = 0.0
        endif
        
        end subroutine getfld_Sol2

        subroutine hedge_Sol2(Pts1,innp,gam0,drange)
        implicit none
        include 'mpif.h'
        double precision, pointer, dimension(:,:) :: Pts1
        integer, intent(in) :: innp
        double precision, intent(in) :: gam0
        double precision, dimension(8), intent(in) :: drange
        double precision :: gam, gambetz0, finv, b0
        double precision, dimension(6) :: ptarry
        integer :: ipt, flg

        b0 = drange(3)
        flg = nint(drange(4))

        do ipt = 1, innp
            gam = gam0 - Pts1(6,ipt)
            gambetz0 = sqrt(gam0*gam0-1.0d0)
            finv = 0.5*Pts1(7,ipt)*Clight/gambetz0*b0

            if (flg.eq.0) then
                ptarry(2) = Pts1(2,ipt)/gambetz0+finv*Pts1(3,ipt)*Scxl
                ptarry(4) = Pts1(4,ipt)/gambetz0-finv*Pts1(1,ipt)*Scxl
            else
                ptarry(2) = Pts1(2,ipt)/gambetz0-finv*Pts1(3,ipt)*Scxl
                ptarry(4) = Pts1(4,ipt)/gambetz0+finv*Pts1(1,ipt)*Scxl
            endif

            ptarry(1) = Pts1(1,ipt)
            ptarry(2) = ptarry(2)*gambetz0
            ptarry(3) = Pts1(3,ipt)
            ptarry(4) = ptarry(4)*gambetz0
            ptarry(5) = Pts1(5,ipt)
            ptarry(6) = Pts1(6,ipt)
            Pts1(1,ipt) = ptarry(1)
            Pts1(2,ipt) = ptarry(2)
            Pts1(3,ipt) = ptarry(3)
            Pts1(4,ipt) = ptarry(4)
            Pts1(5,ipt) = ptarry(5)
            Pts1(6,ipt) = ptarry(6)
        enddo

        end subroutine hedge_Sol2

      end module Sol2class
