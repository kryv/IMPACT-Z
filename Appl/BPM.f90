!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! BPMclass: Beam position monitor class in Lattice module of APPLICATION 
!           layer.
! Version: 2.0
! Author: Ji Qiang, LBNL, 2/21/06
! Description: This class defines the different beam diagnostics at given
!              beam position.
! Comments:
!  1) Itype = -1, shift the transverse centroid position to 0.
!  2) Itype = -2, shift the transverse centroid position and angle to 0.
!                 (this one not work yet due to conflict of definition)
!  3) Itype = -10, mismatch the beam distribution by the amount given in
!                  Param(3) - Param(10).
!  4) Itype = -21, shift the beam centroid in 6D phase space by the amount
!                  given in Param(3) - Param(10).
!  5) Itype = -22, rotate in x-y and px-py plane by param(3) and param(4) degrees.
!QZ
!QZ	05/18/2007
!QZ	modified Itype = -21 by taking account of velocity factor
!QZ add Itype = -25 same amount 6D shifts for all particles (no magnetic effect) 
!QZ	Itype = -25 is for misalignment purpose --> shift beam instead of element
!QZ
!QZ	07/23/2007
!QZ add Itype = -26 to fulfil the function of "shift both the transverse centroid
!QZ positions and angles to 0, and to resolve the conflict with Itype = -2
!QZ
!----------------------------------------------------------------
      module BPMclass
        use PhysConstclass
        use Dataclass
        use MTrndclass
        integer, private, parameter :: Nparam = 8
        character(256) :: strip_dir = ''
        integer :: iostrpfile
        type BPM
          !Itype < 0
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : radius
          !      (3) : xmax
          !      (4) : pxmax
          !      (5) : ymax
          !      (6) : pymax
          !      (7) : zmax
          !      (8) : pzmax
        end type BPM
        interface getparam_BPM
          module procedure getparam1_BPM,  &
                          getparam2_BPM,   &
                          getparam3_BPM
        end interface
        interface setparam_BPM
          module procedure setparam1_BPM,  &
                           setparam2_BPM, setparam3_BPM
        end interface
      contains
        subroutine construct_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_BPM
   
        subroutine setparam1_BPM(this,i,value)
        implicit none
        type (BPM), intent(out) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_BPM

        subroutine setparam2_BPM(this,values)
        implicit none
        type (BPM), intent(out) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_BPM

        subroutine setparam3_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_BPM
   
        subroutine getparam1_BPM(this,i,blparam) 
        implicit none 
        type (BPM), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_BPM
  
        subroutine getparam2_BPM(this,blparams)
        implicit none
        type (BPM), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_BPM

        subroutine getparam3_BPM(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (BPM), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_BPM

        subroutine shift_BPM(Pts1,itype,innp,nptot)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: itype,innp,nptot
        double precision, pointer, dimension(:,:) :: Pts1
        double precision:: x0lc,px0lc,y0lc,py0lc
        double precision, dimension(4) :: tmplc,tmpgl
        integer :: i,j,ierr

        tmplc = 0.0
        tmpgl = 0.0
!QZ		if(itype.eq.(-2)) then
!QZ 07/23/07
        if(itype.eq.(-26)) then
          x0lc = 0.0
          px0lc = 0.0
          y0lc = 0.0
          py0lc = 0.0
          do i = 1, innp
            x0lc = x0lc + Pts1(1,i)
            px0lc = px0lc + Pts1(2,i)
            y0lc = y0lc + Pts1(3,i)
            py0lc = py0lc + Pts1(4,i)
          enddo

          tmplc(1) = x0lc
          tmplc(2) = px0lc
          tmplc(3) = y0lc
          tmplc(4) = py0lc
        
          call MPI_ALLREDUCE(tmplc,tmpgl,4,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

          tmpgl(1) = tmpgl(1)/nptot
          tmpgl(2) = tmpgl(2)/nptot
          tmpgl(3) = tmpgl(3)/nptot
          tmpgl(4) = tmpgl(4)/nptot

          do i = 1, innp
            Pts1(1,i) = Pts1(1,i) - tmpgl(1)
            Pts1(2,i) = Pts1(2,i) - tmpgl(2)
            Pts1(3,i) = Pts1(3,i) - tmpgl(3)
            Pts1(4,i) = Pts1(4,i) - tmpgl(4)
          enddo
        else if(itype.eq.(-1)) then
          x0lc = 0.0
          y0lc = 0.0
          do i = 1, innp
            x0lc = x0lc + Pts1(1,i)
            y0lc = y0lc + Pts1(3,i)
          enddo

          tmplc(1) = x0lc
          tmplc(3) = y0lc
        
          call MPI_ALLREDUCE(tmplc,tmpgl,4,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

          tmpgl(1) = tmpgl(1)/nptot
          tmpgl(3) = tmpgl(3)/nptot

          do i = 1, innp
            Pts1(1,i) = Pts1(1,i) - tmpgl(1)
            Pts1(3,i) = Pts1(3,i) - tmpgl(3)
          enddo
        else
        endif

        end subroutine shift_BPM

        !mismatch the beam at given location.
        !Here, the storage Param(3:8) is used to store the mismatch factors
        subroutine scale_BPM(Pts1,innp,xmis,pxmis,ymis,pymis,zmis,pzmis)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: xmis,pxmis,ymis,pymis,zmis,pzmis
        integer :: i
 
        do i = 1, innp
            Pts1(1,i) = xmis*Pts1(1,i)
            Pts1(2,i) = pxmis*Pts1(2,i)
            Pts1(3,i) = ymis*Pts1(3,i)
            Pts1(4,i) = pymis*Pts1(4,i)
            Pts1(5,i) = zmis*Pts1(5,i)
            Pts1(6,i) = pzmis*Pts1(6,i)
        enddo
 
        end subroutine scale_BPM

        !mismatch the beam at given location.
        !Here, the storage Param(3:8) is used to store the mismatch factors
        !rotate the beam in the x-y plane, px-py plane by a t1,t2
        !degree.t1 = Param(3),t2 = Param(4)
        subroutine xypxpyrot_BPM(Pts1,innp,t1,t2)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: t1,t2
        real*8 :: theta1,theta2,tmp1,tmp2
        integer :: i
 
        theta1 = t1*2*asin(1.0)/180
        theta2 = t2*2*asin(1.0)/180

        do i = 1, innp
            tmp1 = Pts1(1,i)*cos(theta1)+Pts1(3,i)*sin(theta1)
            tmp2 = -Pts1(1,i)*sin(theta1)+Pts1(3,i)*cos(theta1)
            Pts1(1,i) = tmp1
            Pts1(3,i) = tmp2
            tmp1 = Pts1(2,i)*cos(theta2)+Pts1(4,i)*sin(theta2)
            tmp2 = -Pts1(2,i)*sin(theta2)+Pts1(4,i)*cos(theta2)
            Pts1(2,i) = tmp1
            Pts1(4,i) = tmp2
        enddo
 
        end subroutine xypxpyrot_BPM

        !! add 20170322 Yue Hao ------    
        subroutine yrotation_BPM(Pts1,innp,t1,gam)

        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: t1
        double precision, intent(in) :: gam
        real*8 :: theta1,pzinmc,tan_pz,tan_theta1,cos_theta1,sin_theta1
        integer :: i
        theta1 = t1*Pi/180.0
        tan_theta1=tan(theta1)
        cos_theta1=cos(theta1)
        sin_theta1=sin(theta1)        
        do i = 1, innp
            pzinmc = sqrt(-Pts1(2,i)**2 - Pts1(4,i)**2 + (gam-Pts1(6,i))**2-1.0)
            tan_pz=tan_theta1/pzinmc
            Pts1(3,i) = Pts1(3,i) + Pts1(4,i)*Pts1(1,i)*tan_pz/(1-pts1(2,i)*tan_pz)
            Pts1(5,i) = Pts1(5,i) + (gam-Pts1(6,i))*Pts1(1,i)*tan_pz/(1-Pts1(2,i)*tan_pz)
            Pts1(1,i) = Pts1(1,i)/cos_theta1/(1-Pts1(2,i)*tan_pz)
            Pts1(2,i) = Pts1(2,i)*cos_theta1+pzinmc*sin_theta1            
        enddo
        
        end subroutine yrotation_BPM
        
        !J. Q. modified for multiple chage state.
        !originally added by M. I. of KEK
        !shift the beam centroid in the 6D phase space.
        !This element can be used to model steering magnet etc.
        !Here, the storage Param(3:8) is used to store the amount of shift.
        !drange(3); shift in x (m)
        !drange(4); shift in Px (rad)
        !drange(5); shift in y (m)
        !drange(6); shift in Py (rad)
        !drange(7); shift in z (deg)
        !drange(8); shift in Pz (MeV)
        !gam; relativistic gamma factor for design particle
        !mass; mass in eV/c^2
        subroutine kick_BPM(Pts1,innp,xshift,pxshift,yshift,pyshift,zshift,&
                   pzshift,gam,mass,chge)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: xshift,pxshift,yshift,pyshift,&
                                      zshift,pzshift
        double precision, intent(in) :: gam,mass,chge
        integer :: i
        double precision :: gambetz,tmp,chgmass
        !q/m of reference particle
        chgmass = chge/mass
        do i = 1, innp
            gambetz = sqrt(-Pts1(2,i)**2-Pts1(4,i)**2+(gam-Pts1(6,i))**2-1.0)
            tmp = Pts1(7,i)/chgmass
            Pts1(1,i) = Pts1(1,i)+xshift/Scxl
!QZ            Pts1(2,i) = Pts1(2,i)+pxshift*gambetz*tmp
!QZ
!QZ	kick of angle proportional to Q/A and also inversely to beta-gamma (velocity)
!QZ	x' of IMPACT = x' in rad * (beta-gamma)zi 
!QZ	so x'i of IMPACT = x'o in rad * (beta-gamma)zi * (beta-gamma)zo/(beta-gamma)zi * (Q/A)i/(Q/A)o
!QZ														velocity factor			      charge factor
			Pts1(2,i) = Pts1(2,i)+pxshift*sqrt(gam*gam-1.0)*tmp
            Pts1(3,i) = Pts1(3,i)+yshift/Scxl
!QZ            Pts1(4,i) = Pts1(4,i)+pyshift*gambetz*tmp
            Pts1(4,i) = Pts1(4,i)+pyshift*sqrt(gam*gam-1.0)*tmp
			Pts1(5,i) = Pts1(5,i)+zshift/Rad2deg
            Pts1(6,i) = Pts1(6,i)+pzshift*1.0e6/mass*tmp
        enddo
 
        end subroutine kick_BPM

!QZ		05/18/2007
!QZ		based on kick_BPM, shift same amount of value for all particles 
!QZ     to simulate misalignment of element.
        !Here, the storage Param(3:8) is used to store the amount of shift.
        !drange(3); shift in x (m)
        !drange(4); shift in Px (rad)
        !drange(5); shift in y (m)
        !drange(6); shift in Py (rad)
        !drange(7); shift in z (deg)
        !drange(8); shift in Pz (MeV)
        !gam; relativistic gamma factor for design particle
        !mass; mass in eV/c^2
        subroutine qzmis_BPM(Pts1,innp,qzxshift,qzpxshift,qzyshift,qzpyshift,&
                   qzzshift,qzpzshift,gam,mass,chge)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: qzxshift,qzpxshift,qzyshift,&
                                      qzpyshift,qzzshift,qzpzshift
        double precision, intent(in) :: gam,mass,chge
        integer :: i
        double precision :: gambetz,tmp,chgmass
        !q/m of reference particle
        chgmass = chge/mass
        do i = 1, innp
            gambetz = sqrt(-Pts1(2,i)**2-Pts1(4,i)**2+(gam-Pts1(6,i))**2-1.0)
            tmp = Pts1(7,i)/chgmass
            Pts1(1,i) = Pts1(1,i)+qzxshift/Scxl
            Pts1(2,i) = Pts1(2,i)+qzpxshift*gambetz
            Pts1(3,i) = Pts1(3,i)+qzyshift/Scxl
            Pts1(4,i) = Pts1(4,i)+qzpyshift*gambetz
            Pts1(5,i) = Pts1(5,i)+qzzshift/Rad2deg
            Pts1(6,i) = Pts1(6,i)+qzpzshift*1.0e6/mass
        enddo
 
        end subroutine qzmis_BPM


        subroutine MStripF_BPM(btd,tcord,nfile,imaa,npoint,rwkinf,rwkinq,&
                   nchrg,nptlist,currlist,qmcclist)
!
! MSU Stripper Routine
!
! Transforms particles through the stripper foil.
!
! In the "foil" subroutine the parameters are:
!
! amf(MCHARS) : A - mass of the ions after the foil
! zmf(MCHARS) : Z - charge of the ions after the foil
! fmf(MCHARS) : Fractions of different ion species after the foil
! f(ang,E)=A*ang*exp(-(ang/sigma_ang)^power_ang)*exp((E-Eo(ang))/sigma)^2),
! where:
!   ang is the angular scattering of the particle,
!   E is the particle energy after the stripper.
!
! The angular dependance assuming energy straggling is described by:
!   E-Eo(ang)=dE(ang)=B*exp((C*ang)^2)-B
!   A is a normalization constant,
!   sigma_ang, power_ang, dEo, B, C, sigma -- an empirical formula
!   Parameters are defined from SRIM simulation results for two cases.
!
! alf(1) : thickness of the foil [um]
!        (0 or negative means "use old LANA model" -- for checking basically)
! alf(2) : thickness variation [%]
! alf(3) : the reference SRIM calculated point for particles energy [MeV/u]
! alf(4) : -"- for thickness [um]
! clf(1) : B (energy loss angular dependance)
! clf(2) : C (energy loss angular dependance)
! clf(3) : not used
! All other parameters are linearly interpolated in 2-parameter space --
!
! input energy and thickness:
!  slf(1) : Eo -- energy loss due to ionization -- constant
!  slf(2) : Eo -- energy loss due to ionization -- energy dependance coefficient
!  slf(3) : Eo -- energy loss due to ionization -- thickness dependance coefficient
!  dlf(1) : sigma -- energy struggling -- constant
!  dlf(2) : sigma -- energy struggling -- energy dependance coefficient
!  dlf(3) : sigma -- energy struggling -- thickness dependance coefficient
!  drf(1) : sigma_ang -- angular scattering -- constant
!  drf(2) : sigma_ang -- angular scattering -- energy dependance coefficient
!  drf(3) : sigma_ang -- angular scattering -- thickness dependance coefficient
!  prf(1) : power_ang -- angular scattering -- constant
!  prf(2) : power_ang -- angular scattering -- energy dependance coefficient
!  prf(3) : power_ang -- angular scattering -- thickness dependance coefficient
!
! In the input data:
!  amfoil -> amf
!  zmfoil -> zmf
!  fmfoil -> fmf
!  alfoil -> alf
!  clfoil -> clf
!  slfoil -> slf
!  dlfoil -> dlf
!  drfoil -> drf
!  prfoil -> prf
!
! Stripper 1: (calculated for carbon foil, thickness 1.78 um ~0.4 mg/cm^2)
!  amfoil= 238.,238.,238.,238.,238.,
!  zmfoil=  73., 71., 72., 74., 75.,
!  fmfoil=  .2, .2, .2, .2, .2,
!  alfoil= 1.78, 5., 12., 1.78,
!  clfoil= 0.78712, 0.0169905,
!  slfoil= 11.78981, 1.004200, -0.1182028,
!  dlfoil= 2.289452d-3, -1.202082d-5, 6.991233d-4,
!  drfoil= 0.4218618, -3.98d-2, 1.724786d-1,
!  prfoil= 1.530484, -8.8d-3, 9.993088d-2,
!
! Stripper 2: (calculated for carbon foil, thickness 64.35 um ~15 mg/cm^2)
!  amfoil= 238.,238.,238.,
!  zmfoil=  88., 87., 89.,
!  fmfoil=  .34, .33, .33,
!  alfoil= 64.35, 5., 90., 64.35,
!  clfoil= 3.1033, 0.0221929,
!  slfoil= 86.93013, 1.024320, -4.841261d-2,
!  dlfoil= 0.01722606, -3.971112d-5, 1.800033d-4,
!  drfoil= 0.4797693, -7.078d-3, 4.236434d-3,
!  prfoil= 1.700679, -1.04d-3, 1.108175d-3,
!
      implicit real (KIND=8) (a-h,o-z)
      implicit INTEGER (KIND=4) (i-n)
      include 'mpif.h'
!      parameter (imaa=100000)
!      common /blcm1/rwkinf,rwkinq
!      COMMON/temp/icolch(imaa)
!      common /bcom/beta(imaa),wbeam
!      common /constcmn/pi,twopi,piovr2,piovr180,brhof,erest,clight
!      common ncell,nt,newtk,nstart,nstop,ntlast,nctotl,irun,nn,&
!       freq,wavel,tzero,pzero,btazro,qstate,tank(30,6),rfqtab(900), &
!       optcon(30),sce(10),vv(900),errl(10),error,xibeam,elimit,ntr1,&
!       ntr2,trans1(20,80),trans2(20,40)
!      common /reference/bfreq,qref,Wref,phiref,betar,bgr   !CKRC added
!      common /charg/ qstate1, current1, qstate2, current2, np1, np2  !CKRC
!      common /dist/ tdist,loc,nelem !accumulate distances in meters    !CKRC
!      double Precision  s_ran,s_gauss
!      common  /data_random/ it
      dimension amf(10),zmf(10),fmf(10), &
                alf(4),clf(2),slf(3),dlf(3),drf(3),prf(3)
      double precision, dimension(8,imaa) :: tcord
      integer, dimension(imaa) :: icolch
      integer :: nchrg
      integer, dimension(100) :: nptlist 
      double precision, dimension(100) :: qmcclist,currlist
      integer :: nptot
      double precision :: currtot
      integer :: myid,ierr
      character(256) ::fileid

      !set random seed on each processor
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      
      nptot = 0
      currtot = 0.0
      do i = 1, nchrg
        nptot = nptot + nptlist(i)
        currtot = currtot + currlist(i)
      enddo

      aunit=931.49432d6
      betatr=btd

      amf = 0.0
      zmf = 0.0
      fmf = 0.0

! readin the amf,zmf,fmf,alf,clf,slf,dlf,drf,prf. (JQ)
!     nfile = 2
!     print*,"before readin:"
     write(fileid,*), nfile
     open(nfile, file=trim(strip_dir)//'fort.'//trim(adjustl(fileid)),status='old',iostat=iostrpfile)
     if (iostrpfile.ne.0) return
     
     read(nfile,*)mchars
     read(nfile,*)(amf(i),i=1,mchars)
!     print*,"amf: ",amf
     read(nfile,*)(zmf(i),i=1,mchars)
!     print*,"zmf: ",zmf
     read(nfile,*)(fmf(i),i=1,mchars)
     read(nfile,*)(alf(i),i=1,4)
     read(nfile,*)(clf(i),i=1,2)
     read(nfile,*)(slf(i),i=1,3)
     read(nfile,*)(dlf(i),i=1,3)
     read(nfile,*)(drf(i),i=1,3)
     read(nfile,*)(prf(i),i=1,3)

     close(nfile)

!
! update the list
     nchrg = mchars
     do i = 1, nchrg
       nptlist(i) = nptot*fmf(i) + 1
       currlist(i) = currtot*fmf(i)
       qmcclist(i) = zmf(i)/amf(i)
     enddo

!
! Added Foil2 here
!
         rwkinq = ( 1.d0 / sqrt( 1.d0 - betatr**2 ) - 1.d0 )*aunit
         nchars= 0
         fmftot= 0.d0
         do ich=1,MCHARS
             if (fmf(ich).eq.0.) goto 5
             nchars= nchars + 1
             fmftot= fmftot + fmf(ich)
             if (fmftot.ge.1.d0) goto 5
         enddo
 5       continue
!        And shift the reference beta (energy) keeping dW(eV/u) for all particles
         do i=1,npoint
            tcord(7,i) = tcord(7,i)*rwkinq !*.01
         enddo
         rwkinf = rwkinq
         if ( alf(1) .gt. 0.d0 ) then
            rwkinq = rwkinq*1.d-6                          ! [eV/u] -> [MeV/u]
            rwkinq = slf(1) + slf(2)*( rwkinq    - alf(3) ) &
                               + slf(3)*( alf(1) - alf(4) )
            rwkinq = rwkinq*1.d6                           ! back [MeV/u] -> [eV/u]
         else
            rwkinq = rwkinq * ( 1.d0 - slf(1)*.01d0 )   ! dE in %
         endif
         betatr = sqrt( 1.d0 - 1.d0/( 1.d0 + rwkinq/aunit )**2 )
         do i=1,npoint
            tcord(7,i) = tcord(7,i)/rwkinq !*100.
         enddo
!
!  End of Foil2
!

      do i=1,npoint
         ach = fmftot*mtrnd()
         fmftit= 0.d0
         do ich=1,nchars
            if (fmftit.le.ach.and.ach.lt.fmftit+fmf(ich)) then
               icolch(i)= ich - 1
               goto 1
            endif
            fmftit= fmftit + fmf(ich)
         enddo
 1       continue
      enddo
      if ( alf(1) .le. 0.d0 ) then
         drfa = drf(1)                                ! [mrad]
         slfa = slf(1)                                ! dE/E[%]
         dlfa = dlf(1)                                ! dE/E[%]
      endif
!
      do i=1,npoint
         tcord(8,i) = zmf(1+icolch(i)) / amf(1+icolch(i))
!
!        "slfa" shift is taken care of in "RdCav" by shifting the reference
!        "slfb" -- from the correlation of the energy loss and angular straggling
!        Avoid "energy gain" in foil (gaussian distr. positive side)
!
         if ( alf(1) .gt. 0.d0 ) then                 ! f(x)=x*exp(-(x/aw)^ae)
            bza = rwkinf + tcord(7,i)*.01*rwkinq      ! dW/W[%] -> W_before[eV/u]
            bza = bza*1.d-6                           ! [eV/u] -> [MeV/u]
            alfa = alf(1)*(1.+alf(2)*.01*(1.-2.*mtrnd())) ! thickness spread [mkm,%]
            drfa = drf(1) + drf(2)*( bza  - alf(3) )  & 
                          + drf(3)*( alfa - alf(4) )  ! [mrad] from thickness [mkm]
            prfa = prf(1) + prf(2)*( bza  - alf(3) )  &
                          + prf(3)*( alfa - alf(4) )  ! [unitless] from thickness [mkm]
            slfa = slf(1) + slf(2)*( bza  - alf(3) )  & 
                          + slf(3)*( alfa - alf(4) )  ! [MeV/u]from thickness [mkm]
            slfa = bza - slfa                         ! dW[MeV/u] (down)
            dlfa = dlf(1) + dlf(2)*( bza  - alf(3) )  & 
                          + dlf(3)*( alfa - alf(4) )  ! [MeV/u]from thickness [mkm]

            dbr = 100.d0                              ! just to enter the loop
            do while ( abs(dbr) .gt. 30.d0 )
               dbr = s_invgammap(2.d0/prfa)           ! "a"=2/ae, "A"=1/(aw^ae)
            enddo
            dbr = (dbr*(drfa**prfa))**(1.d0/prfa)     ! (A*x^ae) -> (x) [mrad]
            slfb = clf(1)*(exp((clf(2)*dbr)**2)-1.)   ! angle-energy corr.-dE[MeV/u]

!
!           "sigma" here is like in "N(x0,sigma)=exp(-((x-x0)/(2*sigma))^2)
!
            if ( dlfa .gt. 0.d0 ) then
               dbz = 1.1*(slfa+slfb)                  ! just to enter the loop
               do while ( dbz.gt.(slfa+slfb) .or. abs(dbz).gt.3.d0*dlfa)
                  dbz = mtrndn()*dlfa                ! *sigma -> N(0,sigma)
               enddo
               bza = bza + dbz - (slfa+slfb)
            elseif ( (slfa+slfb) .gt. 0.d0 ) then
               bza = bza - (slfa+slfb)
            endif
            tcord(7,i) = ( bza*1.d6/rwkinq - 1. )*100.     ! W[MeV/u]->dW[eV/u]->dW/W[%]
         else !if ( alf(1) .le. 0.d0 ) then           ! Gauss
            if ( drfa .gt. 0.d0 ) then
               dbr = 100.d0                           ! just to enter the loop
               do while ( abs(dbr) .gt. 3.d0 )
                  dbr = mtrndn()                     ! N(0,1)
               enddo
               dbr = dbr*drfa                         ! *sigma -> N(0,sigma) [mrad]
            else
               dbr = 0.d0
            endif
            slfb = 0.d0
            if ( dlfa .gt. 0.d0 ) then
               dbz = 1.1*(slfa+slfb)                  ! just to enter the loop
               do while ( dbz.gt.(slfa+slfb) .or. abs(dbz).gt.3.d0*dlfa)
                  dbz = mtrndn()*dlfa                ! *sigma -> N(0,sigma)
               enddo
               tcord(7,i) = tcord(7,i) + dbz - slfb
            elseif ( slfb .gt. 0.d0 ) then
               tcord(7,i) = tcord(7,i) - slfb
            endif
         endif
!
         dba = 2.d0*PI*mtrnd()
         tcord(5,i) = tcord(5,i) + dbr*cos(dba)
         tcord(6,i) = tcord(6,i) + dbr*sin(dba)
      enddo

!      close(ifile)
!      rewind(nfile)

! 
! update list

      return
      end subroutine MStripF_BPM

!=======================================================================
!
!      double precision function s_ran()
!      implicit double precision(a-h,o-z)
!
! Generate uniform random number (0,1).
! 
!      double precision ran1
!      external         ran1
!      integer               it
!      common  /data_random/ it
!
!      s_ran = ran1(it)
!      return
!      end function s_ran
!
!=======================================================================
!
!      double precision function s_gauss()
!      implicit double precision(a-h,o-z)
!
! Generate Gaussian pair of random numbers.
!
!      double precision avr,sig,x,y,x2
!      double precision x,y,x2
!
!      integer               it
!      common  /data_random/ it
!      logical          iflg
!      data             iflg /.true./
!      save iflg, x2
!
!      if ( iflg ) then
!         call rannor( it, x, y )
!         x2 = y
!         s_gauss = x
!         iflg=.false.
!      else
!         s_gauss = x2
!         iflg=.true.
!      endif
!      return
!      end function s_gauss
!
!=======================================================================

      DOUBLE PRECISION FUNCTION gammaln( xvar )
      implicit double precision(a-h,o-z)
      INTEGER j
      DOUBLE PRECISION xvar, ser, stp, tmp, x, y, cof(6)
      SAVE stp, cof
      DATA stp, cof / 2.5066282746310005d+0 &
                   ,76.1800917294714600d+0,-86.505320329416770d+0 &
                   ,24.0140982408309100d+0, -1.231739572450155d+0 &
                   , 0.1208650973866179d-2, -0.539523938495300d-5/

      x = xvar
      y = x
      tmp = x + 5.5d0
      tmp = ( x + 0.5d0 )*log(tmp) - tmp
      ser = 1.000000000190015d+0
      do j=1,6
        y = y + 1.d0
        ser = ser + cof(j)/y
      enddo
      gammaln = tmp + log( stp * ser/x )
      return
      END FUNCTION gammaln

!=======================================================================

      DOUBLE PRECISION FUNCTION gammaser( a, x, gln, ifail )
      implicit double precision(a-h,o-z)
      INTEGER ITMAX,ifail
      DOUBLE PRECISION a, x, gln, EPS
      PARAMETER ( ITMAX=100, EPS=3.d-16 )
!U    USES gammaln
      INTEGER n
      DOUBLE PRECISION ap, del, sum
!
      ifail=0
      gln = gammaln(a)
      if ( x .le. 0.d0 ) then
        ifail=1
        if ( x .lt. 0.d0 ) then
         ifail=-1
        endif
        gammaser = 0.d0
        return
      endif
      ap = a
      sum = 1.d0/a
      del = sum
      do n=1,ITMAX
        ap = ap + 1.d0
        del = del*x/ap
        sum = sum + del
        if ( abs(del) .lt. abs(sum)*EPS ) then
          gammaser = sum*exp(-x + a*log(x) - gln)
          return
        endif
      enddo
      ifail=-2
      return
      END FUNCTION gammaser

!=======================================================================

      DOUBLE PRECISION FUNCTION gammacf( a, x, gln, ifail )
      implicit double precision(a-h,o-z)
      INTEGER ITMAX,ifail
      DOUBLE PRECISION a, x, gln, EPS,FPMIN
      PARAMETER ( ITMAX=100, EPS=3.d-16, FPMIN=1.d-30 )
!U    USES gammaln
      INTEGER i
      DOUBLE PRECISION an, b, c, d, del, h

      ifail=0
      gln = gammaln(a)
      b = x + 1.d0 - a
      c = 1.d0/FPMIN
      d = 1.d0/b
      h = d
      do i=1,ITMAX
        an = -i*(i-a)
        b = b + 2.d0
        d = an*d+b
        if ( abs(d) .lt. FPMIN ) d = FPMIN
        c = b + an/c
        if ( abs(c) .lt. FPMIN ) c = FPMIN
        d = 1.d0/d
        del = d*c
        h = h*del
        if ( abs(del-1.d0) .lt. EPS ) then
          ggg =-x + a*log(x) - gln
          gammacf = exp(ggg)*h
          return
        endif
      enddo
      ifail=-3
      return
      END FUNCTION gammacf
!=====================================================================
      DOUBLE PRECISION FUNCTION gammap( a, x, ifail )
      implicit double precision(a-h,o-z)
      INTEGER ifail
      DOUBLE PRECISION a,x
!U    USES gammacf, gammaser
      DOUBLE PRECISION gln

      ifail=0
      if ( x.lt.0.d0 .or. a.le.0.d0 ) then
         ifail=10
         return
      endif
      if ( x .lt. a+1.d0 ) then
        gammap = gammaser( a, x, gln, ifail )
      else
        gammap = 1.d0 - gammacf( a, x, gln, ifail )
      endif
      return
      END FUNCTION gammap

!=======================================================================

      SUBROUTINE gaussleg(x1,x2,x,w,n)
      implicit double precision(a-h,o-z)

      INTEGER n
      DOUBLE PRECISION x1, x2, x(n), w(n)
      INTEGER i, j, m
      DOUBLE PRECISION p1, p2, p3, pp, xl, xm, z, z1, EPS, PI
      PARAMETER ( EPS=1.d-14, PI=3.14159265358d0 )
      m = (n+1)/2
      xm = .5d0*(x2+x1)
      xl = .5d0*(x2-x1)
      do i=1,m
        z = cos(PI*(i-.25d0)/(n+.5d0))
 1      continue
          p1 = 1.d0
          p2 = 0.d0
          do j=1,n
            p3 = p2
            p2 = p1
            p1 = ( (2.d0*j-1.d0)*z*p2 - (j-1.d0)*p3 )/j
          enddo
          pp = n*(z*p1-p2)/(z*z-1.d0)
          z1 = z
          z = z1 - p1/pp
        if ( abs(z-z1) .gt. EPS ) goto 1
        x(i) = xm - xl*z
        x(n+1-i) = xm + xl*z
        w(i) = 2.d0*xl/( (1.d0-z*z)*pp*pp )
        w(n+1-i) = w(i)
      enddo
      return
      END SUBROUTINE gaussleg

!=======================================================================

      double precision function s_invgammap( a )
      implicit double precision(a-h,o-z)
!     +-----------+ Generate f(x)=x*exp(-(x/xw)^xe) using inverse
! fnc.|s_invgammap| Incomplete Gamma Function distribution
!     +-----------+ "P(a,x)=gamma(a,x)/Gamma(a)" -> "x(P)|_(a)".

      double precision a,ff,x1,x2,eps,eta,x,fx
      integer ifail
!      double precision ran1,s_gammap,gammap,zbrent,aa,ffix
!      external         ran1,s_gammap,gammap,zbrent
      integer               it,ifaila
      common  /data_random/ it
      common  /data_gamma_random/ aa,ffix,ifaila
!
!      call random_seed (it)
!      it = 2
!      it = -200
      ff = mtrnd()  ! ff F(x) value, now calculate x(F)

      x1 = 0.d0
      x2 = 1.d+6
      ifail=0
      do while ( gammap( a, x2, ifail ) .lt. ff )
         x2 = x2 * 2.d0
         if (ifail.ne.0) then
            s_invgammap = 0.d0
            return
         endif
      enddo
      ifaila=0
      aa=a
      ffix=ff
      eps=1.d-6
      eta=1.d-6
      x = zbrent( s_gammap, x1, x2, eps, eta, ifail )
      ifail=0
      fx = gammap( a, x, ifail )
      fx=fx
      s_invgammap = x
      return
      end function s_invgammap

!=======================================================================

      double precision function s_gammap( x )
      implicit double precision(a-h,o-z)
!
! fnc.|s_gammap| Define the single-argument-function for the root-finder.
!

!      double precision gammap,x,a,ffix
      double precision x,a,ffix
!      external         gammap
      integer it,ifail
      common  /data_random/ it
      common  /data_gamma_random/ a,ffix,ifail
      s_gammap = gammap( a, x, ifail ) - ffix
      return
      END function s_gammap
!=======================================================================
!
!      DOUBLE PRECISION FUNCTION ran1( idum )
!      implicit double precision(a-h,o-z)
!      INTEGER idum, IA,IM,IQ, IR, NTAB, NDIV
!      DOUBLE PRECISION   AM, EPS, RNMX
!      PARAMETER (IA=16807, IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,&
!                 NTAB=32, NDIV=1+(IM-1)/NTAB, EPS=3.d-16, RNMX=1.d0-EPS)
!      INTEGER j, k, iv(NTAB), iy
!      SAVE iv, iy
!      DATA iv /NTAB*0/, iy /0/
!      if ( idum.le.0 .or. iy.eq.0 ) then
!         idum = max( -idum, 1 )
!         do j=NTAB+8,1,-1
!           k = idum/IQ
!           idum = IA*(idum - k*IQ) - IR*k
!           if (idum.lt.0) idum = idum + IM
!           if (j.le.NTAB) iv(j) = idum
!         enddo
!         iy = iv(1)
!       endif
!       k = idum/IQ
!       idum = IA*(idum - k*IQ) - IR*k
!       if ( idum .lt. 0 ) idum = idum + IM
!       j = 1 + iy/NDIV
!       iy = iv(j)
!       iv(j) = idum
!       ran1 = min( AM*iy, RNMX )
!      return
!      END function ran1
!========================================================================
!      SUBROUTINE rannor( idum, x, y )
!      implicit double precision(a-h,o-z)
!     calls ran1
!      INTEGER idum
!      DOUBLE PRECISION x, y, r, p, a, ran1, DPI
!      DOUBLE PRECISION x, y, r, p, a,DPI
!      PARAMETER( DPI=2.d0*3.14159265358d0 )
!      r = ran1(idum)
!      p = DPI * ran1(idum)
!      a = sqrt( -2.d0 * log(r) )
!      x = a * cos(p)
!      y = a * sin(p)
!      return
!      end subroutine rannor 
!=======================================================================

      DOUBLE PRECISION FUNCTION zbrent( func, x1, x2, eps, eta, ifail )       
      implicit double precision(a-h,o-z)
      EXTERNAL func       
      INTEGER  ifail, iter, ITMAX       
      DOUBLE PRECISION func, x1, x2, eps, eta, EPS0
      PARAMETER(ITMAX=100, EPS0=3.d-16 )
      DOUBLE PRECISION a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
      a = x1
      b = x2
      fa = func(a)       
      fb = func(b)
      ifail = 0
      if ( (fa.gt.0.d0 .and. fb.gt.0.d0)     & 
       .or.(fa.lt.0.d0 .and. fb.lt.0.d0) ) then
         ifail = -1          
         return
      endif       
      c = b
      fc = fb
      do iter=1,ITMAX
         if ( (fb.gt.0.d0 .and. fc.gt.0.d0) &
         .or.(fb.lt.0.d0 .and. fc.lt.0.d0) ) then           
           c = a            ! Rename a, b, c and adjust bounding interval d
           fc = fa
           d = b - a
           e = d
         endif         
         if (abs(fc).lt.abs(fb)) then
           a = b
           b = c
           c = a
           fa = fb
           fb = fc
           fc = fa
         endif         
         tol1 = 2.d0*EPS0*abs(b) + 0.5d0*eps
         xm = .5d0*(c-b)
         if (abs(xm).le.tol1 .or. abs(fb).le.eta) then ! Convergence check
           zbrent = b
           return
         endif
         if (abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
           s = fb/fa    ! Attempt inverse quadratic interpolation
           if(a-c.eq.0.0d0) then
             p = 2.d0*xm*s
             q = 1.d0 - s
           else
             q = fa/fc
             r = fb/fc
             p = s*(2.d0*xm*q*(q-r) - (b-a)*(r-1.d0))
             q = (q-1.d0)*(r-1.d0)*(s-1.d0)
           endif           
           if (p.gt.0.d0) q = -q      ! Check whether in bounds
           p = abs(p)
           if (2.d0*p .lt. min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
             e = d                    ! Accept interpolation
             d = p/q
           else
             d = xm      !Interpolation failed, use bisection
             e = d
           endif         
         else                  ! Bounds decreesing too slowly, use bisection
           d = xm
           e = d
         endif
         a = b  ! Move last best guess to "a"
         fa = fb
         if (abs(d) .gt. tol1) then   ! Evaluate new trial root
           b = b + d
         else           
           b = b + sign(tol1,xm)
         endif
         fb = func(b)
       enddo       
      ifail = 1
      zbrent = b
      return

      end function zbrent


      end module BPMclass
