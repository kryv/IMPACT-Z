!----------------------------------------------------------------
! (c) Copyright, 2003 by the Regents of the University of California.
! Distributionclass: Initial distribution of charged beam bunch class in
!                    Beam module of APPLICATION layer.
! Version: 2.0
! Author: Ji Qiang, Robert Ryne, LBNL, 7/25/03
! Description: This class defines initial distributions for the charged
!              particle beam bunch information in the accelerator.
! Comments: we have added three attributes to each particle:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
      module Distributionclass
        use Pgrid2dclass
        use CompDomclass
        use BeamBunchclass
        use Timerclass
        use NumConstclass
        use PhysConstclass
        use MTrndclass

        character(256) :: particle_dname = ''
        character(256) :: particle_fname = 'partcl.data'
        integer :: iodistfile

      contains
        ! sample the particles with intial distribution.
        subroutine sample_Dist(this,distparam,nparam,flagdist,geom,grid,Flagbc,&
                               nchrg,nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nparam,Flagbc
        integer, intent(inout) :: nchrg
        double precision, dimension(nparam) :: distparam
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer, dimension(nchrg) :: nptlist
        type (BeamBunch), intent(inout) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: flagdist
        integer :: myid, myidx, myidy,seedsize,i,isize
        !integer seedarray(1)
        integer, allocatable, dimension(:) :: seedarray
        double precision rancheck

        iodistfile = 0

        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        if(flagdist.eq.1) then
          call Uniform_Dist(this,nparam,distparam,geom,grid)
        else if(flagdist.eq.2) then
!          call Gauss1_Dist(this,nparam,distparam,geom,grid)
          call Gauss3_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.3) then
          call Waterbag_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.4) then
          call Semigauss_Dist(this,nparam,distparam,geom,grid)
        else if(flagdist.eq.5) then
          call KV3d_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.6) then
          call regen_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.7) then
          call GaussGamma_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.8) then
          call WaterGamma_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.9) then
          call KVGamma_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.10) then
          call regen2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.11) then
          call Gauss4new_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.12) then
          call Gauss7_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.14) then
          call Regen7_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.15) then
          call GaussDouble_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.16) then
          call WaterbagMC_Dist(this,nparam,distparam,grid,0,nchrg,&
                               nptlist,qmcclist,currlist)
        else if(flagdist.eq.17) then
          call GaussMC_Dist(this,nparam,distparam,grid,0,nchrg,&
                               nptlist,qmcclist,currlist)
        else if(flagdist.eq.18) then
          call regenMSU_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.19) then
          call readin_Dist(this,nchrg,nptlist,qmcclist)
        else
          print*,"Initial distribution not available!!"
          iodistfile = 1
          return
          !stop
        endif

        !deallocate(seedarray)

        end subroutine sample_Dist

        subroutine Gauss1_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange
        double precision, dimension(6,2) :: a
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2
        double precision,allocatable,dimension(:,:) :: ptstmp
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj,intvsamp,pid
        double precision :: t0
        double precision, allocatable, dimension(:,:) :: x1,x2,x3,x4,x5,x6

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        avgpts = this%Npt/(npx*npy)

        call getlcrange_CompDom(geom,lcrange)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = 0

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
! The performance inside this loop might be improved due to
! a lot of subroutine call in this loop.
        allocate(x1(2,this%Npt/2))
        allocate(x2(2,this%Npt/2))
        allocate(x3(2,this%Npt/2))
        allocate(x4(2,this%Npt/2))
        allocate(x5(2,this%Npt/2))
        allocate(x6(2,this%Npt/2))
        call normVec(x1,this%Npt/2)
        call normVec(x2,this%Npt/2)
        call normVec(x3,this%Npt/2)
        call normVec(x4,this%Npt/2)
        call normVec(x5,this%Npt/2)
        call normVec(x6,this%Npt/2)
        !print*,"sig1: ",sig1,sig2,sig3,sig4,sig5,sig6
        do ii = 1, this%Npt/2
          !x-px:
!          call normdv(x1)
!          call normdv(x2)
!         Correct Gaussian distribution.
          a(1,1) = xmu1 + sig1*x1(1,ii)/sq12
          a(2,1) = xmu2 + sig2*(-muxpx*x1(1,ii)/sq12+x1(2,ii))
          a(1,2) = xmu1 + sig1*x2(1,ii)/sq12
          a(2,2) = xmu2 + sig2*(-muxpx*x2(1,ii)/sq12+x2(2,ii))
!         Rob's Gaussian distribution.
          !a(1,1) = xmu1 + sig1*x1(1)
          !a(2,1) = xmu2 + sig2*(muxpx*x1(1)+sq12*x1(2))
          !a(1,2) = xmu1 + sig1*x2(1)
          !a(2,2) = xmu2 + sig2*(muxpx*x2(1)+sq12*x2(2))
          !y-py
!          call normdv(x1)
!          call normdv(x2)
!         Correct Gaussian distribution.
          a(3,1) = xmu3 + sig3*x3(1,ii)/sq34
          a(4,1) = xmu4 + sig4*(-muypy*x3(1,ii)/sq34+x3(2,ii))
          a(3,2) = xmu3 + sig3*x4(1,ii)/sq34
          a(4,2) = xmu4 + sig4*(-muypy*x4(1,ii)/sq34+x4(2,ii))
!         Rob's Gaussian distribution.
          !a(3,1) = xmu3 + sig3*x1(1)
          !a(4,1) = xmu4 + sig4*(muypy*x1(1)+sq34*x1(2))
          !a(3,2) = xmu3 + sig3*x2(1)
          !a(4,2) = xmu4 + sig4*(muypy*x2(1)+sq34*x2(2))
          !z-pz
!          call normdv(x1)
!          call normdv(x2)
!         Correct Gaussian distribution.
          a(5,1) = xmu5 + sig5*x5(1,ii)/sq56
          a(6,1) = xmu6 + sig6*(-muzpz*x5(1,ii)/sq56+x5(2,ii))
          a(5,2) = xmu5 + sig5*x6(1,ii)/sq56
          a(6,2) = xmu6 + sig6*(-muzpz*x6(1,ii)/sq56+x6(2,ii))
!         Rob's Gaussian distribution.
          !a(5,1) = xmu5 + sig5*x1(1)
          !a(6,1) = xmu6 + sig6*(muzpz*x1(1)+sq56*x1(2))
          !a(5,2) = xmu5 + sig5*x2(1)
          !a(6,2) = xmu6 + sig6*(muzpz*x2(1)+sq56*x2(2))

          do jj = 1, 2

            pid = (ii-1)*2+jj
          if((npx.gt.1).and.(npy.gt.1)) then

          if((myidx.eq.0).and.(myidy.eq.0)) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.0) then
            if((a(5,jj).le.lcrange(6)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.(npx-1)) then
            if((a(5,jj).gt.lcrange(5)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.0) then
            if((a(3,jj).le.lcrange(4)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.(npy-1)) then
            if((a(3,jj).gt.lcrange(3)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          else
            if( ((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
               .and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          endif

          else if(npx.gt.1) then
            if(myidx.eq.0) then
              if(a(5,jj).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            else if(myidx.eq.(npx-1)) then
              if(a(5,jj).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            else
              if((a(5,jj).gt.lcrange(5)).and.(a(5,jj).le.lcrange(6))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            endif
          else if(npy.gt.1) then
            if(myidy.eq.0) then
              if(a(3,jj).le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            else if(myidy.eq.(npy-1)) then
              if(a(3,jj).gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            else
              if((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            endif
          else
            numpts0 = numpts0 + 1
            do j = 1, 6
              this%Pts1(j,numpts0) = a(j,jj)
            enddo
          endif

          enddo
        enddo

        deallocate(x1)
        deallocate(x2)
        deallocate(x3)
        deallocate(x4)
        deallocate(x5)
        deallocate(x6)

!        call MPI_BARRIER(comm2d,ierr)
        allocate(ptstmp(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 6
            ptstmp(j,i0) = this%Pts1(j,i0)
          enddo
        enddo
        deallocate(this%Pts1)
        allocate(this%Pts1(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 6
            this%Pts1(j,i0) = ptstmp(j,i0)
          enddo
!            this%Pts1(5,i0) = -ptstmp(5,i0)
!            this%Pts1(6,i0) = -ptstmp(6,i0)
        enddo
        deallocate(ptstmp)

        this%Nptlocal = numpts0
!        print*,"numpts0: ",numpts0

!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss1_Dist

        ! sample the particles with intial TRUNCATED GAUSS distribution
        ! using rejection method.
        subroutine Gauss2_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r,r1,r2,x0,y0,x1,x2
        double precision :: zlcmin,zlcmax,ylcmin,ylcmax,fvalue,z0
        double precision :: tmp1,tmp2,xrange,pxrange,yrange,pyrange,zrange,pzrange
        integer :: totnp,npy,npx
        integer :: avgpts,numpts,isamz,isamy
        integer :: myid,myidx,myidy
!        integer seedarray(1)
        double precision :: t0,x11

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        seedarray(1)=(1021+myid)*(myid+7)
!        seedarray(1)=(1121+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        xrange = 4.0*sig1/sqrt(1.0-muxpx*muxpx)
        pxrange = 4.0*sig2/sqrt(1.0-muxpx*muxpx)
        yrange = 4.0*sig3/sqrt(1.0-muypy*muypy)
        pyrange = 4.0*sig4/sqrt(1.0-muypy*muypy)
        zrange = 4.0*sig5/sqrt(1.0-muzpz*muzpz)
        pzrange = 4.0*sig6/sqrt(1.0-muzpz*muzpz)

        call getlcrange_CompDom(geom,lcrange)
        zlcmin = lcrange(5)
        zlcmax = lcrange(6)
        ylcmin = lcrange(3)
        ylcmax = lcrange(4)

        ! scale z and y for sampling purpose.
        zlcmin = zlcmin*sqrt(1.0-muzpz*muzpz)/sig5
        zlcmax = zlcmax*sqrt(1.0-muzpz*muzpz)/sig5
        ylcmin = ylcmin*sqrt(1.0-muypy*muypy)/sig3
        ylcmax = ylcmax*sqrt(1.0-muypy*muypy)/sig3
        if(zlcmax.le.0) then
          z0 = zlcmax
        else if(zlcmin.ge.0) then
          z0 = zlcmin
        else if((zlcmax.gt.0).and.(zlcmin.lt.0)) then
          z0 = 0.0
        else
        endif
        if(ylcmax.le.0) then
          y0 = ylcmax
        else if(ylcmin.ge.0) then
          y0 = ylcmin
        else if((ylcmax.gt.0).and.(ylcmin.lt.0)) then
          y0 = 0.0
        else
        endif
!        write(6,1234)myid,myidx,myidy,zlcmin,zlcmax,ylcmin,ylcmax,z0,y0
 1234   format(3i4,1x,6(1pe10.3,1x))

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        numpts = 0
        isamz = 0
        isamy = 0
        do
          ! rejection sample.
10        call mt_random(r)
          r1 = zlcmin + r*(zlcmax-zlcmin)
          fvalue = exp(-0.5*(r1*r1-z0*z0))
          call mt_random(r2)
          isamz = isamz + 1
          if(r2.gt.fvalue) goto 10
          x1 = r1
20        call mt_random(r)
          r1 = ylcmin + r*(ylcmax-ylcmin)
          fvalue = exp(-0.5*(r1*r1-y0*y0))
          call mt_random(r2)
          isamy = isamy + 1
          if(r2.gt.fvalue) goto 20
          x2 = r1

          numpts = numpts + 1
          if(numpts.gt.avgpts) exit
!          tmp1 = sig5*x1
!          tmp2 = sig3*x2
!          if((abs(tmp1).gt.zrange).or.(abs(tmp2).gt.yrange)) goto 10
          !z-y
          this%Pts1(5,numpts) = xmu5 + sig5*x1/rootz
          this%Pts1(3,numpts) = xmu3 + sig3*x2/rooty
          !rob's distribution:
          !this%Pts1(5,numpts) = xmu5 + sig5*x1
          !this%Pts1(3,numpts) = xmu3 + sig3*x2
          !pz-py
30        call normdv(gs)
          tmp1 = sig6*(-muzpz*x1/rootz+gs(1))
          tmp2 = sig4*(-muypy*x2/rooty+gs(2))
          if((abs(tmp1).gt.pzrange).or.(abs(tmp2).gt.pyrange)) goto 30
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x1/rootz+gs(1))
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x2/rooty+gs(2))
          !rob's distribution:
          !this%Pts1(6,numpts) = xmu6 + sig6*(muzpz*x1+rootz*gs(1))
          !this%Pts1(4,numpts) = xmu4 + sig4*(muypy*x2+rooty*gs(2))
          !x-px
40        call normdv(gs)
          tmp1 = sig1*gs(1)/rootx
          tmp2 = sig2*(-muxpx*gs(1)/rootx+gs(2))
          if((abs(tmp1).gt.xrange).or.(abs(tmp2).gt.pxrange)) goto 40
          this%Pts1(1,numpts) = xmu1 + sig1*gs(1)/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*gs(1)/rootx+gs(2))
          !rob's distribution.
          !this%Pts1(1,numpts) = xmu1 + sig1*gs(1)
          !this%Pts1(2,numpts) = xmu2 + sig2*(muxpx*gs(1)+rootx*gs(2))
        enddo

        this%Nptlocal = avgpts

!        print*,avgpts,isamz,isamy

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss2_Dist

        subroutine Gauss3_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp,pid
!        integer seedarray(1)
        double precision :: t0,x11

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        if(mod(avgpts,10).ne.0) then
          print*,"The number of particles has to be an integer multiple of 10Nprocs"
          stop
        endif

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)

        intvsamp = 10
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
        enddo

        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts

        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss3_Dist

        subroutine Gauss4_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, dimension(2) :: x1
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
!        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1,&
                            ddx,ddy,ddz,ddx1,ddy1,ddz1

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        !dd = 0.000005
        !dd = 0.000010
        ddx = 0.00002
        !ddy = 2.0 ! for test2
        ddy = 0.0003 ! for viz
        !ddy = 0.0005
        !ddx = 0.00002
        !ddy = 0.00003
        !ddx = 0.00001
        !ddy = 0.00001
        ddz = 0.0
        !ddx1 = 0.000001
        !ddy1 = 0.000001
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        numpts = 0
        do
          !x-px:
10        call mt_random(r)
!          r1 = -10*sig1 + r*20*sig1
          r1 = -20*sig1 + r*40*sig1
          call mt_random(r)
!          r2 = -10*sig2 + r*20*sig2
          r2 = -20*sig2 + r*40*sig2
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*r2*muxpx/sig1/sig2+r2*r2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+0.1*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+0.1*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+ddx*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*(r1/sig1)**3+&
                   ddx1*(r1/sig1)**5)*muxpx/sig1/sig2+&
                   (r2+ddx*(r1/sig1)**3+ddx1*(r1/sig1)**5)**2/sig2/sig2)/2)
          call mt_random(r3)
          if(r3.gt.fvalue) goto 10
          xr1 = r1
          px1 = r2
          !y-py
20        call mt_random(r)
          r1 = -20*sig3 + r*40*sig3
          call mt_random(r)
          r2 = -20*sig4 + r*40*sig4
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*r2*muypy/sig3/sig4+r2*r2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-0.0001*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-0.0001*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-ddy*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
! used for previous 3 case studies.
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*(r1/sig3)**3-&
!                   ddy1*(r1/sig3)**5)*muypy/sig3/sig4+&
!                   (r2-ddy*(r1/sig3)**3-ddy1*(r1/sig3)**5)**2/sig4/sig4)/2)
!test 1
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3-r1 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3-r1)**2/sig4/sig4)/2)
!test 2
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3+r2 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3+r2)**2/sig4/sig4)/2)
!test 3
          fvalue = exp(-((r1-ddy*(r2/sig4)**3)**2/sig3/sig3+2*&
                   (r1-ddy*(r2/sig4)**3)*r2*muypy/sig3/sig4+&
                   r2**2/sig4/sig4)/2)
          call mt_random(r3)
          if(r3.gt.fvalue) goto 20
          y1 = r1
          py1 = r2
          numpts = numpts + 1
!          this%Pts1(1,numpts) = xmu1 + xr1*16.8*1.75
!          this%Pts1(2,numpts) = xmu2 + px1/2.5/1.75
!          this%Pts1(3,numpts) = xmu3 + y1*18.7*1.75
!          this%Pts1(4,numpts) = xmu4 + py1/3.9/1.75
          this%Pts1(1,numpts) = xmu1 + xr1
          this%Pts1(2,numpts) = xmu2 + px1
          this%Pts1(3,numpts) = xmu3 + y1
          this%Pts1(4,numpts) = xmu4 + py1
          !z-pz
          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(5,numpts) = xmu5 + sig5*x1(1)/sq56
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x1(1)/sq56+x1(2))
          if(numpts.ge.avgpts) exit
        enddo
        print*,"numpts: ",numpts

        this%Nptlocal = avgpts

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss4_Dist

        subroutine Gauss4new_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,pid
        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1
        double precision, allocatable, dimension(:,:) :: x1
        double precision, allocatable, dimension(:) :: ranumx,ranumy,ranum2x,ranum2y
        integer :: isamx,isamy,iranx,irany,intvsamp

        call starttime_Timer(t0)

        frac = 1.0
        rr = 1.0
        rr(3,:) = 0.0

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)
        if(mod(avgpts,10).ne.0) then
          print*,"The number of particles has to be an integer multiple of 10Nprocs"
          stop
        endif

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        !ddx = 0.0002
        !ddy = 0.0002
        !ddx = 0.00002 !for test2 and 3
        !ddx = 0.0000 !from T3E
        ddx = 0.00004
        !ddy = 2.0 !for test2
        !ddy = 0.0003 !for test3
        !ddy = 0.0001 !from T3E
        ddy = 0.0003
        ddz = 0.0
        !ddx1 = 0.0002 ! for emitdistort2
        !ddy1 = 0.0002 ! for emitdistort2
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        numpts = 0


        !intvsamp = avgpts
        intvsamp = 10
        allocate(ranum2x(2*intvsamp))
        allocate(ranum2y(2*intvsamp))
        allocate(ranumx(intvsamp))
        allocate(ranumy(intvsamp))
        allocate(x1(2,intvsamp))
        call normVec(x1,intvsamp)

        isamx = 0
        isamy = 0
        do
          !x-px:
10        continue
          isamx = isamx + 1
          if(mod(isamx-1,intvsamp).eq.0) then
            call mt_random(ranum2x)
            call mt_random(ranumx)
          endif
          iranx = 2*mod(isamx-1,intvsamp)
!          call mt_random(r)
          r1 = -20*sig1 + ranum2x(iranx+1)*40*sig1
!          call mt_random(r)
          r2 = -20*sig2 + ranum2x(iranx+2)*40*sig2
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*r2*muxpx/sig1/sig2+r2*r2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+0.1*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+0.1*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+ddx*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*(r1/sig1)**3+&
                   ddx1*(r1/sig1)**5)*muxpx/sig1/sig2+&
                   (r2+ddx*(r1/sig1)**3+ddx1*(r1/sig1)**5)**2/sig2/sig2)/2)
!          call mt_random(r3)
          r3 = ranumx(iranx/2+1)
          if(r3.gt.fvalue) goto 10
          xr1 = r1
          px1 = r2
          !y-py
20        continue
          isamy = isamy + 1
          if(mod(isamy-1,intvsamp).eq.0) then
            call mt_random(ranum2y)
            call mt_random(ranumy)
          endif
          irany = 2*mod(isamy-1,intvsamp)
!          call mt_random(r)
          r1 = -20*sig3 + ranum2y(irany+1)*40*sig3
!          call mt_random(r)
          r2 = -20*sig4 + ranum2y(irany+2)*40*sig4
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*r2*muypy/sig3/sig4+r2*r2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-0.0001*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-0.0001*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-ddy*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!used for the previous 3 case studies.
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*(r1/sig3)**3-&
!                   ddy1*(r1/sig3)**5)*muypy/sig3/sig4+&
!                   (r2-ddy*(r1/sig3)**3-ddy1*(r1/sig3)**5)**2/sig4/sig4)/2)
!test 1
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3-r1 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3-r1)**2/sig4/sig4)/2)
!test 2
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3+r2 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3+r2)**2/sig4/sig4)/2)
!test 3
          fvalue = exp(-((r1-ddy*(r2/sig4)**3-ddy1*(r2/sig4)**5)**2/sig3/sig3+2*&
                   (r1-ddy*(r2/sig4)**3-ddy1*(r2/sig4)**5)*r2*muypy/sig3/sig4+&
                   r2**2/sig4/sig4)/2)

          r3 = ranumy(irany/2+1)
!          call mt_random(r3)
          if(r3.gt.fvalue) goto 20
          y1 = r1
          py1 = r2
          numpts = numpts + 1
          this%Pts1(1,numpts) = xmu1 + xr1*factx
          this%Pts1(2,numpts) = xmu2 + px1*factx
          this%Pts1(3,numpts) = xmu3 + y1*facty
          this%Pts1(4,numpts) = xmu4 + py1*facty
          !z-pz
!          call normdv(x1)
          this%Pts1(5,numpts) = xmu5 + sig5*x1(1,numpts)/sq56*factz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x1(1,numpts)/sq56+x1(2,numpts))*factz
          if(numpts.ge.avgpts) exit
        enddo

        deallocate(ranum2x)
        deallocate(ranum2y)
        deallocate(ranumx)
        deallocate(ranumy)
        deallocate(x1)

        this%Nptlocal = avgpts

        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1)
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2)

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xmu1 + xxx
          this%Pts1(2,i) = xmu2 + yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xmu3 + xxx
          this%Pts1(4,i) = xmu4 + yyy
        enddo
        print*,"numpts: ",numpts

        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        end subroutine Gauss4new_Dist

        subroutine Gauss5_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, dimension(2) :: x1
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
!        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: al0x,al1x,ga0x,ga1x,b0x,b1x,r11x,r22x,r12x
        double precision :: al0y,al1y,ga0y,ga1y,b0y,b1y,r11y,r22y,r12y
        double precision :: al0z,al1z,ga0z,ga1z,b0z,b1z,r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        !used for the original comparison at beginning
        factx = 0.93894
        facty = 0.72253
        !factx = 0.971047
        !facty = 0.793620
        !used for 3rd order comparison at quad 20
        !factx = 0.95776371
        !facty = 0.85342924
        ! used for 5th order comparison at quad 20
        !factx = 0.8466795
        !facty = 0.8491518
        factz = 1.0
        !dd = 0.000005
        !dd = 0.000010
        !used for the original comparison at beginning
        ddx = 0.0002
        ddy = 0.0002
        !used for 3rd order comparison at quad 20
        !ddx = 0.00002
        !ddy = 0.00003
        ! used for 5th order comparison at quad 20
        !ddx = 0.00001
        !ddy = 0.00001
        !ddx = 0.0
        !ddy = 0.0
        ddz = 0.0
        !ddx1 = 0.000001
        !ddy1 = 0.000001
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        numpts = 0

!test 2
        ddx = 0.00002
        ddy = 2.0
        factx = 0.95618
        facty = 6.55394
        do
          !x-px:
10        call mt_random(r)
          !used for the original comparison at beginning
!          r1 = -10*sig1 + r*20*sig1
          ! used for 5th order comparison at quad 20
          r1 = -20*sig1 + r*40*sig1
          call mt_random(r)
          !used for the original comparison at beginning
!          r2 = -10*sig2 + r*20*sig2
          ! used for 5th order comparison at quad 20
          r2 = -20*sig2 + r*40*sig2
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*r2*muxpx/sig1/sig2+r2*r2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+0.1*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+0.1*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+ddx*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*(r1/sig1)**3+&
                   ddx1*(r1/sig1)**5)*muxpx/sig1/sig2+&
                   (r2+ddx*(r1/sig1)**3+ddx1*(r1/sig1)**5)**2/sig2/sig2)/2)
          call mt_random(r3)
          if(r3.gt.fvalue) goto 10
          xr1 = r1
          px1 = r2
          !y-py
20        call mt_random(r)
          r1 = -20*sig3 + r*40*sig3
          call mt_random(r)
          r2 = -20*sig4 + r*40*sig4
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*r2*muypy/sig3/sig4+r2*r2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-0.0001*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-0.0001*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-ddy*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!used for the previous 3 case studies.
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*(r1/sig3)**3-&
!                   ddy1*(r1/sig3)**5)*muypy/sig3/sig4+&
!                   (r2-ddy*(r1/sig3)**3-ddy1*(r1/sig3)**5)**2/sig4/sig4)/2)
!test 1
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3-r1 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3-r1)**2/sig4/sig4)/2)
!test 2
          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3+r2 &
                   )*muypy/sig3/sig4+&
                   (ddy*(r2/sig4)**3+r2)**2/sig4/sig4)/2)

          call mt_random(r3)
          if(r3.gt.fvalue) goto 20
          y1 = r1
          py1 = r2
          numpts = numpts + 1
!          this%Pts1(1,numpts) = xmu1 + xr1*16.8*1.75
!          this%Pts1(2,numpts) = xmu2 + px1/2.5/1.75
!          this%Pts1(3,numpts) = xmu3 + y1*18.7*1.75
!          this%Pts1(4,numpts) = xmu4 + py1/3.9/1.75
          this%Pts1(1,numpts) = xmu1 + xr1*factx
          this%Pts1(2,numpts) = xmu2 + px1*factx
          this%Pts1(3,numpts) = xmu3 + y1*facty
          this%Pts1(4,numpts) = xmu4 + py1*facty
          !z-pz
          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(5,numpts) = xmu5 + sig5*x1(1)/sq56*factz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x1(1)/sq56+x1(2))*factz
          if(numpts.ge.avgpts) exit
        enddo
        print*,"numpts: ",numpts

        this%Nptlocal = avgpts
        !al0x = 1.16489
        !al1x = 1.56003
        !b0x = 2.123973
        !b1x = 1.977976
        al0x = 1.41627
        !al1x = 1.91224
        al1x = 1.66392
        b0x = 3.215056
        !b1x = 2.779062
        b1x = 2.237308
! test 2
        al0x = 1.16954
        al1x = 1.59331
        b0x = 2.12649
        b1x = 1.96820
        ga0x = (1.0+al0x*al0x)/b0x*0.1363
        ga1x = (1.0+al1x*al1x)/b1x*0.1363
        r11x = sqrt(ga1x/ga0x)
        r22x = 1.0/r11x
        r12x = (al1x-al0x)/sqrt(ga1x*ga0x)
        al0y = -1.41952
        !al1y = -1.94419
        al1y = -1.65727
        b0y = 3.236013
        !b1y = 2.318362
        b1y = 2.2658332
! test 2
        al0y = -1.63889
        b0y = 4.58248
        al1y = -0.77426
        b1y = 101.5378
        ga0y = (1.0+al0y*al0y)/b0y*0.1363
        ga1y = (1.0+al1y*al1y)/b1y*0.1363
        r11y = sqrt(ga1y/ga0y)
        r22y = 1.0/r11y
        r12y = (al1y-al0y)/sqrt(ga1y*ga0y)
!        r11x = 1.0
!        r12x = 0.0
!        r22x = 1.0
!        r11y = 1.0
!        r12y = 0.0
!        r22y = 1.0
        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xmu1 + xxx
          this%Pts1(2,i) = xmu2 + yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xmu3 + xxx
          this%Pts1(4,i) = xmu4 + yyy
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss5_Dist

        subroutine Gauss6_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, dimension(2) :: x1
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
!        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: al0x,al1x,ga0x,ga1x,b0x,b1x,r11x,r22x,r12x
        double precision :: al0y,al1y,ga0y,ga1y,b0y,b1y,r11y,r22y,r12y
        double precision :: al0z,al1z,ga0z,ga1z,b0z,b1z,r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1
        double precision, dimension(3) :: frac,al0,ga0,epson0,al1,ga1,epson1
        double precision, dimension(3,3) :: rr

        call starttime_Timer(t0)

        frac = 1.0
        rr = 1.0
        rr(3,:) = 0.0

        call Waterbag_Dist(this,nparam,distparam,grid,0)
        call gammaepson_Dist(this,al0,ga0,epson0)
        !print*,"a10: ",al0
        !print*,"epson0: ",epson0
        call Gaussdistort_Dist(this,nparam,distparam,grid,rr,frac,0)
        call gammaepson_Dist(this,al1,ga1,epson1)
        !print*,"a11: ",al1
        !print*,"epson1: ",epson1

        do i = 1, 3
          frac(i) = sqrt(epson0(i)/epson1(i))
          rr(1,i) = sqrt(ga1(i)/ga0(i))
          rr(2,i) = 1.0/rr(1,i)
          rr(3,i) = (al1(i)-al0(i))/sqrt(ga1(i)*ga0(i))
        enddo
        !no distort for z-pz yet.
        !frac(3) = 1.0

        call Gaussdistort_Dist(this,nparam,distparam,grid,rr,frac,1)

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss6_Dist

        subroutine Gaussdistort_Dist(this,nparam,distparam,grid,rr,frac,resamp)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,resamp
        double precision, dimension(nparam) :: distparam
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1,z1,pz1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1,xx

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call mt_random(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        if(resamp.eq.1) then
          call Waterbag_Dist(this,nparam,distparam,grid,1)
        endif

        ddx = -4.0e2 !for alphax > 0.0
        ddx = -8.0e2 !for alphax > 0.0
        !ddy = -3.0e2
        ddy = 1.0e7  !for alphay < 0.0
        ddy = 3.0e7  !for alphay < 0.0
        !ddy = 6.0e7  !for alphay < 0.0 for the purpose of display
        !ddz = -5.0e-3
        !ddz = 5.0e11 !for alphaz < 0.0
        ddz = 0.0
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        print*,"prepare:"
        do i = 1, avgpts
          xx = this%Pts1(1,i) - xmu1
          this%Pts1(2,i) = this%Pts1(2,i) + ddx*xx**3
          !for alpha > 0
          !xx = this%Pts1(3,i) - xmu3
          !this%Pts1(4,i) = this%Pts1(4,i) + ddy*xx**3
          !for alpha < 0
          xx = this%Pts1(4,i) - xmu4
          this%Pts1(3,i) = this%Pts1(3,i) + ddy*xx**3
          !xx = this%Pts1(5,i) - xmu5
          !this%Pts1(6,i) = this%Pts1(6,i) + ddz*xx**3
          xx = this%Pts1(6,i) - xmu6
          this%Pts1(5,i) = this%Pts1(5,i) + ddz*xx**3
        enddo

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        do i = 1, avgpts
          this%Pts1(1,i) = this%Pts1(1,i)*factx
          this%Pts1(2,i) = this%Pts1(2,i)*factx
          this%Pts1(3,i) = this%Pts1(3,i)*facty
          this%Pts1(4,i) = this%Pts1(4,i)*facty
          this%Pts1(5,i) = this%Pts1(5,i)*factz
          this%Pts1(6,i) = this%Pts1(6,i)*factz
        enddo

        this%Nptlocal = avgpts

        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1)
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2)
        r11z = rr(1,3)
        r22z = rr(2,3)
        r12z = rr(3,3)

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xxx
          this%Pts1(2,i) = yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xxx
          this%Pts1(4,i) = yyy
          xxx = this%Pts1(5,i)*r11z + this%Pts1(6,i)*r12z
          yyy = this%Pts1(6,i)*r22z
          this%Pts1(5,i) = xxx
          this%Pts1(6,i) = yyy
        enddo

        end subroutine Gaussdistort_Dist

        subroutine Gaussdistortold_Dist(this,nparam,distparam,grid,rr,frac)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
!        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1,z1,pz1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1
        double precision, allocatable, dimension(:) :: ranumx,ranumy,&
                          ranum2x,ranum2y,ranumz,ranum2z
        integer :: isamx,isamy,iranx,irany,intvsamp,isamz,iranz

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)
        if(mod(avgpts,10).ne.0) then
          print*,"The number of particles has to be an integer multiple of 10Nprocs"
          stop
        endif

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        !allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        !ddx = 0.0002
        !ddy = 0.0002
        !ddx = 0.00002 !for test2 and 3
        !ddx = 0.0000 !from T3E
        ddx = 0.00004
        !ddy = 2.0 !for test2
        !ddy = 0.0003 !for test3
        !ddy = 0.0001 !from T3E
        ddy = 0.0003
        !ddz = 0.3
        !ddz = 0.0005
        ddz = 3.0005
        !ddx1 = 0.0002 ! for emitdistort2
        !ddy1 = 0.0002 ! for emitdistort2
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        numpts = 0

        !intvsamp = avgpts
        intvsamp = 10
        allocate(ranum2x(2*intvsamp))
        allocate(ranum2y(2*intvsamp))
        allocate(ranum2z(2*intvsamp))
        allocate(ranumx(intvsamp))
        allocate(ranumy(intvsamp))
        allocate(ranumz(intvsamp))

        isamx = 0
        isamy = 0
        isamz = 0
        do
          !x-px:
10        continue
          isamx = isamx + 1
          if(mod(isamx-1,intvsamp).eq.0) then
            call mt_random(ranum2x)
            call mt_random(ranumx)
          endif
          iranx = 2*mod(isamx-1,intvsamp)
          r1 = -20*sig1 + ranum2x(iranx+1)*40*sig1
          r2 = -20*sig2 + ranum2x(iranx+2)*40*sig2
          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*(r1/sig1)**3+&
                   ddx1*(r1/sig1)**5)*muxpx/sig1/sig2+&
                   (r2+ddx*(r1/sig1)**3+ddx1*(r1/sig1)**5)**2/sig2/sig2)/2)
          r3 = ranumx(iranx/2+1)
          if(r3.gt.fvalue) goto 10
          xr1 = r1
          px1 = r2
          !y-py
20        continue
          isamy = isamy + 1
          if(mod(isamy-1,intvsamp).eq.0) then
            call mt_random(ranum2y)
            call mt_random(ranumy)
          endif
          irany = 2*mod(isamy-1,intvsamp)
          r1 = -20*sig3 + ranum2y(irany+1)*40*sig3
          r2 = -20*sig4 + ranum2y(irany+2)*40*sig4
          fvalue = exp(-((r1-ddy*(r2/sig4)**3-ddy1*(r2/sig4)**5)**2/sig3/sig3+2*&
                   (r1-ddy*(r2/sig4)**3-ddy1*(r2/sig4)**5)*r2*muypy/sig3/sig4+&
                   r2**2/sig4/sig4)/2)

          r3 = ranumy(irany/2+1)
          if(r3.gt.fvalue) goto 20
          y1 = r1
          py1 = r2
          !z-pz
30        continue
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call mt_random(ranum2z)
            call mt_random(ranumz)
          endif
          iranz = 2*mod(isamz-1,intvsamp)
          r1 = -20*sig5 + ranum2z(iranz+1)*40*sig5
          r2 = -20*sig6 + ranum2z(iranz+2)*40*sig6
          fvalue = exp(-((r1-ddz*(r2/sig5)**3-ddz1*(r2/sig6)**5)**2/sig5/sig5+2*&
                   (r1-ddz*(r2/sig6)**3-ddz1*(r2/sig6)**5)*r2*muzpz/sig5/sig6+&
                   r2**2/sig6/sig6)/2)
          r3 = ranumz(iranz/2+1)
          if(r3.gt.fvalue) goto 30
          z1 = r1
          pz1 = r2

          numpts = numpts + 1
          this%Pts1(1,numpts) = xmu1 + xr1*factx
          this%Pts1(2,numpts) = xmu2 + px1*factx
          this%Pts1(3,numpts) = xmu3 + y1*facty
          this%Pts1(4,numpts) = xmu4 + py1*facty
          this%Pts1(5,numpts) = xmu5 + z1*factz
          this%Pts1(6,numpts) = xmu6 + pz1*factz

          if(numpts.ge.avgpts) exit
        enddo

        deallocate(ranum2x)
        deallocate(ranum2y)
        deallocate(ranum2z)
        deallocate(ranumx)
        deallocate(ranumy)
        deallocate(ranumz)

        this%Nptlocal = avgpts

        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1)
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2)
        r11z = rr(1,3)
        r22z = rr(2,3)
        r12z = rr(3,3)

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xmu1 + xxx
          this%Pts1(2,i) = xmu2 + yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xmu3 + xxx
          this%Pts1(4,i) = xmu4 + yyy
          xxx = this%Pts1(5,i)*r11z + this%Pts1(6,i)*r12z
          yyy = this%Pts1(6,i)*r22z
          this%Pts1(5,i) = xmu5 + xxx
          this%Pts1(6,i) = xmu6 + yyy
        enddo
        print*,"numpts: ",numpts,isamx,isamy,isamz

        end subroutine Gaussdistortold_Dist

        subroutine gammaepson_Dist(this,CSalpha,CSgamma,CSepson)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: this
        double precision, dimension(3), intent(out) :: CSalpha,CSgamma,CSepson
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,yrms
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0
        double precision ::sqx,sqpx,sqy,sqpy,sqz,sqpz
        integer :: i,ierr
        double precision:: qmc,xl,xt,pi
        double precision, dimension(15) :: tmplc,tmpgl
        double precision, dimension(3)  :: CSbeta

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg

        innp = this%Nptlocal
        nptot = this%Npt

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
        pi = 2*asin(1.0)

        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*this%Pts1(2,i)
          px0lc = px0lc + this%Pts1(2,i)
          sqsum2local = sqsum2local + this%Pts1(2,i)*this%Pts1(2,i)
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          ypylocal = ypylocal + this%Pts1(3,i)*this%Pts1(4,i)
          py0lc = py0lc + this%Pts1(4,i)
          sqsum4local = sqsum4local + this%Pts1(4,i)*this%Pts1(4,i)
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)

          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
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

        call MPI_ALLREDUCE(tmplc,tmpgl,15,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)

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

        epsx2 = (sqsum1*sqsum2-xpx*xpx)
        epsy2 = (sqsum3*sqsum4-ypy*ypy)
        epsz2 = (sqsum5*sqsum6-zpz*zpz)
        epx = sqrt(max(epsx2,0.0d0))
        epy = sqrt(max(epsy2,0.0d0))
        epz = sqrt(max(epsz2,0.0d0))
        xrms = sqrt(sqsum1)
        yrms = sqrt(sqsum3)
        zrms = sqrt(sqsum5)

        CSalpha(1) = -xpx/epx
        CSalpha(2) = -ypy/epy
        CSalpha(3) = -zpz/epz
        CSbeta(1) = (xrms*xl)**2/(epx*xl)
        CSbeta(2) = (yrms*xl)**2/(epy*xl)
        CSbeta(3) = (zrms*xt)**2/(epz*qmc*xt)

        CSgamma(:) = (1.0 + CSalpha(:)*CSalpha(:))/CSbeta(:)*xl
        CSgamma(3) = (1.0 + CSalpha(3)*CSalpha(3))/CSbeta(3)/qmc*xt

        CSepson(1) = epx*xl
        CSepson(2) = epy*xl
        CSepson(3) = epz*qmc*xt

        end subroutine gammaepson_Dist

        subroutine normdv(y)
        implicit none
        include 'mpif.h'
        double precision, dimension(2), intent(out) :: y
        double precision :: twopi,x1,x2,epsilon

        epsilon = 1.0e-18

        twopi = 4.0*asin(1.0)
        call mt_random(x2)
10      call mt_random(x1)
!        x1 = 0.5
!10      x2 = 0.6
        if(x1.eq.0.0) goto 10
!        if(x1.eq.0.0) x1 = epsilon
        y(1) = sqrt(-2.0*log(x1))*cos(twopi*x2)
        y(2) = sqrt(-2.0*log(x1))*sin(twopi*x2)

        end subroutine normdv

        subroutine normVec(y,num)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: num
        double precision, dimension(2,num), intent(out) :: y
        double precision :: twopi,epsilon
        double precision, dimension(num) :: x1,x2
        integer :: i

        epsilon = 1.0e-18

        twopi = 4.0*asin(1.0)
        call mt_random(x2)
        call mt_random(x1)
        do i = 1, num
          if(x1(i).eq.0.0) x1(i) = epsilon
          y(1,i) = sqrt(-2.0*log(x1(i)))*cos(twopi*x2(i))
          y(2,i) = sqrt(-2.0*log(x1(i)))*sin(twopi*x2(i))
        enddo

        end subroutine normVec

        ! sample the particles with intial distribution
        ! using rejection method.
        subroutine Waterbag_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts,isamz,isamy
        integer :: myid,myidx,myidy,iran,intvsamp,pid,j
!        integer seedarray(2)
        double precision :: t0,x11
        double precision, allocatable, dimension(:) :: ranum6

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        seedarray(2)=(101+2*myid)*(myid+4)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray)
        call mt_random(x11)
        !print*,"x11: ",myid,x11

        avgpts = this%Npt/(npx*npy)
        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs"
        !  stop
        !endif

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif
        numpts = 0
        isamz = 0
        isamy = 0
        intvsamp = avgpts
        !intvsamp = 10
        allocate(ranum6(6*intvsamp))

        do
          ! rejection sample.
10        continue
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call mt_random(ranum6)
!            call random_number(ranum6)
          endif
          iran = 6*mod(isamz-1,intvsamp)
          r1 = 2.0*ranum6(iran+1)-1.0
          r2 = 2.0*ranum6(iran+2)-1.0
          r3 = 2.0*ranum6(iran+3)-1.0
          r4 = 2.0*ranum6(iran+4)-1.0
          r5 = 2.0*ranum6(iran+5)-1.0
          r6 = 2.0*ranum6(iran+6)-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 10
          isamy = isamy + 1
          numpts = numpts + 1
          if(numpts.gt.avgpts) exit
!x-px:
          x1 = r1*sqrt(8.0)
          x2 = r2*sqrt(8.0)
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          x3 = r3*sqrt(8.0)
          x4 = r4*sqrt(8.0)
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          x5 = r5*sqrt(8.0)
          x6 = r6*sqrt(8.0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
        enddo

        deallocate(ranum6)

        this%Nptlocal = avgpts
        !print*,"avgpts: ",avgpts

!        print*,avgpts,isamz,isamy
        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Waterbag_Dist

        subroutine KV3d_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,j,pid
!        integer seedarray(1)
        double precision :: t0,x11,twopi

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        twopi = 4*asin(1.0)

        do numpts = 1, avgpts
          call mt_random(r1)
          call mt_random(r2)
          call mt_random(r3)
          r4 = sqrt(r1)
          r5 = sqrt(1.0-r1)
          r2 = r2*twopi
          r3 = r3*twopi
          x1 = 2*r4*cos(r2)
          x2 = 2*r4*sin(r2)
          x3 = 2*r5*cos(r3)
          x4 = 2*r5*sin(r3)
!x-px:
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation.
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          call mt_random(r5)
          r5 = 2*r5 - 1.0
          call mt_random(r6)
          r6 = 2*r6 - 1.0
          x5 = r5*sqrt(3.0)
          x6 = r6*sqrt(3.0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
        enddo

        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine KV3d_Dist

! This sampling does not work properly if one-dimensional PE is 1
        subroutine Uniform_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange
        double precision, dimension(6,1) :: a
        double precision, dimension(2) :: x1, x2
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj,pid
        double precision :: t0

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        avgpts = this%Npt/(npx*npy)

        call getlcrange_CompDom(geom,lcrange)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = 0

        ! initial allocate 'avgpts' particles on each processor.
        this%Pts1 = 0.0

        do ii = 1, this%Npt
          call mt_random(r1)
          r1 = (2*r1 - 1.0)*sqrt(3.0)
          call mt_random(r2)
          r2 = (2*r2 - 1.0)*sqrt(3.0)
          a(1,1) = xmu1 + sig1*r1/sq12
          a(2,1) = xmu2 + sig2*(-muxpx*r2/sq12+r2)
          call mt_random(r1)
          r1 = (2*r1 - 1.0)*sqrt(3.0)
          call mt_random(r2)
          r2 = (2*r2 - 1.0)*sqrt(3.0)
          a(3,1) = xmu3 + sig3*r1/sq34
          a(4,1) = xmu4 + sig4*(-muypy*r2/sq34+r2)
          call mt_random(r1)
          r1 = (2*r1 - 1.0)*sqrt(3.0)
          call mt_random(r2)
          r2 = (2*r2 - 1.0)*sqrt(3.0)
          a(5,1) = xmu5 + sig5*r1/sq56
          a(6,1) = xmu6 + sig6*(-muzpz*r2/sq56+r2)

          do jj = 1, 1

          if(totnp.ne.1) then

          if((myidx.eq.0).and.(myidy.eq.0)) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidx.eq.0) then
            if((a(5,jj).le.lcrange(6)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidx.eq.(npx-1)) then
            if((a(5,jj).gt.lcrange(5)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidy.eq.0) then
            if((a(3,jj).le.lcrange(4)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidy.eq.(npy-1)) then
            if((a(3,jj).gt.lcrange(3)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          else
            if( ((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
               .and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          endif

          else
            numpts0 = numpts0 + 1
            do j = 1, 6
              this%Pts1(j,numpts0) = a(j,jj)
            enddo
          endif

          enddo
        enddo

!        call MPI_BARRIER(comm2d,ierr)

        this%Nptlocal = numpts0

!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)
        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo


        end subroutine Uniform_Dist

! This sampling does not work properly if one-dimensional PE is 1
        subroutine Semigauss_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange
        double precision, dimension(6,2) :: a
        double precision, dimension(3) :: x1, x2
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2, r, r3
        double precision,allocatable,dimension(:,:) :: ptstmp
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj,pid
        double precision :: t0

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        avgpts = this%Npt/(npx*npy)

        call getlcrange_CompDom(geom,lcrange)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = 0

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
! The performance inside this loop might be improved due to
! a lot of subroutine call in this loop.
        do ii = 1, this%Npt
          ! rejection sample.
10        call mt_random(r1)
          call mt_random(r2)
          call mt_random(r3)
          r1 = 2.0*r1-1.0
          r2 = 2.0*r2-1.0
          r3 = 2.0*r3-1.0
          if(r1**2+r2**2+r3**2.gt.1.0) goto 10
          x2(1) = r1
          x2(2) = r2
          x2(3) = r3
          call normdv2(x1)

          !x-px:
!         Correct Gaussian distribution.
          a(1,1) = xmu1 + sig1*x2(1)/sq12*sqrt(5.0)
          a(2,1) = xmu2 + sig2*(-muxpx*x2(1)/sq12+x1(1))
!         Rob's Gaussian distribution.
          !a(1,1) = xmu1 + sig1*x2(1)*sqrt(5.0)
          !a(2,1) = xmu2 + sig2*(muxpx*x2(1)+sq12*x1(1))
          !y-py
!         Correct Gaussian distribution.
          a(3,1) = xmu3 + sig3*x2(2)/sq34*sqrt(5.0)
          a(4,1) = xmu4 + sig4*(-muypy*x2(2)/sq34+x1(2))
!         Rob's Gaussian distribution.
          !a(3,1) = xmu3 + sig3*x2(2)*sqrt(5.0)
          !a(4,1) = xmu4 + sig4*(muypy*x2(2)+sq34*x1(2))
          !z-pz
!         Correct Gaussian distribution.
          a(5,1) = xmu5 + sig5*x2(3)/sq56*sqrt(5.0)
          a(6,1) = xmu6 + sig6*(-muzpz*x2(3)/sq56+x1(3))
!         Rob's Gaussian distribution.
          !a(5,1) = xmu5 + sig5*x2(3)*sqrt(5.0)
          !a(6,1) = xmu6 + sig6*(muzpz*x2(3)+sq56*x1(3))

          do jj = 1, 1

          if(totnp.ne.1) then

          if((myidx.eq.0).and.(myidy.eq.0)) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.0) then
            if((a(5,jj).le.lcrange(6)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.(npx-1)) then
            if((a(5,jj).gt.lcrange(5)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.0) then
            if((a(3,jj).le.lcrange(4)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.(npy-1)) then
            if((a(3,jj).gt.lcrange(3)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          else
            if( ((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
               .and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          endif

          else
            numpts0 = numpts0 + 1
            do j = 1, 6
              this%Pts1(j,numpts0) = a(j,jj)
            enddo
            this%Pts1(7,numpts0) = this%Charge/this%mass
            this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
            this%Pts1(9,numpts0) = ii

          endif

          enddo
        enddo

!        call MPI_BARRIER(comm2d,ierr)
        allocate(ptstmp(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 9
            ptstmp(j,i0) = this%Pts1(j,i0)
          enddo
        enddo
        deallocate(this%Pts1)
        allocate(this%Pts1(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 9
            this%Pts1(j,i0) = ptstmp(j,i0)
          enddo
        enddo
        deallocate(ptstmp)

        this%Nptlocal = numpts0
!        print*,"numpts0: ",numpts0

!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Semigauss_Dist

        subroutine normdv2(y)
        implicit none
        include 'mpif.h'
        double precision, dimension(3), intent(out) :: y
        double precision :: sumtmp,x
        integer :: i

        sumtmp = 0.0
        do i = 1, 12
          call mt_random(x)
          sumtmp = sumtmp + x
        enddo
        y(1) = sumtmp - 6.0

        sumtmp = 0.0
        do i = 1, 12
          call mt_random(x)
          sumtmp = sumtmp + x
        enddo
        y(2) = sumtmp - 6.0

        sumtmp = 0.0
        do i = 1, 12
          call mt_random(x)
          sumtmp = sumtmp + x
        enddo
        y(3) = sumtmp - 6.0

        end subroutine normdv2

        subroutine regen_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:6,i)
          sumx = sumx + Ptcl(1,i)
          sumx2 = sumx2 + Ptcl(1,i)*Ptcl(1,i)
          sumy = sumy + Ptcl(3,i)
          sumy2 = sumy2 + Ptcl(3,i)*Ptcl(3,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumy2
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        do j = 1, inipts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
          gammabet = sqrt(gamma*gamma-1.0)
          Ptcl(2,j) = (Ptcl(2,j)*gammabet)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gammabet)*pyscale + xmu4
!          betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
!          Ptcl(2,j) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
!          Ptcl(4,j) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j))
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j))
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j))
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j))
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j))
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j))
        enddo

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then

          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call mt_random(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call mt_random(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.04*xmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(1,ii) = Ptcl(3,i)+r
!              this%Pts1(1,ii) = Ptcl(1,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.04*pxmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(2,ii) = Ptcl(4,i)+r
!              this%Pts1(2,ii) = Ptcl(2,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                !r = (2.0*r-1.0)*0.04*ymax
                r = (2.0*r-1.0)*0.02*ymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(3,ii) = Ptcl(1,i)+r
!              this%Pts1(3,ii) = Ptcl(3,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                !r = (2.0*r-1.0)*0.04*pymax
                r = (2.0*r-1.0)*0.02*pymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(4,ii) = Ptcl(2,i)+r
!              this%Pts1(4,ii) = Ptcl(4,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                !r = (2.0*r-1.0)*0.04*zmax
                r = (2.0*r-1.0)*0.005*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                !r = (2.0*r-1.0)*0.04*pzmax
                r = (2.0*r-1.0)*0.002*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(2,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        deallocate(Ptcl)

        do j = 1, this%Nptlocal
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        end subroutine regen_Dist

        subroutine GaussGamma_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4
        double precision :: sq12,sq34
        double precision, dimension(2) :: x1
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,inipts,j,pid
!        integer seedarray(1)
        double precision :: t0,x11
        double precision :: alphaz,alphapz,lambdaz,lambdapz,pi,synangle,&
        kenergy,mccq,tmp1,tmp2,Emin,Emax,phimin,phimax,gamma0,fe,femax,fphi,&
        fphimax,e1,phi1,r1,r2,r3,r4

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        alphaz = distparam(15)
        lambdaz = distparam(16)
        alphapz = distparam(17)
        lambdapz = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        print*,myid,x11

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
        close(12)
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq
        pi = 2*asin(1.0)
        synangle = synangle*pi/180.0

        Emin = 0.0
        !Emax = 6.9
        Emax = 7.5
        phimin = -1.8
        phimax = 0.6

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        numpts = 0
        print*,"kenergy: ",kenergy,synangle,avgpts,alphaz,lambdaz,&
               alphapz,lambdapz
        tmp1 = Emax - alphaz/lambdaz
        tmp2 = phimin + alphapz/lambdapz
        do
          ! rejection sample.
10        call mt_random(r1)
          r1 = Emin + (Emax-Emin)*r1
          femax = (Emax - tmp1)**alphaz*exp(-lambdaz*(Emax-tmp1))
          fe = (Emax - r1)**alphaz*exp(-lambdaz*(Emax-r1))/femax
          call mt_random(r3)
          if(r3.gt.fe) goto 10
          e1 = r1

20        call mt_random(r2)
          r2 = phimin + (phimax-phimin)*r2
          fphimax = (tmp2-phimin)**alphapz*exp(-lambdapz*(tmp2-phimin))
          fphi = (r2-phimin)**alphapz*exp(-lambdapz*(r2-phimin))/fphimax
          call mt_random(r4)
          if(r4.gt.fphi) goto 20
          phi1 = r2

          numpts = numpts + 1
          if(numpts.gt.avgpts) exit

          !x-px:
          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(1,numpts) = xmu1 + sig1*x1(1)/sq12
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1(1)/sq12+x1(2))
          !y-py
          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(3,numpts) = xmu3 + sig3*x1(1)/sq34
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x1(1)/sq34+x1(2))
          !z-pz
          this%Pts1(5,numpts) = xmu5 + (phi1-synangle)
          this%Pts1(6,numpts) = xmu6 + gamma0 - (1+e1*1.0e6/mccq)
        enddo

        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine GaussGamma_Dist

        subroutine WaterGamma_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4
        double precision :: sq12,sq34
        double precision, dimension(2) :: x1
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,inipts,j,pid
!        integer seedarray(1)
        double precision :: t0,x11
        double precision :: alphaz,alphapz,lambdaz,lambdapz,pi,synangle,&
        kenergy,mccq,tmp1,tmp2,Emin,Emax,phimin,phimax,gamma0,fe,femax,fphi,&
        fphimax,e1,phi1,r1,r2,r3,r4,r5,r6,xx1,xx2,xx3,xx4

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        alphaz = distparam(15)
        lambdaz = distparam(16)
        alphapz = distparam(17)
        lambdapz = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        print*,myid,x11

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
        close(12)
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq
        pi = 2*asin(1.0)
        synangle = synangle*pi/180.0

        Emin = 0.0
        !Emax = 6.9
        Emax = 7.5
        phimin = -1.8
        phimax = 0.6

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        numpts = 0
        print*,"kenergy: ",kenergy,synangle,avgpts,alphaz,lambdaz,&
               alphapz,lambdapz
        tmp1 = Emax - alphaz/lambdaz
        tmp2 = phimin + alphapz/lambdapz
        do
          ! rejection sample.
10        call mt_random(r1)
          r1 = Emin + (Emax-Emin)*r1
          femax = (Emax - tmp1)**alphaz*exp(-lambdaz*(Emax-tmp1))
          fe = (Emax - r1)**alphaz*exp(-lambdaz*(Emax-r1))/femax
          call mt_random(r3)
          if(r3.gt.fe) goto 10
          e1 = r1

20        call mt_random(r2)
          r2 = phimin + (phimax-phimin)*r2
          fphimax = (tmp2-phimin)**alphapz*exp(-lambdapz*(tmp2-phimin))
          fphi = (r2-phimin)**alphapz*exp(-lambdapz*(r2-phimin))/fphimax
          call mt_random(r4)
          if(r4.gt.fphi) goto 20
          phi1 = r2

30        call mt_random(r1)
          call mt_random(r2)
          call mt_random(r3)
          call mt_random(r4)
          call mt_random(r5)
          call mt_random(r6)
          r1 = 2.0*r1-1.0
          r2 = 2.0*r2-1.0
          r3 = 2.0*r3-1.0
          r4 = 2.0*r4-1.0
          r5 = 2.0*r5-1.0
          r6 = 2.0*r6-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 30
!x-px:
          xx1 = r1*sqrt(8.0)
          xx2 = r2*sqrt(8.0)
          xx3 = r3*sqrt(8.0)
          xx4 = r4*sqrt(8.0)

          numpts = numpts + 1
          if(numpts.gt.avgpts) exit

          !x-px:
!          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(1,numpts) = xmu1 + sig1*xx1/sq12
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*xx1/sq12+xx2)
          !y-py
!          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(3,numpts) = xmu3 + sig3*xx3/sq34
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*xx3/sq34+xx4)
          !z-pz
          this%Pts1(5,numpts) = xmu5 + (phi1-synangle)
          this%Pts1(6,numpts) = xmu6 + gamma0 - (1+e1*1.0e6/mccq)
        enddo

        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine WaterGamma_Dist

        subroutine KVGamma_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4
        double precision :: sq12,sq34
        double precision, dimension(2) :: x1
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,inipts,j,pid
!        integer seedarray(1)
        double precision :: t0,x11
        double precision :: alphaz,alphapz,lambdaz,lambdapz,pi,synangle,&
        kenergy,mccq,tmp1,tmp2,Emin,Emax,phimin,phimax,gamma0,fe,femax,fphi,&
        fphimax,e1,phi1,r1,r2,r3,r4,r5,xx1,xx2,xx3,xx4,twopi

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        alphaz = distparam(15)
        lambdaz = distparam(16)
        alphapz = distparam(17)
        lambdapz = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        print*,myid,x11

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
        close(12)
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq
        pi = 2*asin(1.0)
        synangle = synangle*pi/180.0

        Emin = 0.0
        !Emax = 6.9
        Emax = 7.5
        phimin = -1.8
        phimax = 0.6

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        numpts = 0
        print*,"kenergy: ",kenergy,synangle,avgpts,alphaz,lambdaz,&
               alphapz,lambdapz
        tmp1 = Emax - alphaz/lambdaz
        tmp2 = phimin + alphapz/lambdapz
        twopi = 4*asin(1.0)
        do
          ! rejection sample.
10        call mt_random(r1)
          r1 = Emin + (Emax-Emin)*r1
          femax = (Emax - tmp1)**alphaz*exp(-lambdaz*(Emax-tmp1))
          fe = (Emax - r1)**alphaz*exp(-lambdaz*(Emax-r1))/femax
          call mt_random(r3)
          if(r3.gt.fe) goto 10
          e1 = r1

20        call mt_random(r2)
          r2 = phimin + (phimax-phimin)*r2
          fphimax = (tmp2-phimin)**alphapz*exp(-lambdapz*(tmp2-phimin))
          fphi = (r2-phimin)**alphapz*exp(-lambdapz*(r2-phimin))/fphimax
          call mt_random(r4)
          if(r4.gt.fphi) goto 20
          phi1 = r2

          numpts = numpts + 1
          if(numpts.gt.avgpts) exit

          call mt_random(r1)
          call mt_random(r2)
          call mt_random(r3)
          r4 = sqrt(r1)
          r5 = sqrt(1.0-r1)
          r2 = r2*twopi
          r3 = r3*twopi
          xx1 = 2*r4*cos(r2)
          xx2 = 2*r4*sin(r2)
          xx3 = 2*r5*cos(r3)
          xx4 = 2*r5*sin(r3)

          !x-px:
!          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(1,numpts) = xmu1 + sig1*xx1/sq12
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*xx1/sq12+xx2)
          !y-py
!          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(3,numpts) = xmu3 + sig3*xx3/sq34
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*xx3/sq34+xx4)
          !z-pz
          this%Pts1(5,numpts) = xmu5 + (phi1-synangle)
          this%Pts1(6,numpts) = xmu6 + gamma0 - (1+e1*1.0e6/mccq)
        enddo

        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine KVGamma_Dist

        subroutine regen2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:6,i)
          sumx = sumx + Ptcl(1,i)
          sumx2 = sumx2 + Ptcl(1,i)*Ptcl(1,i)
          sumy = sumy + Ptcl(3,i)
          sumy2 = sumy2 + Ptcl(3,i)*Ptcl(3,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumy2
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        do j = 1, inipts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
          gammabet = sqrt(gamma*gamma-1.0)
          Ptcl(2,j) = (Ptcl(2,j)*gammabet)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gammabet)*pyscale + xmu4
!          betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
!          Ptcl(2,j) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
!          Ptcl(4,j) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j))
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j))
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j))
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j))
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j))
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j))
        enddo

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then

          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call mt_random(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call mt_random(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                !r = (2.0*r-1.0)*0.015*xmax !old one
                r = (2.0*r-1.0)*0.04*xmax !tt59
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(1,ii) = Ptcl(3,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax !old one
                r = (2.0*r-1.0)*0.04*pxmax !tt59
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(2,ii) = Ptcl(4,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                r = (2.0*r-1.0)*0.04*ymax !tt59
                !r = (2.0*r-1.0)*0.02*ymax !old one
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(3,ii) = Ptcl(1,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                r = (2.0*r-1.0)*0.04*pymax !tt59
                !r = (2.0*r-1.0)*0.02*pymax !old one
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(4,ii) = Ptcl(2,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                r = (2.0*r-1.0)*0.04*zmax !tt59
                !r = (2.0*r-1.0)*0.005*zmax !old one
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                r = (2.0*r-1.0)*0.04*pzmax !tt59
                !r = (2.0*r-1.0)*0.002*pzmax !old one
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo

          do j = 1,1
              ii = j+nset*(i-1)
              !tt35
              !this%Pts1(1,ii) = this%Pts1(2,ii)*1.2
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(2,ii)*1.2
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt37
              !this%Pts1(1,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(2,ii)*0.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt38
              !this%Pts1(1,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(3,ii) = this%Pts1(2,ii)*0.0
              !this%Pts1(3,ii) = 0.0
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.5
              !tt39
              !this%Pts1(1,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt32,tt40, nominal
              !tt41
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.0
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.0
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt42
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt43
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt44 nominal case with centroid offset for off-energy particles
              !this%Pts1(3,ii) = this%Pts1(3,ii) + 0.00025/xl
              !tt45
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 + 0.001/xl
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt46
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5 + 0.001/xl
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 + 0.001/xl
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt47
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5 + 0.001/xl
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 + 0.001/xl
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt48
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.1 + 0.001/xl
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.1 + 0.001/xl
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt49
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5 + 0.00025/xl
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt50 new partcl.data, different percentage off-energy particle
              !tt51,tt52 new partcl.data, Quad scan
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt53 new partcl.data
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.5
              call mt_random(r)
              this%Pts1(5,ii) = 2*pi*r - synangle*pi/180.0
              !this%Pts1(6,ii) = gamma0 - (1+0.075e6/mccq)
              call mt_random(r)
              this%Pts1(6,ii) = gamma0 - (1+(6.3+0.3*r)*1.0e6/mccq)
              !this%Pts1(6,ii) = gamma0 - (1+(6.3+0.4*r)*1.0e6/mccq)
              !this%Pts1(6,ii) = gamma0 - (1+(6.0+0.7*r)*1.0e6/mccq)
              !ttmatch using nominal output only, no off-energy particles
              !tt53 new partcl.data, no off-energy particle
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt54 new partcl.data, no off-energy particle
              !this%Pts1(1,ii) = this%Pts1(1,ii)*0.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(3,ii) = this%Pts1(3,ii)*0.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.5
              !tt55 new partcl.data, no off-energy particle,conserve RMS
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.17
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.17
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.17
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.17
              !tt56 new partcl.data, no off-energy particle,conserve RMS
              !this%Pts1(1,ii) = this%Pts1(1,ii)*0.56
              !this%Pts1(2,ii) = this%Pts1(2,ii)*0.56
              !this%Pts1(3,ii) = this%Pts1(3,ii)*0.56
              !this%Pts1(4,ii) = this%Pts1(4,ii)*0.56
              !tt57 test steering magnets, the other same as tt56
              !tt58 test steering magnets, no off-energy pt, matched beam.
              !tt59, 1% and 5% off-energy, big mismatch,new repopulation radii
!              this%Pts1(1,ii) = this%Pts1(1,ii)*2.00
!              this%Pts1(2,ii) = this%Pts1(2,ii)*2.00
!              this%Pts1(3,ii) = this%Pts1(3,ii)*2.00
!              this%Pts1(4,ii) = this%Pts1(4,ii)*2.00
              !tt60, 5% off-energy, big mismatch,new repopulation radii:400kev
          enddo

          sumx = sumx + this%Pts1(2,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        deallocate(Ptcl)

        do j = 1, this%Nptlocal
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = j
        enddo

        end subroutine regen2_Dist

        subroutine Gauss7_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(3) :: al0,ga0,epson0
        double precision :: t0

        call starttime_Timer(t0)

        call Waterbag_Dist(this,nparam,distparam,grid,0)
        call gammaepson_Dist(this,al0,ga0,epson0)
        call Distort_Dist(this,al0,ga0,epson0,grid)

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss7_Dist

        subroutine Distort_Dist(this,al0,ga0,epson0,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(3) :: al0,ga0,epson0
        integer :: avgpts
        integer :: i,ierr,nptot,j,pid
        double precision :: t0
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy,vtmp
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1,xx
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac,CSbeta,CSalpha,CSgamma,CSepson
        double precision, dimension(6) :: ptctmp
        double precision, dimension(15) :: tmplc,tmpgl
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,yrms
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0
        double precision ::sqx,sqpx,sqy,sqpy,sqz,sqpz
        double precision:: qmc,xl,xt,pi
        integer :: totnp,npy,npx,myidy,myid,myidx

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg
        nptot = this%Npt
        avgpts = this%Nptlocal

        ddx = -4.0e2 !for alphax > 0.0
        !ddx = -8.0e2 !for alphax > 0.0
        ddx = -4.0e2 !for alphax > 0.0
        !ddy = -3.0e2
        ddy = 1.0e7  !for alphay < 0.0
        ddy = 3.0e7  !for alphay < 0.0
        ddy = 2.0e7  !for alphay < 0.0 for vis
        !ddy = -4.0e2 !for alphay > 0.0 for visulization
        !ddy = 6.0e7  !for alphay < 0.0 for the purpose of display
        !ddz = -5.0e-3
        !ddz = 5.0e11 !for alphaz < 0.0
        ddz = 0.0
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0

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
        pi = 2*asin(1.0)

        do i = 1, avgpts
          ptctmp(1) = this%Pts1(1,i)
          vtmp = this%Pts1(1,i)
          ptctmp(2) = this%Pts1(2,i) + ddx*vtmp**3

          !for alpha > 0
          !ptctmp(3) = this%Pts1(3,i)
          !vtmp = this%Pts1(3,i)
          !ptctmp(4) = this%Pts1(4,i) + ddy*vtmp**3
          !for alpha < 0
          vtmp = this%Pts1(4,i)
          ptctmp(3) = this%Pts1(3,i) + ddy*vtmp**3
          ptctmp(4) = this%Pts1(4,i)

          !ptctmp(5) = this%Pts1(5,i)
          !vtmp = this%Pts1(5,i)
          !ptctmp(6) = this%Pts1(6,i) + ddz*vtmp**3
          vtmp = this%Pts1(6,i)
          ptctmp(5) = this%Pts1(5,i) + ddz*vtmp**3
          ptctmp(6) = this%Pts1(6,i)

          x0lc = x0lc + ptctmp(1)
          sqsum1local = sqsum1local + ptctmp(1)*ptctmp(1)
          xpxlocal = xpxlocal + ptctmp(1)*ptctmp(2)
          px0lc = px0lc + ptctmp(2)
          sqsum2local = sqsum2local + ptctmp(2)*ptctmp(2)
          y0lc = y0lc + ptctmp(3)
          sqsum3local = sqsum3local + ptctmp(3)*ptctmp(3)
          ypylocal = ypylocal + ptctmp(3)*ptctmp(4)
          py0lc = py0lc + ptctmp(4)
          sqsum4local = sqsum4local + ptctmp(4)*ptctmp(4)
          z0lc = z0lc + ptctmp(5)
          sqsum5local = sqsum5local + ptctmp(5)*ptctmp(5)
          zpzlocal = zpzlocal + ptctmp(5)*ptctmp(6)
          pz0lc = pz0lc + ptctmp(6)
          sqsum6local = sqsum6local + ptctmp(6)*ptctmp(6)
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

        call MPI_ALLREDUCE(tmplc,tmpgl,15,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)

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

        epsx2 = (sqsum1*sqsum2-xpx*xpx)
        epsy2 = (sqsum3*sqsum4-ypy*ypy)
        epsz2 = (sqsum5*sqsum6-zpz*zpz)
        epx = sqrt(max(epsx2,0.0d0))
        epy = sqrt(max(epsy2,0.0d0))
        epz = sqrt(max(epsz2,0.0d0))
        xrms = sqrt(sqsum1)
        yrms = sqrt(sqsum3)
        zrms = sqrt(sqsum5)

        CSalpha(1) = -xpx/epx
        CSalpha(2) = -ypy/epy
        CSalpha(3) = -zpz/epz
        CSbeta(1) = (xrms*xl)**2/(epx*xl)
        CSbeta(2) = (yrms*xl)**2/(epy*xl)
        CSbeta(3) = (zrms*xt)**2/(epz*qmc*xt)
        CSgamma(:) = (1.0 + CSalpha(:)*CSalpha(:))/CSbeta(:)*xl
        CSgamma(3) = (1.0 + CSalpha(3)*CSalpha(3))/CSbeta(3)/qmc*xt

        CSepson(1) = epx*xl
        CSepson(2) = epy*xl
        CSepson(3) = epz*qmc*xt

        do i = 1, 3
          frac(i) = sqrt(epson0(i)/CSepson(i))
          rr(1,i) = sqrt(CSgamma(i)/ga0(i))
          rr(2,i) = 1.0/rr(1,i)
          rr(3,i) = (CSalpha(i)-al0(i))/sqrt(CSgamma(i)*ga0(i))
        enddo

        do i = 1, avgpts
          vtmp = this%Pts1(1,i)
          this%Pts1(2,i) = this%Pts1(2,i) + ddx*vtmp**3
          !for alpha > 0
          !vtmp = this%Pts1(3,i)
          !this%Pts1(4,i) = this%Pts1(4,i) + ddy*vtmp**3
          !for alpha < 0
          vtmp = this%Pts1(4,i)
          this%Pts1(3,i) = this%Pts1(3,i) + ddy*vtmp**3
          !vtmp = this%Pts1(5,i)
          !this%Pts1(6,i) = this%Pts1(6,i) + ddz*vtmp**3
          vtmp = this%Pts1(6,i)
          this%Pts1(5,i) = this%Pts1(5,i) + ddz*vtmp**3
        enddo

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        do i = 1, avgpts
          this%Pts1(1,i) = this%Pts1(1,i)*factx
          this%Pts1(2,i) = this%Pts1(2,i)*factx
          this%Pts1(3,i) = this%Pts1(3,i)*facty
          this%Pts1(4,i) = this%Pts1(4,i)*facty
          this%Pts1(5,i) = this%Pts1(5,i)*factz
          this%Pts1(6,i) = this%Pts1(6,i)*factz
        enddo

        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1)
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2)
        r11z = rr(1,3)
        r22z = rr(2,3)
        r12z = rr(3,3)

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xxx
          this%Pts1(2,i) = yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xxx
          this%Pts1(4,i) = yyy
          xxx = this%Pts1(5,i)*r11z + this%Pts1(6,i)*r12z
          yyy = this%Pts1(6,i)*r22z
          this%Pts1(5,i) = xxx
          this%Pts1(6,i) = yyy
        enddo

        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        end subroutine Distort_Dist

        subroutine Regen7_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(3) :: al0,ga0,epson0
        double precision :: t0

        call starttime_Timer(t0)

        call Waterbag_Dist(this,nparam,distparam,grid,0)
        call gammaepson_Dist(this,al0,ga0,epson0)
        call regendstort_Dist(this,nparam,distparam,geom,grid,Flagbc)
        call DistReg_Dist(this,al0,ga0,epson0,grid)

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Regen7_Dist

        subroutine DistReg_Dist(this,al0,ga0,epson0,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, dimension(3) :: al0,ga0,epson0
        type (Pgrid2d), intent(in) :: grid
        integer :: avgpts
        integer :: i,ierr,nptot,j,pid
        double precision :: t0
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy,vtmp
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1,xx
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac,CSbeta,CSalpha,CSgamma,CSepson
        double precision, dimension(6) :: ptctmp
        double precision, dimension(15) :: tmplc,tmpgl
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,yrms
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0
        double precision ::sqx,sqpx,sqy,sqpy,sqz,sqpz
        double precision:: qmc,xl,xt,pi
        integer :: totnp,npx,npy,myid,myidx,myidy

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg
        nptot = this%Npt
        avgpts = this%Nptlocal

        ddx = 0.0
        ddy = 0.0
        ddz = 0.0
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0

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
        pi = 2*asin(1.0)

        do i = 1, avgpts
          ptctmp(1) = this%Pts1(1,i)
          vtmp = this%Pts1(1,i)
          ptctmp(2) = this%Pts1(2,i) + ddx*vtmp**3

          !for alpha > 0
          !ptctmp(3) = this%Pts1(3,i)
          !vtmp = this%Pts1(3,i)
          !ptctmp(4) = this%Pts1(4,i) + ddy*vtmp**3
          !for alpha < 0
          vtmp = this%Pts1(4,i)
          ptctmp(3) = this%Pts1(3,i) + ddy*vtmp**3
          ptctmp(4) = this%Pts1(4,i)

          !ptctmp(5) = this%Pts1(5,i)
          !vtmp = this%Pts1(5,i)
          !ptctmp(6) = this%Pts1(6,i) + ddz*vtmp**3
          vtmp = this%Pts1(6,i)
          ptctmp(5) = this%Pts1(5,i) + ddz*vtmp**3
          ptctmp(6) = this%Pts1(6,i)

          x0lc = x0lc + ptctmp(1)
          sqsum1local = sqsum1local + ptctmp(1)*ptctmp(1)
          xpxlocal = xpxlocal + ptctmp(1)*ptctmp(2)
          px0lc = px0lc + ptctmp(2)
          sqsum2local = sqsum2local + ptctmp(2)*ptctmp(2)
          y0lc = y0lc + ptctmp(3)
          sqsum3local = sqsum3local + ptctmp(3)*ptctmp(3)
          ypylocal = ypylocal + ptctmp(3)*ptctmp(4)
          py0lc = py0lc + ptctmp(4)
          sqsum4local = sqsum4local + ptctmp(4)*ptctmp(4)
          z0lc = z0lc + ptctmp(5)
          sqsum5local = sqsum5local + ptctmp(5)*ptctmp(5)
          zpzlocal = zpzlocal + ptctmp(5)*ptctmp(6)
          pz0lc = pz0lc + ptctmp(6)
          sqsum6local = sqsum6local + ptctmp(6)*ptctmp(6)
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

        call MPI_ALLREDUCE(tmplc,tmpgl,15,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)

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

        epsx2 = (sqsum1*sqsum2-xpx*xpx)
        epsy2 = (sqsum3*sqsum4-ypy*ypy)
        epsz2 = (sqsum5*sqsum6-zpz*zpz)
        epx = sqrt(max(epsx2,0.0d0))
        epy = sqrt(max(epsy2,0.0d0))
        epz = sqrt(max(epsz2,0.0d0))
        xrms = sqrt(sqsum1)
        yrms = sqrt(sqsum3)
        zrms = sqrt(sqsum5)

        CSalpha(1) = -xpx/epx
        CSalpha(2) = -ypy/epy
        CSalpha(3) = -zpz/epz
        CSbeta(1) = (xrms*xl)**2/(epx*xl)
        CSbeta(2) = (yrms*xl)**2/(epy*xl)
        CSbeta(3) = (zrms*xt)**2/(epz*qmc*xt)
        CSgamma(:) = (1.0 + CSalpha(:)*CSalpha(:))/CSbeta(:)*xl
        CSgamma(3) = (1.0 + CSalpha(3)*CSalpha(3))/CSbeta(3)/qmc*xt

        CSepson(1) = epx*xl
        CSepson(2) = epy*xl
        CSepson(3) = epz*qmc*xt

        do i = 1, 3
          frac(i) = sqrt(epson0(i)/CSepson(i))
          rr(1,i) = sqrt(CSgamma(i)/ga0(i))
          rr(2,i) = 1.0/rr(1,i)
          rr(3,i) = (CSalpha(i)-al0(i))/sqrt(CSgamma(i)*ga0(i))
        enddo

        do i = 1, avgpts
          vtmp = this%Pts1(1,i)
          this%Pts1(2,i) = this%Pts1(2,i) + ddx*vtmp**3
          !for alpha > 0
          !vtmp = this%Pts1(3,i)
          !this%Pts1(4,i) = this%Pts1(4,i) + ddy*vtmp**3
          !for alpha < 0
          vtmp = this%Pts1(4,i)
          this%Pts1(3,i) = this%Pts1(3,i) + ddy*vtmp**3
          !vtmp = this%Pts1(5,i)
          !this%Pts1(6,i) = this%Pts1(6,i) + ddz*vtmp**3
          vtmp = this%Pts1(6,i)
          this%Pts1(5,i) = this%Pts1(5,i) + ddz*vtmp**3
        enddo

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        do i = 1, avgpts
          this%Pts1(1,i) = this%Pts1(1,i)*factx
          this%Pts1(2,i) = this%Pts1(2,i)*factx
          this%Pts1(3,i) = this%Pts1(3,i)*facty
          this%Pts1(4,i) = this%Pts1(4,i)*facty
          this%Pts1(5,i) = this%Pts1(5,i)*factz
          this%Pts1(6,i) = this%Pts1(6,i)*factz
        enddo

        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1)
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2)
        r11z = rr(1,3)
        r22z = rr(2,3)
        r12z = rr(3,3)

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xxx
          this%Pts1(2,i) = yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xxx
          this%Pts1(4,i) = yyy
          xxx = this%Pts1(5,i)*r11z + this%Pts1(6,i)*r12z
          yyy = this%Pts1(6,i)*r22z
          this%Pts1(5,i) = xxx
          this%Pts1(6,i) = yyy
        enddo

        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        end subroutine DistReg_Dist

        subroutine regendstort_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,ikeep
        double precision, dimension(6) :: lcrange,a,tmptcl
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        ikeep = 0
        do i = 1, inipts
          read(12,*)tmptcl(1:6)
          if(abs(tmptcl(6)-kenergy/1.0e6).lt.1.0) then
            ikeep = ikeep + 1
            Ptcl(1:6,ikeep) = tmptcl(1:6)
            sumx = sumx + Ptcl(1,ikeep)
            sumx2 = sumx2 + Ptcl(1,ikeep)*Ptcl(1,ikeep)
            sumy = sumy + Ptcl(3,ikeep)
            sumy2 = sumy2 + Ptcl(3,ikeep)*Ptcl(3,ikeep)
          else
          endif
        enddo
        inipts = ikeep
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumy2,inipts
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        do j = 1, inipts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
          gammabet = sqrt(gamma*gamma-1.0)
          Ptcl(2,j) = (Ptcl(2,j)*gammabet)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gammabet)*pyscale + xmu4
!          betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
!          Ptcl(2,j) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
!          Ptcl(4,j) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j))
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j))
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j))
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j))
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j))
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j))
        enddo

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then

          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        deallocate(this%Pts1)
        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call mt_random(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call mt_random(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.04*xmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(1,ii) = Ptcl(3,i)+r
!              this%Pts1(1,ii) = Ptcl(1,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.04*pxmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(2,ii) = Ptcl(4,i)+r
!              this%Pts1(2,ii) = Ptcl(2,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                !r = (2.0*r-1.0)*0.04*ymax
                r = (2.0*r-1.0)*0.02*ymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(3,ii) = Ptcl(1,i)+r
!              this%Pts1(3,ii) = Ptcl(3,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                !r = (2.0*r-1.0)*0.04*pymax
                r = (2.0*r-1.0)*0.02*pymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(4,ii) = Ptcl(2,i)+r
!              this%Pts1(4,ii) = Ptcl(4,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                !r = (2.0*r-1.0)*0.04*zmax
                r = (2.0*r-1.0)*0.005*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                !r = (2.0*r-1.0)*0.04*pzmax
                r = (2.0*r-1.0)*0.002*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(2,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        deallocate(Ptcl)

        end subroutine regendstort_Dist

        subroutine GaussDouble_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp,pid
!        integer seedarray(1)
        double precision :: t0,x11,Xfrac,Yfrac
        double precision :: sig1a,sig2a,sig3a,sig4a,sig5a,sig6a
        double precision :: sig1b,sig2b,sig3b,sig4b,sig5b,sig6b
        integer :: num1,num2

        call starttime_Timer(t0)

        !Xfrac = 5.0
        !Yfrac = 0.01
        !Xfrac = 4.0
        !Yfrac = 0.0
        !Xfrac = 4.0
        !Yfrac = 0.001
        !Xfrac = 4.0
        !Yfrac = 0.0002
        !Xfrac = 4.0
        !Yfrac = 0.01
        !Xfrac = 4.0
        !Yfrac = 0.05
        !Xfrac = 4.0
        !Yfrac = 0.005
        !Xfrac = 2.0
        !Yfrac = 0.01
        Xfrac = 2.0
        Yfrac = 0.05

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(x11)
        print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        if(mod(avgpts,10).ne.0) then
          print*,"The number of particles has to be an integer multiple of 10Nprocs"
          stop
        endif

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sig1a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig1
        sig2a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig2
        sig3a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig3
        sig4a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig4
        sig5a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig5
        sig6a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig6

        sig1b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig1
        sig2b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig2
        sig3b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig3
        sig4b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig4
        sig5b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig5
        sig6b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig6

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)

        intvsamp = 10000 !used in the previous simulation
!        intvsamp = avgpts
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

!        num1 = intvsamp*(1-Yfrac) + 1
        num2 = intvsamp*(Yfrac/(1.0+Yfrac))
        num1 = intvsamp - num2

        print*,"num1: ",num1,num2
        print*,"sig1-6a: ",sig1a,sig2a,sig3a,sig4a,sig5a,sig6a
        print*,"sig1-6b: ",sig1b,sig2b,sig3b,sig4b,sig5b,sig6b

        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, num1
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1a*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2a*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3a*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4a*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5a*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6a*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
          do k = num1+1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1b*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2b*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3b*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4b*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5b*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6b*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
        enddo

        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)
        print*,"after Gauss Double: ",avgpts

        end subroutine GaussDouble_Dist

        ! sample the particles with intial distribution
        ! using rejection method for multi-charge state.
        subroutine WaterbagMC_Dist(this,nparam,distparam,grid,flagalloc,nchrg,&
                                   nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc,nchrg
        double precision, dimension(nparam) :: distparam
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer, dimension(nchrg) :: nptlist, nptlistsum
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx,restnpt
        integer :: avgpts,numpts,isamz,isamy, tmpnpt
        integer :: myid,myidx,myidy,iran,i,j,intvsamp,ierr
        integer, allocatable, dimension(:) :: shead, scount
        double precision :: t0,x11
        double precision, allocatable, dimension(:) :: ranum6
        double precision, allocatable, dimension(:,:) :: bpt9d

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call mt_random(x11)
        !print*,"x11: ",myid,x11

        avgpts = this%Npt/totnp
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
        !avgpts = this%Npt/(npx*npy)
        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs"
        !  stop
        !endif

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

      if (myid.eq.0) then
        allocate(bpt9d(9,this%Npt))
        numpts = 0
        isamz = 0
        isamy = 0
        !intvsamp = avgpts
        intvsamp = this%Npt
        allocate(ranum6(6*intvsamp))

        tmpnpt = 0
        do i = 1, nchrg
            tmpnpt = tmpnpt + nptlist(i)
            nptlistsum(i) = tmpnpt
        enddo


        do
          ! rejection sample.
10        continue
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call mt_random(ranum6)
          endif
          iran = 6*mod(isamz-1,intvsamp)
          r1 = 2.0*ranum6(iran+1)-1.0
          r2 = 2.0*ranum6(iran+2)-1.0
          r3 = 2.0*ranum6(iran+3)-1.0
          r4 = 2.0*ranum6(iran+4)-1.0
          r5 = 2.0*ranum6(iran+5)-1.0
          r6 = 2.0*ranum6(iran+6)-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 10
          isamy = isamy + 1
          numpts = numpts + 1
          !if(numpts.gt.avgpts) exit
          if(numpts.gt.this%Npt) exit
!x-px:
          x1 = r1*sqrt(8.0)
          x2 = r2*sqrt(8.0)
          !Correct transformation.
          !this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          !this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          bpt9d(1,numpts) = xmu1 + sig1*x1/rootx
          bpt9d(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          x3 = r3*sqrt(8.0)
          x4 = r4*sqrt(8.0)
          !correct transformation
          !this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          !this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          bpt9d(3,numpts) = xmu3 + sig3*x3/rooty
          bpt9d(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          x5 = r5*sqrt(8.0)
          x6 = r6*sqrt(8.0)
          !correct transformation
          !this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          !this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          bpt9d(5,numpts) = xmu5 + sig5*x5/rootz
          bpt9d(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale

          do j=1, nchrg
              if (numpts.le.nptlistsum(j)) then
                bpt9d(7,numpts) = qmcclist(j)
                bpt9d(8,numpts) = currlist(j)/Scfreq/nptlist(j)*qmcclist(j)/abs(qmcclist(j))
                bpt9d(9,numpts) = numpts
                exit
            endif
          enddo

        enddo

        deallocate(ranum6)
      endif

        call MPI_SCATTERV(bpt9d(1,1), scount, shead, MPI_DOUBLE_PRECISION,&
                          this%Pts1(1,1), 9*avgpts, MPI_DOUBLE_PRECISION, &
                          0, MPI_COMM_WORLD, ierr)

        this%Nptlocal = avgpts

        !ii = 0
        !do i = 1, nchrg
        !  kk = nptlist(i)/totnp
        !  do j = 1, kk
        !    ii = ii + 1
        !    this%Pts1(7,ii) = qmcclist(i)
        !    this%Pts1(8,ii) = currlist(i)/Scfreq/nptlist(i)*qmcclist(i)/abs(qmcclist(i))
        !  enddo
        !enddo

!        do j = 1, avgpts
!          pid = j + myid*totnp
!          this%Pts1(7,j) = this%Charge/this%mass
!          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
!          this%Pts1(9,j) = pid
!        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine WaterbagMC_Dist

        subroutine GaussMC_Dist(this,nparam,distparam,grid,flagalloc,nchrg,&
                                   nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc,nchrg
        double precision, dimension(nparam) :: distparam
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer, dimension(nchrg) :: nptlist, nptlistsum
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3
        integer :: totnp,npy,npx,restnpt
        integer :: avgpts,numpts, tmpnpt
        integer :: myid,myidx,myidy,i,j,intvsamp,ierr
        integer, allocatable, dimension(:) :: shead, scount
        double precision :: t0,x11
        double precision, allocatable, dimension(:,:) :: bpt9d

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call mt_random(x11)
!        print*,myid,x11


        avgpts = this%Npt/totnp
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
!        avgpts = this%Npt/(npx*npy)
!        if(mod(avgpts,10).ne.0) then
!          print*,"The number of particles has to be an integer multiple of 10Nprocs"
!          stop
!        endif

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

        if (myid.eq.0) then
          allocate(bpt9d(9,this%Npt))
          intvsamp = this%Npt

          tmpnpt = 0
          do i = 1, nchrg
            tmpnpt = tmpnpt + nptlist(i)
            nptlistsum(i) = tmpnpt
          enddo

          allocate(x1(2,intvsamp))
          allocate(x2(2,intvsamp))
          allocate(x3(2,intvsamp))

          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do i = 1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            !this%Pts1(1,i) = xmu1 + sig1*x1(1,i)/sq12
            !this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,i)/sq12+x1(2,i))
            bpt9d(1,i) = xmu1 + sig1*x1(1,i)/sq12
            bpt9d(2,i) = xmu2 + sig2*(-muxpx*x1(1,i)/sq12+x1(2,i))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            !this%Pts1(3,i) = xmu3 + sig3*x2(1,i)/sq34
            !this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,i)/sq34+x2(2,i))
            bpt9d(3,i) = xmu3 + sig3*x2(1,i)/sq34
            bpt9d(4,i) = xmu4 + sig4*(-muypy*x2(1,i)/sq34+x2(2,i))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            !this%Pts1(5,i) = xmu5 + sig5*x3(1,i)/sq56
            !this%Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,i)/sq56+x3(2,i))
            bpt9d(5,i) = xmu5 + sig5*x3(1,i)/sq56
            bpt9d(6,i) = xmu6 + sig6*(-muzpz*x3(1,i)/sq56+x3(2,i))

            do j=1, nchrg
              if (i.le.nptlistsum(j)) then
                bpt9d(7,i) = qmcclist(j)
                bpt9d(8,i) = currlist(j)/Scfreq/nptlist(j)*qmcclist(j)/abs(qmcclist(j))
                bpt9d(9,i) = i
                exit
              endif
            enddo
          enddo

          deallocate(x1)
          deallocate(x2)
          deallocate(x3)
        endif

        call MPI_SCATTERV(bpt9d(1,1), scount, shead, MPI_DOUBLE_PRECISION,&
                          this%Pts1(1,1), 9*avgpts, MPI_DOUBLE_PRECISION, &
                          0, MPI_COMM_WORLD, ierr)

        this%Nptlocal = avgpts

        !ii = 0
        !do i = 1, nchrg
        !  kk = nptlist(i)/totnp
        !  do j = 1, kk
        !    ii = ii + 1
        !    this%Pts1(7,ii) = qmcclist(i)
        !    this%Pts1(8,ii) = currlist(i)/Scfreq/nptlist(i)*qmcclist(i)/abs(qmcclist(i))
        !  enddo
        !enddo

        !do j = 1, avgpts
        !  pid = j + myid*totnp
!          this%Pts1(7,j) = this%Charge/this%mass
!          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
        !  this%Pts1(9,j) = pid
        !enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine GaussMC_Dist

        subroutine regenMCold_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabetz
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi,tmp3
        real*8 :: gammabet

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:7,i)
          sumx = sumx + Ptcl(1,i)
          sumx2 = sumx2 + Ptcl(1,i)*Ptcl(1,i)
          sumy = sumy + Ptcl(3,i)
          sumy2 = sumy2 + Ptcl(3,i)*Ptcl(3,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumy2,sumx,sumy
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        print*,"gamma0: ",gamma0
!once used for LANL and MUS
!        do j = 1, inipts
!!          print*,"j: ",j,Ptcl(1,j)
!          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
!          tmp3 = Ptcl(3,j)
!          Ptcl(3,j) = (Ptcl(2,j)/100.0/xl)*yscale + xmu3
!          !this is for MSU input only
!          !gamma = Ptcl(6,j)*(gamma0-1)/100+gamma0
!          !this is for Lanl input only
!          gamma = Ptcl(6,j)*1.0e6/mccq+1.0
!          gammabetz = sqrt((gamma**2-1.0)/(1+(Ptcl(4,j)/1000)**2 &
!             +(Ptcl(5,j)/1000))**2)
!          Ptcl(2,j) = (Ptcl(4,j)*gammabetz/1000)*pxscale + xmu2
!          Ptcl(4,j) = (Ptcl(5,j)*gammabetz/1000)*pyscale + xmu4
!          !unit in rad
!          Ptcl(5,j) = (tmp3*2.0*asin(1.0)/180.0)*zscale+xmu5
!          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j))
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j))
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j))
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j))
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j))
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j))
!        enddo

!for ANL input
        do j = 1, inipts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
!for ANL data
!          gammabet = sqrt(gamma*gamma-1.0)
!          Ptcl(2,j) = (Ptcl(2,j)*gammabet)*pxscale + xmu2
!          Ptcl(4,j) = (Ptcl(4,j)*gammabet)*pyscale + xmu4
!for LANL data from readrfq.f
          betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
          Ptcl(2,j) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j))
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j))
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j))
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j))
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j))
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j))
        enddo

        print *,"pass 1: ",inipts

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then

          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call mt_random(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call mt_random(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.04*xmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(1,ii) = Ptcl(3,i)+r
              this%Pts1(1,ii) = Ptcl(1,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.04*pxmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(2,ii) = Ptcl(4,i)+r
              this%Pts1(2,ii) = Ptcl(2,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                !r = (2.0*r-1.0)*0.04*ymax
                r = (2.0*r-1.0)*0.02*ymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(3,ii) = Ptcl(1,i)+r
              this%Pts1(3,ii) = Ptcl(3,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                !r = (2.0*r-1.0)*0.04*pymax
                r = (2.0*r-1.0)*0.02*pymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(4,ii) = Ptcl(2,i)+r
              this%Pts1(4,ii) = Ptcl(4,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                !r = (2.0*r-1.0)*0.04*zmax
                r = (2.0*r-1.0)*0.005*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                !r = (2.0*r-1.0)*0.04*pzmax
                r = (2.0*r-1.0)*0.002*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(2,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        do j = 1, this%Nptlocal
          pid = j + myid*totnp
          !for ANL input only
          this%Pts1(7,j) = Ptcl(7,j)
          !this%Pts1(7,j) = Ptcl(7,j)/mccq
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        deallocate(Ptcl)

        end subroutine regenMCold_Dist

        subroutine regenMSU_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabetz
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi,tmp3
        real*8 :: gammabet

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:7,i)
          sumx = sumx + Ptcl(1,i)
          sumx2 = sumx2 + Ptcl(1,i)*Ptcl(1,i)
          sumy = sumy + Ptcl(3,i)
          sumy2 = sumy2 + Ptcl(3,i)*Ptcl(3,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumy2,sumx,sumy
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        print*,"gamma0: ",gamma0
        do j = 1, inipts
!          print*,"j: ",j,Ptcl(1,j)
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          tmp3 = Ptcl(3,j)
          Ptcl(3,j) = (Ptcl(2,j)/100.0/xl)*yscale + xmu3
          !this is for MSU input only
          gamma = Ptcl(6,j)*(gamma0-1)/100+gamma0
          gammabetz = sqrt((gamma**2-1.0)/(1+(Ptcl(4,j)/1000)**2 &
             +(Ptcl(5,j)/1000))**2)
          Ptcl(2,j) = (Ptcl(4,j)*gammabetz/1000)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(5,j)*gammabetz/1000)*pyscale + xmu4
          !unit in rad
          Ptcl(5,j) = (tmp3*2.0*asin(1.0)/180.0)*zscale+xmu5
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j))
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j))
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j))
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j))
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j))
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j))
        enddo

        print *,"pass 1: ",inipts

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then

          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call mt_random(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call mt_random(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.04*xmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(1,ii) = Ptcl(3,i)+r
              this%Pts1(1,ii) = Ptcl(1,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.04*pxmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(2,ii) = Ptcl(4,i)+r
              this%Pts1(2,ii) = Ptcl(2,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                !r = (2.0*r-1.0)*0.04*ymax
                r = (2.0*r-1.0)*0.02*ymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(3,ii) = Ptcl(1,i)+r
              this%Pts1(3,ii) = Ptcl(3,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                !r = (2.0*r-1.0)*0.04*pymax
                r = (2.0*r-1.0)*0.02*pymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(4,ii) = Ptcl(2,i)+r
              this%Pts1(4,ii) = Ptcl(4,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                !r = (2.0*r-1.0)*0.04*zmax
                r = (2.0*r-1.0)*0.005*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                !r = (2.0*r-1.0)*0.04*pzmax
                r = (2.0*r-1.0)*0.002*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(2,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        do j = 1, this%Nptlocal
          pid = j + myid*totnp
          this%Pts1(7,j) = Ptcl(7,j)/mccq
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        deallocate(Ptcl)

        end subroutine regenMSU_Dist

        subroutine regenANL_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabetz
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi,tmp3
        real*8 :: gammabet

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:7,i)
          sumx = sumx + Ptcl(1,i)
          sumx2 = sumx2 + Ptcl(1,i)*Ptcl(1,i)
          sumy = sumy + Ptcl(3,i)
          sumy2 = sumy2 + Ptcl(3,i)*Ptcl(3,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumy2,sumx,sumy
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        print*,"gamma0: ",gamma0

!for ANL input
        do j = 1, inipts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
          gammabet = sqrt(gamma*gamma-1.0)
          Ptcl(2,j) = (Ptcl(2,j)*gammabet)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gammabet)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j))
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j))
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j))
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j))
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j))
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j))
        enddo

        print *,"pass 1: ",inipts

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then

          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call mt_random(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call mt_random(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.04*xmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(1,ii) = Ptcl(3,i)+r
              this%Pts1(1,ii) = Ptcl(1,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.04*pxmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(2,ii) = Ptcl(4,i)+r
              this%Pts1(2,ii) = Ptcl(2,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                !r = (2.0*r-1.0)*0.04*ymax
                r = (2.0*r-1.0)*0.02*ymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(3,ii) = Ptcl(1,i)+r
              this%Pts1(3,ii) = Ptcl(3,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                !r = (2.0*r-1.0)*0.04*pymax
                r = (2.0*r-1.0)*0.02*pymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(4,ii) = Ptcl(2,i)+r
              this%Pts1(4,ii) = Ptcl(4,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                !r = (2.0*r-1.0)*0.04*zmax
                r = (2.0*r-1.0)*0.005*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                !r = (2.0*r-1.0)*0.04*pzmax
                r = (2.0*r-1.0)*0.002*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(2,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        do j = 1, this%Nptlocal
          pid = j + myid*totnp
          !for ANL input only
          this%Pts1(7,j) = Ptcl(7,j)
          !this%Pts1(7,j) = Ptcl(7,j)/mccq
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        deallocate(Ptcl)

        end subroutine regenANL_Dist

        subroutine regenLANL_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabetz
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi,tmp3
        real*8 :: gammabet

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call mt_random(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        read(12,*)inipts,kenergy,synangle
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:7,i)
          sumx = sumx + Ptcl(1,i)
          sumx2 = sumx2 + Ptcl(1,i)*Ptcl(1,i)
          sumy = sumy + Ptcl(3,i)
          sumy2 = sumy2 + Ptcl(3,i)*Ptcl(3,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumy2,sumx,sumy
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        print*,"gamma0: ",gamma0

        do j = 1, inipts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
          betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
          Ptcl(2,j) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j))
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j))
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j))
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j))
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j))
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j))
        enddo

        print *,"pass 1: ",inipts

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then

          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call mt_random(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call mt_random(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.04*xmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(1,ii) = Ptcl(3,i)+r
              this%Pts1(1,ii) = Ptcl(1,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.04*pxmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(2,ii) = Ptcl(4,i)+r
              this%Pts1(2,ii) = Ptcl(2,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                !r = (2.0*r-1.0)*0.04*ymax
                r = (2.0*r-1.0)*0.02*ymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(3,ii) = Ptcl(1,i)+r
              this%Pts1(3,ii) = Ptcl(3,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                !r = (2.0*r-1.0)*0.04*pymax
                r = (2.0*r-1.0)*0.02*pymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(4,ii) = Ptcl(2,i)+r
              this%Pts1(4,ii) = Ptcl(4,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                !r = (2.0*r-1.0)*0.04*zmax
                r = (2.0*r-1.0)*0.005*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call mt_random(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                !r = (2.0*r-1.0)*0.04*pzmax
                r = (2.0*r-1.0)*0.002*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(2,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        do j = 1, this%Nptlocal
          pid = j + myid*totnp
          this%Pts1(7,j) = Ptcl(7,j)
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        deallocate(Ptcl)

        end subroutine regenLANL_Dist

        subroutine readin_Dist(this,nchrg,nptlist,qmcclist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(inout) :: nchrg
        double precision, dimension(nchrg), intent(inout) :: qmcclist
        integer, dimension(nchrg), intent(inout) :: nptlist
        integer :: i,j,avgpts,myid,totnp,ierr,npsum,restnpt, itt
        integer, allocatable, dimension(:) :: shead, scount
        double precision :: x1,x2
        double precision, allocatable, dimension(:,:) :: bpt9d

        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,totnp,ierr)

        if (myid.eq.0) then
          open(unit=12,file=trim(particle_dname)//trim(particle_fname),status='old',iostat=iodistfile)
        endif

        call MPI_BCAST(iodistfile,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (iodistfile.ne.0) then
            if(myid.eq.0) print*, 'Distribution.f90 open(): Distribution file (',&
                           trim(particle_dname)//trim(particle_fname),&
                           ') does not found.'
            return
        endif

        if (myid.eq.0) read(12,*)npsum,x1,x2

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

        if (myid.eq.0) then
            allocate(bpt9d(9,npsum))
            do j = 1, npsum
              read(12,*) bpt9d(1:9, j)
            enddo
            close(12)

            qmcclist = 0.0
            nptlist = 0
            nchrg = 0
            do i = 1, npsum
                if (any(bpt9d(7,i) .eq. qmcclist)) then
                    do j = 1, nchrg
                        if (bpt9d(7,i) .eq. qmcclist(j)) then
                            nptlist(j) = nptlist(j) + 1
                            exit
                        endif
                    enddo
                else
                    nchrg = nchrg + 1
                    qmcclist(nchrg) = bpt9d(7,i)
                    nptlist(nchrg) = 1
                endif
            enddo
        endif

        call MPI_SCATTERV(bpt9d(1,1), scount, shead, MPI_DOUBLE_PRECISION,&
                          this%Pts1(1,1), 9*avgpts, MPI_DOUBLE_PRECISION, &
                          0, MPI_COMM_WORLD, ierr)

        this%Nptlocal = avgpts

        call MPI_BCAST(nchrg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qmcclist,100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nptlist,100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        end subroutine readin_Dist

      end module Distributionclass
