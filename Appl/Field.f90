!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! FieldQuantclass: 3D field quantity class in Field module of APPLICATION
!                 layer.
! Version: 2.0
! Author: Ji Qiang, Robert Ryne, LBNL, 1/05/02
! Description: This class defines a 3-D field quantity in the accelerator.
!              The field quantity can be updated at each step.
! Comments:
!----------------------------------------------------------------
      module FieldQuantclass
        use Timerclass
        use CompDomclass
        use Pgrid2dclass
        use FFTclass
        use Besselclass
        use Transposeclass
        type FieldQuant
!          private
          !# of mesh points in x and y directions.
          integer :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal
          !# field quantity array.
          double precision, pointer, dimension(:,:,:) :: FieldQ
        end type FieldQuant
      contains
        !Initialize field class.
        subroutine construct_FieldQuant(this,innx,inny,innz,geom,grid)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz
        integer :: myid, myidx, myidy, nptot,nproccol,nprocrow
        integer, allocatable, dimension(:,:,:) :: LocalTable

        call getsize_Pgrid2d(grid,nptot,nproccol,nprocrow)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        allocate(LocalTable(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(nproccol.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprocrow.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif

        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal))
        this%FieldQ = 0.0

        deallocate(LocalTable)

        end subroutine construct_FieldQuant

        ! set field quantity.
        subroutine set_FieldQuant(this,innx,inny,innz,geom,grid,&
                                  nprx,npry)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz, nprx, npry
        integer :: myid, myidx, myidy
        integer, dimension(2,0:nprx-1,0:npry-1)::LocalTable

        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(npry.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprx.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif
        this%Nxlocal = innx

        deallocate(this%FieldQ)
        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal))
        this%FieldQ = 0.0

        end subroutine set_FieldQuant

!----------------------------------------------------------------------
! update potential (solving Possion's equation) with 3D isolated
! boundary conditions.
        subroutine update3O_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        type (FieldQuant), intent(inout) :: this
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc2,nypzlc2
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)

        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1)
        inyglb = glmshnm(2)
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              rho(i,j,k) = source(i,j+jadd,k+kadd)
            enddo
          enddo
        enddo

        ! +1 is from the real to complex fft.
        nsxy1 = (inxglb+1)/nproccol
        nsxy2 = (inxglb+1) - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo

        nsyz1 = 2*inyglb/nprocrow
        nsyz2 = 2*inyglb - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo

        nxpylc2 = xpystable(myidy)
        nypzlc2 = ypzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        ! Open boundary conditions!
        call openBC3D(innx,inny,innz,rho,hx,hy,hz, &
        nxpylc2,nypzlc2,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,&
        inxglb,inyglb,inzglb)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)
            enddo
          enddo
        enddo

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update3O_FieldQuant

        ! Solving Poisson's equation with open BCs.
        subroutine openBC3D(innx,inny,innz,rho,hx,hy,hz,&
           nxpylc2,nypzlc2,myidz,myidy,npz,npy,commrow,commcol,comm2d, &
           pztable,pytable,ypzstable,xpystable,inxglb,&
           inyglb,inzglb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc2,nypzlc2
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz
        double precision :: scalex,scaley,scalez
        integer :: i,j,k,n1,n2,n3,nylc22,nzlc22
        double complex, allocatable, dimension(:,:,:) :: rho2out
        double complex, allocatable, dimension(:,:,:) :: grn
        integer :: ginny,ginnz,ierr

        n1 = 2*inxglb
        n2 = 2*inyglb
        n3 = 2*inzglb

        nylc22 = nxpylc2
        nzlc22 = nypzlc2

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        allocate(rho2out(n3,nylc22,nzlc22))

        call fft3d1_FFT(n1,n2,n3,innz,inny,nylc22,nzlc22,&
                1,scalex,scaley,scalez,rho,ypzstable,pztable,&
        xpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy, &
        rho2out)

        !c compute FFT of the Green function on the grid:
        ! here the +1 is from the unsymmetry of green function
        ! on double-sized grid.
        if(myidz.eq.(npz-1)) then
           ginnz = innz + 1
        else
           ginnz = innz
        endif
        if(myidy.eq.(npy-1)) then
           ginny = inny + 1
        else
           ginny = inny
        endif
        allocate(grn(n3,nylc22,nzlc22))
        call greenf1(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
               hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
                  ypzstable,pztable,xpystable,pytable,grn)

        ! multiply transformed charge density and transformed Green
        ! function:
        do k = 1, nzlc22
          do j = 1, nylc22
            do i = 1, n3
              rho2out(i,j,k) = rho2out(i,j,k)*grn(i,j,k)
            enddo
          enddo
        enddo

        deallocate(grn)

        ! inverse FFT:
        scalex = 1.0/float(n1)
        scaley = 1.0/float(n2)
        scalez = 1.0/float(n3)
        call invfft3d1_FFT(n3,n2,n1,nylc22,nzlc22,inny,innz,&
               -1,scalex,scaley,scalez,rho2out,pztable,ypzstable,&
        pytable,xpystable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)

        deallocate(rho2out)

        end subroutine openBC3D

        ! green function for extended array.
        subroutine greenf1(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (2*nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(2*nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        gxrtable = xrtable
        gxrtable(npx-1) = xrtable(npx-1) + 1
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo

        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, nx+1
              jj = j + nblocky
              kk = k + nblockx
                kkk = kk - 1
                jjj = jj - 1
                iii = i - 1
              if((iii*iii+jjj*jjj+kkk*kkk) .ne. 0) then
               grn(i,j,k)=1.0/sqrt((hx*iii)**2+(hy*jjj)**2+(hz*kkk)**2)
              endif
            enddo
          enddo
        enddo
        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along y and z dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3/2+1,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2+1)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3/2+1
              tmp3(i,j) = x1(i,j,k)
            enddo
            do i = n3/2+2,n3
              tmp3(i,j) = x1(n3-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf1

!--------------------------------------------------------------------
! Field potential from 2D transverse open, 1D longidutinal periodic
        !update potential (solving Possion's equation) with isolated
        ! boundary conditions.
        subroutine update2O1P_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        type (FieldQuant), intent(inout) :: this
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc2,nypzlc2
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)

        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        pztable(nprocrow-1) = pztable(nprocrow-1) - 1 !periodic BC.
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1)
        inyglb = glmshnm(2)
        inzglb = glmshnm(3)-1  !periodic BC.

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              rho(i,j,k) = source(i,j+jadd,k+kadd)
            enddo
          enddo
        enddo

        ! +1 is from the double precision to complex fft.
        nsxy1 = (inxglb+1)/nproccol
        nsxy2 = (inxglb+1) - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo

        nsyz1 = 2*inyglb/nprocrow
        nsyz2 = 2*inyglb - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo

        nxpylc2 = xpystable(myidy)
        nypzlc2 = ypzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        ! Open boundary conditions!
        call open2perd1(innx,inny,innz,rho,hx,hy,hz, &
        nxpylc2,nypzlc2,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,&
        inxglb,inyglb,inzglb)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)
            enddo
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
        !if(myidz.eq.(nprocrow-1)) then
        !  do j = 1, inny
        !    do i = 1, innx
        !      this%FieldQ(i,j+jadd,innz+1+kadd) = rho(i,j,1+kadd)
        !    enddo
        !  enddo
        !endif

        !open(8,file="phir",status="unknown")
        !do i = 1, innx
        !    write(8,1000)(i-1)*hx,this%FieldQ(i,1,1),&
        !    this%FieldQ(i,1,innz/2),this%FieldQ(i,1,innz)
        !enddo
        !close(8)
        !open(9,file="phithe",status="unknown")
        !do j = 1, inny+1
        !  write(9,1000)(j-1)*hy,this%FieldQ(2,j,2),&
        !  this%FieldQ(3,j,2),this%FieldQ(4,j,2)
        !enddo
        !close(9)
        !open(10,file="phiz",status="unknown")
        !do k = 1, innz+1
        !  write(10,1000)(k-1)*hz,this%FieldQ(1,1,k),&
        !  this%FieldQ(innx/2,1,k),this%FieldQ(innx,1,k)
        !enddo
        !close(10)

1000    format(4(1x,e13.6))

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update2O1P_FieldQuant

        ! Solving Poisson's equation with open BCs.
        subroutine open2perd1(innx,inny,innz,rho,hx,hy,hz,&
           nxpylc2,nypzlc2,myidz,myidy,npz,npy,commrow,commcol,comm2d, &
           pztable,pytable,ypzstable,xpystable,inxglb,&
           inyglb,inzglb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc2,nypzlc2
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz
        !double precision, dimension(innx,inny,innz) :: rhosum
        double precision :: scalex,scaley,scalez,sumtest,sumtest2
        integer :: i,j,k,n1,n2,n3,nylc22,nzlc22
        double complex, allocatable, dimension(:,:,:) :: rho2out
        double complex, allocatable, dimension(:,:,:) :: grn
        integer :: ginny,ginnz,ierr

!        rhosum = rho

        n1 = 2*inxglb
        n2 = 2*inyglb
        n3 = inzglb    ! periodic boundary condition.

        nylc22 = nxpylc2
        nzlc22 = nypzlc2

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        allocate(rho2out(n3,nylc22,nzlc22))

        call fft3d2_FFT(n1,n2,n3,innz,inny,nylc22,nzlc22,&
                1,scalex,scaley,scalez,rho,ypzstable,pztable,&
        xpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy, &
        rho2out)

        !c compute FFT of the Green function on the grid:
        ! here the +1 is from the unsymmetry of green function
        ! on double-sized grid.
        ginnz = innz  !periodic boundary condition.
        if(myidy.eq.(npy-1)) then
           ginny = inny + 1
        else
           ginny = inny
        endif
        allocate(grn(n3,nylc22,nzlc22))
        call greenf2(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
               hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
                  ypzstable,pztable,xpystable,pytable,grn)

        ! multiply transformed charge density and transformed Green
        ! function:
        do k = 1, nzlc22
          do j = 1, nylc22
            do i = 1, n3
              rho2out(i,j,k) = rho2out(i,j,k)*grn(i,j,k)
            enddo
          enddo
        enddo

        deallocate(grn)

        ! inverse FFT:
        scalex = 1.0/float(n1)
        scaley = 1.0/float(n2)
        scalez = 1.0/float(n3)
        call invfft3d2_FFT(n3,n2,n1,nylc22,nzlc22,inny,innz,&
               -1,scalex,scaley,scalez,rho2out,pztable,ypzstable,&
        pytable,xpystable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)

!        sumtest = 0.0
!        sumtest2 = 0.0
!        do k = 1, innz
!          do j = 1, inny
!            do i = 1, innx
!              sumtest2 = sumtest2 + rhosum(i,j,k)
!              sumtest = sumtest + abs(rho(i,j,k)-rhosum(i,j,k))
!            enddo
!          enddo
!        enddo
!        print*,"sumtest-fft: ",sumtest,sumtest2

        deallocate(rho2out)

        end subroutine open2perd1

        ! green function for extended array.
        subroutine greenf2(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr,kz
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0
        double precision :: btlmda,sumtmp,tmptmp

        call starttime_Timer(t0)

        btlmda = nz*hz
!        print*,"btlmda: ",btlmda
        gxrtable = xrtable  !periodic BC
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo
        sumtmp = 0.0
        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, nx+1
              jj = j + nblocky
              kk = k + nblockx
              if(kk.le.nz/2) then
                kkk = kk - 1
              else
                kkk = kk - nz - 1
              endif
              jjj = jj - 1
              iii = i - 1
              grn(i,j,k) = 0.0
              do kz = -8, 8
                tmptmp=(hx*iii)**2+(hy*jjj)**2+(hz*kkk+kz*btlmda)**2
                if(tmptmp .gt. 1.0e-16) then
                 sumtmp=1.0/sqrt(tmptmp)
                endif
                grn(i,j,k) = grn(i,j,k) + sumtmp
              enddo
            enddo
          enddo
        enddo
        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0
        n1 = 2*nx
        n2 = 2*ny
        n3 = nz   !periodic boundary condition.
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along x and y dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along x dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3
              tmp3(i,j) = x1(i,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf2
!-------------------------------------------------------------------------

        subroutine update3_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr,&
          besscoef,bessnorm,gml,modth,nmod)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        type (FieldQuant), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        double precision, dimension(nxlc,nxlc,nmod), intent(in) :: besscoef
        double precision, dimension(nxlc,nmod), intent(in) :: bessnorm,gml
        integer, dimension(nylcr), intent(in) :: modth
        integer, intent(in) :: nxlc,nylc,nzlc,nmod,&
                               nprocrow,nproccol,nylcr,nzlcr
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: pztable,xpzstable,pzdisp
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,comm2d,commcol,commrow
        integer :: nxpylc,nxpzlc,nsxy1,nsxy2,nsxz1,nsxz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr,ksign,kk
        double precision :: t0,pi,fourpi,testsum,totsum,testrho,totrho,rad,scaley
        double precision :: sec0,sec1,sec,tval,sumtest
        double precision :: sumrho0,sumrho0tot,sumrho1,sumrho1tot,sumrho2,&
        sumrho2tot,sumrho3,sumrho3tot,sumrho4,sumrho4tot,sumrho5,sumrho5tot
        double precision :: sumbess,sumbesstot,sumbn,sumbntot,sumod,sumodtot

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)

        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        pztable(nprocrow-1) = pztable(nprocrow-1)
        pzdisp(0) = 0
        do i = 1, nprocrow-1
          pzdisp(i) = pzdisp(i-1)+pztable(i-1)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo
        pytable(nproccol-1) = pytable(nproccol-1) - 1 !periodic BC.

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1)
        inyglb = glmshnm(2)-1  !periodic BC.
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        !used for transpose between serial x and parallel y.
        nsxy1 = inxglb/nproccol
        nsxy2 = inxglb - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo
        nxpylc = xpystable(myidy)

        !used for transpose between serial x and parallel z.
        nsxz1 = inxglb/nprocrow  !periodic BC
        nsxz2 = inxglb - nprocrow*nsxz1
        do i = 0, nprocrow-1
          if(i.le.(nsxz2-1)) then
            xpzstable(i) = nsxz1+1
          else
            xpzstable(i) = nsxz1
          endif
        enddo
        nxpzlc = xpzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        do k = 1, innz
          do j = 1, inny
            rho(1,j,k) = source(1,j+jadd,k+kadd)/(0.125*hx*hx*hy*hz)
            do i = 2, innx-1
              rho(i,j,k) = source(i,j+jadd,k+kadd)/((i-1)*hx*hx*hy*hz)*&
                           1.0
            enddo
            rho(innx,j,k) = source(innx,j+jadd,k+kadd)/&
            (0.25*(2*innx-2.5)*hx*hx*hy*hz)
          enddo
        enddo

        rad = (innx-1)*hx
        pi = 2*asin(1.0)
        !do k = 1, innz
        !  kk = k + pzdisp(myidz)
        !  do j = 1, inny
        !    do i = 1, innx
        !      if((i-1)*hx.le.0.5*rad) then
        !         !rho(i,j,k) = (i-1)*hx*(cos((j-1)*hy))**2*24.0/&
        !         !             (pi*rad**3)
        !         !rho(i,j,k) = 10.0
        !         rho(i,j,k) = exp(-10*((i-1)*hx+(kk-1)*hz)**2)
        !      else
        !         rho(i,j,k) = 0.0
        !      endif
        !    enddo
        !  enddo
        !enddo

!        sumrho0 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho0 = sumrho0 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho0,sumrho0tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho0: ",sumrho0tot,innx,inny,innz,myidz,myidy

        !sec0 = MPI_WTIME()
        ksign = 1
        scaley = 1.0
        call fft3d3_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        ksign,scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

!        sumrho1 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho1 = sumrho1 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho1,sumrho1tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho1: ",sumrho1tot

        call Bessr_Bessel(inxglb,inny,innz,rho,besscoef,&
                         bessnorm,hx,modth,myidy,nmod)

!        sumrho2 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho2 = sumrho2 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho2,sumrho2tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho2: ",sumrho2tot
!        sumbess = 0.0
!        sumbn = 0.0
!        if(myidy.eq.0) then
!
!        do k = 1,nmod
!          do j = 1, inxglb
!            sumbn = sumbn + bessnorm(j,k)
!            do i = 1, inxglb
!              sumbess = sumbess + besscoef(i,j,k)
!            enddo
!          enddo
!        enddo
!
!        else
!
!        do k = 1,nmod
!          do j = 1, inxglb
!            sumbn = sumbn + bessnorm(j,k)
!            do i = 1, inxglb
!              sumbess = sumbess + besscoef(i,j,k)
!            enddo
!          enddo
!        enddo
!
!        endif
!
!        sumod = 0.0
!        do i = 1, inny
!          sumod = sumod + modth(i)
!        enddo
!        call MPI_ALLREDUCE(sumbess,sumbesstot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(sumbn,sumbntot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(sumod,sumodtot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumbess: ",sumbesstot,sumbntot,sumodtot,nmod

        call Gaussz3(inxglb,inzglb,innz,nxpzlc,inny,rho,&
        xpzstable,pztable,nprocrow,commrow,myidz,myidy,gml,modth,hz,&
        nmod)

!        sumrho3 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho3 = sumrho3 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho3,sumrho3tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho3: ",sumrho3tot

        call BessrI_Bessel(inxglb,inny,innz,rho,besscoef,modth,&
                               myidy,nmod)

!        sumrho4 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho4 = sumrho4 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho4,sumrho4tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho4: ",sumrho4tot

        ksign = -1
        scaley = 1.0/float(inyglb)
        call fft3d3_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        ksign,scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)
        !sec1 = MPI_WTIME()
        !sec = sec1 - sec0
        !call MPI_ALLREDUCE(sec,tval,1,MPI_REAL,MPI_MAX, &
        !                     MPI_COMM_WORLD,ierr)
        !print*,"Field time: ",tval

!        sumrho5 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho5 = sumrho5 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho5,sumrho5tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho5: ",sumrho5tot

        !testrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      testrho = testrho + rho(i,j,k)
        !    enddo
        !  enddo
        !enddo
        !call MPI_ALLREDUCE(testrho,totrho,1,MPI_REAL,MPI_SUM, &
        !                     MPI_COMM_WORLD,ierr)
        !print*,"totrho: ",totrho

        fourpi = 4*pi
        ! 1/(4*pi) is due to the multiplication of 4pi in Rob's code.
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*fourpi
            enddo
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
!        if(myidy.eq.(nproccol-1)) then
!          do k = 1, innz
!            do i = 1, innx
!              this%FieldQ(i,inny+1+jadd,k+kadd)=rho(i,1,k)*fourpi
!            enddo
!          enddo
!        endif
!        sumtest = 0.0
!        do k = 1, innz
!          do j = 1, inny
!            do i = 1, innx
!              sumtest = sumtest + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        print*,"sumtest, potential: ",sumtest

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update3_FieldQuant

        subroutine Gaussz3(nx,nz,nsizez,nsizexz,nsizey,x,&
        xstable,xrtable,nprocrow,commrow,myidx,myidy,gm,modth,hz,nmod)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,nz,nsizez,nsizexz,nsizey,nmod
        integer,intent(in) :: nprocrow,commrow
        integer,intent(in) :: myidx,myidy
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        double precision, dimension(nx,nsizey,nsizez), intent(inout) :: x
        double precision, dimension(nx,nmod), intent(in) :: gm
        double precision, intent(in) :: hz
        integer, dimension(nsizey), intent(in) :: modth
        integer,dimension(0:nprocrow-1) :: rdisp
        integer :: i,j,k,jm
        double precision :: t0,tmp1,tmp2
        double precision, dimension(nz) :: aa,cc
        double precision, dimension(nz,nsizey,nsizexz) :: bb
        double precision, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        rdisp(0) = 0
        do i = 1, nprocrow-1
          rdisp(i) = rdisp(i-1) + xstable(i-1)
        enddo

        allocate(x0(nz,nsizey,nsizexz))
        !transpose z direction into serial direction:
        call trans3d3r_TRANSP(nx,nsizey,nsizez,nsizexz,x,x0,nprocrow,&
                      xstable,xrtable,commrow,myidx,nz)

        do i = 1, nz -1
          aa(i) = 1.0
        enddo
        aa(nz) = 0.0
        cc(1) = 0.0
        do i = 2, nz
          cc(i) = 1.0
        enddo
        do j = 1, nsizey
          if(myidy.ne.0) then
            jm = modth(j)-modth(1) + 1
          else
!            jm = modth(j)-modth(1) + 1
!            if(j.eq.2) jm = j
            if(j.le.2) then
              jm = j
            else
              jm = modth(j)-modth(1) + 2
            endif
          endif
          do k = 1, nsizexz
            tmp1 = hz*hz*gm(k+rdisp(myidx),jm)*gm(k+rdisp(myidx),jm)
            tmp2 = exp(-hz*gm(k+rdisp(myidx),jm))
            bb(1,j,k) = -2 - tmp1 + tmp2
            do i = 2, nz-1
              bb(i,j,k) = -2 - tmp1
            enddo
            bb(nz,j,k) = -2 - tmp1 + tmp2
          enddo
        enddo

        !eliminate the lower off-diagnal terms.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 2, nz
              bb(i,j,k) = bb(i,j,k) - aa(i-1)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo
        !the negative sign due to the "-" in the Poisson's equation.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 1, nz
              x0(i,j,k) = -x0(i,j,k)*hz*hz
            enddo
          enddo
        enddo

        !forward substitute for Gauss Elemination.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 2, nz
              x0(i,j,k) = x0(i,j,k) - x0(i-1,j,k)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute
        do k = 1, nsizexz
          do j = 1, nsizey
            x0(nz,j,k) = x0(nz,j,k)/bb(nz,j,k)
            do i = nz-1, 1, -1
              x0(i,j,k) = (x0(i,j,k) - x0(i+1,j,k)*aa(i))/bb(i,j,k)
            enddo
          enddo
        enddo

        call trans3d3r_TRANSP(nz,nsizey,nsizexz,nsizez,x0,x,nprocrow,&
                      xrtable,xstable,commrow,myidx,nx)

        deallocate(x0)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine Gaussz3

!-------------------------------------------------------------------------

        subroutine update4_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr,&
          gblam)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        type (FieldQuant), intent(inout) :: this
        double precision, intent(in) :: gblam
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable,xpzstable
        integer, dimension(0:nproccol-1) :: xpystable,pytable,zpystable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc,nypzlc,nxpzlc,nzpylc
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0,pi,fourpi,testsum,totsum,testrho,totrho

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)

        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        pztable(nprocrow-1) = pztable(nprocrow-1) - 1 !periodic BC.
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo
        pytable(nproccol-1) = pytable(nproccol-1) - 1 !periodic BC.

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1)
        inyglb = glmshnm(2)-1  !periodic BC.
        inzglb = glmshnm(3)-1  !periodic BC.

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        nsxy1 = inxglb/nproccol
        nsxy2 = inxglb - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo
        nxpylc = xpystable(myidy)

        ! /2 + 1 is from the double precision to complex fft.
        nsyz1 = (inyglb/2+1)/nprocrow  !periodic BC
        nsyz2 = (inyglb/2+1) - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo
        nypzlc = ypzstable(myidz)

        nsxy1 = inxglb/nprocrow
        nsxy2 = inxglb - nprocrow*nsxy1
        do i = 0, nprocrow-1
          if(i.le.(nsxy2-1)) then
            xpzstable(i) = nsxy1+1
          else
            xpzstable(i) = nsxy1
          endif
        enddo
        nxpzlc = xpzstable(myidz)

        nsyz1 = inzglb/nproccol  !periodic BC
        nsyz2 = inzglb - nproccol*nsyz1
        do i = 0, nproccol-1
          if(i.le.(nsyz2-1)) then
            zpystable(i) = nsyz1+1
          else
            zpystable(i) = nsyz1
          endif
        enddo
        nzpylc = zpystable(myidy)

        call getmsize_CompDom(fldgeom,msize)
        !hx = msize(1)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)
        !if(myid.eq.0) then
        !  print*,"hx: ",hx,hy,hz,(inxglb-1)*hx
        !endif

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        !from the "source(i,j,k): # of particles on the grid" to get the
        !density on the grid.
        do k = 1, innz
          do j = 1, inny
            rho(1,j,k) = source(1,j+jadd,k+kadd)/(0.25*hx*hx*hy*hz) !@- 0.25 -> 0.125
            do i = 2, innx-1
              rho(i,j,k) = source(i,j+jadd,k+kadd)/((i-1)*hx*hx*hy*hz)*1.0
            enddo
            rho(innx,j,k) = source(innx,j+jadd,k+kadd)/&
            (0.25*(2*innx-2.5)*hx*hx*hy*hz)
          enddo
        enddo

        ! Finite boundary conditions!
        call finiteBCnew(innx,inny,innz,rho,hx,hy,hz,nxpylc,nypzlc,&
        nxpzlc,nzpylc,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,xpzstable,zpystable,&
        inxglb,inyglb,inzglb,gblam)

        pi = 2*asin(1.0)
        fourpi = 4*pi
        ! 1/(4*pi) is due to the multiplication of 4pi in Rob's code.
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*fourpi
            enddo
          enddo
        enddo

        !wrong peiodic condition, works only if nprocrow = 1
!        if(myidz.eq.(nprocrow-1)) then
!          do j = 1, inny
!            do i = 1, innx
!              this%FieldQ(i,j+jadd,innz+1+kadd)=rho(i,j,1)*fourpi
!            enddo
!          enddo
!        endif
!        if(myidy.eq.(nproccol-1)) then
!         do k = 1, innz
!            do i = 1, innx
!              this%FieldQ(i,inny+1+jadd,k+kadd)=rho(i,1,k)*fourpi
!            enddo
!          enddo
!        endif
!        if((nprocrow.eq.1).and.(nproccol.eq.1)) then
!          do i = 1, innx
!            this%FieldQ(i,inny+1,innz+1) = rho(i,1,1)*fourpi
!          enddo
!        endif

        !open(8,file="phir",status="unknown")
        !do i = 1, innx
        !    write(8,1000)(i-1)*hx,this%FieldQ(i,1,1),&
        !    this%FieldQ(i,1,innz/2),this%FieldQ(i,1,innz)
        !enddo
        !close(8)
        !open(9,file="phithe",status="unknown")
        !do j = 1, inny+1
        !  write(9,1000)(j-1)*hy,this%FieldQ(2,j,2),&
        !  this%FieldQ(3,j,2),this%FieldQ(4,j,2)
        !enddo
        !close(9)
        !open(10,file="phiz",status="unknown")
        !do k = 1, innz+1
        !  write(10,1000)(k-1)*hz,this%FieldQ(1,1,k),&
        !  this%FieldQ(innx/2,1,k),this%FieldQ(innx,1,k)
        !enddo
        !close(10)

1000    format(4(1x,e13.6))

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update4_FieldQuant

        ! Solving Poisson's equation with open BCs.
        subroutine finiteBC(innx,inny,innz,rho,hx,hy,hz,nxpylc,nypzlc,&
        nxpzlc,nzpylc,myidz,myidy,npz,npy,commrow,commcol,comm2d,&
        pztable,pytable,ypzstable,xpystable,xpzstable,zpystable,inxglb,&
        inyglb,inzglb,gblam)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc,nypzlc,nxpzlc,nzpylc
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable,&
                                                   xpzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable,&
                                                   zpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz, gblam
        double precision :: scalex,scaley,scalez,pi
        integer :: i,j,k,kk,nzz
        double complex, dimension(innx,inny,innz) :: rhofft
        integer, dimension(0:npz-1) :: zdisp
        integer, dimension(0:npy-1) :: ydisp
        double precision, dimension(innx-1) :: aa,cc
        double precision, dimension(innx-1,inny,innz) :: bb
!        double precision, dimension(innx,inny,innz) :: testrho
!        double precision :: diffrho

        zdisp(0) = 0
        do i = 1, npz-1
          zdisp(i) = zdisp(i-1) + pztable(i-1)
        enddo
        ydisp(0) = 0
        do i = 1, npy-1
          ydisp(i) = ydisp(i-1) + pytable(i-1)
        enddo

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        !testrho = rho
        call fft3d4_FFT(inxglb,inyglb,inzglb,innz,inny,nxpylc,&
        nypzlc,1,scalex,scaley,scalez,rho,ypzstable,pztable,xpystable,&
        pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy,rhofft)

        aa(1) = 2.0/(0.5*hx*hx)
        cc(1) = 0.0
        do i = 2, innx-1
          aa(i) = 1.0/(hx*hx) + 0.5/((i-1)*hx*hx)
          cc(i) = 1.0/(hx*hx) - 0.5/((i-1)*hx*hx)
        enddo

        nzz = inzglb/2 + 1
        pi = 2*asin(1.0)
        do k = 1, innz
          kk = k + zdisp(myidz)
          do j = 1, inny
            !bb(1,j,k) = -2.0/(0.5*hx*hx) - &
            !            (2*pi*(k+zdisp(myidz)-1)/gblam)**2
            if(kk.le.nzz) then
              bb(1,j,k) = -2.0/(0.5*hx*hx) - &
                        (2*pi*(kk-1)/gblam)**2
            else
              bb(1,j,k) = -2.0/(0.5*hx*hx) - &
                        (2*pi*(kk-1-inzglb)/gblam)**2
            endif
            do i = 2, innx-1
              !bb(i,j,k) = -2.0/(hx*hx)-(j+ydisp(myidy)-1)**2/&
              !((i-1)*(i-1)*hx*hx)-(2*pi*(k+zdisp(myidz)-1)/gblam)**2
              if(kk.le.nzz) then
                bb(i,j,k) = -2.0/(hx*hx)-(j+ydisp(myidy)-1)**2/&
                ((i-1)*(i-1)*hx*hx)-(2*pi*(kk-1)&
                /gblam)**2
              else
                bb(i,j,k) = -2.0/(hx*hx)-(j+ydisp(myidy)-1)**2/&
                ((i-1)*(i-1)*hx*hx)-(2*pi*(kk-1-inzglb)&
                /gblam)**2
              endif
            enddo
          enddo
        enddo

        do k = 1, innz
          do j = 1, inny
            do i = 2, innx-1
              bb(i,j,k) = bb(i,j,k) - aa(i-1)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo
        !the negative sign due to the "-" in the Poisson's equation.
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx-2
              rhofft(i,j,k) = - rhofft(i,j,k)
            enddo
            rhofft(innx,j,k) = 0.0
            rhofft(innx-1,j,k) = -rhofft(innx-1,j,k)-aa(innx-1)*&
                                 rhofft(innx,j,k)
          enddo
        enddo

        !forward substitute for Gauss Elemination.
        do k = 1, innz
          do j = 1, inny
            rhofft(1,j,k) = rhofft(1,j,k)
            do i = 2, innx-1
              rhofft(i,j,k) = rhofft(i,j,k) &
                               -rhofft(i-1,j,k)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute.
        do k = 1, innz
          do j = 1, inny
            rhofft(innx-1,j,k) = rhofft(innx-1,j,k)/bb(innx-1,j,k)
            do i = innx-2, 1, -1
              rhofft(i,j,k) = (rhofft(i,j,k)-aa(i)*rhofft(i+1,j,k))/&
                               bb(i,j,k)
            enddo
          enddo
        enddo

        ! inverse FFT:
        scalex = 1.0/dble(inxglb)
        scaley = 1.0/dble(inyglb)
        scalez = 1.0/dble(inzglb)
        call invfft3d4_FFT(inxglb,inyglb,inzglb,innz,inny,nxpzlc,&
        nzpylc,-1,scalex,scaley,scalez,rhofft,xpzstable,pztable,&
        zpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)

        !diffrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      diffrho = diffrho + abs(testrho(i,j,k)-rho(i,j,k))
        !    enddo
        !  enddo
        !enddo
        !print*,"diffrho: ",diffrho

        end subroutine finiteBC

        ! Solving Poisson's equation with open BCs.
        subroutine finiteBCnew(innx,inny,innz,rho,hx,hy,hz,nxpylc,nypzlc,&
        nxpzlc,nzpylc,myidz,myidy,npz,npy,commrow,commcol,comm2d,&
        pztable,pytable,ypzstable,xpystable,xpzstable,zpystable,inxglb,&
        inyglb,inzglb,gblam)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc,nypzlc,nxpzlc,nzpylc
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable,&
                                                   xpzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable,&
                                                   zpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz, gblam
        double precision :: scalex,scaley,scalez,pi,hx2,fac2
        integer :: i,j,k,kk,nzz
        double complex, dimension(innx,inny,innz) :: rhofft
        integer, dimension(0:npz-1) :: zdisp
        integer, dimension(0:npy-1) :: ydisp
        double precision, dimension(innx-1) :: aa,cc
        double precision, dimension(innx-1,inny,innz) :: bb
!        double precision, dimension(innx,inny,innz) :: testrho
!        double precision :: diffrho

        zdisp(0) = 0
        do i = 1, npz-1
          zdisp(i) = zdisp(i-1) + pztable(i-1)
        enddo
        ydisp(0) = 0
        do i = 1, npy-1
          ydisp(i) = ydisp(i-1) + pytable(i-1)
        enddo

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        !testrho = rho
        call fft3d4_FFT(inxglb,inyglb,inzglb,innz,inny,nxpylc,&
        nypzlc,1,scalex,scaley,scalez,rho,ypzstable,pztable,xpystable,&
        pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy,rhofft)

        hx2 = hx*hx

        aa(1) = 2.0/(0.5*hx2)
        cc(1) = 0.0
        do i = 2, innx-1
          aa(i) = 1.0/(hx2) + 0.5/((i-1)*hx2)
          cc(i) = 1.0/(hx2) - 0.5/((i-1)*hx2)
        enddo

        nzz = inzglb/2 + 1
        pi = 2*asin(1.0)
        fac2 = (2.0*pi/gblam)**2
        do k = 1, innz
          kk = k + zdisp(myidz)
          do j = 1, inny
            if(kk.le.nzz) then
              bb(1,j,k) = -2.0/(0.5*hx2) - fac2*int(kk-1)**2
            else
              bb(1,j,k) = -2.0/(0.5*hx2) - fac2*int(kk-1-inzglb)**2
            endif
            do i = 2, innx-1
              if(kk.le.nzz) then
                bb(i,j,k) = -2.0/hx2-int(j+ydisp(myidy)-1)**2/&
                ((i-1)*(i-1)*hx2)-fac2*int(kk-1)**2
              else
                bb(i,j,k) = -2.0/hx2-int(j+ydisp(myidy)-1)**2/&
                ((i-1)*(i-1)*hx2)-fac2*int(kk-1-inzglb)**2
              endif
            enddo
          enddo
        enddo

        do k = 1, innz
          do j = 1, inny
            do i = 2, innx-1
              bb(i,j,k) = bb(i,j,k) - aa(i-1)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        bb = 1.0/bb
        !the negative sign due to the "-" in the Poisson's equation.
        rhofft(innx,:,:) = 0.0
        rhofft = -rhofft

        !forward substitute for Gauss Elemination.
        do k = 1, innz
          do j = 1, inny
            do i = 2, innx-1
              rhofft(i,j,k) = rhofft(i,j,k) &
                               -rhofft(i-1,j,k)*cc(i)*bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute.
        do k = 1, innz
          do j = 1, inny
            rhofft(innx-1,j,k) = rhofft(innx-1,j,k)*bb(innx-1,j,k)
            do i = innx-2, 1, -1
              rhofft(i,j,k) = (rhofft(i,j,k)-aa(i)*rhofft(i+1,j,k))*bb(i,j,k)
            enddo
          enddo
        enddo

        ! inverse FFT:
        scalex = 1.0/dble(inxglb)
        scaley = 1.0/dble(inyglb)
        scalez = 1.0/dble(inzglb)
        call invfft3d4_FFT(inxglb,inyglb,inzglb,innz,inny,nxpzlc,&
        nzpylc,-1,scalex,scaley,scalez,rhofft,xpzstable,pztable,&
        zpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)

        end subroutine finiteBCnew

!-------------------------------------------------------------------------
        subroutine update5_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        type (FieldQuant), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: pztable,xpzstable,zdisp
        integer, dimension(0:nproccol-1) :: xpystable,pytable,ydisp
        double precision, dimension(nxlc-1,nylcr,nzlcr) :: rho
        double precision, dimension(nxlc-1) :: tmprho
        integer :: myid,myidz,myidy,comm2d,commcol,commrow
        integer :: nxpylc,nxpzlc,nsxy1,nsxy2,nsxz1,nsxz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr,ksign,kk,jj
        double precision :: t0,pi,fourpi,testsum,totsum,testrho,totrho,&
                xrad,yrad,scaley
        !double precision , dimension(nxlc-1,nylcr,nzlcr) :: testt
        !double precision, dimension(nxlc) :: betan
        double precision :: length

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)

        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        zdisp(0) = 0
        do i = 0, nprocrow-1
          zdisp(i) = zdisp(i-1)+pztable(i-1)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo
        pytable(nproccol-1) = pytable(nproccol-1) - 1 !periodic BC.
        ydisp(0) = 0
        do i = 1, nproccol-1
          ydisp(i) = ydisp(i-1) + pytable(i-1)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1)-1  !periodic BC
        inyglb = glmshnm(2)-1  !periodic BC.
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc-1  !periodic BC along x.

        !used for transpose between serial x and parallel y.
        nsxy1 = inxglb/nproccol
        nsxy2 = inxglb - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo
        nxpylc = xpystable(myidy)

        !used for transpose between serial x and parallel z.
        nsxz1 = inxglb/nprocrow  !periodic BC
        nsxz2 = inxglb - nprocrow*nsxz1
        do i = 0, nprocrow-1
          if(i.le.(nsxz2-1)) then
            xpzstable(i) = nsxz1+1
          else
            xpzstable(i) = nsxz1
          endif
        enddo
        nxpzlc = xpzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)
        xrad = inxglb*hx
        yrad = inyglb*hy

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        do k = 1, innz
          do j = 1, inny
            rho(1,j,k) = source(1,j+jadd,k+kadd)/(hx*hy*hz)
            do i = 2, innx-1
              rho(i,j,k) = source(i,j+jadd,k+kadd)/(hx*hy*hz)
            enddo
            rho(innx,j,k) = source(innx,j+jadd,k+kadd)/&
            (hx*hy*hz)
          enddo
        enddo
        pi = 2*asin(1.0)
        !do k = 1, innz
        !  kk = k + zdisp(myidz)
        !  do j = 1, inny
        !    jj = j + ydisp(myidy)
        !    do i = 1, innx
        !      if((i.gt.innx/4).and.(i.lt.3*innx/4).and. &
        !         (jj.gt.inyglb/4).and.(jj.lt.3*inyglb/4) ) then
        !         !rho(i,j,k) = (i-1)*hx*(cos((j-1)*hy))**2*24.0/&
        !         !             (pi*rad**3)
        !         !rho(i,j,k) = 10.0
        !         rho(i,j,k) = exp(-10*((i-1)*hx+(kk-1)*hz)**2)
        !!      else
        !         rho(i,j,k) = 0.0
        !      endif
        !    enddo
        !  enddo
        !enddo
!
        !testt = rho

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              tmprho(i) = rho(i,j,k)
            enddo
            call sinft(tmprho,innx)
            do i = 1, innx
              rho(i,j,k) = tmprho(i)
            enddo
          enddo
        enddo

        scaley = 1.0
        call fft3d5_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

        call Gaussz5(inxglb,inzglb,innz,nxpzlc,inny,rho,&
        xpzstable,pztable,nprocrow,nproccol,commrow,myidz,myidy,ydisp,&
        hz,xrad,yrad)

        scaley = 2.0/float(inyglb)
        call fft3d5_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              tmprho(i) = rho(i,j,k)
            enddo
            call sinft(tmprho,innx)
            do i = 1, innx
              rho(i,j,k) = tmprho(i)*2/float(innx)
            enddo
          enddo
        enddo

        !testrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      !testrho = testrho + abs(rho(i,j,k)-testt(i,j,k))
        !      testrho = testrho + rho(i,j,k)
        !    enddo
        !  enddo
        !enddo
        !call MPI_ALLREDUCE(testrho,totrho,1,MPI_REAL,MPI_SUM, &
        !                     MPI_COMM_WORLD,ierr)
        !!print*,"totrho: ",totrho
        !print*,"testrho: ",testrho

        fourpi = 4*pi
        ! 1/(4*pi) is due to the multiplication of 4pi in Rob's code.
        do k = 1, innz
          do j = 1, inny
            this%FieldQ(1,j+jadd,k+kadd) = 0.0
            do i = 2, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*fourpi
            enddo
            this%FieldQ(innx+1,j+jadd,k+kadd) = 0.0
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
        if(myidy.eq.0) then
          do k = 1, innz
            do i = 1, innx+1
              this%FieldQ(i,1+jadd,k+kadd)= 0.0
            enddo
          enddo
        endif
        if(myidy.eq.(nproccol-1)) then
          do k = 1, innz
            do i = 1, innx+1
              this%FieldQ(i,inny+1+jadd,k+kadd)= 0.0
            enddo
          enddo
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update5_FieldQuant

        subroutine Gaussz5(nx,nz,nsizez,nsizexz,nsizey,x,&
        xstable,xrtable,nprocrow,nproccol,commrow,myidx,myidy,ydisp,&
        hz,xa,ya)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,nz,nsizez,nsizexz,nsizey
        integer,intent(in) :: nprocrow,commrow,nproccol
        integer,intent(in) :: myidx,myidy
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        double precision, dimension(nx,nsizey,nsizez), intent(inout) :: x
        double precision, intent(in) :: hz,xa,ya
        integer, dimension(0:nproccol-1), intent(in) :: ydisp
        integer,dimension(0:nprocrow-1) :: xdisp
        integer :: i,j,k,jm,kk,jj
        double precision :: t0,tmp1,tmp2,pi
        double precision, dimension(nz) :: aa,cc
        double precision, dimension(nz,nsizey,nsizexz) :: bb
        double precision, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        pi = 2.0*asin(1.0)
        xdisp(0) = 0
        do i = 1, nprocrow-1
          xdisp(i) = xdisp(i-1) + xstable(i-1)
        enddo

        allocate(x0(nz,nsizey,nsizexz))
        !transpose z direction into serial direction:
        call trans3d3r_TRANSP(nx,nsizey,nsizez,nsizexz,x,x0,nprocrow,&
                      xstable,xrtable,commrow,myidx,nz)

        do i = 1, nz -1
          aa(i) = 1.0
        enddo
        aa(nz) = 0.0
        cc(1) = 0.0
        do i = 2, nz
          cc(i) = 1.0
        enddo
        do k = 1, nsizexz
          kk = xdisp(myidx) + k
          do j = 1, nsizey
            jj = ydisp(myidy) + j
            tmp1 = hz*hz*(((jj-1)*pi/ya)**2+((kk-1)*pi/xa)**2)
            tmp2 = exp(-sqrt(tmp1))
            bb(1,j,k) = -2 - tmp1 + tmp2
            do i = 2, nz-1
              bb(i,j,k) = -2 - tmp1
            enddo
            bb(nz,j,k) = -2 - tmp1 + tmp2
          enddo
        enddo

        !eliminate the lower off-diagnal terms.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 2, nz
              bb(i,j,k) = bb(i,j,k) - aa(i-1)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo
        !the negative sign due to the "-" in the Poisson's equation.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 1, nz
              x0(i,j,k) = -x0(i,j,k)*hz*hz
            enddo
          enddo
        enddo

        if(myidx.ne.0 .or. myidy.ne.0) then
        !forward substitute for Gauss Elemination.
          do k = 1, nsizexz
          do j = 1, nsizey
            do i = 2, nz
              x0(i,j,k) = x0(i,j,k) - x0(i-1,j,k)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute
        do k = 1, nsizexz
          do j = 1, nsizey
            x0(nz,j,k) = x0(nz,j,k)/bb(nz,j,k)
            do i = nz-1, 1, -1
              x0(i,j,k) = (x0(i,j,k) - x0(i+1,j,k)*aa(i))/bb(i,j,k)
            enddo
          enddo
        enddo

        else

        !forward substitute for Gauss Elemination.
        do k = 2, nsizexz
          do j = 2, nsizey
            do i = 2, nz
              x0(i,j,k) = x0(i,j,k) - x0(i-1,j,k)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute
        do k = 2, nsizexz
          do j = 2, nsizey
            x0(nz,j,k) = x0(nz,j,k)/bb(nz,j,k)
            do i = nz-1, 1, -1
              x0(i,j,k) = (x0(i,j,k) - x0(i+1,j,k)*aa(i))/bb(i,j,k)
            enddo
          enddo
        enddo

        endif

        call trans3d3r_TRANSP(nz,nsizey,nsizexz,nsizez,x0,x,nprocrow,&
                      xrtable,xstable,commrow,myidx,nx)

        deallocate(x0)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return

        end subroutine Gaussz5
!-------------------------------------------------------------------------

        subroutine update6_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        type (FieldQuant), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: pztable,xpzstable,zdisp
        integer, dimension(0:nproccol-1) :: xpystable,pytable,ydisp
        double precision, dimension(nxlc-1,nylcr,nzlcr) :: rho
        double precision, dimension(nxlc-1) :: tmprho
        integer :: myid,myidz,myidy,comm2d,commcol,commrow
        integer :: nxpylc,nxpzlc,nsxy1,nsxy2,nsxz1,nsxz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr,ksign,jj,kk
        double precision :: t0,pi,fourpi,testsum,totsum,testrho,totrho,&
                xrad,yrad,scaley,zrad
        !double precision , dimension(nxlc-1,nylcr,nzlcr) :: testt
        !double precision, dimension(nxlc) :: betan
        double precision :: length

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)

        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        pztable(nprocrow-1) = pztable(nprocrow-1) - 1 !periodic BC.
        zdisp(0) = 0
        do i = 0, nprocrow-1
          zdisp(i) = zdisp(i-1)+pztable(i-1)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo
        pytable(nproccol-1) = pytable(nproccol-1) - 1 !periodic BC.
        ydisp(0) = 0
        do i = 1, nproccol-1
          ydisp(i) = ydisp(i-1) + pytable(i-1)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1)-1  !periodic BC
        inyglb = glmshnm(2)-1  !periodic BC.
        inzglb = glmshnm(3)-1  !periodic BC.

        innz = nzlcr
        inny = nylcr
        innx = nxlc-1  !periodic BC along x.

        !used for transpose between serial x and parallel y.
        nsxy1 = inxglb/nproccol
        nsxy2 = inxglb - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo
        nxpylc = xpystable(myidy)

        !used for transpose between serial x and parallel z.
        nsxz1 = inxglb/nprocrow  !periodic BC
        nsxz2 = inxglb - nprocrow*nsxz1
        do i = 0, nprocrow-1
          if(i.le.(nsxz2-1)) then
            xpzstable(i) = nsxz1+1
          else
            xpzstable(i) = nsxz1
          endif
        enddo
        nxpzlc = xpzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)
        xrad = inxglb*hx
        yrad = inyglb*hy
        zrad = inzglb*hz

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        do k = 1, innz
          do j = 1, inny
            rho(1,j,k) = source(1,j+jadd,k+kadd)/(hx*hy*hz)
            do i = 2, innx-1
              rho(i,j,k) = source(i,j+jadd,k+kadd)/(hx*hy*hz)
            enddo
            rho(innx,j,k) = source(innx,j+jadd,k+kadd)/&
            (hx*hy*hz)
          enddo
        enddo
        pi = 2*asin(1.0)
        !do k = 1, innz
        !  kk = k + zdisp(myidz)
        !  do j = 1, inny
        !    jj = j + ydisp(myidy)
        !    do i = 1, innx
        !      if((i.gt.innx/4).and.(i.lt.3*innx/4).and. &
        !         (jj.gt.inyglb/4).and.(jj.lt.3*inyglb/4) ) then
        !         !rho(i,j,k) = (i-1)*hx*(cos((j-1)*hy))**2*24.0/&
        !         !             (pi*rad**3)
        !         rho(i,j,k) = 10.0
        !         !rho(i,j,k) = exp(-10*((i-1)*hx+(k-1)*hz)**2)
        !      else
        !         rho(i,j,k) = 0.0
        !      endif
        !    enddo
        !  enddo
        !enddo
        !testrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      testrho = testrho + rho(i,j,k)*hx*hy*hz
        !    enddo
        !  enddo
        !enddo
        !print*,"testrho-first: ",testrho

        !testt = rho

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              tmprho(i) = rho(i,j,k)
            enddo
            call sinft(tmprho,innx)
            do i = 1, innx
              rho(i,j,k) = tmprho(i)
            enddo
          enddo
        enddo

        scaley = 1.0
        call fft3d6_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

        call Gaussz6(inxglb,inzglb,innz,nxpzlc,inny,rho,&
        xpzstable,pztable,nprocrow,nproccol,commrow,myidz,myidy,ydisp,&
        xrad,yrad,zrad)

        scaley = 2.0/float(inyglb)
        call fft3d6_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              tmprho(i) = rho(i,j,k)
            enddo
            call sinft(tmprho,innx)
            do i = 1, innx
              rho(i,j,k) = tmprho(i)*2/float(innx)
            enddo
          enddo
        enddo

        !testrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      !testrho = testrho + abs(rho(i,j,k)-testt(i,j,k))
        !      testrho = testrho + rho(i,j,k)
        !    enddo
        !  enddo
        !enddo
        !print*,"testrho: ",testrho
        !call MPI_ALLREDUCE(testrho,totrho,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
        !                     MPI_COMM_WORLD,ierr)
        !print*,"totrho: ",totrho,testrho

        fourpi = 4*pi
        ! 1/(4*pi) is due to the multiplication of 4pi in Rob's code.
        do k = 1, innz
          do j = 1, inny
            this%FieldQ(1,j+jadd,k+kadd) = 0.0
            do i = 2, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*fourpi
            enddo
            this%FieldQ(innx+1,j+jadd,k+kadd) = 0.0
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
        if(myidy.eq.0) then
          do k = 1, innz
            do i = 1, innx+1
              this%FieldQ(i,1+jadd,k+kadd)= 0.0
            enddo
          enddo
        endif
        if(myidy.eq.(nproccol-1)) then
          do k = 1, innz
            do i = 1, innx+1
              this%FieldQ(i,inny+1+jadd,k+kadd)= 0.0
            enddo
          enddo
        endif
        if(myidz.eq.(nprocrow-1)) then
          if(myidy.ne.(nproccol-1)) then
            do j = 1, inny
              do i = 1, innx+1
                this%FieldQ(i,j+jadd,innz+1+kadd)=rho(i,j,1)*fourpi
              enddo
            enddo
          else
            do j = 1, inny+1
              do i = 1, innx+1
                this%FieldQ(i,j+jadd,innz+1+kadd)=rho(i,j,1)*fourpi
              enddo
            enddo
          endif
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update6_FieldQuant

        subroutine Gaussz6(nx,nz,nsizez,nsizexz,nsizey,x,&
        xstable,xrtable,nprocrow,nproccol,commrow,myidx,myidy,ydisp,&
        xa,ya,za)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,nz,nsizez,nsizexz,nsizey
        integer,intent(in) :: nprocrow,commrow,nproccol
        integer,intent(in) :: myidx,myidy
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        double precision, dimension(nx,nsizey,nsizez), intent(inout) :: x
        double precision, intent(in) :: xa,ya,za
        integer, dimension(0:nproccol-1), intent(in) :: ydisp
        integer,dimension(0:nprocrow-1) :: xdisp
        integer :: i,j,k,jm,kk,jj,ii,ksign
        double precision :: t0,tmp2,pi,scalez
        double precision, dimension(nz,nsizey) :: tmp1,tmp10
        double precision, allocatable, dimension(:,:,:) :: x0
        integer :: jadd,kadd

        call starttime_Timer(t0)

        pi = 2.0*asin(1.0)
        xdisp(0) = 0
        do i = 1, nprocrow-1
          xdisp(i) = xdisp(i-1) + xstable(i-1)
        enddo

        allocate(x0(nz,nsizey,nsizexz))
        !transpose z direction into serial direction:
        call trans3d3r_TRANSP(nx,nsizey,nsizez,nsizexz,x,x0,nprocrow,&
                      xstable,xrtable,commrow,myidx,nz)

        if(myidx.ne.0) then
          kadd = 0
        else
          kadd = 1
        endif
        if(myidy.ne.0) then
          jadd = 0
        else
          jadd = 1
        endif

        ksign = 1
        scalez = 1.0
        pi = 2*asin(1.0)
        do k = 1+kadd, nsizexz
          kk = xdisp(myidx) + k - 1
          do j = 1, nsizey
            do i = 1, nz
              tmp1(i,j) = x0(i,j,k)
            enddo
          enddo

          call fftrclocal_FFT(ksign,scalez,tmp1,nz,nsizey,tmp10)

          do j = 1+jadd, nsizey
            jj = ydisp(myidy) + j - 1

            x0(1,j,k) = tmp10(1,j)/((kk*pi/xa)**2+(jj*pi/ya)**2)
            x0(2,j,k) = tmp10(2,j)/((kk*pi/xa)**2+(jj*pi/ya)**2+&
                          (2*(nz/2)*pi/za)**2)
            !x0(1,j,k) = tmp10(1,j)
            !x0(2,j,k) = tmp10(2,j)
            do i = 3, nz
              ii = (i+1)/2 - 1
              !no negative sign is due to -rho in Poisson's equation.
              x0(i,j,k) = tmp10(i,j)/((kk*pi/xa)**2+(jj*pi/ya)**2+&
                          (2*ii*pi/za)**2)
              !x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        ksign = -1
        scalez = 1.0/float(nz)
        !scalez = 2.0/float(nz)
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 1, nz
              tmp1(i,j) = x0(i,j,k)
            enddo
          enddo

          call fftcrlocal_FFT(ksign,scalez,tmp1,nz,nsizey,tmp10)

          do j = 1, nsizey
            do i = 1, nz
              x0(i,j,k) = tmp10(i,j) !using machine lib mfftcrlocal
              !x0(i,j,k) = tmp10(i,j)*scalez*2 !using NR lib
            enddo
          enddo
        enddo

        call trans3d3r_TRANSP(nz,nsizey,nsizexz,nsizez,x0,x,nprocrow,&
                      xrtable,xstable,commrow,myidx,nx)

        deallocate(x0)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine Gaussz6

!-------------------------------------------------------------------------
        subroutine setval_FieldQuant(this,i,j,k,value)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        integer, intent(in) :: i, j, k
        double precision, intent(in) :: value

        this%FieldQ(i,j,k) = value

        end subroutine setval_FieldQuant

        double precision function get_FieldQuant(this,i,j,k)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(in) :: i, j, k

        get_FieldQuant = this%FieldQ(i,j,k)

        end function get_FieldQuant

        subroutine getglb_FieldQuant(this,temp)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        type (FieldQuant), intent(out) :: temp
        integer :: i, j, k, lcnz,lcny,lcnx
        double precision :: value

        lcnx = this%Nxlocal
        lcny = this%Nylocal
        lcnz = this%Nzlocal

        do k = 1, lcnz
          do j = 1, lcny
            do i = 1, lcnx
              value = get_FieldQuant(this,i,j,k)
              call setval_FieldQuant(temp,i,j,k,value)
            enddo
          enddo
        enddo

        end subroutine getglb_FieldQuant

        subroutine getlcgrid_FieldQuant(this,nxlc,nylc,nzlc)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(out) :: nxlc,nylc,nzlc

        nxlc = this%Nxlocal
        nylc = this%Nylocal
        nzlc = this%Nzlocal

        end subroutine

        subroutine destruct_FieldQuant(this)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this

        deallocate(this%FieldQ)

        end subroutine destruct_FieldQuant

      end module FieldQuantclass
