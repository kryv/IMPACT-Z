module Messaging

    !-----------------!
    ! input variables !
    !-----------------!

    integer :: Npcol, Nprow, Dm, Npin, Nx, Ny, Nz

    integer :: Flagdist, Flagerr, Flagsubstep, Flagbc, Flagmap, Flagdiag, Rstartflg
    double precision :: Bcurr, Bkenergy, Bmass, Bcharge, Bfreq, Perdlen, &
                        Xwallrad, Ywallrad, Phsini
    double precision :: Zposini = 0.0
    double precision, dimension(21) :: Distparam

    integer :: Nchrgin, Nchrg
    double precision :: Bchargein
    double precision, dimension(100) :: Currlistin, Currlist, Qmcclistin, Qmcclist
    integer, dimension(100) :: Nptlistin, Nptlist

    double precision :: Outputflag

    integer, allocatable, dimension(:) :: Blnseg, Blnstp, Blntyp
    double precision, allocatable, dimension(:) :: Blnlen
    double precision, allocatable, dimension(:,:) :: Blnparams

    !-----------------!
    ! inner variables !
    !-----------------!

    character :: version_fsrc*256 = '1.0.3'
    
    integer :: Np, Nplocal, Nxlocal, Nylocal, Nzlocal, Nblem

    integer :: iqr, idr, ibpm, iccl, iccdtl, idtl, isc, icf, islrf,&
               isl, isl2, idipole, iemfld, imultpole, irfq, Idiag, Idiagbpm, Idstout
    integer, allocatable, dimension(:) :: Ielems, Blnidg
    double precision, allocatable, dimension(:) :: Blnpos

    integer :: Lcini = 1
    integer :: Lcfin = 1
    integer :: runseq = 0
    integer :: dstseq = 0
    integer :: flagmcbin = 0

    integer :: isrd = 1010      ! Seed for random number generators

    integer, allocatable, dimension(:) :: searchres
    integer, allocatable, dimension(:,:) :: hchrgpos

    double precision, allocatable, dimension(:) :: qmlabel

    double precision, allocatable, dimension(:) :: &
        hrefz, hrefphi, hrefgma, hrefeng, hrefbeta
    double precision, allocatable, dimension(:,:) :: &
        hxcen, hycen, hzcen, hxrms, hyrms, hzrms, &
        hxpcen, hypcen, hzpcen, hxprms, hyprms, hzprms, &
        hxtwsa, hytwsa, hztwsa, hxtwsb, hytwsb, hztwsb, &
        hxepsn, hyepsn, hzepsn, hxmax, hymax, hphimax, &
        hxpmax, hypmax, hdemax, hrefr, &
        hx3rd, hy3rd, hz3rd, hxp3rd, hyp3rd, hzp3rd, &
        hx4th, hy4th, hz4th, hxp4th, hyp4th, hzp4th, &
        hxepsnp, hyepsnp, hzepsnp, hxepsnf, hyepsnf, hzepsnf, &
        hrrms, hrnpr, hrmax

    integer, allocatable, dimension(:) :: hlpmin, hlpmax, hnchg
    integer, allocatable, dimension(:,:) :: hnpchg, hptot

    integer, allocatable, dimension(:) :: hbunchid, hbunchsmp

    double precision, allocatable, dimension(:,:) :: ninedata, idistd
    double precision :: ninedata_bg, ninedata_mass, ninedata_scxl
    integer :: ninedata_id, idistnp

end module Messaging