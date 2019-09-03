module PyFunction
    use AdvAccelclass
    use Messaging
    implicit none
    include 'mpif.h'

    contains

    subroutine Setup()
        implicit none

        call init_Data()
        runseq = 1
    end subroutine Setup


    subroutine Teardown()
        implicit none
        double precision :: time

        call destruct_Accel(time)
        call mpi_Teardown(time)
    end subroutine Teardown


    subroutine input(filename)
        implicit none
        character, intent(in) :: filename*256

        call load_Input(filename)
        runseq = 2
    end subroutine


    subroutine load_data()
        implicit none

        call load_AllData()
    end subroutine load_data


    subroutine reset_data()
        implicit none

        call destruct_AdvData()
    end subroutine reset_data


    subroutine init_dist(usrbcurr, usrbmass, usrbkenergy, usrphsini)
        implicit none
        double precision, intent(in) :: usrbcurr, usrbmass, usrbkenergy, usrphsini

        call init_Particle(usrbcurr, usrbmass, usrbkenergy, usrphsini)
        runseq = 3
        dstseq = 1
    end subroutine init_dist


    subroutine run(ini,fin)
        implicit none
        integer, intent(in) :: ini, fin

        call run_Accel(ini, fin)
        runseq = 4
        dstseq = 1
    end subroutine run

    subroutine construct()
        implicit none

        call construct_Beamline(0)
        call construct_PhysConst(Bfreq)
    end subroutine

    subroutine configure0(readindex, rtlength, rtseg, rtstep, rtitype, rtvarr)
        implicit none
        integer, intent(in) :: readindex
        double precision, intent(out) :: rtlength
        integer, intent(out) :: rtseg, rtstep, rtitype
        double precision, dimension(27), intent(out) :: rtvarr

        call get_Beamln(readindex, rtlength, rtseg, rtstep, rtitype, rtvarr)
    end subroutine configure0


    subroutine configure1(readindex, readrow, newvalue, ierr)
        implicit none
        integer, intent(in) :: readindex, readrow
        double precision, intent(in) :: newvalue
        integer, intent(out) :: ierr

        if (readrow.eq.4) then
            ierr = 1
        else
            call configure(readindex, readrow, newvalue)
            ierr = 0
        endif
    end subroutine configure1


    subroutine search(readtype)
        implicit none
        integer, intent(in) :: readtype

        call search_Index(readtype)
    end subroutine search


    subroutine get_distdata(idin, sflag, force, ierr)
        implicit none
        integer, intent(in) :: idin, sflag, force
        integer, intent(out) :: ierr
        integer :: i, rid, tmpid

        if (sflag.eq.1 .or. idin.eq.0) then
            rid = idin
        else
            rid = -5
            tmpid = 0
            do i = 1, size(hbunchid)
                tmpid = tmpid + 1
                if (idin.eq.hbunchid(i)) then
                    rid = tmpid
                    exit
                endif
            enddo
        endif

        if (rid.le.-2) then
            ierr = 1
        else if ((rid.ge.-1) .and. (rid.ne.ninedata_id .or. force.eq.1)) then
            call get_Dist(rid)
            ierr = 0
        endif
    end subroutine get_distdata


    subroutine set_datadir(dname)
        implicit none
        character(256), intent(in) :: dname
        dataclass_dir = dname
    end subroutine set_datadir


    subroutine set_partdir(dname)
        implicit none
        character(256), intent(in) :: dname
        particle_dname = dname
    end subroutine set_partdir


    subroutine set_partfile(fname)
        implicit none
        character(256), intent(in) :: fname
        particle_fname = fname
    end subroutine set_partfile


    subroutine set_strpdir(dname)
        implicit none
        character(256), intent(in) :: dname
        strip_dir = dname
    end subroutine set_strpdir


    subroutine outhist(outname,qid)
        implicit none
        character, intent(in) :: outname*256
        integer, intent(in) :: qid

        call Output_Hist(outname,Flagdiag,qid)
    end subroutine outhist


    subroutine outdist(outname)
        implicit none
        character, intent(in) :: outname*256

        call Output_Dist(outname)
    end subroutine outdist


    subroutine export(outname)
        implicit none
        character, intent(in) :: outname*256

        call export_Setting(outname)
    end subroutine export


    subroutine get_iostat(ioinp, iodist, iodata, iostrp)
        implicit none
        integer, intent(out) :: ioinp, iodist, iodata, iostrp
        ioinp = ioinpfile
        iodist = iodistfile
        iodata = iodatafile
        iostrp = iostrpfile
    end subroutine get_iostat

end module