program main
use AdvAccelclass
implicit none
include 'mpif.h'
integer :: myrank, ierr
character :: filename*256, outname*256, keys*256, expname*256
integer :: i, run_mode, readA, readindex, readrow, readtype
double precision :: time,time0,time1,time2, time3, newvalue
integer, allocatable, dimension(:) :: srlst

double precision :: rtlength
integer :: rtseg, rtstep, rtitype, parasize
double precision, dimension(28) :: rtvarr

call mpi_Setup(time)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

Zposini = 0.0
idistnp = -1.0

run_mode = 0
filename = 'test.in'
do i=1, iargc()
    call getarg(i,keys)
    if (keys(1:1).eq.'-') then
        if (index(keys,'i').ne.0) run_mode = 1
        if (index(keys,'h').ne.0) run_mode = 2
        if (index(keys,'t').ne.0) run_mode = 3
    else
        filename = keys
    endif
enddo

if (run_mode.eq.1) then
    call load_Input(filename)
    if (ioinpfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
    
    if (myrank.eq.0) print*, "start load_AllData()"
    call load_AllData()
    if (iodatafile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)

    do while(readA.ne.9)
        readA = 0
        readindex = 0
        readrow = 0
        newvalue = 0.0
        outname = 'result'

        if (myrank.eq.0) then
            print*, '-------------------------------'
            print*, '--- IMPACT interactive mode ---'
            print*, '1 : Run simulation and output results'
            print*, '2 : Configure element parameter'
            print*, '3 : Reload input file'
            print*, '4 : Reload 3D grid data'
            print*, '5 : Index search by element type'
            print*, '6 : Check element parameter'
            print*, '7 : Export input-file'
            print*, '9 : Exit IMPACT'

            write (*,fmt='(a)', advance='no') '>>> '
            read(*,*), readA
        endif

        call MPI_BCAST(readA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if (readA.eq.-1) then
            call init_Particle()
            if (iodistfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
            call run_Accel()
            if (iostrpfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
            if (myrank.eq.0) print*, 'Finish simulation and results are stored in memory.'
        else if ((readA.eq.-2).or.(readA.eq.1)) then
            if (myrank.eq.0) then
                print*, 'Please input output-file header name. (default name is "result")'
                write (*,fmt='(a)', advance='no') '>>> '
                read(*,'(A)'), outname
                if (outname.eq.'') outname='result'
            endif
            call MPI_BCAST(outname,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

            if (readA.eq.1) then
                call init_Particle()
                if (iodistfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
                call run_Accel()
                if (iostrpfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
                if (myrank.eq.0) print*, 'Finish simulation.'
            endif

            call Output_Hist(outname,Flagdiag,1)
            call Output_Dist(outname)
            if (myrank.eq.0) print*, 'Finish results output to file.'
        else if (readA.eq.2) then
            if (myrank.eq.0) then
                print*, 'Pleas input element index, parameter row, and new value.'
                write (*,fmt='(a)', advance='no') '>>> '
                read(*,*) readindex, readrow, newvalue
            endif
            call MPI_BCAST(readindex,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(readrow,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(newvalue,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            if (readrow.le.4) then
               if (myrank.eq.0) print*, 'Underconstruction. Configuration is skipped.'
            else
                call configure(readindex, readrow-4, newvalue)
                if (myrank.eq.0) print*, 'Finish configuration.'
            endif
        else if (readA.eq.3) then
            if (myrank.eq.0) then
                print*, 'Please input input-file name. (default name is "test.in")'
                write (*,fmt='(a)', advance='no') '>>> '
                read(*,'(A)'), filename
                if (filename.eq.'') filename='test.in'
            endif
            call MPI_BCAST(filename,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
            call load_Input(filename)
            if (ioinpfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
        else if (readA.eq.4) then
            call destruct_AdvData()
            call load_AllData()
            if (iodatafile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
        else if (readA.eq.5) then
            if (myrank.eq.0) then
                print*, 'Please input type value for index seacrh.'
                write (*,fmt='(a)', advance='no') '>>> '
                read(*,*), readtype
                call search_Index(readtype)
                print*, 'Indices of element type :', readtype, ' are '
                print*, searchres
            endif
        else if (readA.eq.6) then
            if (myrank.eq.0) then
                print*, 'Please input index of element.'
                write (*,fmt='(a)', advance='no') '>>> '
                read(*,*), readindex
                call check_Index(readindex, ierr)
                if (ierr.eq.0) then
                    call get_Beamln(readindex, rtlength, rtseg, rtstep, rtitype, rtvarr)
                    parasize = 1
                    if (rtitype.lt.0) then
                        parasize = 7
                    else if (rtitype.eq.0) then
                        parasize = 1
                    else if (rtitype.eq.1) then
                        parasize = 8
                    else if (rtitype.eq.2) then
                        parasize = 4
                    else if (rtitype.eq.3) then
                        parasize = 8
                    else if (rtitype.eq.4) then
                        parasize = 10
                    else if (rtitype.eq.5) then
                        parasize = 13
                    else if (rtitype.eq.101) then
                        parasize = 24
                    else if (rtitype.eq.102) then
                        parasize = 10
                    else if (rtitype.eq.103) then
                        parasize = 10
                    else if (rtitype.eq.104) then
                        parasize = 10
                    else if (rtitype.eq.105) then
                        parasize = 11
                    else if (rtitype.eq.106) then
                        parasize = 11
                    else if (rtitype.eq.110) then
                        parasize = 19
                    else
                    endif

                    print*, 'Parameters of element id :', readindex, ' are '
                    print*, rtlength, rtseg, rtstep, rtitype, rtvarr(1:parasize)
                endif
            endif
        else if (readA.eq.7) then
            if (myrank.eq.0) then
                print*, 'Please input export-file name. (default name is "test.in.new")'
                write (*,fmt='(a)', advance='no') '>>> '
                read(*,'(A)'), expname
                if (expname.eq.'') expname='test.in.new'
                call export_Setting(expname)
                print*, 'Finish exporting.'
            endif
        endif
    enddo
elseif (run_mode.eq.2) then
    print*, "usage : AdvImpact [option] [file]"
    print*, "Options and arguments :"
    print*, "-h    : print this message without running simulation"
    print*, "-i    : interactive mode"
    print*, "file  : user-specified input file name (default is 'test.in')"

elseif (run_mode.eq.3) then

    if (myrank.eq.0) time0 = MPI_WTIME()

    call load_Input(filename)
    if (ioinpfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)

    if (myrank.eq.0) print*, "start load_AllData()"
    call load_AllData()
    if (iodatafile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)

    if (myrank.eq.0) print*, "start run_Accel()"
    call init_Particle()
    if (iodistfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
    
    call run_Accel()
    if (iostrpfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
    
    filename = 'result1'
    call Output_Hist(filename,Flagdiag,1)
    call Output_Dist(filename)

    if (myrank.eq.0) time1 = MPI_WTIME()

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (myrank.eq.0) print*, "start 2nd run"

    call init_Particle()
    if (iodistfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
    
    if (myrank.eq.0) time3 = MPI_WTIME()
    
    call run_Accel()
    if (iostrpfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)

    if (myrank.eq.0) time2 = MPI_WTIME()

    filename = 'result2'
    call Output_Hist(filename,Flagdiag,1)
    call Output_Dist(filename)

    if (myrank.eq.0) then
        print*,"--- 1st run time ------"
        print*,time1-time0
        print*,"--- 2nd run time ------"
        print*,time2-time1
    endif

else

    if (myrank.eq.0) time0 = MPI_WTIME()

    call load_Input(filename)
    if (ioinpfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)

    call load_AllData()
    if (iodatafile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)

    call init_Particle()
    if (iodistfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)
    
    call run_Accel()
    if (iostrpfile.ne.0) call MPI_ABORT(MPI_COMM_WORLD, 1)

    if (myrank.eq.0) time1 = MPI_WTIME()

    filename = 'fort'
    call Output_Hist(filename,Flagdiag,1)
    call Output_Dist(filename)

endif


!if (myrank.eq.0) print*, "start destruct_Accel()"
call destruct_Accel(time)
call mpi_Teardown(time)

end program main
