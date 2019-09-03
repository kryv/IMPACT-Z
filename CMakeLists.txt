cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)

project (IMPACT)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

option(BLD_ORIGINAL "Build original version IMPACT" OFF)
option(BLD_PYAPI "Build python interface" ON)

message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message (STATUS "Build original version IMPACT: ${BLD_ORIGINAL}")
message (STATUS "Build python interface: ${BLD_PYAPI}")

set (CMAKE_Fortran_FLAGS "-fPIC -pipe")
set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
set (CMAKE_Fortran_FLAGS_DEBUG "-O1 -g -pg")

find_package (MPI REQUIRED)
set (CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
include_directories (${MPI_Fortran_INCLUDE_PATH})

if(BLD_ORIGINAL)
    set (AppNames
        Impact
        AdvImpact
    )
else()
    set (AppNames
        AdvImpact
    )
endif()

set (main_Impact Contrl/main.f90)

set (Control_Impact
    Contrl/Input.f90
    Contrl/Output.f90
    Contrl/AccSimulator.f90
)

set (main_AdvImpact AdvContrl/AdvMain.f90)

set (Control_AdvImpact
    AdvContrl/AdvData.f90
    AdvContrl/AdvInput.f90
    AdvContrl/AdvOutput.f90
    AdvContrl/AdvAccel.f90
    AdvContrl/Messaging.f90
)

set (DataStruct
    DataStruct/NumConst.f90
    DataStruct/PhysConst.f90
    DataStruct/Pgrid.f90
    DataStruct/Data.f90
)

set (Func
    Func/Timer.f90
    Func/Transpose.f90
    Func/Fldmger.f90
    Func/Ptclmger.f90
    Func/FFT.f90
    Func/Bessel.f90
    Func/MTrnd.f90
)

set (Appl
    Appl/BPM.f90
    Appl/CCL.f90
    Appl/CCDTL.f90
    Appl/DTL.f90
    Appl/SC.f90
    Appl/DriftTube.f90
    Appl/Quadrupole.f90
    Appl/ConstFoc.f90
    Appl/SolRF.f90
    Appl/Sol.f90
    Appl/Sol2.f90
    Appl/Dipole.f90
    Appl/EMfld.f90
    Appl/Multipole.f90
    Appl/RFQ.f90
    Appl/BeamLineElem.f90
    Appl/CompDom.f90
    Appl/BeamBunch.f90
    Appl/Field.f90
    Appl/Distribution.f90
)

add_library(impact-core STATIC
        ${DataStruct}
        ${Func}
        ${Appl}
)
target_link_libraries(impact-core ${MPI_Fortran_LIBRARIES})

set_target_properties (impact-core PROPERTIES
    Fortran_MODULE_DIRECTORY modules
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

foreach (appname ${AppNames})

    string(TOLOWER ${appname} libname)

    add_library(${libname}-ctrl STATIC
        ${Control_${appname}}
    )

    target_link_libraries(${libname}-ctrl
        impact-core
        ${MPI_Fortran_LIBRARIES}
    )

    set_target_properties (${libname}-ctrl PROPERTIES
        Fortran_MODULE_DIRECTORY modules
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    )

    add_executable (${appname}
        ${main_${appname}}
    )

    target_link_libraries(${appname}
        impact-core
        ${libname}-ctrl
        ${MPI_Fortran_LIBRARIES}
    )

    set_target_properties (${appname} PROPERTIES
        Fortran_MODULE_DIRECTORY modules
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

    install(TARGETS ${appname}
        RUNTIME DESTINATION bin
    )
endforeach ()

if (BLD_PYAPI)
    find_package(PythonInterp 2.6)

    if(PYTHONINTERP_FOUND)
        execute_process(COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/pyconfig.py
            ${CMAKE_CURRENT_BINARY_DIR}/pyconfig.cmake
            RESULT_VARIABLE _PY_CONFIG
        )
        if(_PY_CONFIG)
            message(STATUS "Failed to fetch build information for ${PYTHON_EXECUTABLE}")
        else()
            include(${CMAKE_CURRENT_BINARY_DIR}/pyconfig.cmake)
            message(STATUS "Python version ${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}.${Python_VERSION_PATCH}")
            message(STATUS "  modules ${Python_MODULE_DIR}")
        endif()

        if (Python_NUMPY_FOUND)
            message(STATUS "  | numpy version ${Python_NUMPY_VERSION}")
        else()
            message(SEND_ERROR " Python module: 'numpy' not found")
        endif()

        if (Python_MPI4PY_FOUND)
            message(STATUS "  | mpi4py version ${Python_MPI4PY_VERSION}")
        else()
            message(SEND_ERROR " Python module: 'mpi4py' not found")
        endif()

        if (Python_NOSE_FOUND)
            message(STATUS "  | nose version ${Python_NOSE_VERSION}")
            enable_testing()
        else()
            message(SEND_ERROR " Python module: 'nose' not found")
        endif()

        include(f2py.cmake)
        add_subdirectory(python/impact)

        add_test(
            NAME IMPACT_Python_API_beamdynamics_test
            COMMAND ${PYTHON_EXECUTABLE} -m nose.core --exe python/impact
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
    else()
        message(STATUS "Python not found")
    endif()
endif()

SET(CPACK_GENERATOR "DEB")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kei Fukushima")
SET(CPACK_PACKAGE_VERSION "0.4.0")

INCLUDE(CPack)
