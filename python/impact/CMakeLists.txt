set (f2py_include_dir
    ${CMAKE_BINARY_DIR}/modules
    ${MPI_Fortran_INCLUDE_PATH}
)

set (f2py_library_dir
    ${CMAKE_BINARY_DIR}/lib
)

set (f2py_library
    advimpact-ctrl
    impact-core
)

f2py_add_package(_impactf
    PySubroutine.f90
    ../../AdvContrl/Messaging.f90
    ../../DataStruct/PhysConst.f90
    ${MPI_Fortran_LIBRARIES}
)

install(
FILES ${CMAKE_BINARY_DIR}/python/impact/_impactf${Python_EXT_SUFFIX}
DESTINATION ${Python_MODULE_DIR}/impact
)

set (py_src
    __init__.py
    control.py
    input.py
    result.py
    const.py
    util.py
)

set (py_test
    test/__init__.py
    test/test_beam_dynamics.py
    test/lattice_sample_solenoid
    test/lattice_sample_quadrupole
    test/lattice_sample_multipole
    test/lattice_sample_hsol
    test/lattice_sample_equad
    test/lattice_sample_dipole
    test/lattice_sample_edipole
    test/lattice_sample_ccl
    test/lattice_sample_rfq
    test/lattice_sample_errors
    test/rfdata222
    test/rfqline111
)

foreach(pyfile IN LISTS py_src py_test)
    get_filename_component(_dir ${pyfile} DIRECTORY)
    add_custom_command(
    TARGET _impactf POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/${pyfile}
            ${CMAKE_CURRENT_BINARY_DIR}/${pyfile}
    )
    install(FILES ${pyfile}
    DESTINATION ${Python_MODULE_DIR}/impact/${_dir}
    )
endforeach()