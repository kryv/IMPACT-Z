find_program(F2PY_EXECUTABLE NAMES
            "f2py${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}"
            "f2py-${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}"
            "f2py${PYTHON_VERSION_MAJOR}"
            "f2py"
            REQUIRED)

if(NOT F2PY_EXECUTABLE)
    message(FATAL_ERROR "Error: f2py not found.")
else()
    message(STATUS "Found f2py: ${F2PY_EXECUTABLE}")
endif()

macro(f2py_add_package _name)

    string(REGEX REPLACE ";?DESTINATION.*" "" _srcs "${ARGN}")
    string(REGEX MATCH "DESTINATION;.*" _dest_dir "${ARGN}")
    string(REGEX REPLACE "^DESTINATION;" "" _dest_dir "${_dest_dir}")

    set(_abs_srcs)
    foreach(_src ${_srcs})
        get_filename_component(_abs_src ${_src} ABSOLUTE)
        list(APPEND _abs_srcs ${_abs_src})
    endforeach(_src ${_srcs})
   
    set(_incs)
    foreach(_dir ${f2py_include_dir})
        list(APPEND _incs "-I${_dir}")
    endforeach(_dir)

    set(_libs)
    foreach(_dir ${f2py_library_dir})
        list(APPEND _libs "-L${_dir}")
    endforeach(_dir)
    
    set(_flibs)
    foreach(_lib ${f2py_library})
        list(APPEND _flibs "-l${_lib}")
    endforeach(_lib)
    
    add_custom_command(
    OUTPUT ${_name}.so
    COMMAND ${F2PY_EXECUTABLE} --quiet --f90exec=${CMAKE_Fortran_COMPILER}
                               --opt=-Wno-unused
                               -c ${_abs_srcs} -m ${_name}
                               ${_incs} ${_libs} ${_flibs}
    DEPENDS ${_srcs} ${f2py_library}
    )

    add_custom_target(${_name} ALL DEPENDS ${_name}.so)
endmacro(f2py_add_package)
