#!/usr/bin/env python
"""
Emit a file suitible to include() in a CMakeLists.txt file
with information from a python interpreter
Compatible for python 2.6 -> 3.4
"""

from __future__ import print_function

import os
import sys

if len(sys.argv)<2:
    out = sys.stdout
else:
    out = open(sys.argv[1], 'w')

from distutils.sysconfig import get_config_var, get_python_inc, get_python_lib

incdirs = [get_python_inc()]
libdirs = [get_config_var('LIBDIR')]
extsuff = get_config_var('EXT_SUFFIX') or get_config_var('SO')

# --------------------- #
# Python module checker #
# --------------------- #

# Check NumPy
have_numpy=['NO', None]
try:
    import numpy
    from numpy.distutils.misc_util import get_numpy_include_dirs
    incdirs += get_numpy_include_dirs()
    have_numpy[0] = 'YES'
    have_numpy[1] = numpy.__version__
except ImportError:
    pass

# Check mpi4py
have_mpi4py=['NO', None]
try:
    import mpi4py
    from mpi4py import MPI
    have_mpi4py[0] = 'YES'
    have_mpi4py[1] = mpi4py.__version__
except ImportError:
    pass

# Check nose
have_nose=['NO', None]
try:
    import nose
    have_nose[0] = 'YES'
    have_nose[1] = nose.__version__
except ImportError:
    pass


libdirs = [get_config_var('LIBDIR')]

# prepend introspected numpy directory so that it is checked before
# system python directory, which may contained a different version
# when virtualenv is used.  Debian helpfully symlinks the numpy headers
# as /usr/include/pythonX.Y/numpy :P
libdirs.reverse()
incdirs.reverse()

# location of extension modules relative to prefix (eg. "lib/python3/dist-packages")
moddir = os.path.relpath(get_python_lib(), get_config_var('exec_prefix'))

print('set(Python_DEFINITIONS, "%s")'%get_config_var('BASECFLAGS'), file=out)

print('set(Python_VERSION "%s")'%get_config_var('VERSION'), file=out)
print('set(Python_VERSION_LD "%s")'%(get_config_var('LDVERSION') or get_config_var('VERSION')), file=out)
print('set(Python_INCLUDE_DIRS "%s")'%';'.join(incdirs), file=out)
print('set(Python_LIBRARY_DIRS "%s")'%';'.join(libdirs), file=out)
print('set(Python_MODULE_DIR "%s")'%moddir, file=out)

print('set(Python_NUMPY_FOUND %s)'%have_numpy[0], file=out)
print('set(Python_NUMPY_VERSION %s)'%have_numpy[1], file=out)
print('set(Python_MPI4PY_FOUND %s)'%have_mpi4py[0], file=out)
print('set(Python_MPI4PY_VERSION %s)'%have_mpi4py[1], file=out)
print('set(Python_NOSE_FOUND %s)'%have_nose[0], file=out)
print('set(Python_NOSE_VERSION %s)'%have_nose[1], file=out)

print('set(Python_VERSION_MAJOR %s)'%sys.version_info[0], file=out)
print('set(Python_VERSION_MINOR %s)'%sys.version_info[1], file=out)
print('set(Python_VERSION_PATCH %s)'%sys.version_info[2], file=out)

print('set(Python_EXT_SUFFIX %s)'%extsuff , file=out)
print('set(Python_FOUND YES)', file=out)