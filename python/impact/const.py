from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from collections import OrderedDict
import numpy as np

__authors__ = "Kei Fukushima"
__copyright__ = "(c) 2017, Facility for Rare Isotope beams, Michigan State University"
__contact__ = "Kei Fukushima <fukushim@frib.msu.edu>"

clight = 299792458.0

ulength = {'m': 1.0, 'mm': 1e3, 'inch': 0.0254}
uradian = {'rad': 1.0, 'mrad': 1e3, 'deg': 180.0/np.pi}
udegree = {'deg': 1.0, 'rad': np.pi/180.0, 'mrad': (np.pi/180.0)*1e3}
uemit = {'m-rad': 1.0, 'mm-mrad': 1e6}
uenergy = {'eV': 1.0, 'MeV': 1e-6}
uMenergy = {'eV': 1e6, 'MeV': 1.0}
uNone = {'1': 1.0}

uelem = {'xmin': '[m]', 'xmax': '[m]', 'ymin': '[m]', 'ymax': '[m]',
         'xshift': '[m]', 'xpshift': '[rad]',
         'yshift': '[m]', 'ypshift': '[rad]',
         'zshift': '[deg]', 'zpshift': '[MeV]', 'phase shift': '[deg]',
         'aperture': '[m]', 'xaperture': '[m]', 'yaperture': '[m]',
         'dx': '[m]', 'dy': '[m]',
         'pitch': '[rad]', 'yaw': '[rad]', 'roll': '[rad]',
         'bz0': '[T]', 'gradient': '[T/m]',
         'driven phase': '[deg]', 'input phase': '[deg]',
         'frequency': '[Hz]', 'frequency2': '[Hz]', 'frequency3': '[Hz]',
         'scaling': '[1]', 'scaling2': '[1]', 'scaling3': '[1]',
         'frequency ratio 2': '[1]', 'frequency ratio 3': '[1]',
         'driven phase offset 2': '[deg]', 'driven phase offset 3': '[deg]',
         'scale error': '[1]', 'scale error 2': '[1]', 'scale error 3': '[1]',
         'phase error': '[deg]', 'phase error 2': '[deg]', 'phase error 3': '[deg]',
         'synchronous phase offset': '[deg]'
         }

ddist = {'Rectangle': 1, 'rectangle': 1, 'rc': 1,
         'Gaussian': 2, 'gaussian': 2, 'gs': 2,
         'Waterbag': 3, 'waterbag': 3, 'wb': 3,
         'Semi-Gaussian': 4, 'semi-gaussian': 4, 'sg': 4,
         'KV': 5, 'kv': 5,
         'mc-Waterbag': 16, 'mc-waterbag': 16, 'mcwb': 16, 'wbmc': 16,
         'mc-Gaussian': 17, 'mc-gaussian': 17, 'mcgs': 17, 'gsmc': 17,
         'User-Input': 19, 'user-input': 19, 'ui': 19
         }

elemtype = OrderedDict((
    ('freq jump', -27),
    ('monitor', -23),
    ('corrector', -21), ('cor', -21),
    ('collimator slit', -13), ('col', -13), ('slit', -13),
    ('drift', 0),
    ('quadrupole', 1), ('quad', 1),
    ('constant focusing', 2), ('constfocus', 2),
    ('solenoid', 3), ('sol', 3),
    ('solenoid2', 13), ('sol2', 13), ('hsolenoid', 13), ('hsol', 13),
    ('dipole', 4), ('bend', 4),
    ('multipole', 5), ('mpole', 5),
    ('drift tube linac', 101), ('dtl', 101),
    ('coupled cavity drift tube linac', 102), ('ccdtl', 102),
    ('coupled cavity linac', 103), ('ccl', 103),
    ('superconducting cavity', 104), ('scc', 104),
    ('solenoid-rf', 105), ('solrf', 105),
    ('rfquadrupole', 106), ('rfquad', 106), ('rfq', 106),
    ('emfield', 110), ('emfld', 110),
    ('em field dipole', 114), ('emfdipole', 114), ('emfbend', 114)
))

dopt = OrderedDict((
    ('option1', 5),
    ('option2', 6),
    ('option3', 7),
    ('option4', 8),
    ('option5', 9),
    ('option6', 10),
    ('option7', 11),
    ('option8', 12),
    ('option9', 13),
    ('option10', 14)
))

dcor = OrderedDict((
    ('aperture', 5), ('aper', 5), ('pipe', 5), ('radius', 5),
    ('xshift', 6),
    ('xpshift', 7), ('pxshift', 7), ('theta_x', 7),
    ('yshift', 8),
    ('ypshift', 9), ('pyshift', 9), ('theta_y', 9),
    ('zshift', 10),
    ('zpshift', 11), ('pzshift', 11)
))

dslit = OrderedDict((
    ('xmin', 5),
    ('xmax', 6),
    ('ymin', 7),
    ('ymax', 8)
))

dbps = OrderedDict((
    ('parameter id', 5), ('pid', 5), ('id', 5),
    ('min', 6),
    ('max', 7),
    ('reference energy offset', 8), ('offset', 8)
))

dfjump = OrderedDict((
    ('dummy', 5),
    ('frequency', 6),
    ('phase shift', 7), ('zshift', 7), ('shift', 7)
))

ddrift = OrderedDict((
    ('aperture', 5), ('aper', 5), ('pipe', 5), ('radius', 5)
))

dquad = OrderedDict((
    ('gradient', 5), ('grad', 5), ('b2', 5), ('k', 5), ('voltage', 5),
    ('switch', 6), ('flag', 6),
    ('aperture', 7), ('aper', 7), ('pipe', 7), ('radius', 7),
    ('dx', 8), ('dx0', 8),
    ('dy', 9), ('dy0', 9),
    ('pitch', 10), ('dx1', 10),
    ('yaw', 11), ('dy1', 11),
    ('roll', 12)
))

dcfcs = OrderedDict((
    ('kx0^2', 5), ('kx2', 5),
    ('ky0^2', 6), ('ky2', 6),
    ('kz0^2', 7), ('kz2', 7),
    ('aperture', 8), ('aper', 8), ('pipe', 8), ('radius', 8),
))

dsol = OrderedDict((
    ('bz0', 5), ('bz', 5), ('b', 5),
    ('fileid', 6), ('file', 6), ('id', 6),
    ('aperture', 7), ('aper', 7), ('pipe', 7), ('radius', 7),
    ('dx', 8), ('dx0', 8),
    ('dy', 9), ('dy0', 9),
    ('pitch', 10), ('dx1', 10),
    ('yaw', 11), ('dy1', 11),
    ('roll', 12)
))

dsol2 = OrderedDict((
    ('bz0', 5), ('bz', 5), ('b', 5),
    ('fileid', 6), ('file', 6), ('id', 6),
    ('aperture', 7), ('aper', 7), ('pipe', 7), ('radius', 7),
    ('by', 8),
    ('bx', 9),
    ('dx', 10), ('dx0', 10),
    ('dy', 11), ('dy0', 11),
    ('pitch', 12), ('dx1', 12),
    ('yaw', 13), ('dy1', 13),
    ('roll', 14)
))

ddip = OrderedDict((
    ('angle', 5), ('phi', 5),
    ('beta*gamma', 6), ('bg', 6),
    ('switch', 7), ('flag', 7),
    ('aperture', 8), ('aper', 8), ('pipe', 8), ('radius', 8),
    ('front anlge', 9), ('angle1', 9), ('phi1', 9),
    ('back angle', 10), ('angle2', 10), ('phi2', 10),
    ('front curvature', 11), ('curv1', 11),
    ('back curvature', 12), ('curv2', 12),
    ('fringe q-term', 13), ('kf', 13),
    ('main q-term', 14), ('k', 14)
))

dmlt = OrderedDict((
    ('gradient', 5), ('grad', 5), ('b2', 5),
    ('b3', 6),
    ('b4', 7),
    ('b5', 8),
    ('b6', 9),
    ('b7', 10),
    ('b8', 11),
    ('aperture', 12), ('aper', 12), ('pipe', 12), ('radius', 12),
    ('dx', 13), ('dx0', 13),
    ('dy', 14), ('dy0', 14),
    ('pitch', 15), ('dx1', 15),
    ('yaw', 16), ('dx1', 16),
    ('roll', 17)
))

ddtl = OrderedDict((
    ('scaling', 5), ('scale', 5), ('scl_fac', 5),
    ('frequency', 6), ('f', 6), ('freq', 6),
    ('driven phase', 7), ('phi0', 7),
    ('fileid', 8), ('file', 8), ('id', 8),
    ('aperture', 9), ('aper', 9), ('pipe', 9), ('radius', 9),
    ('quad 1 length', 10), ('q1len', 10),
    ('quad 1 gradient', 11), ('q1grad', 11),
    ('quad 2 length', 12), ('q2len', 12),
    ('quad 2 gradient', 13), ('q2grad', 13),
    ('dx', 14), ('dx0', 14),
    ('dy', 15), ('dy0', 15),
    ('pitch', 16), ('dx1', 16),
    ('yaw', 17), ('dx1', 17),
    ('roll', 18)
))

dccdtl = OrderedDict((
    ('scaling', 5), ('scale', 5), ('scl_fac', 5),
    ('frequency', 6), ('f', 6), ('freq', 6),
    ('driven phase', 7), ('phi0', 7),
    ('fileid', 8), ('file', 8), ('id', 8),
    ('aperture', 9), ('aper', 9), ('pipe', 9), ('radius', 9),
    ('dx', 10), ('dx0', 10),
    ('dy', 11), ('dy0', 11),
    ('pitch', 12), ('dx1', 12),
    ('yaw', 13), ('dy1', 13),
    ('roll', 14)
))

dccl = OrderedDict((
    ('scaling', 5), ('scale', 5), ('scl_fac', 5),
    ('frequency', 6), ('f', 6), ('freq', 6),
    ('input phase', 7), ('phi0', 7),
    ('fileid', 8), ('file', 8), ('id', 8),
    ('aperture', 9), ('aper', 9), ('pipe', 9), ('radius', 9),
    ('dx', 10), ('dx0', 10),
    ('dy', 11), ('dy0', 11),
    ('pitch', 12), ('dx1', 12),
    ('yaw', 13), ('dy1', 13),
    ('roll', 14),
    ('synchronous phase flag', 15), ('sync_flag', 15), ('syncflag', 15),
    ('scale error', 16), ('dscl', 16), ('er_scl', 16),
    ('phase error', 17), ('dphi', 17), ('er_phi', 17),
    ('scaling2', 18), ('scale2', 18), ('scl_fac2', 18),
    ('frequency ratio 2', 19), ('f2', 19), ('freq2', 19),
    ('driven phase offset 2', 20), ('phi0_2', 20),
    ('scaling3', 21), ('scale3', 21), ('scl_fac3', 21),
    ('frequency ratio 3', 22), ('f3', 22), ('freq3', 22),
    ('driven phase offset 3', 23), ('phi0_3', 23),
    ('scale error 2', 24), ('dscl2', 24), ('er_scl2', 24),
    ('phase error 2', 25), ('dphi2', 25), ('er_phi2', 25),
    ('scale error 3', 26), ('dscl3', 26), ('er_scl3', 26),
    ('phase error 3', 27), ('dphi3', 27), ('er_phi3', 27),
    ('synchronous phase offset', 28), ('sync_offset', 28)
))

dsc = OrderedDict((
    ('scaling', 5), ('scale', 5), ('scl_fac', 5),
    ('frequency', 6), ('f', 6), ('freq', 6),
    ('driven phase', 7), ('phi0', 7),
    ('fileid', 8), ('file', 8), ('id', 8),
    ('aperture', 9), ('aper', 9), ('pipe', 9), ('radius', 9),
    ('dx', 10), ('dx0', 10),
    ('dy', 11), ('dy0', 11),
    ('pitch', 12), ('dx1', 12),
    ('yaw', 13), ('dy1', 13),
    ('roll', 14)
))

dsolrf = OrderedDict((
    ('scaling', 5), ('scale', 5), ('scl_fac', 5),
    ('frequency', 6), ('f', 6), ('freq', 6),
    ('driven phase', 7), ('phi0', 7),
    ('fileid', 8), ('file', 8), ('id', 8),
    ('aperture', 9), ('aper', 9), ('pipe', 9), ('radius', 9),
    ('dx', 10), ('dx0', 10),
    ('dy', 11), ('dy0', 11),
    ('pitch', 12), ('dx1', 12),
    ('yaw', 13), ('dy1', 13),
    ('roll', 14),
    ('bz0', 15), ('bz', 15), ('b', 15)
))

drfq = OrderedDict((
    ('scaling', 5), ('scale', 5), ('scl_fac', 5),
    ('frequency', 6), ('f', 6), ('freq', 6),
    ('driven phase', 7), ('phi0', 7),
    ('fileid', 8), ('file', 8), ('id', 8),
    ('aperture', 9), ('aper', 9), ('pipe', 9), ('radius', 9),
    ('modulation', 10), ('mod', 10),
    ('dx', 11), ('dx0', 11),
    ('dy', 12), ('dy0', 12),
    ('pitch', 13), ('dx1', 13),
    ('yaw', 14), ('dy1', 14),
    ('roll', 15)
))


demfld = OrderedDict((
    ('scaling', 5), ('scale', 5), ('scl_fac', 5),
    ('frequency', 6), ('f', 6), ('freq', 6),
    ('input phase', 7), ('phi0', 7),
    ('fileid', 8), ('file', 8), ('id', 8),
    ('xaperture', 9), ('xaper', 9), ('xpipe', 9), ('xradius', 9),
    ('yaperture', 10), ('yaper', 10), ('ypipe', 10), ('yradius', 10),
    ('dx', 11), ('dx0', 11),
    ('dy', 12), ('dy0', 12),
    ('pitch', 13), ('dx1', 13),
    ('yaw', 14), ('dy1', 14),
    ('roll', 15),
    ('data type', 16), ('data', 16),
    ('coordinate', 17), ('coor', 17),
    ('synchronous phase flag', 18), ('sync_flag', 18), ('syncflag', 18),
    ('scale error', 19), ('dscl', 19), ('er_scl', 19),
    ('phase error', 20), ('dphi', 20), ('er_phi', 20),
    ('scaling2', 21), ('scale2', 21), ('scl_fac2', 21),
    ('frequency2', 22), ('f2', 22), ('freq2', 22),
    ('driven phase offset 2', 23), ('phi0_2', 23),
    ('scaling3', 24), ('scale3', 24), ('scl_fac3', 24),
    ('frequency3', 25), ('f3', 25), ('freq3', 25),
    ('driven phase offset 3', 26), ('phi0_3', 26),
    ('scale error 2', 27), ('dscl2', 27), ('er_scl2', 27),
    ('phase error 2', 28), ('dphi2', 28), ('er_phi2', 28),
    ('scale error 3', 29), ('dscl3', 29), ('er_scl3', 29),
    ('phase error 3', 30), ('dphi3', 30), ('er_phi3', 30),
    ('synchronous phase offset', 31), ('sync_offset', 31)
))

dbndfld = OrderedDict((
    ('scaling', 5), ('scale', 5), ('scl_fac', 5),
    ('angle', 6), ('phi', 6),
    ('fringe length', 7), ('flen', 7),
    ('fileid', 8), ('file', 8), ('id', 8),
    ('xaperture', 9), ('xaper', 9), ('xpipe', 9), ('xwidth', 9),
    ('yaperture', 10), ('yaper', 10), ('ypipe', 10), ('ywidth', 10),
    ('dx', 11), ('dx0', 11),
    ('dy', 12), ('dy0', 12),
    ('pitch', 13), ('dx1', 13),
    ('yaw', 14), ('dy1', 14),
    ('roll', 15),
    ('xoffset', 16), ('xoff', 16),
    ('yoffset', 17), ('yoff', 17),
))


def vrdct(typ):
    if typ == -27:
        ret = dfjump
    elif typ == -21:
        ret = dcor
    elif typ == -14:
        ret = dbps
    elif typ == -13:
        ret = dslit
    elif typ < 0:
        ret = dopt
    elif typ == 0:
        ret = ddrift
    elif typ == 1:
        ret = dquad
    elif typ == 2:
        ret = dcfcs
    elif typ == 3:
        ret = dsol
    elif typ == 4:
        ret = ddip
    elif typ == 5:
        ret = dmlt
    elif typ == 13:
        ret = dsol2
    elif typ == 101:
        ret = ddtl
    elif typ == 102:
        ret = dccdtl
    elif typ == 103:
        ret = dccl
    elif typ == 104:
        ret = dsc
    elif typ == 105:
        ret = dsolrf
    elif typ == 106:
        ret = drfq
    elif typ == 110:
        ret = demfld
    elif typ == 114:
        ret = dbndfld
    else:
        ret = -1
    return ret


def vrlen(typ):
    if typ == -27:
        ret = 3
    if typ == -14:
        ret = 4
    elif typ == -13:
        ret = 4
    elif typ < 0:
        ret = 7
    elif typ == 0:
        ret = 1
    elif typ == 1:
        ret = 8
    elif typ == 2:
        ret = 4
    elif typ == 3:
        ret = 8
    elif typ == 4:
        ret = 10
    elif typ == 5:
        ret = 13
    elif typ == 13:
        ret = 10
    elif typ == 101:
        ret = 24
    elif typ == 102:
        ret = 10
    elif typ == 103:
        ret = 24
    elif typ == 104:
        ret = 10
    elif typ == 105:
        ret = 11
    elif typ == 106:
        ret = 11
    elif typ == 110:
        ret = 27
    elif typ == 114:
        ret = 13

    return ret
