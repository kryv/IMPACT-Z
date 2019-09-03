from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import impact._impactf as im
import re
import numpy as np

from mpi4py import MPI

from collections import OrderedDict

__authors__ = "Kei Fukushima"
__copyright__ = "(c) 2017, Facility for Rare Isotope beams, Michigan State University"
__contact__ = "Kei Fukushima <fukushim@frib.msu.edu>"

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

try:
    basestring
except NameError:
    basestring = str


def _usup2d(x, f, unit0, dcnv, abc, *args, **kws):
    unit = kws.get('unit', unit0)
    csid = kws.get('qid', -1) + 1
    idx = None
    value = None
    zstt = im.messaging.zposini

    if rank == 0:
        csx = 0 if csid >= x.shape[0] else csid
        stt, end = im.messaging.hchrgpos[csx]

        if len(x.shape) == 1:
            xtmp = x[stt-1:end] + zstt
        else:
            xtmp = x[csx, stt-1:end] + zstt

        if len(f.shape) == 1:
            ftmp = f[stt-1:end]
        else:
            ftmp = f[csx, stt-1:end]

        for elem in args:
            if isinstance(elem, str) and (elem in dcnv):
                unit = elem
            elif isinstance(elem, dict):
                idx = elem.get('id')
            elif elem is not None:
                value = np.asarray([elem]).flatten() + zstt

        if value is None and idx is None:
            ret = ftmp*dcnv[unit]
        elif value is None and idx is not None:
            ret = ftmp[idx]*dcnv[unit]
        else:
            ret = np.interp(value, xtmp, ftmp)*dcnv[unit]
    else:
        ret = None

    if abc and size != 1:
        ret = comm.bcast(ret, root=0)

    return ret


def _dpath(dname):
    dname = dname.strip()
    if isinstance(dname, basestring):
        if (len(dname) == 0):
            path = ''
        elif (dname[-1] == '/'):
            path = dname
        else:
            path = dname+'/'
    else:
        path = ''
    return path


def invdict(dct, val):
    """Inverse search from dictionary. Returns first key."""
    return list(dct.keys())[list(dct.values()).index(val)]


def invdictall(dct, val):
    """Inverse search from dictionary. Returns all keys."""
    return list(map(str, np.asarray(list(dct.keys()))[np.asarray(list(dct.values())) == val]))


def invdictls(dct, val, enc):
    """Inverse search from complex dictionary. Returns all keys."""
    ntag = []
    for (j, vl) in enumerate(dct.values()):
        if isinstance(vl, (int, np.integer)):
            vl = [vl]
        if val in vl:
            ntag.append(enc(list(dct.keys())[j]))
    return ntag


def get_tag(fname, skip=11):
    """Get tag dictionary from IMPACT input file."""
    dtag = OrderedDict({})
    idx = 1 - skip

    if rank == 0:
        with open(fname, 'r') as f:
            lines = f.readlines()
        for ln in lines:
            fln = ln.replace(' ', '')
            if fln[0] != '!':
                loc = fln.find('/tag=')
                if loc != -1:
                    for blc in ['{}', '[]', '()']:
                        if fln[loc+5] == blc[0]:
                            stt = loc+6
                            end = stt+fln[stt:].find(blc[1])
                            tag = None if end < stt else eval(fln[stt:end])
                            if isinstance(tag, basestring):
                                dtag = add_tag(dtag, {tag: idx})
                            elif tag is not None:
                                for tg in tag:
                                    dtag = add_tag(dtag, {tg: idx})

                idx += 1

    dtag = comm.bcast(dtag, root=0)
    return dtag


def add_tag(dtag, newtag):
    for tag in newtag.keys():
        if tag in dtag:
            if isinstance(dtag[tag], (int, np.integer)):
                if not newtag[tag] == dtag[tag]:
                    dtag[tag] = [dtag[tag], newtag[tag]]
            else:
                if not newtag[tag] in dtag[tag]:
                    dtag[tag].append(newtag[tag])
        else:
            dtag[tag] = newtag[tag]

    return dtag


def put_tag(fname, dtag, skip=11):
    """Put tag for IMPACT input file."""
    if rank == 0:
        with open(fname, 'r') as f:
            lines = f.readlines()

        f2 = open(fname, 'w')
        for i, ln in enumerate(lines):
            ln = re.sub(r'0+E', '0E', ln)
            ln = re.sub(r'E\+00', '', ln)
            if i >= skip:
                ntag = invdictls(dtag, i-skip+1, str)
                if len(ntag) != 0:
                    f2.write(ln[:-1]+' tag='+str(ntag)+ln[-1:])
                else:
                    f2.write(ln)
            else:
                f2.write(ln)
        f2.close()
