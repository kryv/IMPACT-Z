from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import copy as cp
import impact._impactf as im
import logging
import numpy as np

from impact import const
from impact import util
from mpi4py import MPI

from collections import OrderedDict
from impact.input import Input
from impact.result import Result

__authors__ = "Kei Fukushima"
__copyright__ = "(c) 2017, Facility for Rare Isotope beams, Michigan State University"
__contact__ = "Kei Fukushima <fukushim@frib.msu.edu>"

logging.basicConfig()
_Logger = logging.getLogger(__name__)

mpicomm = MPI.COMM_WORLD
mpisize = mpicomm.Get_size()
mpirank = mpicomm.Get_rank()

try:
    basestring
except NameError:
    basestring = str

try:
    im.pyfunction.setup()
    version_fsrc = np.asarray([im.messaging.version_fsrc]).tostring().decode().strip()
except ImportError:
    _Logger.error('cannot setup from _impact.so')


class Sequence(Input, Result):
    """Main class of python API for Advanced IMPACT code.

    - **Attribute**

    .. autosummary::
        subdir
        seed

    - **Methods**

    .. autosummary::
        read
        load
        distribute
        run
        construct
        insert
        conf
        search
        output
        save
        checker

    Parameters
    ----------
    input : str (optional)
        Input file path.

    subdir : str (optional)
        Directory path for data files.
    """

    def __init__(self, input=None, subdir=None):
        if isinstance(input, basestring):
            self.read(input)
        else:
            self._filein = None

        if isinstance(subdir, basestring):
            self.subdir = subdir
        else:
            self._subdir = None

        self._zstart = None

    def __len__(self):
        return im.messaging.nblem

    @property
    def subdir(self):
        """str: Directory path for data files."""
        return self._subdir

    @subdir.setter
    def subdir(self, v):
        _set_subdir(v)
        self._subdir = util._dpath(v)

    @property
    def seed(self):
        """int: Seed value of pseudo random number generator."""
        return im.messaging.isrd

    @seed.setter
    def seed(self, v):
        im.messaging.isrd = int(v)

    def read(self, fname=None, debug=False):
        """Read input file and set simulation parameters.

        Parameters
        ----------
        fname : str
            Input file path.

        debug : bool
            Flag for debug mode. If True, skip over input parameter correctness checker.

        Notes
        -----
        This function is automatically called in case of ``input`` is defined
        in :py:class:`Sequence <impact.control.Sequence>`.
        """
        if isinstance(fname, basestring):
            self._filein = fname

        _read_input(self._filein)
        if not debug:
            self.checker()
        self.dtag = util.get_tag(fname)

    def load(self, dname=None, reload=False):
        """Load 3D field data

        Parameters
        ----------
        dname : str (optional)
            Directory path for data files. Default path is defined by :py:func:`subdir <impact.control.Sequence.subdir>`.

        reload : bool
            If true, force to reload all 3D field data.
        """
        if reload:
            _reset_data()
        _load_data(dname)

    def distribute(self, particles=None, unit='physical', dname=None, fname='partcl.data', **kws):
        """Load 3D field data

        Parameters
        ----------
        particles : (9, n) shape ndarray (optional)
            All particle information for the simulation. *n* is the number of particles.

            * row 1 : x positions [m]

            * row 2 : xp momentums [rad]

            * row 3 : y positions [m]

            * row 4 : yp momentums [rad]

            * row 5 : z positions [deg]

            * row 6 : zp momentums [MeV]

            * row 7 : charge to mass ratios [c^2/eV]

            * row 8 : charge per macro-particle weights [1]

            * row 9 : ID numbers [1]

        unit : str (optional)
            Unit information for *particles*.

            * 'physical' : Use physical units [m, rad, m, rad, deg, MeV]. (default)

            * 'impact' : Use impact internal units.

        dname : str (optional)
            Directory path for particle data file. Default path is defined by :py:func:`subdir <impact.control.Sequence.subdir>`.

        fname : str (optional)
            File name or path for particle data file. Default is 'partcl.data'

        Mass : float (optional)
            Mass of the reference particle (unit is [eV/c^2]).

        Energy : float (optional)
            Initial kinetic energy (unit is [eV]).

        Frequency : float (optional)
            Scaling frequency (unit is [Hz]).

        zstart : float (optional)
            Starting z point of the simulation.

        Phaseini : float (optional)
            Initial phase of the reference particle (unit is [rad]).

        Notes
        -----
        ``dname`` and ``fname`` are used in case of :py:func:`Disttype <impact.input.Input.Disttype>` equals to 19.
        """

        if 'Mass' in kws:
            self.Mass = kws.get('Mass')
        if 'Energy' in kws:
            self.Energy = kws.get('Energy')
        if 'Frequency' in kws:
            self.Frequency = kws.get('Frequency')
        if 'Phaseini' in kws:
            self.Phaseini = kws.get('Phaseini')

        self._zstart = kws.get('zstart', None)

        if isinstance(particles, np.ndarray):
            if unit == 'physical':
                if im.messaging.runseq < 2:
                    if not self._flagfreq:
                        raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")
                    if not self._flagmass:
                        raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
                    if not self._flagengy:
                        raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

                g0 = self.Energy/self.Mass + 1.0
                b0g0 = np.sqrt(g0*g0 - 1.0)
                scxl = const.clight/(2.0*np.pi*self.Frequency)
                cnv = np.ones(9, dtype=np.float64)
                cnv[0] = 1.0/scxl
                cnv[1] = b0g0
                cnv[2] = 1.0/scxl
                cnv[3] = b0g0
                cnv[4] = np.pi/180.0
                cnv[5] = 1.0e6/self.Mass
                part = particles*cnv.reshape(9, 1)
            elif unit.lower() == 'impact':
                part = particles
            else:
                KeyError('Unknown unit type: ' + str(unit))
        else:
            part = None

        _init_dist(part, dname=dname, fname=fname, **kws)

    def run(self, start=0, end=-1, zstart=None, Phaseini=None):
        """Run particle tracking simulation.

        Parameters
        ----------
        start : int or str (optional)
            Starting element of the particle tracking simulation.

        end : int or str (optional)
            Ending element of the particle tracking simulation.

        zstart : float (optional)
            Initial z position of the reference particle. Default is 0.0 .

        Phaseini : float (optional)
            Initial phase of the reference particle (unit is [rad]).

        """

        if isinstance(start, basestring):
            ini = self.dtag[start]
            if not isinstance(ini, (int, np.integer)):
                ini = min(ini)
        else:
            ini = int(start)

        if isinstance(end, basestring):
            fin = self.dtag[end]
            if not isinstance(fin, (int, np.integer)):
                fin = max(fin)
        else:
            fin = int(end)

        if zstart is None:
            im.messaging.zposini = float(self._zstart) if self._zstart is not None else float(0.0)
        else:
            im.messaging.zposini = float(zstart)

        if Phaseini is not None:
            self.Phaseini = float(Phaseini)

        _run_accel(ini, fin)

    def _cnv_element(self, i, el):
        if isinstance(el, dict):
            if 'length' in el:
                lng = el['length']
            elif 'l' in el:
                lng = el['l']
            elif 'L' in el:
                lng = el['L']
            elif 'start' in el and 'end' in el:
                lng = el['end'] - el['start']
            else:
                raise TypeError("Lattice element must have 'length' (or 'l', 'L') parameter.")

            if 'segment' in el:
                seg = el['segment']
            elif 'seg' in el:
                seg = el['seg']
            else:
                raise TypeError("Lattice element must have 'segment' (or 'seg') parameter.")

            if 'step' in el:
                stp = el['step']
            elif 'pid' in el:
                stp = el['pid']
            else:
                raise TypeError("Lattice element must have 'step' (or 'pid') parameter.")

            if 'type' in el:
                ttyp = el['type']
            else:
                raise TypeError("Lattice element must have 'type' parameter.")
            typ = const.elemtype[ttyp] if isinstance(ttyp, basestring) else ttyp

            idd = {}
            param = np.zeros(28)
            tdct = const.vrdct(typ)
            for kw in el:
                if kw.isdigit():
                    ikw = int(kw)
                    if ikw == 1:
                        lng = el[kw]
                    elif ikw == 2:
                        seg = el[kw]
                    elif ikw == 3:
                        stp = el[kw]
                    elif ikw == 4:
                        typ = el[kw]
                    else:
                        idd[kw] = el[kw]
                else:
                    if kw.lower() == 'tag':
                        if isinstance(el[kw], basestring):
                            self.dtag = util.add_tag(self.dtag, {el[kw]: i+1})
                        else:
                            for tag in el[kw]:
                                self.dtag = util.add_tag(self.dtag, {tag: i+1})
                    elif kw.lower() in tdct:
                        kid = tdct[kw.lower()]
                        pkw = util.invdict(tdct, kid)
                        param[kid-4] = el[pkw] if pkw in el else el[kw]

            # overwrite by index input
            for id in idd:
                param[int(id)] = idd[id]

        elif isinstance(el, (list, tuple)):
            lng = el[0]
            seg = el[1]
            stp = el[2]
            typ = el[3]
            esize = len(el[4:])
            param = np.zeros(28)
            param[0] = -1.0
            if isinstance(el[-1], basestring):
                self.dtag = util.add_tag(self.dtag, {el[-1]: i+1})
                param[1:esize] = el[4:-1]
            elif isinstance(el[-1], (list, tuple)):
                for tag in el[-1]:
                    self.dtag = util.add_tag(self.dtag, {tag: i+1})
                param[1:esize] = el[4:-1]
            else:
                param[1:esize+1] = el[4:]
        else:
            raise TypeError('Element information must be a dict or list.')

        return lng, seg, stp, typ, param

    def construct(self, lattice):
        """Construct beam transport line.

        Parameters
        ----------
        lattice : 2D array or list of dict
            List of lattice elements.

        Examples
        --------
        >>> drift = [0.5, 10, 10, 0, 0.2]
        # define element by using list of numbers (format is the same as the input file)
        >>> quad1 = {'length':0.2, 'seg':4, 'step':4, 'type':'quadrupole', 'B2': 10.0, 'aper':0.2}
        # define element by using python dictionary
        >>> quad2 = quad1.copy()
        >>> quad2['B2'] = -10.0
        >>> FODO = [drift, quad1, drift, quad2, drift]
        >>> sq.construct(FODO)
        """

        nl = len(lattice)
        lngs = np.zeros(nl)
        segs = np.zeros(nl, dtype=np.int32)
        stps = np.zeros(nl, dtype=np.int32)
        typs = np.zeros(nl, dtype=np.int32)
        params = np.zeros((28, nl))
        params[0, :] = -1
        self.dtag = OrderedDict()

        for i, el in enumerate(lattice):
            lng, seg, stp, typ, param = self._cnv_element(i, el)

            lngs[i] = lng
            segs[i] = seg
            stps[i] = stp
            typs[i] = typ
            params[:, i] = param

        _construct(lngs, segs, stps, typs, params)

    def insert(self, index, element):
        """Insert new lattice element.

        Parameters
        ----------
        index : int
            Insert new lattice element before index. If index equals to 0, *append* the new element.

        element : list or dict
            Lattice element information.

        Examples
        --------
        >>> quad1 = {'length':0.2, 'seg':4, 'step':4, 'type':'quadrupole', 'B2': 10.0, 'aper':0.2}
        # define element by using python dictionary
        >>> sq.construct(3, quad1) # insert new element
        """
        if index == 0:
            idx = im.messaging.nblem
        else:
            idx = index - 1

        def addnum(idx, tmp):
            new = tmp + 1 if tmp > idx else tmp
            return new

        for ky in self.dtag:
            tmpids = self.dtag[ky]
            if isinstance(tmpids, (int, np.integer)):
                newids = addnum(idx, tmpids)
            else:
                newids = [addnum(idx, tmpid) for tmpid in tmpids]

            self.dtag[ky] = newids

        lng, seg, stp, typ, param = self._cnv_element(idx, element)

        tlngs = list(im.messaging.blnlen)
        tlngs.insert(idx, lng)
        lngs = np.asarray(tlngs)

        tsegs = list(im.messaging.blnseg)
        tsegs.insert(idx, seg)
        segs = np.asarray(tsegs)

        tstps = list(im.messaging.blnstp)
        tstps.insert(idx, stp)
        stps = np.asarray(tstps)

        ttyps = list(im.messaging.blntyp)
        ttyps.insert(idx, typ)
        typs = np.asarray(ttyps)

        tparams = list(map(list, list(im.messaging.blnparams)))
        for tp, pp in zip(tparams, param):
            tp.insert(idx, pp)
        params = np.asarray(tparams)

        _construct(lngs, segs, stps, typs, params)

    def conf(self, index=None, col=None, value=None, syntax=True):
        """Configure lattice parameters.

        Parameters
        ----------
        index : int or str
            Index number of the element or tag name of the element.

        col : int or str
            Column number of the target parameter or name of the target parameter.

        col : dict
            dict key for the target parameter and dict value for target value.

        value : float
            New value for the target parameter.

        syntax : bool
            Print flag for full parameter of the element.

        Returns
        -------
        dict
            Dictionary style parameter list of the element

        Examples
        --------
        >>> sq = Sequence('test.in')
        >>> tmp = sq.conf(3) # show parameter of 3rd element
        index :  3
        length :  0.24 [m]
        segment :  60
        step :  20
        type :  110  (emfld)
        scaling ('scale', 'scl_fac') :  0.64
        frequency ('f', 'freq') :  80500000.0 [Hz]
        driven phase ('phi0',) :  349.74086645 [deg]
        fileid ('file', 'id') :  421.0
        xaperture ('xaper', 'xpipe', 'xradius') :  0.017 [m]
        yaperture ('yaper', 'ypipe', 'yradius') :  0.017 [m]
        dx  :  0.0 [m]
        dy  :  0.0 [m]
        pitch  :  0.0 [rad]
        yaw  :  0.0 [rad]
        roll  :  0.0 [rad]
        data type ('data',) :  1.0
        coordinate ('coor',) :  2.0
        tag :  ['cav1']
        >>> sq.conf(3, {'phi0':360.0}) # set new value for the driven phase
        >>> tmp = sq.conf(3, 'phi0') # confirm new parameter
        driven phase ('phi0',) :  360.0 [deg]

        If you have set the tag, followings return same results.

        >>> sq.conf('cav1') # show parameter of 3rd element
        >>> sq.conf('cav1', {'phi0':360.0}) # set new value for the driven phase
        """

        if index is None:
            return _conf_para(range(1, im.messaging.nblem+1), None, None, self.dtag, False)
        else:
            if isinstance(col, dict):
                if ('tag' in col):
                    if isinstance(col['tag'], basestring):
                        self.dtag = util.add_tag(self.dtag, {col['tag']: index})
                    else:
                        for tag in col['tag']:
                            self.dtag = util.add_tag(self.dtag, {tag: index})
                    return None

            return _conf_para(index, col, value, self.dtag, syntax)

    def search(self, type=None, tag=None):
        """Search lattice element.

        Parameters
        ----------
        type : int or str
            Element type ID or type name for _search.

        tag : str
            Tag name of the element.

        Returns
        -------
        ndarray
            List of element ID

        Notes
        -----
        In case of input both type and tag, it returns `AND` _search result.
        """
        if type is not None and tag is None:
            return _search(type)
        elif type is None and isinstance(tag, basestring):
            return np.asarray([self.dtag[tag]]).flatten()
        elif type is not None and isinstance(tag, basestring):
            lstyp = _search(type)
            lstag = np.asarray([self.dtag[tag]]).flatten()
            return np.intersect1d(lstyp, lstag)

    def output(self, fname='fort', qid=-1):
        """Output original IMPACT results to file.

        Parameters
        ----------
        fname : str
            Header name of output file. Default is 'fort'.

        qid : int
            Charge state index for return value. Default is -1 (Total of the all charge states).
            Index is compatible with :py:func:`qmlabel <impact.result.Result.qmlabel>`.
        """
        if qid > len(self.qmlabel):
            raise ValueError('Input qid is out of range')
        elif qid == 0:
            iid = 1
        else:
            iid = qid + 2

        _output_hist(fname=fname, iqid=iid)
        _output_dist(fname)

    def save(self, fname='test.in.new'):
        """Save current setting to file.

        Parameters
        ----------
        fname : str
            File name for saved settings. Default is 'test.in.new'.
        """
        _save_testin(fname)
        util.put_tag(fname, self.dtag)

    def checker(self):
        """Check input parameter correctness."""
        if mpisize != self.Mpisize:
            _Logger.warning('Number of proccessors does not match to the input value. Adjust to actual size.')
            self.Mpisize = mpisize

        self.Flaginteg = self.Flaginteg
        self.Flagdiag = self.Flagdiag
        self.Flagbc = self.Flagbc
        self.Nxgrid = self.Nxgrid
        self.Nygrid = self.Nygrid
        self.Nzgrid = self.Nzgrid

        if self.Nptot <= 0.0:
            raise ValueError('Nptot must be larger than zero.')

        if self.Nptot != np.sum(self.Nplist):
            raise ValueError('Nptot does not match to the sum of Nplist.')

        if self.Frequency <= 0.0:
            raise ValueError('Frequency must be larger than zero.')

        if self.Mass <= 0.0:
            raise ValueError('Mass must be larger than zero.')

        if self.Disttype not in [1, 2, 3, 4, 5, 16, 17, 19]:
            raise ValueError('Wrong input value for Disttype.')

        self.Distxquadratic = self.Distxquadratic
        self.Distyquadratic = self.Distyquadratic
        self.Distzquadratic = self.Distzquadratic


def _read_input(fname):
    """Read input file"""
    if im.messaging.runseq < 1:
        raise Exception('setup() must be called before _read_input()')
    elif isinstance(fname, basestring):
        im.pyfunction.input(fname)
        _iostat(['inp'])
    else:
        raise TypeError('build() takes a str or unicode argument for input file')


def _set_subdir(dname):
    """Set all sub directory for IMPACT"""
    if isinstance(dname, basestring):
        im.pyfunction.set_datadir(util._dpath(dname))
        im.pyfunction.set_partdir(util._dpath(dname))
        im.pyfunction.set_strpdir(util._dpath(dname))
    elif (dname is not None):
        _Logger.warning('Input argument must be a str or unicode. Input argument is ignored.')


def _load_data(dname=None):
    """Load 3D field data"""
    if isinstance(dname, basestring):
        im.pyfunction.set_datadir(util._dpath(dname))
    elif (dname is not None):
        _Logger.warning('Input argument must be a str or unicode. Input argument is ignored.')

    im.pyfunction.load_data()
    _iostat(['data'])


def _reset_data():
    im.pyfunction.reset_data()


def _init_dist(particles=None, dname=None, fname='partcl.data', **kws):
    """Initialize particle distribution"""
    if im.messaging.runseq < 2:
        raise Exception('input parameters must be set before _init_dist()')

    if isinstance(particles, np.ndarray):
        i, j = np.shape(particles)
        if i != 9:
            raise ValueError('shape mismatch: input particles must have (9, n) shape')
        im.messaging.idistnp = j
        im.messaging.idistd = particles
        unq = np.unique(particles[6, :])
        im.messaging.nchrg = len(unq)
        tmp1 = np.zeros(100)
        tmp1[0:len(unq)] = np.asarray([np.sum(particles[6, :] == v) for v in unq])
        im.messaging.nptlist = tmp1
        tmp2 = np.zeros(100)
        tmp2[0:len(unq)] = unq
        im.messaging.qmcclist = tmp2
    else:
        im.messaging.idistnp = -1

    if isinstance(dname, basestring):
        im.pyfunction.set_partdir(util._dpath(dname))
    elif (dname is not None):
        _Logger.warning('Input argument must be a str or unicode. Input argument is ignored.')

    if isinstance(fname, basestring):
        im.pyfunction.set_partfile(fname)
    elif (fname is not None):
        _Logger.warning('Input argument must be a str or unicode. Input argument is ignored.')

    flagsc = kws.get('Flagsc', im.messaging.bcurr)
    mass = kws.get('Mass', im.messaging.bmass)
    eng = kws.get('Energy', im.messaging.bkenergy)
    phsini = kws.get('Phaseini', im.messaging.phsini)
    im.pyfunction.init_dist(flagsc, mass, eng, phsini)
    _iostat(['dist'])


def _run_accel(ini=0, fin=0):
    """Run particle tracking simulation"""
    if im.messaging.runseq < 2:
        raise Exception('Input parameters must be set before _init_dist()')
    if im.messaging.dstseq == 0:
        raise Exception('Initial distribution must be set before _run_accel()')

    im.pyfunction.run(int(ini), int(fin))
    _iostat(['data', 'strp'])


def _construct(lngs, segs, stps, typs, params):
    """Construct Beam line"""

    n, m = params.shape
    if list(map(len, [lngs, segs, stps, typs])) != [m]*4:
        raise TypeError('Input data size does not match.')
    if m <= 0:
        raise ValueError('Number of Lattice elements must be greater than 0.')
    if n != 28:
        raise TypeError('Parameter array size does not match.')

    im.messaging.nblem = np.asarray(m, dtype=np.int32)
    im.messaging.blnlen = np.asarray(lngs)
    im.messaging.blnseg = np.asarray(segs, dtype=np.int32)
    im.messaging.blnstp = np.asarray(stps, dtype=np.int32)
    im.messaging.blntyp = np.asarray(typs, dtype=np.int32)
    im.messaging.blnparams = np.asarray(params)
    im.pyfunction.construct()


def _conf_para(index, column=None, value=None, idic=None, syntax=True):
    """Configure lattice parameters"""
    if isinstance(index, basestring):
        if isinstance(idic, dict):
            index = idic[index]
        else:
            raise TypeError('wrong input arguments for _conf_para():', index, idic)

    if isinstance(index, (int, np.integer)):
        index = [index]

    if isinstance(column, dict):
        col = list(column.keys())
        val = list(column.values())
    else:
        col = [column]
        val = [value]

    maxid = im.messaging.nblem
    if (val[0] is None):
        ret = []
        for idx in index:
            for c in col:
                if idx > maxid:
                    raise IndexError('Input index is out of lattice range (max '+str(im.messaging.nblem)+').')

                lng, seg, stp, typ, vrr = im.pyfunction.configure0(idx)
                zs, ze = im.messaging.blnpos[idx-1:idx+1]

                dline = list([vrr[v] for v in range(const.vrlen(typ))])

                if typ == -2:
                    keystp = 'pid'
                else:
                    keystp = 'step'

                ddct = OrderedDict((('index', idx), ('length', lng), ('segment', seg),
                                    (keystp, stp), ('type', typ), ('start', zs), ('end', ze)))

                if syntax:
                    if c in [None, 0, 'index']:
                        print('index : ', idx, '(start : ', zs, ', end : ', ze, ' [m])')
                    if c in [None, 1, 'length']:
                        print('length : ', lng, '[m]')
                    if c in [None, 2, 'segment']:
                        print('segment : ', seg)
                    if c in [None, 3, keystp]:
                        print(keystp+' : ', stp)

                    if typ in const.elemtype.values():
                        if c in [None, 4, 'type']:
                            print('type : ', typ, ' ('+util.invdict(const.elemtype, typ)+')')
                    else:
                        if c in [None, 4, 'type']:
                            print('type : ', typ)

                delem = const.vrdct(typ)
                if delem != -1:
                    for i, vr in enumerate(dline):
                        vname = util.invdict(const.vrdct(typ), i+5)
                        alias = tuple(util.invdictall(const.vrdct(typ), i+5)[1:])
                        info = '' if len(alias) == 0 else alias
                        unit = const.uelem[vname] if vname in const.uelem else ''
                        ddct.update(((vname, vr),))
                        for vn in alias:
                            ddct.update(((str(vn), vr),))
                        if syntax and (c in [None, i+5, vname] or c in alias):
                            print(vname, info, ': ', vr, unit)

                if idic is not None:
                    tagname = util.invdictls(idic, idx, str)
                    ddct.update((('tag', tagname),))
                    if syntax and c in [None, 0, 'tag']:
                        print('tag : ', list(map(str, tagname)))

                if syntax:
                    print('')
                ret.append(ddct)
        ret = ret[0] if len(ret) == 1 else ret

    elif (col is not None and val is not None):
        ret = None
        for idx in index:
            for c, v in zip(col, val):
                if idx > maxid:
                    raise IndexError('Input index is out of lattice range (max '+str(im.messaging.nblem)+').')

                lng, seg, stp, typ, vrr = im.pyfunction.configure0(idx)

                if isinstance(c, basestring):
                    ccol = const.vrdct(typ)[c.lower()]
                elif isinstance(c, (int, np.integer)):
                    ccol = c
                else:
                    raise TypeError('wrong input arguments for _conf_para() :', c)

                if ccol == 4:
                    raise IndexError('Can not change element type by conf().')
                if ccol > 31:
                    raise IndexError('Input col must be less than or equal to 31.')

                im.pyfunction.configure1(idx, ccol, v)
    else:
        raise TypeError('wrong input arguments for _conf_para()')

    return ret


def _search(typ):
    """Search lattice elements"""
    if isinstance(typ, basestring):
        index = const.elemtype[typ]
    elif isinstance(typ, (int, np.integer)):
        index = typ
    else:
        raise TypeError('wrong input arguments for _search() :', str(typ))

    im.pyfunction.search(index)
    mpicomm.barrier()
    return cp.copy(im.messaging.searchres)


def _save_testin(fname='test.in.new'):
    """Save current setting to file"""
    if im.messaging.runseq < 2:
        raise Exception('Input parameters must be set before exporting.')
    elif isinstance(fname, basestring):
        im.pyfunction.export(fname)
    else:
        raise TypeError('wrong input arguments for _save_testin()')


def _output_hist(fname='fort', iqid=1):
    """Output history data to file"""
    if im.messaging.hrefz is None:
        raise Exception('no exportable data')
    elif isinstance(fname, basestring):
        im.pyfunction.outhist(fname, iqid)
    else:
        raise TypeError('wrong input arguments for _output_hist()')


def _output_dist(fname='fort'):
    """Output distribution data to file"""
    if (im.messaging.hbunchid is None) or (len(im.messaging.hbunchid) == 0):
        raise Exception('no exportable data')
    elif isinstance(fname, basestring):
        im.pyfunction.outdist(fname)
    else:
        raise TypeError('wrong input arguments for _output_dist()')


def _iostat(key):
    """Return io status in IMPACT"""
    iolst = im.pyfunction.get_iostat()
    iolst = mpicomm.bcast(iolst, root=0)

    if iolst[0] != 0 and 'inp' in key:
        raise Exception('Input file does not found.')
    if iolst[1] != 0 and 'dist' in key:
        raise Exception('Initial distribution file does not found.')
    if iolst[2] != 0 and 'data' in key:
        raise Exception('Field data file does not found.')
    if iolst[3] != 0 and 'strp' in key:
        raise Exception('Stripper data file does not found.')
