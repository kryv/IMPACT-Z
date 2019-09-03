from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import impact._impactf as im
import logging
import numpy as np

from impact import const

__authors__ = "Kei Fukushima"
__copyright__ = "(c) 2017, Facility for Rare Isotope beams, Michigan State University"
__contact__ = "Kei Fukushima <fukushim@frib.msu.edu>"

logging.basicConfig()
_Logger = logging.getLogger(__name__)

try:
    basestring
except NameError:
    basestring = str


class Input(object):
    """Base class for input parameters of Advanced IMPACT.

    - **Attributes for simulation parameter**

    .. autosummary::
        Mpisize
        Ndim
        Nxgrid
        Nygrid
        Nzgrid
        Xbound
        Ybound
        Zperiod
        Flaginteg
        Flagerror
        Flagdiag
        Flagbc
        Flagsubstep
        Flagrestart
        Flagsc
        Outputemit

    - **Attributes for beam parameter**

    .. autosummary::
        Nptot
        Ncharge
        Nplist
        Curlist
        Ctmlist
        Energy
        Mass
        Frequency
        Chargestate
        Phaseini
        Disttype
        Distxsigma
        Distysigma
        Distzsigma
        Distxlambda
        Distylambda
        Distzlambda
        Distxmu
        Distymu
        Distzmu
        Distxtwiss
        Distytwiss
        Distztwiss
        Distxquadratic
        Distyquadratic
        Distzquadratic
        Distxoffset
        Distyoffset
        Distzoffset
        Distxpoffset
        Distypoffset
        Distzpoffset
        Distxoffset_m
        Distyoffset_m
        Distzoffset_deg
        Distxpoffset_rad
        Distypoffset_rad
        Distzpoffset_eV
        Distxmismatch
        Distymismatch
        Distzmismatch
        Distxpmismatch
        Distypmismatch
        Distzpmismatch
    """

    _flagmass = False
    _flagfreq = False
    _flagengy = False
    _xtwiss = np.zeros(3)
    _ytwiss = np.zeros(3)
    _ztwiss = np.zeros(3)

    @property
    def Mpisize(self):
        """int: Number of processors in MPI communicator."""
        return im.messaging.npcol*im.messaging.nprow

    @Mpisize.setter
    def Mpisize(self, v):
        im.messaging.npcol = int(v)
        im.messaging.nprow = 1

    @property
    def Ndim(self):
        """int: Dimension of the simulation space (dummy parameter)."""
        return im.messaging.dm

    @Ndim.setter
    def Ndim(self, v):
        im.messaging.dm = int(v)

    @property
    def Nptot(self):
        """int: Total number of macro particles."""
        return im.messaging.npin

    @Nptot.setter
    def Nptot(self, v):
        im.messaging.npin = int(v)

    @property
    def Flaginteg(self):
        """int: Flag for integration method. 1 for linear map, 2 for nonlinear Lorentz integrator."""
        return im.messaging.flagmap

    @Flaginteg.setter
    def Flaginteg(self, v):
        if not int(v) in [1, 2]:
            raise ValueError("Input value must be 1 or 2.")
        else:
            im.messaging.flagmap = int(v)

    @property
    def Flagerror(self):
        """int: Flag for error study.

        Notes
        -----
        - 0 : No error simulation

        - 1 : Use ``dx``, ``dy``, ``pitch``, ``yaw`` and ``roll`` errors to the lattice element.

        - 2 : Use entrance offsets (``dx0``, ``dy0``), exit offsets (``dx1``, ``dy1``), and ``roll`` errors to the lattice element.
        """
        return im.messaging.flagerr

    @Flagerror.setter
    def Flagerror(self, v):
        im.messaging.flagerr = bool(v)

    @property
    def Flagdiag(self):
        """int: Flag for diagnostic routine. 1 or 3 for standard, 2 or 4 for `n` % emittance, radius output.

        Notes
        -----
        - 1 or 2 : Diagnostic routine is called for each segment.

        - 3 or 4 : Diagnostic routine is called in case of element type = -23 (monitor).
        """
        return im.messaging.flagdiag

    @Flagdiag.setter
    def Flagdiag(self, v):
        if not int(v) in [1, 2, 3, 4]:
            raise ValueError("Input value must be 1, 2, 3, or 4.")
        else:
            im.messaging.flagdiag = int(v)

    @property
    def Nxgrid(self):
        """int: Number of x mesh grid for space charge (PIC) calculation.

        Notes
        -----
        The x mesh number must be :math:`2^n` (for :py:func:`Flagbc <impact.input.Input.Flagbc>` equals to 1, 2, 3, or 4)
        or :math:`2^n+1` (for :py:func:`Flagbc <impact.input.Input.Flagbc>` equals to 5 or 6).
        """
        return im.messaging.nx

    @Nxgrid.setter
    def Nxgrid(self, v):
        if (not np.log2(v).is_integer()) and (not np.log2(v-1).is_integer()):
            raise ValueError("Input value must be 2^n or 2^n + 1.")
        else:
            im.messaging.nx = int(v)

    @property
    def Nygrid(self):
        """int: Number of y mesh grid for space charge (PIC) calculation.

        Notes
        -----
        The y mesh number must be :math:`2^n` (for :py:func:`Flagbc <impact.input.Input.Flagbc>` equals to 1 or 2)
        or :math:`2^n+1` (for :py:func:`Flagbc <impact.input.Input.Flagbc>` equals to 3, 4, 5, or 6).
        """
        return im.messaging.ny

    @Nygrid.setter
    def Nygrid(self, v):
        if (not np.log2(v).is_integer()) and (not np.log2(v-1).is_integer()):
            raise ValueError("Input value must be 2^n or 2^n + 1.")
        else:
            im.messaging.ny = int(v)

    @property
    def Nzgrid(self):
        """int: Number of z mesh grid for space charge (PIC) calculation.

        Notes
        -----
        The z mesh number must be :math:`2^n` (for :py:func:`Flagbc <impact.input.Input.Flagbc>` equals to 1, 3, or 5)
        or :math:`2^n+1` (for :py:func:`Flagbc <impact.input.Input.Flagbc>` equals to 2, 4, or 6).
        """
        return im.messaging.nz

    @Nzgrid.setter
    def Nzgrid(self, v):
        if (not np.log2(v).is_integer()) and (not np.log2(v-1).is_integer()):
            raise ValueError("Input value must be 2^n or 2^n + 1.")
        else:
            im.messaging.nz = int(v)

    @property
    def Flagbc(self):
        """int: Flag of boundary condition for space charge (PIC) calculation.

        Notes
        -----
        - 1 : 3D open. In this case, all mesh number (:py:func:`Nxgrid <impact.input.Input.Nxgrid>`,
          :py:func:`Nygrid <impact.input.Input.Nygrid>`, :py:func:`Nzgrid <impact.input.Input.Nzgrid>`) must be :math:`2^n`.

        - 2 : Transverse open, longitudinal periodic. In this case,
          :py:func:`Nxgrid <impact.input.Input.Nxgrid>` and :py:func:`Nygrid <impact.input.Input.Nygrid>` must be :math:`2^n`,
          :py:func:`Nzgrid <impact.input.Input.Nzgrid>` must be :math:`2^n+1`.

        - 3 : Transverse finite, longitudinal open round pipe. In this case,
          :py:func:`Nxgrid <impact.input.Input.Nxgrid>` and :py:func:`Nzgrid <impact.input.Input.Nzgrid>` must be :math:`2^n`,
          :py:func:`Nygrid <impact.input.Input.Nygrid>` must be :math:`2^n+1`.

        - 4 : Transverse finite, longitudinal periodic round pipe. In this case
          :py:func:`Nxgrid <impact.input.Input.Nxgrid>` must be :math:`2^n`,
          :py:func:`Nygrid <impact.input.Input.Nygrid>` and :py:func:`Nzgrid <impact.input.Input.Nzgrid>` must be :math:`2^n+1`.

        - 5 : Transverse finite, longitudinal open rectangular pipe. In this case,
          :py:func:`Nzgrid <impact.input.Input.Nzgrid>` must be :math:`2^n`,
          :py:func:`Nxgrid <impact.input.Input.Nxgrid>` and :py:func:`Nygrid <impact.input.Input.Nygrid>` must be :math:`2^n+1`.

        - 6 : Transverse finite, longitudinal periodic rectangular pipe. In this case,
          all mesh number ( :py:func:`Nxgrid <impact.input.Input.Nxgrid>`, :py:func:`Nygrid <impact.input.Input.Nygrid>`,
          :py:func:`Nzgrid <impact.input.Input.Nzgrid>`) must be :math:`2^n+1`.
        """
        return im.messaging.flagbc

    @Flagbc.setter
    def Flagbc(self, v):
        if not int(v) in [1, 2, 3, 4, 5, 6]:
            raise ValueError("Input value must be in 1 to 6.")
        else:
            im.messaging.flagbc = int(v)

    @property
    def Xbound(self):
        """float: x (horizontal) pipe width for space charge (PIC) calculation (unit is [m])."""
        return im.messaging.xwallrad

    @Xbound.setter
    def Xbound(self, v):
        im.messaging.xwallrad = float(v)

    @property
    def Ybound(self):
        """float: y (vertical) pipe width for space charge (PIC) calculation (unit is [m])."""
        return im.messaging.ywallrad

    @Ybound.setter
    def Ybound(self, v):
        im.messaging.ywallrad = float(v)

    @property
    def Zperiod(self):
        """float: z periodic length for space charge (PIC) calculation (unit is [m])."""
        return im.messaging.perdlen

    @Zperiod.setter
    def Zperiod(self, v):
        im.messaging.perdlen = float(v)

    @property
    def Flagrestart(self):
        """bool: Flag for restarting the simulation (dummy parameter)."""
        return im.messaging.rstartflg

    @Flagrestart.setter
    def Flagrestart(self, v):
        im.messaging.rstartflg = bool(v)

    @property
    def Flagsubstep(self):
        """bool: Flag of sub-cycle for space charge calculation."""
        return im.messaging.flagsubstep

    @Flagsubstep.setter
    def Flagsubstep(self, v):
        im.messaging.flagsubstep = bool(v)

    @property
    def Ncharge(self):
        """int: Number of charge states"""
        return im.messaging.nchrgin

    @Ncharge.setter
    def Ncharge(self, v):
        im.messaging.nchrgin = int(v)

    @property
    def Nplist(self):
        """list of int: List of particle number for each charge state."""
        return im.messaging.nptlistin[0:self.Ncharge]

    @Nplist.setter
    def Nplist(self, v):
        if self.Ncharge != len(v):
            _Logger.warning('List length does not match Ncharge.')
        ls = np.zeros(100, dtype=int)
        ls[0:len(v)] = v
        im.messaging.nptlistin = ls

    @property
    def Curlist(self):
        """list of float: List of beam current for each charge state. (unit is [A])"""
        return im.messaging.currlistin[0:self.Ncharge]

    @Curlist.setter
    def Curlist(self, v):
        if self.Ncharge != len(v):
            _Logger.warning('List length does not match Ncharge.')
        ls = np.zeros(100)
        ls[0:len(v)] = v
        im.messaging.currlistin = ls

    @property
    def Ctmlist(self):
        """list of float: List of charge to mass ratio for each charge state. (unit is [c^2/eV])"""
        return im.messaging.qmcclistin[0:self.Ncharge]

    @Ctmlist.setter
    def Ctmlist(self, v):
        if self.Ncharge != len(v):
            _Logger.warning('List length does not match Ncharge.')
        ls = np.zeros(100)
        ls[0:len(v)] = v
        im.messaging.qmcclistin = ls

    @property
    def Flagsc(self):
        """bool: Flag of space-charge calculation."""
        return bool(im.messaging.bcurr)

    @Flagsc.setter
    def Flagsc(self, v):
        im.messaging.bcurr = float(v)

    @property
    def Energy(self):
        """float: Initial kinetic energy (unit is [eV])."""
        return im.messaging.bkenergy

    @Energy.setter
    def Energy(self, v):
        im.messaging.bkenergy = float(v)
        self._flagengy = True

    @property
    def Mass(self):
        """float: Mass of the reference particle (unit is [eV/c^2])."""
        return im.messaging.bmass

    @Mass.setter
    def Mass(self, v):
        im.messaging.bmass = float(v)
        self._flagmass = True

    @property
    def Chargestate(self):
        """float: Charge state of the reference particle (unit is [1])."""
        return im.messaging.bchargein

    @Chargestate.setter
    def Chargestate(self, v):
        im.messaging.bchargein = float(v)

    @property
    def Frequency(self):
        """float: Scaling frequency (unit is [Hz])."""
        return im.messaging.bfreq

    @Frequency.setter
    def Frequency(self, v):
        im.messaging.bfreq = float(v)
        im.physconstclass.construct_physconst(im.messaging.bfreq)
        self._flagfreq = True

    @property
    def Phaseini(self):
        """float: Initial phase of the reference particle (unit is [rad])."""
        return im.messaging.phsini

    @Phaseini.setter
    def Phaseini(self, v):
        im.messaging.phsini = float(v)

    @property
    def Outputemit(self):
        """float: Percentage for emittance and radius output (unit is [%]).

        Notes
        -----
        This value is used in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.
        """
        return im.messaging.outputflag

    @Outputemit.setter
    def Outputemit(self, v):
        im.messaging.outputflag = float(v)

    # -------------------------------- #
    # Initial Distribution Parameters  #
    # -------------------------------- #

    @property
    def Disttype(self):
        """int: Initial distribution type.

        Notes
        -----
        - 1 : 6D rectangle in phase space for single charge state.

        - 2 : 6D Gaussian distribution for single charge state.

        - 3 : 6D Waterbag distribution for single charge state.

        - 4 : Semi-Gaussian distribution for single charge state.

        - 5 : KV distribution for single charge state.

        - 16 : 6D Waterbag distribution for multi charge states (recommended).

        - 17 : 6D Gaussian distribution for multi charge states (recommended).

        - 19 : User input distribution (recommended). File path is defined in :py:func:`distribute <impact.control.Sequence.distribute>`.
        """
        return im.messaging.flagdist

    @Disttype.setter
    def Disttype(self, v):
        if isinstance(v, basestring):
            v = const.ddist[v.lower()]

        im.messaging.flagdist = int(v)

    @property
    def Distxsigma(self):
        """float: x standard deviation of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[0]

    @Distxsigma.setter
    def Distxsigma(self, v):
        im.messaging.distparam[0] = float(v)

    @property
    def Distxlambda(self):
        """float: x' (xp, px) standard deviation of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[1]

    @Distxlambda.setter
    def Distxlambda(self, v):
        im.messaging.distparam[1] = float(v)

    @property
    def Distxmu(self):
        """float: x-x' coupling term of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[2]

    @Distxmu.setter
    def Distxmu(self, v):
        im.messaging.distparam[2] = float(v)

    @property
    def Distysigma(self):
        """float: y standard deviation of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[7]

    @Distysigma.setter
    def Distysigma(self, v):
        im.messaging.distparam[7] = float(v)

    @property
    def Distylambda(self):
        """float: y' (yp, py) standard deviation of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[8]

    @Distylambda.setter
    def Distylambda(self, v):
        im.messaging.distparam[8] = float(v)

    @property
    def Distymu(self):
        """float: y-y' coupling term of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[9]

    @Distymu.setter
    def Distymu(self, v):
        im.messaging.distparam[9] = float(v)

    @property
    def Distzsigma(self):
        """float: z standard deviation of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[14]

    @Distzsigma.setter
    def Distzsigma(self, v):
        im.messaging.distparam[14] = float(v)

    @property
    def Distzlambda(self):
        """float: z' (zp, pz) standard deviation of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[15]

    @Distzlambda.setter
    def Distzlambda(self, v):
        im.messaging.distparam[15] = float(v)

    @property
    def Distzmu(self):
        """float: z-z' coupling term of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[16]

    @Distzmu.setter
    def Distzmu(self, v):
        im.messaging.distparam[16] = float(v)

    @property
    def Distxmismatch(self):
        """float: x mismatch factor of the initial distribution (unit is [1])."""
        return im.messaging.distparam[3]

    @Distxmismatch.setter
    def Distxmismatch(self, v):
        im.messaging.distparam[3] = float(v)

    @property
    def Distxpmismatch(self):
        """float: x' (xp, px) mismatch factor of the initial distribution (unit is [1])."""
        return im.messaging.distparam[4]

    @Distxpmismatch.setter
    def Distxpmismatch(self, v):
        im.messaging.distparam[4] = float(v)

    @property
    def Distxoffset(self):
        """float: x position offset of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[5]

    @Distxoffset.setter
    def Distxoffset(self, v):
        im.messaging.distparam[5] = float(v)

    @property
    def Distxpoffset(self):
        """float: x' (xp, px) momentum offset of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[6]

    @Distxpoffset.setter
    def Distxpoffset(self, v):
        im.messaging.distparam[6] = float(v)

    @property
    def Distymismatch(self):
        """float: y mismatch factor of the initial distribution (unit is [1])."""
        return im.messaging.distparam[10]

    @Distymismatch.setter
    def Distymismatch(self, v):
        im.messaging.distparam[10] = float(v)

    @property
    def Distypmismatch(self):
        """float: x' (xp, px) mismatch factor of the initial distribution (unit is [1])."""
        return im.messaging.distparam[11]

    @Distypmismatch.setter
    def Distypmismatch(self, v):
        im.messaging.distparam[11] = float(v)

    @property
    def Distyoffset(self):
        """float: y position offset of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[12]

    @Distyoffset.setter
    def Distyoffset(self, v):
        im.messaging.distparam[12] = float(v)

    @property
    def Distypoffset(self):
        """float: y' (yp, py) momentum offset of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[13]

    @Distypoffset.setter
    def Distypoffset(self, v):
        im.messaging.distparam[13] = float(v)

    @property
    def Distzmismatch(self):
        """float: z mismatch factor of the initial distribution (unit is [1])."""
        return im.messaging.distparam[17]

    @Distzmismatch.setter
    def Distzmismatch(self, v):
        im.messaging.distparam[17] = float(v)

    @property
    def Distzpmismatch(self):
        """float: z' (zp, pz) mismatch factor of the initial distribution (unit is [1])."""
        return im.messaging.distparam[18]

    @Distzpmismatch.setter
    def Distzpmismatch(self, v):
        im.messaging.distparam[18] = float(v)

    @property
    def Distzoffset(self):
        """float: z position offset of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[19]

    @Distzoffset.setter
    def Distzoffset(self, v):
        im.messaging.distparam[19] = float(v)

    @property
    def Distzpoffset(self):
        """float: z' (zp, pz) offset of the initial distribution (IMPACT internal unit)."""
        return im.messaging.distparam[20]

    @Distzpoffset.setter
    def Distzpoffset(self, v):
        im.messaging.distparam[20] = float(v)

    # --------------------------------------------- #
    # Special Input Format for Initial Distribution #
    # --------------------------------------------- #

    @property
    def Distxtwiss(self):
        """list of float: x twiss parameter for the initial distribution.

        List of input parameters and units is [alpha(1), beta(m), normalized emittance(m-rad)].

        Notes
        -----
        Input list length must be 3.

        keywords: Courant Snyder parameters, horizontal
        """
        return self._xtwiss

    @Distxtwiss.setter
    def Distxtwiss(self, v):
        if len(v) != 3:
            raise TypeError('List length must be 3 (Twiss parameter [alpha(1), beta(m), normalized emittance(m-rad)]).')

        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        alpha, beta, emit = v
        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)
        self.Distxsigma = np.sqrt(emit/b0g0*beta/(1.0 + alpha*alpha))/im.physconstclass.scxl
        self.Distxlambda = b0g0*np.sqrt(emit/b0g0/beta)
        self.Distxmu = alpha/np.sqrt(1.0 + alpha*alpha)
        self._xtwiss[:] = alpha, beta, emit

    @property
    def Distxquadratic(self):
        """list of float: x quadratic parameter (IMPACT default) for the initial distribution.

        List of input parameters is [sigma, lambda, mu], IMPACT internal units.

        Notes
        -----
        Input list length must be 3.

        keywords: horizontal
        """
        return np.asarray(im.messaging.distparam[0:3])

    @Distxquadratic.setter
    def Distxquadratic(self, v):
        if len(v) != 3:
            raise TypeError('List length must be 3 (Quadratic parameter [sigma, lambda, mu]).')

        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        sigma, lambda_, mu = v
        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)
        alpha = np.sqrt(mu*mu/(1-mu*mu))
        fac = im.physconstclass.scxl*sigma/np.sqrt(1.0 - mu*mu)
        if fac == 0.0:
            beta = 0.0
            emit = 0.0
        else:
            beta = fac*b0g0/lambda_
            emit = fac*lambda_
        im.messaging.distparam[0:3] = sigma, lambda_, mu
        self._xtwiss[:] = alpha, beta, emit

    @property
    def Distytwiss(self):
        """list of float: y twiss parameter for the initial distribution.

        List of input parameters and units is [alpha(1), beta(m), normalized emittance(m-rad)].

        Notes
        -----
        Input list length must be 3.

        keywords: Courant Snyder parameters, vertical
        """
        return self._ytwiss

    @Distytwiss.setter
    def Distytwiss(self, v):
        if len(v) != 3:
            raise TypeError('List length must be 3 (Twiss parameter [alpha(1), beta(m), normalized emittance(m-rad)]).')

        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        alpha, beta, emit = v
        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)
        self.Distysigma = np.sqrt(emit/b0g0*beta/(1.0 + alpha*alpha))/im.physconstclass.scxl
        self.Distylambda = b0g0*np.sqrt(emit/b0g0/beta)
        self.Distymu = alpha/np.sqrt(1.0 + alpha*alpha)
        self._ytwiss[:] = alpha, beta, emit

    @property
    def Distyquadratic(self):
        """list of float: y quadratic parameter (IMPACT default) for the initial distribution.

        List of input parameters is [sigma, lambda, mu], IMPACT internal units.

        Notes
        -----
        Input list length must be 3.

        keywords: vertical
        """
        return np.asarray(im.messaging.distparam[7:10])

    @Distyquadratic.setter
    def Distyquadratic(self, v):
        if len(v) != 3:
            raise TypeError('List length must be 3 (Quadratic parameter [sigma, lambda, mu]).')

        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        sigma, lambda_, mu = v
        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)
        alpha = np.sqrt(mu*mu/(1.0 - mu*mu))
        fac = im.physconstclass.scxl*sigma/np.sqrt(1.0 - mu*mu)
        if fac == 0.0:
            beta = 0.0
            emit = 0.0
        else:
            beta = fac*b0g0/lambda_
            emit = fac*lambda_
        im.messaging.distparam[7:10] = sigma, lambda_, mu
        self._ytwiss[:] = alpha, beta, emit

    @property
    def Distztwiss(self):
        """list of float: z twiss parameter for the initial distribution.

        List of input parameters and units is [alpha(1), beta(m), normalized emittance(m-rad)].

        Notes
        -----
        Input list length must be 3.

        keywords: Courant Snyder parameters, longitudinal
        """
        return self._ztwiss

    @Distztwiss.setter
    def Distztwiss(self, v):
        if len(v) != 3:
            raise TypeError('List length must be 3 (Twiss parameter [alpha(1), beta(m), normalized emittance(m-rad)]).')

        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        alpha, beta, emit = v
        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)

        m2 = 1.0/self.Energy*g0/(1.0+g0)
        d2 = 1.0/360.0/self.Frequency*const.clight*b0g0/g0

        self.Distzsigma = const.udegree['rad']*np.sqrt(emit*beta/b0g0/d2/d2/(1.0 + alpha*alpha))
        self.Distzlambda = np.sqrt(emit/beta/b0g0/m2/m2)/self.Mass
        self.Distzmu = alpha/np.sqrt(1.0 + alpha*alpha)
        self._ztwiss[:] = alpha, beta, emit

    @property
    def Distzquadratic(self):
        """list of float: z quadratic parameter (IMPACT default) for the initial distribution.

        List of input parameters is [sigma, lambda, mu], IMPACT internal units.

        Notes
        -----
        Input list length must be 3.

        keywords: vertical
        """
        return np.asarray(im.messaging.distparam[14:17])

    @Distzquadratic.setter
    def Distzquadratic(self, v):
        if len(v) != 3:
            raise TypeError('List length must be 3 (Quadratic parameter [sigma, lambda, mu]).')

        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        sigma, lambda_, mu = v
        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)

        m2 = 1.0/self.Energy*g0/(1.0+g0)
        d2 = 1.0/360.0/self.Frequency*const.clight*b0g0/g0

        alpha = np.sqrt(mu*mu/(1-mu*mu))
        fac = const.uradian['deg']*sigma*np.sqrt(1 - mu*mu)
        if fac == 0.0:
            beta = 0.0
            emit = 0.0
        else:
            beta = fac/self.Mass/lambda_*d2/m2
            emit = fac*self.Mass*lambda_*b0g0*d2*m2
        im.messaging.distparam[14:17] = sigma, lambda_, mu
        self._ztwiss[:] = alpha, beta, emit

    @property
    def Distxoffset_m(self):
        """float: x position offset of the initial distribution (unit is [m])."""
        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")

        return im.messaging.distparam[5]*im.physconstclass.scxl

    @Distxoffset_m.setter
    def Distxoffset_m(self, v):
        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")

        im.messaging.distparam[5] = float(v)/im.physconstclass.scxl

    @property
    def Distxpoffset_rad(self):
        """float: x' (xp, px) momentum offset of the initial distribution (unit is [rad])."""
        if im.messaging.runseq < 2:
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)

        return im.messaging.distparam[6]/b0g0

    @Distxpoffset_rad.setter
    def Distxpoffset_rad(self, v):
        if im.messaging.runseq < 2:
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)

        im.messaging.distparam[6] = float(v)*b0g0

    @property
    def Distyoffset_m(self):
        """float: y position offset of the initial distribution (unit is [m])."""
        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")

        return im.messaging.distparam[12]*im.physconstclass.scxl

    @Distyoffset_m.setter
    def Distyoffset_m(self, v):
        if im.messaging.runseq < 2:
            if not self._flagfreq:
                raise ValueError("Scaling frequency is not defined. Please input 'Frequency' at first.")

        im.messaging.distparam[12] = float(v)/im.physconstclass.scxl

    @property
    def Distypoffset_rad(self):
        """float: y' (yp, py) momentum offset of the initial distribution (unit is [rad])."""
        if im.messaging.runseq < 2:
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)

        return im.messaging.distparam[13]/b0g0

    @Distypoffset_rad.setter
    def Distypoffset_rad(self, v):
        if im.messaging.runseq < 2:
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")
            if not self._flagengy:
                raise ValueError("Kinetic energy is not defined. Please input 'Energy' at first.")

        g0 = self.Energy/self.Mass + 1.0
        b0g0 = np.sqrt(g0*g0 - 1.0)

        im.messaging.distparam[13] = float(v)*b0g0

    @property
    def Distzoffset_deg(self):
        """float: z position offset of the initial distribution (unit is [deg])."""
        return im.messaging.distparam[12]*const.uradian['deg']

    @Distzoffset_deg.setter
    def Distzoffset_deg(self, v):
        im.messaging.distparam[12] = float(v)/const.uradian['deg']

    @property
    def Distzpoffset_eV(self):
        """float: z' (zp, pz) offset of the initial distribution (unit is [eV])."""
        if im.messaging.runseq < 2:
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")

        return im.messaging.distparam[13]*self.Mass

    @Distzpoffset_eV.setter
    def Distzpoffset_eV(self, v):
        if im.messaging.runseq < 2:
            if not self._flagmass:
                raise ValueError("Particle mass is not defined. Please input 'Mass' at first.")

        im.messaging.distparam[13] = float(v)/self.Mass
