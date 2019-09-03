from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import impact._impactf as im
import logging
import numpy as np

from impact import const
from impact import util
from mpi4py import MPI

from collections import OrderedDict

__authors__ = "Kei Fukushima"
__copyright__ = "(c) 2017, Facility for Rare Isotope beams, Michigan State University"
__contact__ = "Kei Fukushima <fukushim@frib.msu.edu>"

logging.basicConfig()
_Logger = logging.getLogger(__name__)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

try:
    basestring
except NameError:
    basestring = str


class Result(object):
    """Base class for simulation results of Advanced IMPACT.

    - **Attributes**

    .. autosummary::
        dtag
        autobcast
        qmlabel

    .. _hist-label:

    - **Methods for history results**

    .. autosummary::
        hrefz
        hrefphi
        hrefgma
        hrefeng
        hrefbeta
        hrefr
        hxcen
        hycen
        hzcen
        hxrms
        hyrms
        hzrms
        hxpcen
        hypcen
        hzpcen
        hxprms
        hyprms
        hzprms
        hxtwsa
        hytwsa
        hztwsa
        hxtwsb
        hytwsb
        hztwsb
        hxepsn
        hyepsn
        hzepsn
        hxepsnp
        hyepsnp
        hzepsnp
        hxepsnf
        hyepsnf
        hzepsnf
        hxmax
        hymax
        hzmax
        hxpmax
        hypmax
        hzpmax
        hrrms
        hrnpr
        hrmax
        hx3rd
        hy3rd
        hz3rd
        hxp3rd
        hyp3rd
        hzp3rd
        hx4th
        hy4th
        hzp4th
        hxp4th
        hyp4th
        hz4th
        hlpmin
        hlpmax
        hptot

    - **Methods for distribution results**

    .. autosummary::
        getx
        getxp
        gety
        getyp
        getz
        getzp
        getctm
        getcpmp
        getid
        getall
    """

    _dtag = OrderedDict()
    _autobcast = True

    @property
    def dtag(self):
        """dict: Tag dictionary of lattice elements."""
        return self._dtag

    @dtag.setter
    def dtag(self, dct):
        if isinstance(dct, dict):
            self._dtag = dct
        else:
            TypeError('Type of tag must be dict.')

    @property
    def autobcast(self):
        """bool: Flag of auto MPI bcast for simulation results."""
        return self._autobcast

    @autobcast.setter
    def autobcast(self, v):
        self._autobcast = bool(v)

    @property
    def qmlabel(self):
        """list: List of charge to mass ratios for history results"""
        return im.messaging.qmlabel

    def _cnv_tagh(self, dcnv, *args, **kws):
        at = kws.get('at', 'exit')

        if len(args) > 2:
            raise TypeError('Number of arguments must less than or equal 2. ('+str(len(args))+' given)')

        name = kws.get('tag', None)
        unit = kws.get('unit', None)
        value = kws.get('value', None)
        for elem in args:
            if isinstance(elem, basestring) and (elem in self.dtag):
                name = elem
            elif isinstance(elem, basestring) and (elem in dcnv):
                unit = elem
            elif isinstance(elem, basestring):
                _Logger.warning('Undefined argument: ' + elem)
            else:
                value = elem

        if name is not None:
            if value is not None:
                raise TypeError("Input arguments are conflicted. Only one of 'tag' or 'value' to be input.")
            zid = np.asarray(self.dtag[name]).flatten()-1
            if at == 'mid':
                value = im.messaging.blnpos[zid] + 0.5*im.messaging.blnlen[zid]
            elif im.messaging.flagdiag == (1 or 2):
                if at == 'entry':
                    value = {'id': im.messaging.blnidg[zid]-im.messaging.blnidg[im.messaging.lcini-1]}
                elif at == 'exit':
                    value = {'id': im.messaging.blnidg[zid+1]-im.messaging.blnidg[im.messaging.lcini-1]}
                else:
                    raise TypeError("undefined keyword for 'at'")

            elif im.messaging.flagdiag == (3 or 4):
                if at == 'entry':
                    value = im.messaging.blnpos[zid]
                elif at == 'exit':
                    value = im.messaging.blnpos[zid] + im.messaging.blnlen[zid]
                else:
                    raise TypeError("undefined keyword for 'at'")

            else:
                raise TypeError("Undefined parameter for Flagdiag")

        return value, unit

    def hrefz(self, *args, **kws):
        """Returns z distance history of the reference particle (default unit is [m]).

        Parameters
        ----------
        *args :
            str -- Tag name or physical unit.

            float or list of float -- Substitution value of interpolation function.

        **kws :
            **tag** : str

                Tag name of the target position.

            **unit** : str

                Physical unit name of the return value.

            **value** : float or list of float

                Substitution value for the interpolation function.

            **qid** : int

                Charge state index for return value. Default is -1 (Total of the all charge states).
                Index is compatible with :py:func:`qmlabel <impact.result.Result.qmlabel>`.

            **at** : ``entry``, ``mid``, or ``exit`` (default)

                Data taking location in the target element. Used when the tag name is given.

        Returns
        -------
        ndarray
            z distance history

        Notes
        -----
        keywords: longitudinal, fort.18

        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hrefz, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hrefphi(self, *args, **kws):
        """Returns absolute phase history of the reference particle (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            absolute phase history

        Notes
        -----
        keywords: fort.18

        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hrefphi, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hrefbeta(self, *args, **kws):
        """Returns Lorentz beta history of the reference particle (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            Lorentz beta history

        Notes
        -----
        keywords: fort.18

        """
        argc = self._cnv_tagh(const.uNone, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hrefbeta, '1', const.uNone, self.autobcast, *argc, **kws)

    def hrefgma(self, *args, **kws):
        """Returns Lorentz gamma history of the reference particle (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            Lorentz gamma history

        Notes
        -----
        keywords: fort.18

        """
        argc = self._cnv_tagh(const.uNone, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hrefgma, '1', const.uNone, self.autobcast, *argc, **kws)

    def hrefeng(self, *args, **kws):
        """Returns kinetic energy history of the reference particle (default unit is [MeV]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            kinetic energy history

        Notes
        -----
        keywords: fort.18
        """
        argc = self._cnv_tagh(const.uMenergy, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hrefeng, 'MeV', const.uMenergy, self.autobcast, *argc, **kws)

    def hrefr(self, *args, **kws):
        """Returns maximum r distance history of particles from pipe center (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            maximum r distance history

        Notes
        -----
        keywords: fort.18
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hrefr, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hxcen(self, *args, **kws):
        """Returns x centroid position history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x centroid position history

        Notes
        -----
        keywords: horizontal, fort.24
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxcen, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hxrms(self, *args, **kws):
        """Returns x rms size history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x rms size history

        Notes
        -----
        keywords: root mean square, envelope, horizontal, fort.24
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxrms, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hxpcen(self, *args, **kws):
        """Returns x'(xp, px) centroid momentum history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x' centroid momentum history

        Notes
        -----
        keywords: horizontal, fort.24
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxpcen, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hxprms(self, *args, **kws):
        """Returns x'(xp, px) rms momentum history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x' rms momentum history

        Notes
        -----
        keywords: root mean square, envelope, horizontal, fort.24
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxprms, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hxtwsa(self, *args, **kws):
        """Returns x twiss parameter alpha history (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x twiss parameter alpha history

        Notes
        -----
        keywords: Courant Snyder parameters, horizontal, fort.24
        """
        argc = self._cnv_tagh(const.uNone, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxtwsa, '1', const.uNone, self.autobcast, *argc, **kws)

    def hxtwsb(self, *args, **kws):
        """Returns x twiss parameter beta history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x twiss parameter beta history

        Notes
        -----
        keywords: Courant Snyder parameters, horizontal, fort.24
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxtwsb, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hxepsn(self, *args, **kws):
        """Returns x normalized emittance history (default unit is [m-rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x normalized emittance history

        Notes
        -----
        keywords: horizontal, fort.24
        """
        argc = self._cnv_tagh(const.uemit, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxepsn, 'm-rad', const.uemit, self.autobcast, *argc, **kws)

    def hxepsnp(self, *args, **kws):
        """Returns `n` % x normalized emittance history (default: `n` = 99.9, unit is [m-rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            `n` % x normalized emittance history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.

        The percentage `n` [%] is defined by :py:func:`Outputemit <impact.input.Input.Outputemit>`.

        keywords: horizontal, fort.24
        """
        argc = self._cnv_tagh(const.uemit, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxepsnp, 'm-rad', const.uemit, self.autobcast, *argc, **kws)

    def hxepsnf(self, *args, **kws):
        """Returns full x normalized emittance history (default unit is [m-rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            full x normalized emittance history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.

        keywords: horizontal, fort.24
        """
        argc = self._cnv_tagh(const.uemit, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxepsnf, 'm-rad', const.uemit, self.autobcast, *argc, **kws)

    def hycen(self, *args, **kws):
        """Returns y centroid position history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y centroid position history

        Notes
        -----
        keywords: vertical, fort.25
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hycen, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hyrms(self, *args, **kws):
        """Returns y rms size history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y rms size history

        Notes
        -----
        keywords: root mean square, envelope, vertical, fort.25
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hyrms, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hypcen(self, *args, **kws):
        """Returns y'(yp, py) centroid momentum history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y' centroid momentum history

        Notes
        -----
        keywords: vertical, fort.25
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hypcen, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hyprms(self, *args, **kws):
        """Returns y'(yp, py) rms momentum history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y' rms momentum history

        Notes
        -----
        keywords: root mean square, envelope, vertical, fort.25
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hyprms, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hytwsa(self, *args, **kws):
        """Returns y twiss parameter alpha history (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y twiss parameter alpha history

        Notes
        -----
        keywords: Courant Snyder parameters, vertical, fort.25
        """
        argc = self._cnv_tagh(const.uNone, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hytwsa, '1', const.uNone, self.autobcast, *argc, **kws)

    def hytwsb(self, *args, **kws):
        """Returns y twiss parameter beta history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y twiss parameter beta history

        Notes
        -----
        keywords: Courant Snyder parameters, vertical, fort.25
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hytwsb, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hyepsn(self, *args, **kws):
        """Returns y normalized emittance history (default unit is [m-rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y normalized emittance history

        Notes
        -----
        keywords: vertical, fort.25
        """
        argc = self._cnv_tagh(const.uemit, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hyepsn, 'm-rad', const.uemit, self.autobcast, *argc, **kws)

    def hyepsnp(self, *args, **kws):
        """Returns `n` % y normalized emittance history (default: `n` = 99.9, unit is [m-rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            `n` % y normalized emittance history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.

        The percentage `n` [%] is defined by :py:func:`Outputemit <impact.input.Input.Outputemit>`.

        keywords: vertical, fort.25
        """
        argc = self._cnv_tagh(const.uemit, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hyepsnp, 'm-rad', const.uemit, self.autobcast, *argc, **kws)

    def hyepsnf(self, *args, **kws):
        """Returns full y normalized emittance history (default unit is [m-rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            full y normalized emittance history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.

        keywords: vertical, fort.25
        """
        argc = self._cnv_tagh(const.uemit, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hyepsnf, 'm-rad', const.uemit, self.autobcast, *argc, **kws)

    def hzcen(self, *args, **kws):
        """Returns z (phase) centroid position history (default unit is [deg]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z centroid position history

        Notes
        -----
        keywords: longitudinal, fort.26
        """
        argc = self._cnv_tagh(const.udegree, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hzcen, 'deg', const.udegree, self.autobcast, *argc, **kws)

    def hzrms(self, *args, **kws):
        """Returns z (phase) rms size history (default unit is [deg]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z rms size history

        Notes
        -----
        keywords: root mean square, envelope, longitudinal, fort.26
        """
        argc = self._cnv_tagh(const.udegree, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hzrms, 'deg', const.udegree, self.autobcast, *argc, **kws)

    def hzpcen(self, *args, **kws):
        """Returns z'(zp, pz) centroid history (default unit is [MeV]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z' centroid history

        Notes
        -----
        keywords: longitudinal, energy deviation, fort.26
        """
        argc = self._cnv_tagh(const.uMenergy, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hzpcen, 'MeV', const.uMenergy, self.autobcast, *argc, **kws)

    def hzprms(self, *args, **kws):
        """Returns z'(zp, pz) rms spread history (default unit is [MeV]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z' rms spread history

        Notes
        -----
        keywords: root mean square, envelope, longitudinal, energy deviation, fort.26
        """
        argc = self._cnv_tagh(const.uMenergy, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hzprms, 'MeV', const.uMenergy, self.autobcast, *argc, **kws)

    def hztwsa(self, *args, **kws):
        """Returns z twiss parameter alpha history (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z twiss parameter alpha history

        Notes
        -----
        keywords: Courant Snyder parameters, longitudinal, fort.26
        """
        argc = self._cnv_tagh(const.uNone, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hztwsa, '1', const.uNone, self.autobcast, *argc, **kws)

    def hztwsb(self, *args, **kws):
        """Returns z twiss parameter beta history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z twiss parameter beta history

        Notes
        -----
        keywords: Courant Snyder parameters, longitudinal, fort.26
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hztwsb, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hzepsn(self, *args, **kws):
        """Returns z normalized emittance history (default unit is [m-rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z normalized emittance history

        Notes
        -----
        keywords: longitudinal, fort.26
        """
        argc = self._cnv_tagh(const.uemit, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hzepsn, 'm-rad', const.uemit, self.autobcast, *argc, **kws)

    def hzepsnp(self, *args, **kws):
        """Returns `n` % z normalized emittance history (default: `n` = 99.9, unit is [m-rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            `n` % z normalized emittance history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.

        The percentage `n` [%] is defined by :py:func:`Outputemit <impact.input.Input.Outputemit>`.

        keywords: longitudinal, fort.26
        """
        argc = self._cnv_tagh(const.uemit, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hzepsnp, 'm-rad', const.uemit, self.autobcast, *argc, **kws)

    def hzepsnf(self, *args, **kws):
        """Returns full z normalized emittance history (default unit is [m-rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            full z normalized emittance history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.

        keywords: longitudinal, fort.26
        """
        argc = self._cnv_tagh(const.uemit, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hzepsnf, 'm-rad', const.uemit, self.autobcast, *argc, **kws)

    def hxmax(self, *args, **kws):
        """Returns maximum x position history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            maximum x position history

        Notes
        -----
        keywords: horizontal, fort.27
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxmax, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hxpmax(self, *args, **kws):
        """Returns maximum x' (xp, px) momentum history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            maximum x' momentum history

        Notes
        -----
        keywords: horizontal, fort.27
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxpmax, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hymax(self, *args, **kws):
        """Returns maximum y position history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            maximum y position history

        Notes
        -----
        keywords: vertical, fort.27
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hymax, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hypmax(self, *args, **kws):
        """Returns maximum y' (yp, py) momentum history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            maximum y' momentum history

        Notes
        -----
        keywords: vertical, fort.27
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hypmax, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hzmax(self, *args, **kws):
        """Returns maximum z position (phase) history (default unit is [deg]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            maximum z position history

        Notes
        -----
        ``hphimax`` is alias of ``hzmax`` (``hphimax`` is named in fortran code).

        keywords: longitudinal, fort.27
        """
        argc = self._cnv_tagh(const.udegree, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hphimax, 'deg', const.udegree, self.autobcast, *argc, **kws)

    # add alias for hzmax
    hphimax = hzmax

    def hzpmax(self, *args, **kws):
        """Returns maximum z' (zp, pz) history (default unit is [MeV]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            maximum z' history

        Notes
        -----
        ``hdemax`` is alias of ``hzpmax`` (``hdemax`` is named in fortran code).

        keywords: longitudinal, energy deviation, fort.27
        """
        argc = self._cnv_tagh(const.uMenergy, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hdemax, 'MeV', const.uMenergy, self.autobcast, *argc, **kws)

    # add alias for hzpmax
    hdemax = hzpmax

    def hlpmin(self, *args, **kws):
        """Returns minimum local particle number history (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            minimum local particle number history

        Notes
        -----
        keywords: fort.28
        """
        argc = self._cnv_tagh(const.uNone, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hlpmin, '1', const.uNone, self.autobcast, *argc, **kws)

    def hlpmax(self, *args, **kws):
        """Returns maximum local particle number history (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            maximum local particle number history

        Notes
        -----
        keywords: fort.28
        """
        argc = self._cnv_tagh(const.uNone, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hlpmax, '1', const.uNone, self.autobcast, *argc, **kws)

    def hptot(self, *args, **kws):
        """Returns total particle number history (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            total particle number history

        Notes
        -----
        keywords: fort.28
        """
        argc = self._cnv_tagh(const.uNone, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hptot, '1', const.uNone, self.autobcast, *argc, **kws)

    def hx3rd(self, *args, **kws):
        """Returns x 3rd root of 3rd moment history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x 3rd moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: horizontal, fort.29
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hx3rd, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hxp3rd(self, *args, **kws):
        """Returns x' (xp, px) 3rd root of 3rd moment history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x' 3rd moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: horizontal, fort.29
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxp3rd, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hy3rd(self, *args, **kws):
        """Returns y 3rd root of 3rd moment history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y 3rd moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: vertical, fort.29
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hy3rd, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hyp3rd(self, *args, **kws):
        """Returns y' (yp, py) 3rd root of 3rd moment history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y' 3rd moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: vertical, fort.29
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hyp3rd, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hz3rd(self, *args, **kws):
        """Returns z (phase) 3rd root of 3rd moment history (default unit is [deg]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z 3rd moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: longitudinal, fort.29
        """
        argc = self._cnv_tagh(const.udegree, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hz3rd, 'deg', const.udegree, self.autobcast, *argc, **kws)

    def hzp3rd(self, *args, **kws):
        """Returns z' (zp, pz) 3rd root of 3rd moment history (default unit is [MeV]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z' 3rd moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: longitudinal, energy deviation, fort.29
        """
        argc = self._cnv_tagh(const.uMenergy, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hzp3rd, 'MeV', const.uMenergy, self.autobcast, *argc, **kws)

    def hx4th(self, *args, **kws):
        """Returns x 4th root of 4th moment history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x 4th moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: horizontal, fort.30
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hx4th, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hxp4th(self, *args, **kws):
        """Returns x' (xp, px) 4th root of 4th moment history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            x' 4th moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: horizontal, fort.30
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hxp4th, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hy4th(self, *args, **kws):
        """Returns y 4th root of 4th moment history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y 4th moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: vertical, fort.30
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hy4th, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hyp4th(self, *args, **kws):
        """Returns y' (yp, py) 4th root of 4th moment history (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            y' 4th moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: vertical, fort.30
        """
        argc = self._cnv_tagh(const.uradian, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hyp4th, 'rad', const.uradian, self.autobcast, *argc, **kws)

    def hz4th(self, *args, **kws):
        """Returns z (phase) 4th root of 4th moment history (default unit is [deg]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z 4th moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: longitudinal, fort.30
        """
        argc = self._cnv_tagh(const.udegree, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hz4th, 'deg', const.udegree, self.autobcast, *argc, **kws)

    def hzp4th(self, *args, **kws):
        """Returns z' (zp, pz) 4th root of 4th moment history (default unit is [MeV]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            z' 4th moment history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 1 or 3.

        keywords: longitudinal, energy deviation, fort.30
        """
        argc = self._cnv_tagh(const.uMenergy, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hzp4th, 'MeV', const.uMenergy, self.autobcast, *argc, **kws)

    def hrrms(self, *args, **kws):
        """Returns r rms history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            r rms size history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.

        keywords: root mean square, envelope, radius, fort.29
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hrrms, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hrnpr(self, *args, **kws):
        """Returns `n` % r size history (default: `n` = 99.9) (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            `n` % r size history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.

        The percentage `n` [%] is defined by :py:func:`Outputemit <impact.input.Input.Outputemit>`.

        keywords: radius, fort.29
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hrnpr, 'm', const.ulength, self.autobcast, *argc, **kws)

    def hrmax(self, *args, **kws):
        """Returns maximum r distance history (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        **kws :
            The same keywords are available with :py:func:`hrefz <impact.result.Result.hrefz>`.

        Returns
        -------
        ndarray
            maximum r distance history

        Notes
        -----
        This value is calculated in case of :py:func:`Flagdiag <impact.input.Input.Flagdiag>` equals to 2 or 4.

        keywords: radius, fort.29
        """
        argc = self._cnv_tagh(const.ulength, *args, **kws)
        return util._usup2d(im.messaging.hrefz, im.messaging.hrmax, 'm', const.ulength, self.autobcast, *argc, **kws)

    def _cnv_tagd(self, dcnv, *args, **kws):
        if len(args) > 2:
            raise TypeError('Number of arguments must less than or equal 2. ('+str(len(args))+' given)')

        name = kws.get('tag', None)
        pidx = kws.get('pid', None)
        for elem in args:
            if isinstance(elem, basestring) and (elem in self.dtag):
                name = elem
            elif isinstance(elem, basestring) and (elem in dcnv):
                pass
            elif isinstance(elem, basestring):
                _Logger.warning('Undefined argument: ' + elem)
            elif isinstance(elem, (int, np.integer)):
                pidx = elem
            else:
                pidx = 0

        if name is not None:
            if pidx is not None:
                raise TypeError("Input arguments are conflicted. Only one of 'tag' or 'pid' to be input.")

            if not isinstance(self.dtag[name], (int, np.integer)):
                _Logger.warning('Multiple target found in tag: ', name, '. First element was chosen.')
                for i in self.dtag[name]:
                    lng, seg, stp, typ, vrr = im.pyfunction.configure0(i)
                    if typ == -2:
                        break

                if typ != -2:
                    raise ValueError('no distribution storage in tag : ', name)

                pidx = stp

            else:
                idx = self.dtag[name]
                lng, seg, stp, typ, vrr = im.pyfunction.configure0(idx)
                pidx = stp
        elif pidx is None:
            pidx = 0

        return pidx

    def getx(self, *args, **kws):
        """Returns x positions of the particles (default unit is [m]).

        Parameters
        ----------
        *args :
            str -- Tag name or physical unit.

            int -- ID number in the flag element. Default is 0 ( = current distribution).

        **kws :
            **tag** : str

                Tag name of the target position.

            **unit** : str

                Physical unit name of the return value.

            **pid** : int

                ID number in the flag element.

        Returns
        -------
        ndarray
            x positions of the particles

        Notes
        -----
        keywords: horizontal, distribution
        """
        pidx = self._cnv_tagd(const.ulength, *args, **kws)
        return _dsup(pidx, 0, im.messaging.ninedata_scxl, 'm', const.ulength, self.autobcast, *args, **kws)

    def getxp(self, *args, **kws):
        """Returns x' (xp, px) momentums of the particles (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`getx <impact.result.Result.getx>`.

        **kws :
            The same keywords are available with :py:func:`getx <impact.result.Result.getx>`.

        Returns
        -------
        ndarray
            x' momentums of the particles

        Notes
        -----
        keywords: horizontal, distribution
        """
        pidx = self._cnv_tagd(const.uradian, *args, **kws)
        return _dsup(pidx, 1, 1.0/im.messaging.ninedata_bg, 'rad', const.uradian, self.autobcast, *args, **kws)

    def gety(self, *args, **kws):
        """Returns y positions of the particles (default unit is [m]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`getx <impact.result.Result.getx>`.

        **kws :
            The same keywords are available with :py:func:`getx <impact.result.Result.getx>`.

        Returns
        -------
        ndarray
            y positions of the particles

        Notes
        -----
        keywords: vertical, distribution
        """
        pidx = self._cnv_tagd(const.ulength, *args, **kws)
        return _dsup(pidx, 2, im.messaging.ninedata_scxl, 'm', const.ulength, self.autobcast, *args, **kws)

    def getyp(self, *args, **kws):
        """Returns y' (yp, py) momentums of the particles (default unit is [rad]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`getx <impact.result.Result.getx>`.

        **kws :
            The same keywords are available with :py:func:`getx <impact.result.Result.getx>`.

        Returns
        -------
        ndarray
            y' momentums of the particles

        Notes
        -----
        keywords: vertical, distribution
        """
        pidx = self._cnv_tagd(const.uradian, *args, **kws)
        return _dsup(pidx, 3, 1.0/im.messaging.ninedata_bg, 'rad', const.uradian, self.autobcast, *args, **kws)

    def getz(self, *args, **kws):
        """Returns z positions (phase) of the particles (default unit is [deg]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`getx <impact.result.Result.getx>`.

        **kws :
            The same keywords are available with :py:func:`getx <impact.result.Result.getx>`.

        Returns
        -------
        ndarray
            z positions of the particles

        Notes
        -----
        keywords: longitudinal, distribution
        """
        pidx = self._cnv_tagd(const.udegree, *args, **kws)
        return _dsup(pidx, 4, 180.0/np.pi, 'deg', const.udegree, self.autobcast, *args, **kws)

    def getzp(self, *args, **kws):
        """Returns z' (zp, pz) of the particles (default unit is [MeV]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`getx <impact.result.Result.getx>`.

        **kws :
            The same keywords are available with :py:func:`getx <impact.result.Result.getx>`.

        Returns
        -------
        ndarray
            z' momentums of the particles

        Notes
        -----
        keywords: energy deviation, longitudinal, distribution
        """
        pidx = self._cnv_tagd(const.uenergy, *args, **kws)
        return _dsup(pidx, 5, im.messaging.ninedata_mass, 'MeV', const.uenergy, self.autobcast, *args, **kws)

    def getctm(self, *args, **kws):
        """Returns charge to mass ratio of the particles (default unit is [c^2/eV]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`getx <impact.result.Result.getx>`.

        **kws :
            The same keywords are available with :py:func:`getx <impact.result.Result.getx>`.

        Returns
        -------
        ndarray
            charge to mass ratio of the particles
        """
        uctm = {'1': im.messaging.ninedata_mass, 'c^2/eV': 1.0}
        pidx = self._cnv_tagd(uctm, *args, **kws)
        return _dsup(pidx, 6, 1.0, 'c^2/eV', uctm, self.autobcast, *args, **kws)

    def getcpmp(self, *args, **kws):
        """Returns charge per macro particle of the particles (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`getx <impact.result.Result.getx>`.

        **kws :
            The same keywords are available with :py:func:`getx <impact.result.Result.getx>`.

        Returns
        -------
        ndarray
            charge per macro-particle of the particles
        """
        pidx = self._cnv_tagd(const.uNone, *args, **kws)
        return _dsup(pidx, 7, 1.0, '1', const.uNone, self.autobcast, *args, **kws)

    def getid(self, *args, **kws):
        """Returns ID number of the particles (default unit is [1]).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`getx <impact.result.Result.getx>`.

        **kws :
            The same keywords are available with :py:func:`getx <impact.result.Result.getx>`.

        Returns
        -------
        ndarray
            ID number of the particles
        """
        pidx = self._cnv_tagd(const.uNone, *args, **kws)
        return _dsup(pidx, 8, 1.0, '1', const.uNone, self.autobcast, *args, **kws)

    def getall(self, *args, **kws):
        """Returns full information of the particles).

        Parameters
        ----------
        *args :
            The same arguments are available with :py:func:`getx <impact.result.Result.getx>`.

        **kws :
            The same keywords are available with :py:func:`getx <impact.result.Result.getx>`.

        Returns
        -------
        (9, n) shape ndarray
            All particle information. *n* is the number of particles.

            * row 1 : x positions [m]

            * row 2 : xp momentums [rad]

            * row 3 : y positions [m]

            * row 4 : yp momentums [rad]

            * row 5 : z positions [deg]

            * row 6 : zp momentums [MeV]

            * row 7 : charge to mass ratios [c^2/eV]

            * row 8 : charge per macro-particle weights [1]

            * row 9 : ID numbers [1]
        """

        unit = kws.get('unit', None)
        if not (unit is None or unit.lower() == 'impact'):
            KeyError("input unit must be None or 'impact'")

        ret = np.array([self.getx(*args, **kws),
                        self.getxp(*args, **kws),
                        self.gety(*args, **kws),
                        self.getyp(*args, **kws),
                        self.getz(*args, **kws),
                        self.getzp(*args, **kws),
                        self.getctm(*args, **kws),
                        self.getcpmp(*args, **kws),
                        self.getid(*args, **kws)])
        return ret


def _dsup(idx, col, cnv0, unit0, dcnv, abc, *args, **kws):
    unit = unit0 if kws.get('unit', unit0) is None else kws.get('unit', unit0)
    order = kws.get('order', False)
    ierr = im.pyfunction.get_distdata(idx, order, 0)

    if ierr:
        raise ValueError('Input ID number does not found.')

    for elem in args:
        if isinstance(elem, basestring) and (elem in dcnv):
            unit = elem

    cnv1 = 1.0 if unit == 'impact' else cnv0*dcnv[unit]

    if rank == 0:
        ret = im.messaging.ninedata[col, :]*cnv1
    else:
        ret = None

    if abc and size != 1:
        ret = comm.bcast(ret, root=0)

    return np.copy(ret)
