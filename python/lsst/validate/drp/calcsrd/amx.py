# LSST Data Management System
# Copyright 2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.

from __future__ import print_function, absolute_import

import numpy as np
import astropy.units as u

from lsst.validate.base import MeasurementBase
from ..util import radiansToMilliarcsec, calcRmsDistances


class AMxMeasurement(MeasurementBase):
    """Measurement of AMx (x=1,2,3): The maximum rms of the astrometric
    distance distribution for stellar pairs with separations of D arcmin
    (repeatability).

    Parameters
    ----------
    metric : `lsst.validate.base.Metric`
        An AM1, AM2 or AM3 `~lsst.validate.base.Metric` instance.
    matchedDataset : lsst.validate.drp.matchreduce.MatchedMultiVisitDataset
    filter_name : `str`
        filter_name (filter name) used in this measurement (e.g., ``'r'``).
    width : `float` or `astropy.units.Quantity`, optional
        Width around fiducial distance to include. [arcmin]
    magRange : 2-element `list`, `tuple`, or `numpy.ndarray`, optional
        brighter, fainter limits of the magnitude range to include.
        Default: ``[17.5, 21.5]`` mag.
    verbose : `bool`, optional
        Output additional information on the analysis steps.
    job : :class:`lsst.validate.drp.base.Job`, optional
        If provided, the measurement will register itself with the Job
        object.
    linkedBlobs : dict, optional
        A `dict` of additional blobs (subclasses of BlobBase) that
        can provide additional context to the measurement, though aren't
        direct dependencies of the computation (e.g., ``matchedDataset``).

    Attributes
    ----------
    rmsDistMas : ndarray
        RMS of distance repeatability between stellar pairs.
    blob : AMxBlob
        Blob with by-products from this measurement.

    Notes
    -----
    This table below is provided ``validate_drp``\ 's :file:`metrics.yaml`.

    LPM-17 dated 2011-07-06

    Specification:
        The rms of the astrometric distance distribution for
        stellar pairs with separation of D arcmin (repeatability)
        will not exceed AMx milliarcsec (median distribution for a large number
        of sources). No more than AFx % of the sample will deviate by more than
        ADx milliarcsec from the median. AMx, AFx, and ADx are specified for
        D=5, 20 and 200 arcmin for x= 1, 2, and 3, in the same order (Table 18).

    The three selected characteristic distances reflect the size of an
    individual sensor, a raft, and the camera. The required median astrometric
    precision is driven by the desire to achieve a proper motion accuracy of
    0.2 mas/yr and parallax accuracy of 1.0 mas over the course of the survey.
    These two requirements correspond to relative astrometric precision for a
    single image of 10 mas (per coordinate).

    ========================= ====== ======= =======
    Astrometric Repeatability          Specification
    ------------------------- ----------------------
                       Metric Design Minimum Stretch
    ========================= ====== ======= =======
            AM1 (milliarcsec)     10      20       5
            AF1 (%)               10      20       5
            AD1 (milliarcsec)     20      40      10
            AM2 (milliarcsec)     10      20       5
            AF2 (%)               10      20       5
            AD2 (milliarcsec)     20      40      10
            AM3 (milliarcsec)     15      30      10
            AF3 (%)               10      20       5
            AD3 (milliarcsec)     30      50      20
    ========================= ====== ======= =======

    Table 18: The specifications for astrometric precision.
    The three blocks of values correspond to D=5, 20 and 200 arcmin,
    and to astrometric measurements performed in the r and i bands.
    """

    metric = None

    def __init__(self, metric, matchedDataset, filter_name, width=2.,
                 magRange=None, linkedBlobs=None, verbose=False):
        MeasurementBase.__init__(self)

        self.metric = metric
        self.filter_name = filter_name

        # Register blob
        self.matchedDataset = matchedDataset

        # Measurement Parameters
        self.register_parameter('D', datum=self.metric.D)

        if not isinstance(width, u.Quantity):
            width = width * u.arcmin
        self.register_parameter('width',
                                quantity=width,
                                label='Width',
                                description='Width of annulus')
        if magRange is None:
            magRange = np.array([17.0, 21.5]) * u.mag
        else:
            assert len(magRange) == 2
            if not isinstance(magRange, u.Quantity):
                magRange = np.array(magRange) * u.mag
        self.register_parameter('magRange',
                                quantity=magRange,
                                description='Stellar magnitude selection '
                                            'range.')

        annulus = self.D + (self.width/2)*np.array([-1, +1])
        self.register_parameter('annulus',
                                quantity=annulus,
                                label='annulus radii',
                                description='Inner and outer radii of '
                                            'selection annulus.')

        # Register measurement extras
        self.register_extra('rmsDistMas', label='RMS')

        # Add external blob so that links will be persisted with
        # the measurement
        if linkedBlobs is not None:
            for name, blob in linkedBlobs.items():
                setattr(self, name, blob)

        matches = matchedDataset.safeMatches
        rmsDistances = calcRmsDistances(
            matches,
            self.annulus.to(u.arcmin).value,
            magRange=self.magRange.to(u.mag).value,
            verbose=verbose)[0]

        if len(rmsDistances) == 0:
            # raise ValidateErrorNoStars(
            #     'No stars found that are %.1f--%.1f arcmin apart.' %
            #     (annulus[0], annulus[1]))
            # FIXME should we still report that this measurement was
            # attempted instead of just crashing.
            print('No stars found that are {0:.1f}--{1:.1f} apart.'.format(
                  self.annulus[0], self.annulus[1]))
            self.rmsDistMas = None
            self.quantity = 5 * u.arcmin  # FIXME
        else:
            self.rmsDistMas = np.asarray(radiansToMilliarcsec(rmsDistances)) \
                * u.milliarcsecond
            self.quantity = np.median(self.rmsDistMas)
