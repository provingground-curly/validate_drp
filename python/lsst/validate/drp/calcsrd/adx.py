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


class ADxMeasurement(MeasurementBase):
    """Measurement of AFx (x=1,2,3): The maximum fraction of astrometric
    distances which deviate by more than ADx milliarcsec (see AMx) (%).

    Parameters
    ----------
    metric : `lsst.validate.base.Metric`
        AD1, AD2, or AD3 `~lsst.validate.base.Metric` instance.
    matchedDataset : lsst.validate.drp.matchreduce.MatchedMultiVisitDataset
    amx : :class:`lsst.validate.drp.calcsrd.AMxMeasurement`
        And AMx measurement, providing the median astrometric scatter in
        the annulus.
    filter_name : str
        filter_name (filter name) used in this measurement (e.g., `'r'`).
    spec_name : str
        Name of a specification level to measure against (e.g., design,
        minimum, stretch).
    verbose : bool, optional
        Output additional information on the analysis steps.
    job : :class:`lsst.validate.drp.base.Job`, optional
        If provided, the measurement will register itself with the Job
        object.
    linkedBlobs : dict, optional
        A `dict` of additional blobs (subclasses of BlobBase) that
        can provide additional context to the measurement, though aren't
        direct dependencies of the computation (e.g., `matchedDataset).

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

    def __init__(self, metric, matchedDataset, amx, filter_name, spec_name,
                 job=None, linkedBlobs=None, verbose=False):
        MeasurementBase.__init__(self)

        self.metric = metric
        self.filter_name = filter_name
        self.spec_name = spec_name

        # register input parameters for serialization
        # note that matchedDataset is treated as a blob, separately
        self.register_parameter('D', datum=amx.parameters['D'])
        self.register_parameter('annulus', datum=amx.parameters['annulus'])
        self.register_parameter('magRange', datum=amx.parameters['magRange'])
        self.register_parameter('AMx', datum=amx.datum)

        self.matchedDataset = matchedDataset

        # Add external blobs so that links will be persisted with
        # the measurement
        if linkedBlobs is not None:
            for name, blob in linkedBlobs.items():
                setattr(self, name, blob)

        # AFx's value at this spec level and filter must be known to measure
        self.register_parameter(
            'AFx',
            datum=self.metric.get_spec_dependency(
                self.spec_name,
                'AF{0:d}'.format(self.metric.x.quantity),
                filter_name=self.filter_name))

        if amx.quantity is not None:
            # No more than AFx of values will deviate by more than the
            # AMx (50th) + AFx percentiles
            # To compute ADx, use measured AMx and spec for AFx.
            afxAtPercentile = np.percentile(
                amx.rmsDistMas.to(u.marcsec),
                100. - self.AFx) * u.marcsec
            self.quantity = afxAtPercentile - amx.quantity
        else:
            # FIXME previously would raise ValidateErrorNoStars
            self.quantity = None

        print('ADx', self.quantity)

        if job:
            job.register_measurement(self)
