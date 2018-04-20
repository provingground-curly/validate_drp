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

from lsst.verify import Measurement, Datum


def measurePF1(metric, pa1, pa2_spec):
    """Measurement of PF1: fraction of samples between median RMS (PA1) and
    PA2 specification.

    Parameters
    ----------
    metric : `lsst.verify.Metric`
        A PF1 `~lsst.verify.Metric` instance.
    pa1 : `lsst.verify.Measurement`
        A PA1 measurement instance.
    pa2_spec : `lsst.verify.Spec`
        An `lsst.verify.Spec` that holds the threshold at which to measure PF1

    Returns
    -------
    measurement : `lsst.verify.Measurement`
        Measurement of PF1 and associated metadata.

    Notes
    -----
    The LSST Science Requirements Document (LPM-17) is commonly referred
    to as the SRD.  The SRD puts a limit that no more than PF1 % of difference
    will vary by more than PA2 millimag.  The design, minimum, and stretch
    goals are PF1 = (10, 20, 5) % at PA2 = (15, 15, 10) millimag following
    LPM-17 as of 2011-07-06, available at http://ls.st/LPM-17.
    """

    datums = {}
    datums['pa2_spec'] = Datum(quantity=pa2_spec.threshold, description="Threshold applied to PA2")
    # Use first random sample from original PA1 measurement
    magDiff = pa1.extras['magDiff'].quantity
    magDiffs = magDiff[0, :]

    quantity = 100 * np.mean(np.abs(magDiffs) > pa2_spec.threshold) * u.Unit('percent')
    return Measurement(metric, quantity, extras=datums)
