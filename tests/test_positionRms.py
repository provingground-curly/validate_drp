#
# LSST Data Management System
# Copyright 2012-2017 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#


import unittest

import numpy as np

from numpy.testing import assert_almost_equal, assert_array_max_ulp

import lsst.utils
from lsst.validate.drp import util


def test_empty_positionRms():
    ra = np.deg2rad(np.array([]))
    dec = np.deg2rad(np.array([]))
    ra_avg = np.mean(ra)
    dec_avg = np.mean(dec)

    obs = util.positionRms(ra_avg, dec_avg, ra, dec)

    assert np.isnan(obs)


def test_single_positionRms():
    ra = np.deg2rad(np.array([10.0010]))
    dec = np.deg2rad(np.array([20.001]))
    ra_avg = np.mean(ra)
    dec_avg = np.mean(dec)

    exp = 0
    obs = util.positionRms(ra_avg, dec_avg, ra, dec)

    assert_array_max_ulp(obs, exp)


def test_basic_positionRms():
    ra = np.deg2rad(np.array([10.0010, 10.0005, 10.0000, 10.0005]))
    dec = np.deg2rad(np.array([20.001, 20.006, 20.002, 20.004]))
    ra_avg = np.mean(ra)
    dec_avg = np.mean(dec)

    exp = 0.0019488135998463505 * 3600 * 1000  # degrees * arcsec/degree * milliarcsec/arcsec
    obs = util.positionRms(ra_avg, dec_avg, ra, dec)

    assert_almost_equal(obs, exp, decimal=7)  # 1e-7 rad == 5.73e-6 deg == 36 milliarcsec


def test_north_pole_positionRms():
    ra = np.deg2rad(np.array([10.0010, 10.0005, 190.0000, 190.0005]))
    dec = np.deg2rad(np.array([89.999, 89.998, 89.999, 89.998]))
    ra_avg = np.mean(ra)
    dec_avg = np.mean(dec)

    exp = 0.0021794464685958273 * 3600 * 1000  # degrees * arcsec/degree * milliarcsec/arcsec
    obs = util.positionRms(ra_avg, dec_avg, ra, dec)

    assert_almost_equal(obs, exp, decimal=7)  # 1e-7 rad == 5.73e-6 deg == 36 milliarcsec


def test_south_pole_positionRms():
    ra = np.deg2rad(np.array([10.0010, 10.0005, 190.0000, 190.0005]))
    dec = -np.deg2rad(np.array([89.999, 89.998, 89.999, 89.998]))
    ra_avg = np.mean(ra)
    dec_avg = np.mean(dec)

    exp = 0.0021794464685958273 * 3600 * 1000  # degrees * arcsec/degree * milliarcsec/arcsec
    obs = util.positionRms(ra_avg, dec_avg, ra, dec)

    assert_almost_equal(obs, exp, decimal=7)  # 1e-7 rad == 5.73e-6 deg == 36 milliarcsec


def test_wrap_positionRms():
    ra = np.deg2rad(np.array([10.0010, 10.0005, 370.0000, 370.0005]))
    dec = np.deg2rad(np.array([20.001, 20.006, 20.002, 20.004]))
    ra_avg = np.mean(ra % (2*np.pi))
    dec_avg = np.mean(dec)

    exp = 0.0019488135998463505 * 3600 * 1000  # degrees * arcsec/degree * milliarcsec/arcsec
    obs = util.positionRms(ra_avg, dec_avg, ra, dec)

    assert_almost_equal(obs, exp, decimal=7)  # 1e-7 rad == 5.73e-6 deg == 36 milliarcsec


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
