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

from numpy.testing import assert_allclose

import lsst.utils
from lsst.validate.drp.calcsrd.amx import matchVisitComputeDistance


def test_basic_matchVisitComputeDistance():
    visit_obj1 = [1, 2, 3, 4]
    visit_obj2 = [1, 2, 3, 4]
    ra_obj1 = np.deg2rad(np.array([10.0010, 10.0005, 10.0000, 10.0005]))
    ra_obj2 = np.deg2rad(np.array([10.0008, 10.0003, 10.0005, 10.0006]))
    dec_obj1 = np.deg2rad(np.array([20.001, 20.006, 20.002, 20.004]))
    dec_obj2 = np.deg2rad(np.array([20.010, 20.030, 20.004, 20.003]))

    exp = [0.0001571138, 0.000418891, 3.58568e-05, 1.75301e-05]
    obs = matchVisitComputeDistance(visit_obj1, ra_obj1, dec_obj1,
                                    visit_obj2, ra_obj2, dec_obj2)
    # matchVisitComputeDistance doesn't preserve order
    # Sort to compare
    exp = np.sort(exp)
    obs = np.sort(obs)

    assert_allclose(exp, obs, atol=1e-7)  # 1e-7 rad == 5.73e-6 deg == 0.00036 arcsec


def test_missing_matchVisitComputeDistance():
    visit_obj1 = [1, 2, 3, 4]
    visit_obj2 = [4, 1]
    ra_obj1 = np.deg2rad(np.array([10.0010, 10.0005, 10.0000, 10.0005]))
    ra_obj2 = np.deg2rad(np.array([10.0006, 10.0008]))
    dec_obj1 = np.deg2rad(np.array([20.001, 20.006, 20.002, 20.004]))
    dec_obj2 = np.deg2rad(np.array([20.003, 20.010]))

    exp = [1.75301e-05, 0.0001571138]
    obs = matchVisitComputeDistance(visit_obj1, ra_obj1, dec_obj1,
                                    visit_obj2, ra_obj2, dec_obj2)

    # matchVisitComputeDistance doesn't preserve order
    # Sort to compare
    exp = np.sort(exp)
    obs = np.sort(obs)
    assert_allclose(exp, obs, atol=1e-7)  # 1e-7 rad == 5.73e-6 deg == 0.00036 arcsec


def test_speed_matchVisitComputeDistance(n=5000):
    # Explicitly convert n to int to ensure consistent behavior in indexing and numpy
    n = int(n)
    visits = np.arange(2*n)

    np.random.shuffle(visits)
    visit_obj1 = visits[:n]
    np.random.shuffle(visits)
    visit_obj2 = visits[:n]

    mu, sigma = 0, 3e-4  # mean 0, ~1 arcsec width
    mean_ra = 10  # degrees
    mean_dec = 20  # degrees

    ra_obj1 = mean_ra + np.random.normal(mu, sigma, n) * np.cos(np.deg2rad(mean_dec))
    ra_obj2 = mean_ra + np.random.normal(mu, sigma, n) * np.cos(np.deg2rad(mean_dec))
    dec_obj1 = mean_dec + np.random.normal(mu, sigma, n)
    dec_obj2 = mean_dec + np.random.normal(mu, sigma, n)

    matchVisitComputeDistance(visit_obj1, ra_obj1, dec_obj1,
                              visit_obj2, ra_obj2, dec_obj2)


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
