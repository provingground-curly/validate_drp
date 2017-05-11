# LSST Data Management System
# Copyright 2017 AURA/LSST.
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

from __future__ import print_function, absolute_import, division

import numpy as np
import os

from lsst.validate.drp.photerrmodel import photErrModel, fitPhotErrModel

sigmaSys, gamma, m5 = 0.01, 0.039, 24.35  # mag, '', mag
# m = randn(100)*2 + 25
# m[m<25]
mag = np.array([20.61924649, 22.38816749, 23.42808267, 23.91382154,
                21.03983081, 22.66172609, 22.99162843, 22.15120047,
                23.81592911, 24.18329084, 24.06935126, 24.56872001,
                22.87464248, 19.65561574, 19.91919161, 22.68216205,
                23.64408218, 23.76364731, 24.37833970, 24.64955185,
                23.97967389, 24.38613033, 23.97724519, 24.13588853,
                21.66608579, 24.05617257, 24.68853882, 24.10951259,
                23.57543272, 23.01713577, 24.85388469, 22.98283602,
                22.37079191, 21.94217066, 23.75503379, 22.49101874,
                22.64036787, 22.01897656, 24.01141995, 21.90216590,
                22.40675442, 24.29205718, 24.54783765, 23.63751763,
                24.34770239, 24.66212856, 24.88785257, 24.49763894])
mag_err = photErrModel(mag, sigmaSys, gamma, m5)
# Resample mag
np.random.seed(23987)
noisy_mag = mag + np.random.randn(len(mag_err))*mag_err


def test_perfect_fit_phot_error_model():
    """Does a simple fit to a small, perfect, data set work."""
    fit_results = fitPhotErrModel(mag, mag_err)
    fit_diff = np.array([fit_results['sigmaSys'].value - sigmaSys,
                         fit_results['gamma'].value - gamma,
                         fit_results['m5'].value - m5])
    eps = 1e-7
    assert np.sum(fit_diff**2) < eps


def test_noisy_fit_phot_error_model():
    """Does a simple fit to a small, perfect, data set work."""
    fit_results = fitPhotErrModel(noisy_mag, mag_err)
    fit_diff = np.array([fit_results['sigmaSys'].value - sigmaSys,
                         fit_results['gamma'].value - gamma,
                         fit_results['m5'].value - m5])
    eps = 1e-2
    assert np.sum(fit_diff**2) < eps


def test_failed_fit_phot_error_model():
    """Does a failed fit recover and return NaN."""
    testDir = os.path.dirname(__file__)
    failing_datafile = os.path.join(testDir, 'mag_magerr_bright.dat')
    failing_mag, failing_mag_err = np.loadtxt(failing_datafile)

    fit_results = fitPhotErrModel(failing_mag, failing_mag_err)

    assert np.isnan(fit_results['sigmaSys'].value)
    assert np.isnan(fit_results['gamma'].value)
    assert np.isnan(fit_results['m5'].value)
