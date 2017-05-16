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

# Set seed for repeatibility
np.random.seed(96701237)
m = np.random.randn(1000)*2 + 25
mag = m[m < 25]
mag_err = photErrModel(mag, sigmaSys, gamma, m5)

# Resample mag
# Set seed for repeatibility.  We explicityly reset it here to a known value
#  in case `photErrModel` called np.random
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
