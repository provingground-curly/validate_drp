#!/usr/bin/env python
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

from __future__ import print_function

import os
import json
import numpy as np
import astropy.units as u

from lsst.utils.tests import TestCase
from lsst.utils import getPackageDir
from lsst.validate.drp.photerrmodel import fitPhotErrModel


class BlobBaseTestCase(TestCase):
    """Test BlobBase functionality."""

    def setUp(self):
        vdrp_path = getPackageDir('validate_drp')
        job_path = os.path.join(vdrp_path, 'tests', 'data',
                                'CfhtQuick_output_r.json')
        with open(job_path) as f:
            self.ref_data = json.load(f)

    def tearDown(self):
        pass

    def test_fit(self):
        """Ensure repeatability of fitPhotErrModel."""
        for b in self.ref_data['blobs']:
            if b['name'] == 'AnalyticPhotometryModel':
                datum = b['data']['sigmaSys']
                sigmaSysExpected = np.array(datum['value']) * u.Unit(datum['units'])

                datum = b['data']['gamma']
                gammaExpected = np.array(datum['value']) * u.Unit(datum['units'])

                datum = b['data']['m5']
                m5Expected = np.array(datum['value']) * u.Unit(datum['units'])

                brightSnr = b['data']['brightSnr']['value']

            if b['name'] == 'MatchedMultiVisitDataset':
                datum = b['data']['mag']
                mag = np.array(datum['value']) * u.Unit(datum['units'])

                datum = b['data']['magerr']
                magErr = np.array(datum['value']) * u.Unit(datum['units'])

                datum = b['data']['snr']
                snr = np.array(datum['value']) * u.Unit(datum['units'])

        print(brightSnr)
        bright = np.where(snr.value > brightSnr)
        print(bright)
        mag = mag[bright]
        magErr = magErr[bright]

        fit_params = fitPhotErrModel(mag, magErr)

        self.assertFloatsAlmostEqual(
            fit_params['sigmaSys'].value,
            sigmaSysExpected.value)

        self.assertFloatsAlmostEqual(
            fit_params['gamma'].value,
            gammaExpected.value)

        self.assertFloatsAlmostEqual(
            fit_params['m5'].value,
            m5Expected.value)
