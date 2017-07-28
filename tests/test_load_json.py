#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2012-2016 LSST Corporation.
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

from __future__ import print_function

import os
import unittest

import lsst.utils

from lsst.validate.base import load_metrics
from lsst.validate.drp.validate import (
    get_filter_name_from_job, load_json_output, plot_metrics, print_metrics)
from lsst.validate.drp.util import repoNameToPrefix


class ParseJsonJob(unittest.TestCase):
    """Testing loading of JSON cache files."""

    def setUp(self):
        testDataDir = os.path.dirname(__file__)
        self.metricsFile = os.path.join(testDataDir, 'metrics.yaml')
        self.jsonFile = os.path.join(testDataDir, 'CfhtQuick_output_r.json')
        self.jsonFile_filter = 'r'

    def testLoadJsonJob(self):
        """Can we load a Job from a JSON file?"""

        # without throwing an error
        job = load_json_output(self.jsonFile)
        # Spot-check a few attributes
        self.assertEqual(len(job._measurements), 28)
        self.assertEqual(set(job.spec_levels), set(['design', 'minimum', 'stretch']))

    def testParseJobFilterName(self):
        """Do we correctly read the filterName from a Job object?"""
        job = load_json_output(self.jsonFile)
        filterName = get_filter_name_from_job(job)
        self.assertEqual(filterName, self.jsonFile_filter)

    def testPrintMetricsFromJsonJobFromMetricsfile(self):
        """Does printing the metrics run without error using metrics.yaml?

        This is in essence half the big test.
        If we load a file do we get printed metrics?
        Next TODO is to actually check for the values printed.
        """
        job = load_json_output(self.jsonFile)
        filterName = get_filter_name_from_job(job)
        metrics = load_metrics(self.metricsFile)
        print_metrics(job, filterName, metrics)

    def testPrintMetricsFromJsonJob(self):
        """Does printing the metrics run without error using itself?

        This is in essence half the big test.
        If we load a file do we get printed metrics?
        Next TODO is to actually check for the values printed.
        """
        job = load_json_output(self.jsonFile)
        filterName = get_filter_name_from_job(job)
        metrics = {meas.metric.name: meas.metric for meas in job.measurements}
        print_metrics(job, filterName, metrics)

    def testPlotMetricsFromJsonJob(self):
        """Does plotting the metrics run and produce the correct filenames?

        This could be loosely seen as a test of the outputPrefix handling
        and thus DM-11410, it's a rather incomplete one because the chain
        of wrapping functions is different.
        """
        job = load_json_output(self.jsonFile)
        filterName = get_filter_name_from_job(job)

        plot_metrics(job, filterName)
        assert os.path.exists('check_astrometry.png')

        outputPrefix = repoNameToPrefix(self.jsonFile)
        plot_metrics(job, filterName, outputPrefix=outputPrefix)
        expectedCheckAstrometryFile = outputPrefix+'check_astrometry.png'
        assert os.path.exists(expectedCheckAstrometryFile)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
