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

from __future__ import absolute_import, division, print_function

import filecmp
import os
import tempfile
import unittest

from lsst.validate.drp import performance_summary


class PerformanceSummaryFromJob(unittest.TestCase):
    """Testing release performance summary."""

    def setUp(self):
        test_data_dir = os.path.dirname(__file__)
        self.report_file = os.path.join(test_data_dir, 'CfhtQuick_output_r_report_file.rst')
        self.json_file = os.path.join(test_data_dir, 'CfhtQuick_output_r.json')
        self.json_file_filter = 'r'

    def test_generate_report_from_json(self):
        """Can we read in JSON results and make a report.

        Reference datasets are from a v13.0 validation_data_hsc run.
        """
        # Manually use temporary directories here,
        #  because I can't figure out how to get py.test tmpdir fixture
        #  to work in the unittest.TestCase context.
        tmp_dir = tempfile.mkdtemp()
        out_file_name = os.path.join(tmp_dir, "performance_summary_test.rst")

        performance_summary.run([self.json_file], out_file_name)

        assert(os.path.exists(out_file_name))
        assert filecmp.cmp(out_file_name, self.report_file)

        # Cleanup our temp files
        os.remove(out_file_name)
        os.removedirs(tmp_dir)
