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

import os
from lsst.validate.drp import performance_summary


def test_generate_report_from_json(tmpdir):
    """Can we read in JSON results and make a report.

    Reference datasets are from a v13.0 validation_data_hsc run.
    """

    reference_files = ['CfhtQuick_output_r.json']
    testDir = os.path.dirname(__file__)
    reference_filepaths = [os.path.join(testDir, f) for f in reference_files]

    outfile = tmpdir.mkdir("write").join(tmpdir, "metrics_report.rst")
    performance_summary.run(reference_filepaths, outfile)

    assert(os.path.exists(outfile.strpath))
