#!/usr/bin/env python

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

from __future__ import division, print_function, absolute_import

import argparse
import os.path

from lsst.utils import getPackageDir
from lsst.validate.drp import report_performance

description = """
Report the performance vs. SRD and specified metrics for a JSON cache of a job.

Produces results to:
OUTFILE
    RST summary table suitable for release report

Example call:
reportPerformance.py tests/CfhtQuick_output_r.json test_report_performance.rst \\
    --release_metrics etc/release_metrics.yaml --release_level FY18
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('json_files', nargs='+',
                        help='Paths to a JSON serializations of Jobs.  Separated with spaces.')
    parser.add_argument('--output_file', default='report_performance.rst',
                        help='Filepath of the output RST file.')
    # 2017-07-21: MWV
    # I'm not enforcing the value of 'srd_level' with a restricted 'choices'
    # because I want to preserve the ability for somewhat arbitrary levels,
    # even in the SRD for increased general utility of the function
    # The acceptable values of 'choices' are defined by what's available
    # in the 'srd_metrics' file.
    parser.add_argument('--srd_level', type=str, default='design',
                        help='Level of srd_metric requirement to meet: ["design", "minimum", "stretch"]')
    parser.add_argument('--release_metrics',
                        default=os.path.join(getPackageDir('validate_drp'),
                                             'etc', 'release_metrics.yaml'),
                        help='Path of YAML file with this release specifications.')
    parser.add_argument('--release_level', type=str, default='FY17',
                        help='Level of release_metric requirement to meet: ["FY17", "FY18", ...]')

    args = parser.parse_args()

    report_performance.run(args.json_files, args.output_file,
                           srd_level=args.srd_level,
                           release_metrics_file=args.release_metrics,
                           release_level=args.release_level)
