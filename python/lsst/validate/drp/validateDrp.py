#!/usr/bin/env python

# LSST Data Management System
# Copyright 2008-2018 AURA/LSST.
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

import argparse
import os.path
import sys

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base.argumentParser import ArgumentParser
from lsst.pipe.supertask.super_task import SuperTask
from lsst.utils import getPackageDir
from lsst.validate.base import load_metrics
from lsst.validate.drp import validate, util

__all__ = ["ValidateDrpTask"]


description = """
Calculate and plot validation Key Project Metrics from the LSST SRD.
http://ls.st/LPM-17

Produces results to:
STDOUT
    Summary of key metrics
REPONAME*.png
    Plots of key metrics.  Generated in current working directory.
REPONAME*.json
    JSON serialization of each KPM.

where REPONAME is based on the repository name but with path separators
replaced with underscores.  E.g., "Cfht/output" -> "Cfht_output_"
"""

class ValidateDrpConfig(pexConfig.Config):
    """`ValidateDrpTask` configuration."""
    brightSnr = pexConfig.Field(
        doc='SNR cut-off for bright stars.',
        dtype=float,
        default=100.
    )
    metricsFilePath = pexConfig.Field(
        doc='Path of YAML file with LPM-17 SRD metric definitions.',
        dtype=str,
        default=os.path.join(getPackageDir('validate_drp'),
                             'etc', 'metrics.yaml')
    )
    outputPrefix = pexConfig.Field(
        doc='Prefix for output file names, including any directories.',
        dtype=str,
        default='validate_drp_'
    )
    metricsPackage = pexConfig.Field(
        doc='Name of the repository with YAML definitions of LPM-17 metrics.',
        dtype=str,
        default='verify_metrics'
    )
    makePlot = pexConfig.Field(
        doc='Render and save plots with matplotlib.',
        dtype=bool,
        default=True
    )
    targetSpecLevel = pexConfig.Field(
        doc='Level of SRD requirement to meet: "minimum", "design", "stretch"',
        dtype=str,
        default='design'
    )
    instrument = pexConfig.Field(
        doc='Name of the instrument.  If None will be extracted from the Butler mapper.'
        dtype=str,
        default=None
    )
    dataset_repo_url = pexConfig.Field(
        doc="""
Location of the dataset used.  If None will be set to the path of the repo.
Intended for use to document provenance in the summary.'
""",
        dtype=str,
        default=None
    )

# Need something more like coadd, that takes a list of data Ids 
# to prepare into single output
class ValidateDrpTaskRunner(pipeBase.TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return pipeBase.TaskRunner.getTargetList(parsedCmd, selectDataList=parsedCmd.selectId.dataList,
                                                 **kwargs)


class ValidateDrpTask(pipeBase.CmdLineTask):
    ConfigClass = ValidateDrpConfig
    RunnerClass = ValidateDrpTaskRunner
    _DefaultName = "validateDrp"

#    def __init__(self, 
#                 **kwargs):

    @pipeBase.timeMethod
    def run(self):
        if not os.path.isdir(args.repo):
            print("Could not find repo %r" % (args.repo,))
            sys.exit(1)

        kwargs = {}
        kwargs['verbose'] = args.verbose
        kwargs['level'] = args.level

        if not os.path.exists(args.metricsFile):
            print('Could not find metric definitions: {0}'.format(args.metricsFile))
            sys.exit(1)
        metrics = load_metrics(args.metricsFile)
        kwargs['metrics'] = metrics

        # Wrap the configs into keyword arguments to pass to exist code.
        # Consider pushing this down further into 'validate' module

        validate.run(args.repo, **kwargs)

    @classmethod
    def _makeArgumentParser(cls):
        """!Create and return an argument parser

        @param[in] cls      the class object
        @return the argument parser for this task.

        This override is used to delay making the data ref list until the dataset type is known;
        this is done in @ref parseAndRun.
        """

        parser = argparse.ArgumentParser(name=cls._DefaultName,
                                         description=description,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_id_argument(name="--id",
                               datasetType="calexp",
                               help="data IDs, e.g. --id visit=12345 ccd=1,2^0,3")

        return parser

