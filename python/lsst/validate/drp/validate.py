# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
"""Main driver functions for metric measurements, plotting, specification
grading, and persistence.
"""

from __future__ import print_function, absolute_import

import os
from textwrap import TextWrapper

import yaml

from lsst.utils import getPackageDir
from lsst.validate.base import Metric, Job

from .util import repoNameToPrefix
from .matchreduce import (MatchedMultiVisitDataset, AnalyticPhotometryModel,
                          AnalyticAstrometryModel)
from .calcsrd import (PA1Measurement, PA2Measurement, PF1Measurement,
                      AMxMeasurement, AFxMeasurement, ADxMeasurement)
from .plot import (plotAMx, plotPA1, plotAnalyticPhotometryModel,
                   plotAnalyticAstrometryModel)


__all__ = ['run', 'runOneFilter']


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def run(repo, dataIds, metrics, outputPrefix=None, level="design", verbose=False, **kwargs):
    """Main executable.

    Runs multiple filters, if necessary, through repeated calls to `runOneFilter`.
    Assesses results against SRD specs at specified `level`.

    Inputs
    ------
    repo : string
        The repository.  This is generally the directory on disk
        that contains the repository and mapper.
    dataIds : list of dict
        List of `butler` data IDs of Image catalogs to compare to reference.
        The `calexp` cpixel image is needed for the photometric calibration.
    metrics : `dict` or `collections.OrderedDict`
        Dictionary of `lsst.validate.base.Metric` instances. Typically this is
        data from ``validate_drp``\ 's ``metrics.yaml`` and loaded with
        `lsst.validate.base.load_metrics`.
    outputPrefix : str, optional
        Specify the beginning filename for output files.
        The name of each filter will be appended to outputPrefix.
    level : str
        The level of the specification to check: "design", "minimum", "stretch"
    verbose : bool
        Provide detailed output.

    Outputs
    -------
    Names of plot files or JSON file are generated based on repository name,
    unless overriden by specifying `ouputPrefix`.
    E.g., Analyzing a repository "CFHT/output"
        will result in filenames that start with "CFHT_output_".
    The filter name is added to this prefix.  If the filter name has spaces,
        there will be annoyance and sadness as those spaces will appear in the filenames.
    """

    allFilters = set([d['filter'] for d in dataIds])

    if outputPrefix is None:
        outputPrefix = repoNameToPrefix(repo)

    jobs = {}
    for filterName in allFilters:
        # Do this here so that each outputPrefix will have a different name for each filter.
        thisOutputPrefix = "%s_%s_" % (outputPrefix.rstrip('_'), filterName)
        theseVisitDataIds = [v for v in dataIds if v['filter'] == filterName]
        job = runOneFilter(repo, theseVisitDataIds, metrics,
                           outputPrefix=thisOutputPrefix,
                           verbose=verbose, filterName=filterName,
                           **kwargs)
        jobs[filterName] = job

    for filterName, job in jobs.items():
        print('')
        print(bcolors.BOLD + bcolors.HEADER + "=" * 65 + bcolors.ENDC)
        print(bcolors.BOLD + bcolors.HEADER + '{0} band summary'.format(filterName) + bcolors.ENDC)
        print(bcolors.BOLD + bcolors.HEADER + "=" * 65 + bcolors.ENDC)

        for specName in job.spec_levels:
            passed = True

            measurementCount = 0
            failCount = 0
            for m in job.measurements:
                if m.value is None:
                    continue
                measurementCount += 1
                if not m.check_spec(specName):
                    passed = False
                    failCount += 1

            if passed:
                print('Passed {level:12s} {count:d} measurements'.format(
                    level=specName, count=measurementCount))
            else:
                msg = 'Failed {level:12s} {failCount} of {count:d} failed'.format(
                    level=specName, failCount=failCount, count=measurementCount)
                print(bcolors.FAIL + msg + bcolors.ENDC)


def runOneFilter(repo, visitDataIds, metrics, brightSnr=100,
                 medianAstromscatterRef=25, medianPhotoscatterRef=25, matchRef=500,
                 makePrint=True, makePlot=True, makeJson=True,
                 filterName=None, outputPrefix=None,
                 verbose=False,
                 **kwargs):
    """Main executable for the case where there is just one filter.

    Plot files and JSON files are generated in the local directory
        prefixed with the repository name (where '_' replace path separators),
    unless overriden by specifying `outputPrefix`.
    E.g., Analyzing a repository "CFHT/output"
        will result in filenames that start with "CFHT_output_".

    Parameters
    ----------
    repo : string
        The repository.  This is generally the directory on disk
        that contains the repository and mapper.
    dataIds : list of dict
        List of `butler` data IDs of Image catalogs to compare to reference.
        The `calexp` cpixel image is needed for the photometric calibration.
    metrics : `dict` or `collections.OrderedDict`
        Dictionary of `lsst.validate.base.Metric` instances. Typically this is
        data from ``validate_drp``\ 's ``metrics.yaml`` and loaded with
        `lsst.validate.base.load_metrics`.
    brightSnr : float, optional
        Minimum SNR for a star to be considered bright
    medianAstromscatterRef : float, optional
        Expected astrometric RMS [mas] across visits.
    medianPhotoscatterRef : float, optional
        Expected photometric RMS [mmag] across visits.
    matchRef : int, optional
        Expectation of the number of stars that should be matched across visits.
    makePrint : bool, optional
        Print calculated quantities (to stdout).
    makePlot : bool, optional
        Create plots for metrics.  Saved to current working directory.
    makeJson : bool, optional
        Create JSON output file for metrics.  Saved to current working directory.
    outputPrefix : str, optional
        Specify the beginning filename for output files.
    filterName : str, optional
        Name of the filter (bandpass).
    verbose : bool, optional
        Output additional information on the analysis steps.

    """
    if outputPrefix is None:
        outputPrefix = repoNameToPrefix(repo)

    matchedDataset = MatchedMultiVisitDataset(repo, visitDataIds,
                                              verbose=verbose)
    photomModel = AnalyticPhotometryModel(matchedDataset)
    astromModel = AnalyticAstrometryModel(matchedDataset)
    linkedBlobs = {'photomModel': photomModel, 'astromModel': astromModel}

    job = Job()

    PA1 = PA1Measurement(metrics['PA1'], matchedDataset, filterName,
                         job=job, linkedBlobs=linkedBlobs,
                         verbose=verbose)

    for specName in metrics['PA2'].get_spec_names(filter_name=filterName):
        PA2Measurement(metrics['PA2'], matchedDataset,
                       pa1=PA1, filter_name=filterName,
                       spec_name=specName, verbose=verbose,
                       job=job, linkedBlobs=linkedBlobs)

    for specName in metrics['PF1'].get_spec_names(filter_name=filterName):
        PF1Measurement(metrics['PF1'], matchedDataset, PA1, filterName,
                       specName, verbose=verbose,
                       job=job, linkedBlobs=linkedBlobs)

    AM1 = AMxMeasurement(metrics['AM1'], matchedDataset, filterName,
                         job=job, linkedBlobs=linkedBlobs, verbose=verbose)

    for specName in metrics['AF1'].get_spec_names(filter_name=filterName):
        AFxMeasurement(metrics['AF1'], matchedDataset, AM1, filterName, specName,
                       job=job, linkedBlobs=linkedBlobs, verbose=verbose)

        ADxMeasurement(metrics['AD1'], matchedDataset, AM1, filterName, specName,
                       job=job, linkedBlobs=linkedBlobs, verbose=verbose)

    AM2 = AMxMeasurement(metrics['AM2'], matchedDataset, filterName,
                         job=job, linkedBlobs=linkedBlobs, verbose=verbose)

    for specName in metrics['AF2'].get_spec_names(filter_name=filterName):
        AFxMeasurement(metrics['AF2'], matchedDataset, AM2, filterName, specName,
                       job=job, linkedBlobs=linkedBlobs, verbose=verbose)

        ADxMeasurement(metrics['AD2'], matchedDataset, AM2, filterName, specName,
                       job=job, linkedBlobs=linkedBlobs, verbose=verbose)

    AM3 = AMxMeasurement(metrics['AM3'], matchedDataset, filterName,
                         job=job, linkedBlobs=linkedBlobs, verbose=verbose)

    for specName in metrics['AF3'].get_spec_names(filter_name=filterName):
        AFxMeasurement(metrics['AF3'], matchedDataset, AM3, filterName, specName,
                       job=job, linkedBlobs=linkedBlobs, verbose=verbose)

        ADxMeasurement(metrics['AD3'], matchedDataset, AM3, filterName, specName,
                       job=job, linkedBlobs=linkedBlobs, verbose=verbose)

    job.write_json(outputPrefix.rstrip('_') + '.json')

    if makePlot:
        if AM1.value:
            plotAMx(job.get_measurement('AM1'),
                    job.get_measurement('AF1', spec_name='design'),
                    filterName, amxSpecName='design',
                    outputPrefix=outputPrefix)
        if AM2.value:
            plotAMx(job.get_measurement('AM2'),
                    job.get_measurement('AF2', spec_name='design'),
                    filterName, amxSpecName='design',
                    outputPrefix=outputPrefix)
        if AM3.value:
            plotAMx(job.get_measurement('AM3'),
                    job.get_measurement('AF3', spec_name='design'),
                    filterName, amxSpecName='design',
                    outputPrefix=outputPrefix)

        plotPA1(job.get_measurement('PA1'), outputPrefix=outputPrefix)

        plotAnalyticPhotometryModel(matchedDataset, photomModel,
                                    outputPrefix=outputPrefix)

        plotAnalyticAstrometryModel(matchedDataset, astromModel,
                                    outputPrefix=outputPrefix)

    if makePrint:
        print(bcolors.BOLD + bcolors.HEADER + "=" * 65 + bcolors.ENDC)
        print(bcolors.BOLD + bcolors.HEADER +
              '{band} band metric measurements'.format(band=filterName) +
              bcolors.ENDC)
        print(bcolors.BOLD + bcolors.HEADER + "=" * 65 + bcolors.ENDC)

        wrapper = TextWrapper(width=65)

        for metricName in metrics:
            metric = metrics[metricName]
            print(bcolors.HEADER + '{name} - {reference}'.format(
                name=metric.name, reference=metric.reference))
            print(wrapper.fill(bcolors.ENDC + '{description}'.format(
                description=metric.description).strip()))

            for specName in metric.get_spec_names(filter_name=filterName):
                try:
                    m = job.get_measurement(metricName,
                                            spec_name=specName,
                                            filter_name=filterName)
                except RuntimeError:
                    print('\tSkipped {specName:12s} no spec'.format(
                        specName=specName))
                    continue

                if m.value is None:
                    print('\tSkipped {specName:12s} no measurement'.format(
                        specName=specName))
                    continue

                spec = metric.get_spec(specName, filter_name=filterName)
                passed = m.check_spec(specName)
                if passed:
                    prefix = bcolors.OKBLUE + '\tPassed '
                else:
                    prefix = bcolors.FAIL + '\tFailed '
                infoStr = '{specName:12s} {meas:.4f} {op} {spec:.4f} {units}'.format(
                    specName=specName,
                    meas=m.value,
                    op=metric.operator_str,
                    spec=spec.value,
                    units=spec.units)
                print(prefix + infoStr + bcolors.ENDC)

    return job
