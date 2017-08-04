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

from astropy.table import Column, Table

from lsst.validate.base import load_metrics
from lsst.validate.drp.validate import get_filter_name_from_job, load_json_output


def run(validation_drp_report_filenames, output_file,
        srd_level=None,
        release_metrics_file=None, release_level=None):
    """
    Parameters
    ---
    validation_drp_report_filenames : [] or str, filepaths for JSON files.
    output_file : str, filepath of output RST file.
    srd_level : str, SRD level to quote.  One of ['design', 'minimum', 'stretch']
    release_metrics_file : str [optional], filepath of metrics YAML file
      While the JSON file itself will store the metrics the was calculated with
        one may wish to compare against an external set of specifications.
      E.g., ${VALIDATE_DRP_DIR}/etc/release_metrics.yaml
        contains the target levels for each fiscal year through to and including
        operational readiness review (ORR).  This file is essentially the same
        as ${VALIDATE_DRP_DIR}/etc/metrics.yaml, but defines levels in terms
        of fiscal year and ORR instead of design/minimum/stretch.
    release_level : str, A specification level in the 'release_metrics_file'
       E.g., 'FY17' or 'ORR'

    Products
    ---
    Writes table of performance metrics to an RST file.
    """
    input_objects = ingest_data(validation_drp_report_filenames)
    input_table = objects_to_table(input_objects, level=srd_level)
    if input_table is None:
        msg = "Table from Job is None.  Returning without writing table"
        print(msg)
        return

    if release_metrics_file is not None and release_level is not None:
        release_metrics = load_metrics(release_metrics_file)
        add_release_metric(input_table, release_metrics, release_level)

    write_report(input_table, output_file)


def ingest_data(filenames):
    """Load JSON files into a list of lsst.validate.base measurement Jobs.

    Parameters
    ----------
    filenames : list of str
        Filenames of JSON files to load.

    Returns
    -------
    list of lsst.validate.base.Job
        Each element is the Job representation of the JSON file.
    """
    jobs = {}
    # Read in JSON output from metrics run
    for file in filenames:
        job = load_json_output(file)
        filter_name = get_filter_name_from_job(job)
        jobs[filter_name] = job

    return jobs


# Identify key data from JSON
def objects_to_table(input_objects, level='design'):
    """Take lsst.validate.base.Job objects and convert to astropy.table.Table

    Parameters
    ----------
    input_objects : list of Job objects
    level : str
        The requirement level to compare to for each metric.
        This is required because there are metrics with dependencies
        E.g., AD1 depends on AF1.  Thus you have to specify a level
        to even get a single AD1 metric.

    Returns
    -------
    astropy.table.Table
        Table with columns needed for final report.
    """
    rows = []
    for filter_name, obj in input_objects.items():
        for meas in obj.measurements:
            # Skip specification levels (called .name in measurement objects)
            #  that are not the 'level' we are looking for.
            if meas.spec_name is not None and meas.spec_name != level:
                continue
            m = meas.metric
            try:
                spec = m.get_spec(level, filter_name=filter_name)
            except:
                msg_format = "Could not load meas.spec_name: '{:s}' at level: '{:s}'"
                print(msg_format.format(meas.spec_name, level))
                continue
            if meas.quantity is None:
                meas_quantity_value = "--"
            else:
                meas_quantity_value = meas.quantity.value
            this_row = [m.name, filter_name, meas_quantity_value, meas.unit,
                        m.operator_str, spec.quantity.value]
            rows.append(this_row)

    if len(rows) == 0:
        msg_format = "No rows loaded from Job at level: '{:s}'"
        print(msg_format.format(level))
        return None

    srd_requirement_col_name = 'SRD Requirement: %s' % level
    col_names = ('Metric', 'Filter', 'Value', 'Unit',
                 'Operator', srd_requirement_col_name)
    output = Table(rows=rows, names=col_names)
    output.add_column(Column(['']*len(output), dtype=str, name='Comments'))
    return output


# Calculate numbers in table
def add_release_metric(data, release_metrics, release_metrics_level):
    release_targets = []
    for row in data:
        metric = release_metrics[row['Metric']]
        spec = metric.get_spec(
            name=release_metrics_level, filter_name=row['Filter'])
        release_targets.append(spec.quantity.value)

    release_targets_col = Column(
        release_targets,
        dtype=float,
        name='Release Target: %s' % release_metrics_level)
    data.add_column(release_targets_col)


def float_or_dash(f, format_string='{:.2f}'):
    """Return string of formatted float, or -- if None.

    Intended use is to provide formatting output for columns
    where where None or non-float value indicates a missing measurement.
    """
    # This try/except handles both None and non-numeric strings.
    try:
        return format_string.format(float(f))
    except:
        return '--'


def blank_none(s):
    """Return a blank for  None or 'None', else return string of input.

    Intended use is to provide formatting output for columns where an empty
    or None value is totally reasonable and expected and should be display
    as a blank ''.
    """
    if s is None:
        return ''
    if s == 'None':
        return ''

    return str(s)


def find_col_name(prefix, colnames):
    """Return the first entry in 'colnames' that starts with 'prefix'."""
    for c in colnames:
        if c.startswith(prefix):
            return c


# Output table
def write_report(data, filename='test.rst', format='ascii.rst'):
    """Write performance report to RST file.

    Parameters
    ----------
    data : astropy.table.Table
    filename : str [optional]
        Filepath of output RST file.
    format : str [optional]
        astropy.table format for output table.

    Creates
    -------
    test.rst
        Output file with RST version of data Table.
    """
    # Find the 'Release Target XYZ' column name
    release_target_col_name = find_col_name('Release Target', data.colnames)
    # Find the 'SRD Requirement XYZ' column name
    srd_requirement_col_name = find_col_name('SRD Requirement', data.colnames)

    col_names = ['Metric', 'Filter', 'Unit', 'Operator',
                 srd_requirement_col_name,
                 release_target_col_name,
                 'Value', 'Comments']
    use_col_names = [c for c in col_names if c in data.colnames]
    # Provide default formats
    for spec_col in (release_target_col_name, srd_requirement_col_name):
        if spec_col in data:
            data[spec_col].info.format = '.1f'
    data['Value'].info.format = float_or_dash
    data['Unit'].info.format = blank_none
    data[use_col_names].write(filename=filename, format=format,
                              include_names=use_col_names,
                              overwrite=True)
