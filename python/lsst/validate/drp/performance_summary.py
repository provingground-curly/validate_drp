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

from lsst.validate.drp.validate import get_filter_name_from_job, load_json_output


def run(validation_drp_report_filenames, output_file, debug=True):
    input_objects = ingest_data(validation_drp_report_filenames)
    input_table = objects_to_table(input_objects)
    output_table = calculate_numbers(input_table)
    write_report(output_table, output_file)


def ingest_data(filenames):
    """"""
    jobs = {}
    # Read in JSON output from metrics run
    for file in filenames:
        job = load_json_output(file)
        filter_name = get_filter_name_from_job(job)
        jobs[filter_name] = job

    return jobs


# Identify key data from JSON
def objects_to_table(input_objects, level='design'):
    """Take objects and convert to table."""
    rows = []
    for filter_name, obj in input_objects.items():
        for meas in obj.measurements:
            m = meas.metric
            spec = m.get_spec(level, filter_name=filter_name)
            if meas.quantity is None:
                meas_quantity_value = "--"
            else:
                meas_quantity_value = meas.quantity.value
            this_row = [m.name, filter_name, meas_quantity_value, meas.unit,
                        m.operator_str, spec.quantity.value, spec.quantity.value]
            rows.append(this_row)

    col_names = ('Metric', 'Filter', 'Value', 'Unit',
                 'Operator', 'Release Target', 'SRD Requirement')
    output = Table(rows=rows, names=col_names)
    output.add_column(Column(['']*len(output), dtype=str, name='Comments'))
    return output


# Calculate numbers in table
def calculate_numbers(input_table):
    updated_table = input_table
    return updated_table


# Output table
def write_report(data, filename='test.rst', format='ascii.rst'):
    col_names = ['Metric', 'Unit', 'SRD Requirement',
                 'Release Target', 'Value', 'Comments']
    data[col_names].write(filename=filename, format=format,
                          include_names=col_names,
                          overwrite=True)
