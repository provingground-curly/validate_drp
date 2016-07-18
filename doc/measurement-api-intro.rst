################################################
Introduction to the validate_drp Measurement API
################################################

``validate_drp`` provides a framework for making and reporting measurements of metrics.
This framework is used internally by ``validate_drp`` and can also be used by other packages that monitor algorithmic performance.

The measurement API is focussed on making metric measurements and metadata available to the SQUASH databoard.
All data containers in the API can serialize themselves in a JSON format that can be submitted to the SQUASH dashboard's web API.
Datasets are also designed to be self-documenting, both from Python and JSON contexts.
Values are annotated with units (astropy.units-compatible), and readable-descriptions.

The measurement API also features a YAML format for defining metrics and specification levels.
Metric objects are constructed from YAML definitions, and make it easy to compare a measurement against a specification.

Main Classes in the Measurement API
===================================

The measurement API consists of classes that can either be used directly, or subclassed.
API classes provide a consistent pattern for defining and measuring metrics.
By using these classes, users get JSON serialization 'for-free', as well the ability to easily grade measurements against specification levels.

These are the key API classes, along with links to further documentation:

- ``Datum`` wraps all numerical data, whether it is a specification level, parameter, measurement value, or member of a 'blob' data class. Datum objects contain the following fields:

  - A ``value`` field, which can be a scalar (int or float) or a sequence (list or numpy.ndarray)
  - A ``units`` field, that annotates the physical units of the ``value``. Units must be astrpy.units-compatible strings. For unitless quantities, units should be an empty string.
  - A ``label`` field, which can be used to decorate a plot axis or legend. The ``label`` should exclude reference of units. This ``label`` field is optional.
  - A ``description`` field that can include free form text that documents the Datum. The ``description`` field is also optional.
  
- ``Metric`` is a class containing metadata about a metric, along with specification levels. Typically Metric objects are constructed from a YAML definition file, though they can also be arbitrarily constructed in Python. See :doc:`metric-dev` for more information.
- ``Specification`` objects are contained inside Metric objects, and define levels of expected measurement performance. Metrics can have multiple specifications indicating different performance goals (e.g., 'minimum', 'design' and 'stretch). Specifications can also be associated with certain observational bandpasses. Typically Specification objects are built automatically by the Metric class from a YAML definition. Specifications can also be added manually to Metric classes. See :doc:`metric-dev` for more inforkmation.
- ``MeasurementBase`` is a base class for making measurements of a metric. All code needed for making a measuerment can be contained by a subclass of `MeasurementBase`. At a minimum, MeasurementBase subclasses store a scalar value, but can also register addition Datum objects that can be persisted to JSON. See :doc:`measurement-dev` for more information.
- ``BlobBase`` is a base class for making Blob classes. Blobs are a way of storing datasets in an object that is both convenient for measurement classes and serializable for JSON. Measurement classes can share a common Blob without duplicating data stored in the SQUASH database. Blobs are linked to the measurements that use them, allowing blobs to power plots that provide context to measurements. See :doc:`blob-dev` for more information.