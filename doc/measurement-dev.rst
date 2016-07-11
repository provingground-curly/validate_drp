############################
Creating Measurement Classes
############################


The ``lsst.validate.drp.base.MeasurmentBase`` abstract base class defines standard interface for writing classes that make measurements of metrics.
The ``MeasurementBase`` baseclass ensures that measurements, along with metadata about the measurement, can be serialized and submitted to the QA database and SQUASH dashboard.
Indeed, this page will focus on the details of properly persisting information about a measurement.

A Minimal Measurement Class
===========================

At a minimum, measurement classes must subclass ``MeasurementBase`` and provide some metadata, such as the name of metric being measured and units of the measurement.
This example is a basic template for a measurement:

.. code-block:: python

   class PA1Measurement(MeasurementBase):

       metric = None
       value = None
       units = 'mmag'
       label = 'PA1'
       
       def __init__(self):
           MeasurementBase.__init__(self)
           
           self.metric = Metric.fromYaml(self.label)
           
           # measurement code
           # ...
           
           self.value = 0  # Scalar value from measurement code

In a measurement class, the following metadata attributes must be specified (their presence is required by the ``MeasurementBase`` abstract base class).

``label``
   The name of the metric.

``metric``
   A ``Metric`` object. In this example, the `Metric.fromYaml` class method constructs a ``Metric`` instance for PA1 from the ``metrics.yaml`` file built into ``validate_drp``.

``units``
   The units of the metric measurement. As in the ``Datum`` class, units should be an ``astropy.units``-compatible string.

The main responsibility of a measurement class is to make a make a measurement.
Measurement should occur during a measurement instance's ``__init__`` method.
Any data required for a measurement should be provided through the measurement class's `__init__` method.

The measurement result stored in ``value``:

``value``
   The value attribute should be a scalar (`float` or `int`), in the same physical units indicated by the ``units`` attribute.
   If a Measurement class is unable to make a measurement, ``value`` should be ``None``.

Storing Measurement Parameters
==============================

Often a measurement code is customized with parameters.
As a means of lightweight provenance, the measurement API provides a way to declare these parameters so that they're persisted to the database using the ``Metric.registerParameter`` method:

.. code-block:: python

   class PA1Measurement(MeasurementBase):

       metric = None
       value = None
       units = 'mmag'
       label = 'PA1'
       schema = 'pa1-1.0.0'
       
       def __init__(self, numRandomShuffles=50):
           MeasurementBase.__init__(self)
           
           self.metric = Metric.fromYaml(self.label)

           self.registerParameter('numRandomShuffles',
                                  value=numRandomShuffles,
                                  units='',
                                  description='Number of random shuffles')
           
           # ... measurement code
                              
In this example, the ``PA1Measurement`` class has a parameter ``numRandomShuffles``.

Accessing parameter values as object attributes
-----------------------------------------------

In addition to registering a parameter for serialization, the ``registerParameters`` method makes the values of parameters available as instance attributes:

.. code-block:: python

   pa1 = PA1Measurment(numRandomShuffles=50)
   pa1.numRandomShuffles # == 50
   
Through attribute access, a parameter's value can be both *read* and *updated*.

Accessing parameters as ``Datum`` objects
-----------------------------------------

Although the values of parameters can be accessed through object attributes, they are stored internally as ``Datum`` objects.
These full ``Datum`` objects can be accessed as items of the ``parameters`` attribute:

.. code-block:: python

   pa1.parameters['numRandomShuffles'].value  # 50
   pa1.parameters['numRandomShuffles'].units  # ''
   pa1.parameters['numRandomShuffles'].label  # numRandomShuffles
   pa1.parameters['numRandomShuffles'].description  # 'Number of random shuffles'

Alternative ways of registering parameters
------------------------------------------

The ``registerParameters`` method is flexible in terms of its arguments.

It's possible to first register a parameter and set its value later:

.. code-block:: python

   self.registerParameter('numRandomShuffles', units='', description=Number of random shuffles')
   # ...
   self.numRandomShuffles = 50

In this example, a label is not set; in this case the ``label`` defaults to the name of the parameter itself.

It's also possible to provide a ``Datum`` to ``registerParameters``:

.. code-block:: python

   self.registerParameter('numRandomShuffles',
                          datum=Datum(50, '', label='shuffles,
                                      description='Number of random shuffles'))

This can be useful when copying a parameter already available as a ``Datum``.

Storing Extra Measurement Outputs
=================================

Although metric measurements are strictly scalar values, it can be useful to store additional measurement by-products.
By registering them, these measurement by-products are automatically serialized with the measurement and available the SQUASH dashboard application.
This allows the dashboard to make rich plots, such as histograms or scatter plots, that help a user understand a scalar metric measurement.

Registering measurement outputs is similar to registering parameters, except that the `registerExtra` method is used.

As an example PA1 measurement code stores the inter-quartile range, RMS of magnitude differences for each random sample, along with the magnitude differences and mean magnitude of each pair of observed stars from each sample.

.. code-block:: python

   class PA1Measurement(MeasurementBase):
   
          metric = None
          value = None
          units = 'mmag'
          label = 'PA1'
          schema = 'pa1-1.0.0'
          
          def __init__(self, numRandomShuffles=50):
              MeasurementBase.__init__(self)
              
              self.metric = Metric.fromYaml(self.label)
              
              # register extras
              self.registerExtra('rms', units='mmag',
                                 description='Photometric repeatability RMS of '
                                             'stellar pairs for each random sampling')
              self.registerExtra('iqr', units='mmag',
                                 description='Photometric repeatability IQR of '
                                             'stellar pairsfor each random sample')
              self.registerExtra('magDiff', units='mmag',
                                 description='Difference magnitudes of stellar source pairs'
                                             'for each random sample')
              self.registerExtra('magMean', units='mag',
                                 description='Mean magnitude of pairs of stellar '
                                             'sources matched across visits, for '
                                             'each random sample.')

              # ... make measurements
              
              # Set values of extras
              self.rms = np.array([pa1.rms for pa1 in pa1Samples])
              self.iqr = np.array([pa1.iqr for pa1 in pa1Samples])
              self.magDiff = np.array([pa1.magDiffs for pa1 in pa1Samples])
              self.magMean = np.array([pa1.magMean for pa1 in pa1Samples])
       
              # The scalar metric measurement
              self.value = np.mean(self.iqr)

The registerExtra method works just like the registerParameter method.
The value of the extra can be set at registration time.
An extra can also be registered with a pre-made ``Datum`` object.

Accessing and updating the values and Datum objects of measurement extras
-------------------------------------------------------------------------

As with parameters, registering an extra allows the value of the extra to be accessed or updated through a measurement object attribute named after the extra itself (see the above example).

Extras are internally stored as ``Datum`` objects, which can be accessed as items of the ``extras`` attribute.
Following the PA1 measurement example:

.. code-block:: python

   pa1 = PA1Measurement()
   pa1.extras['rms'].value  # == pa1.rms
   pa1.extras['rms'].units  # 'mmag'
   pa1.extras['rms'].label  # 'rms'
   pa1.extras['rms'].decription  # 'Photometric repeatability RMS ...'
