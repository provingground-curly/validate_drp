################################################
Providing Datasets to Measurements through Blobs
################################################

A common pattern is to reduce a raw dataset, and share that dataset between several measurements.
The measurement API allows such datasets to be expressed as 'blobs,'
In the context of the measurement API, a blob is an object that contains ``Datum`` objects.

There are several advantages with storing input datasets in blob objects:

- Any pre-reduction of a raw dataset can be done in the blob object, which helps to keep an organized codebase.
- Blobs can be passed to to measurement objects, which simplifies the construction of measurements.
- Blobs are automatically serialized alongside measurements and are available to the SQUASH dashboard. Blobs can be shared among several measurements, with the blob data only be stored once.

Template for a Blob Class
=========================

Blobs are subclasses of `BlobBase` that register one or more datums.

.. code-block:: python

   class SimpleBlob(BlobBase):
   
   def __init__(self, gMags, iMags):
       self.registerDatum('g', value=gMags, units='mag',
                          description='g-band magnitudes')
       self.registerDatum('i', value=iMags, units='mag',
                          description='r-band magnitudes')
       self.registerDatum('gi', units='mag',
                          description='g-i colour')
       
       self.gi = self.g - self.r

In this example, the g and i attributes are initially registered with values.
A third blob attribute, ``gi`` is also declared and it's value is computed afterwards.

Notice that, like parameters and extras of measurement classes, the values of Datum objects can be accessed and updated directly through class attributes.

Accessing datum objects
-----------------------

Internally, blob attributes are stored as Datum objects, which can be accessed as items of the `datums` attribute.

.. code-block:: python

   blob = SimpleBlob(g, i)
   blob.datums['gi'].value  # == blob.gi
   blob.datums['gi'].units  # 'mag'
   blob.datums['gi'].label  # 'gi', this was automatically set from the name
   blob.datums['gi'].description  # g-i colour

Linking measurements to blobs
-----------------------------

When a blob is used by a measurement, the measurement class should declare that usage so that the SQUASH dashboard can provide rich context to measurements.
Measurement classes can accomplish this by making the blob an instance attribute.
For example:

.. code-block:: python

   class MeanColor(MeasurementBase):
   
   label = 'MeanColour'
   units = 'mag'
   
   def __init__(self, simpleBlob):
       self.metric = Metric.fromYaml(self.label)
       self.simpleBlob = simpleBlob
       self.value = np.mean(self.simpleBlob.gi)

Accessing blobs in measurements
-------------------------------

In addition to simply accessing blobs associated with a measurment through the instance attribute, blobs are also available as items of the measurement's `blobs` attribute:

.. code-block:: python

   color = SimpleBlob(g, i)
   meanColor = MeanColor(color)
   meanColor.blobs['simpleBlob'].gi  # array of g-i colours
