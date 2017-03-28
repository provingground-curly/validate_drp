# Configure to force using astrometry_net instead of pre-loaded reference catalogs

from lsst.pipe.tasks.setConfigFromEups import setPhotocalConfigFromEups, setAstrometryConfigFromEups


from lsst.meas.astrom import LoadAstrometryNetObjectsTask
config.processCcd.calibrate.astromRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
config.processCcd.calibrate.photoRefObjLoader.retarget(LoadAstrometryNetObjectsTask)

setPhotocalConfigFromEups(config.processCcd.calibrate.photoCal)

menu = { "ps1*": {}, # No changes
         "sdss*": { "filterMap": {"y": "z"} }, # No y-band, use z instead
         "2mass*": { "filterMap": {ff:"J" for ff in 'grizy'} }, # No optical; use J 
       }

setAstrometryConfigFromEups(config.processCcd.calibrate.astromRefObjLoader, menu)


