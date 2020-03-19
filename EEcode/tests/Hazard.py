# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 11:23:17 2019

@author: MEvans
This is a test equivalent to the Hazard, KY study area
https://code.earthengine.google.com/2a1ed3f6bcfd360aae96f9d788ff0985
"""
import numpy
import ee
import ee.mapclient
import analyze
import dictionaries

ee.Initialize()

aoi = ee.Geometry.Polygon(
        [[[-83.37017153264765, 37.48081395204879],
          [-83.37486536622822, 37.31933288374584],
          [-83.05319468739128, 37.30974497135589],
          [-83.05556035288436, 37.47934591201635]]])

projdate = ee.Date('2018-08-01')
dictionary = dictionaries.forest
cvz = ee.Number(dictionary.get('cvz'))
rcvz = ee.Number(dictionary.get('rcvz'))
ndviz = ee.Number(dictionary.get('ndviz'))
ndsiz = ee.Number(dictionary.get('ndsiz'))
ndwiz = ee.Number(dictionary.get('ndwiz'))
nbrz = ee.Number(dictionary.get('nbrz'))
lda = ee.Number(dictionary.get('lda'))
intercept = ee.Number(dictionary.get('int'))
cd_id = 'test'

output = analyze.analyze(aoi, projdate, cvz, nbrz, ndsiz, ndviz, ndwiz, rcvz, 0, intercept, lda, cd_id)

print(output)

ee.mapclient.addToMap(output[3])

ee.Export.table.toDrive({
  collection: output[3],
  description: studyarea + '_pythonPolys',
  fileFormat: 'GeoJSON'
})

arr = numpy.random.rand(49, 49)
red = arr[6:-6, 6:-6]
print (red.shape)
xlist = [5,2,3]
ylist = [4,6,7,1]
ylen = len(ylist)
cor = [x*y for x in xlist for y in ylist]
print(cor)
nested = [cor[i:i+ylen] for i in range(0, len(cor), ylen)]
print(chr(65))
print(['training' + chr(i) + 'tfrecord.gz' for i in range(65,74)])
