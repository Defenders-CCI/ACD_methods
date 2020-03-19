# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 11:23:17 2019

@author: MEvans
This is a test equivalent to the Hazard, KY study area
https://code.earthengine.google.com/12179c0c5593b6be557b1bd43a28eba5
"""
import numpy
import ee
import ee.mapclient
import analyze
import dictionaries

ee.Initialize()

aoi = ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                    [[[-80.81991758407679, 33.88836406786576],
                      [-80.82417447706393, 33.81222611807673],
                      [-80.71508881286002, 33.80919082945717],
                      [-80.72153527034902, 33.881822895147785]]]
                    ),
            {
              "system:index": "0"
            })])

projdate = ee.Date('2017-08-20')
dictionary = dictionaries.wetland
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
