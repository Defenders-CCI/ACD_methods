# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 11:23:17 2019

@author: MEvans
This is a test equivalent to the dune sagebrush lizard study area
https://code.earthengine.google.com/12179c0c5593b6be557b1bd43a28eba5
"""

import ee
import ee.mapclient
import analyze
import dictionaries

ee.Initialize()
    
aoi = ee.Geometry.Polygon(
        [[[-103.12765923269012, 32.06015197325721],
          [-103.08508721120575, 31.91980036812035],
          [-102.94089165456512, 31.98097580591193],
          [-102.98689690358856, 32.10262271189585]]])

projdate = ee.Date('2017-08-20')

dictionary = dictionaries.scrub
cvz = ee.Number(dictionary.get('cvz'))
rcvz = ee.Number(dictionary.get('rcvz'))
ndviz = ee.Number(dictionary.get('ndviz'))
ndsiz = ee.Number(dictionary.get('ndsiz'))
ndwiz = ee.Number(dictionary.get('ndwiz'))
nbrz = ee.Number(dictionary.get('nbrz'))
lda = ee.Number(dictionary.get('lda'))
intercept = ee.Number(dictionary.get('int'))
cd_id = 'test'

output = analyze.analyze_iw(aoi, projdate, cvz, nbrz, ndsiz, ndviz, ndwiz, rcvz, 0, intercept, lda)

print(output)

ee.mapclient.addToMap(output[3])

task = ee.batch.Export.table.toDrive(
  collection= output[3],
  description= studyarea + '_pythonPolys',
  fileFormat= 'GeoJSON'
)
