# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:10:00 2020

@author: MEvans
"""

import ee, analyze, dictionaries

ee.Initialize()

# give this aoi a name
testId = 'HughesMillCA'

# TODO: split up analyze functions so we don't need these
# grab the relevant dictionary of lda coefficients
dictionary = dictionaries.forest

aoi = ee.Geometry.Polygon(
        [[[-120.83127262496578, 39.10457008576222],
          [-120.83127262496578, 39.06952752960459],
          [-120.76518299483882, 39.06952752960459],
          [-120.76518299483882, 39.10457008576222]]], None, False);

doi = '2017-07-01' 

landcover = 'forest'

output = analyze.analyze_iw(
    ee.Feature(aoi, {'mode':landcover}),
     doi, dictionary, 0, testId)

print(output[3].size().getInfo())