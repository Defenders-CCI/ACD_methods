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
        [[[-85.38630927322764, 29.715555535604025],
          [-85.3816744160499, 29.705715285963432],
          [-85.37897074936289, 29.70668444426876],
          [-85.37824118851084, 29.710113698554586],
          [-85.38115943191904, 29.717120940944707]]]);

doi = '2018-06-01' 

landcover = 'forest'

output = analyze.analyze_iw(
    ee.Feature(aoi, {'mode':landcover}),
     doi, dictionary, 0, testId)

print(output[3].size().getInfo())