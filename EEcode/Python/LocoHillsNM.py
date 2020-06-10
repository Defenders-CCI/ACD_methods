# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 10:33:39 2020

@author: MEvans
"""

import ee, analyze, dictionaries

#credentials = ee.ServiceAccountCredentials('change-detection-1@change-detection-dow.iam.gserviceaccount.com', 'change-detection-dow-e29ba5eb4f51.json')
#ee.Initialize(credentials)
ee.Authenticate()
ee.Initialize()


# This test case was generated 2020-06-08. No changes appear in the imagery. We want to confirm.

# give this aoi a name
testId = 'LocoHillsNM'

# TODO: split up analyze functions so we don't need these
# grab the relevant dictionary of lda coefficients
dictionary = dictionaries.scrub

# This was the GEE javascript aoi
#aoi = ee.Geometry.Polygon(
#        [[[-104.07540229751181, 32.785947073773244],
#          [-104.07540229751181, 32.719248277130205],
#          [-103.9593592066915, 32.719248277130205],
#          [-103.9593592066915, 32.785947073773244]]])

# THis is the aoi from deployed habitat patrol
aoi = ee.Geometry.Polygon(
        [[[-104.084129,32.786986],
          [-104.085503,32.740938],
          [-104.011688,32.741515],
          [-104.012718,32.789728],
          [-104.084129,32.786986]]])

doi = '2020-01-01' 

landcover = 'forest'

output = analyze.analyze_iw(
    ee.Feature(aoi, {'mode':landcover}),
    doi,
    dictionary,
    0,
    testId)

print(output[3].size().getInfo())