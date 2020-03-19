# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 09:21:46 2019

@author: MEvans
"""
#we probably don't need to import ee here, could create these as python dictionaries...?
import ee

#Collection of dictionaries storing parameter values for various habitat types
forest = ee.Dictionary({
  'int': 0,
  'cvz': 0.020213447,
  'rcvz': -0.241673763,
  'ndviz': 0.171443325,
  'ndsiz': 0.242237481,
  'ndwiz': -0.021797543,
  'nbrz': -0.004327796,
  'lda': 5
});

scrub = ee.Dictionary({
  'int': 0,
  'cvz': -0.004616377,
  'rcvz': -0.756532635,
  'ndviz': -0.034184409,
  'ndsiz': -0.005930105,
  'ndwiz': -0.208912147,
  'nbrz': -0.065727911,
  'lda': 5
});
desert = ee.Dictionary({
  'int': 0,
  'cvz': 0.03535703,
  'rcvz': 0.37298094,
  'ndviz': 0.81159062,
  'ndsiz': 0.31934502,
  'ndwiz': -0.04074573,
  'nbrz': 0.03053187,
  'lda': 1.6
});
wetland = ee.Dictionary({
  'int': 0,
  'cvz': 0.01539318,
  'rcvz': -0.53627410,
  'ndviz': -0.05276547,
  'ndsiz': 0.04784782,
  'ndwiz': -0.21047313,
  'nbrz': 0.28553036,
  'lda': 2.5
});

grassland = ee.Dictionary({
  'int': 0,
  'cvz': 0.003335746,
  'nbrz': 0.200660811,
  'ndsiz': -0.048249043,
  'ndviz': 0.199691164,
  'ndwiz': 0.012553511,
  'rcvmaxz': -0.473307836,
  'lda': 5
});
