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
  'cv_z': 0.020213447,
  'rcv_z': -0.241673763,
  'ndvi_z': 0.171443325,
  'ndsi_z': 0.242237481,
  'ndwi_z': -0.021797543,
  'nbr_z': -0.004327796,
  'lda': 5
});

scrub = ee.Dictionary({
  'int': 0,
  'cv_z': -0.004616377,
  'rcv_z': -0.756532635,
  'ndvi_z': -0.034184409,
  'ndsi_z': -0.005930105,
  'ndwi_z': -0.208912147,
  'nbr_z': -0.065727911,
  'lda': 5
});
desert = ee.Dictionary({
  'int': 0,
  'cv_z': 0.03535703,
  'rcv_z': 0.37298094,
  'ndvi_z': 0.81159062,
  'ndsi_z': 0.31934502,
  'ndwi_z': -0.04074573,
  'nbr_z': 0.03053187,
  'lda': 1.6
});
wetland = ee.Dictionary({
  'int': 0,
  'cv_z': 0.01539318,
  'rcv_z': -0.53627410,
  'ndvi_z': -0.05276547,
  'ndsi_z': 0.04784782,
  'ndwi_z': -0.21047313,
  'nbr_z': 0.28553036,
  'lda': 2.5
});

grassland = ee.Dictionary({
  'int': 0,
  'cvz': 0.003335746,
  'nbr_z': 0.200660811,
  'ndsi_z': -0.048249043,
  'ndvi_z': 0.199691164,
  'ndwi_z': 0.012553511,
  'rcvmax_z': -0.473307836,
  'lda': 5
});
