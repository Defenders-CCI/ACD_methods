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
  'rcvmax_z': -0.241673763,
  'ndvi_z': 0.171443325,
  'ndsi_z': 0.242237481,
  'ndwi_z': -0.021797543,
  'nbr_z': -0.004327796,
  'lda': 5
});

scrub = ee.Dictionary({
  'int': 0,
  'cv_z': 0.001523492,#-0.004616377,
  'rcvmax_z': -0.228383542,#-0.756532635,
  'ndvi_z': 0.374194520,#-0.034184409,
  'ndsi_z': -0.177040740,#-0.005930105,
  'ndwi_z': -0.135408066,#-0.208912147,
  'nbr_z': -0.241886898,#-0.065727911,
  'lda': 0.8
});
desert = ee.Dictionary({
  'int': 0,
  'cv_z': 0.03535703,
  'rcvmax_z': 0.37298094,
  'ndvi_z': 0.81159062,
  'ndsi_z': 0.31934502,
  'ndwi_z': -0.04074573,
  'nbr_z': 0.03053187,
  'lda': 1.6
});
wetland = ee.Dictionary({
  'int': 0,
  'cv_z': 0.01539318,
  'rcvmax_z': -0.53627410,
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
  chi  0.0195215749
v1   0.0001151189
v2  -0.0006378971
v3  -0.0004738995
v4  -0.0007036161
v5   0.0004734193
v6  -0.0011570211
});

iwwetland
cv_z      0.02608609
nbr_z     0.15379691
ndsi_z    0.09082271
ndvi_z    0.09873193
ndwi_z   -0.12424875
rcvmax_z -0.11168147

madwetland
chi  3.183977e-02
v1  -5.495463e-04
v2  -1.935902e-05
v3  -5.771866e-05
v4   3.543529e-04
v5   8.904052e-04
v6  -2.347725e-03

iwgrassland
cv_z      0.01501676
nbr_z     0.09690559
ndsi_z   -0.08688555
ndvi_z   -0.14835658
ndwi_z   -0.19044280
rcvmax_z -0.20351300

madgrassland


madScrub
chi  0.0112356429
v1  -0.0019951636
v2   0.0002209054
v3   0.0004750912
v4  -0.0015628258
v5  -0.0017628213
v6   0.0005914612

