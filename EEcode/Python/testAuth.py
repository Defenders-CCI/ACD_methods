# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 10:48:53 2020

@author: MEvans
"""

import ee

## INITIALIZE WITH SERVICE ACCOUNT
#credentials = ee.ServiceAccountCredentials('change-detection-1@change-detection-dow.iam.gserviceaccount.com', 'change-detection-dow-e29ba5eb4f51.json')
#ee.Initialize(credentials)

## INITIALIZE WITH USER ACCOUNT
ee.Initialize()

img = ee.Image('JRC/GSW1_1/YearlyHistory/2018')
print(img.bandNames().getInfo())

imgCol = ee.ImageCollection("JRC/GSW1_1/YearlyHistory")
print(imgCol.size().getInfo())