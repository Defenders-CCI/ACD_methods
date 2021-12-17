# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 10:48:53 2020

@author: MEvans
"""

import ee

## INITIALIZE WITH SERVICE ACCOUNT
key = 'C:/Users/mevans/OneDrive - Defenders of Wildlife/Analyses/ACD/acd-app-04a10fe38611.json'
credentials = ee.ServiceAccountCredentials('acd-app@appspot.gserviceaccount.com', key)
ee.Initialize(credentials)

## INITIALIZE WITH USER ACCOUNT
#ee.Initialize()

img = ee.Image('JRC/GSW1_1/YearlyHistory/2018')
print(img.bandNames().getInfo())

imgCol = ee.ImageCollection("JRC/GSW1_1/YearlyHistory")
print(imgCol.size().getInfo())