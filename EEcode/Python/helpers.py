# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 12:41:37 2020

@author: MEvans
"""

import ee
from random import randint
ee.Initialize()

def buffer(ft):
    """
    Buffer an input feature by 10km
    """
    return ft.buffer(10000).bounds(5)

def set_id(ft):
    return ft.set('id', ft.id())

def buffer_and_id(ft):
    return buffer(set_id(ft))

def create_acd_aois(n):
    """
    Generate random 20km2 boxes within the USA with NLCD landcover data
    Parameters:
        n (int): number of boxes to create
    Returns:
        Tuple <ee.FeatureCollection, ee.List>: output boxes with properties
        'id' and 'landcover', and a list of centroid coordinates.
    """
    # get landcover data
    nlcd = ee.Image("USGS/NLCD/NLCD2016");
    # get world country boundaries
    states = ee.FeatureCollection("TIGER/2018/States")
    # create list of excluded us states
    exclusions = ee.List(
            ['Alaska', 'Hawaii', 'Guam', 'American Samoa', 'Puerto Rico', 'United States Virgin Islands', 'Commonwealth of the Northern Mariana Islands']
            )
    # create a single polygon representing the continental USA
    usa = ee.Feature(
            states.filter(
                    ee.Filter.inList('NAME', exclusions).Not())\
                    .union()\
                    .first()
                    )
            
    seed = randint(0, 100)
    
    # generate 1000 random points within the us
    random = ee.FeatureCollection.randomPoints(
      region = usa.geometry(),
      points = n,
      seed = seed
    )
    
    # get list of aoi centroid coordinates for metadata
    coords = random.geometry(5).coordinates()

    # create 10k radius boxes around random points
    boxes = random.map(buffer_and_id)
    
    # add the most common landcover within aoi as a property
    aois = nlcd.select('landcover').reduceRegions(
            collection = boxes,
            reducer = ee.Reducer.mode(),
            scale = 30)
    
    # filter polygons that primarily include ag or pasture
    aois = aois.filter(
            ee.Filter.And(
                    ee.Filter.neq('mode', 81),
                    ee.Filter.neq('mode', 82)
                    )
            )
    
    return aois, coords

def sample_output(output, polygons, aoi_id):
    """
    Collect algorithm output metrics within polgons and export as csv
    Parameters:
        output (ee.Image): n-band image output of either iw or mad
        polgyons (ee.FeatureCollection): polygons within which to sample
        aoi_id (str): unique identifier for aoi to be used as output file name
    Returns:
        export task sending a csv file to google drive
    """

    data = output.sampleRegions(
            collection = polygons,
            properties = [],
            scale = 10,
            tileScale = 8)
    
    task = ee.batch.Export.table.toDrive(
            collection = data,
            #folder = {GDrive subdirectory path},
            description = aoi_id,
            fileFormat = 'CSV'
            )
    
    task.start()
    
if __name__ == '__main__':
    print('How many aois should we create?')
    size = input()
    create_acd_aois(size)