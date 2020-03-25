# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 12:41:37 2020

@author: MEvans
"""

import ee
ee.Initialize()

def buffer(ft):
    """
    Buffer an input feature by 10km
    """
    return ft.buffer(10000)

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
    countries = ee.FeatureCollection("USDOS/LSIB/2013")
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

    # generate 1000 random points within the us
    random = ee.FeatureCollection.randomPoints(
      region = usa,
      points = n
    )
    
    coords = random.geometry(5).coordinates(5)

    # create 10k radius boxes around random points
    boxes = random.map(buffer_and_id)
    
    # add the 
    aois = nlcd.select('landcover').reduceRegions(
            collection = boxes,
            reducer = ee.Reducer.mode(),
            scale = 30)
    
    # filter polygons that primarily include ag or pasture
    aois = aois.filter(
            ee.Filter.And(
                    ee.Filter.neq('landcover', 81),
                    ee.Filter.neq('landcover', 82)
                    )
            )
    
    return aois, coords

if __name__ == '__main__':
    print('How many aois should we create?')
    size = input()
    create_acd_aois(size)