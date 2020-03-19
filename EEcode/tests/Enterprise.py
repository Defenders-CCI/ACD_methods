# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 17:14:54 2019

@author: MEvans
"""

import folium
import ee
import ee.mapclient
import analyze

ee.Initialize()

def add_ee_layer(self, ee_image_object, vis_params, name):
  map_id_dict = ee.Image(ee_image_object).getMapId(vis_params)
  folium.raster_layers.TileLayer(
    tiles = map_id_dict['tile_fetcher'].url_format,
    attr = "Map Data © Google Earth Engine",
    name = name,
    overlay = True,
    control = True
  ).add_to(self)

# Add EE drawing method to folium.
folium.Map.add_ee_layer = add_ee_layer

aoi = ee.Geometry.Polygon(
        [[[-115.31969349704991, 35.98044115914671],
          [-115.2175549777628, 35.977662882104646],
          [-115.2175549777628, 36.044453329376914],
          [-115.31677525364171, 36.04473092617846]]])
    

doi = ee.Date('2018-05-12')

output = analyze.analyze_mad(aoi, doi, 0, 5)

vis_params = {
        'palette': ['ffffff', 'd63000'],
        'min': 80,
        'max': 200}

#ee.mapclient.addToMap(ee.Image(output.get('chi2')))

#ALTERNATIVE LEAFLET VISUALIZATION


my_map = folium.Map(location=[36.025, -115.228], zoom_start=12)

# Add the elevation model to the map object.
my_map.add_ee_layer(output.get('chi2'), vis_params, 'chi2')

# Add a layer control panel to the map.
my_map.add_child(folium.LayerControl())

# Display the map.
my_map.save(outfile = 'mad5.html')