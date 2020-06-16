#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 14:27:01 2020

@author: emilywilliams
"""
### trying a thing where you make a grid
import numpy as np
import pandas as pd
loc = '/Users/emilywilliams/Documents/DrivingData/ColDatShort/ProcessedData/SCcar_20200406_dat.csv'
df = pd.read_csv(loc)


bottomLeft = (df.LAT.min(), df.LONG.min())
bottomRight = (df.LAT.min(), df.LONG.max())
topLeft = (df.LAT.min(), df.LONG.max())
topRight = (df.LAT.max(), df.LONG.max())


cols = np.linspace(bottomLeft[1], bottomRight[1], num=18)
rows = np.linspace(bottomLeft[0], topLeft[0], num=15)
df['col'] = np.searchsorted(cols, df['LONG'])
df['row'] = np.searchsorted(rows, df['LAT'])

## trying another thing

listthing = os.listdir('/Users/emilywilliams/Documents/DrivingData/ColDat/ProcessedData/').copy() 
index = 0
for file in listthing:
    if file.startswith(s1) and file.endswith("dat.csv"):
        if index == 0:
            first = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/ColDat/ProcessedData/' + file)
            index += 1
        elif index != 0:
            second = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/ColDat/ProcessedData/' + file)
            first = pd.concat([first, second], axis=0)
            index += 1

firstgeo = makeGPD(first,'LAT','LONG')
firstgeo.crs = {'init': 'epsg:32610'}
firstgeo = firstgeo.to_crs(epsg=32610)
firstgeo.to_file('allPTS.JSON',driver='GeoJSON')

import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np
#points = gpd.read_file('points.shp')
points = firstgeo.copy()
xmin,ymin,xmax,ymax =  points.total_bounds


topdist = haversine(ymax,xmin,ymax,xmax)
botdist = haversine(ymin,xmin,ymin,xmax)

leftdist = haversine(ymin,xmin,ymax,xmin)
rightdist = haversine(ymin,xmax,ymax,xmax)

cols_top = np.ceil(topdist/100)
cols_bot = np.ceil(botdist/100)
cols_num = max(cols_top,cols_bot)
cols = int(cols_num)

rows_top = np.ceil(leftdist/100)
rows_bot = np.ceil(rightdist/100)
rows_num = max(rows_top,rows_bot)

height = (ymax - ymin)/rows_num
width = (xmax -xmin)/cols_num
rows = int(rows_num)

#width = 4e-3
#height = 4e-3
#rows = int(np.ceil((ymax-ymin) /  height))
#cols = int(np.ceil((xmax-xmin) / width))

XleftOrigin = xmin
XrightOrigin = xmin + width
YtopOrigin = ymax
YbottomOrigin = ymax- height
polygons = []
for i in range(cols):
    Ytop = YtopOrigin
    Ybottom =YbottomOrigin
    for j in range(rows):
        polygons.append(Polygon([(XleftOrigin, Ytop), (XrightOrigin, Ytop), (XrightOrigin, Ybottom), (XleftOrigin, Ybottom)])) 
        Ytop = Ytop - height
        Ybottom = Ybottom - height
    XleftOrigin = XleftOrigin + width
    XrightOrigin = XrightOrigin + width

grid = gpd.GeoDataFrame({'geometry':polygons})
grid.crs = {'init': 'epsg:32610'}
grid = grid.to_crs(epsg=32610)

grid.to_file("grid.shp")
grid.to_file('grid.JSON', driver="GeoJSON")


##### 
import geopandas as gpd

points.crs = {'init': 'epsg:32610'}
points = points.to_crs(epsg=32610)


joined = gpd.sjoin(points, polygons, how="inner", op='intersects')

# Read the data.
polygons = grid.copy()
points2 = points.iloc[0:5000,:]
# Make a copy because I'm going to drop points as I
# assign them to polys, to speed up subsequent search.
pts = points2.copy() 

# We're going to keep a list of how many points we find.
pts_in_polys = []

# Loop over polygons with index i.
for i, poly in polygons.iterrows():

    # Keep a list of points in this poly
    pts_in_this_poly = []

    # Now loop over all points with index j.
    for j, pt in pts.iterrows():
        if poly.geometry.contains(pt.geometry):
            # Then it's a hit! Add it to the list,
            # and drop it so we have less hunting.
            pts_in_this_poly.append(pt.geometry)
            pts = pts.drop([j])

    # We could do all sorts, like grab a property of the
    # points, but let's just append the number of them.
    pts_in_polys.append(len(pts_in_this_poly))

# Add the number of points for each poly to the dataframe.
polygons['number of points'] = gpd.GeoSeries(pts_in_polys)

#############
import geopandas as gpd
import pandas as pd

polys = grid.copy()
polys['PolyID'] = polys.index
points2 = points.copy()

dfsjoin = gpd.sjoin(polys,points2) #Spatial join Points to polygons
dfpivot = pd.pivot_table(dfsjoin,index='PolyID',columns='Food',aggfunc={'Food':len})
dfpivot.columns = dfpivot.columns.droplevel()

dfpolynew = polys.merge(dfpivot, how='left',on='PolyID')

polys = gpd.read_file('zip:///Users/emilywilliams/Downloads/TGWgas_Polygon16/Main16.zip')
polys.crs = {'init': 'epsg:32610'}
polys = polys.to_crs(epsg=32610)
polys.to_file('main16.JSON',driver = 'GeoJSON')

#####

grid = glt.gridify_data(gdf_points, 8000,'time', method=np.min)  # Choose a method to aggregate point values in each grid cell
    
m_plot_dataframe(grid, column='time', contour_poly_width=0.5, edgecolor="grey")
