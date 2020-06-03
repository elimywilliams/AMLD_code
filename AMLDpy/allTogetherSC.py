#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 09:21:43 2020

@author: emilywilliams
"""

#### ALL THAT NEEDS TO BE CHANGED
functionFileLoc = '/Users/emilywilliams/Documents/GitHub/AMLD_CODE/AMLDpy/'
xCar = 'ColLi'
inDir = "/Users/emilywilliams/Documents/DrivingData/ColDat/"
outFolder = "/Users/emilywilliams/Documents/DrivingData/ColDat/Processed/"
rawDatLoc = "/Users/emilywilliams/Documents/DrivingData/ColDat" ##where the .txt files are located
shpFileLocName = "/Users/emilywilliams/Documents/DrivingData/ColDat/compFinal.json"
processedFileLoc = "/Users/emilywilliams/Documents/DrivingData/ColDat/"
OPshpFileLocName = "/Users/emilywilliams/Documents/DrivingData/ColDat/OP_Final.json"


####### starting the 
s1 = "COLLi"
s2 = "Peaks_" + str(s1)
s3 = "Filtered" + str()


##################################
### IMPORTING NECESSARY MODULES
##################################

import contextily as ctx
import os, sys, datetime, time, math, csv, numpy,gzip
sys.path.insert(1, functionFileLoc) ##change this to location of the "allFunctions.py" file
from math import radians, sin, cos, sqrt, asin
import geopandas as gpd
import pandas as pd #
from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
import matplotlib.pyplot as plt
from datetime import datetime
from allFunctionsSC import IdentifyPeaks,filterPeak,unique,unIfInt,IsInPK,intersect,passCombine, weightedLoc,verPk,estEmissions,haversine,ProcessRawData,wt_time_Locs,sumthing


###### POTENTIALLY COULD CHANGE I SUPPOSE
s1 = "COLLi" ## the start of the file that you used
s2 = "Peaks_" + str(s1)
s3 = "Filtered" + str()

##### START OF THE ALGORITHM
xCar = s1
start = time.time()

if __name__ == '__main__':

    inDirRaw = rawDatLoc
    outDir = inDir

    dateList = []
    x1=""
    count = 0
    listthing = os.listdir(inDirRaw)
    for file in listthing:
        #print (file)
        if file.endswith(".txt"):
            dateF = file[11:17]
            if dateF in set(dateList):
                bFirst = False
            
            if dateF not in set(dateList):
                bFirst = True
                dateList.append(dateF)

            
            #if file[:10] == x1:
             #   bFirst = False
            #else:
             #   bFirst = True
            x1 = file[:10]
            
            #xCar = "CSULi"
            xDate = file[:10]
            theResult = ProcessRawData(xCar, xDate, inDirRaw, file, bFirst, 1, outDir)
            #del(bFirst)
            count = count + 1


if __name__ == '__main__':
    listthing = os.listdir(inDir)
    
    #for file in os.listdir(inDir):
    for file in listthing:
        if file.startswith(s1) and file.endswith("dat.csv"):
            xCar = s1
            xDate = file[6:14]
            print(file)
            theResult = IdentifyPeaks(xCar, xDate, inDir, file)

index = 0
numproc = 0
listthing = os.listdir(inDir)

for file in listthing:
    if file.startswith(s2) and file.endswith('.csv'):
        file_loc = inDir + file
        nonempt = False
        if pd.read_csv(file_loc).size != 0:
            index += 1
            nonempt = True
        #print(index)
        #print(file)
        #print(file)
        woo = file
        xCar = s1
        xDate = file[12:20]
        if index == 1 and nonempt:
            mainThing = filterPeak(xCar,xDate,inDir,file,outFolder,whichpass = index )
            first = mainThing
            firstfile = file
            #print(file)
            numproc += 1

        if index != 1 and nonempt:
            secondThing = filterPeak(xCar,xDate,inDir,file,outFolder,whichpass = index )
            mainThing = passCombine(mainThing,secondThing)
            #print(file)
            #print('i combined things')
           # print(numproc)
        #print(index)
        
        
print("Hurray I processed "+ str(index) + ' files. They are now stored in the folder: ' + str(outFolder))

###########
combinedDat = mainThing.copy()



## plotwith the thing
combinedDatPk = combinedDat[['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']].drop_duplicates().reset_index()[['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']]

wcomb = weightedLoc(combinedDatPk,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB')
wcomb = wcomb[wcomb.pk_LAT.isna()==False]
wcomb = wcomb.rename(columns={'pk_LAT': 'overall_LAT','pk_LON':'overall_LON'})



###### CREATE MAP THING (AND HOPE IT WORKS) ##########################################
together = pd.merge(wcomb,combinedDatPk,on = ['min_read'])
 

totmaxCH4 = together.groupby('min_read',as_index = False).pk_maxCH4_AB.max().rename(columns = {'pk_maxCH4_AB':'VP_maxch4ab'})


togetherMore = pd.merge(together,totmaxCH4,on = ['min_read'])
togetherMore = togetherMore[['min_read','numtimes','overall_LAT','overall_LON','VP_maxch4ab']].drop_duplicates()


together = together[['min_read','overall_LAT','overall_LON','PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes']]
together['em'] = together.apply(lambda y: estEmissions(y['pk_maxCH4_AB']),axis=1)
together = together[together.numtimes!=1]

# creating a geometry column 
geometry = [Point(xy) for xy in zip(together['overall_LON'], together['overall_LAT'])]

# Coordinate reference system : WGS84
crs = {'init': 'epsg:4326'}

# Creating a Geographic data frame 
gdf = gpd.GeoDataFrame(together, crs=crs, geometry=geometry)
gdf = gdf.to_crs(epsg = 3857)


###################################################################
################## SAVE THE FILE TO A SHAPE FILE (FOR EXPORT TO ARCMAP PERHAPS)######################################
gdf.to_file(shpFileLocName, driver="GeoJSON")
end = time.time()
print("The processing of the data is complete. I processed " + str(index) + " files and converted them to a GeoJSON file, which is located here: " + str(shpFileLocName) + ". The processing took " + str(round((end-start)/60,3)) + str(" minutes."))





### SAVING ALL OBSERVED PEAKS (NOT JUST VERIFIED)
together = pd.merge(wcomb,combinedDatPk,on = ['min_read'])
 

totmaxCH4 = together.groupby('min_read',as_index = False).pk_maxCH4_AB.max().rename(columns = {'pk_maxCH4_AB':'VP_maxch4ab'})


togetherMore = pd.merge(together,totmaxCH4,on = ['min_read'])
togetherMore = togetherMore[['min_read','numtimes','overall_LAT','overall_LON','VP_maxch4ab']].drop_duplicates()


together = together[['min_read','overall_LAT','overall_LON','PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes']]
together['em'] = together.apply(lambda y: estEmissions(y['pk_maxCH4_AB']),axis=1)

# creating a geometry column 
geometry = [Point(xy) for xy in zip(together['overall_LON'], together['overall_LAT'])]

# Coordinate reference system : WGS84
crs = {'init': 'epsg:4326'}

# Creating a Geographic data frame 
gdf = gpd.GeoDataFrame(together, crs=crs, geometry=geometry)
gdf = gdf.to_crs(epsg = 3857)


###################################################################
################## SAVE THE FILE TO A SHAPE FILE (FOR EXPORT TO ARCMAP PERHAPS)######################################
gdf.to_file(OPshpFileLocName, driver="GeoJSON")
end = time.time()
print("The processing of the data is complete. I processed " + str(index) + " files and converted them to a GeoJSON file, which is located here: " + str(shpFileLocName) + ". The processing took " + str(round((end-start)/60,3)) + str(" minutes."))


