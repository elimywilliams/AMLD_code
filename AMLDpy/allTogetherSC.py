#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 09:21:43 2020
@author: emilywilliams
"""

############## ALL THAT NEEDS TO BE CHANGED ##############################################

## WHERE THE allFunctionsSC.py file is located
functionFileLoc = '/Users/emilywilliams/Documents/GitHub/AMLD_CODE/AMLDpy/'
## what you want the car to be named
xCar = 'SCcar' # might need to be 5 letters? need to check into that
## Folder with .txt Data
rawDatLoc = "/Users/emilywilliams/Documents/DrivingData/ColDat" 

## Folder to put results in (will make subfolders later)
resFolder = "/Users/emilywilliams/Documents/DrivingData/ColDat/"

############################################################################################

################### ASSIGNING FOLDERS FOR RESULTS
## WHERE TO PUT LEAKS
rawDir =  resFolder + 'RawData/'
#inDir = resFolder + 'ObservedPeaks/'
opDir = resFolder + 'ObservedPeaks/'

#outFolder = resFolder + 'FilteredObservedPeaks/'
filtopDir = resFolder + 'FilteredObservedPeaks/'
finRes = resFolder + 'FinalShpFiles/'
shpFileLocName = finRes + 'verifiedPKs.json'
processedFileLoc = resFolder + 'ProcessedData/'
OPshpFileLocName = finRes + "OP_Final.json"

## replace inDir with opDir

#inDir = "/Users/emilywilliams/Documents/DrivingData/ColDat/"

# where you want the 
#outFolder = "/Users/emilywilliams/Documents/DrivingData/ColDat/Processed/"
#shpFileLocName = "/Users/emilywilliams/Documents/DrivingData/ColDat/compFinal.json"
#processedFileLoc = "/Users/emilywilliams/Documents/DrivingData/ColDat/"
#OPshpFileLocName = "/Users/emilywilliams/Documents/DrivingData/ColDat/OP_Final.json"


####### POTENTIALLY COULD CHANGE
s1 = xCar
s2 = "Peaks_" + str(s1)
s3 = "Filtered" + str()


##################################
### IMPORTING NECESSARY MODULES
##################################

import contextily as ctx
import pandas as pd
import os, sys, datetime, time, math, csv, numpy,gzip,shutil
sys.path.insert(1, functionFileLoc) ##change this to location of the "allFunctions.py" file
from math import radians, sin, cos, sqrt, asin
from numpy import log
import geopandas as gpd
import pandas as pd #
from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
import matplotlib.pyplot as plt
from datetime import datetime
from allFunctionsSC import IdentifyPeaks,filterPeak,unique,unIfInt,IsInPK,intersect,passCombine, weightedLoc,verPk,estEmissions,haversine,ProcessRawData,wt_time_Locs,sumthing


### create paths

#### move files to new path

#### CREATING NECESSARY FOLDERS

foldList = [rawDir,resFolder,opDir,filtopDir,finRes,processedFileLoc]
for x in foldList:
    if os.path.isdir(x) == False:
        try:
            os.mkdir(x)
        except OSError:
            print ("Creation of the directory %s failed" % x)
        else:
            print ("Successfully created the directory %s " % x)

### MOVING RAW FILES TO THE RAW DATA FILE FOLDER
listthing = os.listdir(rawDatLoc)

for file in listthing:
    if file.endswith(".txt"):
        shutil.move(rawDatLoc+'/' + file,rawDir)

########################################################################################



##### START OF THE ALGORITHM
start = time.time()

if __name__ == '__main__':  

    dateList = []
    x1=""
    count = 0
    listthing = os.listdir(rawDir)
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
            theResult = ProcessRawData(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc)
            #del(bFirst)
            count = count + 1


if __name__ == '__main__':
    listthing = os.listdir(processedFileLoc).copy() 
    for file in listthing:
        if file.startswith(s1) and file.endswith("dat.csv"):
            xDate = file[6:14]
            theResult = IdentifyPeaks(xCar, xDate, processedFileLoc, file,opDir)

index = 0
numproc = 0
listthing = os.listdir(opDir).copy()

for file in listthing:
    if file.startswith(s2) and file.endswith('.csv'):
        file_loc = opDir + file
        nonempt = False
        if pd.read_csv(file_loc).size != 0:
            index += 1
            nonempt = True
        xDate = file[12:20]
        if index == 1 and nonempt:
            mainThing = filterPeak(xCar,xDate,opDir,file,filtopDir,whichpass = index )
            numproc += 1

        if index != 1 and nonempt:
            secondThing = filterPeak(xCar,xDate,opDir,file,filtopDir,whichpass = index )
            mainThing = passCombine(mainThing,secondThing)

        
print("I processed "+ str(index) + ' days of driving. The processed files are now stored in the folder: ' + str(filtopDir))

###########
#combinedDat = mainThing.copy()

## plotwith the thing
combinedDatPk = mainThing.loc[:,['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']].drop_duplicates().reset_index().loc[:,['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']]
wcomb = weightedLoc(combinedDatPk,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB')
wcomb2 = wcomb.loc[wcomb.pk_LAT.isna()==False,:]
del(wcomb)
wcomb = wcomb2.loc[:,:].rename(columns={'pk_LAT': 'overall_LAT','pk_LON':'overall_LON'})


###### CREATE MAP THING (AND HOPE IT WORKS) ##########################################
together = pd.merge(wcomb,combinedDatPk,on = ['min_read'])
totmaxCH4 = together.groupby('min_read',as_index = False).pk_maxCH4_AB.max().rename(columns = {'pk_maxCH4_AB':'VP_maxch4ab'})


togetherMore2 = pd.merge(together,totmaxCH4,on = ['min_read'])
togetherMore = togetherMore2.loc[:,['min_read','numtimes','overall_LAT','overall_LON','VP_maxch4ab']].drop_duplicates()


together = together.loc[:,['min_read','overall_LAT','overall_LON','PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes']]
together['em'] = together.apply(lambda y: estEmissions(y['pk_maxCH4_AB']),axis=1)

### verified
togetherVer = together.loc[together.numtimes!=1,:]

# creating a geometry column 
geometry = [Point(xy) for xy in zip(togetherVer['overall_LON'], togetherVer['overall_LAT'])]
crs = {'init': 'epsg:4326'}

# Creating a Geographic data frame 
gdf = gpd.GeoDataFrame(togetherVer, crs=crs, geometry=geometry)
gdf = gdf.to_crs(epsg = 3857)
gdf.to_file(shpFileLocName, driver="GeoJSON")


### SAVING ALL OBSERVED PEAKS (NOT JUST VERIFIED)

# creating a geometry column 
geometry = [Point(xy) for xy in zip(together['overall_LON'], together['overall_LAT'])]
gdf = gpd.GeoDataFrame(together, crs=crs, geometry=geometry)
gdf = gdf.to_crs(epsg = 3857)

gdf.to_file(OPshpFileLocName, driver="GeoJSON")
end = time.time()
print("I created a GeoJSON file with the verified peak information, as well as the observed peak information, and it is located here: " + str(shpFileLocName) + ". The processing took " + str(round((end-start)/60,3)) + str(" minutes."))


