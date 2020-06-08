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
rawDatLoc = "/Users/emilywilliams/Documents/DrivingData/ColDatShort" 

## Folder to put results in (will make subfolders later)
resFolder = "/Users/emilywilliams/Documents/DrivingData/ColDatShort/"

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
allPksCSVLoc = finRes + 'overallPeaks.csv'
finalInfoLoc = finRes + 'summaryInfo.csv'
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

threshold = '0.1'
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
from allFunctionsSC import IdentifyPeaks,filterPeak,unique,unIfInt,IsInPK,intersect,\
                            passCombine, weightedLoc,verPk,estEmissions,haversine,\
                           ProcessRawData,wt_time_Locs,sumthing,summarizeDat,\
                           makeGEO, makeGPD


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
            theResult = IdentifyPeaks(xCar, xDate, processedFileLoc, file,opDir,processedFileLoc,threshold)

index = 0
numproc = 0
listthing = os.listdir(opDir).copy()

for file in listthing:
    if file.startswith(s2) and file.endswith('.csv') and not file.endswith('info.csv'):
        file_loc = opDir + file
        csv_loc = opDir + "/" + file[:-4] + '_info.csv'
        
        nonempt = False
        if pd.read_csv(file_loc).size != 0:
            index += 1
            nonempt = True
        xDate = file[12:20]
        if index == 1 and nonempt:
            mainThing = filterPeak(xCar,xDate,opDir,file,filtopDir,whichpass = index )
            numproc += 1
            mainInfo = pd.read_csv(csv_loc)
            #mainInfo = mainInfo1.copy()

        if index != 1 and nonempt:
            secondThing = filterPeak(xCar,xDate,opDir,file,filtopDir,whichpass = index )
            mainThing = passCombine(mainThing,secondThing)
            tempInfo = pd.read_csv(csv_loc)
            mainInfo = pd.concat([mainInfo, tempInfo], axis=0)
            
mainInfo.reset_index().FILENAME.to_csv(finalInfoLoc)
        
print("I processed "+ str(index) + ' days of driving. The processed files are now stored in the folder: ' + str(filtopDir))

################
combined = summarizeDat(mainThing) ## finds locations and mean log ch4 for each peak (either verified or non yet)

## combined so only with the same overall peak
uniquePk = combined.loc[:,'min_read'].drop_duplicates()
uniqueList = combined.loc[uniquePk.index,['min_read','recombine']]
uniqueOther = combined.loc[:,['min_read','overallLON','overallLAT','mnlogCH4',
                             'verified','numtimes']].drop_duplicates()

unique_gdf = makeGPD(combined.loc[:,['min_read','pk_LAT','pk_LON']],'pk_LAT','pk_LON')
combinedGeo = unique_gdf.dissolve(by='min_read',as_index = False)


allTog = pd.merge(combinedGeo,uniqueOther,on=['min_read'])
allTog = pd.merge(allTog,uniqueList,on=['min_read'])
allTog['em'] = allTog.apply(lambda y: estEmissions(y['mnlogCH4']),axis=1)


##### SPLITTING IF THE PEAKS WERE VERIFIED OR NOT
verTog = allTog.loc[allTog.numtimes!= 1,:]

verTog.drop(columns=['recombine']).to_file(shpFileLocName, driver="GeoJSON")
allTog.drop(columns=['recombine']).to_file(OPshpFileLocName, driver="GeoJSON")
allTog.to_csv(allPksCSVLoc)

end = time.time()
print("I created three summary files located here: " + str(finRes) + ". The processing took " + str(round((end-start)/60,3)) + str(" minutes."))
