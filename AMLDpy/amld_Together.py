#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday July 28 
@author: emilywilliams
"""

############## ALL THAT NEEDS TO BE CHANGED ##############################################

## WHERE THE amld_Functions.py file is located
functionFileLoc = '/Users/emilywilliams/Documents/GitHub/AMLD_CODE/AMLDpy/'
## Folder with .txt Data
rawDatLoc = "/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/allData"

## Folder to put results in (will make subfolders later)
resFolder = "/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/allData/"


## CarID 
xCar = 'SoCCar' # might need to be 5 letters? Need to check that!

## What Proportion above Baseline to flag as elevated (i.e. 0.1 = 10% higher)
threshold = '0.05'

## How many minutes to include in background calculation (minutes)
#timethresh = '1.7'
timethresh = '5.0'
#timethresh = '2.5'

## How many minutes to skip at the beginning of the dataset (i.e. if Collin is at his house)
initialTimeIgnore = '0'

# minimum number of elevated readings required for an observed peak
minElevated = '1'

## Lag time for CH4 to reach sensor (in seconds)
#shift = -4
shift = -4

## Is this an engineering file?
#engineering = True
engineering = False
aeris = False

# Not super sure what timePush is but thats cool
timePush = 5 #min
timePush = 0 

###
backObs = '102'
#backObs = '10'

maxCarSpeed = '45'
minCarSpeed = '2'

baseLinePerc = '50' ## median
#baseLinePerc = '25' #Q1

buff = '30'

###############################################################################
###### DON'T CHANGE ANYTHING BELOW THIS (UNLESS YOU CAN FIX IT) ###############
###############################################################################

# STARTING ALGORITHM (NAMING FOLDERS AND SUCH) 
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
finalMain = finRes + 'mainThing.csv'
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

import sys
sys.path.insert(0, functionFileLoc) # FINDING FUNCTIONS FOLDER TO IMPORT FROM
from amld_Functions import unique,unIfInt,\
                            intersect,weightedLoc,verPk,estEmissions,\
                            haversine,wt_time_Locs,sumthing,makeGEO,\
                            makeGPD,summarizeDat,getQuad,calcTheta,\
                            calcBearing,ProcessRawDataEng,strList,\
                            countTimes,IdentifyPeaks,filterPeak,\
                            passCombine,sumData2,addOdometer,ProcessRawData,ProcessRawDataAeris
import rtree
import pygeos
import numpy as np
import os
import os, sys, datetime, time, math, csv, numpy,gzip,shutil
from math import radians, sin, cos, sqrt, asin
from numpy import log
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import matplotlib.pyplot as plt
import ast
from datetime import datetime
import time, swifter


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

##### THIS PORTION OF THE CODE ALLOWS US TO ITERATIVELY ADD IN MORE DATA
#        PUT THE NEW TEXT FILES INTO THE OVERALL FOLDER AND IT WILL DO THE REST
rawTexts = pd.DataFrame(os.listdir(rawDir)).loc[pd.DataFrame(os.listdir(rawDir))[0].str.endswith('.txt')]


### DONT ADD NEW FILE WITH PRE-EXISTING DATE [NEED TO WRITE CODE TO DEAL WITH THAT]
toAnalyse = []
toIdentify = []
toFilter = []

if os.path.exists(finalInfoLoc):
    analysed = pd.read_csv(finalInfoLoc)
    for index,row in rawTexts.reset_index().iterrows():
        text = rawTexts[0].iloc[index]
        if analysed[analysed['FILENAME'].astype(str).str.contains(text)].shape[0] < 1:
            toAnalyse.append(text)
            toIdentify.append(s1 + '_20' + text[11:17] + '_dat.csv')
            toFilter.append(s2 + '_20' + text[11:17] + '.csv')
elif not os.path.exists(finalInfoLoc):
    for index,row in rawTexts.reset_index().iterrows():
        text = rawTexts[0].iloc[index]
        toAnalyse.append(text)
        toIdentify.append(s1 + '_20' + text[11:17] + '_dat.csv')
        toFilter.append(s2 + '_20' + text[11:17] + '.csv')
        

##### START OF THE ALGORITHM
start = time.time()

if __name__ == '__main__':
    dateList = []
    x1=""
    count = 0
    #listthing = os.listdir(rawDir)
    listthing = toAnalyse ## changing to include the files we have to analyse?
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
            #theResult = ProcessRawData(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc)
            if engineering:
                theResult = ProcessRawDataEng(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc,initialTimeIgnore,shift,maxCarSpeed,minCarSpeed)
            elif not engineering:
                if aeris:
                    theResult = ProcessRawDataAeris(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc,initialTimeIgnore,shift,maxCarSpeed,minCarSpeed)
                elif not aeris:
                    theResult = ProcessRawData(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc,initialTimeIgnore,shift,maxCarSpeed,minCarSpeed)
            #del(bFirst)
            count = count + 1


if __name__ == '__main__':
    #listthing = os.listdir(processedFileLoc).copy() 
    listthing = toIdentify.copy()
    for file in listthing:
        if file.startswith(s1) and file.endswith("dat.csv"):
            xDate = file[len(xCar)+1:len(xCar) + 9]
            theResult = IdentifyPeaks(xCar, xDate, processedFileLoc, file,opDir,processedFileLoc,engineering,threshold,timethresh,minElevated,backObs,baseLinePerc)
if __name__ == '__main__':
    index = 0
    numproc = 0
    listthing = os.listdir(opDir)
    for file in listthing:
        if file.startswith(s2) and (file.endswith('.csv') and not file.endswith('info.csv')):
            print("About to Filter Peaks in the file: " + str(file))
            file_loc = opDir + file
            csv_loc = opDir  + file[:-4]+ '_info.csv'
            nonempt = False
            if pd.read_csv(file_loc).size != 0:
                index += 1
                nonempt = True
                #xDate = file[12:20]
                xDate = file[len(xCar) + 7: len(xCar) + 15]
                filterPeak(xCar,xDate,opDir,file,filtopDir,buffer = buff,whichpass = index )
        
end = time.time()
#print(end - start)

    
if not os.path.exists(finalMain):
    toCombine = os.listdir(filtopDir)
    toCombineList = []
    index = 0
    for file in toCombine:
        if file.startswith('Filtered') and file.endswith('.csv') and not file.endswith('info.csv'):
            toCombineList.append(file)
    index = 1
    for file in toCombineList:
        if index == 1:
            loc = filtopDir + file
            mainThing = (pd.read_csv(loc))
        elif index != 1:
            loc2 = filtopDir + file
            secThing = pd.read_csv(loc2)
            woo = passCombine(mainThing,secThing,xCar, buff)
            mainThing = woo.copy()
        index = index + 1
        print(file)
    mainThing = mainThing.copy().reset_index(drop = True)
    mainThing['numtimes']  = mainThing.apply(lambda x: countTimes(x.recombine,xCar),axis=1)
    #mainThing['numtimes']  = mainThing.recombine.swifter.apply(lambda x: countTimes(x))

# =============================================================================
# if not os.path.exists(finalMain):
#     index = 0
#     numproc = 0
#     listthing = os.listdir(opDir)
#     for file in listthing:
#         if file.startswith(s2) and file.endswith('.csv') and not file.endswith('info.csv'):
#             file_loc = opDir + file
#             csv_loc = opDir + "/" + file[:-4] + '_info.csv'
#             
#             nonempt = False
#             if pd.read_csv(file_loc).size != 0:
#                 index += 1
#                 nonempt = True
#             xDate = file[12:20]
#             if index == 1 and nonempt:
#                 mainThing = filterPeak(xCar,xDate,opDir,file,filtopDir,whichpass = index )
#                 numproc += 1
#                 mainInfo = pd.read_csv(csv_loc)
#                 #mainInfo = mainInfo1.copy()
#     
#             if index != 1 and nonempt:
#                 secondThing = filterPeak(xCar,xDate,opDir,file,filtopDir,whichpass = index )
#                 mainThing = passCombine(mainThing,secondThing)
#                 tempInfo = pd.read_csv(csv_loc)
#                 mainInfo = pd.concat([mainInfo, tempInfo], axis=0)
#     toCombine = os.listdir(filtopDir)
#     toCombineList = []
#     for file in toCombine:
#         if file.startswith(s3) and file.endswith('.csv') and not file.endswith('info.csv'):
#             toCombineList.append(file)
#     for x in range(len(toCombineList)):
#         if x == 0: ## first thing
#             mainThing = toCombineList[x]
#             mainLoc = filtopDir + mainThing
# 
#             print('first mainthing = ' + str(toCombineList[x]))
#         elif x!= 0:
#             secondThing = toCombineList[x]
#             secLoc = filtopDir + secondThing
#             print('combining = ' + 'main thing and ' + str(toCombineList[x]))
#             #mainThing = passCombine(mainThing,secondThing)
# =============================================================================
                
## figure out the mainInfo thing
#mainInfo.reset_index().FILENAME.to_csv(finalInfoLoc)
mainThing.reset_index(drop=True).to_csv(finalMain)
        
print("I processed "+ str(index) + ' days of driving. The processed files are now stored in the folder: ' + str(filtopDir))

################

#combined = summarizeDat(mainThing) ## finds locations and mean log ch4 for each peak (either verified or non yet)
combined = sumData2(mainThing) ## finds locations and mean log ch4 for each peak (either verified or non yet)


## combined so only with the same overall peak
uniquePk = combined.loc[:,['min_read']].drop_duplicates()
uniqueList = combined.loc[uniquePk.index,['min_read','recombine']]
uniqueOther = combined.loc[:,['min_read','overallLON','overallLAT','mnlogCH4',\
                             'verified','numtimes','minDist','maxDist']].drop_duplicates()

unique_gdf = makeGPD(uniqueOther,'overallLAT','overallLON')
unique_gdf2 = pd.merge(unique_gdf,uniqueList,on = ['min_read'])

#unique_gdf = makeGPD(combined.loc[:,['min_read','pk_LAT','pk_LON']],'pk_LAT','pk_LON')
#combinedGeo = unique_gdf.dissolve(by='min_read',as_index = False)


#allTog = pd.merge(combinedGeo,uniqueOther,on=['min_read'])
#allTog = pd.merge(allTog,uniqueList,on=['min_read'])
allTog = unique_gdf2.copy()
#allTog['em'] = allTog.apply(lambda y: estEmissions(y['mnlogCH4']),axis=1)
allTog['em'] = allTog['mnlogCH4'].swifter.apply(lambda y: estEmissions(y))



#allTog['threshold'] = allTog.apply(lambda x: threshold,axis =1 )

allTog['threshold'] = allTog['em'].swifter.apply(lambda x: threshold)


##### SPLITTING IF THE PEAKS WERE VERIFIED OR NOT
verTog = allTog.loc[allTog.numtimes!= 1,:]

if verTog.size > 0:
    verTog.drop(columns=['recombine']).to_file(shpFileLocName, driver="GeoJSON")
    print('I found ' + str(len(verTog.min_read.unique()))+" verified peaks")

if verTog.size ==0:
    print("Sorry, no verified peaks were found.")  
if allTog.size> 0:
    allTog.drop(columns=['recombine']).to_file(OPshpFileLocName, driver="GeoJSON")
    allTog.to_csv(allPksCSVLoc)

if allTog.size == 0:
    print("Sorry, no observed peaks were found in the given data")
    
end = time.time()
print("I analysed the data using a threshold of " + str(float(threshold)*100 + 100) + "% for an elevated reading" )
print("where the threshold was calculated using the " + str(baseLinePerc) + 'th percentile over ' + str(backObs) + ' observations')
print("I filtered the speed of the car to be between " + str(minCarSpeed) + 'mph and ' + str(maxCarSpeed) + 'mph')
print("To create an observed peak, I required there to be a minimum of " + str(minElevated) + " observations within 30 seconds")
print("I created three summary files located here: " + str(finRes) + ".")
print("The processing took " + str(round((end-start)/60,3)) + str(" minutes."))
print("I found " + str(len(mainThing.min_read.unique()))+ " Observed Peaks")
