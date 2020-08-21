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
rawDatLoc = "/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/trussville_dat"

## Folder to put results in (will make subfolders later)
resFolder = "/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/trussville_dat/"

xCar = 'TrussCar' #CAR NAME TO APPEAR IN FILENAMES OBSERVED PEAK NAMES
threshold = '0.05'  #What Proportion above Baseline to flag as elevated (i.e. 0.1 = 10% higher)
timethresh = '5.0'  ## How many minutes to include in background calculation (minutes)
initialTimeIgnore = '0' ## How many minutes to skip at the beginning of the dataset (i.e. if Collin is at his house)
minElevated = '1' # minimum number of elevated readings required for an observed peak
shift = 0  ## Lag time for CH4 to reach sensor (in seconds)
engineering = False #is this an engineering file
aeris = True # is this from the aeris instrument
timePush = 0 #not sure what this is
backObs = '102' ### NUMBER OF OBSERVATIONS TO INCLUDE IN THE BACKGROUND
maxCarSpeed = '45' #maximum car speed to allow (mph)
minCarSpeed = '2' # minimum car speed to allow (mph)

baseLinePerc = '50' #what percentile to use as a backgorund calculation
buff = '30' # distance of buffer (m) to use

###############################################################################
###### DON'T CHANGE ANYTHING BELOW THIS (UNLESS YOU CAN FIX IT) ###############
###############################################################################

# STARTING ALGORITHM (NAMING FOLDERS AND SUCH) 
## WHERE TO PUT LEAKS
rawDir =  resFolder + 'RawData/'
opDir = resFolder + 'ObservedPeaks/'

filtopDir = resFolder + 'FilteredObservedPeaks/'
finRes = resFolder + 'FinalShpFiles/'
shpFileLocName = finRes + 'verifiedPKs.json'
processedFileLoc = resFolder + 'ProcessedData/'
OPshpFileLocName = finRes + "OP_Final.json"
allPksCSVLoc = finRes + 'overallPeaks.csv'
finalInfoLoc = finRes + 'summaryInfo.csv'
finalMain = finRes + 'mainThing.csv'

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

import rtree, pygeos,os, sys, datetime, time, math, numpy, csv, gzip,shutil,ast,swifter
from math import radians, sin, cos, sqrt, asin
import numpy as np
from numpy import log
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import matplotlib.pyplot as plt
from datetime import datetime


#### CREATING NECESSARY FOLDERS
addingFiles = False
foldList = [rawDir,resFolder,opDir,filtopDir,finRes,processedFileLoc]
for x in foldList:
    if os.path.isdir(x) == False:
        try:
            os.mkdir(x)
        except OSError:
            print("Creation of the directory %s failed" % x)
        else:
            print("Successfully created the directory %s " % x)
            
### MOVING RAW FILES TO THE RAW DATA FILE FOLDER
for file in os.listdir(rawDatLoc):
    if file.endswith(".txt"):
        shutil.move(rawDatLoc+'/' + file,rawDir)
########################################################################################

##### THIS PORTION OF THE CODE ALLOWS US TO ITERATIVELY ADD IN MORE DATA
#        PUT THE NEW TEXT FILES INTO THE OVERALL FOLDER AND IT WILL DO THE REST
rawTexts = pd.DataFrame(os.listdir(rawDir)).loc[pd.DataFrame(os.listdir(rawDir))[0].str.endswith('.txt')]

### DON'T ADD NEW FILE WITH PRE-EXISTING DATE [NEED TO WRITE CODE TO DEAL WITH THAT]
toAnalyse, toIdentify, toFilter = [[] for _ in range(3)]


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
    for file in toAnalyse:
        #print (file)
        if file.endswith(".txt"):
            dateF = file[11:17]
            if dateF in set(dateList):
                bFirst = False
            
            if dateF not in set(dateList):
                bFirst = True
                dateList.append(dateF)
            x1 = file[:10]
            xDate = file[:10]
            if engineering:
                theResult = ProcessRawDataEng(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc,initialTimeIgnore,shift,maxCarSpeed,minCarSpeed)
            elif not engineering:
                if aeris:
                    theResult = ProcessRawDataAeris(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc,initialTimeIgnore,shift,maxCarSpeed,minCarSpeed)
                elif not aeris:
                    theResult = ProcessRawData(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc,initialTimeIgnore,shift,maxCarSpeed,minCarSpeed)
            count = count + 1

if __name__ == '__main__':
    for file in toIdentify:
        if file.startswith(s1) and file.endswith("dat.csv"):
            xDate = file[len(xCar)+1:len(xCar) + 9]
            theResult = IdentifyPeaks(xCar, xDate, processedFileLoc, file,opDir,processedFileLoc,engineering,threshold,timethresh,minElevated,backObs,baseLinePerc)
if __name__ == '__main__':
    index = 0
    numproc = 0
    for file in os.listdir(opDir):
        if file.startswith(s2) and (file.endswith('.csv') and not file.endswith('info.csv')):
            file_loc = opDir + file
            csv_loc = opDir  + file[:-4]+ '_info.csv'
            nonempt = False
            if pd.read_csv(file_loc).size != 0:
                index += 1
                nonempt = True
                xDate = file[len(xCar) + 7: len(xCar) + 15]
                filterPeak(xCar,xDate,opDir,file,filtopDir,buffer = buff,whichpass = index )

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
            mainInfo = pd.read_csv(filtopDir + file[:-4] + '_info.csv')
            mainThing = (pd.read_csv(loc))
        elif index != 1:
            loc2 = filtopDir + file
            secThing = pd.read_csv(loc2)
            secInfo = pd.read_csv(filtopDir + file[:-4] + '_info.csv')
            woo = passCombine(mainThing,secThing,xCar, buff)
            mainInfo = pd.concat([secInfo,mainInfo])
            mainThing = woo.copy()
        index = index + 1
        print(file)
    mainThing = mainThing.copy().reset_index(drop = True)
    mainThing['numtimes']  = mainThing.apply(lambda x: countTimes(x.recombine,xCar),axis=1)

elif os.path.exists(finalMain):
    addingFiles = True
    toCombine = list(map(lambda x: 'FilteredPeaks_' + x[6:],toFilter))
    mainThing = pd.read_csv(finalMain)
    mainInfo = pd.read_csv(finalInfoLoc)
    curVP = mainThing.loc[mainThing.numtimes > 1, :].min_read.drop_duplicates().shape[0]
    curOP = mainThing.min_read.drop_duplicates().shape[0]
    for file in toCombine:
        loc = filtopDir + file
        if os.path.exists(loc):
            secThing = pd.read_csv(loc)
            mainThing = passCombine(mainThing,secThing,xCar, buff)
            mainInfo = pd.concat([mainInfo,pd.read_csv(filtopDir + file[:-4] + '_info.csv')])
            print(file)
    mainThing = mainThing.copy().reset_index(drop = True)
    mainThing['numtimes']  = mainThing.apply(lambda x: countTimes(x.recombine,xCar),axis=1)

mainInfo.drop_duplicates().reset_index(drop=True).FILENAME.to_csv(finalInfoLoc)
mainThing.reset_index(drop=True).to_csv(finalMain)

combined = sumData2(mainThing) ## finds locations and mean log ch4 for each peak (either verified or non yet)

## combined so only with the same overall peak
uniquePk = combined.loc[:,['min_read']].drop_duplicates()
uniqueList = combined.loc[uniquePk.index,['min_read','recombine']]
uniqueOther = combined.loc[:,['min_read','overallLON','overallLAT','mnlogCH4',
                             'verified','numtimes','minDist','maxDist']].drop_duplicates()
allTog = pd.merge(makeGPD(uniqueOther,'overallLAT','overallLON'),uniqueList,on = ['min_read'])
allTog['em'] = allTog['mnlogCH4'].swifter.apply(lambda y: estEmissions(y))
allTog['threshold'] = allTog['em'].swifter.apply(lambda x: threshold)


##### SPLITTING IF THE PEAKS WERE VERIFIED OR NOT
verTog = allTog.loc[allTog.numtimes!= 1,:]

if verTog.size > 0:
    verTog.drop(columns=['recombine']).to_file(shpFileLocName, driver="GeoJSON")
    print(f'I found {len(verTog.min_read.unique())} verified peaks')
    vpNew = len(verTog.min_read.unique())
if verTog.size ==0:
    print("Sorry, no verified peaks were found.")
    vpNew = 0
if allTog.size> 0:
    allTog.drop(columns=['recombine']).to_file(OPshpFileLocName, driver="GeoJSON")
    allTog.to_csv(allPksCSVLoc)

if allTog.size == 0:
    print("Sorry, no observed peaks were found in the given data")

if not addingFiles:
    print(f"I processed {len(toFilter)} days of driving. I analysed the data using a threshold of {100 + float(threshold)*100}% for an elevated reading, \n \
    where the threshold was calculated using the {baseLinePerc}th percentile over {backObs} observations. \n \
    I filtered the speed of the car to be between {minCarSpeed}mph and {maxCarSpeed}mph.\n \
    I created 3 summary files located here:{finRes}.\n \
    The processing took {round((time.time()-start)/60,3)} minutes. \n \
    I found {len(mainThing.min_read.unique())} observed peaks.")

elif addingFiles:
    print(f"I processed an additional {len(toFilter)} days of driving. I analysed the data using a threshold of {100 + float(threshold) * 100}% for an elevated reading, \n \
    where the threshold was calculated using the {baseLinePerc}th percentile over {backObs} observations. \n \
    I filtered the speed of the car to be between {minCarSpeed}mph and {maxCarSpeed}mph.\n \
    I created 3 summary files located here:{finRes}.\n \
    The processing took {round((time.time() - start) / 60, 3)} minutes. \n \
    I found {len(mainThing.min_read.unique()) - curOP} additional observed peaks, and {vpNew - curVP} VPs.")