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
xCar = 'SCCar' # might need to be 5 letters? need to check into that

## Folder with .txt Data
rawDatLoc = "/Users/emilywilliams/Documents/DrivingData/ColDatEngShort" 

## Folder to put results in (will make subfolders later)
resFolder = "/Users/emilywilliams/Documents/DrivingData/ColDatEngShort/"
initialTimeIgnore = '5'
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

threshold = '0.1'
timethresh = '1.7' ## time to include in background calculation (minutes)
##################################
### IMPORTING NECESSARY MODULES
##################################

#import contextily as ctx
import pandas as pd
import os, sys, datetime, time, math, csv, numpy,gzip,shutil
from math import radians, sin, cos, sqrt, asin
from numpy import log
import geopandas as gpd
import pandas as pd #
from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
import matplotlib.pyplot as plt
from datetime import datetime
sys.path.insert(0, functionFileLoc) ##change this to location of the "allFunctions.py" file
from allFunctionsSC import IdentifyPeaks,filterPeak,unique,unIfInt,IsInPK,intersect,\
                            passCombine, weightedLoc,verPk,estEmissions,haversine,\
                           ProcessRawData,wt_time_Locs,sumthing,makeGEO, makeGPD, summarizeDat
from allFunctionsSC_engineering import ProcessRawDataEng,getQuad,calcTheta,calcBearing


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
            
            theResult = ProcessRawData(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc)
            
            #del(bFirst)
            count = count + 1


if __name__ == '__main__':
    #listthing = os.listdir(processedFileLoc).copy() 
    listthing = toIdentify.copy()
    for file in listthing:
        if file.startswith(s1) and file.endswith("dat.csv"):
            xDate = file[6:14]
            theResult = IdentifyPeaks(xCar, xDate, processedFileLoc, file,opDir,processedFileLoc,threshold,timethresh)

if not os.path.exists(finalMain):
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
elif os.path.exists(finalMain):
    index = 0
    numproc = 0
    listthing = toFilter
    main = pd.read_csv(finalMain)
    mainGeo = main.loc[:,'geometry']
    mainThing = makeGPD(main.drop(columns = ['geometry']).copy(),'pk_LAT','pk_LON')
    
    mainInfo = pd.read_csv(finalInfoLoc)
    for file in listthing:
        if file.startswith(s2) and file.endswith('.csv') and not file.endswith('info.csv'):
            file_loc = opDir + file
            csv_loc = opDir + "/" + file[:-4] + '_info.csv'
            nonempt = False
            if pd.read_csv(file_loc).size != 0:
                index += 1
                nonempt = True
            xDate = file[12:20]
            if  nonempt:
                secondThing = filterPeak(xCar,xDate,opDir,file,filtopDir,whichpass = index )
                mainThing = passCombine(mainThing,secondThing)
                tempInfo = pd.read_csv(csv_loc)
                mainInfo = pd.concat([mainInfo, tempInfo], axis=0)
                
          
mainInfo.reset_index().FILENAME.to_csv(finalInfoLoc)
mainThing.to_csv(finalMain)
        
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
if verTog.size > 0:
    verTog.drop(columns=['recombine']).to_file(shpFileLocName, driver="GeoJSON")
if verTog.size ==0:
    print("Sorry, no verified peaks were found.")  
if allTog.size> 0:
    allTog.drop(columns=['recombine']).to_file(OPshpFileLocName, driver="GeoJSON")
    allTog.to_csv(allPksCSVLoc)

if allTog.size == 0:
    print("Sorry, no observed peaks were found in the given data")
    
end = time.time()
print("I created three summary files located here: " + str(finRes) + ". The processing took " + str(round((end-start)/60,3)) + str(" minutes."))
