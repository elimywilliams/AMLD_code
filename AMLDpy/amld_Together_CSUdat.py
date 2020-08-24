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
rawDatLoc = "/Users/emilywilliams/Documents/GitHub/CSU_SC/SC_Algorithm_Alex_45mph/SC_Algorithm_Alex_15pc_45mph_median2"
## Folder to put results in (will make subfolders later)
resFolder = "/Users/emilywilliams/Documents/GitHub/CSU_SC/SC_Algorithm_Alex_45mph/SC_Algorithm_Alex_15pc_45mph_median2/"

## CarID 
xCar = 'CSULi' # might need to be 5 letters? Need to check that!

## What Proportion above Baseline to flag as elevated (i.e. 0.1 = 10% higher)
threshold = '0.15'

## How many minutes to include in background calculation (minutes)
#timethresh = '1.7' 
time_thresh= '5'

## How many minutes to skip at the beginning of the dataset (i.e. if Colin is at his house)
initial_Time_Ignore = '0'

# minimum number of elevated readings required for an observed peak
minElevated = '2'

## Lag time for CH4 to reach sensor (in seconds)
shift = 0

## Is this an engineering file?
engineering = False

# Not super sure what timePush is but thats cool
timePush = 5 #min
timePush = 0 

###
backObs = '1020'
maxCarSpeed = '45'
minCarSpeed = '2'

baseLinePerc = '50' 



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
                            passCombine,sumData2,addOdometer
import numpy as np
import os
import pandas as pd
import os, sys, datetime, time, math, csv, numpy,gzip,shutil
from math import radians, sin, cos, sqrt, asin
from numpy import log
import geopandas as gpd
import pandas as pd #
from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
import matplotlib.pyplot as plt
import ast
from datetime import datetime
import time

### SPECIAL IDENTIFY PEAKS
########################################################################
#### IDENTIFY PEAKS  
# Input: a .csv file with processed data (already have gone through 'processRawDataEng')
# Output: saves many files, but finds elevated readings

def IdentifyPeaksCSU( xCar, xDate, xDir, xFilename,outDir,processedFileLoc,threshold = '.1',xTimeThreshold = '5.0',minElevated = '2',xB = '1020',basePerc = '50'):
    import csv, numpy    
    import geopandas as gpd
    import shutil 
    try:
        baseCalc = float(basePerc)
        xABThreshold = float(threshold)
        minElevated = float(minElevated)
        #xABThreshold = 0.1                 # above baseline threshold above the mean value
        xDistThreshold = 160.0                 # find the maximum CH4 reading of observations within street segments of this grouping distance in meters
        xSDF = 4                    # multiplier times standard deviation for floating baseline added to mean
        #xB = 1020       # the number of records that constitutes the floating baseline time -- 7200 = 1 hour (assuming average of 0.5 seconds per record)
        #xB = 102 # since it is 1 record/second
        
        xB = int(xB)
        #xB = 300 #five min?
        xTimeThreshold = float(xTimeThreshold)
        
        fn = xDir + "/" + xFilename      #set raw text file to read in
        fnOut = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".csv"       #set CSV format output for observed peaks for a given car, day, city
        fnShape = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".shp"
        fnLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".log"       #set CSV output for observed peaks for a given car, day, city
        pkLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + "_info.csv"       #set CSV output for observed peaks for a given car, day, city
        
        jsonOut =  outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".json"       #set CSV format output for observed peaks for a given car, day, city

        infOut = processedFileLoc + xCar + "_" + xDate.replace("-","") + "_info.csv"
        print(str(outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + "_info.csv"))
        
        fLog = open(fnLog, 'w')        
        shutil.copy(infOut,pkLog)

        
        #field column indices for various variables
        fDate = 0; fTime = 1; fEpochTime = 2; fNanoSeconds = 3; fLat = 4; fLon = 5; fVelocity = 6; fU = 7; fV = 8; fW = 9; fBCH4 = 10; fTCH4 = 12

        #read data in from text file and extract desired fields into a list, padding with 5 minute and hourly average
        x1 = []; x2 = []; x3 = []; x4 = []; x5 = []; x6 = []; x7 = []; x8 = []
        
        count = -1
        with open(fn, 'r') as f:
            t = csv.reader(f)
            for row in t:
                if count < 0:
                    count += 1
                    continue
                
                datet= row[fDate].replace("-","")+row[fTime].replace(":","")

                #x1.append(float(epoch)); 
                x1.append(float(str(row[fEpochTime]) + '.' + str(row[fNanoSeconds])));
                x2.append(float(int(datet)));
                x3.append(float(row[fLat]));
                x4.append(float(row[fLon])); 
                x5.append(float(row[fBCH4]));
                x6.append(float(row[fTCH4]))
                x7.append(0.0); 
                x8.append(0.0)
                
                #print (str(row[fLat])+ str(row[1]))
                count += 1
        print ("Number of observations processed: " + str(count))

        #convert lists to numpy arrays
        aEpochTime = numpy.array(x1); aDateTime = numpy.array(x2); aLat = numpy.array(x3); aLon = numpy.array(x4); aCH4 = numpy.array(x5); aTCH4 = numpy.array(x6)
        aMean = numpy.array(x7); aThreshold = numpy.array(x8)

        xLatMean = numpy.mean(aLat)
        xLonMean = numpy.mean(aLon)
        
        fLog.write ( "Day CH4_mean = " + str(numpy.mean(aCH4)) + ", Day CH4_SD = " + str(numpy.std(aCH4)) + "\n")
        fLog.write( "Center lon/lat = " + str(xLonMean) + ", " + str(xLatMean) + "\n")
        #pkLog.write('hi')
        lstCH4_AB = []

        #generate list of the index for observations that were above the threshold
        for i in range(0,count-2):
            if ((count-2)>xB):
                topBound = min((i+xB), (count-2))
                botBound = max((i-xB), 0)

                for t in range(min((i+xB), (count-2)), i, -1):
                    if aEpochTime[t] < (aEpochTime[i]+(xB/2)):
                        topBound = t
                        break
                for b in range(max((i-xB), 0), i):
                    if aEpochTime[b] > (aEpochTime[i]-(xB/2)):
                        botBound = b
                        break

                xCH4Mean = numpy.percentile(aCH4[botBound:topBound],baseCalc)
               # xCH4SD = numpy.std(aCH4[botBound:topBound])
            else:
                xCH4Mean = numpy.percentile(aCH4[0:(count-2)],baseCalc)
                #xCH4SD = numpy.std(aCH4[0:(count-2)])
            xThreshold = xCH4Mean + (xCH4Mean * xABThreshold)
            
            if (aCH4[i] > xThreshold):
                lstCH4_AB.append(i)
                aMean[i] = xCH4Mean    #insert mean + SD as upper quartile CH4 value into the array to later retreive into the peak calculation
                aThreshold[i] = xThreshold

        # now group the above baseline threshold observations into groups based on distance threshold
        lstCH4_ABP = []
        xDistPeak = 0.0
        xCH4Peak = 0.0
        xTime = 0.0
        cntPeak = 0
        cnt = 0
        sID = ""
        sPeriod5Min = ""
        prevIndex = 0
        for i in lstCH4_AB:   
            if (cnt == 0):
                xLon1 = aLon[i]; xLat1 = aLat[i]
            else:
                # calculate distance between points
                xDist = haversine(xLat1, xLon1, aLat[i], aLon[i])
                xDistPeak += xDist
                xCH4Peak += (xDist * (aCH4[i] - aMean[i]))
                xLon1 = aLon[i]; xLat1 = aLat[i]
                if (sID == ""):
                    xTime = aEpochTime[i]
                    sID = str(xCar) + "_" + str(xTime)
                    sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1))) # 30 sec
                if ((aEpochTime[i]-aEpochTime[prevIndex]) > xTimeThreshold):       #initial start of a observed peak
                    cntPeak += 1
                    xTime = aEpochTime[i]
                    xDistPeak = 0.0
                    xCH4Peak = 0.0
                    sID = str(xCar) + "_" + str(xTime)
                    sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1))) # 30 sec
                    #print str(i) +", " + str(xDist) + "," + str(cntPeak) +"," + str(xDistPeak)         
                lstCH4_ABP.append([sID, xTime, aEpochTime[i], aDateTime[i], aCH4[i], aLon[i], aLat[i], aMean[i] ,aThreshold[i], xDistPeak, xCH4Peak, aTCH4[i], sPeriod5Min])
            cnt += 1
            prevIndex = i
    
        # Finding peak_id larger than 160.0 m
        tmpsidlist = []
        for r in lstCH4_ABP:
            if (float(r[9])>160.0) and (r[0] not in tmpsidlist):
                tmpsidlist.append(r[0])
        cntPeak-=len(tmpsidlist)

        fLog.write ( "Number of peaks found: " + str(cntPeak) + "\n")
        print (xCar + "\t" + xDate + "\t" + xFilename + "\t" + str(count) + "\t" + str(len(lstCH4_ABP)))
        #### calculate attribute for the area under the curve -- PPM
        
        #write out the observed peaks to a csv to be read into a GIS
        fOut = open(fnOut, 'w')
        #s = "PEAK_NUM,EPOCHSTART,EPOCH,DATETIME,CH4,LON,LAT,CH4_BASELINE,CH4_THRESHOLD,PEAK_DIST_M,PEAK_CH4,TCH4,PERIOD5MIN\n"
        s = "OP_NUM,OP_EPOCHSTART,OB_EPOCH,OB_DATETIME,OB_CH4,OB_LON,OB_LAT,OB_CH4_BASELINE,OB_CH4_THRESHOLD,OP_PEAK_DIST_M,OP_PEAK_CH4,OB_TCH4,OB_PERIOD5MIN\n"

        fOut.write(s)

        truecount = 0
        for r in lstCH4_ABP:
            if r[0] not in tmpsidlist:
                s = ''
                for rr in r:
                    s += str(rr) + ','
                s = s[:-1]
                s += '\n'
                fOut.write(s)
                truecount += 1
        fOut.close()
        fLog.close()
        import pandas as pd
        openFile = pd.read_csv(fnOut)
        from shapely.geometry import Point
        if openFile.shape[0] != 0:
            tempCount = openFile.groupby('OP_NUM',as_index=False).OP_EPOCHSTART.count().rename(columns={'OP_EPOCHSTART':'Frequency'})
            tempCount = tempCount.loc[tempCount.Frequency>=minElevated,:]
            if tempCount.shape[0]==0:
             print("No Observed Peaks with enough Elevated Readings Found in the file: " + str(xFilename) )
            elif tempCount.shape[0]!=0:
                oFile = pd.merge(openFile,tempCount,on=['OP_NUM'])
                openFile = oFile.copy()
                del(oFile)
                openFile['minElevated'] = openFile.apply(lambda x: int(minElevated),axis=1)
                openFile.to_csv(fnOut,index=False)
                openFile['OB_CH4_AB'] = openFile.loc[:,'OB_CH4'].sub(openFile.loc[:,'OB_CH4_BASELINE'], axis = 0) 
                
                fileWt = weightedLoc(openFile,'OB_LAT','OB_LON','OP_NUM','OB_CH4_AB').loc[:,:].rename(columns = {'OB_LAT':'pk_LAT','OB_LON':'pk_LON'}).reset_index(drop = True)
                geometry_temp = [Point(xy) for xy in zip(fileWt['pk_LON'], fileWt['pk_LAT'])]
                crs = {'init': 'epsg:4326'}
                    
                    #geometry is the point of the lat/lon
                #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)
                
                ## BUFFER AROUND EACH 'OP_NUM' OF 30 M
                gdf_buff = gpd.GeoDataFrame(fileWt, crs=crs, geometry=geometry_temp)
                #gdf_buff = makeGPD(datFram,'LON','LAT')
                gdf_buff = gdf_buff.to_crs(epsg=32610)
                gdf_buff['geometry'] = gdf_buff.loc[:,'geometry'].buffer(30) 
                gdf_buff.to_file(jsonOut, driver="GeoJSON")
        elif openFile.shape[0] == 0:
            print("No Observed Peaks Found in the file: " + str(xFilename) )

        if truecount > 0:
            #arcpy.MakeXYEventLayer_management(fnOut,"LON","LAT",xCar + 
            #"L","GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',
            #SPHEROID['WGS_1984',6378137.0,298.257223563]],
            #PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];
            #-400 -400 1000000000;-100000 10000;-100000 10000;
            #8.98315284119522E-09;0.001;0.001;IsHighPrecision","#")
            #arcpy.FeatureToPoint_management(xCar + "L",fnShape,"CENTROID")
            #arcpy.Delete_management(xCar+"L")
            return True
    except ValueError:
            print ("Error in Identify Peaks")
            return False


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

#rawTexts = pd.DataFrame(os.listdir(rawDir)).loc[pd.DataFrame(os.listdir(rawDir))[0].str.endswith('.txt')]

procTexts = pd.DataFrame(os.listdir(processedFileLoc)).loc[pd.DataFrame(os.listdir(processedFileLoc))[0].str.endswith('dat.csv')]

### DONT ADD NEW FILE WITH PRE-EXISTING DATE [NEED TO WRITE CODE TO DEAL WITH THAT]
toAnalyse = []
toIdentify = []
toFilter = []

#if os.path.exists(finalInfoLoc):
#    analysed = pd.read_csv(finalInfoLoc)
#    for index,row in rawTexts.reset_index().iterrows():
#        text = rawTexts[0].iloc[index]
#        if analysed[analysed['FILENAME'].astype(str).str.contains(text)].shape[0] < 1:
#            toAnalyse.append(text)
#            toIdentify.append(s1 + '_20' + text[11:17] + '_dat.csv')
#            toFilter.append(s2 + '_20' + text[11:17] + '.csv')
#elif not os.path.exists(finalInfoLoc):
#    
    
    
for index,row in procTexts.reset_index().iterrows():
    text = procTexts[0].iloc[index]
    toAnalyse.append(text)
    #toIdentify.append(text[:-7] + '_dat.csv')
    toIdentify.append(text)
    toFilter.append(s2 + text[5:-8] + '.csv')
    

##### START OF THE ALGORITHM
start = time.time()

# =============================================================================
# if __name__ == '__main__':  
# 
#     dateList = []
#     x1=""
#     count = 0
#     #listthing = os.listdir(rawDir)
#     listthing = toAnalyse ## changing to include the files we have to analyse?
#     for file in listthing:
#         #print (file)
#         if file.endswith(".txt"):
#             dateF = file[11:17]
#             if dateF in set(dateList):
#                 bFirst = False
#             
#             if dateF not in set(dateList):
#                 bFirst = True
#                 dateList.append(dateF)
# 
#             
#             #if file[:10] == x1:
#              #   bFirst = False
#             #else:
#              #   bFirst = True
#             x1 = file[:10]
#             
#             #xCar = "CSULi"
#             xDate = file[:10]
#             #theResult = ProcessRawData(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc)
#             if engineering:
#                 theResult = ProcessRawDataEng(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc,initialTimeIgnore,shift,maxCarSpeed,minCarSpeed)
#             elif not engineering:
#                 theResult = ProcessRawData(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc)
#             #del(bFirst)
#             count = count + 1
# =============================================================================


if __name__ == '__main__':
    #listthing = os.listdir(processedFileLoc).copy() 
    listthing = toIdentify.copy()
    for file in listthing:
        if file.startswith(s1) and file.endswith("dat.csv"):
            xDate = file[6:14]
            theResult = IdentifyPeaksCSU(xCar, xDate, processedFileLoc, file,opDir,processedFileLoc,threshold,time_thresh,minElevated,backObs)

if __name__ == '__main__':
    index = 0
    numproc = 0
    listthing = os.listdir(opDir)
    for file in listthing:
        if file.startswith(s2) and file.endswith('.csv') and not file.endswith('info.csv'):
            print(file)
            file_loc = opDir + file
            csv_loc = opDir  + file[:-4]+ '_info.csv'
            nonempt = False
            if pd.read_csv(file_loc).size != 0:
                index += 1
                nonempt = True
                xDate = file[12:20]
                filterPeak(xCar,xDate,opDir,file,filtopDir,whichpass = index )
        
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
            woo = passCombine(mainThing,secThing)
            mainThing = woo.copy()
        index = index + 1
        print(file)
    mainThing['numtimes']  = mainThing.apply(lambda x: countTimes(x.recombine),axis=1)

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
mainThing.to_csv(finalMain)
        
print("I processed "+ str(index) + ' days of driving. The processed files are now stored in the folder: ' + str(filtopDir))

################

#combined = summarizeDat(mainThing) ## finds locations and mean log ch4 for each peak (either verified or non yet)
combined = sumData2(mainThing) ## finds locations and mean log ch4 for each peak (either verified or non yet)


## combined so only with the same overall peak
uniquePk = combined.loc[:,['min_read']].drop_duplicates()
uniqueList = combined.loc[uniquePk.index,['min_read','recombine']]
uniqueOther = combined.loc[:,['min_read','overallLON','overallLAT','mnlogCH4',\
                             'verified','numtimes']].drop_duplicates()

unique_gdf = makeGPD(uniqueOther,'overallLAT','overallLON')
unique_gdf2 = pd.merge(unique_gdf,uniqueList,on = ['min_read'])

#unique_gdf = makeGPD(combined.loc[:,['min_read','pk_LAT','pk_LON']],'pk_LAT','pk_LON')
#combinedGeo = unique_gdf.dissolve(by='min_read',as_index = False)


#allTog = pd.merge(combinedGeo,uniqueOther,on=['min_read'])
#allTog = pd.merge(allTog,uniqueList,on=['min_read'])
allTog = unique_gdf2.copy()
allTog['em'] = allTog.apply(lambda y: estEmissions(y['mnlogCH4']),axis=1)
allTog['threshold'] = allTog.apply(lambda x: threshold,axis =1 )


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
#print(end - start)