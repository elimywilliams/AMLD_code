#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
## all functions colin
Created on Mon Apr  6 14:51:23 2020

@author: emilywilliams
"""

# /Users/emilywilliams/Documents/mixedmods
import pandas as pd

#xDir =  "/Users/emilywilliams/Documents/mixedmods/"
#outFolder = "/Users/emilywilliams/Documents/mixedmods/"
#xCar = 'ColLi'
#xFilename = 'col_redo.csv'

#fileloc = "/Users/emilywilliams/Desktop//COL_drive1.txt"
#outFolder = "/Users/emilywilliams/Desktop/COl_Results/"

#infile = pd.read_csv(r"/Users/emilywilliams/Desktop//COL_drive1.txt")
#infile.to_csv(r"/Users/emilywilliams/Desktop//COL_drive1.csv",index=None)
#fileloc = "/Users/emilywilliams/Desktop//COL_drive1.csv"


#xDir = "/Users/emilywilliams/Desktop/"
#xDate = "04012020"
#xCar = "COLLi"
#xFilename= "COL_drive1.csv"

def haversine(lat1, lon1, lat2, lon2, radius=6371): # 6372.8 = earth radius in kilometers
    from math import radians, sin, cos, sqrt, asin

    dLat = radians(lat2 - lat1)
    dLon = radians(lon2 - lon1)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    c = 2*asin(sqrt(sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2))

    return radius*c*1000 #return in meters


def ProcessRawData( xCar, xDate, xDir, xFilename, bFirst, gZIP, xOut):
    from datetime import datetime
    import os
    try:
        
        #xOutDir = xDir 
        
        xMaxCarSpeed = 45.0 / 2.23694     # assumes 45 mph as max, convert to meters per second 
        xMinCarSpeed = 2.0 / 2.23694      # minimum 2 miles per hour
        xMinRSSI = 50  #if RSSI is below this we don't like it
        

        # reading in the data with specific headers
        #          0     1    2    3       4           5    6       7        8        9          10                 11              12           13            14      15      16      17        18         19         20         21         22         23        24   25  26       27           28       29           30       31       32       33  34        35   36   37  38   39       40       41   42       43   44   45   46   47   48   49   50   51     52     53     54
        sHeader = "Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude"
        sOutHeader = "DATE,TIME,SECONDS,NANOSECONDS,VELOCITY,U,V,W,BCH4,BRSSI,TCH4,TRSSI,PRESS_MBAR,INLET,TEMPC,CH4,H20,C2H6,R,C2C1,BATTV,POWMV,CURRMA,SOCPER,LAT,LONG\n"

        # somehow gZIP is indicating if  it is the first file name (I think if it is 0 then it is the first file)
        if gZIP == 0:
            f = gzip.open(xDir + "/" + xFilename, 'r') #if in python 3, change this to "r" or just "b" can't remember but something about a bit not a string
        else:
            f = open(xDir + "/" + xFilename, 'r')
        
        # process    
        #if first time on this car/date, then write header out
        
        xdat = str('20') + xFilename[11:17]
        
        #fnOut = xOutDir + xCar + "_" + xDate.replace("-", "") + "_dat.csv"       #set CSV output for raw data
        #fnLog = xOutDir + xCar + "_" + xDate.replace("-", "") + "_log.csv"       #output for logfile
        
        fnOut = xOut  + xCar + "_" + xdat + "_dat.csv"       #set CSV output for raw data
        fnLog = xOut  + xCar + "_" + xdat + "_log.csv"       #output for logfile
       
        if bFirst:
            fOut = open(fnOut, 'w')
            fOut.write(sOutHeader)
            fLog = open(fnLog, 'w')
            print ("fnLog: "+fnOut) 
        if not bFirst:
            fOut = open(fnOut, 'a')
            fLog = open(fnLog, 'a')
            print('hi')
        
        #read all lines
        xCntObs = -1
        xCntGoodValues = 0
        for row in f:
            #print(row)
            woo = row
            bGood = True
            if xCntObs < 0:
                bGood = False
                xCntObs += 1
            if bGood:
                lstS = row.split(",")
                dtime = lstS[0]
                dateob = datetime(int(dtime[6:10]),int(dtime[0:2]),int(dtime[3:5]),int(dtime[11:13]),int(dtime[14:16]),int(dtime[17:19]),int(float(dtime[19:23])*1000000))
                epoch = dateob.strftime('%s.%f')
                dtime = int(dateob.strftime('%Y%m%d%H%M%S'))
  
                # change this once we have QA/QC stuff
                
#                # if RSSI of bottome sensor is below 50
#                if float(lstS[28]) < xMinRSSI:
#                    fLog.write("RSSI (Bottom) value less than 50: "+ str(lstS[28]) + "\n")
#                    continue
#                # Car Speed
#                if float(lstS[12]) > xMaxCarSpeed:
#                    fLog.write("Car speed of " + str(float(lstS[12])) + " exceeds max threshold of: " + str(xMaxCarSpeed) + "\n")
#                    continue
#                if float(lstS[12]) < xMinCarSpeed:
#                    fLog.write("Car speed of " + str(float(lstS[12])) + " less than min threshold of: " + str(xMinCarSpeed) + "\n")
#                    continue

                # For some reason it is producing its longitude in positive number while USA is located at negative longitude
                # thats why we do -1 * float(lstS[7])
                
                # fix this when we have stuffs
                
#                s1 = str(lstS[1])+","+str(lstS[2])+","+str(lstS[3])+","+str(lstS[4])+","+str(lstS[6])+","
#                s1 += str(-1 * float(lstS[7]))+","+str(lstS[12])+","+str(lstS[14])+","+str(lstS[15])+","+str(lstS[16])+","+str(lstS[25])+","
#                s1 += str(lstS[28])+","+str(lstS[38])+","+str(lstS[41])+"\n"
                
                ## choosing what to write in the .csv


                csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(dateob.strftime('%H:%M:%S'))  + ',' + str(int(pd.to_numeric(dateob.strftime('%s.%f')))) + ',' + str(pd.to_numeric(dateob.strftime('%f')) *1000) + str(',')
                csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(lstS[4]) + ',' + str('0') + ','+  str(lstS[4]) + ','
                csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(lstS[4]) + ',' + str(lstS[5]) +',' +  str(lstS[6]) + ','
                csvWrite += str(lstS[7]) + ',' + str(lstS[8]) + ',' + str(lstS[9]) + ',' + str(lstS[10]) + ','+ str(lstS[11]) + ',' + str(lstS[12]) + ',' + str(lstS[13]) + str(',') + str(lstS[14]) 
                fOut.write(csvWrite)
#                xCntGoodValues += 1

            xCntObs += 1

        #sOut = str(gZIP) + "," + str(f) + "," + str(xCntObs) + "," + str(xCntGoodValues) + "\n"
        #fLog.write(sOut)
               
        fOut.close()
        fLog.close()
        
        #xDate = dateob.strftime("%Y%m%d")
        
        #newfnOut = xOutDir + xCar + "_" + xDate + "_dat.csv"       #set CSV output for raw data
        #newfnLog = xOutDir + xCar + "_" + xDate + "_log.csv"  
        
       # os.rename(str(fnOut),str(newfnOut))
       # os.rename(str(fnLog),str(newfnLog))


        #print (xCar + "\t" + xDate + "\t" + xCar + '_' + xDate + '.csv' + "\t" + str(xCntObs) + "\t" + str(xCntGoodValues) + "\t" + str(gZIP))
        print (xCar + "\t" + xdat + "\t" + fnOut[-22:] + "\t" + str(xCntObs) + "\t" + str(xCntGoodValues) + "\t" + str(gZIP))

        return True
    except ValueError:
        return False
    
    
## POTENTIALLY GOOD ONE    
def IdentifyPeaks( xCar, xDate, xDir, xFilename):
    import contextily as ctx
    import os, sys, datetime, time, math, csv, numpy,gzip
    from math import radians, sin, cos, sqrt, asin
    import geopandas as gpd
    import pandas as pd #
    from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
    import matplotlib.pyplot as plt
    from datetime import datetime
    
    try:

        xABThreshold = 0.1                 # above baseline threshold above the mean value
        xDistThreshold = 160.0                 # find the maximum CH4 reading of observations within street segments of this grouping distance in meters
        xSDF = 4                    # multiplier times standard deviation for floating baseline added to mean
        xB = 1020       # the number of records that constitutes the floating baseline time -- 7200 = 1 hour (assuming average of 0.5 seconds per record)
        #xB = 300 #five min?
        xTimeThreshold = 5.0
        
        fn = xDir + "/" + xFilename      #set raw text file to read in
        #fn = fileloc
        fnOut = xDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".csv"       #set CSV format output for observed peaks for a given car, day, city
        fnShape = xDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".shp"
        fnLog = xDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".log"       #set CSV output for observed peaks for a given car, day, city

        fLog = open(fnLog, 'w')

        #field column indices for various variables
        fDate = 0; 
        fTime = 1; 
        fEpochTime = 2; 
        fNanoSeconds = 3; 
        fVelocity = 4; 
        fU = 5; 
        fV = 6; 
        fW = 7;
        fBCH4 = 10;#fBCH4 = 8; 
        #fBRSSI = 9;
        fTCH4 = 10;
        TRSSI = 11;
        PRESS = 12;
        INLET = 13;
        TEMP = 14;
        CH4 = 15;
        H20 = 16;
        C2H6 = 17;
        R = 18;
        C2C1 = 19;
        BATT = 20;
        POWER = 21;
        CURR = 22;
        SOCPER = 23;
        fLat = 24; 
        fLon = 25; 

        
        #fDateTime = 0; 
        #fTime = 16; 
        #fEpochTime = 3;  fLat = 13; fLon = 14; fBCH4 = 4; fTCH4 =4
        #fNanoSeconds = 18;
        #fVelocity = 6;
        
        #read data in from text file and extract desired fields into a list, padding with 5 minute and hourly average
        x1 = []; x2 = []; x3 = []; x4 = []; x5 = []; x6 = []; x7 = []; x8 = []
        
        
    
        count = -1
        with open(fn, 'r') as f:
            t = csv.reader(f)
            for row in t:
                if count < 0:
                    count += 1
                    continue
                #woo = row
                #epoch = float(row[fEpochTime]+"."+row[fNanoSeconds][0])
                #print(row[fLat])
                #dtime = row[fTime]
                #dateob = datetime(int(dtime[6:10]),int(dtime[0:2]),int(dtime[3:5]),int(dtime[11:13]),int(dtime[14:16]),int(dtime[17:19]),int(float(dtime[19:23])*1000000))
                #epoch = dateob.strftime('%s.%f')
                #dtime = int(dateob.strftime('%Y%m%d%H%M%S'))
  
                
                #epoch = float(row[fEpochTime])
                #woo = row
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

                xCH4Mean = numpy.percentile(aCH4[botBound:topBound],50)
                xCH4SD = numpy.std(aCH4[botBound:topBound])
            else:
                xCH4Mean = numpy.percentile(aCH4[0:(count-2)],50)
                xCH4SD = numpy.std(aCH4[0:(count-2)])
            xThreshold = xCH4Mean + (xCH4Mean * 0.1)
            
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
        s = "PEAK_NUM,EPOCHSTART,EPOCH,DATETIME,CH4,LON,LAT,CH4_BASELINE,CH4_THRESHOLD,PEAK_DIST_M,PEAK_CH4,TCH4,PERIOD5MIN\n"
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

        if truecount > 0:
            #arcpy.MakeXYEventLayer_management(fnOut,"LON","LAT",xCar + "L","GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;0.001;0.001;IsHighPrecision","#")
            #arcpy.FeatureToPoint_management(xCar + "L",fnShape,"CENTROID")
            #arcpy.Delete_management(xCar+"L")
            return True
    except ValueError:
            print ("Error in Identify Peaks")
            return False

        
        
#def IdentifyPeaks( xCar, xDate, xDir, xFilename):
#    try:
#
#        xABThreshold = 0.1                 # above baseline threshold above the mean value
#        xDistThreshold = 160.0                 # find the maximum CH4 reading of observations within street segments of this grouping distance in meters
#        xSDF = 4                    # multiplier times standard deviation for floating baseline added to mean
#        xB = 1020       # the number of records that constitutes the floating baseline time -- 7200 = 1 hour (assuming average of 0.5 seconds per record)
#        #xB = 300 #five min?
#        xTimeThreshold = 5.0
#        
#        fn = xDir + "/" + xFilename      #set raw text file to read in
#        #fn = fileloc
#        fnOut = xDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".csv"       #set CSV format output for observed peaks for a given car, day, city
#        fnShape = xDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".shp"
#        fnLog = xDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".log"       #set CSV output for observed peaks for a given car, day, city
#
#        fLog = open(fnLog, 'w')
#
#        #field column indices for various variables
#        #fDate = 0; fTime = 1; fEpochTime = 2; fNanoSeconds = 3; fLat = 4; fLon = 5; fVelocity = 6; fU = 7; fV = 8; fW = 9; fBCH4 = 10; fTCH4 = 12
#        
#        fDateTime = 0; 
#        #fTime = 16; 
#        fEpochTime = 3;  fLat = 13; fLon = 14; fBCH4 = 4; fTCH4 =4
#        #fNanoSeconds = 18;
#        #fVelocity = 6;
#        
#        #read data in from text file and extract desired fields into a list, padding with 5 minute and hourly average
#        x1 = []; x2 = []; x3 = []; x4 = []; x5 = []; x6 = []; x7 = []; x8 = []
#        
#        
#    
#        count = -1
#        with open(fn, 'r') as f:
#            t = csv.reader(f)
#            for row in t:
#                if count < 0:
#                    count += 1
#                    continue
#                #epoch = float(row[fEpochTime]+"."+row[fNanoSeconds][0])
#                #print(row[fLat])
#                dtime = row[fDateTime]
#                dateob = datetime(int(dtime[6:10]),int(dtime[0:2]),int(dtime[3:5]),int(dtime[11:13]),int(dtime[14:16]),int(dtime[17:19]),int(float(dtime[19:23])*1000000))
#                epoch = dateob.strftime('%s.%f')
#                dtime = int(dateob.strftime('%Y%m%d%H%M%S'))
#  
#                
#                #epoch = float(row[fEpochTime])
#                woo = row
#                #print(row[1])
#
#                x1.append(float(epoch)); x2.append(float(dtime));
#                x3.append(float(row[fLat]));
#                x4.append(float(row[fLon])); 
#                x5.append(float(row[fBCH4]));
#                x6.append(float(row[fTCH4]))
#                x7.append(0.0); x8.append(0.0)
#                
#                count += 1
#        print ("Number of observations processed: " + str(count))
#
#        #convert lists to numpy arrays
#        aEpochTime = numpy.array(x1); aDateTime = numpy.array(x2); aLat = numpy.array(x3); aLon = numpy.array(x4); aCH4 = numpy.array(x5); aTCH4 = numpy.array(x6)
#        aMean = numpy.array(x7); aThreshold = numpy.array(x8)
#
#        xLatMean = numpy.mean(aLat)
#        xLonMean = numpy.mean(aLon)
#        
#        fLog.write ( "Day CH4_mean = " + str(numpy.mean(aCH4)) + ", Day CH4_SD = " + str(numpy.std(aCH4)) + "\n")
#        fLog.write( "Center lon/lat = " + str(xLonMean) + ", " + str(xLatMean) + "\n")
#        
#        lstCH4_AB = []
#
#        #generate list of the index for observations that were above the threshold
#        for i in range(0,count-2):
#            if ((count-2)>xB):
#                topBound = min((i+xB), (count-2))
#                botBound = max((i-xB), 0)
#
#                for t in range(min((i+xB), (count-2)), i, -1):
#                    if aEpochTime[t] < (aEpochTime[i]+(xB/2)):
#                        topBound = t
#                        break
#                for b in range(max((i-xB), 0), i):
#                    if aEpochTime[b] > (aEpochTime[i]-(xB/2)):
#                        botBound = b
#                        break
#
#                xCH4Mean = numpy.percentile(aCH4[botBound:topBound],50)
#                xCH4SD = numpy.std(aCH4[botBound:topBound])
#            else:
#                xCH4Mean = numpy.percentile(aCH4[0:(count-2)],50)
#                xCH4SD = numpy.std(aCH4[0:(count-2)])
#            xThreshold = xCH4Mean + (xCH4Mean * 0.1)
#            
#            if (aCH4[i] > xThreshold):
#                lstCH4_AB.append(i)
#                aMean[i] = xCH4Mean    #insert mean + SD as upper quartile CH4 value into the array to later retreive into the peak calculation
#                aThreshold[i] = xThreshold
#
#        # now group the above baseline threshold observations into groups based on distance threshold
#        lstCH4_ABP = []
#        xDistPeak = 0.0
#        xCH4Peak = 0.0
#        xTime = 0.0
#        cntPeak = 0
#        cnt = 0
#        sID = ""
#        sPeriod5Min = ""
#        prevIndex = 0
#        for i in lstCH4_AB:   
#            if (cnt == 0):
#                xLon1 = aLon[i]; xLat1 = aLat[i]
#            else:
#                # calculate distance between points
#                xDist = haversine(xLat1, xLon1, aLat[i], aLon[i])
#                xDistPeak += xDist
#                xCH4Peak += (xDist * (aCH4[i] - aMean[i]))
#                xLon1 = aLon[i]; xLat1 = aLat[i]
#                if (sID == ""):
#                    xTime = aEpochTime[i]
#                    sID = str(xCar) + "_" + str(xTime)
#                    sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1))) # 30 sec
#                if ((aEpochTime[i]-aEpochTime[prevIndex]) > xTimeThreshold):       #initial start of a observed peak
#                    cntPeak += 1
#                    xTime = aEpochTime[i]
#                    xDistPeak = 0.0
#                    xCH4Peak = 0.0
#                    sID = str(xCar) + "_" + str(xTime)
#                    sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1))) # 30 sec
#                    #print str(i) +", " + str(xDist) + "," + str(cntPeak) +"," + str(xDistPeak)         
#                lstCH4_ABP.append([sID, xTime, aEpochTime[i], aDateTime[i], aCH4[i], aLon[i], aLat[i], aMean[i] ,aThreshold[i], xDistPeak, xCH4Peak, aTCH4[i], sPeriod5Min])
#            cnt += 1
#            prevIndex = i
#    
#        # Finding peak_id larger than 160.0 m
#        tmpsidlist = []
#        for r in lstCH4_ABP:
#            if (float(r[9])>160.0) and (r[0] not in tmpsidlist):
#                tmpsidlist.append(r[0])
#        cntPeak-=len(tmpsidlist)
#
#        fLog.write ( "Number of peaks found: " + str(cntPeak) + "\n")
#        print (xCar + "\t" + xDate + "\t" + xFilename + "\t" + str(count) + "\t" + str(len(lstCH4_ABP)))
#        #### calculate attribute for the area under the curve -- PPM
#        
#        #write out the observed peaks to a csv to be read into a GIS
#        fOut = open(fnOut, 'w')
#        s = "PEAK_NUM,EPOCHSTART,EPOCH,DATETIME,CH4,LON,LAT,CH4_BASELINE,CH4_THRESHOLD,PEAK_DIST_M,PEAK_CH4,TCH4,PERIOD5MIN\n"
#        fOut.write(s)
#
#        truecount = 0
#        for r in lstCH4_ABP:
#            if r[0] not in tmpsidlist:
#                s = ''
#                for rr in r:
#                    s += str(rr) + ','
#                s = s[:-1]
#                s += '\n'
#                fOut.write(s)
#                truecount += 1
#        fOut.close()
#        fLog.close()
#
#        if truecount > 0:
#            #arcpy.MakeXYEventLayer_management(fnOut,"LON","LAT",xCar + "L","GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;0.001;0.001;IsHighPrecision","#")
#            #arcpy.FeatureToPoint_management(xCar + "L",fnShape,"CENTROID")
#            #arcpy.Delete_management(xCar+"L")
#            return True
#    except ValueError:
#            print ("Error in Identify Peaks")
#            return False
            
        
#def IdentifyPeaks( xCar, xDate, xDir, xFilename):
#    try:
#
#        xABThreshold = 0.1                 # above baseline threshold above the mean value
#        xDistThreshold = 160.0                 # find the maximum CH4 reading of observations within street segments of this grouping distance in meters
#        xSDF = 4                    # multiplier times standard deviation for floating baseline added to mean
#        xB = 1020       # the number of records that constitutes the floating baseline time -- 7200 = 1 hour (assuming average of 0.5 seconds per record)
#        #xB = 300 #five min?
#        xTimeThreshold = 5.0
#        
#        fn = xDir + "/" + xFilename      #set raw text file to read in
#        fn = fileloc
#        fnOut = xDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".csv"       #set CSV format output for observed peaks for a given car, day, city
#        fnShape = xDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".shp"
#        fnLog = xDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".log"       #set CSV output for observed peaks for a given car, day, city
#
#        fLog = open(fnLog, 'w')
#
#        #field column indices for various variables
#        #fDate = 0; fTime = 1; fEpochTime = 2; fNanoSeconds = 3; fLat = 4; fLon = 5; fVelocity = 6; fU = 7; fV = 8; fW = 9; fBCH4 = 10; fTCH4 = 12
#        
#        fDateTime = 4; 
#        #fTime = 16; 
#        fEpochTime = 3;  fLat = 5; fLon = 6; fBCH4 = 7; fTCH4 =7
#        #fNanoSeconds = 18;
#        #fVelocity = 6;
#        
#        #read data in from text file and extract desired fields into a list, padding with 5 minute and hourly average
#        x1 = []; 
#        x2 = []; 
#        x3 = []; 
#        x4 = []; 
#        x5 = []; 
#        x6 = []; 
#        x7 = []; 
#        x8 = [];
#        x9 = []; #inlet number
#        x10 = []; #pressure
#        x11 = []; #temp
#        
#    
#    
#        count = -1
#        with open(fn, 'r') as f:
#            
#            t = f.readlines()
#            #t = csv.reader(f)
#            for row in t:
#                if count < 0:
#                    count += 1
#                    continue
#                #epoch = float(row[fEpochTime]+"."+row[fNanoSeconds][0])
#                #print(row[fLat])
#                
#                #epoch = float(row[fEpochTime])
#               # print(row[0:1])
#                #date = (row[6:10] + row[0:2]+row[3:5])
#                time = str((row[10:19].replace(":","")))
#                woo = row
#
#                dateob = datetime(int(row[6:10]),int(row[0:2]),int(row[3:5]),int(row[11:13]),int(row[14:16]),int(row[17:19]),int(float(row[19:23])*1000000))
#                epoch = dateob.strftime('%s.%f')
#                dtime = dateob.strftime('%Y%m%d%H%M%S')
#                #print(month)
#                #print(row[1])
#                #print(row[0:6])
#                if count < 10:
#                    date = str(row[6:10] + row[0:2]+row[3:5])
#                    #print((row[11:20]))
#                    time = str((row[10:19].replace(":","")))
#                    datetime = (date + time).replace(" ","")
#                    print(datetime)
#                    #print((row[0:10].replace("/","")))
#
#                #datetime = row[fDate].replace("-","")+row[fTime].replace(":","")
#                x1.append(float(epoch));
#                x2.append(float(dtime));
#                x3.append(float(row[fLat]));
#                x4.append(float(row[fLon])); 
#                x5.append(float(row[fBCH4]));
#                x6.append(float(row[fTCH4]))
#                x7.append(0.0); x8.append(0.0)
#                x8.append()
#                x9.append(float(row[24:25]))
#                x10.append(float(row[26:33]))
#                x11.append(float(row[34:41]))
##                
#                count += 1
#        print ("Number of observations processed: " + str(count))
#
#        #convert lists to numpy arrays
#        aEpochTime = numpy.array(x1); aDateTime = numpy.array(x2); aLat = numpy.array(x3); aLon = numpy.array(x4); aCH4 = numpy.array(x5); aTCH4 = numpy.array(x6)
#        aMean = numpy.array(x7); aThreshold = numpy.array(x8)
#
#        xLatMean = numpy.mean(aLat)
#        xLonMean = numpy.mean(aLon)
#        
#        fLog.write ( "Day CH4_mean = " + str(numpy.mean(aCH4)) + ", Day CH4_SD = " + str(numpy.std(aCH4)) + "\n")
#        fLog.write( "Center lon/lat = " + str(xLonMean) + ", " + str(xLatMean) + "\n")
#        
#        lstCH4_AB = []
#
#        #generate list of the index for observations that were above the threshold
#        for i in range(0,count-2):
#            if ((count-2)>xB):
#                topBound = min((i+xB), (count-2))
#                botBound = max((i-xB), 0)
#
#                for t in range(min((i+xB), (count-2)), i, -1):
#                    if aEpochTime[t] < (aEpochTime[i]+(xB/2)):
#                        topBound = t
#                        break
#                for b in range(max((i-xB), 0), i):
#                    if aEpochTime[b] > (aEpochTime[i]-(xB/2)):
#                        botBound = b
#                        break
#
#                xCH4Mean = numpy.percentile(aCH4[botBound:topBound],50)
#                xCH4SD = numpy.std(aCH4[botBound:topBound])
#            else:
#                xCH4Mean = numpy.percentile(aCH4[0:(count-2)],50)
#                xCH4SD = numpy.std(aCH4[0:(count-2)])
#            xThreshold = xCH4Mean + (xCH4Mean * 0.1)
#            
#            if (aCH4[i] > xThreshold):
#                lstCH4_AB.append(i)
#                aMean[i] = xCH4Mean    #insert mean + SD as upper quartile CH4 value into the array to later retreive into the peak calculation
#                aThreshold[i] = xThreshold
#
#        # now group the above baseline threshold observations into groups based on distance threshold
#        lstCH4_ABP = []
#        xDistPeak = 0.0
#        xCH4Peak = 0.0
#        xTime = 0.0
#        cntPeak = 0
#        cnt = 0
#        sID = ""
#        sPeriod5Min = ""
#        prevIndex = 0
#        for i in lstCH4_AB:   
#            if (cnt == 0):
#                xLon1 = aLon[i]; xLat1 = aLat[i]
#            else:
#                # calculate distance between points
#                xDist = haversine(xLat1, xLon1, aLat[i], aLon[i])
#                xDistPeak += xDist
#                xCH4Peak += (xDist * (aCH4[i] - aMean[i]))
#                xLon1 = aLon[i]; xLat1 = aLat[i]
#                if (sID == ""):
#                    xTime = aEpochTime[i]
#                    sID = str(xCar) + "_" + str(xTime)
#                    sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1))) # 30 sec
#                if ((aEpochTime[i]-aEpochTime[prevIndex]) > xTimeThreshold):       #initial start of a observed peak
#                    cntPeak += 1
#                    xTime = aEpochTime[i]
#                    xDistPeak = 0.0
#                    xCH4Peak = 0.0
#                    sID = str(xCar) + "_" + str(xTime)
#                    sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1))) # 30 sec
#                    #print str(i) +", " + str(xDist) + "," + str(cntPeak) +"," + str(xDistPeak)         
#                lstCH4_ABP.append([sID, xTime, aEpochTime[i], aDateTime[i], aCH4[i], aLon[i], aLat[i], aMean[i] ,aThreshold[i], xDistPeak, xCH4Peak, aTCH4[i], sPeriod5Min])
#            cnt += 1
#            prevIndex = i
#    
#        # Finding peak_id larger than 160.0 m
#        tmpsidlist = []
#        for r in lstCH4_ABP:
#            if (float(r[9])>160.0) and (r[0] not in tmpsidlist):
#                tmpsidlist.append(r[0])
#        cntPeak-=len(tmpsidlist)
#
#        fLog.write ( "Number of peaks found: " + str(cntPeak) + "\n")
#        print (xCar + "\t" + xDate + "\t" + xFilename + "\t" + str(count) + "\t" + str(len(lstCH4_ABP)))
#        #### calculate attribute for the area under the curve -- PPM
#        
#        #write out the observed peaks to a csv to be read into a GIS
#        fOut = open(fnOut, 'w')
#        s = "PEAK_NUM,EPOCHSTART,EPOCH,DATETIME,CH4,LON,LAT,CH4_BASELINE,CH4_THRESHOLD,PEAK_DIST_M,PEAK_CH4,TCH4,PERIOD5MIN\n"
#        fOut.write(s)
#
#        truecount = 0
#        for r in lstCH4_ABP:
#            if r[0] not in tmpsidlist:
#                s = ''
#                for rr in r:
#                    s += str(rr) + ','
#                s = s[:-1]
#                s += '\n'
#                fOut.write(s)
#                truecount += 1
#        fOut.close()
#        fLog.close()
#
#        if truecount > 0:
#            #arcpy.MakeXYEventLayer_management(fnOut,"LON","LAT",xCar + "L","GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;0.001;0.001;IsHighPrecision","#")
#            #arcpy.FeatureToPoint_management(xCar + "L",fnShape,"CENTROID")
#            #arcpy.Delete_management(xCar+"L")
#            return True
#    except ValueError:
#            print ("Error in Identify Peaks")
#            return False
#
def filterPeak(xCar,xDate,xDir,xFilename, outFolder,whichpass = 0):
    import pandas as pd #
    import geopandas as gpd
    from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
#
    file_loc = xDir + xFilename
    new_loc = outFolder + "Filtered" + xFilename
    new_loc_json = new_loc[:-3] + 'json'

    datFram = pd.read_csv(file_loc)
    
    datFram_cent =  datFram.copy()
    datFram_cent['CH4_AB'] = datFram_cent['CH4'].sub(datFram_cent['CH4_BASELINE'], axis = 0) 
    maxch4 = datFram_cent.groupby('PEAK_NUM',as_index = False).CH4_AB.max().rename(columns = {'CH4_AB':'pk_maxCH4_AB'})
    datFram_cent_wmax = pd.merge(datFram_cent,maxch4,on = ['PEAK_NUM'])  
    datFram_wtLoc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB').rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'})
    datFram_wtLocMax = pd.merge(datFram_wtLoc,maxch4,on = ['PEAK_NUM'])
    
    pass_info = datFram.copy()
    geometry_temp = [Point(xy) for xy in zip(datFram['LON'], datFram['LAT'])]
    crs = {'init': 'epsg:4326'}
    gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)
    gdf_buff = gdf_buff.to_crs(epsg=32610)
    gdf_buff['geometry'] = gdf_buff.geometry.buffer(30) 
    gdf_tog = pd.merge(gdf_buff,pass_info,on = ['PEAK_NUM', 'EPOCHSTART', 'EPOCH', 'DATETIME', 'CH4', 'LON', 'LAT',
       'CH4_BASELINE', 'CH4_THRESHOLD', 'PEAK_DIST_M', 'PEAK_CH4', 'TCH4',
       'PERIOD5MIN'])
    gdf_bind_pks = gdf_tog.dissolve(by = 'PEAK_NUM',as_index=False)[['PEAK_NUM','geometry']]
    
    if gdf_bind_pks.shape[0] > 1:
        data_overlap = gpd.GeoDataFrame(crs=gdf_bind_pks.crs)
        data_temp = gdf_bind_pks.copy()
        for index, row in data_temp.iterrows():
            data_temp1=data_temp.loc[data_temp.PEAK_NUM!=row.PEAK_NUM,]
            # check if intersection occured
            overlaps=data_temp1[data_temp1.geometry.overlaps(row.geometry)]['PEAK_NUM'].tolist()
            #print(len(overlaps))
            if len(overlaps)>0:
                #print(len(overlaps))
                temp_list=[]
                #print('checking')
                # compare the area with threshold 
                for y in overlaps:
                    temp_area=gpd.overlay(data_temp.loc[data_temp.PEAK_NUM==y,],data_temp.loc[data_temp.PEAK_NUM==row.PEAK_NUM,],how='intersection')
                    temp_area=temp_area.loc[temp_area.geometry.area>=0.001]
                    if temp_area.shape[0]>0:
                        temp_union = gpd.overlay(data_temp.loc[data_temp.PEAK_NUM==y,],data_temp.loc[data_temp.PEAK_NUM==row.PEAK_NUM,],how='union')
                        data_overlap=gpd.GeoDataFrame(pd.concat([temp_union,data_overlap],ignore_index=True),crs=data_temp.crs)
      
        if data_overlap.size > 0: 
                firstnull = data_overlap[data_overlap.PEAK_NUM_1.isnull()]
                firstnull['PEAK_NUM_1'] = firstnull['PEAK_NUM_2']
                
                secnull = data_overlap[data_overlap.PEAK_NUM_2.isnull()]
                secnull['PEAK_NUM_2'] = secnull['PEAK_NUM_1']
                
                withoutNA = data_overlap.dropna()
                allTog = pd.concat([firstnull,secnull,withoutNA]).reset_index()
                
                over = allTog.copy()
                over['sorted']=over.apply(lambda y: sorted([y['PEAK_NUM_1'],y['PEAK_NUM_2']]),axis=1)
                over['sorted']=over.sorted.apply(lambda y: ''.join(y))
                over = over.drop_duplicates('sorted')
                over['combined']= [list(x) for x in list(over.loc[:,['PEAK_NUM_1','PEAK_NUM_2']].to_numpy())]
                over['VER_NUM'] = over.apply(lambda y: y.combined,axis=1)
                over['min_val']=over.apply(lambda y: min(y.combined),axis=1)
                over=over.reset_index()[['PEAK_NUM_1','PEAK_NUM_2','geometry','combined','min_val']]
                
                
                #over = data_overlap.copy()
                #over['sorted']=over.apply(lambda y: sorted([y['PEAK_NUM_1'],y['PEAK_NUM_2']]),axis=1)
                #over['sorted']=over.sorted.apply(lambda y: ''.join(y))
                #over = over.drop_duplicates('sorted')
                #over['combined']= [list(x) for x in list(over.loc[:,['PEAK_NUM_1','PEAK_NUM_2']].to_numpy())]
                #over['VER_NUM'] = over.apply(lambda y: y.combined,axis=1)
                #over['min_val']=over.apply(lambda y: min(y.combined),axis=1)
                #over=over.reset_index()[['PEAK_NUM_1','PEAK_NUM_2','geometry','combined','min_val']]
                        
                overcop = over.copy()
                overcop['recombine'] = overcop.combined
                for index, row in overcop.iterrows():
                    united = row.recombine
                    #print(united)
                    for index2, row2 in overcop.iterrows():
                        united_temp = unIfInt(united,row2.recombine)
                        if united_temp != None:
                            united = united_temp
                            #print(united)
                    overcop.recombine[index] = united
                    #print(united)
                    del(united)   
                overcop['recombine']= overcop.apply(lambda y: sorted(y.recombine),axis=1)
                overcop['min_read'] = overcop.apply(lambda y: min(y.recombine),axis=1)
                newOverlap = overcop.dissolve(by='min_read',as_index=False).loc[:,['min_read','geometry','recombine']]

  
                combined = gdf_bind_pks.copy()
                combined['recombine'] = [list(x) for x in list(combined.loc[:,['PEAK_NUM']].to_numpy())]
                combined['numtimes'] = 1
                combined['newgeo'] = combined.geometry
                combined['min_read'] = combined.PEAK_NUM
                for index,row in combined.iterrows():
                    for index2,row2 in newOverlap.iterrows():
                        if row.PEAK_NUM in row2.recombine:
                            combined.recombine[index] = row2.recombine
                            combined.newgeo[index] = row2.geometry
                            combined.min_read[index] = row2.min_read
                combined['numtimes'] = combined.apply(lambda y: len(y.recombine),axis = 1)
                combined_reduced = combined[['PEAK_NUM','newgeo','recombine','numtimes','min_read']]
                gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['PEAK_NUM'])
                gdf_pass_pks['verified'] = gdf_pass_pks.apply(lambda y: (True if y.numtimes > 1 else False),axis=1 )
        if data_overlap.size == 0:
           gdf_pass_pks = gdf_bind_pks.copy()
           gdf_pass_pks['min_read']= gdf_pass_pks['PEAK_NUM']
           gdf_pass_pks['numtimes'] = 1
           gdf_pass_pks['newgeo'] = gdf_pass_pks.geometry
           gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['PEAK_NUM']].to_numpy())]
           gdf_pass_pks['verified'] = False
           gdf_pass_pks['oldgeo'] = gdf_pass_pks.geometry.copy()
           gdf_pass_pks['geometry'] = gdf_pass_pks.newgeo.copy()
    if gdf_bind_pks.shape[0] == 1:
        gdf_pass_pks = gdf_bind_pks.copy()
        gdf_pass_pks['min_read']= gdf_pass_pks['PEAK_NUM']
        gdf_pass_pks['numtimes'] = 1
        gdf_pass_pks['newgeo'] = gdf_pass_pks.geometry
        gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['PEAK_NUM']].to_numpy())]
        gdf_pass_pks['verified'] = False
        epdat = pass_info[['PEAK_NUM','EPOCHSTART']]
        gdf_pass_pks = pd.merge(gdf_pass_pks,epdat,on = ['PEAK_NUM']) 
    gdf_pass_pks['oldgeo'] = gdf_pass_pks.geometry.copy()
    gdf_pass_pks['geometry'] = gdf_pass_pks.newgeo.copy()
    del(gdf_pass_pks['newgeo'])
    gdf_pass_pks['pass'] = whichpass
    gdf_tot = pd.merge(gdf_pass_pks,datFram_wtLocMax,on = ['PEAK_NUM'])
    #gdf_pass_pks.to_csv(new_loc, index = False)
    #gdf_pass_pks.to_file(str(file[0:20]) + '.shp')
    #return(gdf_pass_pks)
        ## condense by peak_num
    gdfcop = gdf_tot.copy()[['PEAK_NUM','geometry','min_read','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    gdfcop = gdfcop.to_crs(epsg=32610)
    gdfcop.to_file(new_loc_json, driver="GeoJSON")


    gdf_tot.to_csv(new_loc, index = False)
    return(gdf_tot)
        
#def filterPeak(xCar,xDate,xDir,xFilename, outFolder,whichpass = 0):
#    import pandas as pd #
#    import contextily as ctx
#    import os, sys, datetime, time, math, csv, numpy,gzip
#    from math import radians, sin, cos, sqrt, asin
#    import geopandas as gpd
#    from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
#    import matplotlib.pyplot as plt
#
#    file_loc = xDir + xFilename
#    new_loc = outFolder + "Filtered" + xFilename
#    
#    datFram = pd.read_csv(file_loc)
#    
#    datFram_cent =  datFram.copy()
#    datFram_cent['CH4_AB'] = datFram_cent['CH4'].sub(datFram_cent['CH4_BASELINE'], axis = 0) 
#    maxch4 = datFram_cent.groupby('PEAK_NUM',as_index = False).CH4_AB.max().rename(columns = {'CH4_AB':'pk_maxCH4_AB'})
#    datFram_cent_wmax = pd.merge(datFram_cent,maxch4,on = ['PEAK_NUM'])
#   
#    datFram_wtLoc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB').rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'})
#    datFram_wtLocMax = pd.merge(datFram_wtLoc,maxch4,on = ['PEAK_NUM'])
#    
#    pass_info = datFram.copy()
#    geometry_temp = [Point(xy) for xy in zip(datFram['LON'], datFram['LAT'])]
#    crs = {'init': 'epsg:4326'}
#    gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)
#    gdf_buff = gdf_buff.to_crs(epsg=32610)
#    gdf_buff['geometry'] = gdf_buff.geometry.buffer(30) 
#    gdf_tog = pd.merge(gdf_buff,pass_info,on = ['PEAK_NUM', 'EPOCHSTART', 'EPOCH', 'DATETIME', 'CH4', 'LON', 'LAT',
#       'CH4_BASELINE', 'CH4_THRESHOLD', 'PEAK_DIST_M', 'PEAK_CH4', 'TCH4',
#       'PERIOD5MIN'])
#    gdf_bind_pks = gdf_tog.dissolve(by = 'PEAK_NUM',as_index=False)[['PEAK_NUM','geometry']]
#    
#    if gdf_bind_pks.shape[0] > 1:
#        data_overlap = gpd.GeoDataFrame(crs=gdf_bind_pks.crs)
#        data_temp = gdf_bind_pks.copy()
#        for index, row in data_temp.iterrows():
#            data_temp1=data_temp.loc[data_temp.PEAK_NUM!=row.PEAK_NUM,]
#            # check if intersection occured
#            overlaps=data_temp1[data_temp1.geometry.overlaps(row.geometry)]['PEAK_NUM'].tolist()
#            #print(len(overlaps))
#            if len(overlaps)>0:
#                #print(len(overlaps))
#                temp_list=[]
#                #print('checking')
#                # compare the area with threshold 
#                for y in overlaps:
#                    temp_area=gpd.overlay(data_temp.loc[data_temp.PEAK_NUM==y,],data_temp.loc[data_temp.PEAK_NUM==row.PEAK_NUM,],how='intersection')
#                    temp_area=temp_area.loc[temp_area.geometry.area>=0.001]
#                    if temp_area.shape[0]>0:
#                        temp_union = gpd.overlay(data_temp.loc[data_temp.PEAK_NUM==y,],data_temp.loc[data_temp.PEAK_NUM==row.PEAK_NUM,],how='union')
#                        data_overlap=gpd.GeoDataFrame(pd.concat([temp_union,data_overlap],ignore_index=True),crs=data_temp.crs)
#      
#        if data_overlap.size > 0: 
#                firstnull = data_overlap[data_overlap.PEAK_NUM_1.isnull()]
#                firstnull['PEAK_NUM_1'] = firstnull['PEAK_NUM_2']
#                
#                secnull = data_overlap[data_overlap.PEAK_NUM_2.isnull()]
#                secnull['PEAK_NUM_2'] = secnull['PEAK_NUM_1']
#                
#                withoutNA = data_overlap.dropna()
#                allTog = pd.concat([firstnull,secnull,withoutNA]).reset_index()
#                
#                over = allTog.copy()
#                over['sorted']=over.apply(lambda y: sorted([y['PEAK_NUM_1'],y['PEAK_NUM_2']]),axis=1)
#                over['sorted']=over.sorted.apply(lambda y: ''.join(y))
#                over = over.drop_duplicates('sorted')
#                over['combined']= [list(x) for x in list(over.loc[:,['PEAK_NUM_1','PEAK_NUM_2']].to_numpy())]
#                over['VER_NUM'] = over.apply(lambda y: y.combined,axis=1)
#                over['min_val']=over.apply(lambda y: min(y.combined),axis=1)
#                over=over.reset_index()[['PEAK_NUM_1','PEAK_NUM_2','geometry','combined','min_val']]
#                
#                
#                #over = data_overlap.copy()
#                #over['sorted']=over.apply(lambda y: sorted([y['PEAK_NUM_1'],y['PEAK_NUM_2']]),axis=1)
#                #over['sorted']=over.sorted.apply(lambda y: ''.join(y))
#                #over = over.drop_duplicates('sorted')
#                #over['combined']= [list(x) for x in list(over.loc[:,['PEAK_NUM_1','PEAK_NUM_2']].to_numpy())]
#                #over['VER_NUM'] = over.apply(lambda y: y.combined,axis=1)
#                #over['min_val']=over.apply(lambda y: min(y.combined),axis=1)
#                #over=over.reset_index()[['PEAK_NUM_1','PEAK_NUM_2','geometry','combined','min_val']]
#                        
#                overcop = over.copy()
#                overcop['recombine'] = overcop.combined
#                for index, row in overcop.iterrows():
#                    united = row.recombine
#                    #print(united)
#                    for index2, row2 in overcop.iterrows():
#                        united_temp = unIfInt(united,row2.recombine)
#                        if united_temp != None:
#                            united = united_temp
#                            #print(united)
#                    overcop.recombine[index] = united
#                    #print(united)
#                    del(united)   
#                overcop['recombine']= overcop.apply(lambda y: sorted(y.recombine),axis=1)
#                overcop['min_read'] = overcop.apply(lambda y: min(y.recombine),axis=1)
#                newOverlap = overcop.dissolve(by='min_read',as_index=False).loc[:,['min_read','geometry','recombine']]
#
#  
#                combined = gdf_bind_pks.copy()
#                combined['recombine'] = [list(x) for x in list(combined.loc[:,['PEAK_NUM']].to_numpy())]
#                combined['numtimes'] = 1
#                combined['newgeo'] = combined.geometry
#                combined['min_read'] = combined.PEAK_NUM
#                for index,row in combined.iterrows():
#                    for index2,row2 in newOverlap.iterrows():
#                        if row.PEAK_NUM in row2.recombine:
#                            combined.recombine[index] = row2.recombine
#                            combined.newgeo[index] = row2.geometry
#                            combined.min_read[index] = row2.min_read
#                combined['numtimes'] = combined.apply(lambda y: len(y.recombine),axis = 1)
#                combined_reduced = combined[['PEAK_NUM','newgeo','recombine','numtimes','min_read']]
#                gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['PEAK_NUM'])
#                gdf_pass_pks['verified'] = gdf_pass_pks.apply(lambda y: (True if y.numtimes > 1 else False),axis=1 )
#        if data_overlap.size == 0:
#           gdf_pass_pks = gdf_bind_pks.copy()
#           gdf_pass_pks['min_read']= gdf_pass_pks['PEAK_NUM']
#           gdf_pass_pks['numtimes'] = 1
#           gdf_pass_pks['newgeo'] = gdf_pass_pks.geometry
#           gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['PEAK_NUM']].to_numpy())]
#           gdf_pass_pks['verified'] = False
#           gdf_pass_pks['oldgeo'] = gdf_pass_pks.geometry.copy()
#           gdf_pass_pks['geometry'] = gdf_pass_pks.newgeo.copy()
#    if gdf_bind_pks.shape[0] == 1:
#        gdf_pass_pks = gdf_bind_pks.copy()
#        gdf_pass_pks['min_read']= gdf_pass_pks['PEAK_NUM']
#        gdf_pass_pks['numtimes'] = 1
#        gdf_pass_pks['newgeo'] = gdf_pass_pks.geometry
#        gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['PEAK_NUM']].to_numpy())]
#        gdf_pass_pks['verified'] = False
#        epdat = pass_info[['PEAK_NUM','EPOCHSTART']]
#        gdf_pass_pks = pd.merge(gdf_pass_pks,epdat,on = ['PEAK_NUM']) 
#    gdf_pass_pks['oldgeo'] = gdf_pass_pks.geometry.copy()
#    gdf_pass_pks['geometry'] = gdf_pass_pks.newgeo.copy()
#    del(gdf_pass_pks['newgeo'])
#    gdf_pass_pks['pass'] = whichpass
#    gdf_tot = pd.merge(gdf_pass_pks,datFram_wtLocMax,on = ['PEAK_NUM'])
#    #gdf_pass_pks.to_csv(new_loc, index = False)
#    #gdf_pass_pks.to_file(str(file[0:20]) + '.shp')
#    #return(gdf_pass_pks)
#    return(gdf_tot)

    
### COMBINE PASS WILL THEN COMBINE THE DIFFERENT PASSES TOGETHER TO BE IN ONE SHAPE FILE

def unique(my_list): 
   return [x for x in my_list if x not in locals()['_[1]']]

def unIfInt(a,b):
    if len(intersect(a,b)) != 0:
        return( list(set(a).union(b)))

def IsInPK (peakNumRow,listPksRow):
    if peakNumRow.PEAK_NUM in listPksRow.combined:
        return(True)
    else:
        return(False)
        
def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

def passCombine (firstgrp, secondgrp):
    import pandas as pd #
    import contextily as ctx
    import os, sys, datetime, time, math, csv, numpy,gzip
    from math import radians, sin, cos, sqrt, asin
    import geopandas as gpd
    from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
    import matplotlib.pyplot as plt

    first_geo = firstgrp.geometry
    crs = {'init': 'epsg:4326'}
    firstgrp = gpd.GeoDataFrame(firstgrp,crs = crs,geometry = first_geo)
    
    first_pks = firstgrp[['PEAK_NUM','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    sec_pks = secondgrp[['PEAK_NUM','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    tot_pks = pd.concat([first_pks,sec_pks])
    
    first_dis = firstgrp.dissolve(by='min_read',as_index=False)[['min_read','geometry','recombine','verified','pass']]
    sec_dis = secondgrp.dissolve(by='min_read',as_index=False)[['min_read','geometry','recombine','verified','pass']]
    gdf_bind_pks = pd.concat([first_dis,sec_dis])[['min_read','geometry','recombine']]
    #gdf_tog = pd.concat([firstgrp,secondgrp])
    gdf_tog = pd.concat([firstgrp.drop(['pk_LAT', 'pk_LON','pk_maxCH4_AB'], axis=1),secondgrp.drop(['pk_LAT', 'pk_LON','pk_maxCH4_AB'], axis=1)])

    gdf_bind_pks['prev_read'] = gdf_bind_pks['min_read'].copy()
    gdf_tog['prev_read'] = gdf_tog['min_read'].copy()
    if gdf_bind_pks.shape[0] > 1:
        data_overlap = gpd.GeoDataFrame(crs=gdf_bind_pks.crs)
        data_temp = gdf_bind_pks.copy()
        for index, row in data_temp.iterrows():
            data_temp1=data_temp.loc[data_temp.min_read!=row.min_read,]
            # check if intersection occured
            overlaps=data_temp1[data_temp1.geometry.overlaps(row.geometry)]['min_read'].tolist()
            if len(overlaps)>0:
                temp_list=[]
                #print('checking')
                # compare the area with threshold 
                for y in overlaps:
                    temp_area=gpd.overlay(data_temp.loc[data_temp.min_read==y,],data_temp.loc[data_temp.min_read==row.min_read,],how='intersection')
                    temp_area=temp_area.loc[temp_area.geometry.area>=0]
                    #temp_union = gpd.overlay(data_temp.loc[data_temp.PEAK_NUM==y,],data_temp.loc[data_temp.PEAK_NUM==row.PEAK_NUM,],how='union')
                    if temp_area.shape[0]>0:
                        temp_union = gpd.overlay(data_temp.loc[data_temp.min_read==y,],data_temp.loc[data_temp.min_read==row.min_read,],how='union')
                        data_overlap=gpd.GeoDataFrame(pd.concat([temp_union,data_overlap],ignore_index=True),crs=data_temp.crs)    
    
        if data_overlap.size != 0: 
            firstnull = data_overlap[data_overlap.min_read_1.isnull()]
            firstnull['min_read_1'] = firstnull['min_read_2']
                
            secnull = data_overlap[data_overlap.min_read_2.isnull()]
            secnull['min_read_2'] = secnull['min_read_1']
                
            withoutNA = data_overlap.dropna()
            allTog = pd.concat([firstnull,secnull,withoutNA]).reset_index()
            
            over = allTog.copy()
            over['sorted']=over.apply(lambda y: sorted([y['min_read_1'],y['min_read_2']]),axis=1)
            over['sorted']=over.sorted.apply(lambda y: ''.join(y))
            over = over.drop_duplicates('sorted')
            over['combined']= [list(x) for x in list(over.loc[:,['min_read_1','min_read_2']].to_numpy())]
            over['VER_NUM'] = over.apply(lambda y: y.combined,axis=1)
            over['min_val']=over.apply(lambda y: min(y.combined),axis=1)
            over=over.reset_index()[['min_read_1','min_read_2','geometry','combined','min_val']]
                
            overcop = over.copy()
            overcop['recombine'] = overcop.combined
            
            for index,row in overcop.iterrows():
                rowwoo = row
                first_thing = first_dis[first_dis['min_read']== row.min_read_1].loc[:,['recombine']]

                #first_thing = firstpass[firstpass['min_read']== row.min_read_1].loc[:,['recombine']]
                firstcomb = first_thing.recombine.explode()
                first_list = firstcomb.reset_index().recombine
                #first_list = firstpass.loc[:,firstpass['min_read']== row.min_read_1].recombine
                #first_list.explode()
                #first_list = firstpass.loc[:,firstpass['min_read']== row.min_read_1].recombine.explode()[0]
               # second_list = secondpass[secondpass['min_read']== row.min_read_2].recombine
                sec_thing = sec_dis[sec_dis['min_read']== row.min_read_1].loc[:,['recombine']]

                #sec_thing = secondpass[secondpass['min_read']== row.min_read_2].loc[:,['recombine']]
                seccomb = sec_thing.recombine.explode()
                sec_list = seccomb.reset_index().recombine
               
                firstdf = pd.DataFrame(first_list)
                secdf = pd.DataFrame(sec_list)
    
                
                tot_df = pd.concat([firstdf,secdf])
                tot_list = tot_df.recombine.tolist()
                
                overcop.recombine[index] = tot_list
    ## this recombines the lists together to have the combined entries together?
       
            overcop['recombine']= overcop.apply(lambda y: sorted(y.recombine),axis=1)
            overcop['min_read'] = overcop.apply(lambda y: min(y.recombine),axis=1)
            newOverlap = overcop.dissolve(by='min_read',as_index=False).loc[:,['min_read','geometry','recombine']]
        
      
            combined = gdf_bind_pks.copy()
            #combined['recombine'] = [list(x) for x in list(combined.loc[:,['min_read']].to_numpy())]
            combined['numtimes'] = 1
            combined['newgeo'] = combined.geometry
            combined['oldgeo'] = combined.geometry.copy()

            #combined['min_read'] = combined.min_rea
            combined = combined.reset_index()
            for index,row in combined.iterrows():
                for index2,row2 in newOverlap.iterrows():
                    if row.min_read in row2.recombine:
                        woo = row2
                        indo = index
                        combined.recombine[index] = row2.recombine
                        #combined.newgeo[index] = row2.geometry
                        combined.newgeo[index] = row2.geometry
                        combined.min_read[index] = row2.min_read
                    
        
            combined['numtimes'] = combined.apply(lambda y: len(y.recombine),axis = 1)
            combined['geometry'] = combined.newgeo.copy()
            del(combined['newgeo'])
            combined_reduced = combined[['min_read','geometry','oldgeo','recombine','numtimes','prev_read']]
            #gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['min_read'])
            gdf_tog = gdf_tog.drop('min_read',axis=1)
            gdf_tog = gdf_tog.drop('numtimes',axis=1)
            gdf_tog = gdf_tog.drop('recombine',axis=1)
           # gdf_tog = gdf_tog.drop('newgeo',axis=1)
            gdf_tog['firstgeo'] = gdf_tog.oldgeo.copy()
            gdf_tog['secondgeo'] = gdf_tog.geometry.copy()
            del(gdf_tog['geometry'])
            del(gdf_tog['oldgeo'])

                                   
            gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['prev_read'])

            gdf_pass_pks['verified'] = gdf_pass_pks.apply(lambda y: (True if y.numtimes > 1 else False),axis=1 )
        if data_overlap.size == 0:
            gdf_pass_pks = gdf_tog
    ## didnt adress if the bind shape was only size only 1
        gdf_tot_pks = pd.merge(gdf_pass_pks,tot_pks,on = ['PEAK_NUM'])
    #return(gdf_pass_pks)
    return(gdf_tot_pks)

    
    
### WEIGHTED AVERAGE OF THE LOCATION

#def weightedLoc(df,lat,lon,by,val2avg):
#    #df_use = df[[str(lat),str(lon),str(by)]]
#    df_use = df[[(lat),(lon),(by),val2avg]]
#    by_val = ""
#    wts = []
#    x_num = []
#    y_num = []
#    results = pd.DataFrame(data= None,columns = [(by),(lat),(lon)])
#
#    for index,x in df_use.iterrows():
#        if by_val and (by_val != df_use.iloc[index][by]):
#            long = sum(x_num)/sum(wts)
#            latitude = sum(y_num)/sum(wts)
#            del x_num[:]
#            del y_num[:]
#            del wts[:]
#            results = results.append({(by):by_val,lon:long,lat:latitude},ignore_index = True)
#            
#        by_val = df_use.iloc[index][(by)]
#        wts.append(float(df_use.iloc[index][val2avg]))
#        x_num.append(float(df_use.iloc[index][(lon)])*float(df_use.iloc[index][val2avg]))
#        y_num.append(float(df_use.iloc[index][(lat)])*float(df_use.iloc[index][val2avg]))
#    if wts:
#        long = sum(x_num)/sum(wts)
#        latitude = sum(y_num)/sum(wts)
#        del x_num[:]
#        del y_num[:]
#        del wts[:]
#        results = results.append({by:by_val,lon:long,lat:latitude},ignore_index = True)
#    return(results)
    
    

def wt_time_Locs(wt,loc):
    return(wt*loc)
    
def sumthing(thing):
    return(sum(thing))


def weightedLoc(df,lat,lon,by,val2avg):
    df_use = df[[(lat),(lon),(by),val2avg]]
    by_val = ""
    wts = []
    x_num = []
    y_num = []
    x_num_wt = []
    y_num_wt = []
    results = pd.DataFrame(data= None,columns = [(by),(lat),(lon)])

    df_use['lat_wt'] = df_use.apply(lambda y: y[lat] * y[val2avg],axis = 1)
    df_use['lon_wt'] = df_use.apply(lambda y: y[lon] * y[val2avg],axis = 1)


    #sumwts =pd.DataFrame( df_use.groupby('min_read').apply(lambda y: sumthing(y['pk_maxCH4_AB'])),columns = {'totwts'})
    sumwts = pd.DataFrame( df_use.groupby(str(by)).apply(lambda y: sumthing(y[str(val2avg)])),columns = {'totwts'})
    sumwts['min_reads'] = sumwts.index.copy()
    sumwts = sumwts.reset_index(drop=True)
    
    sumwts = sumwts.rename(columns={"min_reads": str(by)})

    #totlats = pd.DataFrame(df_use.groupby('min_read').apply(lambda y: sumthing(y['lat_wt'])),columns = ['totlats'])
    totlats = pd.DataFrame(df_use.groupby(str(by)).apply(lambda y: sumthing(y['lat_wt'])),columns = ['totlats'])

    totlats['min_reads'] = totlats.index.copy()
    totlats = totlats.reset_index(drop=True)
    totlats = totlats.rename(columns={"min_reads": str(by)})

    #totlons = pd.DataFrame(df_use.groupby('min_read').apply(lambda y: sumthing(y['lon_wt'])),columns = ['totlons'])
    totlons = pd.DataFrame(df_use.groupby(str(by)).apply(lambda y: sumthing(y['lon_wt'])),columns = ['totlons'])

    totlons['min_reads'] = totlons.index.copy()
    totlons = totlons.reset_index(drop=True)
    totlons = totlons.rename(columns={"min_reads": str(by)})

    #df_use = pd.merge(totlats,df_use,on = 'min_read')
    #df_use = pd.merge(totlons,df_use,on = 'min_read')
    #df_use = pd.merge(sumwts,df_use,on = 'min_read')
    
    df_use = pd.merge(totlats,df_use,on = str(by))
    df_use = pd.merge(totlons,df_use,on = str(by))
    
    df_use = pd.merge(sumwts,df_use,on = str(by))
    

    df_use['overall_LON'] = df_use.apply(lambda y: y['totlons']/y['totwts'],axis = 1)
    df_use['overall_LAT'] = df_use.apply(lambda y: y['totlats']/y['totwts'],axis = 1)

    toreturn = df_use[[str(by),'overall_LON','overall_LAT']].drop_duplicates()
    toreturn = toreturn.rename(columns = {'overall_LON':str(lon),'overall_LAT':str(lat)})

    return(toreturn)   


 
#def weightedLoc(df,lat,lon,by,val2avg):
#    #df_use = df[[str(lat),str(lon),str(by)]]
#    df_use = df[[(lat),(lon),(by),val2avg]]
#    by_val = ""
#    wts = []
#    x_num = []
#    y_num = []
#    x_num_wt = []
#    y_num_wt = []
#    results = pd.DataFrame(data= None,columns = [(by),(lat),(lon)])
#    add1 = 0
#    if  df_use.shape[0] == 1:
#        long = float(df_use.iloc[0][(lon)])
#        latitude =  float(df_use.iloc[0][(lat)])
#        by_val = df_use.iloc[0][(by)]
#        results = results.append({(by):by_val,lon:long,lat:latitude},ignore_index = True)
#    if  df_use.shape[0] > 1:
#        for index,x in df_use.iterrows():
#            #print(index)
#            if df_use.iloc[index][by] == 'CSULi_1568812598.4':
#                add1 += 1
#            #if index == 12:
#             #   print(df_use.iloc[index][by])
#            if by_val and (by_val != df_use.iloc[index][by]) and index != 1 :
#                #print(index)
#                if sum(wts)!= 0:
#                    long = sum(x_num)/sum(wts)
#                    latitude = sum(y_num)/sum(wts)
#                
#                
#                del x_num[:]
#                del y_num[:]
#                del x_num_wt[:]
#                del y_num_wt[:]
#                del wts[:]
#                results = results.append({(by):by_val,lon:long,lat:latitude},ignore_index = True)
#                
#            by_val = df_use.iloc[index][(by)]
#            wts.append(float(df_use.iloc[index][val2avg]))
#            x_num.append(float(df_use.iloc[index][(lon)])*float(df_use.iloc[index][val2avg]))
#            y_num.append(float(df_use.iloc[index][(lat)])*float(df_use.iloc[index][val2avg]))
#            
#            if by_val and (by_val == df_use.iloc[index][by]) and index == len(df_use.index)-1 :
#                #print(index)
#                long = sum(x_num)/sum(wts)
#                latitude = sum(y_num)/sum(wts)
#                del x_num[:]
#                del y_num[:]
#                del x_num_wt[:]
#                del y_num_wt[:]
#                del wts[:]
#                results = results.append({(by):by_val,lon:long,lat:latitude},ignore_index = True)
#            if by_val and len(df_use[df_use[by]==by_val][by])==1:
#                #print(df_use[df_use[by]==by_val][by])
#                #long = sum(x_num)/sum(wts)
#                #latitude = sum(y_num)/sum(wts)
#                #print("issue")
#                #print(index)
#                long = float(df_use.iloc[index][(lon)])
#                #print(long)
#                latitude = float(df_use.iloc[index][(lat)])
#                #print(latitude)
#                #print(long)
#                #print(latitude)
#                
#                del x_num[:]
#                del y_num[:]
#                del x_num_wt[:]
#                del y_num_wt[:]
#                del wts[:]
#                results = results.append({(by):by_val,lon:long,lat:latitude},ignore_index = True)
#                 
#    
#    #if wts:
#     #   long = sum(x_num)/sum(wts)
#      #  latitude = sum(y_num)/sum(wts)
#      #  del x_num[:]
#      #  del y_num[:]
#      #  del wts[:]
#       # results = results.append({by:by_val,lon:long,lat:latitude},ignore_index = True)
#    results = results[results[lat].isna()== False].reset_index()[[(by),(lat),(lon)]]
#    results = results.drop_duplicates()
#    return(results)
#    
    
#def weightedLoc(df,lat,lon,by,val2avg):
#    import pandas as pd
#    #df_use = df[[str(lat),str(lon),str(by)]]
#    df_use = df[[(lat),(lon),(by),val2avg]]
#    by_val = ""
#    wts = []
#    x_num = []
#    y_num = []
#    x_num_wt = []
#    y_num_wt = []
#    results = pd.DataFrame(data= None,columns = [(by),(lat),(lon)])
#    
#    if  df_use.shape[0] == 1:
#        long = float(df_use.iloc[0][(lon)])
#        latitude =  float(df_use.iloc[0][(lat)])
#        by_val = df_use.iloc[0][(by)]
#        results = results.append({(by):by_val,lon:long,lat:latitude},ignore_index = True)
#    if  df_use.shape[0] > 1:
#        for index,x in df_use.iterrows():
#            if by_val and (by_val != df_use.iloc[index][by]) :
#                long = sum(x_num)/sum(wts)
#                latitude = sum(y_num)/sum(wts)
#                del x_num[:]
#                del y_num[:]
#                del x_num_wt[:]
#                del y_num_wt[:]
#                del wts[:]
#                results = results.append({(by):by_val,lon:long,lat:latitude},ignore_index = True)
#                
#            by_val = df_use.iloc[index][(by)]
#            wts.append(float(df_use.iloc[index][val2avg]))
#            x_num.append(float(df_use.iloc[index][(lon)])*float(df_use.iloc[index][val2avg]))
#            y_num.append(float(df_use.iloc[index][(lat)])*float(df_use.iloc[index][val2avg]))
#            
#            if by_val and (by_val == df_use.iloc[index][by]) and index == len(df_use.index)-1 :
#                print(index)
#                long = sum(x_num)/sum(wts)
#                latitude = sum(y_num)/sum(wts)
#                del x_num[:]
#                del y_num[:]
#                del x_num_wt[:]
#                del y_num_wt[:]
#                del wts[:]
#                results = results.append({(by):by_val,lon:long,lat:latitude},ignore_index = True)
#            if by_val and len(df_use[df_use[by]==by_val][by])==1:
#                print(by)
#                long = sum(x_num)/sum(wts)
#                latitude = sum(y_num)/sum(wts)
#                del x_num[:]
#                del y_num[:]
#                del x_num_wt[:]
#                del y_num_wt[:]
#                del wts[:]
#                results = results.append({(by):by_val,lon:long,lat:latitude},ignore_index = True)
#                 
#
#    #if wts:
#     #   long = sum(x_num)/sum(wts)
#      #  latitude = sum(y_num)/sum(wts)
#      #  del x_num[:]
#      #  del y_num[:]
#      #  del wts[:]
#       # results = results.append({by:by_val,lon:long,lat:latitude},ignore_index = True)
#    return(results)


####   
    
def verPk(totalData):
    import pandas as pd #
    import contextily as ctx
    import os, sys, datetime, time, math, csv, numpy,gzip
    from math import radians, sin, cos, sqrt, asin
    import geopandas as gpd
    from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
    import matplotlib.pyplot as plt

    totalData = totalData[totalData.numtimes != 1]
    pkRed = totalData[['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']].drop_duplicates().reset_index()
    verLoc = weightedLoc(pkRed,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB')
    pkRed['logCH4']=pkRed.apply(lambda y: log(y.pk_maxCH4_AB),axis=1)
    mnVals = pkRed.groupby('min_read',as_index = False).logCH4.mean()
    together = pd.merge(verLoc,mnVals,on = ['min_read'])
    geometry_temp = [Point(xy) for xy in zip(together['pk_LON'], together['pk_LAT'])]
    crs = {'init': 'epsg:32610'}
    tog_dat = gpd.GeoDataFrame(together,crs = crs,geometry=geometry_temp)
    tog_dat = tog_dat.to_crs(epsg = 3857)


def estEmissions(excessCH4):
    import math
    a = 0.4630664
    b = 0.7443749
    a1 = 1.2889
    b1 = 0.35232
    a2 = 1.755891
    b2 = 0.4438203
    
    m = math.exp((excessCH4 - a)/b)
   # if m < math.exp(3.157):
    #    if m < math.exp(2):
     #       m = math.exp((np.log(m) - a1)/b1)
      #  if m > math.exp(2):
       #     m = math.exp((np.log(m) - a2)/b2)
    return(m)
        
    
        