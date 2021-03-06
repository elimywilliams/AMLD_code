#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
## all functions colin
Created on Mon Apr  6 14:51:23 2020

@author: emilywilliams
"""
import pandas as pd

#mport contextily as ctx
#import pandas as pd
import datetime
from numpy import log
import geopandas as gpd
from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
from datetime import datetime


def haversine(lat1, lon1, lat2, lon2, radius=6371): # 6372.8 = earth radius in kilometers
    from math import radians, sin, cos, sqrt, asin

    dLat = radians(lat2 - lat1)
    dLon = radians(lon2 - lon1)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    c = 2*asin(sqrt(sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2))

    return radius*c*1000 #return in meters


def getQuad(x,y):
    try:
        x = int(x)
        y = int(y)
    except ValueError:
        return(0)
    
    if y>=0 and x>0:
        return(1)
    elif y>=0 and x<0:
        return(2)
    elif y<0 and x<0:
        return(3)
    else:
        return(4)
        


def calcBearing(lat1,lat2,long1,long2,radians):
    from math import atan2
    from numpy import pi
    from math import radians, sin, cos, sqrt, asin

    lat1r = lat1* (pi/180)
    lat2r = lat2*(pi/180)
    long1r = long1*(pi/180)
    long2r = long2*(pi/180)
    X =cos(lat2r)*sin(long2r-long1r)
    Y = cos(lat1r)*sin(lat2r) - (sin(lat1r)*cos(lat2r)*cos(long2r-long1r))
  
    theta = atan2(X,Y)
    theta = theta % (2*pi)
    
    if not radians:
        return (theta * 180/pi)
    elif radians:
        return(theta)
        
def calcTheta(U,V,quad,h_length,radians):
    theta = np.arcsin(U/h_length)
    import numpy as np
    if quad == 1:
        theta = theta
    elif quad == 2:
        theta = -theta + np.pi/2
    elif quad -- 3:
        theta = np.pi/2 + theta + np.pi
    elif quad == 4:
        theta = 3*np.pi/2
    theta = 2*np.pi - theta
    if not radians:
        theta = theta * 180/np.pi
        
    return(theta)



def ProcessRawDataEng( xCar, xDate, xDir, xFilename, bFirst, gZIP, xOut,initialTimeBack,shift):
    from datetime import datetime
    import os
    #import os
    import gzip
    import csv
    from math import floor
    try:
        #shift = -4
        #shift = 0

        #xOutDir = xDir 
        
        xMaxCarSpeed = 45.0 / 2.23694     # assumes 45 mph as max, convert to meters per second 
        xMinCarSpeed = 2.0 / 2.23694      # minimum 2 miles per hour
        xMinRSSI = 50  #if RSSI is below this we don't like it
        

        # reading in the data with specific headers
        #          0     1    2    3       4           5    6       7        8        9          10                 11              12           13            14      15      16      17        18         19         20         21         22         23        24   25  26       27           28       29           30       31       32       33  34        35   36   37  38   39       40       41   42       43   44   45   46   47   48   49   50   51     52     53     54
        #sHeader = "Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude"
        #sHeader = "Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude"
        #sHeader = "Time Stamp,Inlet Number,P (mbars),T0 (degC),T5 (degC), Laser PID Readout,Det PID Readout,win0Fit0,win0Fit1,win0Fit3,win1Fit4,win0Fit5,win0Fit6,win0Fit7,win0Fit8,win0Fit9,win1Fit0,win1Fit1,win1Fit2,win1Fit3,win1Fit4,win1Fit5,win1Fit6,Det Bkgd,Ramp Ampl,CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Battery T (degC),FET T (degC),GPS Time,Latitude,Longitude"
        sHeader = "Time Stamp,Inlet Number,P (mbars),T0 (degC),T5 (degC),Laser PID Readout,Det PID Readout,win0Fit0,win0Fit1,win0Fit2,win0Fit3,win0Fit4,win0Fit5,win0Fit6,win0Fit7,win0Fit8,win0Fit9,win1Fit0,win1Fit1,win1Fit2,win1Fit3,win1Fit4,win1Fit5,win1Fit6,Det Bkgd,Ramp Ampl,CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Battery T (degC),FET T (degC),GPS Time,Latitude,Longitude"
        sOutHeader = "DATE,TIME,SECONDS,NANOSECONDS,VELOCITY,U,V,W,BCH4,BRSSI,TCH4,TRSSI,PRESS_MBAR,INLET,TEMPC,CH4,H20,C2H6,R,C2C1,BATTV,POWMV,CURRMA,SOCPER,LAT,LONG\n"
        
        headerNames = sHeader.split(',')
        GPS_loc = 37
        
        
        infoHeader = "FILENAME\n"
        # somehow gZIP is indicating if  it is the first file name (I think if it is 0 then it is the first file)
        if gZIP == 0:
            f = gzip.open(xDir + "/" + xFilename, 'r') #if in python 3, change this to "r" or just "b" can't remember but something about a bit not a string
        else:
            f = open(xDir + "/" + xFilename, 'r')
            f = open(xDir  + xFilename, 'r')

        
        # process    
        #if first time on this car/date, then write header out
        
        xdat = str('20') + xFilename[11:17]
        
        #fnOut = xOutDir + xCar + "_" + xDate.replace("-", "") + "_dat.csv"       #set CSV output for raw data
        #fnLog = xOutDir + xCar + "_" + xDate.replace("-", "") + "_log.csv"       #output for logfile
        
        fnOut = xOut  + xCar + "_" + xdat + "_dat.csv"       #set CSV output for raw data
        fnLog =  xOut  + xCar + "_" + xdat + "_log.csv"       #output for logfile
        infOut = xOut + xCar + "_" + xdat + "_info.csv"
        
        firsttime = int(float(open(xDir + xFilename).readlines().pop(1).split(',')[37][:-4]))
        fnOutTemp = xOut  + xCar + "_" + xdat + "temp_dat.csv"       #
        
        if bFirst:
            #3fOut = open(fnOutTemp, 'w')
            #fOut.write(sOutHeader)
            fLog = open(fnLog, 'w')
            infOut = open(infOut,'w')
            infOut.write(infoHeader)
            print ("fnLog: "+fnOut) 
        if not bFirst:
            fOut = open(fnOut, 'a')
            fLog = open(fnLog, 'a')
            infOut = open(infOut,'a')
        
        fOut = open(fnOutTemp, 'w')
        fOut.write(sOutHeader)

        
        
        #read all lines
        xCntObs = -1
        xCntGoodValues = 0
        rownum = 0 
        for row in f:
            woo = row
            #print(row)
            bGood = True
            if xCntObs < 0:
                bGood = False
                xCntObs += 1
                #lstSfirst = woo.split(',')
                #gpstimefirst = lstSfirst[GPS_loc]
                #firstseconds = floor(float(gpstimefirst))

            if bGood:
                #lstS = row.split(",")
                lstS = row.split(',')
                gpstime = lstS[GPS_loc]
                dtime = lstS[0]
                dt = lstS[1]
                time_dt = lstS[2]
                epoch = lstS[3]
                #nano = lstS[4]
                
                gps_time = lstS[37]
                dateob = datetime.fromtimestamp(int(gps_time[:-4]))
                nano= gps_time[-4:]
                
                #dateob = datetime(int(dt[0:4]),int(dt[5:7]),int(dt[8:10]),int(time_dt[0:2]),int(time_dt[3:5]),int(time_dt[6:8]),int(float(nano)*1e-9))
                
                dtime  = int(dateob.strftime('%Y%m%d%H%M%S'))
                #Date = dateob.strftime('%m%/%d/%Y')
                Date = dateob.strftime('%Y-%m-%d')

                GPS_Time = dateob.strftime('%H%:%M:%S')
                seconds = floor(float(gpstime))
                nano = dateob.strftime('%f')
                
                #dateob = datetime(int(dtime[6:10]),int(dtime[0:2]),int(dtime[3:5]),int(dtime[11:13]),int(dtime[14:16]),int(dtime[17:19]),int(float(dtime[19:23])*1000000))
                #epoch = dateob.strftime('%s.%f')
  
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
                import sys
                if sys.platform.startswith('win'):
                    ## DATE, TIME, SECONDS,NANOSECONDS
                    csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(dateob.strftime('%H:%M:%S'))  + ',' + str(float(pd.to_numeric(dateob.strftime('%S.%f')))) + ',' + str(pd.to_numeric(dateob.strftime('%f')) *1000) + str(',')
                    ## VELOCITY, U,V,W,BCH4,BRSSI,TCH4
                    csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(lstS[26]) + ',' + str('0') + ','+  str(lstS[26]) + ','
                    ## TRSSI, PRESS_MBAR, INLET, TEMPC, CH4, H20,C2H6
                    csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(lstS[26]) + ',' + str(lstS[27]) +',' +  str(lstS[28]) + ','
                    # R, C2C1, BATTV, POWMV,CURRMA, SOCPER,LAT,LONG
                    csvWrite += str(lstS[29]) + ',' + str(lstS[30]) + ',' + str(lstS[31]) + ',' + str(lstS[32]) + ','+ str(lstS[33]) + ',' + str(lstS[34]) + ',' + str(lstS[38]) + str(',') + str(lstS[39])
                
# =============================================================================
#                 if not sys.platform.startswith('win'):
#                     ## DATE, TIME, SECONDS,NANOSECONDS
#                     csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(dateob.strftime('%H:%M:%S'))  + ',' + str((int(floor(pd.to_numeric(dateob.strftime('%s.%f')))))) + ',' + str((pd.to_numeric(dateob.strftime('%f')) *1000)) + str(',')
#                     ## VELOCITY, U,V,W,BCH4,BRSSI,TCH4
#                     csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(lstS[26]) + ',' + str('0') + ','+  str(lstS[26]) + ','
#                     ## TRSSI, PRESS_MBAR, INLET, TEMPC, CH4, H20,C2H6
#                     csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(lstS[26]) + ',' + str(lstS[27]) +',' +  str(lstS[28]) + ','
#                     # R, C2C1, BATTV, POWMV,CURRMA, SOCPER,LAT,LONG
#                     csvWrite += str(lstS[29]) + ',' + str(lstS[30]) + ',' + str(lstS[31]) + ',' + str(lstS[32]) + ','+ str(lstS[33]) + ',' + str(lstS[34]) + ',' + str(lstS[38]) + str(',') + str(lstS[39][:-1]) + str('\n')                   
#                 #fOut.write('\n')
#                 fOut.write(csvWrite)
#                 #fOut.write('\n')
#                 
# =============================================================================
                if not sys.platform.startswith('win'):
                    ## DATE, TIME, SECONDS,NANOSECONDS
                    csvWrite = str(Date) + ',' + str(GPS_Time)  + ',' + str(seconds) + ',' + str(nano) + str(',')
                    ## VELOCITY, U,V,W,BCH4,BRSSI,TCH4
                    csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(lstS[26]) + ',' + str('0') + ','+  str(lstS[26]) + ','
                    ## TRSSI, PRESS_MBAR, INLET, TEMPC, CH4, H20,C2H6
                    csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(lstS[26]) + ',' + str(lstS[27]) +',' +  str(lstS[28]) + ','
                    # R, C2C1, BATTV, POWMV,CURRMA, SOCPER,LAT,LONG
                    csvWrite += str(lstS[29]) + ',' + str(lstS[30]) + ',' + str(lstS[31]) + ',' + str(lstS[32]) + ','+ str(lstS[33]) + ',' + str(lstS[34]) + ',' + str(lstS[38]) + str(',') + str(lstS[39])    
                #fOut.write('\n')
                
                #print( seconds > firstseconds + (60*5))
                if seconds >= (firsttime + (60*float(initialTimeBack))):
                    fOut.write(csvWrite)
                
                
                del(csvWrite)
#                xCntGoodValues += 1
                

            xCntObs += 1

        #sOut = str(gZIP) + "," + str(f) + "," + str(xCntObs) + "," + str(xCntGoodValues) + "\n"
        #fLog.write(sOut)

        infOut.write(str(xFilename)+'\n')

        fOut.close()
        fLog.close()
        infOut.close()
        
        #xDate = dateob.strftime("%Y%m%d")
        
        #newfnOut = xOutDir + xCar + "_" + xDate + "_dat.csv"       #set CSV output for raw data
        #newfnLog = xOutDir + xCar + "_" + xDate + "_log.csv"  
        

        print (xCar + "\t" + xdat + "\t" + fnOut[-22:] + "\t" + str(xCntObs) + "\t" + str(xCntGoodValues) + "\t" + str(gZIP))
        from numpy import pi
        import numpy as np
        def calcVel(timediff,distance):
            if timediff == 0:
                return(0)
            elif timediff != 0:
                return(distance/timediff)
                
        radians = False
        wind_df = pd.read_csv(fnOutTemp)        
        wind_df['QUADRANT'] = wind_df.apply(lambda row: getQuad(row['U'],row['V']),axis=1)
        wind_df['secnan'] = wind_df.apply(lambda row: row['SECONDS'],axis=1) # + row['NANOSECONDS']*1e-9,axis=1)
        wind_df['prev_LAT'] = wind_df.LAT.shift(periods = 1)
        wind_df['next_LAT'] = wind_df.LAT.shift(periods = -1)
        wind_df['prev_LONG'] = wind_df.LONG.shift(periods = 1)
        wind_df['next_LONG'] = wind_df.LONG.shift(periods = -1)
        wind_df['prev_TIME'] = wind_df.secnan.shift(periods = 1)
        wind_df['next_TIME'] = wind_df.secnan.shift(periods = -1)
        wind_df['distance'] = wind_df.apply(lambda row: haversine(row['prev_LAT'],row['prev_LONG'],row['next_LAT'],row['next_LONG']),axis=1)
        wind_df['bearing'] = wind_df.apply(lambda row: calcBearing(row['prev_LAT'],row['next_LAT'],row['prev_LONG'],row['next_LONG'],radians),axis=1)
        wind_df['timediff'] = wind_df.apply(lambda row: row['next_TIME'] - row['prev_TIME'],axis = 1)
        wind_df['VELOCITY'] = wind_df.apply(lambda row:calcVel(row['timediff'],row['distance']),axis=1)
        wind_df['U_cor'] = wind_df.apply(lambda row:row['U'] + row['VELOCITY'],axis = 1)
        wind_df['horz_length'] = wind_df.apply(lambda row: np.sqrt(row['U_cor']**2 + row['V']**2),axis=1)
        wind_df['uncor_theta'] = wind_df.apply(lambda row :calcBearing(row['U_cor'],row['V'],row['QUADRANT'],row['horz_length'],radians),axis = 1)
        wind_df['adj_theta'] = wind_df.apply(lambda row: (row['uncor_theta'] + row['bearing'])%360,axis =1)
        wind_df['totalWind'] = wind_df.apply(lambda row: np.sqrt(row['horz_length']**2 + row['W']**2),axis = 1)
        wind_df['phi'] = wind_df.apply(lambda row: np.arctan(row['horz_length']),axis=1)
        wind_df['shift_CH4'] = wind_df.CH4.shift(periods = int(float(shift)))
        wind_df['raw_CH4'] = wind_df.apply(lambda row: row['BCH4'],axis=1)
        wind_df['BCH4']= wind_df.loc[:,['shift_CH4']]
        wind_df['CH4']= wind_df.loc[:,['shift_CH4']]
        wind_df['TCH4']= wind_df.loc[:,['shift_CH4']]
        
        wind_df2 = wind_df[wind_df.CH4.notnull()]
        wind_df3 = wind_df2.drop(['QUADRANT', 'secnan','prev_LAT','next_LAT','prev_LONG','next_LONG','prev_TIME','next_TIME','distance','timediff','uncor_theta','CH4'],axis = 1)
        wind_df3['CH4'] = wind_df3.loc[:,'shift_CH4']
        wind_df3 = wind_df3.drop(['shift_CH4'],axis = 1)
        wind_df3 = wind_df3.loc[:,['DATE','TIME','SECONDS','NANOSECONDS','VELOCITY','U','V','W','BCH4','BRSSI','TCH4','TRSSI','PRESS_MBAR','INLET' \
                                   , 'TEMPC','CH4','H20','C2H6','R','C2C1','BATTV','POWMV','CURRMA','SOCPER','LAT','LONG','bearing','U_cor', \
                                   'horz_length','adj_theta','totalWind','phi','raw_CH4']]
        wind_df4 = wind_df3.loc[wind_df3.totalWind.notnull(),:]
        #firstTime = wind_df3.SECONDS.min() + 60 *(initialTimeBack)                           
        #wind_df4 = wind_df3.loc[wind_df3.SECONDS > firstTime,:]                      
       # wind_df3.to_csv(fnOutTemp,index=False)
        
        if bFirst:
            wind_df4.to_csv(fnOut,index=False)
        elif not bFirst:
            norm = pd.read_csv(fnOut)
            pd.concat([norm,wind_df4]).sort_values(by='SECONDS').reset_index(drop=True).to_csv(fnOut,index=False)
        os.remove(fnOutTemp)
        return True
    except ValueError:
        return False
    
def IdentifyPeaks(xCar, xDate, xDir, xFilename,outDir,processedFileLoc,threshold = '.1',xTimeThreshold = '5.0',engineer='False'):
    import csv, numpy    
    import shutil 
    try:
        xABThreshold = float(threshold)

        #xABThreshold = 0.1                 # above baseline threshold above the mean value
        xDistThreshold = 160.0                 # find the maximum CH4 reading of observations within street segments of this grouping distance in meters
        xSDF = 4                    # multiplier times standard deviation for floating baseline added to mean
        xB = 1020       # the number of records that constitutes the floating baseline time -- 7200 = 1 hour (assuming average of 0.5 seconds per record)
        #xB = 300 #five min?
        xTimeThreshold = float(xTimeThreshold)
        
        fn = xDir + "/" + xFilename      #set raw text file to read in
        fnOut = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".csv"       #set CSV format output for observed peaks for a given car, day, city
        fnOutjson = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".JSON"       #set CSV format output for observed peaks for a given car, day, city

        fnShape = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".shp"
        fnLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".log"       #set CSV output for observed peaks for a given car, day, city
        pkLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + "_info.csv"       #set CSV output for observed peaks for a given car, day, city
        
        infOut = processedFileLoc + xCar + "_" + xDate.replace("-","") + "_info.csv"
        print(str(outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + "_info.csv"))
        
        fLog = open(fnLog, 'w')        
        shutil.copy(infOut,pkLog)


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
                #x1.append(float(str(row[fEpochTime]) + '.' + str(row[fNanoSeconds])));
                if engineer:
                    x1.append(float(str(row[fEpochTime])));# + '.' + str(row[fNanoSeconds])));
                elif not engineer:
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

                xCH4Mean = numpy.percentile(aCH4[botBound:topBound],50)
               # xCH4SD = numpy.std(aCH4[botBound:topBound])
            else:
                xCH4Mean = numpy.percentile(aCH4[0:(count-2)],50)
                #xCH4SD = numpy.std(aCH4[0:(count-2)])
            xThreshold = xCH4Mean + (xCH4Mean * xABThreshold)
            
            if (aCH4[i] > xThreshold):
                lstCH4_AB.append(i)
                aMean[i] = xCH4Mean.copy()    #insert mean + SD as upper quartile CH4 value into the array to later retreive into the peak calculation
                aThreshold[i] = xThreshold.copy()

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
        
        
        # EXTRA INFO
        identPks = pd.read_csv(fnOut)
        identPks['OB_CH4_AB'] = identPks.loc[:,'OB_CH4'].sub(identPks.loc[:,'OB_CH4_BASELINE'], axis = 0) 
        locpks = weightedLoc(identPks,'OB_LAT','OB_LON','OP_NUM','OB_CH4_AB').rename(columns = {'OB_LAT':'pk_LAT','OB_LON':'pk_LON'})

        #maxch4 = identPks.groupby('OP_NUM',as_index = False).OB_CH4_AB.max().rename(columns = {'OB_CH4_AB':'pk_maxCH4_AB'})
        #identPksWt = pd.merge(identPks,maxch4,on = ['OP_NUM'])

        geo = [Point(xy) for xy in zip(locpks[('pk_LON')], locpks[('pk_LAT')])]
        crs = {'init': 'epsg:4326'}
        savegeo = gpd.GeoDataFrame(locpks, crs=crs, geometry=geo)
        savegeo = savegeo.to_crs(epsg=32610).copy()
        savegeo.to_file(fnOutjson, driver="GeoJSON")
        
        geo = [Point(xy) for xy in zip(identPks[('OB_LON')], identPks[('OB_LAT')])]
        crs = {'init': 'epsg:4326'}
        savegeo = gpd.GeoDataFrame(identPks, crs=crs, geometry=geo)
        savegeo = savegeo.to_crs(epsg=32610).copy()
        savegeo.to_file(fnOutjson, driver="GeoJSON")

        
        if truecount > 0:
            #arcpy.MakeXYEventLayer_management(fnOut,"LON","LAT",xCar + "L","GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;0.001;0.001;IsHighPrecision","#")
            #arcpy.FeatureToPoint_management(xCar + "L",fnShape,"CENTROID")
            #arcpy.Delete_management(xCar+"L")
            return True
        
    except ValueError:
            print ("Error in Identify Peaks")
            return False
        
def ProcessRawDataEngOLD( xCar, xDate, xDir, xFilename, bFirst, gZIP, xOut,initialTimeBack):
    from datetime import datetime
    #import os
    import gzip
    import csv
    from math import floor
    try:
        
        #xOutDir = xDir 
        
        xMaxCarSpeed = 45.0 / 2.23694     # assumes 45 mph as max, convert to meters per second 
        xMinCarSpeed = 2.0 / 2.23694      # minimum 2 miles per hour
        xMinRSSI = 50  #if RSSI is below this we don't like it
        

        # reading in the data with specific headers
        #          0     1    2    3       4           5    6       7        8        9          10                 11              12           13            14      15      16      17        18         19         20         21         22         23        24   25  26       27           28       29           30       31       32       33  34        35   36   37  38   39       40       41   42       43   44   45   46   47   48   49   50   51     52     53     54
        #sHeader = "Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude"
        #sHeader = "Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude"
        #sHeader = "Time Stamp,Inlet Number,P (mbars),T0 (degC),T5 (degC), Laser PID Readout,Det PID Readout,win0Fit0,win0Fit1,win0Fit3,win1Fit4,win0Fit5,win0Fit6,win0Fit7,win0Fit8,win0Fit9,win1Fit0,win1Fit1,win1Fit2,win1Fit3,win1Fit4,win1Fit5,win1Fit6,Det Bkgd,Ramp Ampl,CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Battery T (degC),FET T (degC),GPS Time,Latitude,Longitude"
        sHeader = "Time Stamp,Inlet Number,P (mbars),T0 (degC),T5 (degC),Laser PID Readout,Det PID Readout,win0Fit0,win0Fit1,win0Fit2,win0Fit3,win0Fit4,win0Fit5,win0Fit6,win0Fit7,win0Fit8,win0Fit9,win1Fit0,win1Fit1,win1Fit2,win1Fit3,win1Fit4,win1Fit5,win1Fit6,Det Bkgd,Ramp Ampl,CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Battery T (degC),FET T (degC),GPS Time,Latitude,Longitude"
        sOutHeader = "DATE,TIME,SECONDS,NANOSECONDS,VELOCITY,U,V,W,BCH4,BRSSI,TCH4,TRSSI,PRESS_MBAR,INLET,TEMPC,CH4,H20,C2H6,R,C2C1,BATTV,POWMV,CURRMA,SOCPER,LAT,LONG\n"
        
        headerNames = sHeader.split(',')
        GPS_loc = 37
        
        
        infoHeader = "FILENAME\n"
        # somehow gZIP is indicating if  it is the first file name (I think if it is 0 then it is the first file)
        if gZIP == 0:
            f = gzip.open(xDir + "/" + xFilename, 'r') #if in python 3, change this to "r" or just "b" can't remember but something about a bit not a string
        else:
            f = open(xDir + "/" + xFilename, 'r')
            f = open(xDir  + xFilename, 'r')

        
        # process    
        #if first time on this car/date, then write header out
        
        xdat = str('20') + xFilename[11:17]
        
        #fnOut = xOutDir + xCar + "_" + xDate.replace("-", "") + "_dat.csv"       #set CSV output for raw data
        #fnLog = xOutDir + xCar + "_" + xDate.replace("-", "") + "_log.csv"       #output for logfile
        
        fnOut = xOut  + xCar + "_" + xdat + "_dat.csv"       #set CSV output for raw data
        fnLog =  xOut  + xCar + "_" + xdat + "_log.csv"       #output for logfile
        infOut = xOut + xCar + "_" + xdat + "_info.csv"
        
        firsttime = int(float(open(xDir + xFilename).readlines().pop(1).split(',')[37][:-4]))
        
        
        if bFirst:
            fOut = open(fnOut, 'w')
            fOut.write(sOutHeader)
            fLog = open(fnLog, 'w')
            infOut = open(infOut,'w')
            infOut.write(infoHeader)
            print ("fnLog: "+fnOut) 
        if not bFirst:
            fOut = open(fnOut, 'a')
            fLog = open(fnLog, 'a')
            infOut = open(infOut,'a')
        
        #read all lines
        xCntObs = -1
        xCntGoodValues = 0
        rownum = 0 
        for row in f:
            woo = row
            #print(row)
            bGood = True
            if xCntObs < 0:
                bGood = False
                xCntObs += 1
                #lstSfirst = woo.split(',')
                #gpstimefirst = lstSfirst[GPS_loc]
                #firstseconds = floor(float(gpstimefirst))

            if bGood:
                #lstS = row.split(",")
                lstS = row.split(',')
                gpstime = lstS[GPS_loc]
                dtime = lstS[0]
                dt = lstS[1]
                time_dt = lstS[2]
                epoch = lstS[3]
                #nano = lstS[4]
                
                gps_time = lstS[37]
                dateob = datetime.fromtimestamp(int(gps_time[:-4]))
                nano= gps_time[-4:]
                
                #dateob = datetime(int(dt[0:4]),int(dt[5:7]),int(dt[8:10]),int(time_dt[0:2]),int(time_dt[3:5]),int(time_dt[6:8]),int(float(nano)*1e-9))
                
                dtime  = int(dateob.strftime('%Y%m%d%H%M%S'))
                #Date = dateob.strftime('%m%/%d/%Y')
                Date = dateob.strftime('%Y-%m-%d')

                GPS_Time = dateob.strftime('%H%:%M:%S')
                seconds = floor(float(gpstime))
                nano = dateob.strftime('%f')
                
                #dateob = datetime(int(dtime[6:10]),int(dtime[0:2]),int(dtime[3:5]),int(dtime[11:13]),int(dtime[14:16]),int(dtime[17:19]),int(float(dtime[19:23])*1000000))
                #epoch = dateob.strftime('%s.%f')
  
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
                import sys
                if sys.platform.startswith('win'):
                    ## DATE, TIME, SECONDS,NANOSECONDS
                    csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(dateob.strftime('%H:%M:%S'))  + ',' + str(float(pd.to_numeric(dateob.strftime('%S.%f')))) + ',' + str(pd.to_numeric(dateob.strftime('%f')) *1000) + str(',')
                    ## VELOCITY, U,V,W,BCH4,BRSSI,TCH4
                    csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(lstS[26]) + ',' + str('0') + ','+  str(lstS[26]) + ','
                    ## TRSSI, PRESS_MBAR, INLET, TEMPC, CH4, H20,C2H6
                    csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(lstS[26]) + ',' + str(lstS[27]) +',' +  str(lstS[28]) + ','
                    # R, C2C1, BATTV, POWMV,CURRMA, SOCPER,LAT,LONG
                    csvWrite += str(lstS[29]) + ',' + str(lstS[30]) + ',' + str(lstS[31]) + ',' + str(lstS[32]) + ','+ str(lstS[33]) + ',' + str(lstS[34]) + ',' + str(lstS[38]) + str(',') + str(lstS[39])
                
# =============================================================================
#                 if not sys.platform.startswith('win'):
#                     ## DATE, TIME, SECONDS,NANOSECONDS
#                     csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(dateob.strftime('%H:%M:%S'))  + ',' + str((int(floor(pd.to_numeric(dateob.strftime('%s.%f')))))) + ',' + str((pd.to_numeric(dateob.strftime('%f')) *1000)) + str(',')
#                     ## VELOCITY, U,V,W,BCH4,BRSSI,TCH4
#                     csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(lstS[26]) + ',' + str('0') + ','+  str(lstS[26]) + ','
#                     ## TRSSI, PRESS_MBAR, INLET, TEMPC, CH4, H20,C2H6
#                     csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(lstS[26]) + ',' + str(lstS[27]) +',' +  str(lstS[28]) + ','
#                     # R, C2C1, BATTV, POWMV,CURRMA, SOCPER,LAT,LONG
#                     csvWrite += str(lstS[29]) + ',' + str(lstS[30]) + ',' + str(lstS[31]) + ',' + str(lstS[32]) + ','+ str(lstS[33]) + ',' + str(lstS[34]) + ',' + str(lstS[38]) + str(',') + str(lstS[39][:-1]) + str('\n')                   
#                 #fOut.write('\n')
#                 fOut.write(csvWrite)
#                 #fOut.write('\n')
#                 
# =============================================================================
                if not sys.platform.startswith('win'):
                    ## DATE, TIME, SECONDS,NANOSECONDS
                    csvWrite = str(Date) + ',' + str(GPS_Time)  + ',' + str(seconds) + ',' + str(nano) + str(',')
                    ## VELOCITY, U,V,W,BCH4,BRSSI,TCH4
                    csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(lstS[26]) + ',' + str('0') + ','+  str(lstS[26]) + ','
                    ## TRSSI, PRESS_MBAR, INLET, TEMPC, CH4, H20,C2H6
                    csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(lstS[26]) + ',' + str(lstS[27]) +',' +  str(lstS[28]) + ','
                    # R, C2C1, BATTV, POWMV,CURRMA, SOCPER,LAT,LONG
                    csvWrite += str(lstS[29]) + ',' + str(lstS[30]) + ',' + str(lstS[31]) + ',' + str(lstS[32]) + ','+ str(lstS[33]) + ',' + str(lstS[34]) + ',' + str(lstS[38]) + str(',') + str(lstS[39])    
                #fOut.write('\n')
                
                #print( seconds > firstseconds + (60*5))
                if seconds >= (firsttime + (60*5)):
                    fOut.write(csvWrite)
                
                
                del(csvWrite)
#                xCntGoodValues += 1
                

            xCntObs += 1

        #sOut = str(gZIP) + "," + str(f) + "," + str(xCntObs) + "," + str(xCntGoodValues) + "\n"
        #fLog.write(sOut)

        infOut.write(str(xFilename)+'\n')

        fOut.close()
        fLog.close()
        infOut.close()
        
        #xDate = dateob.strftime("%Y%m%d")
        
        #newfnOut = xOutDir + xCar + "_" + xDate + "_dat.csv"       #set CSV output for raw data
        #newfnLog = xOutDir + xCar + "_" + xDate + "_log.csv"  
        

        print (xCar + "\t" + xdat + "\t" + fnOut[-22:] + "\t" + str(xCntObs) + "\t" + str(xCntGoodValues) + "\t" + str(gZIP))
        from numpy import pi
        import numpy as np
        def calcVel(timediff,distance):
            if timediff == 0:
                return(0)
            elif timediff != 0:
                return(distance/timediff)
        radians = False
        wind_df = pd.read_csv(fnOut)        
        wind_df['QUADRANT'] = wind_df.apply(lambda row: getQuad(row['U'],row['V']),axis=1)
        wind_df['secnan'] = wind_df.apply(lambda row: row['SECONDS'],axis=1) # + row['NANOSECONDS']*1e-9,axis=1)
        wind_df['prev_LAT'] = wind_df.LAT.shift(periods = 1)
        wind_df['next_LAT'] = wind_df.LAT.shift(periods = -1)
        wind_df['prev_LONG'] = wind_df.LONG.shift(periods = 1)
        wind_df['next_LONG'] = wind_df.LONG.shift(periods = -1)
        wind_df['prev_TIME'] = wind_df.secnan.shift(periods = 1)
        wind_df['next_TIME'] = wind_df.secnan.shift(periods = -1)
        wind_df['distance'] = wind_df.apply(lambda row: haversine(row['prev_LAT'],row['prev_LONG'],row['next_LAT'],row['next_LONG']),axis=1)
        wind_df['bearing'] = wind_df.apply(lambda row: calcBearing(row['prev_LAT'],row['next_LAT'],row['prev_LONG'],row['next_LONG'],radians),axis=1)
        wind_df['timediff'] = wind_df.apply(lambda row: row['next_TIME'] - row['prev_TIME'],axis = 1)
        wind_df['VELOCITY'] = wind_df.apply(lambda row:calcVel(row['timediff'],row['distance']),axis=1)
        wind_df['U_cor'] = wind_df.apply(lambda row:row['U'] + row['VELOCITY'],axis = 1)
        wind_df['horz_length'] = wind_df.apply(lambda row: np.sqrt(row['U_cor']**2 + row['V']**2),axis=1)
        wind_df['uncor_theta'] = wind_df.apply(lambda row :calcBearing(row['U_cor'],row['V'],row['QUADRANT'],row['horz_length'],radians),axis = 1)
        wind_df['adj_theta'] = wind_df.apply(lambda row: (row['uncor_theta'] + row['bearing'])%360,axis =1)
        wind_df['totalWind'] = wind_df.apply(lambda row: np.sqrt(row['horz_length']**2 + row['W']**2),axis = 1)
        wind_df['phi'] = wind_df.apply(lambda row: np.arctan(row['horz_length']),axis=1)
        wind_df['shift_CH4'] = wind_df.CH4.shift(periods = shift)
        wind_df['raw_CH4'] = wind_df.apply(lambda row: row['BCH4'],axis=1)
        wind_df['BCH4']= wind_df.loc[:,['shift_CH4']]
        wind_df['CH4']= wind_df.loc[:,['shift_CH4']]
        wind_df['TCH4']= wind_df.loc[:,['shift_CH4']]
        
        wind_df2 = wind_df[wind_df.CH4.notnull()]
        wind_df3 = wind_df2.drop(['QUADRANT', 'secnan','prev_LAT','next_LAT','prev_LONG','next_LONG','prev_TIME','next_TIME','distance','timediff','uncor_theta','CH4'],axis = 1)
        wind_df3['CH4'] = wind_df3.loc[:,'shift_CH4']
        wind_df3 = wind_df3.drop(['shift_CH4'],axis = 1)
        wind_df3 = wind_df3.loc[:,['DATE','TIME','SECONDS','NANOSECONDS','VELOCITY','U','V','W','BCH4','BRSSI','TCH4','TRSSI','PRESS_MBAR','INLET' \
                                   , 'TEMPC','CH4','H20','C2H6','R','C2C1','BATTV','POWMV','CURRMA','SOCPER','LAT','LONG','bearing','U_cor', \
                                   'horz_length','adj_theta','totalWind','phi','raw_CH4']]
        #firstTime = wind_df3.SECONDS.min() + 60 *(initialTimeBack)                           
        #wind_df4 = wind_df3.loc[wind_df3.SECONDS > firstTime,:]                      
        wind_df3.to_csv(fnOut,index=False)
        

        return True
    except ValueError:
        return False

    
## POTENTIALLY GOOD ONE    
def IdentifyPeaksOLD( xCar, xDate, xDir, xFilename,outDir,processedFileLoc,threshold = '.1',xTimeThreshold = '5.0'):
    import csv, numpy    
    import shutil 
    try:
        xABThreshold = float(threshold)

        #xABThreshold = 0.1                 # above baseline threshold above the mean value
        xDistThreshold = 160.0                 # find the maximum CH4 reading of observations within street segments of this grouping distance in meters
        xSDF = 4                    # multiplier times standard deviation for floating baseline added to mean
        xB = 1020       # the number of records that constitutes the floating baseline time -- 7200 = 1 hour (assuming average of 0.5 seconds per record)
        #xB = 300 #five min?
        xTimeThreshold = float(xTimeThreshold)
        
        #fn = xDir + "/" + xFilename      #set raw text file to read in
        fn =  xDir  + xFilename
        fnOut = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".csv"       #set CSV format output for observed peaks for a given car, day, city
        fnShape = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".shp"
        fnLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".log"       #set CSV output for observed peaks for a given car, day, city
        pkLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + "_info.csv"       #set CSV output for observed peaks for a given car, day, city
        
        infOut = processedFileLoc + xCar + "_" + xDate.replace("-","") + "_info.csv"
        print(str(outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + "_info.csv"))
        
        fLog = open(fnLog, 'w')        
        shutil.copy(infOut,pkLog)


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
        
        


        #read data in from text file and extract desired fields into a list, padding with 5 minute and hourly average
        x1 = []; x2 = []; x3 = []; x4 = []; x5 = []; x6 = []; x7 = []; x8 = []
        
        count = -1
        with open(fn, 'r') as f:
            t = csv.reader(f)
            for row in t:
                woo = row
                if count < 0:
                    count += 1
                    continue
                
                datet= row[fDate].replace("-","")+row[fTime].replace(":","")

                #x1.append(float(epoch)); 
                x1.append(float(str(row[fEpochTime]))); #+ '.' + str(row[fNanoSeconds])));
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

                xCH4Mean = numpy.percentile(aCH4[botBound:topBound],50)
               # xCH4SD = numpy.std(aCH4[botBound:topBound])
            else:
                xCH4Mean = numpy.percentile(aCH4[0:(count-2)],50)
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

def makeGEO(df,lat,lon):
    geo = [Point(xy) for xy in zip(df[(lon)], df[(lat)])]
    return (geo)

def makeGPD(df,lat,lon,cps = 'epsg:4326'):
    gdf = gpd.GeoDataFrame(df, crs=cps, geometry=makeGEO(df,lat,lon))
    return(gdf)
    
        
def filterPeak(xCar,xDate,xDir,xFilename, outFolder,whichpass = 0):
    import pandas as pd #
    import geopandas as gpd
    import shutil
    from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry

    file_loc = xDir + xFilename
    new_loc = outFolder + "Filtered" + xFilename
    new_loc_json = new_loc[:-3] + 'json'

    oldInfo = xDir + 'Peaks_' + xCar + "_" + xDate.replace("-","") + "_info.csv"
    newInfo = outFolder + 'FilteredPeaks_' + xCar + "_" + xDate.replace("-","") + "_info.csv"
    
    shutil.copy(oldInfo,newInfo)


    datFram = pd.read_csv(file_loc)
    datFram_cent =  datFram.loc[:,:]  
    datFram_cent['CH4_AB'] = datFram.loc[:,'CH4'].sub(datFram.loc[:,'CH4_BASELINE'], axis = 0) 
    
 
    ### MAXCH4 is a df with the max methane (above baseline) in the given observed peak
    maxch4 = datFram_cent.groupby('PEAK_NUM',as_index = False).CH4_AB.max().rename(columns = {'CH4_AB':'pk_maxCH4_AB'})
        
    ### FINDING WEIGHTED LOCATION OF THE OP, BY THE ABOVE BASELINE CH4 LEVEL
    #wtloc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB')
    #datFram_wtLoca =  wtloc.copy()
    #datFram_wtLoc = datFram_wtLoca.rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'})
    
    datFram_wtLoc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB').loc[:,:].rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'})
    #datFram_wtLoc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB').rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'}).copy()
    datFram_wtLocMax = pd.merge(datFram_wtLoc,maxch4,on = ['PEAK_NUM'])
    
    pass_info = datFram.copy()
    
    ## MIGHT NEED TO CHANGE BACK
    geometry_temp = [Point(xy) for xy in zip(datFram['LON'], datFram['LAT'])]
    crs = {'init': 'epsg:4326'}
        
        #geometry is the point of the lat/lon
    gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)
    
    #gdf_buff = makeGPD(datFram,'LON','LAT')
    gdf_buff = gdf_buff.to_crs(epsg=32610)
    gdf_buff['geometry'] = gdf_buff.loc[:,'geometry'].buffer(30) 
    
    pass_info_new = pass_info.rename(columns={"geometry": 'pk_geo'})
    
    gdf_tog = pd.merge(gdf_buff,pass_info_new,on = ['PEAK_NUM', 'EPOCHSTART', 'EPOCH', 'DATETIME', 'CH4', 'LON', 'LAT',
       'CH4_BASELINE', 'CH4_THRESHOLD', 'PEAK_DIST_M', 'PEAK_CH4', 'TCH4',
       'PERIOD5MIN'])
    
    gdf_bind_pks = gdf_tog.dissolve(by = 'PEAK_NUM',as_index=False).loc[:,['PEAK_NUM','geometry']]
                  
    if gdf_bind_pks.shape[0] > 1:
        data_overlap = gpd.GeoDataFrame(crs=gdf_bind_pks.crs)
        data_temp = gdf_bind_pks.copy()
        for index, row in data_temp.iterrows():
            data_temp1=data_temp.loc[data_temp.PEAK_NUM!=row.PEAK_NUM,]
            # check if intersection occured
            overlaps=data_temp1[data_temp1.geometry.overlaps(row.geometry)]['PEAK_NUM'].tolist()
            if len(overlaps)>0:
                
                # compare the area with threshold 
                for y in overlaps:
                    temp_area=gpd.overlay(data_temp.loc[data_temp.PEAK_NUM==y,],data_temp.loc[data_temp.PEAK_NUM==row.PEAK_NUM,],how='intersection')
                    temp_area=temp_area.loc[temp_area.geometry.area>=0.001]
                    if temp_area.shape[0]>0:
                        temp_union = gpd.overlay(data_temp.loc[data_temp.PEAK_NUM==y,],data_temp.loc[data_temp.PEAK_NUM==row.PEAK_NUM,],how='union')
                        data_overlap=gpd.GeoDataFrame(pd.concat([temp_union,data_overlap],ignore_index=True),crs=data_temp.crs)
        if data_overlap.size > 0: 
                
                firstnull2 = data_overlap.loc[data_overlap.PEAK_NUM_1.isnull(),:]
                firstnull = firstnull2.copy()
                firstnull.loc[:,'PEAK_NUM_1'] = firstnull2.loc[:,'PEAK_NUM_2']
                
                secnull2 = data_overlap.loc[data_overlap.PEAK_NUM_2.isnull(),:]
                
                secnull = secnull2.copy()
                secnull.loc[:,'PEAK_NUM_2'] = secnull2.loc[:,'PEAK_NUM_1']
                
                withoutNA = data_overlap.copy().dropna()
                allTog = pd.concat([firstnull,secnull,withoutNA]).reset_index().copy()
                
                over = allTog.copy()
                over['sorted']=over.apply(lambda y: sorted([y['PEAK_NUM_1'],y['PEAK_NUM_2']]),axis=1)
                over['sorted']=over.sorted.apply(lambda y: ''.join(y))
                over = over.drop_duplicates('sorted')
                over['combined']= [list(x) for x in list(over.loc[:,['PEAK_NUM_1','PEAK_NUM_2']].to_numpy())]
                over['VER_NUM'] = over.apply(lambda y: y.combined,axis=1)
                over['min_val']=over.apply(lambda y: min(y.combined),axis=1)
                over=over.reset_index().loc[:,['PEAK_NUM_1','PEAK_NUM_2','geometry','combined','min_val']]
                
                overcop = over.copy()
                overcop.loc[:,'recombine'] = overcop.loc[:,'combined']
                
                for index, row in overcop.iterrows():
                    united = row.recombine
                    for index2, row2 in overcop.iterrows():
                        united_temp = unIfInt(united,row2.recombine)
                        if united_temp != None:
                            united = united_temp
                    overcop.at[index, 'recombine']  = united.copy()
                    del(united)   
                overcop['recombine']= overcop.apply(lambda y: sorted(y.recombine),axis=1).copy()
                overcop['min_read'] = overcop.apply(lambda y: min(y.recombine),axis=1).copy()
                newOverlap = overcop.dissolve(by='min_read',as_index=False).loc[:,['min_read','geometry','recombine']].copy()

  
                combined = gdf_bind_pks.copy()
                combined['recombine'] = [list(x) for x in list(combined.loc[:,['PEAK_NUM']].to_numpy())]
                combined['numtimes'] = 1
                combined['newgeo'] = combined.loc[:,'geometry']
                combined['min_read'] = combined.loc[:,"PEAK_NUM"]
                for index,row in combined.iterrows():
                    for index2,row2 in newOverlap.iterrows():
                        if row.PEAK_NUM in row2.recombine:
                            combined.at[index, 'recombine']  = row2.recombine.copy()
                            combined.at[index, 'newgeo']  = row2.copy().geometry
                            combined.at[index,'min_read'] = row2.copy().min_read
                            
                combined['numtimes'] = combined.apply(lambda y: len(y.recombine),axis = 1).copy()
                combined_reduced = combined.loc[:,['PEAK_NUM','newgeo','recombine','numtimes','min_read']]
                gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['PEAK_NUM']).copy()
                gdf_pass_pks['verified'] = gdf_pass_pks.apply(lambda y: (True if y.numtimes > 1 else False),axis=1 ).copy()
        if data_overlap.size == 0:
           gdf_pass_pks = gdf_bind_pks.copy()
           gdf_pass_pks['min_read']= gdf_pass_pks.loc[:,'PEAK_NUM']
           gdf_pass_pks['numtimes'] = 1
           gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:,'geometry']
           gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['PEAK_NUM']].to_numpy())].copy()
           gdf_pass_pks['verified'] = False
           gdf_pass_pks['oldgeo'] = gdf_pass_pks.loc[:,'geometry']
           gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
    if gdf_bind_pks.shape[0] == 1:
        gdf_pass_pks = gdf_bind_pks.copy()
        gdf_pass_pks['min_read']= gdf_pass_pks.loc[:,'PEAK_NUM']
        gdf_pass_pks['numtimes'] = 1
        gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:,'geometry']
        gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['PEAK_NUM']].to_numpy())].copy()
        gdf_pass_pks['verified'] = False
        epdat = pass_info.loc[:,['PEAK_NUM','EPOCHSTART']]
        gdf_pass_pks = pd.merge(gdf_pass_pks,epdat,on = ['PEAK_NUM']).copy()
    gdf_pass_pks['oldgeo'] = gdf_pass_pks.loc[:,"geometry"]
    gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
    del(gdf_pass_pks['newgeo'])
    gdf_pass_pks['pass'] = whichpass
    gdf_tot = pd.merge(gdf_pass_pks,datFram_wtLocMax,on = ['PEAK_NUM']).copy()
    ## condense by peak_num
    gdfcop = gdf_tot.loc[:,['PEAK_NUM','geometry','min_read','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    gdfcop = gdfcop.to_crs(epsg=32610).copy()
    gdfcop.to_file(new_loc_json, driver="GeoJSON")


    gdf_tot.to_csv(new_loc, index = False) 
    return(gdf_tot)
        

    
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
    import geopandas as gpd

    first_geo = firstgrp.copy().geometry
    crs = {'init': 'epsg:4326'}
    firstgrp = gpd.GeoDataFrame(firstgrp.copy(),crs = crs,geometry = first_geo)
    
    first_pks = firstgrp.loc[:,['PEAK_NUM','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    sec_pks = secondgrp.loc[:,['PEAK_NUM','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    tot_pks = pd.concat([first_pks,sec_pks])
    
    first_dis = firstgrp.dissolve(by='min_read',as_index=False)[['min_read','geometry','recombine','verified','pass']].copy()
    sec_dis = secondgrp.dissolve(by='min_read',as_index=False)[['min_read','geometry','recombine','verified','pass']].copy()
    gdf_bind_pks = pd.concat([first_dis,sec_dis])[['min_read','geometry','recombine']].copy()
    gdf_tog = pd.concat([firstgrp.drop(['pk_LAT', 'pk_LON','pk_maxCH4_AB'], axis=1),secondgrp.drop(['pk_LAT', 'pk_LON','pk_maxCH4_AB'], axis=1)]).copy()

    gdf_bind_pks['prev_read'] = gdf_bind_pks.loc[:,'min_read']
    gdf_tog['prev_read'] = gdf_tog.loc[:,'min_read']
    if gdf_bind_pks.shape[0] > 1:
        data_overlap = gpd.GeoDataFrame(crs=gdf_bind_pks.crs).copy()
        data_temp = gdf_bind_pks.copy()
        for index, row in data_temp.iterrows():
            data_temp1=data_temp.loc[data_temp.min_read!=row.min_read,]
            # check if intersection occured
            overlaps=data_temp1[data_temp1.geometry.overlaps(row.geometry)]['min_read'].tolist()
            if len(overlaps)>0:
                # compare the area with threshold 
                for y in overlaps:
                    temp_area=gpd.overlay(data_temp.loc[data_temp.min_read==y,],data_temp.loc[data_temp.min_read==row.min_read,],how='intersection')
                    temp_area=temp_area.loc[temp_area.geometry.area>=0]
                    #temp_union = gpd.overlay(data_temp.loc[data_temp.PEAK_NUM==y,],data_temp.loc[data_temp.PEAK_NUM==row.PEAK_NUM,],how='union')
                    if temp_area.shape[0]>0:
                        temp_union = gpd.overlay(data_temp.loc[data_temp.min_read==y,],data_temp.loc[data_temp.min_read==row.min_read,],how='union')
                        data_overlap=gpd.GeoDataFrame(pd.concat([temp_union,data_overlap],ignore_index=True),crs=data_temp.crs)    
    
        if data_overlap.size != 0: 
            firstnull = data_overlap[data_overlap.min_read_1.isnull()].copy()
            firstnull.loc[:,'min_read_1'] = firstnull.loc[:,'min_read_2']
                
            secnull = data_overlap[data_overlap.min_read_2.isnull()].copy()
            secnull['min_read_2'] = secnull['min_read_1'].copy()
                
            withoutNA = data_overlap.dropna().copy()
            allTog = pd.concat([firstnull,secnull,withoutNA]).reset_index().copy()
            
            
            over = allTog.copy()
            over['sorted']=over.apply(lambda y: sorted([y['min_read_1'],y['min_read_2']]),axis=1).copy()
            over['sorted']=over.sorted.apply(lambda y: ''.join(y)).copy()
            over = over.drop_duplicates('sorted').copy()
            over['combined']= [list(x) for x in list(over.loc[:,['min_read_1','min_read_2']].to_numpy())].copy()
            over['VER_NUM'] = over.apply(lambda y: y.combined,axis=1).copy()
            over['min_val']=over.apply(lambda y: min(y.combined),axis=1).copy()
            over=over.reset_index()[['min_read_1','min_read_2','geometry','combined','min_val']].copy()
                
            overcop = over.copy()
            overcop['recombine'] = overcop.combined.copy()
            
            for index,row in overcop.iterrows():
                #rowwoo = row
                
                ## fixthis
                
                first_thing = first_dis[first_dis['min_read']== row.min_read_1].loc[:,['recombine']]
                #first_thing = first_dis.loc[first_dis.min_read == row.min_read_1,'recombine']
                #first_thing = first_dis.loc[:,first_dis.min_read == row.min_read_1].loc[:,'recombine']

                #first_thing = firstpass[firstpass['min_read']== row.min_read_1].loc[:,['recombine']]
                firstcomb = first_thing.recombine.explode().copy()
                first_list = firstcomb.reset_index().recombine.copy()
                #first_list = firstpass.loc[:,firstpass['min_read']== row.min_read_1].recombine
                #first_list.explode()
                #first_list = firstpass.loc[:,firstpass['min_read']== row.min_read_1].recombine.explode()[0]
               # second_list = secondpass[secondpass['min_read']== row.min_read_2].recombine
               
               
               ## fix this 
               
                sec_thing = sec_dis[sec_dis['min_read']== row.min_read_1].loc[:,['recombine']]

                #sec_thing = secondpass[secondpass['min_read']== row.min_read_2].loc[:,['recombine']]
                seccomb = sec_thing.recombine.explode()
                sec_list = seccomb.reset_index().recombine
               
                firstdf = pd.DataFrame(first_list)
                secdf = pd.DataFrame(sec_list)
    
                
                tot_df = pd.concat([firstdf,secdf])
                tot_list = tot_df.recombine.tolist()
                
                overcop.at[index, 'recombine']  = tot_list.copy()




    ## this recombines the lists together to have the combined entries together?
       
            overcop['recombine']= overcop.apply(lambda y: sorted(y.recombine),axis=1).copy()
            overcop['min_read'] = overcop.apply(lambda y: min(y.recombine),axis=1).copy()
            newOverlap = overcop.dissolve(by='min_read',as_index=False).loc[:,['min_read','geometry','recombine']].copy()
        
      
            combined = gdf_bind_pks.copy()
            #combined['recombine'] = [list(x) for x in list(combined.loc[:,['min_read']].to_numpy())]
            combined['numtimes'] = 1
            combined['newgeo'] = combined.copy().geometry
            combined['oldgeo'] = combined.copy().geometry

            #combined['min_read'] = combined.min_rea
            combined = combined.reset_index()
            for index,row in combined.iterrows():
                for index2,row2 in newOverlap.iterrows():
                    if row.min_read in row2.recombine:                      
                        combined.at[index, 'recombine']  = row2.recombine.copy()
                        combined.at[index, 'min_read']  = row2.copy().min_read
        
            combined['numtimes'] = combined.copy().apply(lambda y: len(y.recombine),axis = 1).copy()
            combined['geometry'] = combined.copy().newgeo
            
            del(combined['newgeo'])
            combined_reduced = combined[['min_read','geometry','oldgeo','recombine','numtimes','prev_read']].copy()
            #gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['min_read'])
            gdf_tog = gdf_tog.drop('min_read',axis=1)
            gdf_tog = gdf_tog.drop('numtimes',axis=1)
            gdf_tog = gdf_tog.drop('recombine',axis=1)
           # gdf_tog = gdf_tog.drop('newgeo',axis=1)
            gdf_tog['firstgeo'] = gdf_tog.copy().oldgeo
            gdf_tog['secondgeo'] = gdf_tog.copy().geometry
            
            del(gdf_tog['geometry'])
            #gdf_tog.drop(columns=['geometry'])
            del(gdf_tog['oldgeo'])

                                   
            gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['prev_read']).copy()

            gdf_pass_pks['verified'] = gdf_pass_pks.copy().apply(lambda y: (True if y.numtimes > 1 else False),axis=1 )
        if data_overlap.size == 0:
            gdf_pass_pks = gdf_tog.copy()
    ## didnt adress if the bind shape was only size only 1
        gdf_tot_pks = pd.merge(gdf_pass_pks,tot_pks,on = ['PEAK_NUM']).copy()
    #return(gdf_pass_pks)
    return(gdf_tot_pks)



def wt_time_Locs(wt,loc):
    return(wt*loc)
    
def sumthing(thing):
    return(sum(thing))


def weightedLoc(df,lat,lon,by,val2avg):
    import pandas as pd
    df_use = df.loc[:,[(lat),(lon),(by),val2avg]]
    df_use.loc[:,'lat_wt'] = df_use.apply(lambda y: y[lat] * y[val2avg],axis = 1).copy()
    df_use.loc[:,'lon_wt'] = df_use.apply(lambda y: y[lon] * y[val2avg],axis = 1).copy()


    #sumwts =pd.DataFrame( df_use.groupby('min_read').apply(lambda y: sumthing(y['pk_maxCH4_AB'])),columns = {'totwts'})
    sumwts = pd.DataFrame(df_use.copy().groupby(str(by)).apply(lambda y: sumthing(y[str(val2avg)])),columns = {'totwts'})
    sumwts.loc[:,'min_reads'] = sumwts.copy().index
    sumwts = sumwts.reset_index(drop=True).rename(columns={"min_reads": str(by)})
    
    #sumwts = sumwts.rename(columns={"min_reads": str(by)})

    totlats = pd.DataFrame(df_use.groupby(str(by)).apply(lambda y: sumthing(y['lat_wt'])),columns = ['totlats'])

    totlats['min_reads'] = totlats.index.copy()
    totlats = totlats.reset_index(drop=True)
    totlats = totlats.rename(columns={"min_reads": str(by)})

    totlons = pd.DataFrame(df_use.groupby(str(by)).apply(lambda y: sumthing(y['lon_wt'])),columns = ['totlons'])

    totlons['min_reads'] = totlons.index.copy()
    totlons = totlons.reset_index(drop=True)
    totlons = totlons.rename(columns={"min_reads": str(by)})

    df_use = pd.merge(totlats,df_use,on = str(by))
    df_use = pd.merge(totlons,df_use,on = str(by))
    df_use = pd.merge(sumwts,df_use,on = str(by))
    

    df_use.loc[:,'overall_LON'] = df_use.apply(lambda y: y['totlons']/y['totwts'],axis = 1)
    df_use.loc[:,'overall_LAT'] = df_use.apply(lambda y: y['totlats']/y['totwts'],axis = 1)

    toreturn = df_use.loc[:,[(str(by)),('overall_LON'),('overall_LAT')]].drop_duplicates()
    toreturn = toreturn.rename(columns = {'overall_LON':str(lon),'overall_LAT':str(lat)})

    return(toreturn)   


    
def verPk(totalData):
    import pandas as pd #
    from numpy import log
    import geopandas as gpd
    from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry

    totalData = totalData[totalData.numtimes != 1]
    pkRed = totalData[['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']].drop_duplicates().reset_index()
    verLoc = weightedLoc(pkRed,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB')
    pkRed.loc[:,('logCH4')]=pkRed.apply(lambda y: log(y.pk_maxCH4_AB),axis=1)
    mnVals = pkRed.groupby('min_read',as_index = False).logCH4.mean()
    together = pd.merge(verLoc,mnVals,on = ['min_read'])
    geometry_temp = [Point(xy) for xy in zip(together['pk_LON'], together['pk_LAT'])]
    crs = {'init': 'epsg:32610'}
    tog_dat = gpd.GeoDataFrame(together,crs = crs,geometry=geometry_temp)
    tog_dat = tog_dat.to_crs(epsg = 3857)


# =============================================================================
# summarizeDat
# Input: Dataframe with columns:
# Output: Dataframe with the average of the log(max_ch4_ab) for each verified peak (or observed if not verified yet)
# 
# =============================================================================

def summarizeDat(totalData):
    pkRed = totalData.loc[:,['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']].drop_duplicates().reset_index().loc[:,['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']]
    verLoc = weightedLoc(pkRed,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').rename(columns = {'pk_LAT':'overallLAT','pk_LON':'overallLON'})
    pkRed['logCH4'] = pkRed.apply(lambda y: log(y.pk_maxCH4_AB),axis = 1)
    mnVals = pkRed.groupby('min_read',as_index=False).logCH4.mean().rename(columns ={'logCH4':'mnlogCH4'}).loc[:,['min_read','mnlogCH4']]
    together = pd.merge(verLoc,mnVals,on = ['min_read'])
    final = pd.merge(together,totalData,on=['min_read'])
    return(final)
    


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
        
    
        