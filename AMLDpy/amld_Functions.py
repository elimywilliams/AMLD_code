#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 10:51:28 2020

@author: emilywilliams
"""


######### ALL FUNCTIONS NECESSARY TO RUN AMLD CODE ####################################

########################################################################
#### unique
# Helper function to find unique items in a list (used later)
# Input: a list
# Output: a list of unique items from original list

def unique(my_list):
    return [x for x in my_list if x not in locals()['_[1]']]


########################################################################
#### unIfInt
# helper function to find intersection of two lists
# Input: two lists
# Output: one list with union of items

def unIfInt(a, b):
    if len(intersect(a, b)) != 0:
        return (list(set(a).union(b)))


########################################################################
#### intersect
# helper function to find the intersection of 2 lists
# Input: two lists (a,b)
# Output: list of intersection of a and b

def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))


########################################################################
#### weightedLoc
# helper function to find weighted Location of a peak
# Input:
#   df = dataframe with data
#   lat = name of column with latitude data
#   lon = name of column with longitude data
#   by = name of column to group by
#   val2avg = name of column to use to weight
# Output: dataframe with weighted location for each grouping variable


def weightedLoc(df, lat, lon, by, val2avg):
    import pandas as pd
    import swifter
    df_use = df.loc[:, [(lat), (lon), (by), val2avg]]
    df_use.loc[:, 'lat_wt'] = df_use.swifter.apply(lambda y: y[lat] * y[val2avg], axis=1).copy()
    df_use.loc[:, 'lon_wt'] = df_use.swifter.apply(lambda y: y[lon] * y[val2avg], axis=1).copy()

    sumwts = pd.DataFrame(df_use.copy().groupby(str(by)).apply(lambda y: sumthing(y[str(val2avg)])), columns={'totwts'})
    sumwts.loc[:, 'min_reads'] = sumwts.copy().index
    sumwts = sumwts.reset_index(drop=True).rename(columns={"min_reads": str(by)})
    totlats = pd.DataFrame(df_use.groupby(str(by)).apply(lambda y: sumthing(y['lat_wt'])), columns=['totlats'])
    totlats['min_reads'] = totlats.index.copy()
    totlats = totlats.reset_index(drop=True)
    totlats = totlats.rename(columns={"min_reads": str(by)})
    totlons = pd.DataFrame(df_use.groupby(str(by)).apply(lambda y: sumthing(y['lon_wt'])), columns=['totlons'])
    totlons['min_reads'] = totlons.index.copy()
    totlons = totlons.reset_index(drop=True)
    totlons = totlons.rename(columns={"min_reads": str(by)})
    df_use = pd.merge(totlats, df_use, on=str(by))
    df_use = pd.merge(totlons, df_use, on=str(by))
    df_use = pd.merge(sumwts, df_use, on=str(by))
    df_use.loc[:, 'overall_LON'] = df_use.swifter.apply(lambda y: y['totlons'] / y['totwts'], axis=1)
    df_use.loc[:, 'overall_LAT'] = df_use.swifter.apply(lambda y: y['totlats'] / y['totwts'], axis=1)
    toreturn = df_use.loc[:, [(str(by)), ('overall_LON'), ('overall_LAT')]].drop_duplicates()
    toreturn = toreturn.rename(columns={'overall_LON': str(lon), 'overall_LAT': str(lat)})
    return (toreturn)


########################################################################
#### verPK
# helper function to make dataframe of verified peaks
# Input: df with all data
# Output: geodataframe of locations

def verPk(totalData):
    import pandas as pd  #
    from numpy import log
    import geopandas as gpd
    import swifter
    from shapely.geometry import Point  # Shapely for converting latitude/longtitude to geometry

    totalData = totalData[totalData.numtimes != 1]
    pkRed = totalData[
        ['PEAK_NUM', 'pk_LON', 'pk_LAT', 'pk_maxCH4_AB', 'numtimes', 'min_read']].drop_duplicates().reset_index()
    verLoc = weightedLoc(pkRed, 'pk_LAT', 'pk_LON', 'min_read', 'pk_maxCH4_AB')
    pkRed.loc[:, ('logCH4')] = pkRed.swifter.apply(lambda y: log(y.pk_maxCH4_AB), axis=1)
    mnVals = pkRed.groupby('min_read', as_index=False).logCH4.mean()
    together = pd.merge(verLoc, mnVals, on=['min_read'])
    geometry_temp = [Point(xy) for xy in zip(together['pk_LON'], together['pk_LAT'])]
    crs = {'init': 'epsg:32610'}
    tog_dat = gpd.GeoDataFrame(together, crs=crs, geometry=geometry_temp)
    tog_dat = tog_dat.to_crs(epsg=3857)


########################################################################
#### estEmissions
# Input: excess CH4 value
# Output: estimated emission

def estEmissions(excessCH4):
    import math
    a = 0.4630664
    b = 0.7443749
    a1 = 1.2889
    b1 = 0.35232
    a2 = 1.755891
    b2 = 0.4438203

    m = math.exp((excessCH4 - a) / b)
    # if m < math.exp(3.157):
    #    if m < math.exp(2):
    #       m = math.exp((np.log(m) - a1)/b1)
    #  if m > math.exp(2):
    #     m = math.exp((np.log(m) - a2)/b2)
    return (m)


########################################################################
#### haversine
# Input: two locations (lat,long), radius of earth
# Output: distance (m) between points

def haversine(lat1, lon1, lat2, lon2, radius=6371):  # 6372.8 = earth radius in kilometers
    from math import radians, sin, cos, sqrt, asin
    dLat = radians(lat2 - lat1)
    dLon = radians(lon2 - lon1)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    c = 2 * asin(sqrt(sin(dLat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dLon / 2) ** 2))
    return radius * c * 1000  # return in meters


########################################################################
#### wt_time_Locs
# Input: weight, location
# Output: wt*location

def wt_time_Locs(wt, loc):
    return (wt * loc)


########################################################################
#### sumthing
# Input: values
# Output: sum of values

def sumthing(thing):
    return (sum(thing))


########################################################################
#### makeGEO
# Input: dataframe, lat,long columns
# Output: converts lat long to points (useful for geodataframes)


def makeGEO(df, lat, lon):
    from shapely.geometry import Point
    geo = [Point(xy) for xy in zip(df[(lon)], df[(lat)])]
    return (geo)


########################################################################
#### makeGPD
# Input: dataframe, lat,long locations, crs to use
# Output: returns a geodataframe with the lat longs as points

def makeGPD(df, lat, lon, cps='epsg:4326'):
    import geopandas as gpd
    gdf = gpd.GeoDataFrame(df, crs=cps, geometry=makeGEO(df, lat, lon))
    return (gdf)


########################################################################
#### summarizeDat
# Input: data from all analyses
# Output: summary information with the weighted locations

def summarizeDat(totalData):
    import pandas as pd
    from numpy import log
    pkRed = totalData.loc[:, ['PEAK_NUM', 'pk_LON', 'pk_LAT', 'pk_maxCH4_AB', 'numtimes',
                              'min_read']].drop_duplicates().reset_index().loc[:,
            ['PEAK_NUM', 'pk_LON', 'pk_LAT', 'pk_maxCH4_AB', 'numtimes', 'min_read']]
    verLoc = weightedLoc(pkRed, 'pk_LAT', 'pk_LON', 'min_read', 'pk_maxCH4_AB').rename(
        columns={'pk_LAT': 'overallLAT', 'pk_LON': 'overallLON'})
    pkRed['logCH4'] = pkRed.apply(lambda y: log(y.pk_maxCH4_AB), axis=1)
    mnVals = pkRed.groupby('min_read', as_index=False).logCH4.mean().rename(columns={'logCH4': 'mnlogCH4'}).loc[:,
             ['min_read', 'mnlogCH4']]
    together = pd.merge(verLoc, mnVals, on=['min_read'])
    final = pd.merge(together, totalData, on=['min_read'])
    return (final)


########################################################################
#### getQuad
# Input: x,y values
# Output: which quadrant of a xy plot it is in (useful for wind info)

def getQuad(x, y):
    try:
        x = int(x)
        y = int(y)
    except ValueError:
        return (0)

    if y >= 0 and x > 0:
        return (1)
    elif y >= 0 and x < 0:
        return (2)
    elif y < 0 and x < 0:
        return (3)
    else:
        return (4)


########################################################################
#### calcTheta
# Input: u,v (wind), quadrant, length of horizontal wind, radians (t/f)
# Output: calculates theta value

def calcTheta(U, V, quad, h_length, radians):
    import numpy as np
    theta = np.arcsin(U / h_length)
    import numpy as np
    if quad == 1:
        theta = theta
    elif quad == 2:
        theta = -theta + np.pi / 2
    elif quad - - 3:
        theta = np.pi / 2 + theta + np.pi
    elif quad == 4:
        theta = 3 * np.pi / 2
    theta = 2 * np.pi - theta
    if not radians:
        theta = theta * 180 / np.pi

    return (theta)


########################################################################
#### calcBearing
# Input: two gps coordinates, radians(t/f)
# Output: direction (degrees/radians) of motion

def calcBearing(lat1, lat2, long1, long2, radians):
    from math import atan2
    from numpy import pi
    from math import radians, sin, cos

    lat1r = lat1 * (pi / 180)
    lat2r = lat2 * (pi / 180)
    long1r = long1 * (pi / 180)
    long2r = long2 * (pi / 180)
    X = cos(lat2r) * sin(long2r - long1r)
    Y = cos(lat1r) * sin(lat2r) - (sin(lat1r) * cos(lat2r) * cos(long2r - long1r))

    theta = atan2(X, Y)
    theta = theta % (2 * pi)

    if not radians:
        return (theta * 180 / pi)
    elif radians:
        return (theta)


########################################################################
#### ProcessRawDataEng
# Input: a .txt file with raw data
# Output: saves a log, and .csv file with processed data (removes unwanted)

def ProcessRawDataEng(xCar, xDate, xDir, xFilename, bFirst, gZIP, xOut, initialTimeBack,
                      shift, maxSpeed='45', minSpeed='2'):
    import pandas as pd
    from datetime import datetime
    import os
    import gzip
    # import csv
    from math import floor
    try:
        xMaxCarSpeed = float(maxSpeed) / 2.23694  # CONVERTED TO M/S (default is 45mph)
        xMinCarSpeed = float(minSpeed) / 2.23694  # CONVERTED TO M/S (default is 2mph)

        ########################################################################
        #### WE DON'T HAVE AN RSSI INPUT
        ### (SO THIS IS A PLACEHOLDER FOR SOME SORT OF QA/QC VARIABLE)
        ##  xMinRSSI = 50  #if RSSI is below this we don't like it
        ##################################################################

        # reading in the (.txt) data with specific headers --> need to change this
        #          0     1    2    3       4           5    6       7        8        9          10                 11              12           13            14      15      16      17        18         19         20         21         22         23        24   25  26       27           28       29           30       31       32       33  34        35   36   37  38   39       40       41   42       43   44   45   46   47   48   49   50   51     52     53     54
        # sHeader = "Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude"
        # sHeader = "Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude"
        # sHeader = "Time Stamp,Inlet Number,P (mbars),T0 (degC),T5 (degC), Laser PID Readout,Det PID Readout,win0Fit0,win0Fit1,win0Fit3,win1Fit4,win0Fit5,win0Fit6,win0Fit7,win0Fit8,win0Fit9,win1Fit0,win1Fit1,win1Fit2,win1Fit3,win1Fit4,win1Fit5,win1Fit6,Det Bkgd,Ramp Ampl,CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Battery T (degC),FET T (degC),GPS Time,Latitude,Longitude"
        sHeader = "Time Stamp,Inlet Number,P (mbars),T0 (degC),T5 (degC),Laser PID Readout,Det PID Readout,win0Fit0,win0Fit1,win0Fit2,win0Fit3,win0Fit4,win0Fit5,win0Fit6,win0Fit7,win0Fit8,win0Fit9,win1Fit0,win1Fit1,win1Fit2,win1Fit3,win1Fit4,win1Fit5,win1Fit6,Det Bkgd,Ramp Ampl,CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Battery T (degC),FET T (degC),GPS Time,Latitude,Longitude"
        sOutHeader = "DATE,TIME,SECONDS,NANOSECONDS,VELOCITY,U,V,W,BCH4,BRSSI,TCH4,TRSSI,PRESS_MBAR,INLET,TEMPC,CH4,H20,C2H6,R,C2C1,BATTV,POWMV,CURRMA,SOCPER,LAT,LONG\n"

        headerNames = sHeader.split(',')
        GPS_loc = 37  # Where the GPS data is located (in the row)

        infoHeader = "FILENAME\n"

        # gZIP is indicating if it is a ZIP file (I don't think I've written this in)
        if gZIP == 0:
            f = gzip.open(xDir + "/" + xFilename,
                          'r')  # if in python 3, change this to "r" or just "b" can't remember but something about a bit not a string
        else:
            f = open(xDir + "/" + xFilename, 'r')
            f = open(xDir + xFilename, 'r')

        ### FIGURING OUT DATE FROM FILENAME (WILL NEED TO CHANGE THIS IF DIFFERENT FILENAME)
        xdat = str('20') + xFilename[11:17]

        # fnOut = xOutDir + xCar + "_" + xDate.replace("-", "") + "_dat.csv"       #set CSV output for raw data
        # fnLog = xOutDir + xCar + "_" + xDate.replace("-", "") + "_log.csv"       #output for logfile

        fnOut = xOut + xCar + "_" + xdat + "_dat.csv"  # set CSV output for raw data
        fnLog = xOut + xCar + "_" + xdat + "_log.csv"  # output for logfile
        infOut = xOut + xCar + "_" + xdat + "_info.csv"

        # FINDING THE FIRST TIME NOTED
        firsttime = int(float(open(xDir + xFilename).readlines().pop(1).split(',')[37][:-4]))

        ## MAKING TEMPORARY FILE (FOR IF LATER YOU HAVE TO ADD A DATE)
        fnOutTemp = xOut + xCar + "_" + xdat + "temp_dat.csv"  #

        if bFirst:
            # 3fOut = open(fnOutTemp, 'w')
            # fOut.write(sOutHeader)
            fLog = open(fnLog, 'w')
            infOut = open(infOut, 'w')
            infOut.write(infoHeader)
            print(f"fnLog:{fnOut}")
        if not bFirst:
            fOut = open(fnOut, 'a')
            fLog = open(fnLog, 'a')
            infOut = open(infOut, 'a')

        fOut = open(fnOutTemp, 'w')
        fOut.write(sOutHeader)

        # READ IN THE LINES
        xCntObs = -1
        xCntGoodValues = 0
        for row in f:
            bGood = True
            if xCntObs < 0:
                bGood = False
                xCntObs += 1

            if bGood:
                lstS = row.split(',')
                gpstime = lstS[GPS_loc]
                dtime = lstS[0]
                dt = lstS[1]
                time_dt = lstS[2]
                epoch = lstS[3]
                # nano = lstS[4]

                gps_time = lstS[37]
                dateob = datetime.fromtimestamp(int(gps_time[:-4]))
                nano = gps_time[-4:]

                # dateob = datetime(int(dt[0:4]),int(dt[5:7]),int(dt[8:10]),int(time_dt[0:2]),int(time_dt[3:5]),int(time_dt[6:8]),int(float(nano)*1e-9))

                dtime = int(dateob.strftime('%Y%m%d%H%M%S'))
                # Date = dateob.strftime('%m%/%d/%Y')
                Date = dateob.strftime('%Y-%m-%d')

                GPS_Time = dateob.strftime('%H%:%M:%S')
                seconds = floor(float(gpstime))
                nano = dateob.strftime('%f')

                # dateob = datetime(int(dtime[6:10]),int(dtime[0:2]),int(dtime[3:5]),int(dtime[11:13]),int(dtime[14:16]),int(dtime[17:19]),int(float(dtime[19:23])*1000000))
                # epoch = dateob.strftime('%s.%f')

                # THIS IS USING THE CSU METHOD. IN OUR METHOD, WE DO THE SPEED LATER IN THE ALGORITHM.

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
                import pandas as pd
                # if sys.platform.startswith('win'):
                #     ## DATE, TIME, SECONDS,NANOSECONDS
                #     csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(dateob.strftime('%H:%M:%S')) + ',' + str(
                #         float(pd.to_numeric(dateob.strftime('%S.%f')))) + ',' + str(
                #         pd.to_numeric(dateob.strftime('%f')) * 1000) + str(',')
                #     ## VELOCITY, U,V,W,BCH4,BRSSI,TCH4
                #     csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(
                #         lstS[26]) + ',' + str('0') + ',' + str(lstS[26]) + ','
                #     ## TRSSI, PRESS_MBAR, INLET, TEMPC, CH4, H20,C2H6
                #     csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(
                #         lstS[26]) + ',' + str(lstS[27]) + ',' + str(lstS[28]) + ','
                #     # R, C2C1, BATTV, POWMV,CURRMA, SOCPER,LAT,LONG
                #     csvWrite += str(lstS[29]) + ',' + str(lstS[30]) + ',' + str(lstS[31]) + ',' + str(
                #         lstS[32]) + ',' + str(lstS[33]) + ',' + str(lstS[34]) + ',' + str(lstS[38]) + str(',') + str(
                #         lstS[39])

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
                #if not sys.platform.startswith('win'):
                if 1 == 1:
                    ## DATE, TIME, SECONDS,NANOSECONDS
                    csvWrite = str(Date) + ',' + str(GPS_Time) + ',' + str(seconds) + ',' + str(nano) + str(',')
                    ## VELOCITY, U,V,W,BCH4,BRSSI,TCH4
                    csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(
                        lstS[26]) + ',' + str('0') + ',' + str(lstS[26]) + ','
                    ## TRSSI, PRESS_MBAR, INLET, TEMPC, CH4, H20,C2H6
                    csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(
                        lstS[26]) + ',' + str(lstS[27]) + ',' + str(lstS[28]) + ','
                    # R, C2C1, BATTV, POWMV,CURRMA, SOCPER,LAT,LONG
                    csvWrite += str(lstS[29]) + ',' + str(lstS[30]) + ',' + str(lstS[31]) + ',' + str(
                        lstS[32]) + ',' + str(lstS[33]) + ',' + str(lstS[34]) + ',' + str(lstS[38]) + str(',') + str(
                        lstS[39])
                # fOut.write('\n')

                #### REMOVING THE FIRST BIT OF DATA (if you need to )
                if seconds >= (firsttime + (60 * float(initialTimeBack))):
                    fOut.write(csvWrite)

                del (csvWrite)
            #                xCntGoodValues += 1

            xCntObs += 1

        # sOut = str(gZIP) + "," + str(f) + "," + str(xCntObs) + "," + str(xCntGoodValues) + "\n"
        # fLog.write(sOut)

        infOut.write(str(xFilename) + '\n')

        fOut.close()
        fLog.close()
        infOut.close()

        # xDate = dateob.strftime("%Y%m%d")

        # newfnOut = xOutDir + xCar + "_" + xDate + "_dat.csv"       #set CSV output for raw data
        # newfnLog = xOutDir + xCar + "_" + xDate + "_log.csv"

        #print(xCar + "\t" + xdat + "\t" + fnOut[-22:] + "\t" + str(xCntObs) + "\t" + str(xCntGoodValues) + "\t" + str(
        #    gZIP))

        print(f"{xCar} \t {xdat} \t {fnOut[-(17 + len(xCar)):]} \t {xCntObs} \t {xCntGoodValues} \t {gZIP}")

        def calcVel(timediff, distance):
            if timediff == 0:
                return (0)
            elif timediff != 0:
                return (distance / timediff)

        import numpy as np
        radians = False
        wind_df = pd.read_csv(fnOutTemp)
        wind_df['QUADRANT'] = wind_df.apply(lambda row: getQuad(row['U'], row['V']), axis=1)

        wind_df['secnan'] = wind_df.apply(lambda row: row['SECONDS'], axis=1)  # + row['NANOSECONDS']*1e-9,axis=1)
        wind_df['prev_LAT'] = wind_df.LAT.shift(periods=1)
        wind_df['next_LAT'] = wind_df.LAT.shift(periods=-1)
        wind_df['prev_LONG'] = wind_df.LONG.shift(periods=1)
        wind_df['next_LONG'] = wind_df.LONG.shift(periods=-1)
        wind_df['prev_TIME'] = wind_df.secnan.shift(periods=1)
        wind_df['next_TIME'] = wind_df.secnan.shift(periods=-1)
        wind_df['distance'] = wind_df.apply(
            lambda row: haversine(row['prev_LAT'], row['prev_LONG'], row['next_LAT'], row['next_LONG']), axis=1)
        wind_df['bearing'] = wind_df.apply(
            lambda row: calcBearing(row['prev_LAT'], row['next_LAT'], row['prev_LONG'], row['next_LONG'], radians),
            axis=1)
        wind_df['timediff'] = wind_df.apply(lambda row: row['next_TIME'] - row['prev_TIME'], axis=1)
        wind_df['VELOCITY'] = wind_df.apply(lambda row: calcVel(row['timediff'], row['distance']), axis=1)
        wind_df['U_cor'] = wind_df.apply(lambda row: row['U'] + row['VELOCITY'], axis=1)
        wind_df['horz_length'] = wind_df.apply(lambda row: np.sqrt(row['U_cor'] ** 2 + row['V'] ** 2), axis=1)
        wind_df['uncor_theta'] = wind_df.apply(
            lambda row: calcBearing(row['U_cor'], row['V'], row['QUADRANT'], row['horz_length'], radians), axis=1)
        wind_df['adj_theta'] = wind_df.apply(lambda row: (row['uncor_theta'] + row['bearing']) % 360, axis=1)
        wind_df['totalWind'] = wind_df.apply(lambda row: np.sqrt(row['horz_length'] ** 2 + row['W'] ** 2), axis=1)
        wind_df['phi'] = wind_df.apply(lambda row: np.arctan(row['horz_length']), axis=1)
        wind_df['shift_CH4'] = wind_df.CH4.shift(periods=int(float(shift)))
        wind_df['raw_CH4'] = wind_df.apply(lambda row: row['BCH4'], axis=1)
        wind_df['BCH4'] = wind_df.loc[:, ['shift_CH4']]
        wind_df['CH4'] = wind_df.loc[:, ['shift_CH4']]
        wind_df['TCH4'] = wind_df.loc[:, ['shift_CH4']]

        wind_df2 = wind_df[wind_df.CH4.notnull()]
        wind_df3 = wind_df2.drop(
            ['QUADRANT', 'secnan', 'prev_LAT', 'next_LAT', 'prev_LONG', 'next_LONG', 'prev_TIME', 'next_TIME',
             'distance', 'timediff', 'uncor_theta', 'CH4'], axis=1)
        wind_df3['CH4'] = wind_df3.loc[:, 'shift_CH4']
        wind_df3 = wind_df3.drop(['shift_CH4'], axis=1)

        wind_df3 = wind_df3.loc[:,
                   ['DATE', 'TIME', 'SECONDS', 'NANOSECONDS', 'VELOCITY', 'U', 'V', 'W', 'BCH4', 'BRSSI', 'TCH4',
                    'TRSSI', 'PRESS_MBAR', 'INLET' \
                       , 'TEMPC', 'CH4', 'H20', 'C2H6', 'R', 'C2C1', 'BATTV', 'POWMV', 'CURRMA', 'SOCPER', 'LAT',
                    'LONG', 'bearing', 'U_cor', \
                    'horz_length', 'adj_theta', 'totalWind', 'phi', 'raw_CH4']]
        wind_df4 = wind_df3.loc[wind_df3.totalWind.notnull(), :]

        wind_df7 = addOdometer(wind_df4, 'LAT', 'LONG')
        wind_df4 = wind_df7.copy()
        wind_df5 = wind_df4.loc[wind_df4.VELOCITY > xMinCarSpeed, :]
        wind_df6 = wind_df5.loc[wind_df5.VELOCITY < xMaxCarSpeed, :]

        del (wind_df4)

        # wind_df7 = addOdometer(wind_df6,'LAT','LONG')
        wind_df4 = wind_df6.copy().drop_duplicates()
        # del(wind_df7)

        # firstTime = wind_df3.SECONDS.min() + 60 *(initialTimeBack)
        # wind_df4 = wind_df3.loc[wind_df3.SECONDS > firstTime,:]
        # wind_df3.to_csv(fnOutTemp,index=False)

        if bFirst:
            wind_df4.to_csv(fnOut, index=False)
        elif not bFirst:
            norm = pd.read_csv(fnOut)
            pd.concat([norm, wind_df4]).sort_values(by='SECONDS').reset_index(drop=True).to_csv(fnOut, index=False)
        os.remove(fnOutTemp)
        return True
    except ValueError:
        return False


############
### processData (not engineering file)
def ProcessRawData(xCar, xDate, xDir, xFilename, bFirst, gZIP, xOut, initialTimeBack, shift, maxSpeed='45',
                   minSpeed='2'):
    import pandas as pd
    from datetime import datetime
    import os
    import gzip
    # import csv
    try:
        xMaxCarSpeed = float(maxSpeed) / 2.23694  # CONVERTED TO M/S (default is 45mph)
        xMinCarSpeed = float(minSpeed) / 2.23694  # CONVERTED TO M/S (default is 2mph)

        ########################################################################
        #### WE DON'T HAVE AN RSSI INPUT
        ### (SO THIS IS A PLACEHOLDER FOR SOME SORT OF QA/QC VARIABLE)
        ##  xMinRSSI = 50  #if RSSI is below this we don't like it
        ##################################################################

        # reading in the data with specific headers
        #          0     1    2    3       4           5    6       7        8        9          10                 11              12           13            14      15      16      17        18         19         20         21         22         23        24   25  26       27           28       29           30       31       32       33  34        35   36   37  38   39       40       41   42       43   44   45   46   47   48   49   50   51     52     53     54
        sHeader = "Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude"
        sOutHeader = "DATE,TIME,SECONDS,NANOSECONDS,VELOCITY,U,V,W,BCH4,BRSSI,TCH4,TRSSI,PRESS_MBAR,INLET,TEMPC,CH4,H20,C2H6,R,C2C1,BATTV,POWMV,CURRMA,SOCPER,LAT,LONG\n"
        infoHeader = "FILENAME\n"
        # somehow gZIP is indicating if  it is the first file name (I think if it is 0 then it is the first file)
        if gZIP == 0:
            f = gzip.open(xDir + "/" + xFilename,
                          'r')  # if in python 3, change this to "r" or just "b" can't remember but something about a bit not a string
        else:
            f = open(xDir + "/" + xFilename, 'r')

        infoHeader = "FILENAME\n"

        # process
        # if first time on this car/date, then write header out
        headerNames = sHeader.split(',')
        xdat = str('20') + xFilename[11:17]
        #xdat = str(xFilename[len(xCar)+1:len(xCar) + 9])

        # fnOut = xOutDir + xCar + "_" + xDate.replace("-", "") + "_dat.csv"       #set CSV output for raw data
        # fnLog = xOutDir + xCar + "_" + xDate.replace("-", "") + "_log.csv"       #output for logfile

        fnOut = xOut + xCar + "_" + xdat + "_dat.csv"  # set CSV output for raw data
        fnLog = xOut + xCar + "_" + xdat + "_log.csv"  # output for logfile
        infOut = xOut + xCar + "_" + xdat + "_info.csv"
        #

        dtime = open(xDir + xFilename).readlines().pop(1).split(',')[0]
        firstdate = datetime(int(dtime[6:10]), int(dtime[0:2]), int(dtime[3:5]), int(dtime[11:13]), int(dtime[14:16]),
                             int(dtime[17:19]), int(float(dtime[19:23]) * 1000000))
        firsttime = firstdate.strftime('%s.%f')

        # firsttime = int(float(open(xDir + xFilename).readlines().pop(1).split(',')[37][:-4]))

        fnOutTemp = xOut + xCar + "_" + xdat + "temp_dat.csv"  #

        if bFirst:
            # fOut = open(fnOut, 'w')
            # fOut.write(sOutHeader)
            fLog = open(fnLog, 'w')
            infOut = open(infOut, 'w')
            infOut.write(infoHeader)
            print(f"fnLog:{fnOut}")
        if not bFirst:
            fOut = open(fnOut, 'a')
            fLog = open(fnLog, 'a')
            infOut = open(infOut, 'a')

        fOut = open(fnOutTemp, 'w')
        fOut.write(sOutHeader)

        # read all lines
        xCntObs = -1
        xCntGoodValues = 0
        for row in f:
            # print(row)
            bGood = True
            if xCntObs < 0:
                bGood = False
                xCntObs += 1
            if bGood:
                lstS = row.split(",")
                dtime = lstS[0]
                dateob = datetime(int(dtime[6:10]), int(dtime[0:2]), int(dtime[3:5]), int(dtime[11:13]),
                                  int(dtime[14:16]), int(dtime[17:19]), int(float(dtime[19:23]) * 1000000))
                # epoch = dateob.strftime('%s.%f')
                # dtime = int(dateob.strftime('%Y%m%d%H%M%S'))

                fdate = datetime(int(dtime[6:10]), int(dtime[0:2]), int(dtime[3:5]), int(dtime[11:13]),
                                 int(dtime[14:16]), int(dtime[17:19]), int(float(dtime[19:23]) * 1000000))
                seconds = fdate.strftime('%s.%f')

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
                    csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(dateob.strftime('%H:%M:%S')) + ',' + str(
                        int(pd.to_numeric(dateob.strftime('%S.%f')))) + ',' + str(
                        pd.to_numeric(dateob.strftime('%f')) * 1000) + str(',')
                    csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(
                        lstS[4]) + ',' + str('0') + ',' + str(lstS[4]) + ','
                    csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(
                        lstS[4]) + ',' + str(lstS[5]) + ',' + str(lstS[6]) + ','
                    csvWrite += str(lstS[7]) + ',' + str(lstS[8]) + ',' + str(lstS[9]) + ',' + str(
                        lstS[10]) + ',' + str(lstS[11]) + ',' + str(lstS[12]) + ',' + str(lstS[13]) + str(',') + str(
                        lstS[14])
                if not sys.platform.startswith('win'):
                    csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(dateob.strftime('%H:%M:%S')) + ',' + str(
                        seconds[:10]) + ',' + str(pd.to_numeric(seconds[11:]) * 1000) + str(',')
                    csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(
                        lstS[4]) + ',' + str('0') + ',' + str(lstS[4]) + ','
                    csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(lstS[3]) + ',' + str(
                        lstS[4]) + ',' + str(lstS[5]) + ',' + str(lstS[6]) + ','
                    csvWrite += str(lstS[7]) + ',' + str(lstS[8]) + ',' + str(lstS[9]) + ',' + str(
                        lstS[10]) + ',' + str(lstS[11]) + ',' + str(lstS[12]) + ',' + str(lstS[13]) + str(',') + str(
                        lstS[14])
                if float(seconds) >= (float(firsttime) + (60 * float(initialTimeBack))):
                    fOut.write(csvWrite)
                    del (seconds)
                del (csvWrite)

            xCntObs += 1

        # sOut = str(gZIP) + "," + str(f) + "," + str(xCntObs) + "," + str(xCntGoodValues) + "\n"
        # fLog.write(sOut)
        infOut.write(str(xFilename) + '\n')

        fOut.close()
        fLog.close()
        infOut.close()

        # xDate = dateob.strftime("%Y%m%d")

        # newfnOut = xOutDir + xCar + "_" + xDate + "_dat.csv"       #set CSV output for raw data
        # newfnLog = xOutDir + xCar + "_" + xDate + "_log.csv"

        #print(xCar + "\t" + xdat + "\t" + fnOut[-22:] + "\t" + str(xCntObs) + "\t" + str(xCntGoodValues) + "\t" + str(
        #    gZIP))

        print(f"{xCar} \t {xdat} \t {fnOut[-(17 + len(xCar)):]} \t {xCntObs} \t {xCntGoodValues} \t  {gZIP}")
        from numpy import pi
        import numpy as np
        def calcVel(timediff, distance):
            if timediff == 0:
                return (0)
            elif timediff != 0:
                return (distance / timediff)

        wind_df = pd.read_csv(fnOutTemp)
        radians = False
        wind_df['QUADRANT'] = wind_df.apply(lambda row: getQuad(row['U'], row['V']), axis=1)
        wind_df['secnan'] = wind_df.apply(lambda row: row['SECONDS'] + row['NANOSECONDS'] * 1e-9,
                                          axis=1)  # + row['NANOSECONDS']*1e-9,axis=1)
        wind_df['prev_LAT'] = wind_df.LAT.shift(periods=1)
        wind_df['next_LAT'] = wind_df.LAT.shift(periods=-1)
        wind_df['prev_LONG'] = wind_df.LONG.shift(periods=1)
        wind_df['next_LONG'] = wind_df.LONG.shift(periods=-1)
        wind_df['prev_TIME'] = wind_df.secnan.shift(periods=1)
        wind_df['next_TIME'] = wind_df.secnan.shift(periods=-1)
        wind_df['distance'] = wind_df.apply(
            lambda row: haversine(row['prev_LAT'], row['prev_LONG'], row['next_LAT'], row['next_LONG']), axis=1)
        wind_df['bearing'] = wind_df.apply(
            lambda row: calcBearing(row['prev_LAT'], row['next_LAT'], row['prev_LONG'], row['next_LONG'], radians),
            axis=1)
        wind_df['timediff'] = wind_df.apply(lambda row: row['next_TIME'] - row['prev_TIME'], axis=1)
        wind_df['VELOCITY'] = wind_df.apply(lambda row: calcVel(row['timediff'], row['distance']), axis=1)
        wind_df['U_cor'] = wind_df.apply(lambda row: row['U'] + row['VELOCITY'], axis=1)
        wind_df['horz_length'] = wind_df.apply(lambda row: np.sqrt(row['U_cor'] ** 2 + row['V'] ** 2), axis=1)
        wind_df['uncor_theta'] = wind_df.apply(
            lambda row: calcBearing(row['U_cor'], row['V'], row['QUADRANT'], row['horz_length'], radians), axis=1)
        wind_df['adj_theta'] = wind_df.apply(lambda row: (row['uncor_theta'] + row['bearing']) % 360, axis=1)
        wind_df['totalWind'] = wind_df.apply(lambda row: np.sqrt(row['horz_length'] ** 2 + row['W'] ** 2), axis=1)
        wind_df['phi'] = wind_df.apply(lambda row: np.arctan(row['horz_length']), axis=1)
        wind_df['shift_CH4'] = wind_df.CH4.shift(periods=int(float(shift)))
        wind_df['raw_CH4'] = wind_df.apply(lambda row: row['BCH4'], axis=1)
        wind_df['BCH4'] = wind_df.loc[:, ['shift_CH4']]
        wind_df['CH4'] = wind_df.loc[:, ['shift_CH4']]
        wind_df['TCH4'] = wind_df.loc[:, ['shift_CH4']]

        wind_df2 = wind_df[wind_df.CH4.notnull()]
        wind_df2 = wind_df.copy()
        wind_df3 = wind_df2.drop(
            ['QUADRANT', 'secnan', 'prev_LAT', 'next_LAT', 'prev_LONG', 'next_LONG', 'prev_TIME', 'next_TIME',
             'timediff', 'uncor_theta', 'CH4'], axis=1)
        wind_df3['CH4'] = wind_df3.loc[:, 'shift_CH4']
        wind_df3 = wind_df3.drop(['shift_CH4'], axis=1)

        wind_df3 = wind_df3.loc[:,
                   ['DATE', 'TIME', 'SECONDS', 'NANOSECONDS', 'VELOCITY', 'U', 'V', 'W', 'BCH4', 'BRSSI', 'TCH4',
                    'TRSSI', 'PRESS_MBAR', 'INLET' \
                       , 'TEMPC', 'CH4', 'H20', 'C2H6', 'R', 'C2C1', 'BATTV', 'POWMV', 'CURRMA', 'SOCPER', 'LAT',
                    'LONG', 'bearing', 'U_cor', \
                    'horz_length', 'adj_theta', 'totalWind', 'phi', 'raw_CH4', 'distance']]
        wind_df3['odometer'] = wind_df3.loc[:, 'distance'].cumsum()
        # wind_df4 = wind_df3.loc[wind_df3.totalWind.notnull(),:]

        wind_df4 = wind_df3.copy()

        # wind_df7 = addOdometer(wind_df4,'LAT','LONG')

        # wind_df4 = wind_df7.copy()
        wind_df5 = wind_df4.loc[wind_df4.VELOCITY > xMinCarSpeed, :]
        wind_df6 = wind_df5.loc[wind_df5.VELOCITY < xMaxCarSpeed, :]

        del (wind_df4)
        wind_df4 = wind_df6.copy().drop_duplicates()
        wind_df5 = wind_df4.loc[wind_df4.CH4.notnull(), :]
        wind_df4 = wind_df5.copy()
        if bFirst:
            wind_df4.to_csv(fnOut, index=False)
        elif not bFirst:
            norm = pd.read_csv(fnOut)
            pd.concat([norm, wind_df4]).sort_values(by='SECONDS').reset_index(drop=True).to_csv(fnOut, index=False)
        os.remove(fnOutTemp)
        return True
    except ValueError:
        return False


def ProcessRawDataAeris(xCar, xDate, xDir, xFilename, bFirst, gZIP, xOut, initialTimeBack, shift, maxSpeed='45',
                        minSpeed='2'):
    import pandas as pd
    from datetime import datetime
    import os
    import gzip
    import numpy as np
    # import csv
    try:
        xMaxCarSpeed = float(maxSpeed) / 2.23694  # CONVERTED TO M/S (default is 45mph)
        xMinCarSpeed = float(minSpeed) / 2.23694  # CONVERTED TO M/S (default is 2mph)

        ########################################################################
        #### WE DON'T HAVE AN RSSI INPUT
        ### (SO THIS IS A PLACEHOLDER FOR SOME SORT OF QA/QC VARIABLE)
        ##  xMinRSSI = 50  #if RSSI is below this we don't like it
        ##################################################################

        # reading in the data with specific headers
        #          0     1    2    3       4           5    6       7        8        9          10                 11              12           13            14      15      16      17        18         19         20         21         22         23        24   25  26       27           28       29           30       31       32       33  34        35   36   37  38   39       40       41   42       43   44   45   46   47   48   49   50   51     52     53     54
        sHeader = "Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude"
        sHeader = 'Time Stamp,Inlet Number,P (mbars),T (degC),CH4 (ppm),H2O (ppm),C2H6 (ppb),R,C2/C1,Battery Charge (V),Power Input (mV),Current (mA),SOC (%),Latitude,Longitude,U (m/sec),V (m/sec),W (m/sec),T (degC),Dir (deg),Speed (m/sec),Compass (deg)'
        sOutHeader = "DATE,TIME,SECONDS,NANOSECONDS,VELOCITY,U,V,W,BCH4,BRSSI,TCH4,TRSSI,PRESS_MBAR,INLET,TEMPC,CH4,H20,C2H6,R,C2C1,BATTV,POWMV,CURRMA,SOCPER,LAT,LONG\n"
        infoHeader = "FILENAME\n"
        # somehow gZIP is indicating if  it is the first file name (I think if it is 0 then it is the first file)
        if gZIP == 0:
            f = gzip.open(xDir + "/" + xFilename,
                          'r')  # if in python 3, change this to "r" or just "b" can't remember but something about a bit not a string
        else:
            f = open(xDir + xFilename, 'r')

        infoHeader = "FILENAME\n"

        # process - if first time on this car/date, then write header out
        headerNames = sHeader.split(',')
        xdat = str('20') + xFilename[11:17]
        fnOut = xOut + xCar + "_" + xdat + "_dat.csv"  # set CSV output for raw data
        fnLog = xOut + xCar + "_" + xdat + "_log.csv"  # output for logfile
        infOut = xOut + xCar + "_" + xdat + "_info.csv"
        #

        dtime = open(xDir + xFilename).readlines().pop(1).split(',')[0]
        firstdate = datetime(int(dtime[6:10]), int(dtime[0:2]), int(dtime[3:5]), int(dtime[11:13]),
                             int(dtime[14:16]), int(dtime[17:19]), int(float(dtime[19:23]) * 1000000))
        firsttime = firstdate.strftime('%s.%f')
        fnOutTemp = xOut + xCar + "_" + xdat + "temp_dat.csv"  #

        if bFirst:
            fLog = open(fnLog, 'w')
            infOut = open(infOut, 'w')
            infOut.write(infoHeader)
            print(f"fnLog:{fnOut}")
        if not bFirst:
            fOut = open(fnOut, 'a')
            fLog = open(fnLog, 'a')
            infOut = open(infOut, 'a')

        fOut = open(fnOutTemp, 'w')
        fOut.write(sOutHeader)

        # read all lines
        xCntObs = -1
        xCntGoodValues = 0
        for row in f:
            # print(row)
            bGood = True
            if xCntObs < 0:
                bGood = False
                xCntObs += 1
            if bGood:
                lstS = row.split(",")
                dtime = lstS[0]
                dateob = datetime(int(dtime[6:10]), int(dtime[0:2]), int(dtime[3:5]), int(dtime[11:13]),
                                  int(dtime[14:16]), int(dtime[17:19]), int(float(dtime[19:23]) * 1000000))
                fdate = datetime(int(dtime[6:10]), int(dtime[0:2]), int(dtime[3:5]), int(dtime[11:13]),
                                 int(dtime[14:16]), int(dtime[17:19]), int(float(dtime[19:23]) * 1000000))
                seconds = fdate.strftime('%s.%f')

                import sys
                if sys.platform.startswith('win'):
                    csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(
                        dateob.strftime('%H:%M:%S')) + ',' + str(
                        int(pd.to_numeric(dateob.strftime('%S.%f')))) + ',' + str(
                        pd.to_numeric(dateob.strftime('%f')) * 1000) + str(',')
                    csvWrite += str('50') + ',' + str('0') + ',' + str('0') + ',' + str('0') + ',' + str(
                        lstS[4]) + ',' + str('0') + ',' + str(lstS[4]) + ','
                    csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(
                        lstS[3]) + ',' + str(lstS[4]) + ',' + str(lstS[5]) + ',' + str(lstS[6]) + ','
                    csvWrite += str(lstS[7]) + ',' + str(lstS[8]) + ',' + str(lstS[9]) + ',' + str(
                        lstS[10]) + ',' + str(lstS[11]) + ',' + str(lstS[12]) + ',' + str(lstS[13]) + str(
                        ',') + str(lstS[14])
                if not sys.platform.startswith('win'):
                    csvWrite = str(dateob.strftime('%Y-%m-%d')) + ',' + str(
                        dateob.strftime('%H:%M:%S')) + ',' + str(seconds[:10]) + ',' + str(
                        pd.to_numeric(seconds[11:]) * 1000) + str(',')
                    csvWrite += str(lstS[20]) + ',' + str(lstS[15]) + ',' + str(lstS[16]) + ',' + str(
                        lstS[17]) + ',' + str(
                        lstS[4]) + ',' + str('0') + ',' + str(lstS[4]) + ','
                    csvWrite += str('0') + ',' + str(lstS[2]) + ',' + str(lstS[1]) + ',' + str(
                        lstS[3]) + ',' + str(lstS[4]) + ',' + str(lstS[5]) + ',' + str(lstS[6]) + ','
                    csvWrite += str(lstS[7]) + ',' + str(lstS[8]) + ',' + str(lstS[9]) + ',' + str(
                        lstS[10]) + ',' + str(lstS[11]) + ',' + str(lstS[12]) + ',' + str(lstS[13]) + str(
                        ',') + str(lstS[14]) + '\n'
                fOut.write(csvWrite)
                xCntObs += 1
            infOut.write(str(xFilename) + '\n')
        fOut.close()
        fLog.close()
        infOut.close()
        #print(xCar + "\t" + xdat + "\t" + fnOut[-22:] + "\t" + str(xCntObs) + "\t" + str(xCntGoodValues) + "\t" + str(
        #    gZIP))
        print(f"{xCar} \t {xdat} \t {fnOut[-(17 + len(xCar)):]} \t  {xCntObs} \t {xCntGoodValues} \t {gZIP}")

        def calcVel(timediff, distance):
            if timediff == 0:
                return (0)
            elif timediff != 0:
                return (distance / timediff)

        wind_df = pd.read_csv(fnOutTemp)
        wind_df_not_null = wind_df.loc[wind_df['LAT'].notnull(),].reset_index(drop=True)
        del (wind_df)
        wind_df = wind_df_not_null.copy()

        radians = False
        wind_df['QUADRANT'] = wind_df.apply(lambda row: getQuad(row['U'], row['V']), axis=1)
        wind_df['secnan'] = wind_df.apply(lambda row: row['SECONDS'] + row['NANOSECONDS'] * 1e-9,
                                          axis=1)  # + row['NANOSECONDS']*1e-9,axis=1)
        wind_df['prev_LAT'] = wind_df.LAT.shift(periods=1)
        wind_df['next_LAT'] = wind_df.LAT.shift(periods=-1)
        wind_df['prev_LONG'] = wind_df.LONG.shift(periods=1)
        wind_df['next_LONG'] = wind_df.LONG.shift(periods=-1)
        wind_df['prev_TIME'] = wind_df.secnan.shift(periods=1)
        wind_df['next_TIME'] = wind_df.secnan.shift(periods=-1)
        wind_df['distance'] = wind_df.apply(
            lambda row: haversine(row['prev_LAT'], row['prev_LONG'], row['next_LAT'], row['next_LONG']), axis=1)
        wind_df['bearing'] = wind_df.apply(
            lambda row: calcBearing(row['prev_LAT'], row['next_LAT'], row['prev_LONG'], row['next_LONG'], radians),
            axis=1)
        wind_df['timediff'] = wind_df.apply(lambda row: row['next_TIME'] - row['prev_TIME'], axis=1)
        # wind_df['VELOCITY_calc'] = wind_df.apply(lambda row:calcVel(row['timediff'],row['distance']),axis=1)
        wind_df['U_cor'] = wind_df.apply(lambda row: row['U'] + row['VELOCITY'], axis=1)
        wind_df['horz_length'] = wind_df.apply(lambda row: np.sqrt(row['U_cor'] ** 2 + row['V'] ** 2), axis=1)
        wind_df['uncor_theta'] = wind_df.apply(
            lambda row: calcBearing(row['U_cor'], row['V'], row['QUADRANT'], row['horz_length'], radians), axis=1)
        wind_df['adj_theta'] = wind_df.apply(lambda row: (row['uncor_theta'] + row['bearing']) % 360, axis=1)
        wind_df['totalWind'] = wind_df.apply(lambda row: np.sqrt(row['horz_length'] ** 2 + row['W'] ** 2), axis=1)
        wind_df['phi'] = wind_df.apply(lambda row: np.arctan(row['horz_length']), axis=1)
        wind_df['shift_CH4'] = wind_df.CH4.shift(periods=int(float(shift)))
        wind_df['raw_CH4'] = wind_df.apply(lambda row: row['BCH4'], axis=1)
        wind_df['BCH4'] = wind_df.loc[:, ['shift_CH4']]
        wind_df['CH4'] = wind_df.loc[:, ['shift_CH4']]
        wind_df['TCH4'] = wind_df.loc[:, ['shift_CH4']]

        wind_df2 = wind_df[wind_df.CH4.notnull()]
        wind_df2 = wind_df.copy()
        wind_df3 = wind_df2.drop(
            ['QUADRANT', 'secnan', 'prev_LAT', 'next_LAT', 'prev_LONG', 'next_LONG', 'prev_TIME', 'next_TIME',
             'distance', 'timediff', 'uncor_theta', 'CH4'], axis=1)
        wind_df3['CH4'] = wind_df3.loc[:, 'shift_CH4']
        wind_df3 = wind_df3.drop(['shift_CH4'], axis=1)

        wind_df3 = wind_df3.loc[:,
                   ['DATE', 'TIME', 'SECONDS', 'NANOSECONDS', 'VELOCITY', 'U', 'V', 'W', 'BCH4', 'BRSSI', 'TCH4',
                    'TRSSI', 'PRESS_MBAR', 'INLET' \
                       , 'TEMPC', 'CH4', 'H20', 'C2H6', 'R', 'C2C1', 'BATTV', 'POWMV', 'CURRMA', 'SOCPER', 'LAT',
                    'LONG', 'bearing', 'U_cor', \
                    'horz_length', 'adj_theta', 'totalWind', 'phi', 'raw_CH4']]
        # wind_df4 = wind_df3.loc[wind_df3.totalWind.notnull(),:]

        wind_df4 = wind_df3.copy()

        # wind_df7 = addOdometer(wind_df4,'LAT','LONG')

        # wind_df4 = wind_df7.copy()
        wind_df5 = wind_df4.loc[wind_df4.VELOCITY > xMinCarSpeed, :]
        wind_df6 = wind_df5.loc[wind_df5.VELOCITY < xMaxCarSpeed, :]

        del (wind_df4)
        wind_df4 = wind_df6.copy().drop_duplicates()
        wind_df5 = wind_df4.loc[wind_df4.CH4.notnull(), :]
        wind_df4 = wind_df5.copy()
        if bFirst:
            wind_df4.to_csv(fnOut, index=False)
        elif not bFirst:
            norm = pd.read_csv(fnOut)
            pd.concat([norm, wind_df4]).sort_values(by='SECONDS').reset_index(drop=True).to_csv(fnOut, index=False)
        os.remove(fnOutTemp)
        return True
    except ValueError:
        return False


########################################################################
#### addOdometer
# function to add column to dataframe with Odometer reading (in kms)
# Input: df, lat, lon
# Output: df with the odometer reading

def addOdometer(df, lat, lon):
    import pandas as pd
    import math
    df_use = df.loc[:, [(lat), (lon)]]
    df_use['prev_LAT'] = df_use.loc[:, (lat)].shift(periods=1)
    df_use['prev_LON'] = df_use.loc[:, (lon)].shift(periods=1)
    df_use['distance2'] = df_use.apply(lambda row: haversine(row['prev_LAT'], row['prev_LON'], row[(lat)], row[(lon)]),
                                       axis=1)
    df_use = df_use.reset_index(drop=True)

    def nanthing(thing):
        if (math.isnan(thing) == True):
            return (0)
        else:
            return (thing)

    df_use.loc[:, 'distance'] = df_use.apply(lambda x: nanthing(x.distance2), axis=1)
    df_use['prev_dist'] = df_use.loc[:, 'distance'].shift(periods=1)
    # df_use['od'] = df
    df_use['odometer'] = df_use['distance'].cumsum()
    df_use['prevod'] = df_use.loc[:, 'odometer'].shift(periods=1)
    df_use['dif'] = df_use.apply(lambda x: x.odometer - x.prevod, axis=1)
    df_use['dif'] = df_use.apply(lambda x: nanthing(x.dif), axis=1)
    return (pd.merge(df, df_use.loc[:, [(lat), (lon), 'odometer', 'distance']], on=[(lat), (lon)]))


########################################################################
#### strList
# helper function to convert a string of a list to just a list
# Input: string of a list thing to
# Output: list

def strList(x):
    import ast
    x = ast.literal_eval(x)
    x = [n.strip() for n in x]
    return (x)


########################################################################
#### COUNTTIMES
# Input: a list of peak times included in a given combined peak
# Output: counts number of times the peak was seen (not in same 5 min period)


def countTimes(opList,xCar):
    if isinstance(opList, str):
        opList = strList(opList)
    if len(opList) == 1:
        numtimes = 1
        return (numtimes)
    else:
        opList.sort()
        numtimes = 1
        index = 1
        for x in opList:
            if index == 1:
                initTime = float(x[len(xCar)+1:])
                initStart = initTime
                initEnd = initTime + 300
                index = index + 1
            if index != 1:
                curTime = float(x[len(xCar)+1:])
                within = curTime < initEnd
                if curTime < initEnd:
                    index = index + 1
                elif curTime >= initEnd:
                    numtimes = numtimes + 1
                    initTime = curTime
                    initStart = initTime
                    initEnd = curTime + 300
                    index = index + 1
        return (numtimes)


def IdentifyPeaks(xCar, xDate, xDir, xFilename, outDir, processedFileLoc, Engineering, threshold='.1',
                  xTimeThreshold='5.0', minElevated='2', xB='1020', basePerc='50'):
    import csv, numpy
    import geopandas as gpd
    import shutil
    try:
        baseCalc = float(basePerc)
        xABThreshold = float(threshold)
        minElevated = float(minElevated)
        # xABThreshold = 0.1                 # above baseline threshold above the mean value
        xDistThreshold = 160.0  # find the maximum CH4 reading of observations within street segments of this grouping distance in meters
        xSDF = 4  # multiplier times standard deviation for floating baseline added to mean
        # xB = 1020       # the number of records that constitutes the floating baseline time -- 7200 = 1 hour (assuming average of 0.5 seconds per record)
        # xB = 102 # since it is 1 record/second

        xB = int(xB)
        # xB = 300 #five min?
        xTimeThreshold = float(xTimeThreshold)

        # fn = xDir + "/" + xFilename      #set raw text file to read in
        fn = xDir + xFilename  # set raw text file to read in

        fnOut = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-",
                                                                    "") + ".csv"  # set CSV format output for observed peaks for a given car, day, city
        fnShape = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-", "") + ".shp"
        fnLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-",
                                                                    "") + ".log"  # set CSV output for observed peaks for a given car, day, city
        pkLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-",
                                                                    "") + "_info.csv"  # set CSV output for observed peaks for a given car, day, city

        jsonOut = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-",
                                                                      "") + ".json"  # set CSV format output for observed peaks for a given car, day, city

        infOut = processedFileLoc + xCar + "_" + xDate.replace("-", "") + "_info.csv"

        ### TEST THING
        fn = xDir + xFilename  # set raw text file to read in

        fnOut = outDir + "Peaks" + "_" + xCar + "_" + xDate +  ".csv"  # set CSV format output for observed peaks for a given car, day, city
        fnShape = outDir + "Peaks" + "_" + xCar + "_" + xDate + ".shp"
        fnLog = outDir + "Peaks" + "_" + xCar + "_" + xDate + ".log"  # set CSV output for observed peaks for a given car, day, city
        pkLog = outDir + "Peaks" + "_" + xCar + "_" + xDate + "_info.csv"  # set CSV output for observed peaks for a given car, day, city

        jsonOut = outDir + "Peaks" + "_" + xCar + "_" + xDate + ".json"  # set CSV format output for observed peaks for a given car, day, city

        infOut = processedFileLoc + xCar + "_" + xDate + "_info.csv"


        #print(str(outDir + "Peaks" + "_" + xCar + "_" + xDate + "_info.csv"))
        print(f"{outDir}Peaks_{xCar}_{xDate}_info.csv")
        fLog = open(fnLog, 'w')
        shutil.copy(infOut, pkLog)

        # field column indices for various variables
        if Engineering == True:
            fDate = 0;
            fTime = 1;
            fEpochTime = 2;
            fNanoSeconds = 3;
            fVelocity = 4;
            fU = 5;
            fV = 6;
            fW = 7;
            fBCH4 = 10;
            fBCH4 = 8;
            fBRSSI = 9;
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
        elif not Engineering:
            fDate = 0;
            fTime = 1;
            fEpochTime = 2;
            fNanoSeconds = 3;
            fVelocity = 4;
            fU = 5;
            fV = 6;
            fW = 7;
            fBCH4 = 8;

            fBRSSI = 9;
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
            fDist = 33;
            fOdometer = 34;

            # read data in from text file and extract desired fields into a list, padding with 5 minute and hourly average
            x1 = [];
            x2 = [];
            x3 = [];
            x4 = [];
            x5 = [];
            x6 = [];
            x7 = [];
            x8 = []
            x9 = []
            x10 = []

            count = -1
            with open(fn, 'r') as f:
                t = csv.reader(f)
                for row in t:
                    woo = row
                    # print(count)
                    if count < 0:
                        count += 1
                        continue

                    datet = row[fDate].replace("-", "") + row[fTime].replace(":", "")
                    ## if not engineering
                    epoch = float(row[fEpochTime] + "." + row[fNanoSeconds][0])

                    # x1.append(float(epoch));

                    # =============================================================================
                    #                 x1.append(float(str(row[fEpochTime]) + '.' + str(row[fNanoSeconds])));
                    #                 x2.append(float(int(datet)));
                    #                 x3.append(float(row[fLat]));
                    #                 x4.append(float(row[fLon]));
                    #                 x5.append(float(row[fBCH4]));
                    #                 x6.append(float(row[fTCH4]))
                    #                 x7.append(0.0);
                    #                 x8.append(0.0)
                    # =============================================================================

                    datetime = row[fDate].replace("-", "") + row[fTime].replace(":", "")

                    x1.append(epoch);
                    x2.append(datetime);
                    if row[fLat] == '':
                        x3.append('')
                    elif row[fLat] != '':
                        x3.append(float(row[fLat]));
                    if row[fLon] == '':
                        x4.append('')
                    elif row[fLon] != '':
                        x4.append(float(row[fLon]));

                    x5.append(float(row[fBCH4]));
                    x6.append(float(row[fTCH4]))
                    x7.append(0.0);
                    x8.append(0.0)
                    x9.append(row[fOdometer])

                    # print (str(row[fLat])+ str(row[1]))
                    count += 1
            print(f"Number of observations processed:{count}")

        # convert lists to numpy arrays
        aEpochTime = numpy.array(x1);
        aDateTime = numpy.array(x2);
        aLat = numpy.array(x3);
        aLon = numpy.array(x4);
        aCH4 = numpy.array(x5);
        aTCH4 = numpy.array(x6)
        aMean = numpy.array(x7);
        aThreshold = numpy.array(x8)
        aOdom = numpy.array(x9)
        xLatMean = numpy.mean(aLat)
        xLonMean = numpy.mean(aLon)

        fLog.write("Day CH4_mean = " + str(numpy.mean(aCH4)) + ", Day CH4_SD = " + str(numpy.std(aCH4)) + "\n")
        fLog.write("Center lon/lat = " + str(xLonMean) + ", " + str(xLatMean) + "\n")
        # pkLog.write('hi')
        lstCH4_AB = []

        # generate list of the index for observations that were above the threshold
        for i in range(0, count - 2):
            if ((count - 2) > xB):
                topBound = min((i + xB), (count - 2))
                botBound = max((i - xB), 0)

                for t in range(min((i + xB), (count - 2)), i, -1):
                    if aEpochTime[t] < (aEpochTime[i] + (xB / 2)):
                        topBound = t
                        break
                for b in range(max((i - xB), 0), i):
                    if aEpochTime[b] > (aEpochTime[i] - (xB / 2)):
                        botBound = b
                        break

                xCH4Mean = numpy.percentile(aCH4[botBound:topBound], baseCalc)
            # xCH4SD = numpy.std(aCH4[botBound:topBound])
            else:
                xCH4Mean = numpy.percentile(aCH4[0:(count - 2)], baseCalc)
                # xCH4SD = numpy.std(aCH4[0:(count-2)])
            xThreshold = xCH4Mean + (xCH4Mean * xABThreshold)

            if (aCH4[i] > xThreshold):
                lstCH4_AB.append(i)
                aMean[
                    i] = xCH4Mean  # insert mean + SD as upper quartile CH4 value into the array to later retreive into the peak calculation
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
                xLon1 = aLon[i];
                xLat1 = aLat[i]
                xOdom = aOdom[i]
            else:
                # calculate distance between points
                xDist = haversine(xLat1, xLon1, aLat[i], aLon[i])
                xDistPeak += xDist
                xCH4Peak += (xDist * (aCH4[i] - aMean[i]))
                xLon1 = aLon[i];
                xLat1 = aLat[i]
                xOdom = aOdom[i]
                if (sID == ""):
                    xTime = aEpochTime[i]
                    sID = str(xCar) + "_" + str(xTime)
                    sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1)))  # 30 sec
                if ((aEpochTime[i] - aEpochTime[prevIndex]) > xTimeThreshold):  # initial start of a observed peak
                    cntPeak += 1
                    xTime = aEpochTime[i]
                    xDistPeak = 0.0
                    xCH4Peak = 0.0
                    sID = str(xCar) + "_" + str(xTime)
                    sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1)))  # 30 sec
                    # print str(i) +", " + str(xDist) + "," + str(cntPeak) +"," + str(xDistPeak)
                lstCH4_ABP.append(
                    [sID, xTime, aEpochTime[i], aDateTime[i], aCH4[i], aLon[i], aLat[i], aMean[i], aThreshold[i],
                     xDistPeak, xCH4Peak, aTCH4[i], sPeriod5Min, xOdom])
            cnt += 1
            prevIndex = i

        # Finding peak_id larger than 160.0 m
        tmpsidlist = []
        for r in lstCH4_ABP:
            if (float(r[9]) > 160.0) and (r[0] not in tmpsidlist):
                tmpsidlist.append(r[0])
        cntPeak -= len(tmpsidlist)

        fLog.write("Number of peaks found: " + str(cntPeak) + "\n")
        #print(xCar + "\t" + xDate + "\t" + xFilename + "\t" + str(count) + "\t" + str(len(lstCH4_ABP)))
        print(f"{xCar} \t {xDate} \t {xFilename} \t {count} \t {len(lstCH4_ABP)}")

        #### calculate attribute for the area under the curve -- PPM

        # write out the observed peaks to a csv to be read into a GIS
        fOut = open(fnOut, 'w')
        # s = "PEAK_NUM,EPOCHSTART,EPOCH,DATETIME,CH4,LON,LAT,CH4_BASELINE,CH4_THRESHOLD,PEAK_DIST_M,PEAK_CH4,TCH4,PERIOD5MIN\n"
        s = "OP_NUM,OP_EPOCHSTART,OB_EPOCH,OB_DATETIME,OB_CH4,OB_LON,OB_LAT,OB_CH4_BASELINE,OB_CH4_THRESHOLD,OP_PEAK_DIST_M,OP_PEAK_CH4,OB_TCH4,OB_PERIOD5MIN,ODOMETER\n"

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
            pkDistDf = openFile.copy().groupby('OP_NUM', as_index=False).apply(
                lambda x: max(x.ODOMETER) - min(x.ODOMETER))
            pkDistDf.columns = ['OP_NUM', 'OP_DISTANCE']
            openFile = pd.merge(openFile.copy(), pkDistDf)

            tempCount = openFile.groupby('OP_NUM', as_index=False).OP_EPOCHSTART.count().rename(
                columns={'OP_EPOCHSTART': 'Frequency'})
            tempCount = tempCount.loc[tempCount.Frequency >= minElevated, :]
            if tempCount.shape[0] == 0:
                print(f"No Observed Peaks with enough Elevated Readings Found in the file: {xFilename}")

            elif tempCount.shape[0] != 0:
                oFile = pd.merge(openFile, tempCount, on=['OP_NUM'])
                openFile = oFile.copy()
                del (oFile)
                openFile['minElevated'] = openFile.apply(lambda x: int(minElevated), axis=1)
                openFile.to_csv(fnOut, index=False)
                openFile['OB_CH4_AB'] = openFile.loc[:, 'OB_CH4'].sub(openFile.loc[:, 'OB_CH4_BASELINE'], axis=0)

                fileWt = weightedLoc(openFile, 'OB_LAT', 'OB_LON', 'OP_NUM', 'OB_CH4_AB').loc[:, :].rename(
                    columns={'OB_LAT': 'pk_LAT', 'OB_LON': 'pk_LON'}).reset_index(drop=True)
                geometry_temp = [Point(xy) for xy in zip(fileWt['pk_LON'], fileWt['pk_LAT'])]
                crs = {'init': 'epsg:4326'}

                # geometry is the point of the lat/lon
                # gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

                ## BUFFER AROUND EACH 'OP_NUM' OF 30 M
                gdf_buff = gpd.GeoDataFrame(fileWt, crs=crs, geometry=geometry_temp)
                # gdf_buff = makeGPD(datFram,'LON','LAT')
                gdf_buff = gdf_buff.to_crs(epsg=32610)
                gdf_buff['geometry'] = gdf_buff.loc[:, 'geometry'].buffer(30)
                gdf_buff.to_file(jsonOut, driver="GeoJSON")
        elif openFile.shape[0] == 0:
            print(f"No Observed Peaks Found in the file:{xFilename}")
    except ValueError:
        print("Error in Identify Peaks")
        return False


########################################################################
#### FILTER PEAK
# Input: a .csv file with peak data (already have gone through
#   ('identifyPeaks')
# Output: checks for overlaps within a day's drive. Finds locations of
#           Observed Peaks

def filterPeak(xCar, xDate, xDir, xFilename, outFolder, buffer='30', whichpass=0):
    ## NECESSARY MODULES
    import pandas as pd  #
    import geopandas as gpd
    import shutil
    from datetime import datetime
    from shapely.geometry import Point
    buffer = float(buffer)

    # MOVING THE FILES NECESSARY & CREATING NEW FILES
    file_loc = xDir + xFilename
    new_loc = outFolder + "Filtered" + xFilename
    new_loc_json = new_loc[:-3] + 'json'

    oldInfo = xDir + 'Peaks_' + xCar + "_" + xDate.replace("-", "") + "_info.csv"
    newInfo = outFolder + 'FilteredPeaks_' + xCar + "_" + xDate.replace("-", "") + "_info.csv"

    shutil.copy(oldInfo, newInfo)

    # identified peaks has the columns:
    # Index(['OP_NUM', 'OP_EPOCHSTART', 'OB_EPOCH', 'OB_DATETIME', 'OB_CH4',
    #   'OB_LON', 'OB_LAT', 'OB_CH4_BASELINE', 'OB_CH4_THRESHOLD',
    #   'OP_PEAK_DIST_M', 'OP_PEAK_CH4', 'OB_TCH4', 'OB_PERIOD5MIN'],
    #  dtype='object')

    datFram = pd.read_csv(file_loc)  # READING IN THE FILE

    if datFram.shape[0] == 0:  # IF THE OBSERVED PEAK FILE WAS EMPTY, MOVE ON
        print("Not filtering this file, no peak in it!")
    elif datFram.shape[0] == 1:  ## IF ONLY HAD ONE OBSERVED PEAK
        datFram_cent = datFram.copy()
        datFram_cent['OB_CH4_AB'] = datFram.loc[:, 'OB_CH4'].sub(datFram.loc[:, 'OB_CH4_BASELINE'], axis=0)
        maxch4 = datFram_cent.groupby('OP_NUM', as_index=False).OB_CH4_AB.max().rename(
            columns={'OB_CH4_AB': 'pk_maxCH4_AB'})
        datFram_wtLoc = weightedLoc(datFram_cent, 'OB_LAT', 'OB_LON', 'OP_NUM', 'OB_CH4_AB').loc[:, :].rename(
            columns={'OB_LAT': 'pk_LAT', 'OB_LON': 'pk_LON'})
        datFram_wtLocMax = pd.merge(datFram_wtLoc, maxch4, on=['OP_NUM'])
        pass_info = datFram.copy()
        geometry_temp = [Point(xy) for xy in zip(datFram_wtLocMax['pk_LON'], datFram_wtLocMax['pk_LAT'])]
        crs = {'init': 'epsg:4326'}
        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)
        gdf_buff = gdf_buff.to_crs(epsg=32610)
        # gdf_buff['geometry'] = gdf_buff.loc[:,'geometry'].buffer(30)
        gdf_buff['geometry'] = gdf_buff.loc[:, 'geometry'].buffer(buffer)
        gdf_tog = pd.merge(gdf_buff, datFram, on=['OP_NUM'])
        gdf_bind_pks = gdf_buff.copy()
        gdf_pass_pks = gdf_bind_pks.copy()
        gdf_pass_pks['min_read'] = gdf_pass_pks.loc[:, 'OP_NUM']
        gdf_pass_pks['numtimes'] = 1
        gdf_pass_pks['numdays'] = 1
        gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:, 'geometry']
        gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:, ['OP_NUM']].to_numpy())].copy()
        gdf_pass_pks['dates'] = gdf_pass_pks.apply(
            lambda x: datetime.fromtimestamp(int(x.OP_NUM[len(xCar) + 1:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'), axis=1)
        gdf_pass_pks['pk_Dates'] = [list(x) for x in list(gdf_pass_pks.loc[:, ['dates']].to_numpy())]
        gdf_pass_pks['min_Date'] = gdf_pass_pks.loc[:, 'dates']
        gdf_pass_pks = gdf_pass_pks.drop(columns=['dates'])

        gdf_pass_pks['verified'] = False
        gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:, "newgeo"]
        together = pd.merge(gdf_pass_pks, gdf_tog, on=['OP_NUM', 'pk_LON', 'pk_LAT', 'pk_maxCH4_AB', 'geometry'])
        together['pass'] = whichpass
        gdf_pass_pks = together.copy()

        gdf_pass_pks['pkGEO'] = gdf_pass_pks.loc[:, "geometry"]
        gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:, "newgeo"]
        del (gdf_pass_pks['newgeo'])
        gdf_pass_pks['pass'] = whichpass

        gdf_op_unique = gdf_pass_pks.loc[:,
                        ['numtimes', 'min_read', 'numdays', 'min_Date', 'verified', 'pass', 'OB_LON',
                         'OB_LAT']].drop_duplicates()
        gdfcop = gdf_pass_pks.loc[:,
                 ['OP_NUM', 'min_read', 'min_Date', 'numtimes', 'verified', 'pass', 'pk_LAT', 'pk_LON',
                  'pk_maxCH4_AB']].drop_duplicates()
        combinedOP = weightedLoc(gdfcop, 'pk_LAT', 'pk_LON', 'min_read', 'pk_maxCH4_AB').loc[:, :].rename(
            columns={'pk_LAT': 'Overall_LAT', 'pk_LON': 'Overall_LON'}).reset_index(drop=True)
        combinedOP1 = pd.merge(combinedOP, gdfcop, on=['min_read'])
        geometry_temp = [Point(xy) for xy in zip(combinedOP1['Overall_LON'], combinedOP1['Overall_LAT'])]
        crs = {'init': 'epsg:4326'}
        gdf_OP = gpd.GeoDataFrame(combinedOP1, crs=crs, geometry=geometry_temp)
        gdf_OP = gdf_OP.to_crs(epsg=32610).copy()
        gdf_OP_reduced = gdf_OP.loc[:, ['min_read', 'geometry', 'numtimes', 'Overall_LON', 'Overall_LAT', 'min_Date',
                                        'verified']].drop_duplicates().reset_index(drop=True)
        gdf_OP_reduced.to_file(new_loc_json, driver="GeoJSON")
        gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']), gdf_pass_pks.drop(columns=['geometry']),
                                     on=['min_read', 'min_Date', 'numtimes', 'pass', 'verified', 'pk_LAT', 'pk_LON',
                                         'OP_NUM', 'pk_maxCH4_AB'])
        gdf_OP_wrecombine.to_csv(new_loc, index=False)

        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)
        unique_peaks = gdf_pass_pks.loc[:, ['OP_NUM', 'pk_LAT', 'pk_LON', 'min_read', 'min_Date']].drop_duplicates()
        unique_peaks['save'] = True
        good_pks = list(unique_peaks.index)

        def getthing(index):
            if index in good_pks:
                return True
            else:
                return False

        gdf_pass_pks['wooind'] = gdf_pass_pks.index
        gdf_pass_pks['save'] = gdf_pass_pks.apply(lambda x: getthing(x.wooind), axis=1)
        unique_pks_tog = gdf_pass_pks.loc[gdf_pass_pks.save == True, :].reset_index(drop=True)
        unique_pks_tog['Latitude'] = unique_pks_tog.loc[:, 'pk_LAT']
        unique_pks_tog['Longitude'] = unique_pks_tog.loc[:, 'pk_LON']
        unique_pks_tog.to_csv(new_loc, index=False)

        # return(gdf_OP_wrecombine)

    elif datFram.shape[0] != 1:
        datFram_cent = datFram.copy()
        datFram_cent['OB_CH4_AB'] = datFram.loc[:, 'OB_CH4'].sub(datFram.loc[:, 'OB_CH4_BASELINE'], axis=0)

        ### MAXCH4 is a df with the max methane (above baseline) in the given observed peak
        # maxch4 = datFram_cent.groupby('PEAK_NUM',as_index = False).CH4_AB.max().rename(columns = {'CH4_AB':'pk_maxCH4_AB'})
        maxch4 = datFram_cent.groupby('OP_NUM', as_index=False).OB_CH4_AB.max().rename(
            columns={'OB_CH4_AB': 'pk_maxCH4_AB'})

        ### FINDING WEIGHTED LOCATION OF THE OP, BY THE ABOVE BASELINE CH4 LEVEL
        # wtloc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB')
        # datFram_wtLoca =  wtloc.copy()
        # datFram_wtLoc = datFram_wtLoca.rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'})

        datFram_wtLoc = weightedLoc(datFram_cent, 'OB_LAT', 'OB_LON', 'OP_NUM', 'OB_CH4_AB').loc[:, :].rename(
            columns={'OB_LAT': 'pk_LAT', 'OB_LON': 'pk_LON'})
        # datFram_wtLoc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB').rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'}).copy()
        datFram_wtLocMax = pd.merge(datFram_wtLoc, maxch4, on=['OP_NUM'])

        pass_info = datFram.copy()

        ## MIGHT NEED TO CHANGE BACK
        # geometry_temp = [Point(xy) for xy in zip(datFram['LON'], datFram['LAT'])]
        # crs = {'init': 'epsg:4326'}

        geometry_temp = [Point(xy) for xy in zip(datFram_wtLocMax['pk_LON'], datFram_wtLocMax['pk_LAT'])]
        crs = {'init': 'epsg:4326'}

        # geometry is the point of the lat/lon
        # gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        ## BUFFER AROUND EACH 'OP_NUM' OF 30 M
        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)
        # gdf_buff = makeGPD(datFram,'LON','LAT')
        gdf_buff = gdf_buff.to_crs(epsg=32610)
        # gdf_buff['geometry'] = gdf_buff.loc[:,'geometry'].buffer(30)
        gdf_buff['geometry'] = gdf_buff.loc[:, 'geometry'].buffer(buffer)
        gdf_tog = pd.merge(gdf_buff, datFram, on=['OP_NUM'])
        gdf_bind_pks = gdf_buff.copy()

        if gdf_bind_pks.shape[0] > 1:
            data_overlap = gpd.GeoDataFrame(crs=gdf_bind_pks.crs)
            data_temp = gdf_bind_pks.copy()
            data_temp = data_temp.to_crs(epsg=32610)

            for index, row in data_temp.iterrows():
                data_temp1 = data_temp.loc[data_temp.OP_NUM != row.OP_NUM, :]
                data_temp1 = data_temp1.to_crs(epsg=32610)

                # check if intersection occured
                overlaps = data_temp1[data_temp1.geometry.overlaps(row.geometry)]['OP_NUM'].tolist()
                if len(overlaps) > 0:

                    # compare the area with threshold
                    for y in overlaps:
                        temp_area = gpd.overlay(data_temp.loc[data_temp.OP_NUM == y,],
                                                data_temp.loc[data_temp.OP_NUM == row.OP_NUM,], how='intersection')
                        temp_area = temp_area.loc[temp_area.geometry.area >= 0.001]
                        if temp_area.shape[0] > 0:
                            temp_union = gpd.overlay(data_temp.loc[data_temp.OP_NUM == y,],
                                                     data_temp.loc[data_temp.OP_NUM == row.OP_NUM,], how='union')
                            data_overlap = gpd.GeoDataFrame(pd.concat([temp_union, data_overlap], ignore_index=True),
                                                            crs=data_temp.crs)
            if data_overlap.size > 0:

                firstnull2 = data_overlap.loc[data_overlap.OP_NUM_1.isnull(), :]
                firstnull = firstnull2.copy()
                firstnull.loc[:, 'OP_NUM_1'] = firstnull2.loc[:, 'OP_NUM_2']

                secnull2 = data_overlap.loc[data_overlap.OP_NUM_2.isnull(), :]

                secnull = secnull2.copy()
                secnull.loc[:, 'OP_NUM_2'] = secnull2.loc[:, 'OP_NUM_1']

                withoutNA = data_overlap.copy().dropna()
                allTog2 = pd.concat([firstnull, secnull, withoutNA]).reset_index().copy()

                allTog2['notsame'] = allTog2.apply(lambda x: x.OP_NUM_1 == x.OP_NUM_2, axis=1)
                allTog = allTog2.loc[allTog2.notsame == False, :].drop(columns=['notsame'])

                over = allTog.copy()
                over['sorted'] = over.apply(lambda y: sorted([y['OP_NUM_1'], y['OP_NUM_2']]), axis=1)
                over['sorted'] = over.sorted.apply(lambda y: ''.join(y))
                over = over.drop_duplicates('sorted')
                over['combined'] = [list(x) for x in list(over.loc[:, ['OP_NUM_1', 'OP_NUM_2']].to_numpy())]
                # over['date1'] = over.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM_1'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                # over['date2'] = over.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM_2'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                over['date1'] = over.apply(
                    lambda x: datetime.fromtimestamp(int(x.OP_NUM_1[len(xCar)+1:x.OP_NUM_1.find('.')])).strftime('%Y-%m-%d'),
                    axis=1)
                over['date2'] = over.apply(
                    lambda x: datetime.fromtimestamp(int(x.OP_NUM_2[len(xCar)+1:x.OP_NUM_2.find('.')])).strftime('%Y-%m-%d'),
                    axis=1)

                def unique(list1):
                    # intilize a null list
                    unique_list = []
                    # traverse for all elements
                    for x in list1:
                        # check if exists in unique_list or not
                        if x not in unique_list:
                            unique_list.append(x)
                    return (unique_list)

                over['dates'] = [list(x) for x in list(over.loc[:, ['date1', 'date2']].to_numpy())]
                over['pk_Dates'] = over.apply(lambda x: unique(x.dates), axis=1)
                over = over.drop(columns=['dates'])

                over['VER_NUM'] = over.apply(lambda y: y.combined, axis=1)
                over['min_val'] = over.apply(lambda y: min(y.combined), axis=1)
                over2 = over.reset_index().loc[:,
                        ['OP_NUM_1', 'OP_NUM_2', 'geometry', 'combined', 'min_val', 'pk_Dates']]

                overcop = over2.copy().rename(columns={'combined': 'recombine'})
                # overcop.loc[:,'recombine'] = overcop.loc[:,'combined']

                for index, row in overcop.iterrows():
                    united = row.recombine
                    undate = row.pk_Dates
                    for index2, row2 in overcop.iterrows():
                        united_temp = unIfInt(united, row2.recombine)
                        undate_temp = unIfInt(undate, row2.pk_Dates)
                        if united_temp != None:
                            united = united_temp
                        if undate_temp != None:
                            undate = undate_temp
                    overcop.at[index, 'recombine'] = united.copy()
                    overcop.at[index, 'pk_Dates'] = undate.copy()

                    del (united)
                    del (undate)

                overcop['recombine'] = overcop.apply(lambda y: sorted(y.recombine), axis=1).copy()
                overcop['pk_Dates'] = overcop.apply(lambda y: sorted(y.pk_Dates), axis=1).copy()
                overcop['min_read'] = overcop.apply(lambda y: min(y.recombine), axis=1).copy()
                overcop['min_Date'] = overcop.apply(lambda y: min(y.pk_Dates), axis=1).copy()

                newOverlap = overcop.dissolve(by='min_read', as_index=False).loc[:,
                             ['min_read', 'geometry', 'recombine', 'min_Date', 'pk_Dates']].copy()

                combined = gdf_bind_pks.copy()
                combined['recombine'] = [list(x) for x in list(combined.loc[:, ['OP_NUM']].to_numpy())]
                # combined['dates'] = combined.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                combined['dates'] = combined.apply(
                    lambda x: datetime.fromtimestamp(int(x.OP_NUM[len(xCar)+1:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'), axis=1)

                combined['pk_Dates'] = [list(x) for x in list(combined.loc[:, ['dates']].to_numpy())]
                combined['min_Date'] = combined.loc[:, 'dates']
                combined['numtimes'] = 1
                combined['newgeo'] = combined.loc[:, 'geometry']
                combined['min_read'] = combined.loc[:, "OP_NUM"]

                for index, row in combined.iterrows():
                    for index2, row2 in newOverlap.iterrows():
                        if row.OP_NUM in row2.recombine:
                            combined.at[index, 'recombine'] = row2.recombine.copy()
                            # combined.at[index, 'newgeo']  = row2.copy().geometry
                            combined.at[index, 'min_read'] = row2.copy().min_read
                            combined.at[index, 'pk_Dates'] = row2.pk_Dates
                            combined.at[index, 'min_Date'] = row2.min_Date

                # combined['numtimes'] = combined.apply(lambda y: len(y.recombine),axis = 1).copy()
                combined['numtimes'] = combined.apply(lambda x: countTimes(x.recombine,xCar), axis=1)

                combined['numdays'] = combined.apply(lambda y: len(y.pk_Dates), axis=1).copy()
                combined_reduced = combined.loc[:,
                                   ['OP_NUM', 'newgeo', 'recombine', 'numtimes', 'min_read', 'numdays', 'pk_Dates',
                                    'min_Date']]
                gdf_pass_pks = pd.merge(gdf_tog, combined_reduced, on=['OP_NUM']).copy()
                gdf_pass_pks['verified'] = gdf_pass_pks.apply(lambda y: (True if y.numtimes > 1 else False),
                                                              axis=1).copy()
            if data_overlap.size == 0:
                gdf_pass_pks = gdf_bind_pks.copy()
                gdf_pass_pks['min_read'] = gdf_pass_pks.loc[:, 'OP_NUM']
                gdf_pass_pks['numtimes'] = 1
                gdf_pass_pks['numdays'] = 1

                gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:, 'geometry']
                gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:, ['OP_NUM']].to_numpy())].copy()
                # gdf_pass_pks['dates'] = gdf_pass_pks.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                gdf_pass_pks['dates'] = gdf_pass_pks.apply(
                    lambda x: datetime.fromtimestamp(int(x.OP_NUM[len(xCar)+1:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'), axis=1)

                gdf_pass_pks['pk_Dates'] = [list(x) for x in list(gdf_pass_pks.loc[:, ['dates']].to_numpy())]
                gdf_pass_pks['min_Date'] = gdf_pass_pks.loc[:, 'dates']
                gdf_pass_pks = gdf_pass_pks.drop(columns=['dates'])

                gdf_pass_pks['verified'] = False
                #           gdf_pass_pks['oldgeo'] = gdf_pass_pks.loc[:,'geometry']
                gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:, "newgeo"]
                together = pd.merge(gdf_pass_pks, gdf_tog,
                                    on=['OP_NUM', 'pk_LON', 'pk_LAT', 'pk_maxCH4_AB', 'geometry'])
                together['pass'] = whichpass
                gdf_pass_pks = together.copy()

        if gdf_bind_pks.shape[0] == 1:
            gdf_pass_pks = gdf_bind_pks.copy()
            gdf_pass_pks['min_read'] = gdf_pass_pks.loc[:, 'OP_NUM']
            gdf_pass_pks['numtimes'] = 1
            gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:, 'geometry']

            gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:, ['OP_NUM']].to_numpy())].copy()
            # gdf_pass_pks['dates'] = gdf_pass_pks.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
            gdf_pass_pks['dates'] = gdf_pass_pks.apply(
                lambda x: datetime.fromtimestamp(int(x.OP_NUM[len(xCar)+1:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'), axis=1)

            gdf_pass_pks['pk_Dates'] = [list(x) for x in list(gdf_pass_pks.loc[:, ['dates']].to_numpy())]
            gdf_pass_pks['min_Date'] = gdf_pass_pks.loc[:, 'dates']
            gdf_pass_pks['numdays'] = 1
            gdf_pass_pks = gdf_pass_pks.drop(columns=['dates'])

            gdf_pass_pks['verified'] = False
            epdat = pass_info.loc[:, ['OP_NUM', 'OP_EPOCHSTART']]
            gdf_pass_pks = pd.merge(gdf_pass_pks, epdat, on=['OP_NUM']).copy()
            data_overlap = pd.DataFrame(columns=['what', 'oh'])

        ### gdf_pass_pks
        #    Index(['OP_NUM', 'pk_LON', 'pk_LAT', 'pk_maxCH4_AB', 'geometry',
        #       'OP_EPOCHSTART', 'OB_EPOCH', 'OB_DATETIME', 'OB_CH4', 'OB_LON',
        #       'OB_LAT', 'OB_CH4_BASELINE', 'OB_CH4_THRESHOLD', 'OP_PEAK_DIST_M',
        #       'OP_PEAK_CH4', 'OB_TCH4', 'OB_PERIOD5MIN', 'OB_CH4_AB', 'newgeo',
        #       'recombine', 'numtimes', 'min_read', 'numdays', 'pk_Dates', 'min_Date',
        #       'verified'],
        #      dtype='object')

        gdf_pass_pks['pkGEO'] = gdf_pass_pks.loc[:, "geometry"]
        gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:, "newgeo"]
        del (gdf_pass_pks['newgeo'])
        gdf_pass_pks['pass'] = whichpass
        gdf_pass_pks['Overall_LON'] = gdf_pass_pks['pk_LON']
        gdf_pass_pks['Overall_LAT'] = gdf_pass_pks['pk_LAT']
        combinedOP1 = gdf_pass_pks.drop(columns=['recombine', 'pk_Dates']).drop_duplicates()

        # gdf_tot = pd.merge(gdf_pass_pks,datFram_wtLocMax.loc[:,['PEAK_NUM','pk_LON','pk_LAT']],on = ['PEAK_NUM','pk_LON','pk_LAT']).copy()
        ## condense by peak_num
        # gdfcop = gdf_tot.loc[:,['PEAK_NUM','geometry','min_read','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()

        #### WANT A DATAFRAME WITH
        # EACH OP SUMMARY
        # COMBINED WITH THE COMBINED SUMMARY
        # gdf_op_unique = gdf_pass_pks.loc[:,['numtimes','min_read','numdays','min_Date','verified','pass','Overall_LON','Overall_LAT']].drop_duplicates()

        if data_overlap.size != 0:
            gdf_op_unique = gdf_pass_pks.loc[:,
                            ['numtimes', 'min_read', 'numdays', 'min_Date', 'verified', 'pass', 'OB_LON',
                             'OB_LAT']].drop_duplicates()
            gdfcop = gdf_pass_pks.loc[:,
                     ['OP_NUM', 'min_read', 'min_Date', 'numtimes', 'verified', 'pass', 'pk_LAT', 'pk_LON',
                      'pk_maxCH4_AB']].drop_duplicates()
            combinedOP = weightedLoc(gdfcop, 'pk_LAT', 'pk_LON', 'min_read', 'pk_maxCH4_AB').loc[:, :].rename(
                columns={'pk_LAT': 'Overall_LAT', 'pk_LON': 'Overall_LON'}).reset_index(drop=True)
            combinedOP1 = pd.merge(combinedOP, gdfcop, on=['min_read'])

        if data_overlap.size == 0 and gdf_bind_pks.shape[0] != 1:
            gdf_op_unique = gdf_pass_pks.loc[:,
                            ['numtimes', 'min_read', 'numdays', 'min_Date', 'verified', 'pass', 'OB_LON',
                             'OB_LAT']].drop_duplicates()
            gdfcop = gdf_pass_pks.loc[:,
                     ['OP_NUM', 'min_read', 'min_Date', 'numtimes', 'verified', 'pass', 'pk_LAT', 'pk_LON',
                      'pk_maxCH4_AB']].drop_duplicates()
            combinedOP = weightedLoc(gdfcop, 'pk_LAT', 'pk_LON', 'min_read', 'pk_maxCH4_AB').loc[:, :].rename(
                columns={'pk_LAT': 'Overall_LAT', 'pk_LON': 'Overall_LON'}).reset_index(drop=True)
            combinedOP1 = pd.merge(combinedOP, gdfcop, on=['min_read'])

        ## TO FINDED WEIGHTED LOCATION OF EACH PK GROUP
        #    gdfcop = gdf_pass_pks.loc[:,['OP_NUM','min_read','min_Date','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
        #    combinedOP = weightedLoc(gdfcop,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').loc[:,:].rename(columns = {'pk_LAT':'Overall_LAT','pk_LON':'Overall_LON'}).reset_index(drop=True)
        #    combinedOP1 = pd.merge(combinedOP,gdfcop,on=['min_read'])

        ## getting the rest of the stuff
        # gdf_justloc = gdfcop.loc[:,['min_read','pk_LAT','pk_LON','min_Date']].reset_index(drop=True)

        # other = gdfcop.loc[:,['min_read','numtimes','min_Date']].reset_index(drop=True)

        geometry_temp = [Point(xy) for xy in zip(combinedOP1['Overall_LON'], combinedOP1['Overall_LAT'])]
        crs = {'init': 'epsg:4326'}
        gdf_OP = gpd.GeoDataFrame(combinedOP1, crs=crs, geometry=geometry_temp)
        gdf_OP = gdf_OP.to_crs(epsg=32610).copy()

        gdf_OP_reduced = gdf_OP.loc[:, ['min_read', 'geometry', 'numtimes', 'Overall_LON', 'Overall_LAT', 'min_Date',
                                        'verified']].drop_duplicates().reset_index(drop=True)
        gdf_OP_reduced.to_file(new_loc_json, driver="GeoJSON")

        # gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']),gdf_pass_pks.drop(columns=['geometry','oldgeo']),on=['min_read','min_Date','numtimes','pass','verified','pk_LAT','pk_LON','OP_NUM','pk_maxCH4_AB'])
        gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']), gdf_pass_pks.drop(columns=['geometry']),
                                     on=['min_read', 'min_Date', 'numtimes', 'pass', 'verified', 'pk_LAT', 'pk_LON',
                                         'OP_NUM', 'pk_maxCH4_AB'])

        # gdf_OP.to_csv(new_loc,index=False)
        gdf_OP_wrecombine.to_csv(new_loc, index=False)

        # geometry is the point of the lat/lon
        # gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)
        unique_peaks = gdf_pass_pks.loc[:, ['OP_NUM', 'pk_LAT', 'pk_LON', 'min_read', 'min_Date']].drop_duplicates()
        unique_peaks['save'] = True
        good_pks = list(unique_peaks.index)

        def getthing(index):
            if index in good_pks:
                return True
            else:
                return False

        gdf_pass_pks['wooind'] = gdf_pass_pks.index
        gdf_pass_pks['save'] = gdf_pass_pks.apply(lambda x: getthing(x.wooind), axis=1)

        unique_pks_tog = gdf_pass_pks.loc[gdf_pass_pks.save == True, :].reset_index(drop=True)
        unique_pks_tog['Latitude'] = unique_pks_tog.loc[:, 'pk_LAT']
        unique_pks_tog['Longitude'] = unique_pks_tog.loc[:, 'pk_LON']
        unique_pks_tog.to_csv(new_loc, index=False)

        return ()


########################################################################
#### sumData2
# Input: the two data groups to combine (already have gone through
#           'filterpeak')
# Output: One dataframe with the combination (combine peaks)

def sumData2(mainDF):
    from numpy import log
    import pandas as pd
    todo = mainDF.loc[:, ['OP_NUM', 'min_Date', 'pk_LON', 'pk_LAT', 'pk_maxCH4_AB', 'numtimes',
                          'min_read','OP_DISTANCE']].drop_duplicates().reset_index(drop=True)
    todo['logCH4'] = todo.apply(lambda y: log(y.pk_maxCH4_AB), axis=1)
    mnVals = todo.groupby('min_read', as_index=False).logCH4.mean().rename(columns={'logCH4': 'mnlogCH4'}).loc[:,
             ['min_read', 'mnlogCH4']]
    opMin = todo.groupby('min_read', as_index=False).OP_DISTANCE.min().rename(columns={'OP_DISTANCE': 'minDist'}).loc[:,
             ['min_read', 'minDist']]
    opMax = todo.groupby('min_read', as_index=False).OP_DISTANCE.max().rename(columns={'OP_DISTANCE': 'maxDist'}).loc[:,
             ['min_read', 'maxDist']]

    verLoc = weightedLoc(todo, 'pk_LAT', 'pk_LON', 'min_read', 'pk_maxCH4_AB').rename(
        columns={'pk_LAT': 'overallLAT', 'pk_LON': 'overallLON'}).reset_index(drop=True)
    together = pd.merge(verLoc, mnVals, on=['min_read'])
    final = pd.merge(together, mainDF, on=['min_read'])
    final = pd.merge(final,opMin,on=['min_read'])
    final = pd.merge(final,opMax,on = ['min_read'])
    return (final)


########################################################################
#### PASS COMBINE
# Input: the two data groups to combine (already have gone through
#           'filterpeak')
# Output: One dataframe with the combination (combine peaks)


def passCombine(firstgroup, secondgroup, xCar, buffer='30'):
    import ast
    def strList(x):
        x = ast.literal_eval(x)
        x = [n.strip() for n in x]
        return (x)

    import pandas as pd  #
    import geopandas as gpd
    from shapely.geometry import Point
    buffer = float(buffer)

    seclist = secondgroup.columns

    if not "prev_read" in seclist:
        sgrp = secondgroup.copy()
        secondgroup['prev_read'] = sgrp.loc[:, 'min_read']
    ### MAKE GEOPANDAS FRAME OF THE FIRST GROUP PEAKS
    first_geo = [Point(xy) for xy in zip(firstgroup['pk_LON'], firstgroup['pk_LAT'])]
    #   crs = {'init': 'epsg:4326'}
    crs = {'init': 'epsg:32610'}

    if 'geometry' in firstgroup.columns:
        firstgrp = gpd.GeoDataFrame(firstgroup.drop(columns=['geometry', 'pkGEO']), crs=crs, geometry=first_geo)
    else:
        firstgrp = gpd.GeoDataFrame(firstgroup.drop(columns=['pkGEO']), crs=crs, geometry=first_geo)

    first_buffg = firstgrp.copy()
    first_buff = first_buffg.copy().drop(columns=['geometry'])
    # first_buff['geometry2'] = first_buff.loc[:,'geometry'].buffer(30)

    # first_buff['geometry2'] = first_buff.apply(lambda x: x.geometry.buffer(0.0001*3),axis =1)

    first_buff['geometry'] = first_buffg.apply(lambda x: x.geometry.buffer(0.00001 * buffer), axis=1)
    # first_buff['geometry'] = first_buffg.apply(lambda x: x.geometry.buffer(0.0001*3),axis =1)

    # first_buff = first_buff.drop(columns = ['geometry'])
    # first_buff['geometry'] = first_buff.loc[:,'geometry2']
    # first_buff = first_buff.drop(columns = ['geometry2'])

    # first_buff = first_buff.to_crs(epsg=32610)
    # first_buff.plot()

    firstgrp = first_buff.copy()

    sec_geo = [Point(xy) for xy in zip(secondgroup['pk_LON'], secondgroup['pk_LAT'])]
    secgrp = gpd.GeoDataFrame(secondgroup.drop(columns=['geometry', 'pkGEO']), crs=crs, geometry=sec_geo)

    sec_buffg = secgrp.copy()
    sec_buff = sec_buffg.copy().drop(columns=['geometry'])
    sec_buff['geometry'] = sec_buffg.apply(lambda x: x.geometry.buffer(0.00001 * buffer), axis=1)
    secgrp = sec_buff.copy()

    ### GATHER THE INDIVIDUAL PEAKS, WITH THEIR LOCATIONS
    # first_pks = firstgrp.loc[:,['min_read','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    # sec_pks = secondgrp.loc[:,['PEAK_NUM','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    # tot_pks = pd.concat([first_pks,sec_pks])

    ## choosing the unique OPs from each group
    first_pks = firstgrp.loc[:, ['OP_NUM', 'min_read', 'pk_LAT', 'pk_LON', 'pk_maxCH4_AB']].drop_duplicates()
    sec_pks = secgrp.loc[:, ['OP_NUM', 'min_read', 'pk_LAT', 'pk_LON', 'pk_maxCH4_AB']].drop_duplicates()

    # combining to make one large dataframe of unique OPs
    tot_pks = pd.concat([first_pks, sec_pks])

    ### COMBINE EACH GROUP (GEOMETRICALLY) BY THEIR OVERALL NAME (COULD BE VERIFIED)
    first_dis = firstgrp.dissolve(by='min_read', as_index=False)[
        ['min_read', 'geometry', 'recombine', 'verified', 'pass']].copy()
    sec_dis = secgrp.dissolve(by='min_read', as_index=False)[
        ['min_read', 'geometry', 'recombine', 'verified', 'pass']].copy()

    gdf_bind_pks = pd.concat([first_dis, sec_dis]).loc[:, ['min_read', 'geometry', 'recombine']]
    # gdf_tog = pd.concat([firstgrp.drop(['pk_LAT', 'pk_LON','pk_maxCH4_AB'], axis=1),secondgrp.drop(['pk_LAT', 'pk_LON','pk_maxCH4_AB'], axis=1)]).copy()
    gdf_tog = pd.concat([firstgroup.drop(['pk_LAT', 'pk_LON', 'pk_maxCH4_AB'], axis=1),
                         secondgroup.drop(['pk_LAT', 'pk_LON', 'pk_maxCH4_AB'], axis=1)]).copy()
    gdf_tog2 = pd.concat([firstgroup, secondgroup]).copy()

    gdf_bind_pks['prev_read'] = gdf_bind_pks.loc[:, 'min_read']
    gdf_tog['prev_read'] = gdf_tog.loc[:, 'min_read']
    if gdf_bind_pks.shape[0] > 1:
        data_overlap = gpd.GeoDataFrame(crs=gdf_bind_pks.crs).copy()
        data_temp = gdf_bind_pks.copy()
        data_temp = data_temp.to_crs(epsg=32610)

        for index, row in data_temp.iterrows():
            data_temp1 = data_temp.loc[data_temp.min_read != row.min_read,]
            data_temp1 = data_temp1.to_crs(epsg=32610)

            # woo = data_temp1
            # what = row

            # check if intersection occured
            overlaps = data_temp1[data_temp1.geometry.overlaps(row.geometry)]['min_read'].tolist()
            if len(overlaps) > 0:
                # compare the area with threshold
                for y in overlaps:
                    temp_area = gpd.overlay(data_temp.loc[data_temp.min_read == y,],
                                            data_temp.loc[data_temp.min_read == row.min_read,], how='intersection')
                    temp_area = temp_area.to_crs(epsg=32610)
                    temp_area = temp_area.loc[temp_area.geometry.area >= 0]
                    # temp_union = gpd.overlay(data_temp.loc[data_temp.PEAK_NUM==y,],data_temp.loc[data_temp.PEAK_NUM==row.PEAK_NUM,],how='union')
                    if temp_area.shape[0] > 0:
                        temp_union = gpd.overlay(data_temp.loc[data_temp.min_read == y,],
                                                 data_temp.loc[data_temp.min_read == row.min_read,], how='union')
                        data_overlap = gpd.GeoDataFrame(pd.concat([temp_union, data_overlap], ignore_index=True),
                                                        crs=data_temp.crs)

        if data_overlap.size != 0:
            firstnull = data_overlap[data_overlap.min_read_1.isnull()].copy()
            firstnull.loc[:, 'min_read_1'] = firstnull.loc[:, 'min_read_2']

            secnull = data_overlap[data_overlap.min_read_2.isnull()].copy()
            secnull['min_read_2'] = secnull['min_read_1'].copy()

            withoutNA = data_overlap.dropna().copy()
            allTog = pd.concat([firstnull, secnull, withoutNA]).reset_index().copy()

            over = allTog.copy().drop(columns=['index'])
            over['sorted'] = over.apply(lambda y: sorted([y['min_read_1'], y['min_read_2']]), axis=1).copy()
            over['sorted'] = over.sorted.apply(lambda y: ','.join(y)).copy()

            over['same'] = over.apply(lambda x: x.min_read_1 != x.min_read_2, axis=1)
            over = over.loc[over.same == True, :].drop(columns=['same'])
            over = over.copy().drop_duplicates('sorted').reset_index(drop=True)

            def checkLst(opList):
                if isinstance(opList, str):
                    opList = strList(opList)
                return (opList)

            # over['bothcombine'] = over.apply(lambda x: sorted(strList(x.recombine_1)+ strList(x.recombine_2)),axis=1)

            over['bothcombine'] = over.apply(lambda x: sorted(checkLst(x.recombine_1) + checkLst(x.recombine_2)),
                                             axis=1)

            over['combined'] = [list(x) for x in list(over.loc[:, ['min_read_1', 'min_read_2']].to_numpy())].copy()
            over['VER_NUM'] = over.apply(lambda y: y.combined, axis=1).copy()
            over['min_val'] = over.apply(lambda y: min(y.combined), axis=1).copy()
            over = over.reset_index().loc[:,
                   ['min_read_1', 'min_read_2', 'geometry', 'combined', 'min_val', 'bothcombine']]

            new_1 = over.copy().drop(columns=['min_read_2']).rename(columns={'min_read_1': 'min_read'})
            new_2 = over.copy().drop(columns=['min_read_1']).rename(columns={'min_read_2': 'min_read'})
            newtog = pd.concat([new_1, new_2])
            minreads = newtog.min_read.unique().tolist()

            if 'prev_read' in gdf_tog2.columns:
                toChange = gdf_tog2[gdf_tog2['min_read'].isin(minreads)].drop(columns=['geometry', 'prev_read'])
            elif 'prev_read' not in gdf_tog2.columns:
                toChange = gdf_tog2[gdf_tog2['min_read'].isin(minreads)].drop(columns=['geometry'])

            toChangecombined = pd.merge(toChange, newtog, on=['min_read']).drop(
                columns=['geometry', 'numtimes', 'verified', 'recombine'])
            toChangecombined = toChangecombined.rename(columns={'min_read': 'prev_read', 'bothcombine': 'recombine'})
            # toChangecombined['numtimes'] = toChangecombined.apply(lambda x: len(x.recombine),axis = 1)
            toChangecombined['numtimes'] = toChangecombined.apply(lambda x: countTimes(x.recombine,xCar), axis=1)

            toChangecombined['verified'] = toChangecombined.apply(lambda x: x.numtimes > 1, axis=1)
            toChangecombined = toChangecombined.rename(columns={'min_val': 'min_read'}).drop(columns=['combined'])

            toNotChange = gdf_tog2[~gdf_tog2['min_read'].isin(minreads)].drop(columns=['geometry'])
            if 'prev_read' in toNotChange.columns:
                toNotChange = toNotChange.drop(columns=['prev_read'])
                toNotChange['prev_read'] = toNotChange.min_read
                # toNotChange.loc[:,'prev_read'] = gdf_tog2[~gdf_tog2['min_read'].isin(minreads)].min_read
            elif 'prev_read' not in toNotChange.columns:
                toNotChange.loc[:, 'prev_read'] = gdf_tog2[~gdf_tog2['min_read'].isin(minreads)].min_read

            newCombined = pd.concat([toChangecombined, toNotChange])
        elif data_overlap.size == 0:
            newCombined = gdf_tog2.copy()
            if 'prev_read' not in newCombined.columns:
                newCombined.loc[:, 'prev_read'] = gdf_tog2.loc[:, 'min_read']
    elif gdf_bind_pks.shape[0] == 1:
        newCombined = gdf_tog2.copy()
        if 'prev_read' not in newCombined.columns:
            newCombined.loc[:, 'prev_read'] = gdf_tog2.loc[:, 'min_read']
    return (newCombined)
    # return(gdf_tot_pks)
