

# =============================================================================
#     first_geo = [Point(xy) for xy in zip(firstgroup['pk_LON'], firstgroup['pk_LAT'])]
#     crs = {'init': 'epsg:32610'}
#     firstgrp = gpd.GeoDataFrame(firstgroup.drop(columns=['geometry','pkGEO']),crs = crs,geometry = first_geo)
#
#     first_buff = firstgrp.copy()
#     first_buff = first_buff.to_crs(epsg=32610)
#     #first_buff['geometry2'] = first_buff.loc[:,'geometry'].buffer(30)
#
#     first_buff['geometry2'] = first_buff.apply(lambda x: x.geometry.buffer(.030),axis =1)
#
#     first_buff = first_buff.drop(columns = ['geometry'])
#     first_buff['geometry'] = first_buff.loc[:,'geometry2']
#     first_buff = first_buff.drop(columns = ['geometry2'])
#
#
#     firstgrp = first_buff.copy()
# #    sec_geo = secondgroup.copy().pkGEO
# #    crs = {'init': 'epsg:4326'}
# #    secondgrp = gpd.GeoDataFrame(secondgroup.copy(),crs = crs,geometry = sec_geo)
# #
#     sec_geo = [Point(xy) for xy in zip(secondgroup['pk_LON'], secondgroup['pk_LAT'])]
#     crs = {'init': 'epsg:32610'}
#     secgrp = gpd.GeoDataFrame(secondgroup.drop(columns=['geometry','pkGEO']),crs = crs,geometry = sec_geo)
#
#     sec_buff = secgrp.copy()
#     sec_buff = sec_buff.to_crs(epsg=32610)
#     sec_buff['geometry'] = sec_buff.loc[:,'geometry'].buffer(30)
#
#     secgrp = sec_buff.copy()
# =============================================================================




# =============================================================================
#
#             over['recombine'] =
#
#             over.apply(lambda x: (x.bothcombine),axis=1)
#
#             overcop = over.copy()
#             overcop['recombine'] = overcop.combined.copy()
#
#
#
#             for index,row in overcop.iterrows():
#                 #rowwoo = row
#
#                 ## fixthis
#
#                 first_thing = first_dis[first_dis['min_read']== row.min_read_1].loc[:,['recombine']]
#                 #first_thing = first_dis.loc[first_dis.min_read == row.min_read_1,'recombine']
#                 #first_thing = first_dis.loc[:,first_dis.min_read == row.min_read_1].loc[:,'recombine']
#
#                 #first_thing = firstpass[firstpass['min_read']== row.min_read_1].loc[:,['recombine']]
#                 firstcomb = first_thing.recombine.explode().copy()
#                 first_list = firstcomb.reset_index().recombine.copy()
#                 #first_list = firstpass.loc[:,firstpass['min_read']== row.min_read_1].recombine
#                 #first_list.explode()
#                 #first_list = firstpass.loc[:,firstpass['min_read']== row.min_read_1].recombine.explode()[0]
#                # second_list = secondpass[secondpass['min_read']== row.min_read_2].recombine
#
#
#                ## fix this
#
#                 sec_thing = sec_dis[sec_dis['min_read']== row.min_read_1].loc[:,['recombine']]
#
#                 #sec_thing = secondpass[secondpass['min_read']== row.min_read_2].loc[:,['recombine']]
#                 seccomb = sec_thing.recombine.explode()
#                 sec_list = seccomb.reset_index().recombine
#
#                 firstdf = pd.DataFrame(first_list)
#                 secdf = pd.DataFrame(sec_list)
#
#
#                 tot_df = pd.concat([firstdf,secdf])
#                 tot_list = tot_df.recombine.tolist()
#
#                 overcop.at[index, 'recombine']  = tot_list.copy()
#
#
#
#
#     ## this recombines the lists together to have the combined entries together?
#
#             overcop['recombine']= overcop.apply(lambda y: sorted(y.recombine),axis=1).copy()
#             overcop['min_read'] = overcop.apply(lambda y: min(y.recombine),axis=1).copy()
#             newOverlap = overcop.dissolve(by='min_read',as_index=False).loc[:,['min_read','geometry','recombine']].copy()
#
#
#             combined = gdf_bind_pks.copy()
#             #combined['recombine'] = [list(x) for x in list(combined.loc[:,['min_read']].to_numpy())]
#             combined['numtimes'] = 1
#             combined['newgeo'] = combined.copy().geometry
#             combined['oldgeo'] = combined.copy().geometry
#
#             #combined['min_read'] = combined.min_rea
#             combined = combined.reset_index()
#             for index,row in combined.iterrows():
#                 for index2,row2 in newOverlap.iterrows():
#                     if row.min_read in row2.recombine:
#                         combined.at[index, 'recombine']  = row2.recombine.copy()
#                         combined.at[index, 'min_read']  = row2.copy().min_read
#
#             combined['numtimes'] = combined.copy().apply(lambda y: len(y.recombine),axis = 1).copy()
#             combined['geometry'] = combined.copy().newgeo
#
#             del(combined['newgeo'])
#             combined_reduced = combined.loc[:,['min_read','geometry','oldgeo','recombine','numtimes','prev_read']]
#             #gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['min_read'])
#             gdf_tog = gdf_tog.drop('min_read',axis=1)
#             gdf_tog = gdf_tog.drop('numtimes',axis=1)
#             gdf_tog = gdf_tog.drop('recombine',axis=1)
#
#             #gdf_tog = gdf_tog.drop(columns = ['min_read','numtimes','recombine'])
#
#             # gdf_tog = gdf_tog.drop('newgeo',axis=1)
#            # gdf_tog['firstgeo'] = gdf_tog.copy().oldgeo
#            # gdf_tog['secondgeo'] = gdf_tog.copy().geometry
#
#            # del(gdf_tog['geometry'])
#             #gdf_tog.drop(columns=['geometry'])
#            # del(gdf_tog['oldgeo'])
#
#
#             gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['prev_read']).copy()
#
#             gdf_pass_pks['verified'] = gdf_pass_pks.copy().apply(lambda y: (True if y.numtimes > 1 else False),axis=1 )
#
# =============================================================================
# =============================================================================
#         if data_overlap.size == 0:
#             gdf_pass_pks = gdf_tog.copy()
#
#     ## didnt adress if the bind shape was only size only 1
#         gdf_tot_pks = pd.merge(gdf_pass_pks,tot_pks,on = ['OP_NUM','min_read']).copy()
#         gdf_tot_pks['numtimes']= gdf_tot_pks.apply(lambda x: len(x.recombine.split(',')),axis=1)
#     #return(gdf_pass_pks)
# =============================================================================


########################################################################
#### IDENTIFY PEAKS
# Input: a .csv file with processed data (already have gone through 'processRawDataEng')
# Output: saves many files, but finds elevated readings

# =============================================================================
# def IdentifyPeaks(xCar, xDate, xDir, xFilename,outDir,processedFileLoc,Engineering,threshold = '.1',xTimeThreshold = '5.0',minElevated = '2',xB = '1020',basePerc = '50'):
#     import csv, numpy
#     import geopandas as gpd
#     import shutil
#     try:
#         baseCalc = float(basePerc)
#         xABThreshold = float(threshold)
#         minElevated = float(minElevated)
#         #xABThreshold = 0.1                 # above baseline threshold above the mean value
#         xDistThreshold = 160.0                 # find the maximum CH4 reading of observations within street segments of this grouping distance in meters
#         xSDF = 4                    # multiplier times standard deviation for floating baseline added to mean
#         #xB = 1020       # the number of records that constitutes the floating baseline time -- 7200 = 1 hour (assuming average of 0.5 seconds per record)
#         #xB = 102 # since it is 1 record/second
#
#         xB = int(xB)
#         #xB = 300 #five min?
#         xTimeThreshold = float(xTimeThreshold)
#
#         fn = xDir + "/" + xFilename      #set raw text file to read in
#         fnOut = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".csv"       #set CSV format output for observed peaks for a given car, day, city
#         fnShape = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".shp"
#         fnLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".log"       #set CSV output for observed peaks for a given car, day, city
#         pkLog = outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + "_info.csv"       #set CSV output for observed peaks for a given car, day, city
#
#         jsonOut =  outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + ".json"       #set CSV format output for observed peaks for a given car, day, city
#
#         infOut = processedFileLoc + xCar + "_" + xDate.replace("-","") + "_info.csv"
#         print(str(outDir + "Peaks" + "_" + xCar + "_" + xDate.replace("-","") + "_info.csv"))
#
#         fLog = open(fnLog, 'w')
#         shutil.copy(infOut,pkLog)
#
#
#         #field column indices for various variables
#         if Engineering == True:
#             fDate = 0; fTime = 1; fEpochTime = 2;
#             fNanoSeconds = 3; fVelocity = 4;  fU = 5;
#             fV = 6; fW = 7;fBCH4 = 10;#fBCH4 = 8;
#             #fBRSSI = 9;
#             fTCH4 = 10; TRSSI = 11;PRESS = 12;
#             INLET = 13; TEMP = 14;CH4 = 15;
#             H20 = 16;C2H6 = 17;R = 18;
#             C2C1 = 19; BATT = 20;POWER = 21;
#             CURR = 22;SOCPER = 23;fLat = 24;
#             fLon = 25;
#         elif not Engineering:
#             fDate = 0; fTime = 1; fEpochTime = 2;
#             fNanoSeconds = 3; fVelocity = 4;  fU = 5;
#             fV = 6; fW = 7;fBCH4 = 8;
#             fBRSSI = 9;
#             fTCH4 = 10; TRSSI = 11;PRESS = 12;
#             INLET = 13; TEMP = 14;CH4 = 15;
#             H20 = 16;C2H6 = 17;R = 18;
#             C2C1 = 19; BATT = 20;POWER = 21;
#             CURR = 22;SOCPER = 23;fLat = 24;
#             fLon = 25;
#
#         #read data in from text file and extract desired fields into a list, padding with 5 minute and hourly average
#         x1 = []; x2 = []; x3 = []; x4 = []; x5 = []; x6 = []; x7 = []; x8 = []
#
#         count = -1
#         with open(fn, 'r') as f:
#             t = csv.reader(f)
#             for row in t:
#                 if count < 0:
#                     count += 1
#                     continue
#
#                 datet= row[fDate].replace("-","")+row[fTime].replace(":","")
#
#                 #x1.append(float(epoch));
#                 x1.append(float(str(row[fEpochTime]) + '.' + str(row[fNanoSeconds])));
#                 x2.append(float(int(datet)));
#                 x3.append(float(row[fLat]));
#                 x4.append(float(row[fLon]));
#                 x5.append(float(row[fBCH4]));
#                 x6.append(float(row[fTCH4]))
#                 x7.append(0.0);
#                 x8.append(0.0)
#
#                 #print (str(row[fLat])+ str(row[1]))
#                 count += 1
#         print ("Number of observations processed: " + str(count))
#
#         #convert lists to numpy arrays
#         aEpochTime = numpy.array(x1); aDateTime = numpy.array(x2); aLat = numpy.array(x3); aLon = numpy.array(x4); aCH4 = numpy.array(x5); aTCH4 = numpy.array(x6)
#         aMean = numpy.array(x7); aThreshold = numpy.array(x8)
#
#         xLatMean = numpy.mean(aLat)
#         xLonMean = numpy.mean(aLon)
#
#         fLog.write ( "Day CH4_mean = " + str(numpy.mean(aCH4)) + ", Day CH4_SD = " + str(numpy.std(aCH4)) + "\n")
#         fLog.write( "Center lon/lat = " + str(xLonMean) + ", " + str(xLatMean) + "\n")
#         #pkLog.write('hi')
#         lstCH4_AB = []
#
#         #generate list of the index for observations that were above the threshold
#         for i in range(0,count-2):
#             if ((count-2)>xB):
#                 topBound = min((i+xB), (count-2))
#                 botBound = max((i-xB), 0)
#
#                 for t in range(min((i+xB), (count-2)), i, -1):
#                     if aEpochTime[t] < (aEpochTime[i]+(xB/2)):
#                         topBound = t
#                         break
#                 for b in range(max((i-xB), 0), i):
#                     if aEpochTime[b] > (aEpochTime[i]-(xB/2)):
#                         botBound = b
#                         break
#
#                 xCH4Mean = numpy.percentile(aCH4[botBound:topBound],baseCalc)
#                # xCH4SD = numpy.std(aCH4[botBound:topBound])
#             else:
#                 xCH4Mean = numpy.percentile(aCH4[0:(count-2)],baseCalc)
#                 #xCH4SD = numpy.std(aCH4[0:(count-2)])
#             xThreshold = xCH4Mean + (xCH4Mean * xABThreshold)
#
#             if (aCH4[i] > xThreshold):
#                 lstCH4_AB.append(i)
#                 aMean[i] = xCH4Mean    #insert mean + SD as upper quartile CH4 value into the array to later retreive into the peak calculation
#                 aThreshold[i] = xThreshold
#
#         # now group the above baseline threshold observations into groups based on distance threshold
#         lstCH4_ABP = []
#         xDistPeak = 0.0
#         xCH4Peak = 0.0
#         xTime = 0.0
#         cntPeak = 0
#         cnt = 0
#         sID = ""
#         sPeriod5Min = ""
#         prevIndex = 0
#         for i in lstCH4_AB:
#             if (cnt == 0):
#                 xLon1 = aLon[i]; xLat1 = aLat[i]
#             else:
#                 # calculate distance between points
#                 xDist = haversine(xLat1, xLon1, aLat[i], aLon[i])
#                 xDistPeak += xDist
#                 xCH4Peak += (xDist * (aCH4[i] - aMean[i]))
#                 xLon1 = aLon[i]; xLat1 = aLat[i]
#                 if (sID == ""):
#                     xTime = aEpochTime[i]
#                     sID = str(xCar) + "_" + str(xTime)
#                     sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1))) # 30 sec
#                 if ((aEpochTime[i]-aEpochTime[prevIndex]) > xTimeThreshold):       #initial start of a observed peak
#                     cntPeak += 1
#                     xTime = aEpochTime[i]
#                     xDistPeak = 0.0
#                     xCH4Peak = 0.0
#                     sID = str(xCar) + "_" + str(xTime)
#                     sPeriod5Min = str(int((aEpochTime[i] - 1350000000) / (30 * 1))) # 30 sec
#                     #print str(i) +", " + str(xDist) + "," + str(cntPeak) +"," + str(xDistPeak)
#                 lstCH4_ABP.append([sID, xTime, aEpochTime[i], aDateTime[i], aCH4[i], aLon[i], aLat[i], aMean[i] ,aThreshold[i], xDistPeak, xCH4Peak, aTCH4[i], sPeriod5Min])
#             cnt += 1
#             prevIndex = i
#
#         # Finding peak_id larger than 160.0 m
#         tmpsidlist = []
#         for r in lstCH4_ABP:
#             if (float(r[9])>160.0) and (r[0] not in tmpsidlist):
#                 tmpsidlist.append(r[0])
#         cntPeak-=len(tmpsidlist)
#
#         fLog.write ( "Number of peaks found: " + str(cntPeak) + "\n")
#         print (xCar + "\t" + xDate + "\t" + xFilename + "\t" + str(count) + "\t" + str(len(lstCH4_ABP)))
#         #### calculate attribute for the area under the curve -- PPM
#
#         #write out the observed peaks to a csv to be read into a GIS
#         fOut = open(fnOut, 'w')
#         #s = "PEAK_NUM,EPOCHSTART,EPOCH,DATETIME,CH4,LON,LAT,CH4_BASELINE,CH4_THRESHOLD,PEAK_DIST_M,PEAK_CH4,TCH4,PERIOD5MIN\n"
#         s = "OP_NUM,OP_EPOCHSTART,OB_EPOCH,OB_DATETIME,OB_CH4,OB_LON,OB_LAT,OB_CH4_BASELINE,OB_CH4_THRESHOLD,OP_PEAK_DIST_M,OP_PEAK_CH4,OB_TCH4,OB_PERIOD5MIN\n"
#
#         fOut.write(s)
#
#         truecount = 0
#         for r in lstCH4_ABP:
#             if r[0] not in tmpsidlist:
#                 s = ''
#                 for rr in r:
#                     s += str(rr) + ','
#                 s = s[:-1]
#                 s += '\n'
#                 fOut.write(s)
#                 truecount += 1
#         fOut.close()
#         fLog.close()
#         import pandas as pd
#         openFile = pd.read_csv(fnOut)
#         from shapely.geometry import Point
#         if openFile.shape[0] != 0:
#             tempCount = openFile.groupby('OP_NUM',as_index=False).OP_EPOCHSTART.count().rename(columns={'OP_EPOCHSTART':'Frequency'})
#             tempCount = tempCount.loc[tempCount.Frequency>=minElevated,:]
#             if tempCount.shape[0]==0:
#              print("No Observed Peaks with enough Elevated Readings Found in the file: " + str(xFilename) )
#             elif tempCount.shape[0]!=0:
#                 oFile = pd.merge(openFile,tempCount,on=['OP_NUM'])
#                 openFile = oFile.copy()
#                 del(oFile)
#                 openFile['minElevated'] = openFile.apply(lambda x: int(minElevated),axis=1)
#                 openFile.to_csv(fnOut,index=False)
#                 openFile['OB_CH4_AB'] = openFile.loc[:,'OB_CH4'].sub(openFile.loc[:,'OB_CH4_BASELINE'], axis = 0)
#
#                 fileWt = weightedLoc(openFile,'OB_LAT','OB_LON','OP_NUM','OB_CH4_AB').loc[:,:].rename(columns = {'OB_LAT':'pk_LAT','OB_LON':'pk_LON'}).reset_index(drop = True)
#                 geometry_temp = [Point(xy) for xy in zip(fileWt['pk_LON'], fileWt['pk_LAT'])]
#                 crs = {'init': 'epsg:4326'}
#
#                     #geometry is the point of the lat/lon
#                 #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)
#
#                 ## BUFFER AROUND EACH 'OP_NUM' OF 30 M
#                 gdf_buff = gpd.GeoDataFrame(fileWt, crs=crs, geometry=geometry_temp)
#                 #gdf_buff = makeGPD(datFram,'LON','LAT')
#                 gdf_buff = gdf_buff.to_crs(epsg=32610)
#                 gdf_buff['geometry'] = gdf_buff.loc[:,'geometry'].buffer(30)
#                 gdf_buff.to_file(jsonOut, driver="GeoJSON")
#         elif openFile.shape[0] == 0:
#             print("No Observed Peaks Found in the file: " + str(xFilename) )
#
#         if truecount > 0:
#             #arcpy.MakeXYEventLayer_management(fnOut,"LON","LAT",xCar +
#             #"L","GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',
#             #SPHEROID['WGS_1984',6378137.0,298.257223563]],
#             #PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];
#             #-400 -400 1000000000;-100000 10000;-100000 10000;
#             #8.98315284119522E-09;0.001;0.001;IsHighPrecision","#")
#             #arcpy.FeatureToPoint_management(xCar + "L",fnShape,"CENTROID")
#             #arcpy.Delete_management(xCar+"L")
#             return True
#     except ValueError:
#             print ("Error in Identify Peaks")
#             return False
#
# =============================================================================


#### extra

# =============================================================================
# def summarizeDat(totalData):
#     pkRed = totalData.loc[:,['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']]. \
#         drop_duplicates().reset_index().loc[:,['PEAK_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','numtimes','min_read']]
#     verLoc = weightedLoc(pkRed,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').rename(columns = {'pk_LAT':'overallLAT','pk_LON':'overallLON'}).reset_index(drop = True)
#     pkRed['logCH4'] = pkRed.apply(lambda y: log(y.pk_maxCH4_AB),axis = 1)
#     mnVals = pkRed.groupby('min_read',as_index=False).logCH4.mean().rename(columns ={'logCH4':'mnlogCH4'}).loc[:,['min_read','mnlogCH4']]
#     together = pd.merge(verLoc,mnVals,on = ['min_read'])
#     final = pd.merge(together,totalData,on=['min_read'])
#     return(final)
# =============================================================================

# =============================================================================
# ########################################################################
# #### IsInPk
# # helper function to find if in thing
# # Input: a list of peak times included in a given combined peak
# # Output: counts number of times the peak was seen (not in same 5 min period)
#
# def IsInPK (peakNumRow,listPksRow):
#     if peakNumRow.PEAK_NUM in listPksRow.combined:
#         return(True)
#     else:
#         return(False)
# =============================================================================


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 13:23:50 2020

@author: emilywilliams
"""

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

#### adding saving stuff

mainInfo.drop_duplicates().reset_index(drop=True).FILENAME.to_csv(final_info_loc)
mainThing.reset_index(drop=True).to_csv(final_main_csv_loc)

combined = summarize_data_2(mainThing) ## finds locations and mean log ch4 for each peak (either verified or non yet)

## combined so only with the same overall peak
uniquePk = combined.loc[:,['min_read']].drop_duplicates()
uniqueList = combined.loc[uniquePk.index,['min_read','recombine']]
uniqueOther = combined.loc[:,['min_read','overallLON','overallLAT','mnlogCH4',
                             'verified','numtimes','minDist','maxDist']].drop_duplicates()
allTog = pd.merge(make_GPD(uniqueOther,'overallLAT','overallLON'),uniqueList,on = ['min_read'])
allTog['em'] = allTog['mnlogCH4'].swifter.apply(lambda y: estimate_emissions(y))
allTog['threshold'] = allTog['em'].swifter.apply(lambda x: threshold)


##### SPLITTING IF THE PEAKS WERE VERIFIED OR NOT
verTog = allTog.loc[allTog.numtimes!= 1,:]

if verTog.size > 0:
    verTog.drop(columns=['recombine']).to_file(shp_file_loc, driver="GeoJSON")
    print(f'I found {len(verTog.min_read.unique())} verified peaks')
    vpNew = len(verTog.min_read.unique())
if verTog.size ==0:
    print("Sorry, no verified peaks were found.")
    vpNew = 0
if allTog.size> 0:
    allTog.drop(columns=['recombine']).to_file(op_shp_file_loc, driver="GeoJSON")
    allTog.to_csv(all_op_csv_loc)

if allTog.size == 0:
    print("Sorry, no observed peaks were found in the given data")


if not addingFiles:
    print(f"I processed {len(to_filter)} days of driving. I analysed the data using a threshold of {100 + float(threshold)*100}% for an elevated reading, \n \
    where the threshold was calculated using the {baseline_percentile}th percentile over {back_obs_num} observations. \n \
    I filtered the speed of the car to be between {min_car_speed}mph and {max_car_speed}mph.\n \
    I created 3 summary files located here:{final_results_dir}.\n \
    The processing took {round((time.time()-start)/60,3)} minutes. \n \
    I found {len(mainThing.min_read.unique())} observed peaks.")

elif addingFiles:
    print(f"I processed an additional {len(to_filter)} days of driving. I analysed the data using a threshold of {100 + float(threshold) * 100}% for an elevated reading, \n \
    where the threshold was calculated using the {baseline_percentile}th percentile over {back_obs_num} observations. \n \
    I filtered the speed of the car to be between {min_car_speed}mph and {max_car_speed}mph.\n \
    I created 3 summary files located here:{final_results_dir}.\n \
    The processing took {round((time.time() - start) / 60, 3)} minutes. \n \
    I found {len(mainThing.min_read.unique()) - curOP} additional observed peaks, and {vpNew - curVP} VPs.")

