### testring the processData aeris

theResult = ProcessRawData(xCar, xDate, rawDir, file, bFirst, 1, processedFileLoc, initialTimeIgnore, shift,
                           maxCarSpeed, minCarSpeed)

def ProcessRawDataAeris(xCar, xDate, xDir, xFilename, bFirst, gZIP, xOut, initialTimeBack, shift, maxSpeed='45', minSpeed='2'):

xCar = 'TrCar'
xDate = '20200630'
xDir = rawDir
xFilename = file
gZIP = 1
xOut = processedFileLoc
initialTimeBack = initialTimeIgnore
maxSpeed = maxCarSpeed
minSpeed = minCarSpeed
filtopDir = '/Users/emilywilliams/Documents/DrivingData/retestData/CoDrive_45_102_5pc_pycharm/FilteredObservedPeaks/'
index = 1


stime = time.time()
filterPeak_normal(xCar, xDate, opDir, file, filtopDir, whichpass=index)
etime = time.time()
print(etime - stime)

stime = time.time()
filterPeak_swift(xCar, xDate, opDir, file, filtopDir, whichpass=index)
etime = time.time()
print(etime - stime)

def filterPeak_normal(xCar,xDate,xDir,xFilename, outFolder,whichpass = 0):
    import pandas as pd #
    import geopandas as gpd
    import shutil
    from datetime import datetime
    from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry

    file_loc = xDir + xFilename
    new_loc = outFolder + "Filtered" + xFilename
    new_loc_json = new_loc[:-3] + 'json'

    oldInfo = xDir + 'Peaks_' + xCar + "_" + xDate.replace("-","") + "_info.csv"
    newInfo = outFolder + 'FilteredPeaks_' + xCar + "_" + xDate.replace("-","") + "_info.csv"

    shutil.copy(oldInfo,newInfo)

    # identified peaks has the columns:
    #Index(['OP_NUM', 'OP_EPOCHSTART', 'OB_EPOCH', 'OB_DATETIME', 'OB_CH4',
    #   'OB_LON', 'OB_LAT', 'OB_CH4_BASELINE', 'OB_CH4_THRESHOLD',
    #   'OP_PEAK_DIST_M', 'OP_PEAK_CH4', 'OB_TCH4', 'OB_PERIOD5MIN'],
    #  dtype='object')

    datFram = pd.read_csv(file_loc)

    if datFram.shape[0] == 0:
        print("Not filtering this file, no peak in it!")
    elif datFram.shape[0] == 1: ## only one thing to begin with
        datFram_cent =  datFram.loc[:,:]
        datFram_cent['OB_CH4_AB'] = datFram.loc[:,'OB_CH4'].sub(datFram.loc[:,'OB_CH4_BASELINE'], axis = 0)
        maxch4 = datFram_cent.groupby('OP_NUM',as_index = False).OB_CH4_AB.max().rename(columns = {'OB_CH4_AB':'pk_maxCH4_AB'})
        datFram_wtLoc = weightedLoc(datFram_cent,'OB_LAT','OB_LON','OP_NUM','OB_CH4_AB').loc[:,:].rename(columns = {'OB_LAT':'pk_LAT','OB_LON':'pk_LON'})
        datFram_wtLocMax = pd.merge(datFram_wtLoc,maxch4,on = ['OP_NUM'])
        pass_info = datFram.copy()
        geometry_temp = [Point(xy) for xy in zip(datFram_wtLocMax['pk_LON'], datFram_wtLocMax['pk_LAT'])]
        crs = {'init': 'epsg:4326'}
        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)
        gdf_buff = gdf_buff.to_crs(epsg=32610)
        gdf_buff['geometry'] = gdf_buff.loc[:,'geometry'].buffer(30)
        gdf_tog = pd.merge(gdf_buff,datFram,on = ['OP_NUM'])
        gdf_bind_pks = gdf_buff.copy()
        gdf_pass_pks = gdf_bind_pks.copy()
        gdf_pass_pks['min_read']= gdf_pass_pks.loc[:,'OP_NUM']
        gdf_pass_pks['numtimes'] = 1
        gdf_pass_pks['numdays'] = 1
        gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:,'geometry']
        gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['OP_NUM']].to_numpy())].copy()
        gdf_pass_pks['dates'] = gdf_pass_pks.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM[6:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'),axis=1)
        gdf_pass_pks['pk_Dates'] = [list(x) for x in list(gdf_pass_pks.loc[:,['dates']].to_numpy())]
        gdf_pass_pks['min_Date'] = gdf_pass_pks.loc[:,'dates']
        gdf_pass_pks = gdf_pass_pks.drop(columns=['dates'])

        gdf_pass_pks['verified'] = False
#           gdf_pass_pks['oldgeo'] = gdf_pass_pks.loc[:,'geometry']
        gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
        together = pd.merge(gdf_pass_pks,gdf_tog,on = ['OP_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','geometry'])
        together['pass'] = whichpass
        gdf_pass_pks = together.copy()

        gdf_pass_pks['pkGEO'] = gdf_pass_pks.loc[:,"geometry"]
        gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
        del(gdf_pass_pks['newgeo'])
        gdf_pass_pks['pass'] = whichpass


        gdf_op_unique = gdf_pass_pks.loc[:,['numtimes','min_read','numdays','min_Date','verified','pass','OB_LON','OB_LAT']].drop_duplicates()
        gdfcop = gdf_pass_pks.loc[:,['OP_NUM','min_read','min_Date','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
        combinedOP = weightedLoc(gdfcop,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').loc[:,:].rename(columns = {'pk_LAT':'Overall_LAT','pk_LON':'Overall_LON'}).reset_index(drop=True)
        combinedOP1 = pd.merge(combinedOP,gdfcop,on=['min_read'])



        geometry_temp = [Point(xy) for xy in zip(combinedOP1['Overall_LON'], combinedOP1['Overall_LAT'])]
        crs = {'init': 'epsg:4326'}
        gdf_OP = gpd.GeoDataFrame(combinedOP1, crs=crs, geometry=geometry_temp)
        gdf_OP = gdf_OP.to_crs(epsg=32610).copy()

        gdf_OP_reduced = gdf_OP.loc[:,['min_read','geometry','numtimes','Overall_LON','Overall_LAT','min_Date','verified']].drop_duplicates().reset_index(drop=True)
        gdf_OP_reduced.to_file(new_loc_json, driver="GeoJSON")

        #gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']),gdf_pass_pks.drop(columns=['geometry','oldgeo']),on=['min_read','min_Date','numtimes','pass','verified','pk_LAT','pk_LON','OP_NUM','pk_maxCH4_AB'])
        gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']),gdf_pass_pks.drop(columns=['geometry']),on=['min_read','min_Date','numtimes','pass','verified','pk_LAT','pk_LON','OP_NUM','pk_maxCH4_AB'])

        #gdf_OP.to_csv(new_loc,index=False)
        gdf_OP_wrecombine.to_csv(new_loc,index=False)

            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)


            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        #gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)



        #gdfcop = gdfcop.to_crs(epsg=32610).copy()
        #gdfcop.to_file(new_loc_json, driver="GeoJSON")


        #gdf_tot.to_csv(new_loc, index = False)
        unique_peaks = gdf_pass_pks.loc[:,['OP_NUM','pk_LAT','pk_LON','min_read','min_Date']].drop_duplicates()
        unique_peaks['save'] = True
        #good_pks = unique_peaks.PEAK_NUM.drop_duplicates().values.tolist()
        good_pks = list(unique_peaks.index)

        def getthing(index):
            if index in good_pks:
                return True
            else:
                return False
        gdf_pass_pks['wooind'] = gdf_pass_pks.index
        gdf_pass_pks['save'] = gdf_pass_pks.apply(lambda x: getthing(x.wooind),axis=1)

       # unique_pks_tog = pd.concat([unique_peaks, gdf_pass_pks.drop(columns=['LON', 'LAT','PEAK_NUM'])], axis=1, join='inner')
       # testa = pd.merge(gdf_pass_pks, unique_peaks, how='left', on=['PEAK_NUM', 'pk_LAT','pk_LON','min_read'])

        unique_pks_tog = gdf_pass_pks.loc[gdf_pass_pks.save == True,:].reset_index(drop=True)
        unique_pks_tog['Latitude'] = unique_pks_tog.loc[:,'pk_LAT']
        unique_pks_tog['Longitude'] = unique_pks_tog.loc[:,'pk_LON']

        #unique_pks_tog.to_csv(new_loc2, index = False)
        unique_pks_tog.to_csv(new_loc, index = False)


        #return(gdf_OP_wrecombine)

    elif datFram.shape[0] != 1:
        datFram_cent =  datFram.copy()
        #datFram_cent['CH4_AB'] = datFram.loc[:,'CH4'].sub(datFram.loc[:,'CH4_BASELINE'], axis = 0)
        datFram_cent['OB_CH4_AB'] = datFram.loc[:,'OB_CH4'].sub(datFram.loc[:,'OB_CH4_BASELINE'], axis = 0)


        ### MAXCH4 is a df with the max methane (above baseline) in the given observed peak
        #maxch4 = datFram_cent.groupby('PEAK_NUM',as_index = False).CH4_AB.max().rename(columns = {'CH4_AB':'pk_maxCH4_AB'})
        maxch4 = datFram_cent.groupby('OP_NUM',as_index = False).OB_CH4_AB.max().rename(columns = {'OB_CH4_AB':'pk_maxCH4_AB'})

        ### FINDING WEIGHTED LOCATION OF THE OP, BY THE ABOVE BASELINE CH4 LEVEL
        #wtloc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB')
        #datFram_wtLoca =  wtloc.copy()
        #datFram_wtLoc = datFram_wtLoca.rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'})

        datFram_wtLoc = weightedLoc(datFram_cent,'OB_LAT','OB_LON','OP_NUM','OB_CH4_AB').loc[:,:].rename(columns = {'OB_LAT':'pk_LAT','OB_LON':'pk_LON'})
        #datFram_wtLoc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB').rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'}).copy()
        datFram_wtLocMax = pd.merge(datFram_wtLoc,maxch4,on = ['OP_NUM'])

        pass_info = datFram.copy()

        ## MIGHT NEED TO CHANGE BACK
        #geometry_temp = [Point(xy) for xy in zip(datFram['LON'], datFram['LAT'])]
        #crs = {'init': 'epsg:4326'}

        geometry_temp = [Point(xy) for xy in zip(datFram_wtLocMax['pk_LON'], datFram_wtLocMax['pk_LAT'])]
        crs = {'init': 'epsg:4326'}

            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        ## BUFFER AROUND EACH 'OP_NUM' OF 30 M
        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)
        #gdf_buff = makeGPD(datFram,'LON','LAT')
        gdf_buff = gdf_buff.to_crs(epsg=32610)
        gdf_buff['geometry'] = gdf_buff.loc[:,'geometry'].buffer(30)

        #pass_info_new = datFram.copy().rename(columns={"geometry": 'pk_geo'})

       # gdf_tog = pd.merge(gdf_buff,pass_info_new,on = ['PEAK_NUM', 'EPOCHSTART', 'EPOCH', 'DATETIME', 'CH4', 'LON', 'LAT',
       #    'CH4_BASELINE', 'CH4_THRESHOLD', 'PEAK_DIST_M', 'PEAK_CH4', 'TCH4',
       #    'PERIOD5MIN'])

        gdf_tog = pd.merge(gdf_buff,datFram,on = ['OP_NUM'])

        #gdf_bind_pks = gdf_tog.dissolve(by = 'PEAK_NUM',as_index=False).loc[:,['PEAK_NUM','geometry']]

        gdf_bind_pks = gdf_buff.copy()


        if gdf_bind_pks.shape[0] > 1:
            data_overlap = gpd.GeoDataFrame(crs=gdf_bind_pks.crs)
            data_temp = gdf_bind_pks.copy()
            for index, row in data_temp.iterrows():
                data_temp1=data_temp.loc[data_temp.OP_NUM!=row.OP_NUM,]
                # check if intersection occured
                overlaps=data_temp1[data_temp1.geometry.overlaps(row.geometry)]['OP_NUM'].tolist()
                if len(overlaps)>0:

                    # compare the area with threshold
                    for y in overlaps:
                        temp_area=gpd.overlay(data_temp.loc[data_temp.OP_NUM==y,],data_temp.loc[data_temp.OP_NUM==row.OP_NUM,],how='intersection')
                        temp_area=temp_area.loc[temp_area.geometry.area>=0.001]
                        if temp_area.shape[0]>0:
                            temp_union = gpd.overlay(data_temp.loc[data_temp.OP_NUM==y,],data_temp.loc[data_temp.OP_NUM==row.OP_NUM,],how='union')
                            data_overlap=gpd.GeoDataFrame(pd.concat([temp_union,data_overlap],ignore_index=True),crs=data_temp.crs)
            if data_overlap.size > 0:

                    firstnull2 = data_overlap.loc[data_overlap.OP_NUM_1.isnull(),:]
                    firstnull = firstnull2.copy()
                    firstnull.loc[:,'OP_NUM_1'] = firstnull2.loc[:,'OP_NUM_2']

                    secnull2 = data_overlap.loc[data_overlap.OP_NUM_2.isnull(),:]

                    secnull = secnull2.copy()
                    secnull.loc[:,'OP_NUM_2'] = secnull2.loc[:,'OP_NUM_1']

                    withoutNA = data_overlap.copy().dropna()
                    allTog2 = pd.concat([firstnull,secnull,withoutNA]).reset_index().copy()


                    allTog2['notsame'] = allTog2.apply(lambda x:x.OP_NUM_1 == x.OP_NUM_2,axis=1)
                    allTog = allTog2.loc[allTog2.notsame == False,:].drop(columns = ['notsame'])


                    over = allTog.copy()
                    over['sorted']=over.apply(lambda y: sorted([y['OP_NUM_1'],y['OP_NUM_2']]),axis=1)
                    over['sorted']=over.sorted.apply(lambda y: ''.join(y))
                    over = over.drop_duplicates('sorted')
                    over['combined']= [list(x) for x in list(over.loc[:,['OP_NUM_1','OP_NUM_2']].to_numpy())]
                    #over['date1'] = over.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM_1'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                    #over['date2'] = over.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM_2'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                    over['date1'] = over.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM_1[6:x.OP_NUM_1.find('.')])).strftime('%Y-%m-%d'),axis=1)
                    over['date2'] = over.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM_2[6:x.OP_NUM_2.find('.')])).strftime('%Y-%m-%d'),axis=1)



                    def unique(list1):
                        # intilize a null list
                        unique_list = []

                        # traverse for all elements
                        for x in list1:
                            # check if exists in unique_list or not
                            if x not in unique_list:
                                unique_list.append(x)
                        return(unique_list)


                    over['dates']= [list(x) for x in list(over.loc[:,['date1','date2']].to_numpy())]
                    over['pk_Dates'] = over.apply(lambda x: unique(x.dates),axis=1)
                    over = over.drop(columns = ['dates'])


                    over['VER_NUM'] = over.apply(lambda y: y.combined,axis=1)
                    over['min_val']=over.apply(lambda y: min(y.combined),axis=1)
                    over2=over.reset_index().loc[:,['OP_NUM_1','OP_NUM_2','geometry','combined','min_val','pk_Dates']]

                    overcop = over2.copy().rename(columns = {'combined':'recombine'})
                    #overcop.loc[:,'recombine'] = overcop.loc[:,'combined']

                    for index, row in overcop.iterrows():
                        united = row.recombine
                        undate = row.pk_Dates
                        for index2, row2 in overcop.iterrows():
                            united_temp = unIfInt(united,row2.recombine)
                            undate_temp = unIfInt(undate,row2.pk_Dates)
                            if united_temp != None:
                                united = united_temp
                            if undate_temp != None:
                                undate = undate_temp
                        overcop.at[index, 'recombine']  = united.copy()
                        overcop.at[index, 'pk_Dates']  = undate.copy()

                        del(united)
                        del(undate)


                    overcop['recombine']= overcop.apply(lambda y: sorted(y.recombine),axis=1).copy()
                    overcop['pk_Dates']= overcop.apply(lambda y: sorted(y.pk_Dates),axis=1).copy()
                    overcop['min_read'] = overcop.apply(lambda y: min(y.recombine),axis=1).copy()
                    overcop['min_Date'] = overcop.apply(lambda y: min(y.pk_Dates),axis=1).copy()

                    newOverlap = overcop.dissolve(by='min_read',as_index=False).loc[:,['min_read','geometry','recombine','min_Date','pk_Dates']].copy()


                    combined = gdf_bind_pks.copy()
                    combined['recombine'] = [list(x) for x in list(combined.loc[:,['OP_NUM']].to_numpy())]
                    #combined['dates'] = combined.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                    combined['dates'] = combined.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM[6:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'),axis=1)



                    combined['pk_Dates'] = [list(x) for x in list(combined.loc[:,['dates']].to_numpy())]

                    combined['min_Date'] = combined.loc[:,'dates']

                    combined['numtimes'] = 1
                    combined['newgeo'] = combined.loc[:,'geometry']
                    combined['min_read'] = combined.loc[:,"OP_NUM"]

                    for index,row in combined.iterrows():
                        for index2,row2 in newOverlap.iterrows():
                            if row.OP_NUM in row2.recombine:
                                combined.at[index, 'recombine']  = row2.recombine.copy()
                                #combined.at[index, 'newgeo']  = row2.copy().geometry
                                combined.at[index,'min_read'] = row2.copy().min_read
                                combined.at[index, 'pk_Dates']  = row2.pk_Dates
                                combined.at[index,'min_Date'] = row2.min_Date

                    #combined['numtimes'] = combined.apply(lambda y: len(y.recombine),axis = 1).copy()
                    combined['numtimes'] = combined.apply(lambda x: countTimes(x.recombine),axis=1)

                    combined['numdays'] = combined.apply(lambda y: len(y.pk_Dates),axis = 1).copy()
                    combined_reduced = combined.loc[:,['OP_NUM','newgeo','recombine','numtimes','min_read','numdays','pk_Dates','min_Date']]
                    gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['OP_NUM']).copy()
                    gdf_pass_pks['verified'] = gdf_pass_pks.apply(lambda y: (True if y.numtimes > 1 else False),axis=1 ).copy()
            if data_overlap.size == 0:
               gdf_pass_pks = gdf_bind_pks.copy()
               gdf_pass_pks['min_read']= gdf_pass_pks.loc[:,'OP_NUM']
               gdf_pass_pks['numtimes'] = 1
               gdf_pass_pks['numdays'] = 1

               gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:,'geometry']
               gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['OP_NUM']].to_numpy())].copy()
               #gdf_pass_pks['dates'] = gdf_pass_pks.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
               gdf_pass_pks['dates'] = gdf_pass_pks.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM[6:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'),axis=1)

               gdf_pass_pks['pk_Dates'] = [list(x) for x in list(gdf_pass_pks.loc[:,['dates']].to_numpy())]
               gdf_pass_pks['min_Date'] = gdf_pass_pks.loc[:,'dates']
               gdf_pass_pks = gdf_pass_pks.drop(columns=['dates'])

               gdf_pass_pks['verified'] = False
    #           gdf_pass_pks['oldgeo'] = gdf_pass_pks.loc[:,'geometry']
               gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
               together = pd.merge(gdf_pass_pks,gdf_tog,on = ['OP_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','geometry'])
               together['pass'] = whichpass
               gdf_pass_pks = together.copy()




        if gdf_bind_pks.shape[0] == 1:
            gdf_pass_pks = gdf_bind_pks.copy()
            gdf_pass_pks['min_read']= gdf_pass_pks.loc[:,'OP_NUM']
            gdf_pass_pks['numtimes'] = 1
            gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:,'geometry']

            gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['OP_NUM']].to_numpy())].copy()
           # gdf_pass_pks['dates'] = gdf_pass_pks.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
            gdf_pass_pks['dates'] = gdf_pass_pks.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM[6:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'),axis=1)

            gdf_pass_pks['pk_Dates'] = [list(x) for x in list(gdf_pass_pks.loc[:,['dates']].to_numpy())]
            gdf_pass_pks['min_Date'] = gdf_pass_pks.loc[:,'dates']
            gdf_pass_pks['numdays'] = 1
            gdf_pass_pks = gdf_pass_pks.drop(columns=['dates'])


            gdf_pass_pks['verified'] = False
            epdat = pass_info.loc[:,['OP_NUM','OP_EPOCHSTART']]
            gdf_pass_pks = pd.merge(gdf_pass_pks,epdat,on = ['OP_NUM']).copy()
            data_overlap = pd.DataFrame(columns = ['what','oh'])

        ### gdf_pass_pks
        #    Index(['OP_NUM', 'pk_LON', 'pk_LAT', 'pk_maxCH4_AB', 'geometry',
        #       'OP_EPOCHSTART', 'OB_EPOCH', 'OB_DATETIME', 'OB_CH4', 'OB_LON',
        #       'OB_LAT', 'OB_CH4_BASELINE', 'OB_CH4_THRESHOLD', 'OP_PEAK_DIST_M',
        #       'OP_PEAK_CH4', 'OB_TCH4', 'OB_PERIOD5MIN', 'OB_CH4_AB', 'newgeo',
        #       'recombine', 'numtimes', 'min_read', 'numdays', 'pk_Dates', 'min_Date',
        #       'verified'],
        #      dtype='object')

        gdf_pass_pks['pkGEO'] = gdf_pass_pks.loc[:,"geometry"]
        gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
        del(gdf_pass_pks['newgeo'])
        gdf_pass_pks['pass'] = whichpass
        gdf_pass_pks['Overall_LON'] = gdf_pass_pks['pk_LON']
        gdf_pass_pks['Overall_LAT'] = gdf_pass_pks['pk_LAT']
        combinedOP1 = gdf_pass_pks.drop(columns = ['recombine','pk_Dates']).drop_duplicates()

        #gdf_tot = pd.merge(gdf_pass_pks,datFram_wtLocMax.loc[:,['PEAK_NUM','pk_LON','pk_LAT']],on = ['PEAK_NUM','pk_LON','pk_LAT']).copy()
        ## condense by peak_num
       # gdfcop = gdf_tot.loc[:,['PEAK_NUM','geometry','min_read','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()

       #### WANT A DATAFRAME WITH
       # EACH OP SUMMARY
       # COMBINED WITH THE COMBINED SUMMARY
       # gdf_op_unique = gdf_pass_pks.loc[:,['numtimes','min_read','numdays','min_Date','verified','pass','Overall_LON','Overall_LAT']].drop_duplicates()

        if data_overlap.size != 0:
            gdf_op_unique = gdf_pass_pks.loc[:,['numtimes','min_read','numdays','min_Date','verified','pass','OB_LON','OB_LAT']].drop_duplicates()
            gdfcop = gdf_pass_pks.loc[:,['OP_NUM','min_read','min_Date','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
            combinedOP = weightedLoc(gdfcop,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').loc[:,:].rename(columns = {'pk_LAT':'Overall_LAT','pk_LON':'Overall_LON'}).reset_index(drop=True)
            combinedOP1 = pd.merge(combinedOP,gdfcop,on=['min_read'])

        if data_overlap.size == 0 and gdf_bind_pks.shape[0] != 1:
            gdf_op_unique = gdf_pass_pks.loc[:,['numtimes','min_read','numdays','min_Date','verified','pass','OB_LON','OB_LAT']].drop_duplicates()
            gdfcop = gdf_pass_pks.loc[:,['OP_NUM','min_read','min_Date','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
            combinedOP = weightedLoc(gdfcop,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').loc[:,:].rename(columns = {'pk_LAT':'Overall_LAT','pk_LON':'Overall_LON'}).reset_index(drop=True)
            combinedOP1 = pd.merge(combinedOP,gdfcop,on=['min_read'])


        ## TO FINDED WEIGHTED LOCATION OF EACH PK GROUP
    #    gdfcop = gdf_pass_pks.loc[:,['OP_NUM','min_read','min_Date','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    #    combinedOP = weightedLoc(gdfcop,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').loc[:,:].rename(columns = {'pk_LAT':'Overall_LAT','pk_LON':'Overall_LON'}).reset_index(drop=True)
    #    combinedOP1 = pd.merge(combinedOP,gdfcop,on=['min_read'])

        ## getting the rest of the stuff
        #gdf_justloc = gdfcop.loc[:,['min_read','pk_LAT','pk_LON','min_Date']].reset_index(drop=True)


        #other = gdfcop.loc[:,['min_read','numtimes','min_Date']].reset_index(drop=True)


        geometry_temp = [Point(xy) for xy in zip(combinedOP1['Overall_LON'], combinedOP1['Overall_LAT'])]
        crs = {'init': 'epsg:4326'}
        gdf_OP = gpd.GeoDataFrame(combinedOP1, crs=crs, geometry=geometry_temp)
        gdf_OP = gdf_OP.to_crs(epsg=32610).copy()

        gdf_OP_reduced = gdf_OP.loc[:,['min_read','geometry','numtimes','Overall_LON','Overall_LAT','min_Date','verified']].drop_duplicates().reset_index(drop=True)
        gdf_OP_reduced.to_file(new_loc_json, driver="GeoJSON")

        #gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']),gdf_pass_pks.drop(columns=['geometry','oldgeo']),on=['min_read','min_Date','numtimes','pass','verified','pk_LAT','pk_LON','OP_NUM','pk_maxCH4_AB'])
        gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']),gdf_pass_pks.drop(columns=['geometry']),on=['min_read','min_Date','numtimes','pass','verified','pk_LAT','pk_LON','OP_NUM','pk_maxCH4_AB'])

        #gdf_OP.to_csv(new_loc,index=False)
        gdf_OP_wrecombine.to_csv(new_loc,index=False)

            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)


            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        #gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)



        #gdfcop = gdfcop.to_crs(epsg=32610).copy()
        #gdfcop.to_file(new_loc_json, driver="GeoJSON")


        #gdf_tot.to_csv(new_loc, index = False)
        unique_peaks = gdf_pass_pks.loc[:,['OP_NUM','pk_LAT','pk_LON','min_read','min_Date']].drop_duplicates()
        unique_peaks['save'] = True
        #good_pks = unique_peaks.PEAK_NUM.drop_duplicates().values.tolist()
        good_pks = list(unique_peaks.index)

        def getthing(index):
            if index in good_pks:
                return True
            else:
                return False
        gdf_pass_pks['wooind'] = gdf_pass_pks.index
        gdf_pass_pks['save'] = gdf_pass_pks.apply(lambda x: getthing(x.wooind),axis=1)

       # unique_pks_tog = pd.concat([unique_peaks, gdf_pass_pks.drop(columns=['LON', 'LAT','PEAK_NUM'])], axis=1, join='inner')
       # testa = pd.merge(gdf_pass_pks, unique_peaks, how='left', on=['PEAK_NUM', 'pk_LAT','pk_LON','min_read'])

        unique_pks_tog = gdf_pass_pks.loc[gdf_pass_pks.save == True,:].reset_index(drop=True)
        unique_pks_tog['Latitude'] = unique_pks_tog.loc[:,'pk_LAT']
        unique_pks_tog['Longitude'] = unique_pks_tog.loc[:,'pk_LON']

        #unique_pks_tog.to_csv(new_loc2, index = False)
        unique_pks_tog.to_csv(new_loc, index = False)

        return()
        #return(gdf_OP_wrecombine)


def filterPeak_swift(xCar,xDate,xDir,xFilename, outFolder,whichpass = 0):
    import pandas as pd #
    import geopandas as gpd
    import shutil
    from datetime import datetime
    from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry

    file_loc = xDir + xFilename
    new_loc = outFolder + "Filtered" + xFilename
    new_loc_json = new_loc[:-3] + 'json'

    oldInfo = xDir + 'Peaks_' + xCar + "_" + xDate.replace("-","") + "_info.csv"
    newInfo = outFolder + 'FilteredPeaks_' + xCar + "_" + xDate.replace("-","") + "_info.csv"

    shutil.copy(oldInfo,newInfo)

    # identified peaks has the columns:
    #Index(['OP_NUM', 'OP_EPOCHSTART', 'OB_EPOCH', 'OB_DATETIME', 'OB_CH4',
    #   'OB_LON', 'OB_LAT', 'OB_CH4_BASELINE', 'OB_CH4_THRESHOLD',
    #   'OP_PEAK_DIST_M', 'OP_PEAK_CH4', 'OB_TCH4', 'OB_PERIOD5MIN'],
    #  dtype='object')

    datFram = pd.read_csv(file_loc)

    if datFram.shape[0] == 0:
        print("Not filtering this file, no peak in it!")
    elif datFram.shape[0] == 1: ## only one thing to begin with
        datFram_cent =  datFram.loc[:,:]
        datFram_cent['OB_CH4_AB'] = datFram.loc[:,'OB_CH4'].sub(datFram.loc[:,'OB_CH4_BASELINE'], axis = 0)
        maxch4 = datFram_cent.groupby('OP_NUM',as_index = False).OB_CH4_AB.max().rename(columns = {'OB_CH4_AB':'pk_maxCH4_AB'})
        datFram_wtLoc = weightedLoc(datFram_cent,'OB_LAT','OB_LON','OP_NUM','OB_CH4_AB').loc[:,:].rename(columns = {'OB_LAT':'pk_LAT','OB_LON':'pk_LON'})
        datFram_wtLocMax = pd.merge(datFram_wtLoc,maxch4,on = ['OP_NUM'])
        pass_info = datFram.copy()
        geometry_temp = [Point(xy) for xy in zip(datFram_wtLocMax['pk_LON'], datFram_wtLocMax['pk_LAT'])]
        crs = {'init': 'epsg:4326'}
        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)
        gdf_buff = gdf_buff.to_crs(epsg=32610)
        gdf_buff['geometry'] = gdf_buff.loc[:,'geometry'].buffer(30)
        gdf_tog = pd.merge(gdf_buff,datFram,on = ['OP_NUM'])
        gdf_bind_pks = gdf_buff.copy()
        gdf_pass_pks = gdf_bind_pks.copy()
        gdf_pass_pks['min_read']= gdf_pass_pks.loc[:,'OP_NUM']
        gdf_pass_pks['numtimes'] = 1
        gdf_pass_pks['numdays'] = 1
        gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:,'geometry']
        gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['OP_NUM']].to_numpy())].copy()
        gdf_pass_pks['dates'] = gdf_pass_pks.swifter.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM[6:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'),axis=1)
        gdf_pass_pks['pk_Dates'] = [list(x) for x in list(gdf_pass_pks.loc[:,['dates']].to_numpy())]
        gdf_pass_pks['min_Date'] = gdf_pass_pks.loc[:,'dates']
        gdf_pass_pks = gdf_pass_pks.drop(columns=['dates'])

        gdf_pass_pks['verified'] = False
#           gdf_pass_pks['oldgeo'] = gdf_pass_pks.loc[:,'geometry']
        gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
        together = pd.merge(gdf_pass_pks,gdf_tog,on = ['OP_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','geometry'])
        together['pass'] = whichpass
        gdf_pass_pks = together.copy()

        gdf_pass_pks['pkGEO'] = gdf_pass_pks.loc[:,"geometry"]
        gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
        del(gdf_pass_pks['newgeo'])
        gdf_pass_pks['pass'] = whichpass


        gdf_op_unique = gdf_pass_pks.loc[:,['numtimes','min_read','numdays','min_Date','verified','pass','OB_LON','OB_LAT']].drop_duplicates()
        gdfcop = gdf_pass_pks.loc[:,['OP_NUM','min_read','min_Date','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
        combinedOP = weightedLoc(gdfcop,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').loc[:,:].rename(columns = {'pk_LAT':'Overall_LAT','pk_LON':'Overall_LON'}).reset_index(drop=True)
        combinedOP1 = pd.merge(combinedOP,gdfcop,on=['min_read'])



        geometry_temp = [Point(xy) for xy in zip(combinedOP1['Overall_LON'], combinedOP1['Overall_LAT'])]
        crs = {'init': 'epsg:4326'}
        gdf_OP = gpd.GeoDataFrame(combinedOP1, crs=crs, geometry=geometry_temp)
        gdf_OP = gdf_OP.to_crs(epsg=32610).copy()

        gdf_OP_reduced = gdf_OP.loc[:,['min_read','geometry','numtimes','Overall_LON','Overall_LAT','min_Date','verified']].drop_duplicates().reset_index(drop=True)
        gdf_OP_reduced.to_file(new_loc_json, driver="GeoJSON")

        #gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']),gdf_pass_pks.drop(columns=['geometry','oldgeo']),on=['min_read','min_Date','numtimes','pass','verified','pk_LAT','pk_LON','OP_NUM','pk_maxCH4_AB'])
        gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']),gdf_pass_pks.drop(columns=['geometry']),on=['min_read','min_Date','numtimes','pass','verified','pk_LAT','pk_LON','OP_NUM','pk_maxCH4_AB'])

        #gdf_OP.to_csv(new_loc,index=False)
        gdf_OP_wrecombine.to_csv(new_loc,index=False)

            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)


            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        #gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)



        #gdfcop = gdfcop.to_crs(epsg=32610).copy()
        #gdfcop.to_file(new_loc_json, driver="GeoJSON")


        #gdf_tot.to_csv(new_loc, index = False)
        unique_peaks = gdf_pass_pks.loc[:,['OP_NUM','pk_LAT','pk_LON','min_read','min_Date']].drop_duplicates()
        unique_peaks['save'] = True
        #good_pks = unique_peaks.PEAK_NUM.drop_duplicates().values.tolist()
        good_pks = list(unique_peaks.index)

        def getthing(index):
            if index in good_pks:
                return True
            else:
                return False
        gdf_pass_pks['wooind'] = gdf_pass_pks.index
        gdf_pass_pks['save'] = gdf_pass_pks.swifter.apply(lambda x: getthing(x.wooind),axis=1)

       # unique_pks_tog = pd.concat([unique_peaks, gdf_pass_pks.drop(columns=['LON', 'LAT','PEAK_NUM'])], axis=1, join='inner')
       # testa = pd.merge(gdf_pass_pks, unique_peaks, how='left', on=['PEAK_NUM', 'pk_LAT','pk_LON','min_read'])

        unique_pks_tog = gdf_pass_pks.loc[gdf_pass_pks.save == True,:].reset_index(drop=True)
        unique_pks_tog['Latitude'] = unique_pks_tog.loc[:,'pk_LAT']
        unique_pks_tog['Longitude'] = unique_pks_tog.loc[:,'pk_LON']

        #unique_pks_tog.to_csv(new_loc2, index = False)
        unique_pks_tog.to_csv(new_loc, index = False)


        #return(gdf_OP_wrecombine)

    elif datFram.shape[0] != 1:
        datFram_cent =  datFram.copy()
        #datFram_cent['CH4_AB'] = datFram.loc[:,'CH4'].sub(datFram.loc[:,'CH4_BASELINE'], axis = 0)
        datFram_cent['OB_CH4_AB'] = datFram.loc[:,'OB_CH4'].sub(datFram.loc[:,'OB_CH4_BASELINE'], axis = 0)


        ### MAXCH4 is a df with the max methane (above baseline) in the given observed peak
        #maxch4 = datFram_cent.groupby('PEAK_NUM',as_index = False).CH4_AB.max().rename(columns = {'CH4_AB':'pk_maxCH4_AB'})
        maxch4 = datFram_cent.groupby('OP_NUM',as_index = False).OB_CH4_AB.max().rename(columns = {'OB_CH4_AB':'pk_maxCH4_AB'})

        ### FINDING WEIGHTED LOCATION OF THE OP, BY THE ABOVE BASELINE CH4 LEVEL
        #wtloc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB')
        #datFram_wtLoca =  wtloc.copy()
        #datFram_wtLoc = datFram_wtLoca.rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'})

        datFram_wtLoc = weightedLoc(datFram_cent,'OB_LAT','OB_LON','OP_NUM','OB_CH4_AB').loc[:,:].rename(columns = {'OB_LAT':'pk_LAT','OB_LON':'pk_LON'})
        #datFram_wtLoc = weightedLoc(datFram_cent,'LAT','LON','PEAK_NUM','CH4_AB').rename(columns = {'LAT':'pk_LAT','LON':'pk_LON'}).copy()
        datFram_wtLocMax = pd.merge(datFram_wtLoc,maxch4,on = ['OP_NUM'])

        pass_info = datFram.copy()

        ## MIGHT NEED TO CHANGE BACK
        #geometry_temp = [Point(xy) for xy in zip(datFram['LON'], datFram['LAT'])]
        #crs = {'init': 'epsg:4326'}

        geometry_temp = [Point(xy) for xy in zip(datFram_wtLocMax['pk_LON'], datFram_wtLocMax['pk_LAT'])]
        crs = {'init': 'epsg:4326'}

            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        ## BUFFER AROUND EACH 'OP_NUM' OF 30 M
        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)
        #gdf_buff = makeGPD(datFram,'LON','LAT')
        gdf_buff = gdf_buff.to_crs(epsg=32610)
        gdf_buff['geometry'] = gdf_buff.loc[:,'geometry'].buffer(30)

        #pass_info_new = datFram.copy().rename(columns={"geometry": 'pk_geo'})

       # gdf_tog = pd.merge(gdf_buff,pass_info_new,on = ['PEAK_NUM', 'EPOCHSTART', 'EPOCH', 'DATETIME', 'CH4', 'LON', 'LAT',
       #    'CH4_BASELINE', 'CH4_THRESHOLD', 'PEAK_DIST_M', 'PEAK_CH4', 'TCH4',
       #    'PERIOD5MIN'])

        gdf_tog = pd.merge(gdf_buff,datFram,on = ['OP_NUM'])

        #gdf_bind_pks = gdf_tog.dissolve(by = 'PEAK_NUM',as_index=False).loc[:,['PEAK_NUM','geometry']]

        gdf_bind_pks = gdf_buff.copy()


        if gdf_bind_pks.shape[0] > 1:
            data_overlap = gpd.GeoDataFrame(crs=gdf_bind_pks.crs)
            data_temp = gdf_bind_pks.copy()
            for index, row in data_temp.iterrows():
                data_temp1=data_temp.loc[data_temp.OP_NUM!=row.OP_NUM,]
                # check if intersection occured
                overlaps=data_temp1[data_temp1.geometry.overlaps(row.geometry)]['OP_NUM'].tolist()
                if len(overlaps)>0:

                    # compare the area with threshold
                    for y in overlaps:
                        temp_area=gpd.overlay(data_temp.loc[data_temp.OP_NUM==y,],data_temp.loc[data_temp.OP_NUM==row.OP_NUM,],how='intersection')
                        temp_area=temp_area.loc[temp_area.geometry.area>=0.001]
                        if temp_area.shape[0]>0:
                            temp_union = gpd.overlay(data_temp.loc[data_temp.OP_NUM==y,],data_temp.loc[data_temp.OP_NUM==row.OP_NUM,],how='union')
                            data_overlap=gpd.GeoDataFrame(pd.concat([temp_union,data_overlap],ignore_index=True),crs=data_temp.crs)
            if data_overlap.size > 0:

                    firstnull2 = data_overlap.loc[data_overlap.OP_NUM_1.isnull(),:]
                    firstnull = firstnull2.copy()
                    firstnull.loc[:,'OP_NUM_1'] = firstnull2.loc[:,'OP_NUM_2']

                    secnull2 = data_overlap.loc[data_overlap.OP_NUM_2.isnull(),:]

                    secnull = secnull2.copy()
                    secnull.loc[:,'OP_NUM_2'] = secnull2.loc[:,'OP_NUM_1']

                    withoutNA = data_overlap.copy().dropna()
                    allTog2 = pd.concat([firstnull,secnull,withoutNA]).reset_index().copy()


                    allTog2['notsame'] = allTog2.swifter.apply(lambda x:x.OP_NUM_1 == x.OP_NUM_2,axis=1)
                    allTog = allTog2.loc[allTog2.notsame == False,:].drop(columns = ['notsame'])


                    over = allTog.copy()
                    over['sorted']=over.swifter.apply(lambda y: sorted([y['OP_NUM_1'],y['OP_NUM_2']]),axis=1)
                    over['sorted']=over.sorted.swifter.apply(lambda y: ''.join(y))
                    over = over.drop_duplicates('sorted')
                    over['combined']= [list(x) for x in list(over.loc[:,['OP_NUM_1','OP_NUM_2']].to_numpy())]
                    #over['date1'] = over.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM_1'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                    #over['date2'] = over.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM_2'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                    over['date1'] = over.swifter.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM_1[6:x.OP_NUM_1.find('.')])).strftime('%Y-%m-%d'),axis=1)
                    over['date2'] = over.swifter.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM_2[6:x.OP_NUM_2.find('.')])).strftime('%Y-%m-%d'),axis=1)



                    def unique(list1):
                        # intilize a null list
                        unique_list = []

                        # traverse for all elements
                        for x in list1:
                            # check if exists in unique_list or not
                            if x not in unique_list:
                                unique_list.append(x)
                        return(unique_list)


                    over['dates']= [list(x) for x in list(over.loc[:,['date1','date2']].to_numpy())]
                    over['pk_Dates'] = over.swifter.apply(lambda x: unique(x.dates),axis=1)
                    over = over.drop(columns = ['dates'])


                    over['VER_NUM'] = over.swifter.apply(lambda y: y.combined,axis=1)
                    over['min_val']=over.swifter.apply(lambda y: min(y.combined),axis=1)
                    over2=over.reset_index().loc[:,['OP_NUM_1','OP_NUM_2','geometry','combined','min_val','pk_Dates']]

                    overcop = over2.copy().rename(columns = {'combined':'recombine'})
                    #overcop.loc[:,'recombine'] = overcop.loc[:,'combined']

                    for index, row in overcop.iterrows():
                        united = row.recombine
                        undate = row.pk_Dates
                        for index2, row2 in overcop.iterrows():
                            united_temp = unIfInt(united,row2.recombine)
                            undate_temp = unIfInt(undate,row2.pk_Dates)
                            if united_temp != None:
                                united = united_temp
                            if undate_temp != None:
                                undate = undate_temp
                        overcop.at[index, 'recombine']  = united.copy()
                        overcop.at[index, 'pk_Dates']  = undate.copy()

                        del(united)
                        del(undate)


                    overcop['recombine']= overcop.swifter.apply(lambda y: sorted(y.recombine),axis=1).copy()
                    overcop['pk_Dates']= overcop.swifter.apply(lambda y: sorted(y.pk_Dates),axis=1).copy()
                    overcop['min_read'] = overcop.swifter.apply(lambda y: min(y.recombine),axis=1).copy()
                    overcop['min_Date'] = overcop.swifter.apply(lambda y: min(y.pk_Dates),axis=1).copy()

                    newOverlap = overcop.dissolve(by='min_read',as_index=False).loc[:,['min_read','geometry','recombine','min_Date','pk_Dates']].copy()


                    combined = gdf_bind_pks.copy()
                    combined['recombine'] = [list(x) for x in list(combined.loc[:,['OP_NUM']].to_numpy())]
                    #combined['dates'] = combined.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                    combined['dates'] = combined.swifter.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM[6:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'),axis=1)



                    combined['pk_Dates'] = [list(x) for x in list(combined.loc[:,['dates']].to_numpy())]

                    combined['min_Date'] = combined.loc[:,'dates']

                    combined['numtimes'] = 1
                    combined['newgeo'] = combined.loc[:,'geometry']
                    combined['min_read'] = combined.loc[:,"OP_NUM"]

                    for index,row in combined.iterrows():
                        for index2,row2 in newOverlap.iterrows():
                            if row.OP_NUM in row2.recombine:
                                combined.at[index, 'recombine']  = row2.recombine.copy()
                                #combined.at[index, 'newgeo']  = row2.copy().geometry
                                combined.at[index,'min_read'] = row2.copy().min_read
                                combined.at[index, 'pk_Dates']  = row2.pk_Dates
                                combined.at[index,'min_Date'] = row2.min_Date

                    #combined['numtimes'] = combined.apply(lambda y: len(y.recombine),axis = 1).copy()
                    combined['numtimes'] = combined.swifter.apply(lambda x: countTimes(x.recombine),axis=1)

                    combined['numdays'] = combined.swifter.apply(lambda y: len(y.pk_Dates),axis = 1).copy()
                    combined_reduced = combined.loc[:,['OP_NUM','newgeo','recombine','numtimes','min_read','numdays','pk_Dates','min_Date']]
                    gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['OP_NUM']).copy()
                    gdf_pass_pks['verified'] = gdf_pass_pks.swifter.apply(lambda y: (True if y.numtimes > 1 else False),axis=1 ).copy()
            if data_overlap.size == 0:
               gdf_pass_pks = gdf_bind_pks.copy()
               gdf_pass_pks['min_read']= gdf_pass_pks.loc[:,'OP_NUM']
               gdf_pass_pks['numtimes'] = 1
               gdf_pass_pks['numdays'] = 1

               gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:,'geometry']
               gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['OP_NUM']].to_numpy())].copy()
               #gdf_pass_pks['dates'] = gdf_pass_pks.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
               gdf_pass_pks['dates'] = gdf_pass_pks.swifter.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM[6:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'),axis=1)

               gdf_pass_pks['pk_Dates'] = [list(x) for x in list(gdf_pass_pks.loc[:,['dates']].to_numpy())]
               gdf_pass_pks['min_Date'] = gdf_pass_pks.loc[:,'dates']
               gdf_pass_pks = gdf_pass_pks.drop(columns=['dates'])

               gdf_pass_pks['verified'] = False
    #           gdf_pass_pks['oldgeo'] = gdf_pass_pks.loc[:,'geometry']
               gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
               together = pd.merge(gdf_pass_pks,gdf_tog,on = ['OP_NUM','pk_LON','pk_LAT','pk_maxCH4_AB','geometry'])
               together['pass'] = whichpass
               gdf_pass_pks = together.copy()




        if gdf_bind_pks.shape[0] == 1:
            gdf_pass_pks = gdf_bind_pks.copy()
            gdf_pass_pks['min_read']= gdf_pass_pks.loc[:,'OP_NUM']
            gdf_pass_pks['numtimes'] = 1
            gdf_pass_pks['newgeo'] = gdf_pass_pks.loc[:,'geometry']

            gdf_pass_pks['recombine'] = [list(x) for x in list(gdf_pass_pks.loc[:,['OP_NUM']].to_numpy())].copy()
           # gdf_pass_pks['dates'] = gdf_pass_pks.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
            gdf_pass_pks['dates'] = gdf_pass_pks.swifter.apply(lambda x: datetime.fromtimestamp(int( x.OP_NUM[6:x.OP_NUM.find('.')])).strftime('%Y-%m-%d'),axis=1)

            gdf_pass_pks['pk_Dates'] = [list(x) for x in list(gdf_pass_pks.loc[:,['dates']].to_numpy())]
            gdf_pass_pks['min_Date'] = gdf_pass_pks.loc[:,'dates']
            gdf_pass_pks['numdays'] = 1
            gdf_pass_pks = gdf_pass_pks.drop(columns=['dates'])


            gdf_pass_pks['verified'] = False
            epdat = pass_info.loc[:,['OP_NUM','OP_EPOCHSTART']]
            gdf_pass_pks = pd.merge(gdf_pass_pks,epdat,on = ['OP_NUM']).copy()
            data_overlap = pd.DataFrame(columns = ['what','oh'])

        ### gdf_pass_pks
        #    Index(['OP_NUM', 'pk_LON', 'pk_LAT', 'pk_maxCH4_AB', 'geometry',
        #       'OP_EPOCHSTART', 'OB_EPOCH', 'OB_DATETIME', 'OB_CH4', 'OB_LON',
        #       'OB_LAT', 'OB_CH4_BASELINE', 'OB_CH4_THRESHOLD', 'OP_PEAK_DIST_M',
        #       'OP_PEAK_CH4', 'OB_TCH4', 'OB_PERIOD5MIN', 'OB_CH4_AB', 'newgeo',
        #       'recombine', 'numtimes', 'min_read', 'numdays', 'pk_Dates', 'min_Date',
        #       'verified'],
        #      dtype='object')

        gdf_pass_pks['pkGEO'] = gdf_pass_pks.loc[:,"geometry"]
        gdf_pass_pks['geometry'] = gdf_pass_pks.loc[:,"newgeo"]
        del(gdf_pass_pks['newgeo'])
        gdf_pass_pks['pass'] = whichpass
        gdf_pass_pks['Overall_LON'] = gdf_pass_pks['pk_LON']
        gdf_pass_pks['Overall_LAT'] = gdf_pass_pks['pk_LAT']
        combinedOP1 = gdf_pass_pks.drop(columns = ['recombine','pk_Dates']).drop_duplicates()

        #gdf_tot = pd.merge(gdf_pass_pks,datFram_wtLocMax.loc[:,['PEAK_NUM','pk_LON','pk_LAT']],on = ['PEAK_NUM','pk_LON','pk_LAT']).copy()
        ## condense by peak_num
       # gdfcop = gdf_tot.loc[:,['PEAK_NUM','geometry','min_read','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()

       #### WANT A DATAFRAME WITH
       # EACH OP SUMMARY
       # COMBINED WITH THE COMBINED SUMMARY
       # gdf_op_unique = gdf_pass_pks.loc[:,['numtimes','min_read','numdays','min_Date','verified','pass','Overall_LON','Overall_LAT']].drop_duplicates()

        if data_overlap.size != 0:
            gdf_op_unique = gdf_pass_pks.loc[:,['numtimes','min_read','numdays','min_Date','verified','pass','OB_LON','OB_LAT']].drop_duplicates()
            gdfcop = gdf_pass_pks.loc[:,['OP_NUM','min_read','min_Date','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
            combinedOP = weightedLoc(gdfcop,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').loc[:,:].rename(columns = {'pk_LAT':'Overall_LAT','pk_LON':'Overall_LON'}).reset_index(drop=True)
            combinedOP1 = pd.merge(combinedOP,gdfcop,on=['min_read'])

        if data_overlap.size == 0 and gdf_bind_pks.shape[0] != 1:
            gdf_op_unique = gdf_pass_pks.loc[:,['numtimes','min_read','numdays','min_Date','verified','pass','OB_LON','OB_LAT']].drop_duplicates()
            gdfcop = gdf_pass_pks.loc[:,['OP_NUM','min_read','min_Date','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
            combinedOP = weightedLoc(gdfcop,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').loc[:,:].rename(columns = {'pk_LAT':'Overall_LAT','pk_LON':'Overall_LON'}).reset_index(drop=True)
            combinedOP1 = pd.merge(combinedOP,gdfcop,on=['min_read'])


        ## TO FINDED WEIGHTED LOCATION OF EACH PK GROUP
    #    gdfcop = gdf_pass_pks.loc[:,['OP_NUM','min_read','min_Date','numtimes','verified','pass','pk_LAT','pk_LON','pk_maxCH4_AB']].drop_duplicates()
    #    combinedOP = weightedLoc(gdfcop,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').loc[:,:].rename(columns = {'pk_LAT':'Overall_LAT','pk_LON':'Overall_LON'}).reset_index(drop=True)
    #    combinedOP1 = pd.merge(combinedOP,gdfcop,on=['min_read'])

        ## getting the rest of the stuff
        #gdf_justloc = gdfcop.loc[:,['min_read','pk_LAT','pk_LON','min_Date']].reset_index(drop=True)


        #other = gdfcop.loc[:,['min_read','numtimes','min_Date']].reset_index(drop=True)


        geometry_temp = [Point(xy) for xy in zip(combinedOP1['Overall_LON'], combinedOP1['Overall_LAT'])]
        crs = {'init': 'epsg:4326'}
        gdf_OP = gpd.GeoDataFrame(combinedOP1, crs=crs, geometry=geometry_temp)
        gdf_OP = gdf_OP.to_crs(epsg=32610).copy()

        gdf_OP_reduced = gdf_OP.loc[:,['min_read','geometry','numtimes','Overall_LON','Overall_LAT','min_Date','verified']].drop_duplicates().reset_index(drop=True)
        gdf_OP_reduced.to_file(new_loc_json, driver="GeoJSON")

        #gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']),gdf_pass_pks.drop(columns=['geometry','oldgeo']),on=['min_read','min_Date','numtimes','pass','verified','pk_LAT','pk_LON','OP_NUM','pk_maxCH4_AB'])
        gdf_OP_wrecombine = pd.merge(gdf_OP.drop(columns=['geometry']),gdf_pass_pks.drop(columns=['geometry']),on=['min_read','min_Date','numtimes','pass','verified','pk_LAT','pk_LON','OP_NUM','pk_maxCH4_AB'])

        #gdf_OP.to_csv(new_loc,index=False)
        gdf_OP_wrecombine.to_csv(new_loc,index=False)

            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)


            #geometry is the point of the lat/lon
        #gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

        #gdf_buff = gpd.GeoDataFrame(datFram_wtLocMax, crs=crs, geometry=geometry_temp)



        #gdfcop = gdfcop.to_crs(epsg=32610).copy()
        #gdfcop.to_file(new_loc_json, driver="GeoJSON")


        #gdf_tot.to_csv(new_loc, index = False)
        unique_peaks = gdf_pass_pks.loc[:,['OP_NUM','pk_LAT','pk_LON','min_read','min_Date']].drop_duplicates()
        unique_peaks['save'] = True
        #good_pks = unique_peaks.PEAK_NUM.drop_duplicates().values.tolist()
        good_pks = list(unique_peaks.index)

        def getthing(index):
            if index in good_pks:
                return True
            else:
                return False
        gdf_pass_pks['wooind'] = gdf_pass_pks.index
        gdf_pass_pks['save'] = gdf_pass_pks.swifter.apply(lambda x: getthing(x.wooind),axis=1)

       # unique_pks_tog = pd.concat([unique_peaks, gdf_pass_pks.drop(columns=['LON', 'LAT','PEAK_NUM'])], axis=1, join='inner')
       # testa = pd.merge(gdf_pass_pks, unique_peaks, how='left', on=['PEAK_NUM', 'pk_LAT','pk_LON','min_read'])

        unique_pks_tog = gdf_pass_pks.loc[gdf_pass_pks.save == True,:].reset_index(drop=True)
        unique_pks_tog['Latitude'] = unique_pks_tog.loc[:,'pk_LAT']
        unique_pks_tog['Longitude'] = unique_pks_tog.loc[:,'pk_LON']

        #unique_pks_tog.to_csv(new_loc2, index = False)
        unique_pks_tog.to_csv(new_loc, index = False)

        return()
        #return(gdf_OP_wrecombine)

