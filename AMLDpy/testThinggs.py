#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:28:18 2020

@author: emilywilliams
"""
mainLoc = '/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/FilteredPeaks_SCCar_20200401.csv'
secLoc = '/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/FilteredPeaks_SCCar_20200410.csv'

mainLoc = '/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/FilteredPeaks_SCCar_20200401.csv'
secLoc = '/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/FilteredPeaks_SCCar_20200417.csv'


firstgroup = pd.read_csv(mainLoc)
secondgroup = pd.read_csv(secLoc)

passCombine(firstgroup,secondgroup)

toCombine = os.listdir('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/')
toCombineList = []
index = 0
for file in toCombine:
    if file.startswith('Filtered') and file.endswith('.csv') and not file.endswith('info.csv'):
        index = index +1
        if index == 1:
            mainThing = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file)
        elif index != 1:
            mainThing = passCombine(mainThing, pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file))


index = 1
for file in toCombineList:
    if index == 1:
        print(file)
        mainThing = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file)
    elif index > 1:
        secThing = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file)
        mainThing = passCombine(mainThing,secThing)
    index = index+1   
    print(index)
        
file1 = 'FilteredPeaks_SCCar_20200410.csv'
file2 = 'FilteredPeaks_SCCar_20200407.csv'
file3 = 'FilteredPeaks_SCCar_20200406.csv'
file4 = 'FilteredPeaks_SCCar_20200417.csv'
file5 = 'FilteredPeaks_SCCar_20200401.csv'
file6= 'FilteredPeaks_SCCar_20200415.csv'
file7 = 'FilteredPeaks_SCCar_20200408.csv'

mainThing = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file1)
secThing = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file2)

woo = passCombine(mainThing,secThing)
mainThing = woo.copy()

secThing = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file3)
woo = passCombine(mainThing,secThing)
mainThing = woo.copy()

secThing = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file4)
woo = passCombine(mainThing,secThing)
mainThing = woo.copy()

secThing = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file4)
woo = passCombine(mainThing,secThing)
mainThing = woo.copy()

secThing = pd.read_csv('/Users/emilywilliams/Documents/DrivingData/422_new/FilteredObservedPeaks/' + file5)
woo = passCombine(mainThing,secThing)
mainThing = woo.copy()

mainLoc = weightedLoc(mainThing,'pk_LAT','pk_LON','min_read','pk_maxCH4_AB').rename(columns={'pk_LON':'Overall_LON','pk_LAT':'Overall_LAT'})

mainTog = pd.merge(mainLoc,mainThing,on=['min_read'])
geometry_temp = [Point(xy) for xy in zip(mainTog['Overall_LON'], mainTog['Overall_LAT'])]
crs = {'init': 'epsg:4326'}
mainGPD = gpd.GeoDataFrame(mainTog,crs = crs,geometry=geometry_temp)
mainGPD = mainGPD.to_crs(epsg = 4326)

woot = mainGPD.loc[:,['min_read','Overall_LON','Overall_LAT','numtimes','verified','geometry']].drop_duplicates().reset_index(drop = True)

woot.to_file('/Users/emilywilliams/Documents/DrivingData/woo.json',driver="GeoJSON")



passCombine(firstgroup,secondgroup)

first_geo = [Point(xy) for xy in zip(firstgroup['pk_LON'], firstgroup['pk_LAT'])]
crs = {'init': 'epsg:4326'}
firstgrp = gpd.GeoDataFrame(firstgroup.drop(columns=['geometry','pkGEO']),crs = crs,geometry = first_geo)


first_buff = firstgrp.copy()

#first_buff['geometry2'] = first_buff.loc[:,'geometry'].buffer(30) 

first_buff['geometry2'] = first_buff.apply(lambda x: x.geometry.buffer(0.0001*3),axis =1)

first_buff = first_buff.drop(columns = ['geometry'])
first_buff['geometry'] = first_buff.loc[:,'geometry2']
first_buff = first_buff.drop(columns = ['geometry2'])
#first_buff = first_buff.to_crs(epsg=32610)
first_buff.plot()

firstgrp = first_buff.copy()


sec_geo = [Point(xy) for xy in zip(secondgroup['pk_LON'], secondgroup['pk_LAT'])]
crs = {'init': 'epsg:4326'}
secgrp = gpd.GeoDataFrame(secondgroup.drop(columns=['geometry','pkGEO']),crs = crs,geometry = sec_geo)


sec_buff = secgrp.copy()

#first_buff['geometry2'] = first_buff.loc[:,'geometry'].buffer(30) 

sec_buff['geometry2'] = sec_buff.apply(lambda x: x.geometry.buffer(0.0001*3),axis =1)

sec_buff = sec_buff.drop(columns = ['geometry'])
sec_buff['geometry'] = sec_buff.loc[:,'geometry2']
sec_buff = sec_buff.drop(columns = ['geometry2'])
#first_buff = first_buff.to_crs(epsg=32610)
sec_buff.plot()

secgrp = sec_buff.copy()


##############



from shapely.geometry import Point # Shapely for converting latitude/longtitude to geometry
geometry_temp = [Point(xy) for xy in zip(together['pk_LON'], together['pk_LAT'])]
crs = {'init': 'epsg:4326'}
firstgroup.columns

import utm
utm.from_latlon(51.2, 7.5)


import pyproj
import shapely.geometry as shpgeo

proj = pyproj.Proj(proj="utm", zone=13, ellps="WGS84", datum="WGS84")

def toFromUTM(shp, proj, inv=False):
    geoInterface = shp.__geo_interface__

    shpType = geoInterface['type']
    coords = geoInterface['coordinates']
    if shpType == 'Polygon':
        newCoord = [[proj(*point, inverse=inv) for point in linring] for linring in coords]
    elif shpType == 'MultiPolygon':
        newCoord = [[[proj(*point, inverse=inv) for point in linring] for linring in poly] for poly in coords]

    return shpgeo.shape({'type': shpType, 'coordinates': tuple(newCoord)})

init_shape_utm = toFromUTM(firstgrp, proj)
buffer_shape_utm = init_shape_utm.buffer(30.48)
buffer_shape_lonlat = toFromUTM(buffer_shape_utm, proj, inv=True)



 if data_overlap.size > 0: 
                    
                    firstnull2 = data_overlap.loc[data_overlap.min_read_1.isnull(),:]
                    firstnull = firstnull2.copy()
                    firstnull.loc[:,'min_read_1'] = firstnull2.loc[:,'min_read_2']
                    
                    secnull2 = data_overlap.loc[data_overlap.min_read_2.isnull(),:]
                    
                    secnull = secnull2.copy()
                    secnull.loc[:,'min_read_2'] = secnull2.loc[:,'min_read_1']
                    
                    withoutNA = data_overlap.copy().dropna()
                    allTog2 = pd.concat([firstnull,secnull,withoutNA]).reset_index().copy()
                    
                    
                    allTog2['notsame'] = allTog2.apply(lambda x:x.min_read_1 == x.min_read_2,axis=1)
                    allTog = allTog2.loc[allTog2.notsame == False,:].drop(columns = ['notsame'])
                    
                    
                    over = allTog.copy()
                    over['sorted']=over.apply(lambda y: sorted([y['OP_NUM_1'],y['OP_NUM_2']]),axis=1)
                    over['sorted']=over.sorted.apply(lambda y: ''.join(y))
                    over = over.drop_duplicates('sorted')
                    over['combined']= [list(x) for x in list(over.loc[:,['OP_NUM_1','OP_NUM_2']].to_numpy())]
                    over['date1'] = over.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM_1'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                    over['date2'] = over.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM_2'][6:-2])).strftime('%Y-%m-%d'),axis=1)
    
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
                    combined['dates'] = combined.apply(lambda x: datetime.fromtimestamp(int(x['OP_NUM'][6:-2])).strftime('%Y-%m-%d'),axis=1)
                    combined['pk_Dates'] = [list(x) for x in list(combined.loc[:,['dates']].to_numpy())]
    
                    combined['min_Date'] = combined.loc[:,'dates']
    
                    combined['numtimes'] = 1
                    combined['newgeo'] = combined.loc[:,'geometry']
                    combined['min_read'] = combined.loc[:,"OP_NUM"]
                    
                    for index,row in combined.iterrows():
                        for index2,row2 in newOverlap.iterrows():
                            if row.OP_NUM in row2.recombine:
                                combined.at[index, 'recombine']  = row2.recombine.copy()
                                combined.at[index, 'newgeo']  = row2.copy().geometry
                                combined.at[index,'min_read'] = row2.copy().min_read
                                combined.at[index, 'pk_Dates']  = row2.pk_Dates
                                combined.at[index,'min_Date'] = row2.min_Date
                                
                    combined['numtimes'] = combined.apply(lambda y: len(y.recombine),axis = 1).copy()
                    combined['numdays'] = combined.apply(lambda y: len(y.pk_Dates),axis = 1).copy()
                    combined_reduced = combined.loc[:,['OP_NUM','newgeo','recombine','numtimes','min_read','numdays','pk_Dates','min_Date']]
                    gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['OP_NUM']).copy()
                    gdf_pass_pks['verified'] = gdf_pass_pks.apply(lambda y: (True if y.numtimes > 1 else False),axis=1 ).copy()


import ast
x = woo
x = ast.literal_eval(x)
x = [n.strip() for n in x]

def strList(x):
    x = ast.literal_eval(x)
    x = [n.strip() for n in x]
    return(x)

uniqueMin1 = overcop.min_read_1.unique()

  for index,row in overcop.iterrows():
                first_thing = first_dis[first_dis['min_read']== row.min_read_1].loc[:,['recombine']]
                firstcomb = first_thing.recombine.explode().copy()
                first_list = firstcomb.reset_index().recombine.copy()
                        
               ## fix this 
               
                sec_thing = sec_dis[sec_dis['min_read']== row.min_read_1].loc[:,['recombine']]
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
            combined_reduced = combined.loc[:,['min_read','geometry','oldgeo','recombine','numtimes','prev_read']]
            #gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['min_read'])
            gdf_tog = gdf_tog.drop('min_read',axis=1)
            gdf_tog = gdf_tog.drop('numtimes',axis=1)
            gdf_tog = gdf_tog.drop('recombine',axis=1)
            
            #gdf_tog = gdf_tog.drop(columns = ['min_read','numtimes','recombine'])
            
            # gdf_tog = gdf_tog.drop('newgeo',axis=1)
           # gdf_tog['firstgeo'] = gdf_tog.copy().oldgeo
           # gdf_tog['secondgeo'] = gdf_tog.copy().geometry
            
           # del(gdf_tog['geometry'])
            #gdf_tog.drop(columns=['geometry'])
           # del(gdf_tog['oldgeo'])

                                   
            gdf_pass_pks = pd.merge(gdf_tog,combined_reduced,on = ['prev_read']).copy()

            gdf_pass_pks['verified'] = gdf_pass_pks.copy().apply(lambda y: (True if y.numtimes > 1 else False),axis=1 )
        if data_overlap.size == 0:
            gdf_pass_pks = gdf_tog.copy()
            
    ## didnt adress if the bind shape was only size only 1
        gdf_tot_pks = pd.merge(gdf_pass_pks,tot_pks,on = ['OP_NUM','min_read']).copy()
        gdf_tot_pks['numtimes']= gdf_tot_pks.apply(lambda x: len(x.recombine.split(',')),axis=1)
    #return(gdf_pass_pks)
    return(gdf_tot_pks)


def findOriginal(rd):
    row = gdf_tog.loc[gdf_tog.min_read == rd,'recombine']
    
    if row.shape[0] == 0:
        return ([])
    elif row.shape[0] != 0:
        row2 = row.drop_duplicates().reset_index(drop=True)[0]
        return(row2)
    
def getBig(recom):
    if type(recom) == list:
    #bigList = strList(recom)
        bigList = recom
        howmany = len(bigList)
        allRecombine = []
    elif type(recom) != list:
        bigList = strList(recom)
        howmany=len(bigList)
    
    allRecombine = []
    for x in bigList:
        woo = findOriginal(x)
        if len(woo) != 0:
            woo2 = strList(findOriginal(x))
            allRecombine = allRecombine + woo2
            #woo = findOriginal(x)
            #allRecombine = allRecombine.append(woo)
    return(allRecombine)
        
new_1 = over.copy().drop(columns= ['min_read_2']).rename(columns={'min_read_1':'min_read'})
getBig(new_1.bothcombine[0])

new_2 = over.copy().drop(columns= ['min_read_1']).rename(columns={'min_read_2':'min_read'})
getBig(new_2.bothcombine[0])

newtog = pd.concat([new_1,new_2])
minreads = newtog.min_read.unique().tolist()


toChange = gdf_tog.loc[(gdf_tog['min_read'] in minreads),: ]

toChange = gdf_tog[gdf_tog['min_read'].isin(minreads)].drop(columns=['geometry'])

toChangecombined = pd.merge(toChange,newtog,on=['min_read']).drop(columns = ['geometry','numtimes','verified','recombine'])
toChangecombined = toChangecombined.rename(columns = {'min_read':'prev_read','bothcombine':'recombine'})
toChangecombined['numtimes'] = toChangecombined.apply(lambda x: len(x.recombine),axis = 1)
toChangecombined['verified'] = toChangecombined.apply(lambda x:x.numtimes>1,axis=1 )
toChangecombined = toChangecombined.rename(columns ={'min_val':'min_read'}).drop(columns = ['combined'])

toNotChange = gdf_tog[~gdf_tog['min_read'].isin(minreads)].drop(columns = ['geometry'])
toNotChange.loc[:,'prev_read'] = gdf_tog[~gdf_tog['min_read'].isin(minreads)].min_read

newCombined = pd.concat([toChangecombined,toNotChange])
