import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

### first thing
relpolys = gpd.read_file('/Users/emilywilliams/Documents/GitHub/Trussville/rel_poly_1983.shp')
techInds = gpd.read_file('/Users/emilywilliams/Documents/GitHub/Trussville/techIndications.shp')
aerisInds = gpd.read_file('/Users/emilywilliams/Documents/GitHub/Trussville/aeris678lks.shp')

totPolys = relpolys.geometry.unary_union
techInds['win'] = techInds.within(totPolys)
techIndsWin = techInds.loc[techInds.win == True,:].reset_index()

aerisInds['within'] = aerisInds.within(totPolys)
aerisIndsWin = aerisInds.loc[aerisInds['within'] == True,:].reset_index()
del(aerisInds)
aerisInds = aerisIndsWin.copy()
del(aerisIndsWin)
#aerisInds = aerisInds.to_crs(epsg=26730)
### add buffer
buffdist = 100
#aeris_buff100 = aerisInds.to_crs(epsg=26730).buffer(0.00001*buffdist)
aeris_buff100 = aerisInds.buffer(0.00001*buffdist)
aeris_buff100_together = aeris_buff100.geometry.unary_union
techwinAeris = techIndsWin.copy()
techwinAeris['win'] = techIndsWin.within(aeris_buff100_together)
techwinAeris = techwinAeris.loc[techwinAeris.win == True,:]

techwinAeris_SC = techwinAeris.copy()
firstThing = scIndsWin.loc[scIndsWin.analysis_num == 0,:]
firstThingbuff = firstThing.buffer(0.00001*buffdist)
firstThingbuff_together = firstThingbuff.geometry.unary_union
techwinAeris_SC['win'] = techwinAeris_SC.within(firstThingbuff_together)
techwinAeris_SC_win = techwinAeris_SC.loc[techwinAeris_SC.win == True,:]


## first


ax = relpolys.plot()
techIndsWin.plot(ax=ax,color = 'black')

## importing the different files
import os
filesLoc = '/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data'
toCheck = os.listdir(filesLoc)
toCheck2 = []
for index,name in enumerate(toCheck):
    if name.startswith('truss_') and 'shift' and 'remove' in name: # and '6pc_102_med_shift_0' not in name:
        toCheck2.append([name])

for index,name in enumerate(toCheck2):
    namesplit = str(name[0]).split('_')
    tperc = float(namesplit[1][:-2])
    tbaselineobs = float(namesplit[2])
    tbackgroundperc = str(namesplit[3])
    tshift = int(namesplit[5])
    tRmin = str(namesplit[6][7:])
    if 'n1' in tRmin:
        tRmin = -1
    elif 'n1' not in tRmin:
        tRmin = int(tRmin)
    tminElevated = int(namesplit[7][0])
    temp_df = gpd.read_file(filesLoc + '/' + str(name[0]) + '/FinalShpFiles/OP_Final.json')
    temp_df.loc[:,'baseline_obs'] = temp_df.apply(lambda x: tbaselineobs,axis = 1)
    temp_df.loc[:,'baseline_perc'] = temp_df.apply(lambda x: tbackgroundperc,axis = 1)
    temp_df.loc[:,'CH4_shift'] = temp_df.apply(lambda x: tshift,axis = 1)
    temp_df.loc[:,'min_Rval'] = temp_df.apply(lambda x: tRmin,axis = 1)
    temp_df.loc[:,'min_elevated'] = temp_df.apply(lambda x: tminElevated,axis = 1)
    temp_df.loc[:,'analysis_num'] = temp_df.apply(lambda x: index,axis = 1)

    if index == 0:
        final_df = temp_df.copy()
    elif index != 0:
        final_df = pd.concat([final_df,temp_df])

scInds = final_df.copy()
scInds['win'] = scInds.within(totPolys)
scIndsWintemp = scInds.loc[scInds.win == True,:].reset_index()
scIndsWin = scIndsWintemp.join(pd.DataFrame(scIndsWintemp.groupby('analysis_num').win.count()).rename(columns = {'win':'indications_win'}),on = ['analysis_num'])

for index,anum in enumerate(scIndsWin.analysis_num.unique()):
    smallDF = scIndsWin.loc[scIndsWin.analysis_num == anum,: ].reset_index()
    smallBuff = smallDF.buffer(0.00001 * buffdist)
    smallBuff_together = smallBuff.geometry.unary_union
    tempTechWinAerisSC = techwinAeris.copy()
    tempTechWinAerisSC['win'] = tempTechWinAerisSC.within(smallBuff_together)
    tempTechWinAerisSC_win = tempTechWinAerisSC.loc[tempTechWinAerisSC.win == True, :]
    numtechwin =tempTechWinAerisSC_win.shape[0]
    numscwin = smallDF.shape[0]
    backobs = smallDF.baseline_obs[0]
    baseperc = smallDF.baseline_perc[0]
    num_shift = smallDF.CH4_shift[0]
    threshold = smallDF.threshold[0]
    minR = smallDF.min_Rval[0]
    minel = smallDF.min_elevated[0]
    ntowa = techwinAeris.shape[0]
    tempData = {'AnalysisNumber': [anum], 'NumberTechObsWinAerisSC':[numtechwin],
                'NumberSCWinPoly':[numscwin],'BackgroundObsNum':[backobs],
                'BaselinePercentile':[baseperc],'CH4Shift':[num_shift],'MinR':[minR],'minElevated':[minel],
                'ElevationThreshold':[threshold],'NumberTechObsWinAeris':[ntowa],
                'PercSC_tot':[numtechwin/ntowa],"Buffer":[buffdist]}

    tempDF = pd.DataFrame.from_dict(tempData)
    if index == 0:
        sumDF = tempDF.copy()
    elif index != 0:
        sumDF = pd.concat([sumDF,tempDF])

sumDF.to_csv('SummaryInfoTrussville.csv')


firstThing = scIndsWin.loc[scIndsWin.analysis_num == 0,:]
firstThingbuff = firstThing.buffer(100)
firstThingbuff_together = firstThingbuff.geometry.unary_union
techwinAeris_SC = techwinAeris.copy()
techwinAeris_SC['win'] = techwinAeris_SC.within(firstThingbuff_together)
techwinAeris_SC = techwinAeris.loc[techwinAeris.win == True,:]

firstThingbuff.to_file('/Users/emilywilliams/Documents/GitHub/Trussville/testFirst.json',driver="GeoJSON")
aeris_buff100.to_file('/Users/emilywilliams/Documents/GitHub/Trussville/testAerisBuff.json',driver="GeoJSON")

techwinAeris.to_file('/Users/emilywilliams/Documents/GitHub/Trussville/techwinaeris.json',driver="GeoJSON")
techwinAeris.to_file('/Users/emilywilliams/Documents/GitHub/Trussville/techwinaeris.shp')
techwinAeris_SC_win.to_file('/Users/emilywilliams/Documents/GitHub/Trussville/techwinaeris_sc.shp')


techIndsWin.to_csv('techIndicationsLocs.csv')
aerisInds.to_csv('aerisIndicationsLocs.csv')