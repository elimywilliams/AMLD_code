import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

testDat = pd.read_csv('/Users/emilywilliams/Downloads/colinSample.csv')

testDat = pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/testAggregate/AggregatedData_1.csv')
testDat['AerisWaitTime'] = testDat.apply(lambda x: 0.5,axis =1)
testDat['ActisenseWaitTime'] = testDat.apply(lambda x: 0.05,axis =1)
testDat['prevTime'] = testDat.timestamp.shift(1)
testDat['timeDiff'] = testDat.apply(lambda x: x.timestamp - x.prevTime,axis = 1)

def analysePD(testDat,awt,aswt):
    testDat['AerisWaitTime'] = testDat.apply(lambda x: str(awt), axis=1)
    testDat['ActisenseWaitTime'] = testDat.apply(lambda x: str(aswt), axis=1)
    testDat['prevTime'] = testDat.timestamp.shift(1)
    testDat['timeDiff'] = testDat.apply(lambda x: x.timestamp - x.prevTime, axis=1)
    return(testDat)



agg1 = analysePD(pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/testAggregate/AggregatedData_1.csv'),.5,.05)
agg2 = analysePD(pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/testAggregate/AggregatedData_2.csv'),.5,.01)
agg3 = analysePD(pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/testAggregate/AggregatedData_3.csv'),.25,.01)
agg4 = analysePD(pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/testAggregate/AggregatedData_4.csv'),.1,.01)
agg5 = analysePD(pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/testAggregate/AggregatedData_5.csv'),.01,.01)
agg5x = analysePD(pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/AggregatedData_x.csv'),.01,.01)

testDat.timeDiff.max()
testDat.timeDiff.min()
testDat.timeDiff.var().sqrt()

np.histogram(testDat.timeDiff.dropna())

n,bins,patches = plt.hist(x=testDat.timeDiff.dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('Time Difference Between Observations')
plt.ylabel('Frequency')
plt.title('Distribution of Time Differences Between Aggregator')
maxfreq = n.max()
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
plt.xlim(xmin = 0.05,xmax = 0.15)

def histDif(pdf):
    testDat = pdf.copy()
    awt = float(testDat['AerisWaitTime'][0])
    aswt = float(testDat['ActisenseWaitTime'][0])
    n,bins,patches = plt.hist(x=testDat.timeDiff.dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
    plt.grid(axis = 'y',alpha = 0.75)
    plt.xlabel('Time Difference Between Observations')
    plt.ylabel('Frequency')
    plt.title('Distribution of Time Differences Between Observations,\n Aeris WT: ' + str(awt) + ', Actisense WT: ' + str(aswt))
    maxfreq = n.max()
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.xlim(xmin = 0.00,xmax = .2)

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(3, 3)
ax1 = plt.subplot(gs[0, :])
ax2 = plt.subplot(gs[1,:-1])
ax3 = plt.subplot(gs[1:, -1])
ax4 = plt.subplot(gs[-1,0])
ax5 = plt.subplot(gs[-1,-2])

histDif(agg1)
histDif(agg2)
histDif(agg3)
histDif(agg4)
histDif(agg5)
histDif(agg5x)

fig, axs = plt.subplots(2, 2)
if axs[0, 0]:



axs[0, 0].plot(x, y)
axs[0, 0].set_title('Axis [0, 0]')
axs[0, 1].plot(x, y, 'tab:orange')
axs[0, 1].set_title('Axis [0, 1]')
axs[1, 0].plot(x, -y, 'tab:green')
axs[1, 0].set_title('Axis [1, 0]')
axs[1, 1].plot(x, -y, 'tab:red')
axs[1, 1].set_title('Axis [1, 1]')



smallDat = testDat.loc[:,['timeDiff','channel']].dropna()
smallDat.boxplot(by = 'channel',column = ['timeDiff'],grid = False)

mn = testDat.timeDiff.dropna().median()
testMn = smallDat.copy()
testMn['absDiff'] = smallDat.apply(lambda x: abs(x.timeDiff - mn), axis = 1)

n,bins,patches = plt.hist(x=testMn.absDiff,bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('Time Difference Between Observations')
plt.ylabel('Frequency')
plt.title('Distribution of abs(x - mn(x)) Differences Between Aggregator')
maxfreq = n.max()
plt.xlim(xmin = 0, xmax = 0.06)

def notNew(df):
    testDat = df.copy()
    testDat['notnewCH4'] = testDat.CH4.eq(testDat.CH4.shift(1))
    testDat['newCH4'] = testDat.apply(lambda x: not x.notnewCH4,axis =1)
    testDat = testDat.drop(columns = ['notnewCH4'])

    testDat['notnewH20'] = testDat.H2O.eq(testDat.H2O.shift(1))
    testDat['newH20'] = testDat.apply(lambda x: not x.notnewH20,axis =1)
    testDat = testDat.drop(columns = ['notnewH20'])
    testDat['newAeris'] = testDat.apply(lambda x: x.newH20 or x.newCH4,axis = 1)


    testDat['notnewLAT'] = testDat.latitude.eq(testDat.latitude.shift(1))
    testDat['newLAT'] = testDat.apply(lambda x: not x.notnewLAT,axis =1)
    testDat = testDat.drop(columns = ['notnewLAT'])

    testDat['notnewLON'] = testDat.longitude.eq(testDat.longitude.shift(1))
    testDat['newLON'] = testDat.apply(lambda x: not x.notnewLON,axis =1)
    testDat = testDat.drop(columns = ['notnewLON'])

    testDat['newGPS'] = testDat.apply(lambda x: x.newLON or x.newLAT,axis = 1)

    ### WIND
    testDat['notnewWS'] = testDat.windSpeed.eq(testDat.windSpeed.shift(1))
    testDat['newWS'] = testDat.apply(lambda x: not x.notnewWS,axis =1)
    testDat = testDat.drop(columns = ['notnewWS'])

    testDat['notnewWA'] = testDat.windAngle.eq(testDat.windAngle.shift(1))
    testDat['newWA'] = testDat.apply(lambda x: not x.notnewWA,axis =1)
    testDat = testDat.drop(columns = ['notnewWA'])

    testDat['notnewWR'] = testDat.windReference.eq(testDat.windReference.shift(1))
    testDat['newWR'] = testDat.apply(lambda x: not x.notnewWR,axis =1)
    testDat = testDat.drop(columns = ['notnewWR'])

    testDat['newWind'] = testDat.apply(lambda x: x.newWA or x.newWR or x.newWS,axis = 1)
    return (testDat)

notNewDat1 = notNew(agg1)
notNewDat2 = notNew(agg2)
notNewDat3 = notNew(agg3)
notNewDat4 = notNew(agg4)
notNewDat5 = notNew(agg5)
notNewDat5x = notNew(agg5x)


testDat.drop(columns = ['newCH4','newH20','newLAT','newLON','newWS','newWA','newWR']).to_csv('checkShift.csv')

def histWind(df):
    newWindDf = df.loc[df.newWind == True, :]
    newWindDf.loc[:, ['wind_timediff']] = newWindDf.timestamp.sub(newWindDf.timestamp.shift(1))

    awt = float(newWindDf['AerisWaitTime'][0])
    aswt = float(newWindDf['ActisenseWaitTime'][0])

    n, bins, patches = plt.hist(x=newWindDf.wind_timediff.dropna(), bins='auto', color='#0504aa', alpha=0.7,
                                rwidth=0.85, density=True)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Time Difference Between Observations (wind)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Time Differences Between Wind Observations,\n Aeris WT: ' + str(awt) + ', Actisense WT: ' + str(aswt))
    maxfreq = n.max()
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.xlim(xmin=0.05, xmax=.4)

def histCH4(df,obs):
    df['obsnum'] = df.index
    awt = float(df['AerisWaitTime'][0])
    aswt = float(df['ActisenseWaitTime'][0])
    df = df.loc[df.obsnum >= 60,:]
    newWindDf = df.loc[df.newCH4 == True, :]
    newWindDf.loc[:, ['wind_timediff']] = newWindDf.timestamp.sub(newWindDf.timestamp.shift(1))
    newWindDf.loc[:, ['wind_obsdiff']] = newWindDf.obsnum.sub(newWindDf.obsnum.shift(1))
    newWindDf = newWindDf.loc[newWindDf.obsnum >= 30,:]

    if obs == True:
        n, bins, patches = plt.hist(x=newWindDf.wind_obsdiff.dropna(), bins='auto', color='#0504aa', alpha=0.7,
                                    rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('Number of repeated observations before new reading (CH4)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Number of Readings between new CH4,\n Aeris WT: ' + str(awt) + ', Actisense WT: ' + str(aswt))
        maxfreq = n.max()
        #plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
        #plt.xlim(xmin=0, xmax=16)
    elif obs == False:
        n, bins, patches = plt.hist(x=newWindDf.wind_timediff.dropna(), bins='auto', color='#0504aa', alpha=0.7,
                                    rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('Time Difference Between Observations (CH4)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Time Differences Between CH4 Observations,\n Aeris WT: ' + str(awt) + ', Actisense WT: ' + str(aswt))
        maxfreq = n.max()
        plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
        #plt.xlim(xmin=0.05, xmax=.4)
def histGPS(df,obs):
    df['obsnum'] = df.index
    awt = float(df['AerisWaitTime'][0])
    aswt = float(df['ActisenseWaitTime'][0])
    df = df.loc[df.obsnum >= 60,:]
    newWindDf = df.loc[df.newGPS == True, :]
    newWindDf.loc[:, ['wind_timediff']] = newWindDf.timestamp.sub(newWindDf.timestamp.shift(1))
    newWindDf.loc[:, ['wind_obsdiff']] = newWindDf.obsnum.sub(newWindDf.obsnum.shift(1))
    newWindDf = newWindDf.loc[newWindDf.obsnum >= 30,:]

    if obs == True:
        n, bins, patches = plt.hist(x=newWindDf.wind_obsdiff.dropna(), bins='auto', color='#0504aa', alpha=0.7,
                                    rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('Number of repeated observations before new reading (GPS)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Number of Readings between new GPS,\n Aeris WT: ' + str(awt) + ', Actisense WT: ' + str(aswt))
        maxfreq = n.max()
        #plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
        #plt.xlim(xmin=0, xmax=16)
    elif obs == False:
        n, bins, patches = plt.hist(x=newWindDf.wind_timediff.dropna(), bins='auto', color='#0504aa', alpha=0.7,
                                    rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('Time Difference Between Observations (GPS)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Time Differences Between GPS Observations,\n Aeris WT: ' + str(awt) + ', Actisense WT: ' + str(aswt))
        maxfreq = n.max()
        plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
        #plt.xlim(xmin=0.05, xmax=.4)

histWind(notNewDat1)
histWind(notNewDat2)
histWind(notNewDat3)
histWind(notNewDat4)
histWind(notNewDat5)
histWind(notNewDat5x)

histCH4(notNewDat1,True)
histCH4(notNewDat2,True)
histCH4(notNewDat3,True)
histCH4(notNewDat4,True)
histCH4(notNewDat5,True)
histCH4(notNewDat5x,True)

histGPS(notNewDat1,True)
histGPS(notNewDat2,True)
histGPS(notNewDat3,True)
histGPS(notNewDat4,True)
histGPS(notNewDat5,True)
histGPS(notNewDat5x,True)

###
df = pd.read_csv('checkShift.csv')
df = testDat.copy()
newWindDf = df.loc[df.newWind == True,:]
newWindDf.loc[:,['wind_timediff']] = newWindDf.timestamp.sub(newWindDf.timestamp.shift(1))

n,bins,patches = plt.hist(x=newWindDf.wind_timediff.dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('Time Difference Between Observations (wind)')
plt.ylabel('Frequency')
plt.title('Distribution of Time Differences Between new wind observations')
maxfreq = n.max()
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
plt.xlim(xmin = 0.05,xmax = 0.15)

newCH4DF = df.loc[df.newCH4 == True,:]
newCH4DF.loc[:,['ch4_timediff']] = newCH4DF.timestamp.sub(newCH4DF.timestamp.shift(1))

n,bins,patches = plt.hist(x=newCH4DF.ch4_timediff.dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('Time Difference Between Observations of New Aeris readings')
plt.ylabel('Frequency')
plt.title('Distribution of Time Differences Between New Aeris Readings')
maxfreq = n.max()
plt.xlim(xmin = 0.85,xmax = 1.2)

newGPSDF = df.loc[df.newGPS == True,:]
newGPSDF.loc[:,['gps_timediff']] = newGPSDF.timestamp.sub(newGPSDF.timestamp.shift(1))

n,bins,patches = plt.hist(x=newGPSDF.gps_timediff.dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('Time diff Between Observations of New GPS readings')
plt.ylabel('Frequency')
plt.title('Distribution of Time Differences Between New GPS Readings')
maxfreq = n.max()
plt.xlim(xmin = 0.05,xmax = .15)

leaksinitial = pd.read_csv('https://raw.githubusercontent.com/elimywilliams/Trussville/master/allLeaksWin.csv')
leaks= leaksinitial.loc[leaksinitial.loc[:,'Time Stamp'] < '09/01/2020 00:00:01']

n,bins,patches = plt.hist(x=leaks.ch4Enhcmt.dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('CH4 Enhancement')
plt.ylabel('Frequency')
plt.title('Distribution of CH4 Enhancement, indications from Aeris (Trussville)')
maxfreq = n.max()
plt.xlim(0,3)
#plt.xlim(xmin = min(leaks.ch4Enhcmt),xmax = max(leaks.ch4Enhcmt))

min(leaks.ch4Enhcmt)
max(leaks.ch4Enhcmt)

n,bins,patches = plt.hist(x=leaks.loc[:,'CH4 (ppm)'].dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('CH4')
plt.ylabel('Frequency')
plt.title('Distribution of CH4 raw ppm, marked as elevated, Aeris (Trussville)')
maxfreq = n.max()
plt.xlim(xmin = 1.8,xmax = 7)

max(leaks.loc[:,'CH4 (ppm)'])

min(leaks.loc[:,'ch4Enhcmt'])
max(leaks.loc[:,'ch4Enhcmt'])

leaks.loc[leaks.ch4Enhcmt == 0,:]
totPD.loc[totPD.raw_CH4 == 2.01460,807]
totPD.loc[totPD.raw_CH4 == 2.00437,:]



import os
list = os.listdir('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/allTruss/ProcessedData')
listuse = []
for index,x in enumerate(list):
    if x.endswith('dat.csv'):
        listuse.append(x)

for index,x in enumerate(listuse):
    tempPD = pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/allTruss/ProcessedData/' + str(x))
    if index == 0:
        totPD = tempPD.copy()
    elif index != 0:
        totPD = pd.concat([totPD,tempPD])

### LOOKING AT ALL CH4

n,bins,patches = plt.hist(x=totPD.loc[:,'raw_CH4'].dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('CH4')
plt.ylabel('Frequency')
plt.title('Distribution of CH4 raw ppm, of all Aeris (Trussville)')
maxfreq = n.max()
plt.xlim(xmin = 1.8,xmax = 2.6)

bigch4 = totPD.loc[totPD.raw_CH4 > 1.8,:]

bigch4.shape[0]/totPD.shape[0]

###
SC_Leaks = pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/allTruss/FinalShpFiles/overallPeaks2_wlog.csv')
n,bins,patches = plt.hist(x=SC_Leaks.loc[:,'ch4_enh'].dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('CH4')
plt.ylabel('Frequency')
plt.title('Distribution of CH4 excess, of all SC Leaks indicated')
maxfreq = n.max()
plt.xlim(xmin = 0,xmax = 3)

smallPD = totPD.loc[totPD.CH4 > 2.01,:]

mainSmall = pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/truss_7pc_1_med_102_ns/FinalShpFiles/mainThing_small.csv')
n,bins,patches = plt.hist(x=mainSmall.loc[:,'OB_CH4'].dropna(),bins = 'auto',color = '#0504aa',alpha = 0.7,rwidth = 0.85,density=True )
plt.grid(axis = 'y',alpha = 0.75)
plt.xlabel('Observed CH4')
plt.ylabel('Frequency')
plt.title('Distribution of Observed CH4, of all SC Leaks indicated')
maxfreq = n.max()
plt.xlim(xmin = 1.8,xmax = 4)