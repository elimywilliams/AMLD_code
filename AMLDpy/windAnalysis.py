#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 14:45:45 2020

@author: emilywilliams
"""
file = '/Users/emilywilliams/SouthernCross/windData.csv'
aWind = pd.read_csv(file)
from numpy import pi
import numpy as np
def calcVel(timediff,distance):
    if timediff == 0:
        return(0)
    elif timediff != 0:
        return(distance/timediff)
        
radians = False
wind_df = aWind.copy()   
wind_df['U'] = wind_df.loc[:,'SONIC_U']
wind_df['V'] = wind_df.loc[:,'SONIC_V']
wind_df['W'] = wind_df.loc[:,'SONIC_W']

wind_df['LAT'] = wind_df.loc[:,'GPS_LAT']
wind_df['LONG'] = wind_df.loc[:,'GPS_LONG']



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
#wind_df['BCH4']= wind_df.loc[:,['shift_CH4']]
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


import windrose

wind_df4.to_csv('/Users/emilywilliams/SouthernCross/windDataAnalysed.csv')
