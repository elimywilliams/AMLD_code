#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 13:35:17 2020

@author: emilywilliams
"""

thisDir = '/Users/emilywilliams/Documents/DrivingData/ColDat/ProcessedData/'

fileLoc = os.listdir(thisDir)

fileLoc[1]

df = pd.read_csv(thisDir + fileLoc[1])

import pandas as pd, numpy, matplotlib.pyplot as plt
from shapely.geometry import LineString
from time import time

coordinates = df.loc[:,['LAT','LONG']].values

line = LineString(coordinates)
tolerance = .00003
simplified_line = line.simplify(tolerance, preserve_topology = False)

print(line.length, 'line length')
print(simplified_line.length, 'simplified line length')
print(len(line.coords), 'coordinate pairs in full data set')
print(len(simplified_line.coords), 'coordinate pairs in simplified data set')
print(round(((1 - float(len(simplified_line.coords)) / float(len(line.coords))) * 100), 1), 'percent compressed')

lon = pd.Series(pd.Series(simplified_line.coords.xy)[1])
lat = pd.Series(pd.Series(simplified_line.coords.xy)[0])
si = pd.DataFrame({'LONG':lon, 'LAT':lat})
si.tail()


start_time = time()

# df_label column will contain the label of the matching row from the original full data set
si['df_label'] = None

# for each coordinate pair in the simplified set
for si_label, si_row in si.iterrows():    
    si_coords = (si_row['LAT'], si_row['LONG'])
    
    # for each coordinate pair in the original full data set
    for df_label, df_row in df.iterrows():
        
        # compare tuples of coordinates, if the points match, save this row's label as the matching one
        if si_coords == (df_row['LAT'], df_row['LONG']):
            si.loc[si_label, 'df_label'] = df_label
            break
            
print('process took %s seconds' % round(time() - start_time, 2))

# select the rows from the original full data set whose labels appear in the df_label column of the simplified data set
rs = df.loc[si['df_label'].dropna().values]

#rs.to_csv('data/summer-travel-gps-simplified.csv', index=False)
rs.tail()


# plot the final simplified set of coordinate points vs the original full set
plt.figure(figsize=(10, 6), dpi=100)
rs_scatter = plt.scatter(rs['LONG'], rs['LAT'], c='m', alpha=0.3, s=150)
df_scatter = plt.scatter(df['LONG'], df['LAT'], c='k', alpha=0.4, s=10)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Simplified set of coordinate points vs original full set')
plt.legend((rs_scatter, df_scatter), ('Simplified', 'Original'), loc='upper left')
plt.show()

rs.to_csv('/Users/emilywilliams/Documents/DrivingData/shortened.csv', index = False)

