#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday July 28 
@author: emilywilliams
"""

############## ALL THAT NEEDS TO BE CHANGED ##############################################

## WHERE THE amld_Functions.py file is located
function_file_Loc = '/Users/emilywilliams/Documents/GitHub/AMLD_CODE/AMLDpy/'

## Folder with .txt Data
raw_data_loc = "/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/trussville_dat"

## Folder to put results in (will make subfolders later)
results_folder_loc = "/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/trussville_dat/"

car_id = 'TrussCar' #CAR NAME TO APPEAR IN FILENAMES OBSERVED PEAK NAMES
threshold = '0.05'  #What Proportion above Baseline to flag as elevated (i.e. 0.1 = 10% higher)
time_thresh = '5.0'  ## How many minutes to include in background calculation (minutes)
initial_time_ignore = '0' ## How many minutes to skip at the beginning of the dataset (i.e. if Collin is at his house)
min_elevated = '1' # minimum number of elevated readings required for an observed peak
shift = 0  ## Lag time for CH4 to reach sensor (in seconds)
engineering = False #is this an engineering file
aeris = True # is this from the aeris instrument
CSU = True
time_push = 0 #not sure what this is
back_obs_num = '102' ### NUMBER OF OBSERVATIONS TO INCLUDE IN THE BACKGROUND
max_car_speed = '45' #maximum car speed to allow (mph)
min_car_speed = '2' # minimum car speed to allow (mph)
baseline_percentile = '50' #what percentile to use as a backgorund calculation
buffer_distance = '30' # distance of buffer (m) to use

###############################################################################
###### DON'T CHANGE ANYTHING BELOW THIS (UNLESS YOU CAN FIX IT) ###############
###############################################################################

# STARTING ALGORITHM (NAMING FOLDERS AND SUCH) 
## WHERE TO PUT LEAKS
raw_data_dir =  results_folder_loc + 'RawData/'
observed_peaks_dir = results_folder_loc + 'ObservedPeaks/'

filtered_observed_peaks_dir = results_folder_loc + 'FilteredObservedPeaks/'
final_results_dir = results_folder_loc + 'FinalShpFiles/'
shp_file_loc = final_results_dir + 'verifiedPKs.json'
processedFileLoc = results_folder_loc + 'ProcessedData/'
op_shp_file_loc = final_results_dir + "OP_Final.json"
all_op_csv_loc = final_results_dir + 'overallPeaks.csv'
final_info_loc = final_results_dir + 'summaryInfo.csv'
final_main_csv_loc = final_results_dir + 'mainThing.csv'

####### POTENTIALLY COULD CHANGE
s1 = car_id
s2 = "Peaks_" + str(s1)
s3 = "Filtered" + str()

##################################
### IMPORTING NECESSARY MODULES
##################################

import sys
sys.path.insert(0, function_file_Loc) # FINDING FUNCTIONS FOLDER TO IMPORT FROM
from amld_Functions import unique,unIfInt,\
                            intersect,weighted_loc,verPk,estimate_emissions,\
                            haversine,wt_time_Locs,sum_values,make_GEO,\
                            make_GPD,summarize_dat,get_quadrant,calc_theta,\
                            calc_bearing,process_raw_data_eng,str_list,\
                            count_times,identify_peaks,filter_peaks,\
                            pass_combine,summarize_data_2,add_odometer,process_raw_data,process_raw_data_aeris,identify_peaks_CSU

import rtree, pygeos,os, sys, datetime, time, math, numpy, csv, gzip,shutil,ast,swifter
from math import radians, sin, cos, sqrt, asin
import numpy as np
from numpy import log
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import matplotlib.pyplot as plt
from datetime import datetime


#### CREATING NECESSARY FOLDERS
addingFiles = False
foldList = [raw_data_dir,results_folder_loc,observed_peaks_dir,filtered_observed_peaks_dir,final_results_dir,processedFileLoc]
for x in foldList:
    if os.path.isdir(x) == False:
        try:
            os.mkdir(x)
        except OSError:
            print("Creation of the directory %s failed" % x)
        else:
            print("Successfully created the directory %s " % x)
            
### MOVING RAW FILES TO THE RAW DATA FILE FOLDER
for file in os.listdir(raw_data_loc):
    if file.endswith(".txt"):
        shutil.move(raw_data_loc+'/' + file,raw_data_dir)
########################################################################################

##### THIS PORTION OF THE CODE ALLOWS US TO ITERATIVELY ADD IN MORE DATA
#        PUT THE NEW TEXT FILES INTO THE OVERALL FOLDER AND IT WILL DO THE REST
raw_txts = pd.DataFrame(os.listdir(raw_data_dir)).loc[pd.DataFrame(os.listdir(raw_data_dir))[0].str.endswith('.txt')]

### DON'T ADD NEW FILE WITH PRE-EXISTING DATE [NEED TO WRITE CODE TO DEAL WITH THAT]
to_analyse, to_identify, to_filter = [[] for _ in range(3)]


if os.path.exists(final_info_loc):
    analysed = pd.read_csv(final_info_loc)
    for index,row in raw_txts.reset_index().iterrows():
        text = raw_txts[0].iloc[index]
        if analysed[analysed['FILENAME'].astype(str).str.contains(text)].shape[0] < 1:
            to_analyse.append(text)
            to_identify.append(s1 + '_20' + text[11:17] + '_dat.csv')
            to_filter.append(s2 + '_20' + text[11:17] + '.csv')
elif not os.path.exists(final_info_loc):
    for index,row in raw_txts.reset_index().iterrows():
        text = raw_txts[0].iloc[index]
        to_analyse.append(text)
        to_identify.append(s1 + '_20' + text[11:17] + '_dat.csv')
        to_filter.append(s2 + '_20' + text[11:17] + '.csv')

##### START OF THE ALGORITHM
start = time.time()

if __name__ == '__main__':
    date_list = []
    x1=""
    count = 0
    for file in to_analyse:
        #print (file)
        if file.endswith(".txt"):
            dateF = file[11:17]
            if dateF in set(date_list):
                bFirst = False
            
            if dateF not in set(date_list):
                bFirst = True
                date_list.append(dateF)
            x1 = file[:10]
            x_date = file[:10]
            if engineering:
                the_result = process_raw_data_eng(car_id, x_date, raw_data_dir, file, bFirst, 1, processedFileLoc,initial_time_ignore,shift,max_car_speed,min_car_speed)
            elif not engineering:
                if aeris:
                    the_result = process_raw_data_aeris(car_id, x_date, raw_data_dir, file, bFirst, 1, processedFileLoc,initial_time_ignore,shift,max_car_speed,min_car_speed)
                elif not aeris:
                    the_result = process_raw_data(car_id, x_date, raw_data_dir, file, bFirst, 1, processedFileLoc,initial_time_ignore,shift,max_car_speed,min_car_speed)
            count = count + 1

if __name__ == '__main__':
    for file in to_identify:
        if file.startswith(s1) and file.endswith("dat.csv"):
            x_date = file[len(car_id)+1:len(car_id) + 9]
            if CSU:
                the_result = identify_peaks_CSU(car_id, x_date, processedFileLoc, file, observed_peaks_dir, processedFileLoc, engineering,
                                          threshold, time_thresh, min_elevated, back_obs_num, baseline_percentile)
            elif not CSU:
                the_result = identify_peaks(car_id, x_date, processedFileLoc, file,observed_peaks_dir,processedFileLoc,engineering,threshold,time_thresh,min_elevated,back_obs_num,baseline_percentile)
if __name__ == '__main__':
    index = 0
    numproc = 0
    for file in os.listdir(observed_peaks_dir):
        if file.startswith(s2) and (file.endswith('.csv') and not file.endswith('info.csv')):
            file_loc = observed_peaks_dir + file
            csv_loc = observed_peaks_dir  + file[:-4]+ '_info.csv'
            non_empty = False
            if pd.read_csv(file_loc).size != 0:
                index += 1
                non_empty = True
                x_date = file[len(car_id) + 7: len(car_id) + 15]
                filter_peaks(car_id,x_date,observed_peaks_dir,file,filtered_observed_peaks_dir,buffer = buffer_distance,whichpass = index )

if not os.path.exists(final_main_csv_loc):
    toCombine = os.listdir(filtered_observed_peaks_dir)
    toCombineList = []
    index = 0
    for file in toCombine:
        if file.startswith('Filtered') and file.endswith('.csv') and not file.endswith('info.csv'):
            toCombineList.append(file)
    index = 1
    for file in toCombineList:
        if index == 1:
            loc = filtered_observed_peaks_dir + file
            mainInfo = pd.read_csv(filtered_observed_peaks_dir + file[:-4] + '_info.csv')
            mainThing = (pd.read_csv(loc))
        elif index != 1:
            loc2 = filtered_observed_peaks_dir + file
            secThing = pd.read_csv(loc2)
            secInfo = pd.read_csv(filtered_observed_peaks_dir + file[:-4] + '_info.csv')
            woo = pass_combine(mainThing,secThing,car_id, buffer_distance)
            mainInfo = pd.concat([secInfo,mainInfo])
            mainThing = woo.copy()
        index = index + 1
        print(file)
    mainThing = mainThing.copy().reset_index(drop = True)
    mainThing['numtimes']  = mainThing.apply(lambda x: count_times(x.recombine,car_id),axis=1)

elif os.path.exists(final_main_csv_loc):
    addingFiles = True
    toCombine = list(map(lambda x: 'FilteredPeaks_' + x[6:],to_filter))
    mainThing = pd.read_csv(final_main_csv_loc)
    mainInfo = pd.read_csv(final_info_loc)
    curVP = mainThing.loc[mainThing.numtimes > 1, :].min_read.drop_duplicates().shape[0]
    curOP = mainThing.min_read.drop_duplicates().shape[0]
    for file in toCombine:
        loc = filtered_observed_peaks_dir + file
        if os.path.exists(loc):
            secThing = pd.read_csv(loc)
            mainThing = pass_combine(mainThing,secThing,car_id, buffer_distance)
            mainInfo = pd.concat([mainInfo,pd.read_csv(filtered_observed_peaks_dir + file[:-4] + '_info.csv')])
            print(file)
    mainThing = mainThing.copy().reset_index(drop = True)
    mainThing['numtimes']  = mainThing.apply(lambda x: count_times(x.recombine,car_id),axis=1)

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