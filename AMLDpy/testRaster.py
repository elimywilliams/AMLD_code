#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 14:27:01 2020

@author: emilywilliams
"""

import numpy as np
import pandas as pd
loc = '/Users/emilywilliams/Documents/DrivingData/ColDatShort/ProcessedData/SCcar_20200406_dat.csv'
df = pd.read_csv(loc)


bottomLeft = (df.LAT.min(), df.LONG.min())
bottomRight = (df.LAT.min(), df.LONG.max())
topLeft = (df.LAT.min(), df.LONG.max())
topRight = (df.LAT.max(), df.LONG.max())


cols = np.linspace(bottomLeft[1], bottomRight[1], num=18)
rows = np.linspace(bottomLeft[0], topLeft[0], num=15)
df['col'] = np.searchsorted(cols, df['LONG'])
df['row'] = np.searchsorted(rows, df['LAT'])