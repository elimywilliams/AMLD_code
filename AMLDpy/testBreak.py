import pandas as pd
import geopandas as gpd
testfile = gpd.read_file('testgpd.geojson')
firstarea = testfile.loc[testfile.OP_NUM == 'TRUSS_1592439600.9',:]
secondarea = testfile.loc[testfile.OP_NUM == 'TRUSS_1592439707.4',:]
gpd.overlay(firstarea,secondarea,how='union')
