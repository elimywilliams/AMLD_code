foundLeaks = pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/truss_short/FinalShpFiles/mainThing.csv')

foundLeaks['avgDeg'] = foundLeaks.apply(lambda x: (x['OB_THETA_AVG']* 180/pi)%360,axis = 1)
lkGeo = [Point(lon, lat) for lon, lat in zip(foundLeaks['OB_LON'], foundLeaks['OB_LAT'])]
crs = 'EPSG:4326'
flGPD = gpd.GeoDataFrame(foundLeaks,crs = crs, geometry = lkGeo)
flGPD.to_file('finalLeaks.shp')
# geometry is the point of the lat/lon
# gdf_buff = gpd.GeoDataFrame(datFram, crs=crs, geometry=geometry_temp)

## BUFFER AROUND EACH 'OP_NUM' WITH BUFFER DISTANCE
gdf_buff = gpd.GeoDataFrame(fileWt, crs=crs, geometry=geometry_temp)
# gdf_buff = makeGPD(datFram,'LON','LAT')

import matplotlib.pyplot as plt
import numpy as np

import iris
import iris.coord_categorisation
import iris.quickplot as qplt

import cartopy
import cartopy.feature as cfeat
import cartopy.crs as ccrs

uwind = foundLeaks.OB_U_AVG
vwind = foundLeaks.OB_V_AVG

import matplotlib.pyplot as plt
import numpy as np
from numpy import ma

X=np.arange(0, 2 * np.pi, .2)
Y=np.ones(X.shape)
U= np.cos(X)
V= np.sin(X)
plt.figure()
plt.plot(X,X,'--')
#Q = plt.quiver(X, Y, Y, Y, units='width')
Q = plt.quiver(X, Y-2, U, V, units='width')
plt.show()

X = foundLeaks.OB_LAT
Y = foundLeaks.OB_LON
U= uwind
V= vwind
plt.figure()
plt.plot(X,Y, '.')
Q = plt.quiver(X,Y,U,V,units = 'width')
plt.show()
