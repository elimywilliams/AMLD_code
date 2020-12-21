import pandas as pd
testWind = pd.read_csv('/Users/emilywilliams/Downloads/September 22 Data 3/AggregatedData.csv').loc[30:,]
testWind2 = testWind.loc[:,['timestamp','latitude','longitude','windSpeed','windAngle']]
testWind2['seconds'] = testWind2.apply(lambda x:math.floor(x.timestamp),axis=1 )
testWind3 = pd.DataFrame(testWind2.groupby('seconds',as_index=False).first())
testWind4 = testWind3.loc[(testWind3.latitude>20) & (testWind3.latitude<50),:]

def haversine(lat1, lon1, lat2, lon2, radius=6371):
    """ calculate the distance between two gps coordinates, using haversine function
    input:
        lat1,lon1: location 1
        lat2,lon2: location 2
        radius: earth's radius at the location used (km). Default is 6371km
    output:
        distance between the points (m)
    """
    from math import radians, sin, cos, sqrt, asin
    dLat = radians(lat2 - lat1)
    dLon = radians(lon2 - lon1)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    c = 2 * asin(sqrt(sin(dLat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dLon / 2) ** 2))
    return radius * c * 1000  # return in meters
def calc_bearing(lat1, lat2, long1, long2, radians):
    """ calculating the direction (bearing) of driving
    input:
        lat1,lon1: location 1
        lat2,lon2: location 2
        radians: T/F
    output:
        direction (degrees/radians) of motion
    """
    from math import atan2
    from numpy import pi
    from math import radians, sin, cos

    lat1r = lat1 * (pi / 180)
    lat2r = lat2 * (pi / 180)
    long1r = long1 * (pi / 180)
    long2r = long2 * (pi / 180)
    X = cos(lat2r) * sin(long2r - long1r)
    Y = cos(lat1r) * sin(lat2r) - (sin(lat1r) * cos(lat2r) * cos(long2r - long1r))

    theta = atan2(X, Y)
    theta = theta % (2 * pi)

    if not radians:
        return (theta * 180 / pi)
    elif radians:
        return (theta)

def calc_velocity(timediff, distance):
    if timediff == 0:
        return (0)
    elif timediff != 0:
        return (abs(distance / timediff))
def calc_car_velocity(df2,lat,long,timecol):
    df = df2.loc[(df2[lat]>30) & (df2[lat]<50),:]
    df['prev_LAT'] = df[lat].shift(periods=1)
    df['next_LAT'] = df[lat].shift(periods=-1)
    df['prev_LONG'] = df[long].shift(periods=1)
    df['next_LONG'] = df[long].shift(periods=-1)
    df['prev_TIME'] = df[timecol].shift(periods=1)
    df['next_TIME'] = df[timecol].shift(periods=-1)
    df['distance'] = df.apply(
        lambda row: haversine(row['prev_LAT'], row['prev_LONG'], row['next_LAT'], row['next_LONG']), axis=1)
    df['bearing'] = df.apply(
        lambda row: calc_bearing(row['prev_LAT'], row['next_LAT'], row['prev_LONG'], row['next_LONG'], True),
        axis=1)
    df['timediff'] = df.apply(lambda row: row['next_TIME'] - row['prev_TIME'], axis=1)
    df['VELOCITY'] = df.apply(lambda row:calc_velocity(row['timediff'],row['distance']),axis=1)
    df['nullVel'] = df.loc[:,'VELOCITY'].isnull()
    df['VELOCITY'] = df.apply(lambda row: 0 if row.nullVel else row.VELOCITY,axis=1)
    return(df.drop(columns = ['nullVel']))


df3 = calc_car_velocity(testWind4,'latitude','longitude','seconds')

df3['delta'] = df3.apply(lambda x: 10*(math.pi/180),axis=1)
df3['theorWindNS'] = df3.apply(lambda x: math.cos(x.bearing)*(x.VELOCITY + x.windSpeed *math.cos(x.windAngle+ x.delta)) -math.sin(x.windAngle+ x.delta)*math.sin(x.bearing)*x.windSpeed,axis=1 )
df3['theorWindEW'] = df3.apply(lambda x: math.sin(x.bearing)*(x.VELOCITY + x.windSpeed *math.cos(x.windAngle+ x.delta)) +math.sin(x.windAngle+ x.delta)*math.cos(x.bearing)*x.windSpeed,axis=1 )
df3['theorWS'] = df3.apply(lambda x: math.sqrt(x.theorWindNS**2 +x.theorWindEW**2 ),axis=1)

bigWS =df3.loc[df3.theorWS>20,:]
smallWS =df3.loc[df3.theorWS<20,:]