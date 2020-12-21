## use true wind to make

# Actual wind speed (m/s)
wind_true = 10

# Actual wind angle (rad, counterclockwise from x axis)
phi_true = 0*pi/2

# Actual Car velocity (m/s)
car_velocity = 10
theta_compass = math.pi/2

# Actual car heading (rad, counterclockwise from x axis)
car_heading = math.pi/2 - theta_compass


wtrue_x = wind_true* math.cos(phi_true)
wtrue_y = wind_true* math.sin(phi_true)

theta_star = - car_heading

wind_car_x = math.cos(theta_star)*wtrue_x - math.sin(theta_star)*wtrue_y
wind_car_y = math.sin(theta_star)*wtrue_x + math.cos(theta_star)*wtrue_y

phi = math.atan((wind_car_y)/(wind_car_x + car_velocity))
w = (wind_car_x + car_velocity)/math.cos(phi)

phi_anem = 3*math.pi/2 - phi
wind_anem = w
theta_compass = math.pi/2 - car_heading

## pick own input
theta_compass = 15 * (math.pi/180)
phi_anem = 30 * (math.pi/180)
car_velocity = 10
wind_anem = 10


theta = (math.pi/2 - theta_compass)%(2*math.pi)
phi = (3*math.pi/2 - phi_anem)%(2*math.pi)

wind_true_x = wind_anem*math.cos(theta)*math.cos(phi) - car_velocity*math.cos(theta) - wind_anem*math.sin(theta)*math.sin(phi)
wind_true_y = wind_anem*math.sin(theta)*math.sin(phi) - car_velocity*math.sin(theta) + wind_anem*math.cos(theta)*math.sin(phi)

if wind_true_x == 0:
    wind_theta_true = 0
elif wind_true_x !=0:
    wind_theta_true = math.atan(wind_true_y/wind_true_x)
wind_speed_true = math.sqrt(wind_true_x**2 + wind_true_y**2)


import pandas as pd
from math import radians, sin, cos, sqrt, asin
wind_true = 10
phi_true = 0

## inputting wind
runDat = pd.read_csv('/Users/emilywilliams/Downloads/September22Data/simulateWindRun.csv')
#d = {'VELOCITY': [15, 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15], 'bearing': [0, 0,0,0,0,0,math.pi/4,math.pi/4,math.pi/4,math.pi/4,math.pi/4,math.pi/4,
#                                                                               95 * (math.pi/180),95 * (math.pi/180),95 * (math.pi/180),95 * (math.pi/180),95 * (math.pi/180)],
#     'timestamp': [1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]}
#runDat = pd.DataFrame(data=d)
runDat['secnan'] = runDat.apply(lambda x: (x.timestamp),axis =1)
runDat2 = calcVehicleSpeed2(runDat)
#runDat2 = runDat.copy()
runDat2['car_heading'] = runDat2.apply(lambda x:math.pi/2 - x.bearing,axis=1)
runDat2['car_velocity'] = runDat2.apply(lambda x: x.VELOCITY,axis=1)
runDat2['wtrue_x']= runDat2.apply(lambda x: wind_true* math.cos(phi_true),axis=1)
runDat2['wtrue_y'] = runDat2.apply(lambda x: wind_true* math.sin(phi_true),axis=1)
runDat2['theta_star'] = runDat2.apply(lambda x: -x.car_heading,axis=1)
runDat2['wind_car_x'] = runDat2.apply(lambda x: math.cos(x.theta_star)*x.wtrue_x - math.sin(x.theta_star)*x.wtrue_y,axis =1)
runDat2['wind_car_y'] = runDat2.apply(lambda x: math.sin(x.theta_star)*x.wtrue_x + math.cos(x.theta_star)*x.wtrue_y,axis =1)

runDat2['phi'] = runDat2.apply(lambda x: math.atan((x.wind_car_y)/(x.wind_car_x + x.car_velocity)),axis=1)
runDat2['w'] = runDat2.apply(lambda x: (x.wind_car_x + x.car_velocity)/math.cos(x.phi),axis=1)
runDat2['phi_anem'] = runDat2.apply(lambda x: 3*math.pi/2 - x.phi, axis =1)
runDat2['wind_anem'] = runDat2.apply(lambda x: x.w,axis = 1)
runDat2['theta_compass'] = runDat2.apply(lambda x: math.pi/2 - x.car_heading,axis=1)

## back calculating wind
runDat2['theta'] = runDat2.apply(lambda x: (math.pi/2 - x.theta_compass)%(2*math.pi),axis=1)
runDat2['phi_calc'] = runDat2.apply(lambda x: (3*math.pi/2 - x.phi_anem)%(2*math.pi),axis=1)
runDat2['wind_true_x'] = runDat2.apply(lambda x: x.wind_anem*math.cos(x.theta)*math.cos(x.phi_calc) - x.car_velocity*math.cos(x.theta) - x.wind_anem*math.sin(x.theta)*math.sin(x.phi_calc),axis=1)
runDat2['wind_true_y'] = runDat2.apply(lambda x: x.wind_anem*math.sin(x.theta)*math.cos(x.phi_calc) - x.car_velocity*math.sin(x.theta) + x.wind_anem*math.cos(x.theta)*math.sin(x.phi_calc),axis=1)
runDat2['wind_theta_true'] = runDat2.apply(lambda x: 0 if x.wind_true_x == 0 else math.atan(x.wind_true_y/x.wind_true_x),axis=1 )
runDat2['wind_speed_true'] = runDat2.apply(lambda x: math.sqrt(x.wind_true_x**2 + x.wind_true_y**2),axis=1)
runDat2.to_csv('runSimulated.csv', index=False)



w = (wind_car_x + car_velocity)/math.cos(phi)






runDat2['U_car'] = runDat2.apply(lambda x: u_fr_rtheta(x.VELOCITY,x.bearing),axis = 1)
runDat2['V_car'] = runDat2.apply(lambda x: v_fr_rtheta(x.VELOCITY,x.bearing),axis = 1)

runDat2['U_twind'] = runDat2.apply(lambda x: u_fr_rtheta(x.TWS,x.TWD),axis = 1)
runDat2['V_twind'] = runDat2.apply(lambda x: v_fr_rtheta(x.TWS,x.TWD),axis = 1)

runDat2.to_csv('runSimulated.csv', index=False)

def u_fr_rtheta(r,theta):
    return(-r * math.cos(theta))
def v_fr_rtheta(r,theta):
    return(r * math.sin(math.pi - theta))

angles = [0,math.pi/4,2*math.pi/4,3*math.pi/4,4*math.pi/4,5*math.pi/4,6*math.pi/4,7*math.pi/4]

for index,angle in enumerate(angles):
    print(angle)
    print('U:' + str(u_fr_rtheta(1*np.sqrt(2),angle)))
    print("V: " + str(v_fr_rtheta(1*np.sqrt(2),angle)))

carDat = pd.read_csv("/Users/emilywilliams/Downloads/September22Data/simulatedWind.csv")

carDat['seconds'] = carDat.apply(lambda x: np.floor(x.timestamp),axis=1)

# convert from polar to UV coordinates of raw wind
carDat['U'] = carDat.apply(lambda x: u_fr_rtheta(x.WS,x.WA_rad),axis = 1)
carDat['V'] = carDat.apply(lambda x: v_fr_rtheta(x.WS,x.WA_rad),axis = 1)

### summarize the file to have just one second readings
summaryCar = carDat.groupby('seconds',as_index=False).mean()
carDat = summaryCar.copy()
carDat = carDat.rename(columns = {'latitude':'LAT','longitude':'LONG'})
carDat['secnan'] = carDat.apply(lambda x: (x.timestamp),axis =1)


def calc_velocity(timediff, distance):
    if timediff == 0:
        return (0)
    elif timediff != 0:
        return (distance / timediff)
def calcVehicleSpeed2(df):
    df['prev_LAT'] = df.LAT.shift(periods=1)
    df['next_LAT'] = df.LAT.shift(periods=-1)
    df['prev_LONG'] = df.LONG.shift(periods=1)
    df['next_LONG'] = df.LONG.shift(periods=-1)
    df['prev_TIME'] = df.secnan.shift(periods=1)
    df['next_TIME'] = df.secnan.shift(periods=-1)
    df['distance'] = df.apply(
        lambda row: haversine(row['prev_LAT'], row['prev_LONG'], row['next_LAT'], row['next_LONG']), axis=1)
    df['bearing'] = df.apply(
        lambda row: calc_bearing(row['prev_LAT'], row['next_LAT'], row['prev_LONG'], row['next_LONG'], radians),
        axis=1)
    df['timediff'] = df.apply(lambda row: row['next_TIME'] - row['prev_TIME'], axis=1)
    df['VELOCITY'] = df.apply(lambda row:calc_velocity(row['timediff'],row['distance']),axis=1)
    df['nullVel'] = df.loc[:,'VELOCITY'].isnull()
    df['VELOCITY'] = df.apply(lambda row: 0 if row.nullVel else row.VELOCITY,axis=1)
    return(df.drop(columns = ['nullVel']))

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
def correctU(df):
    df['U_cor'] = df.apply(lambda x: x.U + x.VELOCITY,axis = 1)
    return(df)
def uv_to_angle(u, v):
    from math import atan2
    from math import pi
    theta = (atan2(v, u) + pi) % (2 * pi)
    if u == 0 and v > 0:
        return (pi / 2)
    if u == 0 and v < 0:
        return (3 * pi / 2)
    if v == 0 and u >= 0:
        return (pi)
    if v == 0 and u < 0:
        return (0)
    quad = get_quadrant(u, v)
    if round(theta, 5) == round(2 * pi, 5):
        return (0)
    if quad == 1:
        return (theta + pi / 2)
    elif quad == 2:
        return (2 * pi - theta)
    elif quad == 3:
        return (pi / 2 + theta)
    elif quad == 4:
        return (pi - theta)
def polar_from_uv(df):
    from math import sqrt
    df['theta_wind'] = df.apply(lambda x: uv_to_angle(x.U_cor,x.V),axis=1)
    df['r_wind'] = df.apply(lambda x: sqrt(x.U**2 + x.V**2),axis=1)
    return(df)
def correct_bearing(df):
    df['bearing_true'] = df.apply(lambda x:x.bearing + (x.theta_wind + math.pi)%(2*(math.pi)),axis=1)
    return(df)
def uv_from_polar(df):
    df['U_final'] = df.apply(lambda x: x.r_wind * math.cos(x.bearing_true),axis = 1)
    df['V_final'] = df.apply(lambda x: -x.r_wind * math.sin(x.bearing_true),axis = 1)
    return(df)

def u_fr_rtheta(r,theta):
    return(r * math.cos(theta))
def v_fr_rtheta(r,theta):
    return(-r * math.sin(theta))
def get_quadrant(x, y):
    """ given an x,y coordinate, return which quadrant it is in
    input:
        x,y values
    output:
        quadrant
    exceptions:

    """
    try:
        x = int(x)
        y = int(y)
    except ValueError:
        return (0)

    if y >= 0 and x > 0:
        return (1)
    elif y >= 0 and x < 0:
        return (2)
    elif y < 0 and x < 0:
        return (3)
    else:
        return (4)


df = calcVehicleSpeed2(carDat)
df = carDat.copy()
df['bearing'] = df.apply(lambda x: x.bearing%(2*math.pi),axis=1)
df['U_cor'] = df.apply(lambda x: -x.U + x.VELOCITY,axis=1)
df['r_final'] = df.apply(lambda x: np.sqrt(x.U_cor**2 + x.V**2),axis=1)
df['theta_semifinal'] = df.apply(lambda x: uv_to_angle(x.U_cor,x.V),axis=1)
df['theta_final'] = df.apply(lambda x: (x.bearing + x.theta_semifinal)%(2*math.pi),axis=1)
df.to_csv('simulatedCorrected', index=False)

df = polar_from_uv(correctU(calcVehicleSpeed2(carDat)))
df2 = correct_bearing(df)
df = correct_bearing(polar_from_uv(correctU(calcVehicleSpeed(carDat))))

