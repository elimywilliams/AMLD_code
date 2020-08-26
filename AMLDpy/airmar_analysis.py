loc = '/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/SouthernCross_JvonFischer_Test.csv'
out_location = 'testWind2.csv'
def process_wind(file_location,out_location):
    import pandas as pd
    from datetime import datetime
    import swifter
    windFile = pd.read_csv(file_location,skiprows=[0,2,3],names = ['TIMESTAMP','RECORD','Ux_1','Uy_1','Uz_1',"Ts_1",
                        'Ux_2', 'Uy_2', 'Uz_2', "Ts_2",
                        'Ux_3', 'Uy_3', 'Uz_3', "Ts_3",
                        'Ux_4', 'Uy_4', 'Uz_4', "Ts_4",
                        'Ux_5', 'Uy_5', 'Uz_5', "Ts_5"])
    windFile["hasDate"]= windFile["TIMESTAMP"].swifter.apply(lambda x: '/' in x)
    firstIndex = windFile.loc[(windFile.hasDate.cumsum()==1),:].index[0]
    windFile2 = windFile.loc[windFile.index >= firstIndex,:].reset_index(drop=True)
    windFile2.loc[:,['date']] = windFile2.swifter.apply(lambda x: x.TIMESTAMP if x.hasDate else False,axis=1)
    windFile2.loc[:,['date2']] = windFile2.swifter.apply(lambda x: x.TIMESTAMP.split(' ')[0] if x.hasDate else False,axis=1)
    windFile2.loc[:,['hours']] = windFile2.swifter.apply(lambda x: x.TIMESTAMP.split(' ')[1].split(':')[0] if x.hasDate else False,axis=1)
    windFile2.loc[:,['mins']] = windFile2.swifter.apply(lambda x: x.TIMESTAMP.split(':')[0] if not x.hasDate else False,axis=1)
    windFile2.loc[:,['sec']] = windFile2.swifter.apply(lambda x: x.TIMESTAMP.split(':')[1].split('.')[0] if not x.hasDate else False,axis=1)
    windFile2.loc[:,['ptsec']] = windFile2.swifter.apply(lambda x: x.TIMESTAMP.split(':')[1].split('.')[1] if not x.hasDate else 0,axis=1)
    for x in range(windFile2.shape[0]):
        if windFile2.date2[x] != False:
            curdate = windFile2.date2[x]
            curtime = str(windFile2.hours[x]) + ":" + str(windFile2.mins[x+1]) + ":" + str(windFile2.sec[x+1]) + '.0'
        elif windFile2.date2[x] == False:
            curdate = windFile2.realdate[x-1]
            curtime = str(windFile2.realtime[x-1].split(':')[0]) + ":" + str(windFile2.mins[x]) + ":" + str(windFile2.sec[x]) + '.' + str(windFile2.ptsec[x])
        windFile2.loc[x,['realdate']] = curdate
        windFile2.loc[x,['realtime']] = curtime
    windFile2.loc[:,['DATETIME']] = windFile2.swifter.apply(lambda x: datetime(int(str(x.realdate).split('/')[2]),int(x.realdate.split('/')[0]),int(x.realdate.split('/')[1]),
             int(x.realtime.split(':')[0]),int(x.realtime.split(':')[1]),int(x.realtime.split(':')[2].split('.')[0]),
             int(x.realtime.split(':')[2].split('.')[1])*100000),axis = 1)
    windFile2 = windFile2.rename(columns= {'realdate':'DATE','realtime':'TIME'})
    windFile2.loc[:,['DATE','TIME','RECORD','Ux_1','Uy_1','Uz_1',"Ts_1",
                        'Ux_2', 'Uy_2', 'Uz_2', "Ts_2",
                        'Ux_3', 'Uy_3', 'Uz_3', "Ts_3",
                        'Ux_4', 'Uy_4', 'Uz_4', "Ts_4",
                        'Ux_5', 'Uy_5', 'Uz_5', "Ts_5"]].to_csv(out_location)
    #windFile2.drop(columns = ['ptsec','hours','mins','sec','date2','date','hasDate']).to_csv(out_location)
    return


st = time.time()
process_wind(loc,out_location)
print(time.time() - st)
