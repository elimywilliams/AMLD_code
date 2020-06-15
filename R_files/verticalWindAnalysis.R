### LIBRARIES
library(readr)

### CHOOSE DF
df <- read_csv("~/Documents/MobileSurveyRaw/CSULi_20191118_dat.csv")

#### LIBRARIES
library(tidyr)
library(dplyr)
library(zoo)
library(readr)
library(openair)
library(ggplot2)
library(grid)
#########

dayDat <- df
hz <- 10
min <- 5
width <- min*60*hz
width <- 1020

#### ADD LAG FOR THE TCH4 (baseline calculation)
## width is how many before/after
df <- df %>%
  mutate(CH4_baseline = rollmean(TCH4, k = width, fill = NA),
         ab_base = -CH4_baseline + TCH4,
         ab_base_perc = TCH4/CH4_baseline) ## how much above baseline!
df2 <- df %>% 
  mutate(secnan = SECONDS + (NANOSECONDS*1e-9)) %>% 
  mutate(prev_lat = as.numeric(lag(LAT,1,default=NA)),
         next_lat = as.numeric(lead(LAT,1,default=NA)),
         prev_lon = as.numeric(lag(LONG,1,default=NA)),
         next_lon = as.numeric(lead(LONG,1,default=NA)),
         prev_time = as.numeric(lag(secnan,1,default=NA)),
         next_time = as.numeric( lead(secnan,1,default=NA))) %>% 
  group_by(DATE,TIME,SECONDS,NANOSECONDS) %>% 
  dplyr::mutate(dist = gcd.hf3(prev_lon,prev_lat,next_lon,next_lat)) %>% 
  dplyr::mutate(bearing = calcBearing(lat1 = prev_lat,lat2 = next_lat,long1 = prev_lon,long2 = next_lon,radians = radians)) %>% 
  filter(dist < 1800) %>% 
  mutate(dist = ifelse(dist > 1e5,NA,dist)) %>% 
  mutate(timediff = next_time - prev_time) %>% 
  mutate(speed = dist/timediff) %>% 
  mutate(U_prime = U + speed) %>% 
  rename(U_old = U) %>% 
  mutate(U = U_prime)

correctedDf <- addWindDat2(df2,"U","V","W",F,T,T,T) %>% 
  mutate(adjtheta = (bearing + theta) %% 360) %>% 
  filter(!is.na(U),!is.na(V)) %>% 
  mutate(U_adj = r*cos(adjtheta*pi/180),
         V_adj = r*sin(adjtheta*pi/180)
  )

###################################

#### LOOK AT HISTOGRAM OF VERTICAL WIND
hist(abs(correctedDf$W* (180/pi)%% 180),main="Vertical wind (degrees) in drive")
correctedDf$w_cor <- abs(correctedDf$W* (180/pi)%% 180)
windRose(correctedDf,ws = 'r',wd='phi')

hist(correctedDf$phi)

correctedDf %>% 
  mutate(num = 1:n()) %>% 
  ggplot(aes(x=num,y=phi,col=r)) + geom_point()









