### LIBRARIES
library(readr)

### CHOOSE DF
df1 <- read_csv("~/Documents/MobileSurveyRaw/CSULi_20191118_dat.csv")
df2 <- read_csv("~/Documents/MobileSurveyRaw/CSULi_20191022_dat.csv")
df3 <- read_csv("~/Documents/MobileSurveyRaw/CSULi_20191115_dat.csv")
df4 <- read_csv("~/Documents/MobileSurveyRaw/CSULi_20191113_dat.csv")
df5 <- read_csv("~/Documents/MobileSurveyRaw/CSULi_20191124_dat.csv")

df <- plyr::join(df1,df2,type='full')
df <- plyr::join(df,df3,type='full')
df <- plyr::join(df,df4,type='full')
df <- plyr::join(df,df5,type='full')

df_all <- df

df <- read_delim("~/Downloads/2020-06-09_12h44m07s.txt", 
                 "\t", escape_double = FALSE, col_types = cols(DATAH = col_skip()), 
                 trim_ws = TRUE) %>% 
  dplyr::rename(LAT = GPS_LAT,LONG = GPS_LONG) %>% 
  mutate(TCH4 = CH4) %>% 
  rename(U = SONIC_U,V=SONIC_V,W=SONIC_W)

##combining all of the datasets
list <- list.files('/Users/emilywilliams/Documents/MobileSurveyRaw')
nlist <- length(list)
index <- 0
for (file in list) {
  if (endsWith(file,'dat.csv')) {
    if(index == 0){
      df <- read_csv(paste('/Users/emilywilliams/Documents/MobileSurveyRaw/',file,sep=''))
    }
    if(index !=0){
      df1 <- read_csv(paste('/Users/emilywilliams/Documents/MobileSurveyRaw/',file,sep=''))
      df <- plyr::join(df,df1,type='full')
    }
    index <-index + 1
  }
}


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
correctedDf_stat <- correctedDf
correctedDf_all <- correctedDf
#### LOOK AT HISTOGRAM OF VERTICAL WIND
hist(abs(correctedDf$W* (180/pi)%% 180),main="Vertical wind (degrees) in drive")
correctedDf$r_vert <- correctedDf$r*sin((correctedDf$phi)*(pi/180))
windRose(correctedDf,ws = 'r_vert',wd='phi')

hist(correctedDf$phi)

correctedDf_stat %>% 
  mutate(num = 1:n()) %>% 
  ggplot(aes(x=num,y=phi,col=log(r))) + geom_point()

hist(correctedDf$phi,prob=T,ylim=c(0,.06),main="Distribution of Vertical Wind \n (0 degrees= horizontal)",col='gray',xlab='Vertical Direction (degrees)')
lines(density(correctedDf2$phi),col='forestgreen',lwd = 2)

quantile(correctedDf$phi,c(0.80))

correctedDf %>% 
  filter(phi > 20) %>% 
  mutate(num = 1:n()) %>% 
  ggplot(aes(x=num,y=phi,col=log(VELOCITY))) + geom_point()


plot(correctedDf$phi,correctedDf$GPS_VELOCITY,main="Car velocity vs. Vertical Wind Component")

hist(correctedDf_stat$r * cos(correctedDf_stat$phi * pi/180),main='Hist of vertical wind speed (stationary)',prob=T,ylim=c(0,.4),col='gray')
hist(correctedDf_all$r * cos(correctedDf_all$phi * pi/180),add=T,prob=T,col='lightblue')

par(mfrow = c(1,2))
boxplot(correctedDf_all$r * cos(correctedDf_all$phi * pi/180),main='All Data',ylab='Vert WS (m/s)',ylim=c(0,25))
boxplot(correctedDf_stat$r * cos(correctedDf_stat$phi * pi/180),main = "Stationary Data",ylab='Vert WS (m/s)',ylim=c(0,25))


correctedDf_stat %>% 
  mutate(num = 1:n()) %>% 
  ggplot(aes(x=num,y = r*cos(phi*pi/180)/r)) + geom_point()

correctedDf_all %>% 
  mutate(num = 1:n()) %>% 
  ggplot(aes(x=num,y = r*cos(phi*pi/180)/r)) + geom_point()

par(mfrow=c(1,2))
hist((cos(correctedDf_all$phi*pi/180)),prob=T,main="All Drives",xlab='Proportion of wind from vertical:total')
hist((cos(correctedDf_stat$phi*pi/180)),prob=T,main="Stationary",xlab='Proportion of wind from vertical:total')

par(mfrow=c(1,2))
boxplot((cos(correctedDf_all$phi*pi/180)),main="All Drives",xlab='Proportion of wind from vertical:total')
boxplot((cos(correctedDf_stat$phi*pi/180)),main="Stationary",xlab='Proportion of wind from vertical:total')

par(mfrow=c(1,2))
boxplot(exp(cos(correctedDf_all$phi*pi/180)),main="All Drives",ylab='Proportion of wind from vertical:total')
boxplot(exp(cos(correctedDf_stat$phi*pi/180)),main="Stationary",ylab='Proportion of wind from vertical:total')


plot(correctedDf$phi)

par(mfrow = c(1,2))
plot(y=correctedDf$r*cos(correctedDf$phi*(pi/180)),
     x=correctedDf$GPS_VELOCITY,main = "Car Velocity vs. Vertical Wind Speed",
     xlab= "Car Velocity (m/s)",ylab="Vertical Wind Component (m/s)")
plot(y=correctedDf$phi,
     x=correctedDf$GPS_VELOCITY,main = "Car Velocity vs. Vertical Wind Direction",
     xlab= "Car Velocity (m/s)",ylab="Vertical Wind Component (degrees)")

fastwind <- correctedDf[correctedDf$r*cos(correctedDf$phi*(pi/180))>median(correctedDf$VELOCITY),]


par(mfrow=c(1,2))
plot(y = fastwind$r*cos(fastwind$phi*(pi/180)),x= fastwind$VELOCITY,xlab = "Velocity (m/s)", ylab="Vertical Wind Speed (m/s)",
     main = "Wind speeds (Vertical and Horizontal) \n vs. Car Velocity")
plot(y=fastwind$r*sin(fastwind$phi*(pi/180)),x= fastwind$VELOCITY, ylab="Horizontal Wind Speed (m/s)",xlab='Velocity (m/s)')


nohorz <- correctedDf[correctedDf$r*cos(correctedDf$phi*(pi/180)) <.75,]
par(mfrow=c(1,2))
plot(y = nohorz$r*cos(nohorz$phi*(pi/180)),x= nohorz$VELOCITY,xlab = "Velocity (m/s)", ylab="Vertical Wind Speed (m/s)",
     main = "Wind speeds (Vertical and Horizontal) \n vs. Car Velocity")
plot(y=nohorz$r*sin(nohorz$phi*(pi/180)),x= nohorz$VELOCITY, ylab="Horizontal Wind Speed (m/s)",xlab='Velocity (m/s)')


plot(fastwind$ab_base_perc)

plot(correctedDf$ab_base_perc)

plot(correctedDf$ab_base_perc,correctedDf$r)

hist(correctedDf$VELOCITY)

correctedDf %>% 
  ggplot(aes(x=r,y=ab_base_perc,col=log(VELOCITY))) + geom_point() +
  ggtitle("CH4 (percentage above baseline) vs. Wind Speed (r)")+
  ylab("CH4 (percentage above baseline)") + xlab("Wind Speed (m/s)") + 
  ggthemes::theme_clean() + geom_hline(yintercept = 1.1) 

###
plot(correctedDf_stat$W)
points(correctedDf_stat$r * cos(correctedDf_stat$phi * (pi/180)),col='blue')



