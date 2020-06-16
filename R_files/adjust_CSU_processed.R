cdat <- read_csv('~/Documents/DrivingData/CoDrive/ProcessedData/SCcar_20200612_dat.csv')

odat <- read_csv("~/Documents/GitHub/CSU_SC/AlexDriving2/CSULi_20200612_dat.csv")

odat <- odat %>% 
  mutate(PRESS_MBAR = 0,
         INLET = 0,
         TEMPC = 0,
         CH4 = TCH4,
         H20 = 0,
         C2H6 = 0,
         R = 0,
         C2C1 = 0,
         BATTV = 0,
         POWMV = 0,
         CURRMA = 0,
         SOCPER = 0) %>% 
  dplyr::select(names(cdat))


list <- list.files("~/Documents/GitHub/CSU_SC/AlexDriving2/")
nlist <- length(list)
index <- 0
for (file in list) {
  if (endsWith(file,'dat.csv')) {
    temp <- read_csv(paste("~/Documents/GitHub/CSU_SC/AlexDriving2/",file,sep='')) %>% 
      mutate(PRESS_MBAR = 0,
             INLET = 0,
             TEMPC = 0,
             CH4 = TCH4,
             H20 = 0,
             C2H6 = 0,
             R = 0,
             C2C1 = 0,
             BATTV = 0,
             POWMV = 0,
             CURRMA = 0,
             SOCPER = 0) %>% 
      dplyr::select(names(cdat))
    
      write.csv(temp,paste('/Users/emilywilliams/Documents/DrivingData/aCoDrive/',file,sep=''))
  }
}
