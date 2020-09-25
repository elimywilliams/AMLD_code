library(readr)
allElevated <- read_csv("~/Documents/GitHub/AMLD_Driving_Data/truss_6pc_102_med_shift_2_new/FinalShpFiles/mainThing.csv")
summary(allElevated$OB_C2H6_AB)
par(mfrow = c(1,2))
hist(allElevated$OB_R,prob=T,
     main="Histogram of R values in Elevated Readings, \n SC Algorithm",
     xlab='Observed R Readings',breaks = 10,xlim=c(-1,1))
abline(v = .7,col = 'blue',lwd=3)
hist(aeris_678_lks$R,prob=T,main="Histogram of R values in Elevated Readings, \n Aeris Algorithm",
     xlab='Observed R Readings',breaks = 10,xlim=c(-1,1))
abline(v = .7,col = 'blue',lwd=3)

par(mfrow = c(1,2))
hist(allElevated$OB_C2H6,prob=T,
     main="Histogram of C2H6 values in Elevated Readings, \n SC Algorithm",
     xlab='Observed C2H6 Readings',breaks = 10)
hist(aeris_678_lks$`C2H6 (ppb)`,prob=T,main="Histogram of C2H6 values in Elevated Readings, \n Aeris Algorithm",
     xlab='Observed C2H6 Readings',breaks = 10)

par(mfrow = c(1,2))
hist(allElevated[allElevated$OB_R >=.7,]$OB_C2C1,prob=T,
     main="Histogram of C2C1 values in Elevated Readings, \n SC Algorithm",
     xlab='Observed C2C1 Readings',breaks = 10)
hist(aeris_678_lks$`C2/C1`,prob=T,main="Histogram of C2C1 values in Elevated Readings, \n Aeris Algorithm",
     xlab='Observed C2C1 Readings',breaks = 10)

smallElevated <- (allElevated %>% 
  filter(OB_R >= .7))

plot(aeris_678_lks$`C2/C1`,aeris_678_lks$`CH4 (ppm)`/aeris_678_lks$`C2H6 (ppb)`)

hist(allElevated[allElevated$OB_R >=.7,]$OB_C2H6)
hist(allElevated$OB_C2H6)

hist(aeris_678_lks$`C2/C1`)
hist(allElevated$OB_C2C1)

hist(aeris_678_lks$`C2/C1`)
summary(allElevated[allElevated$OB_R >=.7,]$OB_C2C1)


allLeaksWin <- read_csv("~/Documents/GitHub/Trussville/allLeaksWin.csv", 
                        col_types = cols(`Time Stamp` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))

datesConsidered <- c(
  '2020-06-16','2020-06-17','2020-06-18','2020-06-19',
  '2020-06-20','2020-06-23','2020-06-24','2020-06-25',
  '2020-06-26','2020-06-27','2020-06-29','2020-06-30',
  '2020-07-01','2020-07-02','2020-07-06','2020-07-07',
  '2020-07-08','2020-07-09','2020-07-10','2020-07-11',
  '2020-07-13','2020-07-14','2020-07-15','2020-07-16','2020-07-18',
  '2020-07-19','2020-07-20','2020-07-21','2020-07-22',
  '2020-07-23','2020-07-24','2020-07-25','2020-07-26',
  '2020-07-27','2020-07-28',
  '2020-08-11','2020-08-12','2020-08-18','2020-08-19',
  '2020-08-20','2020-08-21','2020-08-24','2020-08-25',
  '2020-08-26','2020-08-27','2020-08-28','2020-08-29',
  '2020-08-30'
)

aeris_678_lks <- allLeaksWin %>% 
  mutate(DateFound = as.Date(`Time Stamp`)) %>% 
  filter(DateFound %in% as.Date(datesConsidered))

write_csv(aeris_678_lks,'aerisLks678.csv')



  filter(DateFound <= as.Date('2020-08-31')) %>% 
  filter(DateFound != as.Date('2020-06-15'))


library(readr)
SC_20200406_dat <- read_csv("~/Documents/GitHub/AMLD_Driving_Data/test_data_ethane/ProcessedData/SC_20200406_dat.csv") %>% 
  mutate(timer = 1:n())
SC_20200406_dat %>% 
  mutate(timer = 1:n()) %>% 
  ggplot(aes(x=timer,y=C2H6)) + geom_point() + geom_line()

SC_20200406_dat %>% 
  mutate(timer = 1:n()) %>% 
  ggplot(aes(x=timer,y=CH4)) + geom_point() + geom_line()


SC_20200406_dat %>% 
  mutate(timer = 1:n()) %>% 
  ggplot(aes(x=timer,y=C2C1)) + geom_point() + geom_line()

test2 <- SC_20200406_dat %>% 
  mutate(timer = 1:n()) %>% 
  tidyr::gather('Aeris_Measurement','Level',c(C2H6,CH4))  

test2 %>% 
  ggplot(aes(x=timer,y=Level,col = Aeris_Measurement))+ geom_point()

plot(SC_20200406_dat$timer,SC_20200406_dat$CH4)
plot(SC_20200406_dat$timer,SC_20200406_dat$CH4)

par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(SC_20200406_dat$timer, SC_20200406_dat$CH4, pch = 16, col = 'black',ylab = "CH4",
     xlab='Observation Number')              # Create first plot
par(new = TRUE)                             # Add new plot
plot(SC_20200406_dat$timer, SC_20200406_dat$C2H6, pch = 17, col = 'blue',              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "",main='Comparison, CH4, C2H6')
axis(side = 4, at = pretty(range( SC_20200406_dat$C2H6)))      # Add second axis
mtext("C2H6", side = 4, line = 3)

cor(SC_20200406_dat$CH4,SC_20200406_dat$C2H6)

plot(SC_20200406_dat$CH4,SC_20200406_dat$C2H6)


movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}

maC2H6 <- movingAverage(SC_20200406_dat$C2H6,n=102)
plot(maC2H6)

maCH4 <- movingAverage(SC_20200406_dat$CH4,n=102)
plot(maCH4)

par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(maCH4, pch = 16, col = 'black',ylab = "Moving Average CH4",
     xlab='Observation Number')              # Create first plot
par(new = TRUE)                             # Add new plot
plot(maC2H6,  pch = 17, col = 'blue',              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "",main='Comparison of Moving Average, CH4, C2H6')
axis(side = 4, at = pretty(range( maC2H6)))      # Add second axis
mtext("C2H6", side = 4, line = 3)

library("zoo")
medch4 <- rollapply(data = SC_20200406_dat$CH4, width = 102, FUN = median,align="center")
medc2h6 <- rollapply(data = SC_20200406_dat$C2H6, width = 102, FUN = median,align="center")

par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(medch4, pch = 16, col = 'black',ylab = "Moving Median CH4",
     xlab='Observation Number')              # Create first plot
par(new = TRUE)                             # Add new plot
plot(medc2h6,  pch = 17, col = 'blue',              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "",main='Comparison of Background Levels, CH4, C2H6')
axis(side = 4, at = pretty(range( maC2H6)))      # Add second axis
mtext("C2H6", side = 4, line = 3)

plot(SC_20200406_dat$C2C1,type='l')

nrow(SC_20200406_dat[SC_20200406_dat$C2C1==-1,])/nrow(SC_20200406_dat)
plot(SC_20200406_dat[SC_20200406_dat$C2C1!=-1,]$C2C1)
summary(SC_20200406_dat$C2C1)
summary(SC_20200406_dat[SC_20200406_dat$C2C1!=-1,]$C2C1)



plot(SC_20200406_dat$R)
summary(SC_20200406_dat$R)

woo <- SC_20200406_dat$CH4/SC_20200406_dat$C2H6
