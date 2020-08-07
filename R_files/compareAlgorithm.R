perc <- c(1,5,10,15)

results <- data.frame(
  perc = c(rep(c(1,2,3,4,5,10,15),4)),
  car = c(rep(c('CSU'),7),
          rep(c('CSU'),7),
          rep(c('SC'),7),
          rep(c('SC'),7)
          ),
  algorithm = c(rep(c("SC"),7),
                rep(c("CSU"),7),
                rep(c("SC"),7),
                rep(c("SC"),7)),
  VP = c(224,NA,NA,NA,16,7,4,
         NA,NA,NA,NA,23,8,3,
         11,2,1,1,0,0,0,
         17,4,1,NA,0,NA,NA
         
         ),
  OP = c(448,NA,NA,NA,70,32,22,
         NA,NA,NA,NA,112,48,29,
         70,39,19,10,11,1,1,
         101,40,20,NA,8,NA,NA
         ),
  backgroundObsSec = c(rep(c(102),7),
                       rep(c(102),7),
                       rep(c(102),7),
                       rep(c(500),7)
                       )
)

resultsSC <- data.frame(
  perc = c(rep(c(1,2,3,4,5,10,15),5)),
  car = c(rep(c('SC'),7),
          rep(c('SC'),7),
          rep(c('SC'),7),
          rep(c('SC'),7),
          rep(c('SC'),7)
          
  ),
  algorithm = c(rep(c("CSU"),7),
                rep(c("CSU"),7),
                rep(c("CSU"),7),
                rep(c("CSU"),7),
                rep(c("CSU"),7)   ),
  
  VP = c(8,2,1,0,0,0,0,
         8,2,1,0,0,0,NA,
         9,3,1,1,0,0,0,
         14,4,1,1,0,0,0,
         19,4,1,1,1,0,0  ),
  OP = c(67,39,19,10,8,1,1,
         73,37,19,10,8,10,NA,
         74,40,20,10,8,1,1,
         83,39,20,10,8,1,1,
         102,43,25,11,8,1,1
  ),
  backgroundObsSec = c(rep(c(102),7),
                       rep(c(51),7),
                       rep(c(204),7),
                       rep(c(500),7),
                       rep(c(1020),7)
  )
)

comparisonResults <- read_csv("~/Documents/GitHub/CSU_SC/comparisonResults.csv")

comparisonResults %>% 
  #filter(MAXSPEED == 45)  %>% 
  filter(CAR == "CSU") %>% 
  filter(BaselinePercentile == 50) %>% 
  #filter(BaselinePercentile == 50) %>% 
  ggplot(aes(x=PERCENTILE,y=VP,linetype=ALGORITHM,col=as.factor(MAXSPEED),shape = as.factor(MAXSPEED))) + 
  geom_point() + geom_line() + ggtitle("# VP (CSU Car) vs. threshold percentile")

comparisonResults %>% 
  #filter(MAXSPEED == 45)  %>% 
 # filter(CAR == "CSU") %>% 
  filter(BaselinePercentile == 50) %>% 
  #filter(VP < 100) %>% 
  #filter(BaselinePercentile == 50) %>% 
  ggplot(aes(x=PERCENTILE,y=VP,col=ALGORITHM,linetype=as.factor(MAXSPEED),shape = as.factor(MAXSPEED))) + 
  geom_point() + geom_line() + ggtitle("# VP vs. threshold percentile") + facet_wrap(~CAR)+
  xlim(c(0,15))+theme_bw()

comparisonResults %>% 
  #filter(MAXSPEED == 45)  %>% 
  filter(CAR == "SC") %>% 
  #filter(BaselinePercentile == 50) %>% 
  #filter(VP < 100) %>% 
  #filter(BaselinePercentile == 50) %>% 
  ggplot(aes(x=PERCENTILE,y=VP,col=ALGORITHM,linetype=as.factor(MAXSPEED),shape = as.factor(BaselinePercentile))) + 
  geom_point() + geom_line() + ggtitle("# VP vs. threshold percentile") + facet_wrap(~BaselinePercentile)+
  xlim(c(0,15))+theme_bw()

comparisonResults %>% 
  #filter(MAXSPEED == 45)  %>% 
  #filter(CAR == "SC") %>% 
  #filter(BaselinePercentile == 50) %>% 
  #filter(VP < 100) %>% 
  #filter(BaselinePercentile == 50) %>% 
  filter(OP < 300) %>% 
  ggplot(aes(x=PERCENTILE,y=OP,col=ALGORITHM,linetype=as.factor(MAXSPEED),shape = as.factor(BaselinePercentile))) + 
  geom_point() + geom_line() + ggtitle("# OP vs. threshold percentile") + facet_grid(CAR~BaselinePercentile)+
  xlim(c(0,15))+theme_bw()

comparisonResults %>% 
  #filter(MAXSPEED == 45)  %>% 
  #filter(CAR == "SC") %>% 
  #filter(BaselinePercentile == 50) %>% 
  #filter(VP < 100) %>% 
  #filter(BaselinePercentile == 50) %>% 
  filter(OP < 300) %>% 
  ggplot(aes(x=PERCENTILE,y=VP,col=ALGORITHM,linetype=as.factor(MAXSPEED),shape = as.factor(BaselinePercentile))) + 
  geom_point() + geom_line() + ggtitle("# VP vs. threshold percentile") + facet_grid(CAR~BaselinePercentile)+
  xlim(c(0,15))+theme_bw()

csu <- comparisonResults %>% 
  filter(CAR == "CSU" ) 
  #filter(PERCENTILE == '15')

sc <- comparisonResults %>% 
  filter(CAR == "SC" ) %>% 
  filter(PERCENTILE < 15) %>% 
  mutate(PERCENTILE = PERCENTILE + 5)

woot <- plyr::join(csu,sc,type='full')

woot %>% 
  ggplot(aes(x=PERCENTILE,y=OP,col = CAR, shape = ALGORITHM,linetype=as.factor(BaselinePercentile))) + geom_point() + geom_line() + facet_wrap(~MAXSPEED)

  mutate(keep = ifelse(CAR == "CSU" && PERCENTILE == 15,T,
                       ifelse(CAR == "SC" && PERCENTILE < 10,T,F
                       ))) %>% 
  filter(keep == T)

comparisonResults %>% 
  #filter(MAXSPEED == 45)  %>% 
  filter(CAR == "SC") %>% 
  filter(BaselinePercentile == 50) %>% 
  filter(PERCENTILE > 4) %>% 
  #filter(BaselinePercentile == 50) %>% 
  ggplot(aes(x=PERCENTILE,y=OP,linetype=as.factor(ALGORITHM),col = as.factor(OBSERVATIONSEC))) + 
  geom_point() + geom_line() + ggtitle("# OP (SC Car) vs. threshold percentile") + facet_wrap(~MAXSPEED)

comparisonResults %>% 
  #filter(MAXSPEED == 45)  %>% 
  filter(CAR == "SC") %>% 
  filter(ALGORITHM == 'SC') %>% 
  #filter(BaselinePercentile == 50) %>% 
  #filter(PERCENTILE > 4) %>% 
  #filter(BaselinePercentile == 50) %>% 
  ggplot(aes(x=PERCENTILE,y=OP,linetype=as.factor(MAXSPEED),col = as.factor(OBSERVATIONSEC))) + 
  geom_point() + geom_line() + 
  ggtitle("# OP (SC Car) vs. threshold percentile") + 
  facet_grid(OBSERVATIONSEC~BaselinePercentile) +
  theme_bw()

comparisonResults %>% 
  #filter(MAXSPEED == 45)  %>% 
  filter(CAR == "SC") %>% 
  filter(ALGORITHM == 'SC') %>% 
  #filter(BaselinePercentile == 50) %>% 
  #filter(PERCENTILE > 4) %>% 
  #filter(BaselinePercentile == 50) %>% 
  ggplot(aes(x=PERCENTILE,y=VP,linetype=as.factor(MAXSPEED),col = as.factor(OBSERVATIONSEC))) + 
  geom_point() + geom_line() + 
  ggtitle("# VP (SC Car) vs. threshold percentile") + 
  facet_grid(OBSERVATIONSEC~BaselinePercentile) +
  theme_bw()



comparisonResults %>% 
  filter(MAXSPEED == 45)  %>% 
  filter(CAR == "CSU") %>% 
  #filter(BaselinePercentile == 50) %>% 
  ggplot(aes(x=PERCENTILE,y=VP,linetype=ALGORITHM,col=CAR,shape=as.factor(OBSERVATIONSEC))) + geom_point() + geom_line()

resultsSC %>% 
  #mutate(perc = ifelse(car == 'SC',perc + 4,perc)) %>% 
  #filter(perc != 1) %>% 
  filter(!is.na(VP)) %>% 
  ggplot(aes(x=perc,y=VP,linetype = algorithm,col=car,shape=as.factor(backgroundObsSec))) + geom_point() + geom_line() + xlab("% Above Baseline") + 
  ylab("# VP")

resultsSC %>% 
  #filter(perc != 1) %>% 
  #mutate(perc = ifelse(car == 'SC',perc + 4,perc)) %>% 
  
  filter(!is.na(VP)) %>% 
  
  ggplot(aes(x=perc,y=OP,linetype = algorithm,col=car,shape=as.factor(backgroundObsSec))) + geom_point() + geom_line() + xlab("% Above Baseline") + 
  ylab("# OP")


resultsSC %>% 
  filter(!is.na(VP)) %>% 
  ggplot(aes(x=perc,y=VP,col=as.factor(backgroundObsSec),shape=as.factor(backgroundObsSec))) + geom_point() + geom_line() + xlab("% Above Baseline") + 
  ylab("# OP") + ggtitle( "VP, with varying threshold and Background")


results %>% 
  filter(!is.na(VP)) %>% 
  filter(car == 'CSU') %>% 
  ggplot(aes(x=perc,y=VP,linetype = algorithm,col=as.factor(backgroundObsSec),shape=as.factor(backgroundObsSec))) + geom_point() + geom_line() + xlab("% Above Baseline") + 
  ylab("# OP") + ggtitle( "VP, with varying threshold and Background")

scRes <- results %>% filter(car!='CSU')
togRes <- plyr::join(results,resultsSC, type="full")


togRes %>% 
  filter(!is.na(VP)) %>% 
  #filter(OP <50) %>% 
  #filter(car == 'CSU') %>% 
  ggplot(aes(x=perc,y=OP,linetype = algorithm,col=as.factor(backgroundObsSec),shape=car)) + geom_point() + geom_line() + xlab("% Above Baseline") + 
  ylab("# OP") + ggtitle( "OP, with varying threshold and Background") + ylim(c(0,120))


