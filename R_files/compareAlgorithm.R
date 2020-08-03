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
         17,NA,NA,NA,0,NA,NA
         
         ),
  OP = c(448,NA,NA,NA,70,32,22,
         NA,NA,NA,NA,112,48,29,
         70,39,19,10,11,1,1,
         101,NA,NA,NA,8,NA,NA
         ),
  backgroundObsSec = c(rep(c(102),7),
                       rep(c(102),7),
                       rep(c(102),7),
                       rep(c(500),7)
                       )
)


results %>% 
  mutate(perc = ifelse(car == 'SC',perc + 4,perc)) %>% 
  #filter(perc != 1) %>% 
  filter(!is.na(VP)) %>% 
  ggplot(aes(x=perc,y=VP,linetype = algorithm,col=car,shape=as.factor(backgroundObsSec))) + geom_point() + geom_line() + xlab("% Above Baseline") + 
  ylab("# VP")

results %>% 
  #filter(perc != 1) %>% 
  mutate(perc = ifelse(car == 'SC',perc + 4,perc)) %>% 
  
  filter(!is.na(VP)) %>% 
  
  ggplot(aes(x=perc,y=OP,linetype = algorithm,col=car,shape=as.factor(backgroundObsSec))) + geom_point() + geom_line() + xlab("% Above Baseline") + 
  ylab("# OP")
