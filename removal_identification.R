library(ggplot2)
library(plyr)
library(tidyverse)


r <- seq(1,20,1)
nFT <- c(2,4)

theta <- tidyr::crossing(r,nFT)


removal <- function(r,nFT,minFT,maxFT){
  
  ft <- round(seq(minFT,maxFT,length.out=nFT),1)
  ft <- append(ft, c(0.1,0.2), after = 1)
  
  data.frame(ft=ft,nFT=nFT, value = exp(-r*ft))
}

all <- mdply(theta, removal,minFT=0,maxFT=2)



ggplot(all, aes(r, value, colour=factor(ft))) + geom_line() +
  geom_point() + theme_bw() + scale_x_continuous(breaks = seq(1,20,by=1)) + scale_y_continuous(breaks=seq(0,1,by=.05))+
  facet_wrap(~nFT) + labs(y="conA",colour="FT",title = "Activation change due removal for different freetime conditions") + 
  geom_vline(xintercept = 8.5,size=.75, color="blue",alpha=0.2) 


