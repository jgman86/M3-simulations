library(data.table) 
library(tidyverse,quietly = T)

# Read in Results 
files <- list.files(pattern = "M3_EE_CS")

df <- files %>% map(readRDS) %>% rbindlist(fill = T)


# Plot Results 

ggplot(df, aes(y=cor_e, x=OtherItems, group=as.factor(NPL), color=as.factor(NPL))) +
  geom_line() + geom_point() + facet_grid(~nFreetime + Retrievals)

       