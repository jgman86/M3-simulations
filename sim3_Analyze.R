library(data.table) 
library(tidyverse,quietly = T)

# Read in Results 
files <- list.files(pattern = "M3_EE_CS",path = "Results/",full.names = T)

df <- files %>% map(readRDS) %>% rbindlist(fill = T)

roi <- df %>%filter(Retrievals==500) 

df_full <- roi %>% filter(fixedf ==0)
df_fixed <- df %>% filter(fixedf ==1)

# Cor ranges
# full model
roi %>%  filter(Retrievals==500, fixedf==0) %>% select("cor_c") %>% range()
roi %>%  filter(Retrievals==500, fixedf==0) %>% select("cor_a") %>% range()
roi %>%  filter(Retrievals==500, fixedf==0) %>% select("cor_f") %>% range()
roi %>%  filter(Retrievals==250, fixedf==0) %>% select("cor_e") %>% range()
roi %>%  filter(Retrievals==250, fixedf==0) %>% select("cor_r") %>% range()

# fixed f

roi %>%  filter(Retrievals==500, fixedf==1) %>% select("cor_c") %>% range()
roi %>%  filter(Retrievals==500, fixedf==1) %>% select("cor_a") %>% range()
roi %>%  filter(Retrievals==500, fixedf==1) %>% select("cor_e") %>% range()
roi %>%  filter(Retrievals==500, fixedf==1) %>% select("cor_r") %>% range()




# Plot Results 
 
ft.labs <- c("FreeTime Conditions:2","Freetime Conditions:4")
names(ft.labs) <- c(2,4)

f.labs <- c("Full Model","Fixed f Parameter")
names(f.labs) <- c(0,1)

# Ran
# New facet label names for supp variable
# c parameter
ggplot(roi, aes(y=cor_c, x=as.factor(NPL), group=as.factor(OtherItems), color=as.factor(OtherItems))) +
  geom_line() + geom_point() + facet_grid(nFreetime~ fixedf,
                                          labeller=labeller(nFreetime=ft.labs, fixedf=f.labs)) + 
  theme_bw() + labs(color="Distractor Items in Other Positions",x= "NPLs", y=expression(italic("r"))) + 
  scale_color_viridis_d()  + theme(legend.position = "top")
# a parameter
ggplot(roi, aes(y=cor_a, x=as.factor(NPL), group=as.factor(OtherItems), color=as.factor(OtherItems))) +
  geom_line() + geom_point() + facet_grid(nFreetime~ fixedf,
                                          labeller=labeller(nFreetime=ft.labs, fixedf=f.labs)) + 
  theme_bw() + labs(color="Distractor Items in Other Positions",x= "NPLs", y=expression(italic("r"))) +
  scale_color_viridis_d() + theme(legend.position = "top")
# e parameter
ggplot(roi, aes(y=cor_e, x=as.factor(NPL), group=as.factor(OtherItems), color=as.factor(OtherItems))) +
  geom_line() + geom_point() + facet_grid(nFreetime~ fixedf,
                                          labeller=labeller(nFreetime=ft.labs, fixedf=f.labs)) + 
  theme_bw() + labs(color="Distractor Items in Other Positions",x= "NPLs", y=expression(italic("r")))  + 
  scale_color_viridis_d() + theme(legend.position = "top")
#r parameter
ggplot(roi, aes(y=cor_r, x=as.factor(NPL), group=as.factor(OtherItems), color=as.factor(OtherItems))) +
  geom_line() + geom_point() + facet_grid(nFreetime~ fixedf,
                                          labeller=labeller(nFreetime=ft.labs, fixedf=f.labs)) + 
  theme_bw() + labs(color="Distractor Items in Other Positions",x= "NPLs", y=expression(italic("r")))  + 
  scale_color_viridis_d() + theme(legend.position = "top")
#r parameter


roi %>% filter(fixedf==0, Retrievals==500) %>%
ggplot(., aes(y=cor_f, x=as.factor(NPL), group=as.factor(OtherItems), color=as.factor(OtherItems))) +
  geom_line() + geom_point() + facet_grid(~nFreetime,labeller=labeller(nFreetime=ft.labs)) + 
  theme_bw() + labs(color="Distractor Items in Other Positions",x= "NPLs", y=expression(italic("r")))  + 
  scale_color_viridis_d() + theme(legend.position = "top")
#r parameter






