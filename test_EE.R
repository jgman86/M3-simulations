suppressPackageStartupMessages({
  library(optparse) # to create a command line interface
  library(SimDesign)
  library(bayesplot)
  library(cmdstanr)
  library(tidyverse)
  library(data.table)
  library(psych)
  library(here)
  library(posterior)
  library(tidybayes)
  library(HDInterval)})

source("Functions/M3_functions.R")

init <- function()
{
  list(hyper_pars=cbind(runif(stan1.dat$J,10,20)),
       subj_pars=cbind(runif(stan1.dat$N,1,10)))
  
}


# Set cmdstan path
cmd_path<-paste0("C:/Coding/cmdstan-2.30.0")
set_cmdstan_path(path=cmd_path)

# Load Models
mod3_EE <- cmdstan_model("Models/M3_ComplexSpan_EE_LKJ_Cholesky_NC.stan")
mod_univ_nc_ff <- cmdstan_model("Models/M3_ComplexSpan_EE_nc_fixedf.stan")

# Simulate Data for Estimation ----

nRetrievals <- 500
minFT <- 0.2
maxFT <- 1
nFT <- c(2,4) # 2,4,10 Conditions between 0.2 and 2
SampleSize <- 100 
con_nFT =4
#conFT <- seq(from = minFT, to = maxFT, length.out = con_nFT) # eventuell log scale 0.2 0.8 2.4 oder so?
conFT <- c(0,0.2,0.4,0.8)
# Set Range for Parameter Means
range_muC <- c(1,30)
range_muA <- c(0,0.5)
range_muF <- c(0,1) # fix to 0.5
range_muE <-c(0,0.5)
range_muR <- c(0,25) # range 0 - 25 empirical derived
eta <- 5 # Simulated N = 10000 with eta = 5, 95 % of all values lie within 0 -+ 0.56

sigC <- c(0.125,0.5)
sigA <- c(0.125,0.5)
sigF <- c(0.1,0.5)
sigE <- c(1,2)
sigR <- c(0.125,0.5) # abhänbgig von removal parameter -> analog zu c und a
sigB <- c(0.0001, 0.1)


# Sample Hyper-Parameter Means with C as fixpoint ----
relCA <- runif(1, min = range_muA[1],max = range_muA[2])
Mean_Cpar <- 8
Mean_Apar <- 4
Mean_Epar <- 0.4
Mean_Rpar <- 25
Mean_Fpar <-  0.6
log_mu_f <- log(Mean_Fpar/(1-Mean_Fpar))
Mean_bpar <- 0.1


# Make Vector with Hyper Pars

hyper_mus <- c(Mean_Cpar,Mean_Apar,log_mu_f,Mean_Epar,Mean_Rpar, Mean_bpar)


# Sample Variances and Set Covariances----

sig_c <- Mean_Cpar*runif(1, min = sigC[1], max= sigC[2])
sig_a <- Mean_Apar*runif(1, min = sigA[1], max= sigA[2])
sig_f <- Mean_Fpar*runif(1, min = sigF[1], max= sigF[2])
sig_e <- Mean_Epar*runif(1, min = sigE[1], max= sigE[2])
sig_r <- Mean_Rpar*runif(1, min = sigR[1], max= sigR[2])
sig_b <- 0.001


sigs <-c(sig_c,sig_a,sig_f,sig_e,sig_r,sig_b)
Sig <- diag(length(hyper_mus))

Sig[1,1] <- (sig_c)^2
Sig[2,2] <- (sig_a)^2
Sig[3,3] <- (sig_f)^2
Sig[4,4] <- (sig_e)^2
Sig[5,5] <- (sig_r)^2
Sig[6,6] <- (sig_b)^2


# Set Correlations for Parameters ----

# Sample Covariance Matrix Sigma

omega <- rlkjcorr(1,length(hyper_mus),eta)

# Little Hack for fixing coveriance of b to zer0

omega[6,1:5] = omega[1:5,6] = 0 #fix cov of b to 0
#omega[3,1:5] = omega[1:5,3] = .00001 #fix cov of f to 0


Sigma <- cor2cov(omega,sigs)


# Sample Parameters from MVN ----

parms <- tmvtnorm::rtmvnorm(n=SampleSize, mean= hyper_mus, sigma=Sigma,
                            lower=c(0,0,-Inf,0,0,0),upper = c(Inf,Inf,Inf,Inf,Inf,Inf))
# Merge Parameters to one Matrix
colnames(parms) <- c("conA","genA","f","e","r","baseA")

parms[,6] <- 0.1
#parms[,3] <- 0
parms[,3] <- 1 / (1+exp(-parms[,3]))

# Simulate Data for Estimation ----

ParmsFT <- matrix(rep(parms,each =con_nFT), nrow = length(parms[,1])*con_nFT, ncol = ncol(parms), byrow = F)
colnames(ParmsFT) <- c("conA","genA","f","e","r","baseA")
FT <- rep(conFT,length.out = nrow(ParmsFT))

data <- simData_CSpanEE(ParmsFT,as.vector(respOpt_Cspan(4,16)),1000,FT)


# Generate Stan Data ----


stan1.dat <- list(count = data[,4:8], 
                 K = 5,
                 R = as.vector(respOpt_Cspan(4,16)),
                 J = length(sigs)-1,
                 N = length(unique(data[,"ID"])),
                 Con = length(unique(data[,"Freetime"])),
                 Freetime = unique(data[,"Freetime"]),
                 retrievals = 1000,
                 scale_b = 0.1)


fit <- mod3_EE$sample(stan1.dat,refresh = 100,
                      init = init,
                      chains = 4,
                      parallel_chains = 4,
                      iter_sampling = 1500,
                      iter_warmup = 1500,
                      show_messages = F)



# summarise results

hyper <- fit$summary(c("hyper_pars","mu_f"), mean,sd,Mode,rhat) %>% pivot_wider(.,names_from = "variable",values_from = c("mean","Mode","sd","rhat")) 

subj <-  fit$summary(c("subj_pars"),mean,sd,Mode,rhat)  %>% mutate(variable = str_remove_all(variable, "subj_pars")) %>% 
  separate(col = variable,into = c("theta","ID"),sep = ",")  %>% 
  mutate(ID = str_remove(ID,pattern = "]"), theta = case_when(theta == "[1" ~ "c", theta == "[2" ~ "a",
                                                              theta == "[3" ~ "log_f",
                                                              theta == "[4" ~ "EE",
                                                              theta== "[5" ~ "r")) %>%
  pivot_wider(.,names_from = "theta",values_from = c("mean","Mode","sd","rhat"))

f <-  fit$summary("f",mean,sd,Mode,rhat) %>% mutate(ID = seq(1:SampleSize), theta = "f") %>% relocate(c("ID","theta"), .before = variable) %>% select(-variable) %>%
  pivot_wider(.,names_from = "theta",values_from = c("mean","Mode","sd","rhat")) 

# Recoveries of Subject Pars
cor(parms[,"conA"],subj$mean_c)
cor(parms[,"genA"],subj$mean_a)
cor(parms[,"f"],f$mean_f)
cor(parms[,"e"],subj$mean_EE)
cor(parms[,"r"],subj$mean_r)

# max rhat for subj pars

max_rhat_c <- subj %>% select(contains("rhat_c")) %>% max()
max_rhat_a <- subj %>% select(contains("rhat_a")) %>% max()
max_rhat_f <- subj %>% select(contains("rhat_f")) %>% max()
max_rhat_e <- subj %>% select(contains("rhat_e")) %>% max()
max_rhat_r <- subj %>% select(contains("rhat_r")) %>% max()




# Stan is noisy, so tell it to be more quiet()
M3 <-  mod3_norm$sample(stan1.dat,
                  refresh = 100,
                  chains = 4,
                  parallel_chains=4, 
                  iter_warmup=500,
                  iter_sampling=1000,
                  adapt_delta=.99,
                  max_treedepth=15,
                  init=init_test,
                  show_messages = F)

fit$summary(c("hyper_pars","mu_f"), mean,sd,Mode,rhat)

fit$summary(c("subj_pars[1,1]"), mean)
    
    
    
    
    [,2],parms[,"conA"])

mcmc_combo(fit$draws(c("hyper_pars","mu_f")))


plot((exp(-seq(10:25)*1.25)*9.2),type="b")

plot(((exp(-seq(1:25)*1))*9.2),type="b")


plot(((exp(-seq(1:10)*c(0,0.2,0.4,0.8)))*1))





g<-Generate_M3(sim3,fo)
r<-Analyze_M3(sim3,g,fo)
s<-Summarise(sim3,r,fo)






