suppressPackageStartupMessages({
  library(optparse) # to create a command line interface
  library(SimDesign)
  library(cmdstanr)
  library(tidyverse)
  library(data.table)
  library(psych)
  library(here)
  library(posterior)
  library(tidybayes)
  library(HDInterval)
  library(bayesplot)})

# Set cmdstan path
#cmd_path<-paste0("C:/Coding/cmdstan-2.30.0")
#set_cmdstan_path(path=cmd_path)

m3_upd <- cmdstan_model("Models/M3_Updating_LKJ_NC.stan")


# Source Functions
source("Functions/M3_functions.R")

# Sim Conditions
SampleSize = 50
N <- 5
K <- 16
nFT <-2
nRetrievals <- 1000

# Timing Conditions
fixtime <- 0.2
enctime <- 0.5
minTime <- 0.1
maxTime <- 0.8

conWCI <- seq(from = minTime, to = maxTime, length.out = nFT)
conCWI <- seq(from = minTime, to = maxTime, length.out = nFT)


# Set Range for Parameter Means
range_muC <- c(1,10)
range_muA <- c(0,0.5)
range_muD <- c(0,0.8)
range_mueU <-c(0,5)
range_muR <- c(0.2,0.8)
eta <- 5 # Simulated N = 10000 with eta = 5, 95 % of all values lie witin 0 -+ 0.56


sigC <- c(0.125,0.5)
sigA <- c(0.125,0.5)
sigD<- c(0.2,0.8)
sigE <- c(1,3)
sigR <- c(0,5)
sigB <- c(0.0001, 0.1)

# Create Test Data 



relCA <- runif(1, min = range_muA[1],max = range_muA[2])
Mean_Cpar <- runif(1, min =range_muC[1], max= range_muC[2])
Mean_Apar <- Mean_Cpar*(relCA)
Mean_Dpar <- runif(1, min =range_muD[1], max = range_muD[2])
Mean_EUpar <- runif(1, min =range_mueU[1], max = range_mueU[2]) #add new vars to fixed objects 
Mean_Rpar <- runif(1, min= range_muR[1], max = range_muR[2])
log_mu_d <- log(Mean_Dpar/(1-Mean_Dpar))
Mean_bpar <- 0.1
n_theta = 5

# Make Vector with Hyper Pars

hyper_mus <- c(Mean_Cpar,Mean_Apar,log_mu_d,Mean_EUpar,Mean_Rpar, Mean_bpar)

# Set Variances 
sig_c <- runif(1, min = sigC[1],max = sigC[2])*Mean_Cpar
sig_a <- runif(1, min = sigA[1],max = sigA[2])*Mean_Apar 
sig_d <- runif(1, min = sigD[1],max = sigD[2])
sig_eu <- runif(1, min = sigE[1], max= sigE[2])
sig_r <- runif(1, min = sigR[1], max= sigR[2])*Mean_Rpar
sig_b <- 0.001

sigs <-c(sig_c,sig_a,sig_d,sig_eu,sig_r,sig_b)
Sig <- diag(length(hyper_mus))

Sig[1,1] <- (sig_c)^2
Sig[2,2] <- (sig_a)^2
Sig[3,3] <- (sig_d)^2
Sig[4,4] <- (sig_eu)^2
Sig[5,5] <- (sig_r)^2
Sig[6,6] <- (sig_b)^2


# Set Correlations for Parameters ----

# Sample Covariance Matrix Sigma

omega <- rlkjcorr(1,length(hyper_mus),eta)

# Little Hack for fixing coveriance of b to zer0

omega[6,1:5] = omega[1:5,6] = 0

Sigma <- cor2cov(omega,sigs)


# Sample Parameters from MVN ----

parms <- tmvtnorm::rtmvnorm(n=SampleSize, mean= hyper_mus, sigma=Sigma,
                            lower=c(0,0,-Inf,0,0,0),
                            upper = c(Inf,Inf,Inf,Inf,Inf,Inf))


# Merge Parameters to one Matrix
colnames(parms) <- c("conA","genA","d","eU","r","baseA")
parms[,3] <- 1 / (1+exp(-parms[,3]))

parms[,6] <- 0.1

# Simulate Data 
ParmsUpd <- matrix(rep(parms,each = length(conCWI)*length(conWCI)), 
                   nrow = nrow(parms)*length(conCWI)*length(conWCI), 
                   ncol = ncol(parms), byrow = F)

colnames(ParmsUpd) <- c("conA","genA","d","eU","r","baseA")

# WordCue <- rep(conWCI,length.out = nrow(parms))
# CueWord <-  rep(conCWI,length.out = nrow(parms))
# fix <- rep(fixtime, length.out = nrow(parms))
# enc <- rep(enctime, length.out = nrow(parms))
t_EU <- conWCI + fixtime + conCWI
t_rm <- fixtime + conCWI + enctime + conWCI

simData_UpdatingModel_test <- function(parmsMMM,respOpts,nRetrievals,CWI,WCI,fixtime,enctime){
  # extract parms
  
  
  if(is.vector(parmsMMM)){
    conA <- parmsMMM["conA"]
    genA <- parmsMMM["genA"]
    EU <- parmsMMM["eU"]
    rm <- parmsMMM["r"]
    d <- parmsMMM["d"]
    baseA <- parmsMMM["baseA"]
    
  }else{
    conA <- parmsMMM[,"conA"]
    genA <- parmsMMM[,"genA"]
    EU <- parmsMMM[,"eU"]
    rm <- parmsMMM[,"r"]
    d <- parmsMMM[,"d"]
    baseA <- parmsMMM[,"baseA"]
  }
  
  Cons <- length(unique(CWI))*length(unique(WCI))
  
  # Compute Extended Updating and Removal Time
  t_EU <- WCI + fixtime + CWI
  t_rm <- fixtime + CWI + enctime + WCI
  
  Con1 <- length(t_EU)
  Con2 <- length(t_rm)
  
  acts <- matrix(NaN, ncol = 7,nrow =nrow(parmsMMM)*Cons)
  
  c <- 1
  
  for(j in 1:nrow(parmsMMM))
    for (i in 1:Con1)
      for(k in 1:Con2)
        
      {
        
        acts[c,1] <- 0.1 + (1+EU[j]*t_EU[i])*(conA[j] + genA[j])
        acts[c,2] <- 0.1 + (1+EU[j]*t_EU[i])*genA[j]
        acts[c,3] <- 0.1 + (exp(-rm[j]*t_rm[k])*d[j]*(1+EU[j]*t_EU[i])*conA[j]) + (1+EU[j]*t_EU[i])*genA[j]
        acts[c,4] <- 0.1 + (1+EU[j]*t_EU[i])*genA[j]
        acts[c,5] <- 0.1
        acts[c,6] <-t_EU[i]
        acts[c,7] <-t_rm[k]
        # 
        c <- c + 1 
      }
  
  
  colnames(acts) <- c("A_IIP","A_IOP","A_DIP","A_DIOP","A_NPL","t_EU","t_rm")
  
  # summarize activations
  # if(length(A_IIP) != 1){
  #   acts <- cbind(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  # }else{
  #   acts <- c(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  # }
  
  # compute summed activation
  sumActs <- apply(t(respOpts*t(acts[,1:5])),1,sum)
  
  Probs <- t(respOpts*t(acts[,1:5]))/sumActs
  
  colnames(Probs) <- c("P_IIP","P_IOP","P_DIP","P_DOP","P_NPL")
  
  simdata <- matrix(NA,ncol = ncol(Probs),nrow = nrow(Probs))
  ID <- matrix(NA, ncol=1, nrow=nrow(Probs))
  i <- 1
  
  for(id in 1:nrow(Probs)){
    
    simdata[id,] <- t(stats::rmultinom(1,nRetrievals,Probs[id,]))
    
    
    ID[id,] <- i
    
    if (id %% Cons == 0)
    {
      
      i <- i+1
      
    }
    
    
  }
  
  

  data <- cbind(ID,simdata,acts[,6:7])
  colnames(data) <- c("ID","IIP","IOP","DIP","DOP","NPL","t_rm","t_EU")
  return(data)
}
data <- simData_UpdatingModel_test(parms,c(1,8,1,4,16),nRetrievals,conCWI,conWCI,fixtime,enctime)

resp <- c()
# Create Stan Data 
stan.dat <- list(count = data[,2:6], 
                 K = 5,
                 R = c(1,8,1,4,16),
                 J = length(sigs)-1,
                 N = length(unique(data[,"ID"])),
                 Con1 = nFT,
                 Con2 = nFT,
                 t_eU = data[,"t_EU"],
                 t_rm = data[,"t_rm"],
                 retrievals = nRetrievals,
                 scale_b = 0.1)


upd_fit <- m3_upd$sample(stan.dat,chains=4,
                         parallel_chains = 4,
                         iter_sampling =1500,
                         iter_warmup=1500,
                         show_messages = F,
                         #adapt_delta=0.95,
                         max_treedepth = 15,
                         init=init_nc)
upd_fit$output()
# Plot some subj theta results
mcmc_combo(upd_fit$draws(c("subj_pars[1,3]","subj_pars[2,3]","subj_pars[4,3]","subj_pars[5,3]","d[3]")))
# Plot hyper pars
mcmc_combo(upd_fit$draws(c("hyper_pars","mu_d")))


hyper <- upd_fit$summary(c("mu_d","hyper_pars"), mean, rhat) 
