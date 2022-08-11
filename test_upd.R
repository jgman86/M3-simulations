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
cmd_path<-paste0("C:/Coding/cmdstan-2.30.0")
set_cmdstan_path(path=cmd_path)

m3_upd <- cmdstan_model("Models/M3_Updating_LKJ_NC.stan")


# Source Functions
source("Functions/M3_functions.R")

# Sim Conditions
SampleSize = 100
N <- 5
K <- 16
nFT <- 4
nRetrievals <- 1000

# Timing Conditions
fixtime <- 0.25
enctime <- 0.5
minTime <- 0.25
maxTime <- 1.75


# Set Range for Parameter Means
range_muC <- c(1,10)
range_muA <- c(0,0.5)
range_muD <- c(0.5,1)
range_mueU <-c(0.2,0.6)
range_muR <- c(0,10)
eta <- 5 # Simulated N = 10000 with eta = 5, 95 % of all values lie witin 0 -+ 0.56


sigC <- c(0.125,0.5)
sigA <- c(0.125,0.5)
sigD<- c(0.2,0.8)
sigE <- c(1,3)
sigR <- c(0,1.5)
sigB <- c(0.0001, 0.1)

# Create Test Data 

conWCI <- seq(from = minTime, to = maxTime, length.out = nFT)
conCWI <- seq(from = minTime, to = maxTime, length.out = nFT)


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
                            lower=c(0,0,-Inf,0,0,0),upper = c(Inf,Inf,Inf,Inf,Inf,Inf))
# Merge Parameters to one Matrix
colnames(parms) <- c("conA","genA","d","eU","r","baseA")
parms[,3] <- 1 / (1+exp(-parms[,3]))

parms[,6] <- 0.1

# Simulate Data 
ParmsUpd <- matrix(rep(parms,each = length(conCWI)*length(conWCI)), nrow = nrow(parms)*length(conCWI)*length(conWCI), ncol = ncol(parms), byrow = F)

colnames(ParmsUpd) <- c("conA","genA","d","eU","r","baseA")

WordCue <- rep(conWCI,length.out = nrow(parms))
CueWord <-  rep(conCWI,length.out = nrow(parms))
fix <- rep(fixtime, length.out = nrow(parms))
enc <- rep(enctime, length.out = nrow(parms))

data <- simData_UpdatingModel(ParmsUpd,as.vector(respOpt_Cspan(N,K)),nRetrievals,conCWI,conWCI,fixtime,enctime)

# Create Stan Data 
stan.dat <- list(count = data[,2:6], 
                 K = 5,
                 R = as.vector(respOpt_Cspan(N,K)),
                 J = length(sigs)-1,
                 N = length(unique(data[,"ID"])),
                 Con1 = length(unique(data[,"t_EU"])),
                 Con2 = length(unique(data[,"t_rm"])),
                 t_eU = data[,"t_EU"],
                 t_rm = data[,"t_rm"],
                 retrievals = nRetrievals,
                 scale_b = 0.1)


m3_upd <- m3_upd$sample(stan.dat,chains=4,
                        parallel_chains = 4,
                        iter_sampling =1500,
                        iter_warmup=1500,
                        show_messages = F,
                        #adapt_delta=0.95,
                        max_treedepth = 15,
                        init=init_nc)

# Plot some subj theta results
mcmc_combo(m3_upd$draws(c("subj_pars[1,3]","subj_pars[2,3]","subj_pars[4,3]","subj_pars[5,3]","f[3]")))
# Plot hyper pars
mcmc_combo(m3_upd$draws(c("hyper_pars","mu_d")))

