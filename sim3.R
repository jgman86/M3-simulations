#####M3 Extended Encoding and Updating Simulations########

### Load Packages ####

library(SimDesign)
library(cmdstanr)
library(tidyverse)
library(data.table)
library(psych)
library(here)
library(posterior)


#### Set Options ####
dir_path <- here()
cmd_path <- paste0("C:/Coding/cmdstan-2.30.0")
set_cmdstan_path(cmd_path)
cmdstan_version()

### Source Functions ####
source(paste0(dir_path,"/Functions/M3_functions.R"))

#### Define Simulation Design and SetUp Model----

###### Varying Simulation Factors ---- 
N <- c(4)
K <- c(8,16)
nRetrievals <- c(500)
nFT<- c(2) # 2,4,10 Conditions between 0.2 and 2

###### Create Simulation Table ---- 
sim3 <- createDesign(OtherItems=N,
                     NPL=K,
                     Retrievals = nRetrievals,
                     nFreetime=nFT)

###### Fixed Simulation Factors ---- 
SampleSize <- 10
reps2con <- 100
minFT <- 0.5
maxFT <- 2

###### Model Path #####
stan_path_M3_EE <- paste0(dir_path,"/Models/M3_ComplexSpan_EE_LKJ_Cholesky_NC.stan")
#stan_path_M3_Upd <- paste0(dir_path,"/Models/M3_") add upadting model 


###### Compile and store within fixed objects for use in runSimulation() #####
fo <- list(M3_CS_EE=cmdstan_model(stan_path_M3_EE),
           SampleSize = SampleSize,
           minFT = minFT,
           maxFT = maxFT,
           # Set Range for Parameter Means
           range_muC  =c(1,100),
           range_muA = c(0,0.5),
           range_muF = c(0,1), # fix to 0.5
           range_muE = c(0,0.5),
           range_muR = c(0,25), # range 0 - 25 empirical derived
           eta = 5, # Simulated N = 10000 with eta = 5, 95 % of all values lie within 0 -+ 0.56
           sigC = c(0.125,0.5),
           sigA = c(0.125,0.5),
           sigF = c(0.1,0.5),
           sigE = c(1,2),
           sigR = c(0.125,0.5), # abhänbgig von removal parameter -> analog zu c und a
           sigB = c(0.0001, 0.1)
           
)


##### Set Up Fitting Function for cmdstan
stan_fit <- function(mod, dat){
  
  set_cmdstan_path(path="C:/Coding/cmdstan-2.30.0/")
  
  
  init <- function()
  {
    list(hyper_pars=cbind(runif(dat$J,10,20)),
         subj_pars=cbind(runif(dat$N,1,10)))
    
  }
  
  
  # Stan is noisy, so tell it to be more quiet()
  M3 <-  quiet(mod$sample(dat,
                         refresh = 0,
                         chains = 4,
                         #parallel_chains=4, 
                         iter_warmup=150,
                         iter_sampling=500,
                         adapt_delta=.95,
                         max_treedepth=15,
                         init = init,
                         show_messages = FALSE))
  
  M3 <- M3$summary("hyper_pars", mean)$mean[1:5]
  
  M3
}

##### Set Up Simulation Functions for Sim Design ----

Generate_M3 <- function(condition, fixed_objects=NULL) {
  Attach(condition)
  
  SampleSize <- fixed_objects$SampleSize
  minFT <- fixed_objects$minFT
  maxFT <- fixed_objects$maxFT
  range_muC  <- fixed_objects$range_muC
  range_muA<- fixed_objects$range_muA
  range_muF <- fixed_objects$range_muF
  range_muE <- fixed_objects$range_muE
  range_muR <- fixed_objects$range_muR
  eta <- fixed_objects$eta
  sigC <- fixed_objects$sigC
  sigA <- fixed_objects$sigA
  sigF <- fixed_objects$sigF
  sigE <- fixed_objects$range_muE
  sigR <-fixed_objects$range_muR
  sigB <- fixed_objects$sigB
  
  
  
  # Generate FreeTime Conditions
  
  conFT <- seq(from = minFT, to = maxFT, length.out = nFreetime) # eventuell log scale 0.2 0.8 2.4 oder so?
  
  
  # Sample Hyper-Parameter Means with C as fixpoint ----
  relCA <- runif(1, min = range_muA[1],max = range_muA[2])
  Mean_Cpar <- runif(1, min =range_muC[1], max= range_muC[2])
  Mean_Apar <- Mean_Cpar*(relCA)
  Mean_Epar <- runif(1, min =range_muE[1], max = range_muE[2])
  Mean_Rpar <- runif(1, min= range_muR[1], max = range_muR[2])
  Mean_Fpar <- runif(1, min= range_muF[1], max = range_muF[2])
  log_mu_f <- log(Mean_Fpar/(1-Mean_Fpar))
  Mean_bpar <- 0.1
  
  # Make Vector with Hyper Pars
  
  hyper_mus <- c(Mean_Cpar,Mean_Apar,log_mu_f,Mean_Epar,Mean_Rpar, Mean_bpar)
  
  
  # Sample Variances and Set Covariances----
  
  sig_c <- runif(1, min = sigC[1], max = sigC[2])*Mean_Cpar
  sig_a <- runif(1, min = sigA[1], max = sigA[2])*Mean_Apar 
  sig_f <- runif(1, min = sigF[1], max= sigF[2])
  sig_e <- runif(1, min = sigE[1], max= sigE[2])
  sig_r <- runif(1, min = sigR[1], max= sigR[2])*Mean_Rpar
  sig_b <- 0.001
  
  
  sigs <-c(sig_c,sig_a,sig_f,sig_e,sig_r,sig_b)
  Sig <- diag(length(hyper_mus))
  
  Sig[1,1] <- (sig_c)^2
  Sig[2,2] <- (sig_a)^2
  Sig[3,3] <- (sig_f)^2
  Sig[4,4] <- (sig_e)^2
  Sig[5,5] <- (sig_r)^2
  Sig[6,6] <- (sig_b)^2
  
  
  # Set Correlations fixed_objectsr Parameters ----
  
  # Sample Covariance Matrix Sigma
  
  omega <- rlkjcorr(1,length(hyper_mus),eta)
  
  # Little Hack fixed_objectsr fixing coveriance of b to zer0
  
  omega[6,1:5] = omega[1:5,6] = 0 #fix cov of b to 0
  #omega[3,1:5] = omega[1:5,3] = 0 #fix cov of f to 0
  
  
  Sigma <- cor2cov(omega,sigs)
  
  
  # Sample Parameters from MVN ----
  
  parms <- tmvtnorm::rtmvnorm(n=SampleSize, mean= hyper_mus, sigma=Sigma,
                              lower=c(0,0,-Inf,0,0,0),upper = c(Inf,Inf,Inf,Inf,Inf,Inf))
  # Merge Parameters to one Matrix
  colnames(parms) <- c("conA","genA","f","e","r","baseA")
  parms[,6] <- 0.1
  parms[,3] <- 1 / (1+exp(-parms[,3]))
  
  #parms[,3] <- 0 if f is fixed to 0.5
  
  
  # Simulate Data fixed_objectsr Estimation ----
  
  ParmsFT <- matrix(rep(parms,each =nFreetime), nrow = length(parms[,1])*nFreetime, ncol = ncol(parms), byrow = F)
  colnames(ParmsFT) <- c("conA","genA","f","e","r","baseA")
  FT <- rep(conFT,length.out = nrow(ParmsFT))
  
  data <- simData_CSpanEE(ParmsFT,as.vector(respOpt_Cspan(OtherItems,NPL)),Retrievals,FT)
  
  
  # Generate Stan Data ----
  
  
  dat <- list(count = data[,4:8], 
              K = 5,
              R = as.vector(respOpt_Cspan(OtherItems,NPL)),
              J = length(sigs)-1,
              N = length(unique(data[,"ID"])),
              Con = length(unique(data[,"Freetime"])),
              Freetime = unique(data[,"Freetime"]),
              retrievals = Retrievals,
              scale_b = 0.1)
  
  dat  
}


Analyze_M3 <- function(condition,dat,fixed_objects=NULL)
{
  
  Attach(condition)
  
  fit3 <- stan_fit(fixed_objects$M3_CS_EE, dat)
  
  ret <- fit3
  
  ret
  
 
  
}



Summarise <- function(condition, results, fixed_objects=NULL) {
  
  
  ret <- c(hyper_mu_c = results[1],hyper_mu_a = results[2], hyper_mu_f = results[3],hyper_mu_e = results[4],hyper_mu_r = results[5])   
  ret
  
}



SimClean()
res <- runSimulation(sim3, replications = 100, generate = Generate_M3, 
                     analyse = Analyze_M3, summarise = Summarise, 
                     fixed_objects = fo, parallel=TRUE, 
                     packages = c("cmdstanr","posterior","tmvtnorm","psych"),ncores =32)


SimExtract(res,what = "summarise")


stan1.dat<-Generate_M3(sim3[1,],fo)









