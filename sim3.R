#!/usr/bin/env Rscript

#####M3 Extended Encoding and Updating Simulations########

### Load Packages ####
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
library(unixtools)})


#### Read in Condtions for Current Job ####

# Create Option List for simulation conditions

option_list <- list(
  make_option(c("-N","--OtherItems"), type="integer",
              help="Number of OtherItems",action = "store",
              metavar="number"),
  make_option(c("-K","--NPL"), type="integer",
              help="Number of NPLs",action = "store",
              metavar="number"),
  make_option(c("-F","--nFreetime"), type="integer",
              help="Number of Freetime Conditions",action = "store",
              metavar="number"),
  make_option(c("-R","--nRetrievals"), type="integer", action="store",
              help="Number of Retrievals",
              metavar="number"),
  make_option(c("-P","--fixedf"), type="integer", action="store",
              help="fixed f condition",
              metavar="number"),
  make_option(c("-I","--JobID"),type="character", action="store",
              help="Slurm JobID",
              metavar="number"),
  make_option(c("-D","--JobDir"),type="character", action="store",
              help="Slurm JobDir",
              metavar="number")
  
)

con <- parse_args(OptionParser(option_list=option_list))

#### Set Up Cmdstan Interface and external Functions###
cmd_path <- paste0("~/R/x86_64-pc-linux-gnu-library/4.1/cmdstan-2.30.0")
quiet(set_cmdstan_path(cmd_path))

ID <- con$JobID

tmp<-paste0("/localscratch/",ID,"/ramdisk")
set.tempdir(tmp)

### Source Functions ####
source("Functions/M3_functions.R")

#### Define Simulation Design and SetUp Model----

###### Varying Simulation Factors ---- 
N <- con$OtherItems
K <- con$NPL
nRetrievals <- con$nRetrievals
nFT<- con$nFreetime

fixed_f <- con$fixedf
default.out <- con$JobDir

# Test Correct Argument Passthrough
# print(N)
# print(K)
# print(nRetrievals)
# print(nFT)
# print(ID)
# print(fixed_f)


# ###### Create Simulation Table ----
sim3 <- createDesign(OtherItems=N,
                     NPL=K,
                     Retrievals = nRetrievals,
                     nFreetime=nFT,
                     fixedf = fixed_f)

###### Fixed Simulation Factors ----
SampleSize <- 100
reps2con <- 100
minFT <- 0.250
maxFT <- 1.75

###### Simulation Options
n_iter = 1500
n_warmup= 1500
adapt_delta = .90
max_treedepth = 15


###### Model Path #####
full_model <- quiet(cmdstan_model("Models/M3_ComplexSpan_EE_LKJ_Cholesky_NC.stan"))
fixedf_model <- quiet(cmdstan_model("Models/M3_ComplexSpan_EE_LKJ_Cholesky_NC_fixedf.stan"))

#stan_path_M3_Upd <- paste0(dir_path,"/Models/M3_") add upadting model


###### Compile and store within fixed objects for use in runSimulation() #####
fo <- list(M3_CS_EE = full_model,
           M3_CS_EE_fixedf = fixedf_model,
           SampleSize = SampleSize,
           minFT = minFT,
           maxFT = maxFT,
           # Set Range for Parameter Means
           range_muC  =c(1,30),
           range_muA = c(0,0.5),
           range_muF = c(0,1), # fix to 0.5
           range_muE = c(0,0.8),
           range_muR = c(0,10), # range 0 - 25 empirical derived
           eta = 5, # Simulated N = 10000 with eta = 5, 95 % of all values lie within 0 -+ 0.56
           sigC = c(0.125,0.5),
           sigA = c(0.125,0.5),
           sigF = c(0.1,0.5),
           sigE = c(1,2),
           sigR = c(0.125,0.5), # abhänbgig von removal parameter -> analog zu c und a
           sigB = c(0.0001, 0.1),
           n_iter = n_iter,
           n_warmup= n_warmup,
           adapt_delta = adapt_delta,
           max_treedepth = max_treedepth
)


##### Set Up Fitting Function for cmdstan
stan_fit <- function(mod,dat,n_warmup,n_iter,adapt_delta,max_treedepth,fixedf){
  
  init <- function()
  {
    list(hyper_pars=cbind(runif(dat$J,10,20)),
         subj_pars=cbind(runif(dat$N,1,10)))
    
  }
  
  
  # Stan is noisy, so tell it to be more quiet()
  M3 <-  quiet(mod$sample(dat,
                          refresh = 0,
                          chains = 4,
                          parallel_chains=4,
                          iter_warmup=n_warmup,
                          iter_sampling=n_iter,
                          adapt_delta=adapt_delta,
                          max_treedepth=max_treedepth,
                          init = init,
                          show_messages = FALSE))
  
  M3_subj <- M3$summary(c("subj_pars"), mean,sd,rhat,Mode,HDInterval::hdi)
  M3_count_rep <- M3$summary(c("count_rep"),mean)
  M3_omega <- M3$summary("cor_mat_lower_tri",mean)
  
  if(fixedf == 0){
    M3_hyper <- M3$summary(c("hyper_pars","mu_f"), mean,Mode,sd,rhat,HDInterval::hdi)
    M3_f <- M3$summary(c("f"), mean,sd,Mode,rhat,HDInterval::hdi)
    M3_ac <- M3$summary(c("ac"), mean,sd,Mode,rhat,HDInterval::hdi)
    
    #Return Summary Object
    M3_sum <- list(M3_hyper,M3_subj,M3_f,M3_ac,M3_count_rep,M3_omega)
  } 
  else {
    M3_hyper <- M3$summary(c("hyper_pars"), mean,Mode,sd,rhat,HDInterval::hdi)
    
    #Return Summary Object
    M3_sum <- list(M3_hyper,M3_subj,M3_count_rep,M3_omega)
  }
  
  rm(M3)

  return(M3_sum)
}

##### Set Up Simulation Functions for Sim Design ----

Generate_M3 <- function(condition, fixed_objects=NULL) {
  Attach(condition)
  
  SampleSize <- fixed_objects$SampleSize
  minFT <- fixed_objects$minFT
  maxFT <- fixed_objects$maxFT
  range_muC  <- fixed_objects$range_muC
  range_muA <- fixed_objects$range_muA
  range_muF <- fixed_objects$range_muF
  range_muE <- fixed_objects$range_muE
  range_muR <- fixed_objects$range_muR
  eta <- fixed_objects$eta
  sigC <- fixed_objects$sigC
  sigA <- fixed_objects$sigA
  sigF <- fixed_objects$sigF
  sigE <- fixed_objects$sigE
  sigR <-fixed_objects$sigR
  sigB <- fixed_objects$sigB
  
  
  # Generate FreeTime Conditions
  
  conFT <- seq(from = minFT, to = maxFT, length.out=nFreetime) # eventuell log scale 0.2 0.8 2.4 oder so?
  
  
  # Sample Hyper-Parameter Means with C as fixpoint ----
  relCA <- runif(1, min = range_muA[1],max = range_muA[2])
  Mean_Cpar <- runif(1, min =range_muC[1], max= range_muC[2])
  Mean_Apar <- Mean_Cpar*(relCA)
  Mean_Epar <- runif(1, min =range_muE[1], max = range_muE[2])
  Mean_Rpar <- runif(1, min= range_muR[1], max = range_muR[2])
  
  
  if(fixedf == 0){
    Mean_Fpar <- runif(1, min= range_muF[1], max = range_muF[2])
    log_mu_f <- log(Mean_Fpar/(1-Mean_Fpar)) 
    n_theta = 5
  } else{
    Mean_Fpar <- 0.5 
    log_mu_f <- 0
    n_theta = 4
  }
  
  
  Mean_bpar <- 0.1
  
  # Make Vector with Hyper Pars
  
  hyper_mus <- c(Mean_Cpar,Mean_Apar,log_mu_f,Mean_Epar,Mean_Rpar, Mean_bpar)
  
  
  # Sample Variances and Set Covariances----
  
  sig_c <- runif(1, min = sigC[1], max = sigC[2])*Mean_Cpar
  sig_a <- runif(1, min = sigA[1], max = sigA[2])*Mean_Apar
  
  
  if(fixedf == 0) {
    sig_f <- runif(1, min = sigF[1], max= sigF[2])
  } else{
    sig_f <- 0.001
  }
  
  
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
  
  
  if(fixedf == 0){
    omega[6,1:5] = omega[1:5,6] = 0 #fix cov of b to 0    
  } else {
    omega[6,1:5] = omega[1:5,6] = 0 #fix cov of b to 0    
    omega[3,1:5] = omega[1:5,3] = 0.0001 #fix cov of f to 0 
  }
  
  
  
  Sigma <- cor2cov(omega,sigs)
  
  
  # Sample Parameters from MVN ----
  
  parms <- tmvtnorm::rtmvnorm(n=SampleSize, mean= hyper_mus, sigma=Sigma,
                              lower=c(0,0,-Inf,0,0,0),upper = c(Inf,Inf,Inf,Inf,Inf,Inf))
  # Merge Parameters to one Matrix
  colnames(parms) <- c("conA","genA","f","e","r","baseA")
  parms[,6] <- 0.1
  
  if(fixedf == 0){
    
    parms[,3] <- 1 / (1+exp(-parms[,3]))
    
    # Create Compound of r and f to counter parameter trade-offs
    
    ac <- exp(parms[,3]-1) * parms[,5]
    parms <- cbind(parms, ac)
    
    
  } else {
    parms[,3] <- 0.5 # 0 if f is fixed to 0.5
    
  }
  
  
  # Simulate Data from theta ----
  
  ParmsFT <- matrix(rep(parms,each =nFreetime), nrow = length(parms[,1])*nFreetime, ncol = ncol(parms), byrow = F)
  
  if(fixedf == 0){
    colnames(ParmsFT) <- c("conA","genA","f","e","r","baseA","ac")
  } else {
    colnames(ParmsFT) <- c("conA","genA","f","e","r","baseA")
  }
  
  FT <- rep(conFT,length.out = nrow(ParmsFT))
  
  data <- simData_CSpanEE(ParmsFT,as.vector(respOpt_Cspan(OtherItems,NPL)),Retrievals,FT)
  
  
  # Generate Stan Data ----
  
  
  stan.dat <- list(count = data[,4:8],
                   K = 5,
                   R = as.vector(respOpt_Cspan(OtherItems,NPL)),
                   J = n_theta,
                   N = length(unique(data[,"ID"])),
                   Con = length(unique(data[,"Freetime"])),
                   Freetime = unique(data[,"Freetime"]),
                   retrievals = Retrievals,
                   scale_b = 0.1)
  
  theta <- parms
  dat  <- list(stan.dat,theta,data,hyper_mus)
  return(dat)
}


Analyze_M3 <- function(condition,dat,fixed_objects=NULL){
  
  Attach(condition)
  
  if (fixedf == 0){
    
    model <-fixed_objects$M3_CS_EE
    
  } else{
    
    model <- fixed_objects$M3_CS_EE_fixedf
    
  }
  
  
  fit3 <- stan_fit(model, dat[[1]],
                   n_warmup = fixed_objects$n_warmup,
                   n_iter=fixed_objects$n_iter,
                   adapt_delta=fixed_objects$adapt_delta,
                   max_treedepth=fixed_objects$max_treedepth,fixedf)
  
  
  ## Calculate Current Repetitions Row
  
  # Merge Data depending on model conditions
  
  if(fixedf == 0){
    
    hyper <- fit3[[1]] %>% mutate(variable = case_when(variable == "hyper_pars[1]" ~ "mu_c",
                                                       variable == "hyper_pars[2]" ~ "mu_a",
                                                       variable == "hyper_pars[3]" ~ "logMu_f",
                                                       variable == "hyper_pars[4]" ~ "mu_e",
                                                       variable == "hyper_pars[5]" ~ "mu_r",
                                                       variable == "mu_f" ~ "mu_f")) %>%
      pivot_wider(.,names_from = "variable",values_from = c("mean","Mode","sd","rhat","upper","lower"))
    
    subj <- fit3[[2]]  %>% mutate(variable = str_remove_all(variable, "subj_pars")) %>%
      separate(col = variable,into = c("theta","ID"),sep = ",")  %>%
      mutate(ID = str_remove(ID,pattern = "]"), theta = case_when(theta == "[1" ~ "c", theta == "[2" ~ "a",
                                                                  theta == "[3" ~ "log_f",
                                                                  theta == "[4" ~ "EE",
                                                                  theta== "[5" ~ "r")) %>%
      pivot_wider(.,names_from = "theta",values_from = c("mean","Mode","sd","rhat","upper","lower"))
    
    
    f <- fit3[[3]] %>% mutate(ID = seq(1:fixed_objects$SampleSize), theta = "f") %>% relocate(c("ID","theta"), .before = variable) %>% select(-variable) %>%
      pivot_wider(.,names_from = "theta",values_from = c("mean","Mode","sd","rhat","upper","lower"))
    
    ac <- fit3[[4]] %>% mutate(ID = seq(1:fixed_objects$SampleSize), theta = "ac") %>% relocate(c("ID","theta"), .before = variable) %>% select(-variable) %>%
      pivot_wider(.,names_from = "theta",values_from = c("mean","Mode","sd","rhat","upper","lower"))
    
    count_rep <- fit3[[5]] %>% separate(variable,into = c("ID","Category"),sep = ",") %>%
      mutate(ID= str_remove_all(ID,pattern="count_rep\\["), Category=str_remove_all(Category,"\\]"), mean=round(mean,0)) %>%
      pivot_wider(., names_from = Category, values_from = mean) %>%
      rename("IIP"=`1`,"IOP" = `2`,"DIP" = `3`,"DIOP" = `4`,"NPL" = `5`)
    
    # RMSE for Hyper Pars
    
    rmse_mu_c <- RMSE(hyper$mean_mu_c,dat[[4]][1])
    rmse_mu_a <- RMSE(hyper$mean_mu_a,dat[[4]][2])
    rmse_mu_f <- RMSE(hyper$mean_mu_f,1/(1+exp(-dat[[4]][3])))
    rmse_mu_e <- RMSE(hyper$mean_mu_e,dat[[4]][4])
    rmse_mu_r <- RMSE(hyper$mean_mu_r,dat[[4]][5])
    
    # Recoveries of Subject theta
    recoveries_c <- cor(dat[[2]][,"conA"],subj$mean_c)
    recoveries_a <- cor(dat[[2]][,"genA"],subj$mean_a)
    recoveries_f <- cor(dat[[2]][,"f"],f$mean_f)
    recoveries_e <- cor(dat[[2]][,"e"],subj$mean_EE)
    recoveries_r <- cor(dat[[2]][,"r"],subj$mean_r)
    recoveries_ac <- cor(dat[[2]][,"ac"],ac$mean_ac)
    
    # RMSE for subject theta
    
    rmse_c <- RMSE(subj$mean_c,dat[[2]][,"conA"])
    rmse_a <- RMSE(subj$mean_a,dat[[2]][,"genA"])
    rmse_f <- RMSE(f$mean_f,dat[[2]][,"f"])
    rmse_e <- RMSE(subj$mean_EE,dat[[2]][,"e"])
    rmse_r <- RMSE(subj$mean_r,dat[[2]][,"r"])
    rmse_ac <- RMSE(ac$mean_ac,dat[[2]][,"ac"])
    
    # max rhat for subject theta
    
    max_rhat_c <- subj %>% select(contains("rhat_c")) %>% max()
    max_rhat_a <- subj %>% select(contains("rhat_a")) %>% max()
    max_rhat_f <- f %>% select(contains("rhat_f")) %>% max()
    max_rhat_e <- subj %>% select(contains("rhat_EE")) %>% max()
    max_rhat_r <- subj %>% select(contains("rhat_r")) %>% max()
    
    # Predictive Correlation for Different Response Categories
    
    cor_rep_IIP <- cor(dat[[3]][,"IIP"],count_rep$IIP)
    cor_rep_IOP <- cor(dat[[3]][,"IOP"],count_rep$IOP)
    cor_rep_DIP <- cor(dat[[3]][,"DIP"],count_rep$DIP)
    cor_rep_DIOP <- cor(dat[[3]][,"DOP"],count_rep$DIOP)
    cor_rep_NPL <- cor(dat[[3]][,"NPL"],count_rep$NPL)
    
    # Merge Results for current iteration
    
    ret <- data.frame(re_c = recoveries_c,
                      rmse_c = rmse_c,
                      rmse_mu_c=rmse_mu_c,
                      rhat_c = max_rhat_c,
                      
                      re_a=recoveries_a,
                      rmse_a = rmse_a,
                      rmse_mu_a=rmse_mu_a,
                      rhat_a=max_rhat_a,
                      
                      re_f= recoveries_f,
                      rmse_f = rmse_f,
                      rmse_mu_f=rmse_mu_f,
                      rhat_f=max_rhat_f,
                      
                      re_e = recoveries_e,
                      rmse_e = rmse_e,
                      rmse_mu_e=rmse_mu_e,
                      rhat_e = max_rhat_e,
                      
                      re_r=recoveries_r,
                      rmse_r = rmse_r,
                      rmse_mu_r=rmse_mu_r,
                      rhat_r=max_rhat_r,
                      
                      re_ac=recoveries_ac,
                      rmse_ac = rmse_ac,
                      
                      cor_IIP = cor_rep_IIP,
                      cor_IOP = cor_rep_IOP,
                      cor_DIP = cor_rep_DIP,
                      cor_DIOP= cor_rep_DIOP,
                      cor_NPL = cor_rep_NPL)
    
    # Return Results from current Iteration 
    
    ret
    
  }
  
  else {
    
    hyper <- fit3[[1]] %>% mutate(variable = case_when(variable == "hyper_pars[1]" ~ "mu_c",
                                                       variable == "hyper_pars[2]" ~ "mu_a",
                                                       variable == "hyper_pars[3]" ~ "mu_e",
                                                       variable == "hyper_pars[4]" ~ "mu_r")) %>%
                                    pivot_wider(.,names_from = "variable",values_from = c("mean","Mode","sd","rhat","upper","lower"))
    
    subj <- fit3[[2]]  %>% mutate(variable = str_remove_all(variable, "subj_pars")) %>%
      separate(col = variable,into = c("theta","ID"),sep = ",")  %>%
      mutate(ID = str_remove(ID,pattern = "]"), theta = case_when(theta == "[1" ~ "c", 
                                                                  theta == "[2" ~ "a",
                                                                  theta == "[3" ~ "EE",
                                                                  theta== "[4" ~ "r")) %>%
      pivot_wider(.,names_from = "theta",values_from = c("mean","Mode","sd","rhat","upper","lower"))
    
    count_rep <- fit3[[3]] %>% separate(variable,into = c("ID","Category"),sep = ",") %>%
      mutate(ID= str_remove_all(ID,pattern="count_rep\\["), Category=str_remove_all(Category,"\\]"), mean=round(mean,0)) %>%
      pivot_wider(., names_from = Category, values_from = mean) %>%
      rename("IIP"=`1`,"IOP" = `2`,"DIP" = `3`,"DIOP" = `4`,"NPL" = `5`)
    
    
    # RMSE for Hyper Pars
    
    rmse_mu_c <- RMSE(hyper$mean_mu_c,dat[[4]][1])
    rmse_mu_a <- RMSE(hyper$mean_mu_a,dat[[4]][2])
    rmse_mu_e <- RMSE(hyper$mean_mu_e,dat[[4]][3])
    rmse_mu_r <- RMSE(hyper$mean_mu_r,dat[[4]][4])
    
    # Recoveries of Subject theta
    
    recoveries_c <- cor(dat[[2]][,"conA"],subj$mean_c)
    recoveries_a <- cor(dat[[2]][,"genA"],subj$mean_a)
    recoveries_e <- cor(dat[[2]][,"e"],subj$mean_EE)
    recoveries_r <- cor(dat[[2]][,"r"],subj$mean_r)

    # RMSE for subject theta
    
    rmse_c <- RMSE(subj$mean_c,dat[[2]][,"conA"])
    rmse_a <- RMSE(subj$mean_a,dat[[2]][,"genA"])
    rmse_e <- RMSE(subj$mean_EE,dat[[2]][,"e"])
    rmse_r <- RMSE(subj$mean_r,dat[[2]][,"r"])

    # max rhat for subject theta
    
    max_rhat_c <- subj %>% select(contains("rhat_c")) %>% max()
    max_rhat_a <- subj %>% select(contains("rhat_a")) %>% max()
    max_rhat_e <- subj %>% select(contains("rhat_EE")) %>% max()
    max_rhat_r <- subj %>% select(contains("rhat_r")) %>% max()
    
    # Predictive Correlation for Different Response Categories
    
    cor_rep_IIP <- cor(dat[[3]][,"IIP"],count_rep$IIP)
    cor_rep_IOP <- cor(dat[[3]][,"IOP"],count_rep$IOP)
    cor_rep_DIP <- cor(dat[[3]][,"DIP"],count_rep$DIP)
    cor_rep_DIOP <- cor(dat[[3]][,"DOP"],count_rep$DIOP)
    cor_rep_NPL <- cor(dat[[3]][,"NPL"],count_rep$NPL)
    
    # Merge Results for current iteration
    
    ret <- data.frame(re_c = recoveries_c,
                      rmse_c = rmse_c,
                      rmse_mu_c=rmse_mu_c,
                      rhat_c = max_rhat_c,
                      re_a=recoveries_a,
                      rmse_a = rmse_a,
                      rmse_mu_a=rmse_mu_a,
                      rhat_a=max_rhat_a,
                      re_e = recoveries_e,
                      rmse_e = rmse_e,
                      rmse_mu_e=rmse_mu_e,
                      rhat_e = max_rhat_e,
                      re_r=recoveries_r,
                      rmse_r = rmse_r,
                      rmse_mu_r=rmse_mu_r,
                      rhat_r=max_rhat_r,
                      cor_IIP = cor_rep_IIP,
                      cor_IOP = cor_rep_IOP,
                      cor_DIP = cor_rep_DIP,
                      cor_DIOP= cor_rep_DIOP,
                      cor_NPL = cor_rep_NPL)
    
    # Return Results from current Iteration 
    
    ret
    
  }
  
  
}


Summarise <- function(condition, results, fixed_objects=NULL) {
  Attach(condition)
  
  if (fixedf == 0){
    
    # Calculate Meta Statistics over all Replications
    ret <- c(cor_c = fisherz2r(mean(fisherz(results$re_c))),
             max_rhat_c = max(results$rhat_c),
             mean_rmse_mu_c = mean(results$rmse_mu_c),
             mean_rmse_c = mean(results$rmse_c),
             cor_a=fisherz2r(mean(fisherz(results$re_a))),
             mean_rmse_mu_a = mean(results$rmse_mu_a),
             mean_rmse_a = mean(results$rmse_a),
             max_rhat_a = max(results$rhat_a),
             cor_f = fisherz2r(mean(fisherz(results$re_f))),
             mean_rmse_mu_f = mean(results$rmse_mu_f),
             mean_rmse_f = mean(results$rmse_f),
             max_rhat_f = max(results$rhat_f),
             cor_e = fisherz2r(mean(fisherz(results$re_e))),
             mean_rmse_mu_e = mean(results$rmse_mu_e),
             mean_rmse_e = mean(results$rmse_e),
             max_rhat_e = max(results$rhat_e),
             cor_r = fisherz2r(mean(fisherz(results$re_r))),
             mean_rmse_mu_r = mean(results$rmse_mu_r),
             mean_rmse_r = mean(results$rmse_r),
             max_rhat_r = max(results$rhat_r),
             cor_ac = fisherz2r(mean(fisherz(results$re_ac))),
             mean_rmse_ac = mean(results$rmse_ac),
             cor_IIP = fisherz2r(mean(fisherz(results$cor_IIP))),
             cor_IOP = fisherz2r(mean(fisherz(results$cor_IOP))),
             cor_DIP = fisherz2r(mean(fisherz(results$cor_DIP))),
             cor_DIOP = fisherz2r(mean(fisherz(results$cor_DIOP))),
             cor_NPL = fisherz2r(mean(fisherz(results$cor_NPL))))
    
     ret
    
  } else {
    
    
    ret <- c(cor_c = fisherz2r(mean(fisherz(results$re_c))),
             max_rhat_c = max(results$rhat_c),
             mean_rmse_mu_c = mean(results$rmse_mu_c),
             mean_rmse_c = mean(results$rmse_c),
             cor_a=fisherz2r(mean(fisherz(results$re_a))),
             mean_rmse_mu_a = mean(results$rmse_mu_a),
             mean_rmse_a = mean(results$rmse_a),
             max_rhat_a = max(results$rhat_a),
             cor_e = fisherz2r(mean(fisherz(results$re_e))),
             mean_rmse_mu_e = mean(results$rmse_mu_e),
             mean_rmse_e = mean(results$rmse_e),
             max_rhat_e = max(results$rhat_e),
             cor_r = fisherz2r(mean(fisherz(results$re_r))),
             mean_rmse_mu_r = mean(results$rmse_mu_r),
             mean_rmse_r = mean(results$rmse_r),
             max_rhat_r = max(results$rhat_r),
             cor_IIP = fisherz2r(mean(fisherz(results$cor_IIP))),
             cor_IOP = fisherz2r(mean(fisherz(results$cor_IOP))),
             cor_DIP = fisherz2r(mean(fisherz(results$cor_DIP))),
             cor_DIOP = fisherz2r(mean(fisherz(results$cor_DIOP))),
             cor_NPL = fisherz2r(mean(fisherz(results$cor_NPL))))
    
    ret
    
  } 
  
}


runSimulation(sim3, replications = reps2con, generate = Generate_M3,
              analyse = Analyze_M3, summarise = Summarise,parallel = T,
              fixed_objects = fo,filename=paste0("M3_EE_CS_FULL_",ID),
              save_details = list(out_rootdir=paste0(default.out,"Results/")),
              packages = c("cmdstanr","posterior","tmvtnorm","psych","tidyverse","tidybayes"))


