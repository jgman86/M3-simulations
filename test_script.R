file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod3_nc <- cmdstan_model(, cpp_options=list(stan_threads=T))
mod3_norm <- cmdstan_model(stan_path_M3_EE)
mod3_univ <- cmdstan_mode()

# Simulate Data for Estimation ----
# 
# 
# 
# 

nRetrievals <- 500
minFT <- 0.2
maxFT <- 2
nFT <- c(2,4) # 2,4,10 Conditions between 0.2 and 2
SampleSize <- 100
con_nFT =4
conFT <- seq(from = minFT, to = maxFT, length.out = con_nFT) # eventuell log scale 0.2 0.8 2.4 oder so?

# Set Range for Parameter Means
range_muC <- c(1,100)
range_muA <- c(0,0.5)
range_muF <- c(0,1) # fix to 0.5
range_muE <-c(0,0.5)
range_muR <- c(0,25) # range 0 - 25 empirical derived
eta <- 5 # Simulated N = 10000 with eta = 5, 95 % of all values lie within 0 -+ 0.56


sigC <- c(0.125,0.5)
sigA <- c(0.125,0.5)
sigF <- c(0.0001,0.1)
sigE <- c(1,2)
sigR <- c(0.125,0.5) # abhänbgig von removal parameter -> analog zu c und a
sigB <- c(0.0001, 0.1)



# Sample Hyper-Parameter Means with C as fixpoint ----
relCA <- runif(1, min = range_muA[1],max = range_muA[2])
Mean_Cpar <- 8
Mean_Apar <- 4
Mean_Epar <- runif(1, min =range_muE[1], max = range_muE[2])
Mean_Rpar <- runif(1, min= range_muR[1], max = range_muR[2])
Mean_Fpar <-  runif(1, min= range_muF[1], max = range_muF[2])
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


# Set Correlations for Parameters ----

# Sample Covariance Matrix Sigma

omega <- rlkjcorr(1,length(hyper_mus),eta)

# Little Hack for fixing coveriance of b to zer0

omega[6,1:5] = omega[1:5,6] = 0 #fix cov of b to 0
#mega[3,1:5] = omega[1:5,3] = .00001 #fix cov of f to 0


Sigma <- cor2cov(omega,sigs)


# Sample Parameters from MVN ----

parms <- tmvtnorm::rtmvnorm(n=SampleSize, mean= hyper_mus, sigma=Sigma,
                            lower=c(0,0,-Inf,0,0,0),upper = c(Inf,Inf,Inf,Inf,Inf,Inf))
# Merge Parameters to one Matrix
colnames(parms) <- c("conA","genA","f","e","r","baseA")
parms[,6] <- 0.1
#parms[,3] <- 0
parms[,3] <- 1 / (1+exp(-parms[,3]))

#parms[,3] <- 0 if f is fixed to 0.5


# Simulate Data for Estimation ----

ParmsFT <- matrix(rep(parms,each =con_nFT), nrow = length(parms[,1])*con_nFT, ncol = ncol(parms), byrow = F)
colnames(ParmsFT) <- c("conA","genA","f","e","r","baseA")
FT <- rep(conFT,length.out = nrow(ParmsFT))

data <- simData_CSpanEE(ParmsFT,as.vector(respOpt_Cspan(4,16)),1000,FT)


# Generate Stan Data ----




init_fun_FT_nc_univ <- function()
{
  list(c=cbind(runif(stan.dat$N,1,10)),
       a=cbind(runif(stan.dat$N,1,10)),
       log_f=cbind(runif(stan.dat$N,-1,1)),
       EE=cbind(runif(stan.dat$N,1,10)),
       r=cbind(runif(stan.dat$N,1,10)),
        mu_c=cbind(runif(stan.dat$N,1,10)),
        mu_a=cbind(runif(stan.dat$N,1,10)),
        mu_r=cbind(runif(stan.dat$N,1,10)),
        mu_e=cbind(runif(stan.dat$N,1,10)))
   
}

init_fun_FT_nc_univ <- function()
{
  list(hyper_pars=cbind(runif(stan.dat$J,1,10)),
        subj_pars=cbind(runif(stan.dat$N,1,10)))
  
}



stan.dat <- list(count = data[,4:8], 
                 K = 5,
                 R = as.vector(respOpt_Cspan(4,16)),
                 J = length(sigs)-1,
                 N = length(unique(data[,"ID"])),
                 Con = length(unique(data[,"Freetime"])),
                 Freetime = unique(data[,"Freetime"]),
                 retrievals = 1000,
                 scale_b = 0.1)
                  #f = .5)

fit3 <- mod3_norm$sample(data = stan1.dat,
                         chains = 4, 
            parallel_chains = 4,iter_warmup = 500,iter_sampling = 1000,
            adapt_delta = .99, max_treedepth = 15,
            init = init_fun_FT_nc_univ,
            show_messages = F)
            
  
  gf<-stan_fit(mod3_norm,dat = stan1.dat)
  rbind(gf)
  bayesplot::mcmc_combo(fit3$draws(c("hyper_pars")),alpha=0.2)
  bayesplot::rhat(fit3$cmdstan_summary())
  
          ret <- as.vector((hyper_mu_c = gf[1,2])),hyper_mu_a = results[2,2], hyper_mu_f = results[3,2],hyper_mu_e = results[4,2],hyper_mu_r = results[4,2])              
  
  post_samples_hyper <- as_draws_matrix(fit3$draws(c("hyper_pars","mu_f")))
  post_samples_subj <- as_draws_matrix(fit3$draws(c("c","a","f","EE","r")))

 
  
 summary_hyper <- posterior::summarise_draws(post_samples_hyper)
 summary_subj <- posterior::summarise_draws(post_samples_subj)
  