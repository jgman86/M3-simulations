################################################################################### Data Simulation ####################################################################################----
# Functions for Response Frequency Simulation

# Complex Span Model ----
simData_CSpan <- function(parmsMMM,respOpts,nRetrievals){
  # extract parms
  if(is.vector(parmsMMM)){
    conA <- parmsMMM["conA"]
    genA <- parmsMMM["genA"]
    filt <- parmsMMM["f"]
    baseA <- parmsMMM["baseA"]
  }else{
    conA <- parmsMMM[,"conA"]
    genA <- parmsMMM[,"genA"]
    filt <- parmsMMM[,"f"]
    baseA <- parmsMMM[,"baseA"]
  }
  
  # compute acts for response categories
  A_IIP <- conA + genA + baseA
  A_IOP <- genA + baseA
  A_DIP <- filt*(conA + genA) + baseA
  A_DOP <- filt*genA + baseA
  A_NPL <- baseA
  
  # summarize activations
  if(length(A_IIP) != 1){
    acts <- cbind(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }else{
    acts <- c(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }
  
  # compute summed activation
  sumActs <- apply(t(respOpts*t(acts)),1,sum)
  
  Probs <- t(respOpts*t(acts))/sumActs
  colnames(Probs) <- c("P_IIP","P_IOP","P_DIP","P_DOP","P_NPL")
  
  simdata <- matrix(NA,ncol = ncol(Probs),nrow = nrow(Probs))
  for(id in 1:nrow(Probs)){
    simdata[id,] <- t(stats::rmultinom(1,nRetrievals,Probs[id,]))
  }
  
  colnames(simdata) <- c("IIP","IOP","DIP","DOP","NPL")
  return(simdata)
}

# Simple Span Model ----
simData_SimpleSpan <- function(parmsMMM,respOpts,nTrials){
  # extract parms
  if(is.vector(parmsMMM)){
    conA <- parmsMMM["conA"]
    genA <- parmsMMM["genA"]
    baseA <- parmsMMM["baseA"]
  }else{
    conA <- parmsMMM[,"conA"]
    genA <- parmsMMM[,"genA"]
    baseA <- parmsMMM[,"baseA"]
  }
  
  # compute acts for response categories
  A_IIP <- conA + genA + baseA
  A_IOP <- genA + baseA
  A_NPL <- baseA
  
  # summarize activations
  if(length(A_IIP) != 1){
    acts <- cbind(A_IIP,A_IOP,A_NPL)
  }else{
    acts <- c(A_IIP,A_IOP,A_NPL)
  }
  
  # compute summed activation
  sumActs <- apply(t(respOpts*t(acts)),1,sum)
  
  Probs <- t(respOpts*t(acts))/sumActs
  colnames(Probs) <- c("P_IIP","P_IOP","P_NPL")
  
  simdata <- matrix(NA,ncol = ncol(Probs),nrow = nrow(Probs))
  for(id in 1:nrow(Probs)){
    simdata[id,] <- t(stats::rmultinom(1,nTrials,Probs[id,]))
  }
  
  colnames(simdata) <- c("IIP","IOP","NPL")
  return(simdata)
}



# Complex Span with extended Encoding ----
simData_CSpanEE <- function(parmsMMM,respOpts,nRetrievals,ft){
  # extract parms
  
  
  if(is.vector(parmsMMM)){
    conA <- parmsMMM["conA"]
    genA <- parmsMMM["genA"]
    filt <- parmsMMM["f"]
    e <- parmsMMM[,"e"]
    rm <- parmsMMM[,"r"]
    baseA <- parmsMMM["baseA"]
    
  }else{
    conA <- parmsMMM[,"conA"]
    genA <- parmsMMM[,"genA"]
    filt <- parmsMMM[,"f"]
    e <- parmsMMM[,"e"]
    rm <- parmsMMM[,"r"]
    baseA <- parmsMMM[,"baseA"]
  }
  
  
  
  # compute acts for response categories
  A_IIP <- baseA + (((1+e*ft)*conA) + genA) 
  A_IOP <- baseA + genA
  #A_IOP <- baseA + (1+e*ft)*genA
  A_DIP <- baseA + filt*((exp(-rm*ft)*conA) + genA) 
  #A_DOP <- baseA + exp(-rm*ft)*filt*genA 
  A_DOP <- baseA + filt*genA 
  A_NPL <- baseA 
  Freetime <- unique(ft)
  
  
  # summarize activations
  if(length(A_IIP) != 1){
    acts <- cbind(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }else{
    acts <- c(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }
  
  # compute summed activation
  sumActs <- apply(t(respOpts*t(acts)),1,sum)
  
  Probs <- t(respOpts*t(acts))/sumActs
  
  colnames(Probs) <- c("P_IIP","P_IOP","P_DIP","P_DOP","P_NPL")
  
  simdata <- matrix(NA,ncol = ncol(Probs),nrow = nrow(Probs))
  ID <- matrix(NA, ncol=1, nrow=nrow(Probs))
  i <- 1 
  
  for(id in 1:nrow(Probs)){
    
    simdata[id,] <- t(stats::rmultinom(1,nRetrievals,Probs[id,]))
    
    
    ID[id,] <- i
    
    if (id %% length(Freetime) == 0)
    {
      i <- i +1
      
    }
    
    
  }
  PC <- simdata[,1] / nRetrievals
  simdata <- cbind(ID,rep(Freetime,length.out = length(ID)),PC,simdata)
  colnames(simdata) <- c("ID","Freetime","PC","IIP","IOP","DIP","DOP","NPL")
  return(simdata)
}


simData_CSpanEE_bestfit <- function(parmsMMM,respOpts,nRetrievals,ft){
  # extract parms
  
  
  if(is.vector(parmsMMM)){
    conA <- parmsMMM["conA"]
    genA <- parmsMMM["genA"]
    filt <- parmsMMM["f"]
    e <- parmsMMM[,"e"]
    rm <- parmsMMM[,"r"]
    baseA <- parmsMMM["baseA"]
    
  }else{
    conA <- parmsMMM[,"conA"]
    genA <- parmsMMM[,"genA"]
    filt <- parmsMMM[,"f"]
    e <- parmsMMM[,"e"]
    rm <- parmsMMM[,"r"]
    baseA <- parmsMMM[,"baseA"]
  }
  
  
  
  # compute acts for response categories
  A_IIP <- baseA + (1+e*ft)*(conA + genA) 
  A_IOP <- baseA + (1+e*ft)*genA
  A_DIP <- baseA + exp(-rm*ft)*filt*(conA) + (filt*genA) 
  A_DOP <- baseA + filt*genA 
  A_NPL <- baseA 
  Freetime <- unique(ft)
  
  
  # summarize activations
  if(length(A_IIP) != 1){
    acts <- cbind(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }else{
    acts <- c(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }
  
  # compute summed activation
  sumActs <- apply(t(respOpts*t(acts)),1,sum)
  
  Probs <- t(respOpts*t(acts))/sumActs
  
  colnames(Probs) <- c("P_IIP","P_IOP","P_DIP","P_DOP","P_NPL")
  
  simdata <- matrix(NA,ncol = ncol(Probs),nrow = nrow(Probs))
  ID <- matrix(NA, ncol=1, nrow=nrow(Probs))
  i <- 1 
  
  for(id in 1:nrow(Probs)){
    
    simdata[id,] <- t(stats::rmultinom(1,nRetrievals,Probs[id,]))
    
    
    ID[id,] <- i
    
    if (id %% length(Freetime) == 0)
    {
      i <- i +1
      
    }
    
    
  }
  PC <- simdata[,1] / nRetrievals
  simdata <- cbind(ID,rep(Freetime,length.out = length(ID)),PC,simdata)
  colnames(simdata) <- c("ID","Freetime","PC","IIP","IOP","DIP","DOP","NPL")
  return(simdata)
}









# Complex Span Updating Variation ----
simData_UpdatingModel <- function(parmsMMM,respOpts,nRetrievals,CWI,WCI,fixtime,enctime){
  # extract parms
  
  
  if(is.vector(parmsMMM)){
    conA <- parmsMMM["conA"]
    genA <- parmsMMM["genA"]
    EU <- parmsMMM_30["eU"]
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
  
  
  # compute acts for response categories
  
  A_IIP <- baseA + (1+EU*t_EU)*(conA + genA) 
  A_IOP <- baseA + (1+EU*t_EU)*genA
  A_DIP <- baseA + exp(-rm*t_rm)*d*(1+EU*t_EU)*(conA + genA) 
  A_DOP <- baseA + exp(-rm*t_rm)*d*(1+EU*t_EU)*genA
  A_NPL <- baseA
  
  
  # summarize activations
  if(length(A_IIP) != 1){
    acts <- cbind(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }else{
    acts <- c(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }
  
  # compute summed activation
  sumActs <- apply(t(respOpts*t(acts)),1,sum)
  
  Probs <- t(respOpts*t(acts))/sumActs
  
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
  
  t_EU<-rep(t_EU,length.out = length(ID))
  t_rm<-rep(t_rm,length.out = length(ID))
  simdata <- cbind(ID,simdata,t_EU,t_rm)
  colnames(simdata) <- c("ID","IIP","IOP","DIP","DOP","NPL","t_EU","t_rm")
  return(simdata)
}

################################################################################### Response Categories ####################################################################################----

# Functions to manipulate the number of Distractors and the numbers of NPLs orthogonally ----

# N = Number of Distractors
# K = Number of NPLs

# IIP = Item in Position
# IOP = Item in other Position
# DIP = Distractor in Position
# DIOP = Distractor in other Position
# NPL = non presented Lure


# Simple Span Response Categories 
respOpt_Sspan <- function(N, K)
{
  k <-1
  respOpts <- matrix(NaN,nrow =length(N)*length(K), ncol = 3)
  
  for (i in N)
  {
    for(j in K)
    {
      
      IIP <- 1 # Item in Position
      IOP <- i # Items in Other Position
      NPL <- j  # Non-presented Lure
      
      
      respOpts[k,] <- cbind(IIP,IOP,NPL)
      
      k <- k+1
      
    }
    
  }
  colnames(respOpts) <- c("IIP","IOP","NPL")
  
  return(respOpts)
}


# For all ComplexSpan Variations (CS,CS_EE,EE) ----

respOpt_Cspan <- function(N, K)
{
  k <-1
  respOpts <- matrix(0,nrow =length(N)*length(K), ncol = 5)
  
  for (i in N)
  {
    for(j in K)
    {
      
      IIP <- 1 # Item in Position
      IOP <- i # Items in Other Position
      DIP <- 1 ## Distractor in Position
      DIOP <- i # Distractor in Other Position
      NPL <- j  # Non-presented Lure
      
      
      respOpts[k,] <- cbind(IIP,IOP,DIP,DIOP,NPL)
      
      k <- k+1
      
    }
    
  }
  
  colnames(respOpts) <- c("IIP","IOP","DIP","DIOP","NPL")
  
  return(respOpts)
}

# For the Simple Span Model ----
respOpt_Sspan <- function(N, K)
{
  k <-1
  respOpts <- matrix(0,nrow =length(N)*length(K), ncol = 3)
  
  for (i in N)
  {
    for(j in K)
    {
      
      IIP <- 1 # Item in Position
      IOP <- i # Items in Other Position
      NPL <- j  # Non-presented Lure
      
      
      respOpts[k,] <- cbind(IIP,IOP,NPL)
      
      k <- k+1
      
    }
    
  }
  
  colnames(respOpts) <- c("IIP","IOP","NPL")
  
  return(respOpts)
}




# Functions to manipulate only Set Size and the maximum response set ----
# ComplexSpan Variations ----

respOpt_CS_SSmaxSet <- function(SetSize, MaxSet)
{
  k <-1
  respOpts <- matrix(nrow =8, ncol = 5)
  for (value in SetSize)
  {
    
    IIP <- 1 # Item in Position
    IOP <- value-1 # Items in Other Position
    DIP <- 1 ## Distractor in Position
    DIOP <- value-1 # Distractor in Other Position
    
    NPL <- MaxSet - (sum(IIP,IOP,DIP,DIOP))  # Non-presented Lure
    
    
    respOpts[k,] <-cbind(IIP,IOP,DIP,DIOP,NPL)
    
    k <- k+1
    
  }
  
  colnames(respOpts) <- c("IIP","IOP","DIP","DIOP","NPL")
  
  return(respOpts)
}

# Simple Span ----
respOpt_SS_SSmaxSet <- function(SetSize, MaxSet)
{
  k <-1
  respOpts <- matrix(nrow =length(SetSize), ncol = 3)
  for (value in SetSize)
  {
    
    IIP <- 1 # Item in Position
    IOP <- value-1 # Items in Other Position
    NPL <- MaxSet - (sum(IIP,IOP))  # Non-presented Lure
    
    
    respOpts[k,] <-cbind(IIP,IOP,NPL)
    
    k <- k+1
    
  }
  
  colnames(respOpts) <- c("IIP","IOP","NPL")
  
  return(respOpts)
}
# Helper Functions ----
# Create Random Correlation Matrix from LKJ Distribution ----
rlkjcorr <- function ( n , K , eta = 1 ) {
  
  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  #if (K == 1) return(matrix(1, 1, 1))
  
  f <- function() {
    alpha <- eta + (K - 2)/2
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
    R[1,1] <- 1
    R[1,2] <- r12
    R[2,2] <- sqrt(1 - r12^2)
    if(K > 2) for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)
      
      # Draw uniformally on a hypersphere
      z <- rnorm(m, 0, 1)
      z <- z / sqrt(crossprod(z)[1])
      
      R[1:m,m+1] <- sqrt(y) * z
      R[m+1,m+1] <- sqrt(1 - y)
    }
    return(crossprod(R))
  }
  R <- replicate( n , f() )
  if ( dim(R)[3]==1 ) {
    R <- R[,,1]
  } else {
    # need to move 3rd dimension to front, so conforms to array structure that Stan uses
    R <- aperm(R,c(3,1,2))
  }
  return(R)
}


# Set Initial Values for STAN Sampler for Multivariate Sampling for ComplexSpan Models (Hardcoded soon)----
init_fun <- function() 
{
  
  list(subj_pars=cbind(runif(stan.dat$N,1,100),
                       runif(stan.dat$N,1,100),
                       runif(stan.dat$N,1,100),
                       runif(stan.dat$N,1,100)))
}
