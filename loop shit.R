
baseA <-rep(0.1,5) 
EU <- rep(0.8,5)
conA <- rep(8,5)
genA <- rep(4,5)
rm <- rep(0.76,5)
d <- rep(0.6,5)
t_EU <- c(0.2,1)
t_rm <- c(0.5,2)
Con1 <- 2
Con2 <- 2
c <- 1
time <- crossing(t_rm,t_EU)

acts <- matrix(NaN, ncol = 7,nrow = 3)
acts_IIP <-c()
colnames(acts) <- c("IIP","IOP","DIP","DIOP","NPL","t_EU","t_rm")
  
for (i in 1:100)
{
  for (j in 1:2)
    {
  

      acts_IIP[j + (i-1)*Con1] <- baseA + (1+EU*t_EU[j])*(conA + genA) 
      # acts[j + (i-1)*Con1*Con2,2] <- baseA + (1+EU*t_EU[j])*genA
      # acts[j + (i-1)*Con1*Con2,3] <- baseA + (exp(-rm*t_rm[j])*d*(1+EU*t_EU[j])*conA) + (1+EU*t_EU[j])*genA 
      # acts[j + (i-1)*Con1*Con2,4] <- baseA + (1+EU*i)*genA
      # acts[j + (i-1)*Con1*Con2,5] <- baseA
      # acts[j + (i-1)*Con1*Con2,6] <- i
      # acts[j + (i-1)*Con1*Con2,7] <- j
    
  }
}

rm(acts)
for (i in 1:100)
{
  for (j in 1:4)
  {
    print(j + (i-1)*4) 
  }
}






A_IIP <- baseA + (1+EU*t_EU)*(conA + genA) 
A_IOP <- baseA + (1+EU*t_EU)*genA
A_DIP <- baseA + (exp(-rm*t_rm[j])*d*(1+EU*t_EU[i])*conA) + (1+EU*t_EU[i])*genA 
A_DOP <- baseA + (1+EU*t_EU)*genA
A_NPL <- baseA


acts <- matrix(NaN, ncol = 7,nrow = nrow(ParmsUpd)*9)


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
  
  # compute acts for response categories
  acts <- matrix(NaN, ncol = 7,nrow = nrow(parmsMMM))
  #colnames(acts) <- c("IIP","IOP","DIP","DIOP","NPL","t_EU","t_rm")
  
  for (h in 1:(nrow(parmsMMM)*Cons))
  {
    for (i in t_EU)
      for (j in t_rm)
      {
        {
          acts[h,1] <- baseA + (1+EU*i)*(conA + genA) 
          acts[h,2] <- baseA + (1+EU*i)*genA
          acts[h,3] <- baseA + (exp(-rm*j)*d*(1+EU*i)*conA) + (1+EU*i)*genA 
          acts[h,4] <- baseA + (1+EU*i)*genA
          acts[h,5] <- baseA
          acts[h,6] <- i
          acts[h,7] <- j
          
        }
      }
  }
  return(acts)
}
data <- simData_UpdatingModel_test(ParmsUpd,as.vector(respOpt_Cspan(N,K)),nRetrievals,conCWI,conWCI,fixtime,enctime)
