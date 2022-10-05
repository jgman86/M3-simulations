
baseA <-c(0.1,0.1,0.1)
EU <- c(0.8,0.4,0.2)
conA <- c(6,9,11)
genA <- c(3,6,7)
rm <-  c(0.7,0.2,0.5)
d <- c(0.6,0.2,0.9)
t_EU <- c(0.2,1)
t_rm <- c(0.5,2)
Con1 <- 2
Con2 <- 2

time <- crossing(t_rm,t_EU)

acts <- matrix(NaN, ncol = 7,nrow =12)
colnames(acts) <- c("IIP","IOP","DIP","DIOP","NPL","t_EU","t_rm")


  c <- 1
  
          for(j in 1:3)
            for (i in 1:2)
              for(k in 1:2)
                
    {
  
      
      acts[c,1] <- j
      acts[c,2] <- 0.1 + (1+EU[j]*t_EU[i])*(conA[j] + genA[j])
      acts[c,2] <- 0.1 + (1+EU[j]*t_EU[i])*genA[j]
      acts[c,3] <- 0.1 + (exp(-rm[j]*t_rm[k])*d[j]*(1+EU[j]*t_EU[i])*conA[j]) + (1+EU[j]*t_EU[i])*genA[j]
      acts[c,4] <- 0.1 + (1+EU[j]*t_EU[i])*genA[j]
      acts[c,5] <- 0.1
      acts[c,6] <-t_EU[i]
      acts[c,7] <-t_rm[k]
      # 
      c <- c + 1 
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
