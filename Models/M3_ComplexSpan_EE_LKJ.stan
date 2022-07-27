// Implemented in Stan by Jan Goetttmann 
// M3 - Model for extended encoding in complex span tasks after Oberauer & Lewandowsky, 2018
// Questions regarding the model code to Jan.Goettmann@psychologie.uni-heidelberg.de

// Best Fit for Complex Span Tasks

data {
  int <lower=0> N;  // number rownumber
  int <lower=0> K;  // categories 
  int <lower=0> J;  // Dims of Cov Matrix
  int <lower=0> Con;  // number of categories
  int R[K];         // number of responses per category
  int count[N*Con,K];   // observed data
  real scale_b;     // set scaling for background noise
  vector[Con] Freetime ;        // freetime conditions after distractor encoding (e.g. 500 ms, 1000 ms or 1500 ms) 
  real f;
}

parameters {
  
  // Defining vector for hyper and subject parameters 
  
  corr_matrix[J] Omega;
  vector<lower=0>[J] sigma;
  vector [J] hyper_pars;
  vector [J] subj_pars[N];
  
  
}


transformed parameters{
  // Transform f Parameter
  
  //real f[N] = inv_logit(subj_pars[,3]);
  //real mu_f = inv_logit(hyper_pars[3]);
  
  // activations
  real acts_IIP[N*Con];
  real acts_IOP[N*Con];
  real acts_DIP[N*Con];
  real acts_DIOP[N*Con];
  real acts_NPL[N*Con];
  
  
  // probabilities
  vector[K] probs[N*Con];
  real SummedActs[N*Con];
  
  
  // loop over subjects and conditions to compute activations and probabilites
  
  for (i in 1:N){ // for each subject
  
  for(j in 1:Con) {
    
    acts_IIP[j + (i-1)*Con] = scale_b + ((1+ subj_pars[i,3]*Freetime[j])*subj_pars[i,1]) + subj_pars[i,2]; // Item in Position                      
    acts_IOP[j + (i-1)*Con] = scale_b + subj_pars[i,2];        // Item in Other Position
    acts_DIP[j + (i-1)*Con] = scale_b + f*((exp(-subj_pars[i,4]*Freetime[j])* subj_pars[i,2]) + subj_pars[i,1]);// Distractor in Position
    acts_DIOP[j + (i-1)*Con] = scale_b +  f*subj_pars[i,2]; // Distractor in other Position
    acts_NPL[j + (i-1)*Con] = scale_b; // non presented Lure
    
    SummedActs[j + (i-1)*Con] = R[1] * acts_IIP[j + (i-1)*Con] + R[2] * acts_IOP[j + (i-1)*Con] + R[3] * acts_DIP[j + (i-1)*Con]+
    R[4] * acts_DIOP[j + (i-1)*Con]+ R[5] * acts_NPL[j + (i-1)*Con];
    
    probs[j + (i-1)*Con,1] = (R[1] * acts_IIP[j + (i-1)*Con]) ./ (SummedActs[j + (i-1)*Con]);  
    probs[j + (i-1)*Con,2] = (R[2] * acts_IOP[j + (i-1)*Con]) ./ (SummedActs[j + (i-1)*Con]);
    probs[j + (i-1)*Con,3] = (R[3] * acts_DIP[j + (i-1)*Con]) ./ (SummedActs[j + (i-1)*Con]);
    probs[j + (i-1)*Con,4] = (R[4] * acts_DIOP[j + (i-1)*Con]) ./ (SummedActs[j + (i-1)*Con]);
    probs[j + (i-1)*Con,5] = (R[5] * acts_NPL[j + (i-1)*Con]) ./ (SummedActs[j + (i-1)*Con]);
  }
  }
}


model {
  
  // priors for hyper parameters
  
  hyper_pars[1] ~ normal(20,10); // c
  hyper_pars[2] ~ normal(2,10); // a
  //hyper_pars[3] ~ normal(0,10); //f
  hyper_pars[3] ~ normal(1,10); // EE
  hyper_pars[4] ~ normal(1,10); // r
  
  Omega ~ lkj_corr(5);
  sigma ~ gamma(1,0.01);
  
  
  // Loop over subjects
  subj_pars[,] ~ multi_normal(hyper_pars,quad_form_diag(Omega, sigma));
  
  
  for (j in 1:Con) {
    for (i in 1:N) {
      // draw data from probabilities determined by MMM parms
      count[j + (i-1)*Con,]  ~ multinomial(probs[j + (i-1)*Con,]);  
    }
  }
}


