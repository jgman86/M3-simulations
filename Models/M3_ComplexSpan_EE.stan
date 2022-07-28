// Implemented in Stan by Jan Goetttmann 
// M3 - Model for extended encoding in complex span tasks after Oberauer & Lewandowsky, 2018
// Questions regarding the model code to Jan.Goettmann@psychologie.uni-heidelberg.de



data {
  int <lower=0> N;  // number rownumber
  int <lower=0> K;  // categories 
  int <lower=0> Con;  // number of categories
  int R[K];         // number of responses per category
  int count[N*Con,K];   // observed data
  real scale_b;     // set scaling for background noise
  vector[Con] Freetime ;        // freetime conditions after distractor encoding (e.g. 500 ms, 1000 ms or 1500 ms) 
}

parameters {
  // subject parameters
  real <lower=0>a[N];
  real <lower=0>c[N];
  real <lower=0>r[N];
  real <lower=0>EE[N];
  // real log_f[N];
  
  // Mu & Sigma for hyper distributions
  real mu_a;
  real <lower=0> sig_a;
  real mu_c;
  real <lower=0> sig_c;
  real mu_r;
  real <lower=0> sig_r;
  real mu_e;
  real <lower=0> sig_e;
  // real  logMu_f;
  // real <lower=0> logSig_f;
  
}

transformed parameters{
  // Transform f Parameter
  
  // real f[N] = inv_logit(log_f);
  // real mu_f = inv_logit(logMu_f);
  // real sig_f = sd(f);
  // 
  
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
    
    acts_IIP[j + (i-1)*Con] = scale_b + ((1+EE[i]*Freetime[j])*c[i]) + a[i]; // Item in Position                      
    acts_IOP[j + (i-1)*Con] = scale_b + a[i];        // Item in Other Position
    acts_DIP[j + (i-1)*Con] = scale_b + 0.5*((exp(-r[i]*Freetime[j])*c[i])+a[i]);// Distractor in Position
    acts_DIOP[j + (i-1)*Con] = scale_b + 0.5 *a[i]; // Distractor in other Position
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
  
  mu_c ~ normal(20,10);
  sig_c ~ gamma(1,0.01);
  
  mu_a ~ normal(2,10);
  sig_a ~ gamma(1,0.01);
  
  mu_r ~ normal(1,10);
  sig_r~ gamma(1,0.01);
  
  mu_e ~ normal(1,10);
  sig_e ~ gamma(1,0.01);
  
  // logMu_f ~ normal(0,1);
  // logSig_f ~ gamma(1,0.01);
  
  // Loop over subjects
  for(i in 1:N){
    // Draw subject parameters from truncated normal
    
    c[i] ~ normal(mu_c, sig_c);
    a[i] ~ normal(mu_a, sig_a);
    r[i] ~ normal (mu_r, sig_r);
    EE[i] ~ normal(mu_e, sig_e);
    // log_f[i] ~ normal(logMu_f, logSig_f);
    
    for (j in 1:Con) {
      // draw data from probabilities determined by MMM parms
      count[j + (i-1)*Con,]  ~ multinomial(probs[j + (i-1)*Con,]);  
      
    }
  }
}

