// Implemented in Stan by Jan Goetttmann 
// M3 - Model for extended encoding in complex span tasks after Oberauer & Lewandowsky, 2018
// Questions regarding the model code to Jan.Goettmann@psychologie.uni-heidelberg.de

// Best Model Equations for for Complex Span Task


data {
  int <lower=0> N;  // number rownumber
  int <lower=0> K;  // categories 
  int <lower=0> Con;  // number of categories
  int R[K];         // number of responses per category
  int count[N*Con,K];   // observed data
  real scale_b;     // set scaling for background noise
  vector[Con] Freetime;        // freetime conditions after distractor encoding (e.g. 500 ms, 1000 ms or 1500 ms) 
  real f;
}

parameters {
  // subject parameters
  real  a_raw[N];
  real  c_raw[N];
  real  r_raw[N];
  real  EE_raw[N];
  //real log_f_raw[N];
  
  // Mu & Sigma for hyper distributions
  real  mu_a;
  real <lower=0> sig_a;
  real  mu_c;
  real <lower=0> sig_c;
  real mu_r;
  real <lower=0> sig_r;
  real  mu_e;
  real <lower=0> sig_e;
  // real  logMu_f;
  // real <lower=0> logSig_f;
  
}

transformed parameters{
  
  real c[N];
  real a[N];
  real EE[N];
  real r[N];
  //real log_f[N];
  //real f[N];
  // real mu_f;
  // real sig_f;
  
  // activations
  real acts_IIP[N*Con];
  real acts_IOP[N*Con];
  real acts_DIP[N*Con];
  real acts_DIOP[N*Con];
  real acts_NPL[N*Con];
  
  
  // probabilities
  vector[K] probs[N*Con];
  real SummedActs[N*Con];
  
  
  for (i in 1:N)
  
  {
    
    c[i] = sig_c * c_raw[i] + mu_c;
    a[i] = sig_a  *a_raw[i] + mu_a;
    //log_f[i] = logSig_f * log_f_raw[i] + logMu_f; 
    EE[i] = mu_e + EE_raw[i] *sig_e;
    r[i] = mu_r + r_raw[i]  * sig_r;
   // f[i] =  inv_logit(log_f[i]);
    
  }
  
  
  // mu_f = inv_logit(logMu_f);
  // sig_f = sd(f);
  
  
  
  // loop over subjects and conditions to compute activations and probabilites
  
  for (i in 1:N){ // for each subject
  for(j in 1:Con) {
    
    // Best Fit for Complex Span Task
    
    acts_IIP[j + (i-1)*Con] = scale_b + ((1+EE[i]*Freetime[j])*c[i]) + a[i]; // Item in Position                      
    acts_IOP[j + (i-1)*Con] = scale_b + a[i];        // Item in Other Position
    acts_DIP[j + (i-1)*Con] = scale_b + f*((exp(-r[i]*Freetime[j])*c[i])+a[i]);// Distractor in Position
    acts_DIOP[j + (i-1)*Con] = scale_b + f *a[i]; // Distractor in other Position
    acts_NPL[j + (i-1)*Con] = scale_b; // non presented Lure
    
    SummedActs[j + (i-1)*Con] = R[1] * acts_IIP[j + (i-1)*Con] + R[2] * acts_IOP[j + (i-1)*Con] + R[3] * acts_DIP[j + (i-1)*Con]+ // account for different numbers auf DIP / DIOP 
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
  
  // logMu_f ~ normal(0,10);
  // logSig_f ~ gamma(1,0.01);
  // 
  
  // Draw subject parameters from truncated normal
  
  
  
  c_raw ~ normal(0,1);
  a_raw ~ normal(0,1);
  r_raw ~ normal(0,1);
  EE_raw ~ normal(0,1);
  // log_f_raw ~ normal(0, 1);
  
  for (j in 1:Con) {
    for (i in 1:N) {
      // draw data from probabilities determined by MMM parms
      count[j + (i-1)*Con,]  ~ multinomial(probs[j + (i-1)*Con,]);  
      
    }
  }
}

