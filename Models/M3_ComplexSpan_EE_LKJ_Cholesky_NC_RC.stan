// Implemented in Stan by Jan Goetttmann 
// M3 - Model for extended encoding in complex span tasks after Oberauer & Lewandowsky, 2018
// Questions regarding the model code to Jan.Goettmann@psychologie.uni-heidelberg.de

// Best Fitting Model for Retro-Cue Experiments (EXP 3 & $ Oberauer et al., 2018)


functions{
  // flatten_lower_tri: function that returns the lower-tri of a matrix, flattened to a vector
  
  vector flatten_lower_tri(matrix mat) {
    int n_cols = cols(mat) ;
    int n_uniq = (n_cols * (n_cols - 1)) %/% 2;
    vector[n_uniq] out ;
    int i = 1;
    for(c in 1:(n_cols-1)){
      for(r in (c+1):n_cols){
        out[i] = mat[r,c];
        i += 1;
      }
    }
    return(out) ;
  }
  
}

data {
  int <lower=0> N;  // number rownumber
  int <lower=0> K;  // categories 
  int <lower=0> J;  // Dims of Cov Matrix
  int <lower=0> Con;  // number of FT categories
  int R[K];         // number of responses per category
  int count[N*Con,K];   // observed data
  real scale_b; // set scaling for background noise
  vector[Con] Freetime;  // freetime conditions after distractor encoding (e.g. 500 ms, 1000 ms or 1500 ms) 
  int retrievals;
  //matrix[N*Con,2] d_weight; // Weights for Response Category inbalances
}

parameters {
  
  // Defining vector for hyper and subject parameters 
  
  cholesky_factor_corr[J] L_Omega;
  vector<lower=0>[J] sigma;
  vector[J] hyper_pars;
  matrix[J,N] theta;
  
}


transformed parameters {
  // non-centered multivariate
  matrix[J,N] subj_pars =  (
    diag_pre_multiply( sigma, L_Omega )
    * theta
    + rep_matrix(hyper_pars,N)
    ) ;
    
    // Transform f Parameter
    real mu_f = inv_logit(hyper_pars[3]);
    row_vector[N] f = inv_logit(subj_pars[3,]);
    
    // build compound of r and f 
    real ac[N];
    
    for (i in 1:N)
    ac[i]= exp(f[i]-1) * subj_pars[5,i]; 
    
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
      
      acts_IIP[j + (i-1)*Con] = scale_b + ((1+subj_pars[4,i]*Freetime[j])* subj_pars[1,i]) + subj_pars[2,i]; // Item in Position                      
      acts_IOP[j + (i-1)*Con] = scale_b + (1+subj_pars[4,i]*Freetime[j])*subj_pars[2,i];        // Item in Other Position
      acts_DIP[j + (i-1)*Con] = scale_b + f[i]*((exp(-subj_pars[5,i]*Freetime[j])*subj_pars[1,i]) + subj_pars[2,i]); // Distractor in Position
      acts_DIOP[j + (i-1)*Con] = scale_b + f[i]*subj_pars[2,i]; // Distractor in other Position
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
  
  hyper_pars[1] ~ normal(20,10); // c
  hyper_pars[2] ~ normal(2,10); // a
  hyper_pars[3] ~ normal(0,10); // f
  hyper_pars[4] ~ normal(1,10); // EE
  hyper_pars[5] ~ normal(1,10); // r
  
  
  // 
  L_Omega ~ lkj_corr_cholesky(2);
  sigma ~ gamma(1,0.01);
  
  
  // Loop over subjects
  
  for (i in 1:N) 
  {
    
    theta[,i] ~ normal(0,1);
    
  }
  
  
  
  for (j in 1:Con) {
    for (i in 1:N) {
      // draw data from probabilities determined by MMM parms
      count[j + (i-1)*Con,]  ~ multinomial(probs[j + (i-1)*Con,]);  
    }
  }
}

generated quantities{
  
  vector[(J*(J-1))%/%2] cor_mat_lower_tri;
  int count_rep[N*Con,K];
  
  
  cor_mat_lower_tri = flatten_lower_tri(multiply_lower_tri_self_transpose(L_Omega));
  
  
  
  
  for (i in 1:N)
  for(j in 1:Con)
  {
    {
      
      count_rep[j + (i-1)*Con,] = multinomial_rng(probs[j + (i-1)*Con,], retrievals);
      
    }
  }
}

