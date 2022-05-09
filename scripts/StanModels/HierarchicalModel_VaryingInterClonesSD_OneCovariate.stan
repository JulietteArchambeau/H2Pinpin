//Model 1:  model with the betaX * Xp on the median of the sigma_clon

data {
  int<lower=1> N;                                                              // Number of observations
  vector[N] y;                                                                 // Response variable
  int<lower=0> nblock;                                                         // Number of blocks
  int<lower=0> nprov;                                                          // Number of provenances
  int<lower=0> nclon;                                                          // Number of clones
  int<lower=0, upper=nblock> bloc[N];                                          // Blocks
  int<lower=0, upper=nprov> prov[N];                                           // Provenances
  int<lower=0, upper=nprov> which_prov[nclon];                                 // Provenances
  int<lower=0, upper=nclon> clon[N];                                           // Clones
  vector[nprov] X;                                                             // Potential predictor of sigma_clon
}


parameters {
  real beta0; // global intercept
  real betaX; // coefficient of Xp
  simplex[4] pi;
  real<lower = 0> sigma_tot;
  real<lower = 0> sigma_K;
  vector[nprov] z_log_sigma_clon;
  
  vector[nblock] z_block;
  vector[nprov] z_prov;
  vector[nclon] z_clon;
  
}


transformed parameters {
  real R_squared;
  
  real<lower = 0>  sigma_r;
  real<lower = 0>  sigma_block;
  real<lower = 0>  sigma_prov;
  
  real mean_sigma_clon;
  vector[nprov] sigma_clon;
  
  vector[nprov] alpha_prov;
  vector[nclon] alpha_clon;
  vector[nblock] alpha_block;
  
  vector[N] mu; // linear predictor
  
  // variance partitioning with the simplex pi
  sigma_r = sqrt(pi[1]) * sigma_tot;
  sigma_block = sqrt(pi[2]) * sigma_tot;
  sigma_prov = sqrt(pi[3]) * sigma_tot;
  
  mean_sigma_clon= sqrt(pi[4]) * sigma_tot;
  sigma_clon = exp(log(mean_sigma_clon) - (square(sigma_K)/2)  + betaX*X + z_log_sigma_clon*sigma_K);
  
  alpha_prov = z_prov*sigma_prov;
  
  for(c in 1:nclon){
    alpha_clon[c] =  z_clon[c]*sigma_clon[which_prov[c]];
  }
  
  alpha_block = z_block*sigma_block;
  
  mu = rep_vector(beta0, N) + alpha_clon[clon] + alpha_block[bloc] + alpha_prov[prov];
  R_squared = 1 - variance(y - mu) / variance(y);
}

model{
  //Priors
  beta0 ~ normal(mean(y),2);
  betaX ~ std_normal();
  sigma_tot ~ student_t(3, 0.0, 1.0);
  
  z_prov ~ std_normal();
  z_block ~ std_normal();
  z_clon ~ std_normal();
  
  z_log_sigma_clon ~ std_normal();
  sigma_K ~ exponential(1);
  
  // Likelihood
  y ~ normal(mu, sigma_r);
}


generated quantities {
  //Variances
  real<lower=0> sigma2_r;
  real<lower=0> sigma2_prov;
  real<lower=0> sigma2_block;
  vector[nprov] sigma2_clon;
  vector[nprov] h2_prov;

  // Posterior predictive check
  vector[N] y_rep;

  // Log-Likelihood (for WAIC/loo computations)
  vector[N] log_lik;

  sigma2_r = square(sigma_r);
  sigma2_prov = square(sigma_prov);
  sigma2_clon = square(sigma_clon);
  for(p in 1:nprov)  h2_prov[p] = sigma2_clon[p]/(sigma2_r+sigma2_clon[p]);
  sigma2_block = square(sigma_block);
  for(i in 1:N)  {
    y_rep[i] = normal_rng(mu[i], sigma_r);
    log_lik[i] = normal_lpdf(y[i]| mu[i], sigma_r); // log probability density function
  }
}

