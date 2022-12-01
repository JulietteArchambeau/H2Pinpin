//Model to estimate the additive genetic variance in the common garden of Asturias

data {
  int<lower=1> N;                                                              // Number of observations
  vector[N] y;                                                                 // Response variable
  int<lower=0> nblock;                                                         // Number of blocks
  int<lower=0> nprov;                                                          // Number of provenances
  int<lower=0> nfam;                                                          // Number of families
  int<lower=0, upper=nblock> bloc[N];                                          // Blocks
  int<lower=0, upper=nprov> prov[N];                                           // Provenances
  int<lower=0, upper=nprov> which_prov[nfam];                                 // Provenances
  int<lower=0, upper=nfam> fam[N];                                           // Families
  vector[nprov] X;                                                             // Potential predictor of the additive genetic variance
}


parameters {
  real beta0; // global intercept
  real betaX; // coefficient of Xp
  simplex[4] pi;
  real<lower = 0> sigma_tot;
  real<lower = 0> sigma_K;
  vector[nprov] z_log_sigma_A;
  
  vector[nblock] z_block;
  vector[nprov] z_prov;
  vector[nfam] z_fam;
  
}


transformed parameters {
  real R_squared;
  
  real<lower = 0>  sigma_r;
  real<lower = 0>  sigma_block;
  real<lower = 0>  sigma_prov;
  
  real mean_sigma_fam;
  vector[nprov] sigma_fam;
  vector[nprov] sigma_A;
  
  vector[nprov] alpha_prov;
  vector[nfam] alpha_fam;
  vector[nblock] alpha_block;
  
  vector[N] mu; // linear predictor
  
  // variance partitioning with the simplex pi
  sigma_r = sqrt(pi[1]) * sigma_tot;
  sigma_block = sqrt(pi[2]) * sigma_tot;
  sigma_prov = sqrt(pi[3]) * sigma_tot;
  
  mean_sigma_fam= sqrt(pi[4]) * sigma_tot;
  sigma_A = exp(log(2*mean_sigma_fam) - (square(sigma_K)/2)  + betaX*X + z_log_sigma_A*sigma_K);
  
  alpha_prov = z_prov*sigma_prov;
  
  for(f in 1:nfam){
    alpha_fam[f] =  z_fam[f]*0.5*sigma_A[which_prov[f]]; 
    sigma_fam[which_prov[f]] = 0.5*sigma_A[which_prov[f]]; // sigma2_fam = 1/4 sigma2_A ==> sigma_fam = 1/2 sigma_A
  }
  
  alpha_block = z_block*sigma_block;
  
  mu = rep_vector(beta0, N) + alpha_fam[fam] + alpha_block[bloc] + alpha_prov[prov];
  R_squared = 1 - variance(y - mu) / variance(y);
}

model{
  //Priors
  beta0 ~ normal(mean(y),2);
  betaX ~ std_normal();
  sigma_tot ~ student_t(3, 0.0, 1.0);
  
  z_prov ~ std_normal();
  z_block ~ std_normal();
  z_fam ~ std_normal();
  
  z_log_sigma_A ~ std_normal();
  sigma_K ~ exponential(1);
  
  // Likelihood
  y ~ normal(mu, sigma_r);
}


generated quantities {
  //Variances
  real<lower=0> sigma2_r;
  real<lower=0> sigma2_prov;
  real<lower=0> sigma2_block;
  vector[nprov] sigma2_fam;
  vector[nprov] sigma2_A;
  vector[nprov] h2_prov;
  
  // Posterior predictive check
  vector[N] y_rep;
  
  // Log-Likelihood (for WAIC/loo computations)
  vector[N] log_lik;
  
  sigma2_r = square(sigma_r);
  sigma2_prov = square(sigma_prov);
  sigma2_fam = square(sigma_fam);
  sigma2_A = square(sigma_A);
  for(p in 1:nprov)  h2_prov[p] = sigma2_A[p]/(sigma2_r+sigma2_A[p]);
  sigma2_block = square(sigma_block);
  for(i in 1:N)  {
    y_rep[i] = normal_rng(mu[i], sigma_r);
    log_lik[i] = normal_lpdf(y[i]| mu[i], sigma_r); // log probability density function
  }
}

