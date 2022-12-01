//Model code to estimate the broad-sense heritability for each trait

data {
  int<lower=1> N;                                                              // Number of observations
  vector[N] y;                                                                 // Response variable
  int<lower=0> nblock;                                                         // Number of blocks
  int<lower=0> nprov;                                                          // Number of provenances
  int<lower=0> nclon;                                                          // Number of clones
  int<lower=0, upper=nblock> bloc[N];                                          // Blocks
  int<lower=0, upper=nprov> prov[N];                                           // Provenances
  int<lower=0, upper=nclon> clon[N];                                           // Clones
}


parameters {
  real beta0; // global intercept
  simplex[4] pi;
  real<lower = 0> sigma_tot;
  
  vector[nblock] z_block;
  vector[nprov] z_prov;
  vector[nclon] z_clon;
  
}


transformed parameters {
  real R_squared;
  
  real<lower = 0>  sigma_r;
  real<lower = 0>  sigma_block;
  real<lower = 0>  sigma_prov;
  real<lower = 0>  sigma_clon;
  
  vector[nprov] alpha_prov;
  vector[nclon] alpha_clon;
  vector[nblock] alpha_block;
  
  vector[N] mu; // linear predictor
  
  // variance partitioning with the simplex pi
  sigma_r = sqrt(pi[1]) * sigma_tot;
  sigma_block = sqrt(pi[2]) * sigma_tot;
  sigma_prov = sqrt(pi[3]) * sigma_tot;
  sigma_clon= sqrt(pi[4]) * sigma_tot;
  
  alpha_prov = z_prov*sigma_prov;
  alpha_clon = z_clon*sigma_clon;
  alpha_block = z_block*sigma_block;
  
  mu = rep_vector(beta0, N) + alpha_clon[clon] + alpha_block[bloc] + alpha_prov[prov];
  R_squared = 1 - variance(y - mu) / variance(y);
}

model{
  //Priors
  beta0 ~ normal(mean(y),2);
  sigma_tot ~ student_t(3, 0.0, 1.0);
  
  z_prov ~ std_normal();
  z_block ~ std_normal();
  z_clon ~ std_normal();
  
  // Likelihood
  y ~ normal(mu, sigma_r);
}


generated quantities {
  //Variances and broad-sense heritability
  real<lower=0> sigma2_r;
  real<lower=0> sigma2_prov;
  real<lower=0> sigma2_block;
  real<lower=0> sigma2_clon;
  real<lower=0> H2;
  
  sigma2_r = square(sigma_r);
  sigma2_prov = square(sigma_prov);
  sigma2_clon = square(sigma_clon);
  H2 = sigma2_clon/(sigma2_r+sigma2_clon);
  sigma2_block = square(sigma_block);

}

