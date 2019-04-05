functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    	real alpha_plus = sum(alpha);
      return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                  - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
}

data {
  int<lower=1> N;
  int<lower=0> counts[N];
}

parameters {
  real<lower=0> alpha_prior;
}

model {
  alpha_prior ~ normal(0,0.25);
  counts ~ dirichlet_multinomial(rep_vector(alpha_prior, N));
}
