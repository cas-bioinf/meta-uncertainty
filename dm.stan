functions {

  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    	real alpha_plus = sum(alpha);

      return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                  - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }

	int[] dirichlet_multinomial_rng(vector alpha, int N) {
	    return multinomial_rng(dirichlet_rng(alpha), N);
	}
}

data {
  int<lower=1> N;
  int<lower=0> counts[N];
}

parameters {
  real<lower=0> alpha_prior;
  //simplex[N] alpha;
}

model {
  alpha_prior ~ normal(0,0.25);
  //alpha ~ dirichlet(rep_vector(alpha_prior, N));
  counts ~ dirichlet_multinomial(rep_vector(alpha_prior, N));
}

generated quantities {
  int counts_gen[N];
  counts_gen = dirichlet_multinomial_rng(to_vector(counts) + alpha_prior, sum(counts));
}
