#Everywhere: rows - observation, columns - genes

default_prior = "empirical"

rdirichlet_multinomial <- function(N, N_draws, alpha) {
  dirichlet_samples <- MCMCpack::rdirichlet(N, alpha)
  apply(dirichlet_samples, MARGIN = 1, FUN = function(x) { rmultinom(1, N_draws, x)[,1] })
}

get_prior <- function(observed, prior) {
  if(prior == "empirical") {
    prior <- alpha_guess_ml(observed)
    if(prior > 0.5) {
      warning("Empirical prior larger than 0.5")
    }
  } else if(!is.numeric(prior) || length(prior) != 1) {
    stop("Prior must be single number or 'empirical'")
  }
  prior
}

# Returns a matrix with N samples (rows) and length(observed) columns
sample_latent_dm_single <- function(N, observed, prior = default_prior) {
  posterior <- observed + get_prior(observed, prior)
  t(MCMCpack::rdirichlet(N, posterior))
}

# Returns a matrix with N samples (rows) and length(observed) columns
sample_posterior_dm_single <- function(N, observed, N_draws = sum(observed), prior = default_prior) {
  posterior <- observed + get_prior(observed, prior)
  t(rdirichlet_multinomial(N, N_draws, posterior))
}


# Returns array of size N * nrow(observed_matrix) * ncol(observed_matrix)
sample_posterior_dm_all <- function(N, observed_matrix, N_draws = max(rowSums(observed_matrix)), 
                                    prior = default_prior) {
  values_raw <- apply(observed_matrix, MARGIN = 1, FUN = function(x) { sample_posterior_dm_single(N, x, N_draws, prior) }) 
  wrong_dim_order <- array(values_raw, c(N, ncol(observed_matrix) , nrow(observed_matrix)), 
        dimnames = list(paste0("S",1:N), colnames(observed_matrix), rownames(observed_matrix)))
  aperm(wrong_dim_order, c(1,3,2))
}

#FUN takes a vector of values for single observation and returns a single number
summarise_posterior_per_observation <- function(posterior, FUN) {
  apply(posterior, MARGIN = c(1, 3), FUN = FUN)
}

#FUN takes a matrix of values for all observation and returns a vector number
summarise_posterior_per_sample <- function(posterior, FUN) {
  apply(posterior, MARGIN = c(1), FUN = FUN)
}


summarise_posterior_richness <- function(posterior) {
  summarise_posterior_per_observation(posterior, function(x) { sum(x > 0) })
}

dirichlet_multinomial_lpmf <- function(y, alpha) {
  if(length(y) != length(alpha)) {
    stop("y and alpha must have equal lengths")
  }
  alpha_plus <-  sum(alpha);
  sum_y <- sum(y)
  
  lgamma(sum_y + 1) + lgamma(alpha_plus) + sum(lgamma(alpha + y))
    - lgamma(alpha_plus+sum_y) - sum(lgamma(alpha)) - sum(lgamma(y + 1));
}

dirichlet_multinomial_lpmf_dalpha <- function(y, alpha) {
  if(length(alpha) != 1) {
    stop("Expecting single alpha")
  }
  k = length(y)
  k * digamma(alpha * k) + sum(digamma(alpha + y)) - k * digamma(alpha * k + sum(y)) -
    k * (digamma(alpha))
}

alpha_guess_ml <- function(observed) {
  lower <- 1e-12
  upper <- 0.5
  val_lower <- dirichlet_multinomial_lpmf_dalpha(observed, lower)
  val_upper <- dirichlet_multinomial_lpmf_dalpha(observed, upper)
  if(sign(val_lower) == sign(val_upper)) {
    warning("Empirical alpha outside of expected interval, using 0.5")
    0.5
  } else {
    uniroot(function(x) { dirichlet_multinomial_lpmf_dalpha(observed, x) },
          lower = lower, upper = upper)$root
  }
}