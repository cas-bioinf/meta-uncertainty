synthetic_convex_mixture <- function(n_otus, n_observations, n_populations = 2,
                                     mean_reads = 10 * n_otus,
                                     mixing_prior = 1 / (n_populations - 1),
                                     include_raw_populations = FALSE) {
  if(mean_reads <= n_otus) {
    stop("mean_reads has to be greater than n_otus")
  }

  background <- rep(1e-3, n_otus)
  #Two complementary vectors
  populations <- matrix(NA_real_, nrow = n_populations, ncol = n_otus)
  for(p in 1:n_populations) {
    code_vector <- rep(0, n_populations)
    code_vector[p] <- 1
    populations[p,] <- rep(code_vector, length.out = n_otus) * rexp(n_otus)
  }

  # Rows are observations
  mixing <- MCMCpack::rdirichlet(n_observations, rep(mixing_prior, n_populations))

  if(include_raw_populations) {
    mixing[1:n_populations,] <- diag(1, nrow = n_populations, ncol = n_populations)
  }

  mapping <- data.frame(Sample = as.character(1:n_observations),
                        mixing = mixing,
                        n_reads = n_otus + rnbinom(n_observations, mu = (mean_reads - n_otus), size = 1),
                        Group = "A"
                        )

  otu_t <- matrix(NA_integer_, nrow = n_observations, ncol = n_otus)
  rownames(otu_t) <- as.character(1:n_observations)
  for(o in 1:n_observations) {
    prob <- background + mixing[o,, drop = FALSE] %*% populations
    #otu_t[o,] <- round(mapping$n_reads[o] * prob)
    otu_t[o,] <- rmultinom(1, mapping$n_reads[o], prob)
  }
  list(otu_t = otu_t, mapping = mapping)
}
