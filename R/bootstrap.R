#Everywhere: rows - observation, columns - OTUs

default_dm_prior = "ml"

dm_model <- NULL

get_dm_model <- function() {
  if(is.null(dm_model)) {
    dm_model <- rstan::stan_model(paste0(path.package("metagenboot"),"/dm.stan"))
  }
  dm_model
}

rdirichlet_multinomial <- function(N, N_reads, alpha) {
  dirichlet_draws <- MCMCpack::rdirichlet(N, alpha)
  apply(dirichlet_draws, MARGIN = 1, FUN = function(x) { rmultinom(1, N_reads, x)[,1] })
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
    warning("ML alpha outside of expected interval, using 0.5")
    0.5
  } else {
    uniroot(function(x) { dirichlet_multinomial_lpmf_dalpha(observed, x) },
            lower = lower, upper = upper)$root
  }
}

alpha_guess_bayes <- function(N, observed) {
  fit <- sampling(get_dm_model(), data = list(N = length(observed), counts = observed),
                  iter = 1000 + max(N, 500), warmup = 1000)
  if(get_num_divergent(fit) > 0) {
    stop("divergences")
  }
  if(get_num_max_treedepth(fit) > 0) {
    stop("treedepth")
  }
  if(length(get_low_bfmi_chains(fit)) > 0) {
    stop("low bfmi")
  }
  s <- summary(fit, pars = "alpha_prior")$summary
  if( s[,"Rhat"] > 1.01) {
    warning("Large Rhat when guessing alpha")
  }
  if( s[,"n_eff"] < N / 2) {
    warning("Low n_eff when guessing alpha")
  }

  rstan::extract(fit, pars = "alpha_prior")$alpha_prior[1:N]
}

prepare_prior <- function(N, observed, prior) {
  if(prior == "ml") {
    prior <- alpha_guess_ml(observed)
    if(prior > 0.5) {
      warning("ML prior larger than 0.5")
    }
  } else if(prior == "bayes") {
    prior <- alpha_guess_bayes(N, observed)
  } else if(!is.numeric(prior) ||
            (length(prior) != 1 && length(prior) != N)) {
    stop("Prior must be single number, vector of length N, 'ml' or 'bayes'")
  }

  if(length(prior) == 1) {
    prior <- rep(prior, N)
  }

  prior
}

get_dirichlet_proportions <- function(N, observed, prior = default_dm_prior) {
  prior_prep <- prepare_prior(N, observed, prior)
  proportions <- matrix(NA_real_, ncol = length(observed), nrow = N)
  for(i in 1:N) {
    proportions[i,] <- MCMCpack::rdirichlet(1, observed + prior_prep[i])
  }
  # reads <- MCMCpack::rdirichlet(N, observed + prior_prep[1])
  proportions
}

# Returns a matrix with N draws (rows) and length(observed) columns
bootstrap_proportions_dm_single <- function(N, observed, prior = default_dm_prior) {
  t(get_dirichlet_proportions(N, observed, prior))
}

# Returns a matrix with N samples (rows) and length(observed) columns
bootstrap_reads_dm_single <- function(N, observed, N_reads = sum(observed), prior = default_dm_prior) {
  if(N_reads == "original") {
    N_reads = sum(observed)
  } else if(!is.numeric(N_reads) || length(N_reads) != 1) {
    stop("N_reads must be single number or 'original'")
  }

  dirichlet_proportions <- get_dirichlet_proportions(N, observed, prior)
  t(apply(dirichlet_proportions, MARGIN = 1, FUN = function(x) { rmultinom(1, N_reads, x)[,1] }))
}


# Returns 3D array with dimensions N, nrow(observed_matrix), ncol(observed_matrix)
bootstrap_reads_dm <- function(N, observed_matrix, N_reads = "original",
                                    prior = default_dm_prior) {
  #TODO check observed_matrix is just integers
  values_raw <- apply(observed_matrix, MARGIN = 1, FUN = function(x) { bootstrap_reads_dm_single(N, x, N_reads, prior) })
  wrong_dim_order <- array(values_raw, c(N, ncol(observed_matrix) , nrow(observed_matrix)),
        dimnames = list(paste0("S",1:N), colnames(observed_matrix), rownames(observed_matrix)))
  aperm(wrong_dim_order, c(1,3,2))
}

# Returns 3D array with dimensions N, nrow(observed_matrix), ncol(observed_matrix)
bootstrap_proportions_dm <- function(N, observed_matrix,
                                    prior = default_dm_prior) {
  values_raw <- apply(observed_matrix, MARGIN = 1, FUN = function(x) { bootstrap_proportions_dm_single(N, x, prior) })
  wrong_dim_order <- array(values_raw, c(N, ncol(observed_matrix) , nrow(observed_matrix)),
                           dimnames = list(paste0("S",1:N), colnames(observed_matrix), rownames(observed_matrix)))
  aperm(wrong_dim_order, c(1,3,2))
}


bootstrap_reads_DESeq2 <- function(N, observed_matrix, mapping, design) {
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = t(observed_matrix), colData = mapping, design = design)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)

  n_otus <- ncol(observed_matrix)
  n_observations <- nrow(observed_matrix)

  predicted_means <- 2^(coef(dds) %*% t(model.matrix(design, mapping ))) * matrix(rep(sizeFactors(dds), each = n_otus), n_otus, n_observations)
  dispersions <- matrix(rep(dispersions(dds), times = n_observations))

  res <- array(NA_integer_, c(N, n_observations, n_otus),
               dimnames = list(paste0("S",1:N), rownames(observed_matrix), colnames(observed_matrix)))
  for(n in 1:N) {
    res[n,,] <- t(matrix(
      rnbinom(n_observations * n_otus, mu = predicted_means, size = dispersions), n_otus, n_observations))
  }
  res
}

bootstrap_reads_rarefy <- function(N, observed_matrix, min_depth = min(rowSums(observed_matrix))) {
  n_otus <- ncol(observed_matrix)
  n_observations <- nrow(observed_matrix)

  res <- array(NA_integer_, c(N, n_observations, n_otus),
               dimnames = list(paste0("S",1:N), rownames(observed_matrix), colnames(observed_matrix)))
  for(n in 1:N) {
    res[n,,] <- rrarefy(observed_matrix, min_depth)
  }
  res
}

bootstrap_reads_jackknife_observations <- function(observed_matrix) {
  n_otus <- ncol(observed_matrix)
  n_observations <- nrow(observed_matrix)

  res <- list()
  for(n in 1:n_observations) {
    res[[n]] <- observed_matrix[-n,, drop = FALSE]
  }
  attr(res, "observation_subset") <- TRUE
  res
}

#Turns value returned by `bootstrap_reads_dm` into a matrix of size N * nrow(observed_matrix), ncol(observed_matrix)
#Useful for clustering etc.
flatten_posterior_draws <- function(posterior_draws, name.delim = "_") {
  if(is.list(posterior_draws)) {
    stop("flatten_posterior_draws not yet implemented for list draws")
  }
  dims <- dim(posterior_draws)
  flat <- matrix(posterior_draws, dims[1] * dims[2], dims[3])
  combined_names <- expand.grid(dimnames(posterior_draws)[[1]], dimnames(posterior_draws)[[2]])
  rownames(flat) <- paste0(combined_names[[2]], name.delim, combined_names[[1]])
  colnames(flat) <- dimnames(posterior_draws)[[3]]
  flat
}

#FUN takes a vector of values for single observation and returns a single number
summarise_bootstrap_per_observation <- function(bootstrap, FUN) {
  if(is.list(posterior)) {
    stop("summarise_posterior_per_observation not yet implemented for list draws")
  }
  t(apply(bootstrap, MARGIN = c(1, 2), FUN = FUN))
}

#FUN takes a matrix of values for all observation and returns a vector number
summarise_bootstrap_per_draw <- function(bootstrap, FUN) {
  if(is.list(bootstrap)) {
    stop("summarise_posterior_per_draw not yet implemented for list draws")
  }
  apply(bootstrap, MARGIN = c(1), FUN = FUN)
}

#FUN takes a mtrix of values and returns any object - a list with result per draw is returned
apply_per_draw <- function(bootstrap, FUN, ...) {
  if(is.list(bootstrap)) {
    stop("apply_per_draw not yet implemented for list draws")
  }
  plyr::alply(bootstrap, .margins = 1, .fun = FUN, ...)
}

summarise_bootstrap_richness <- function(bootstrap) {
  summarise_bootstrap_per_observation(bootstrap, function(x) { sum(x > 0) })
}

get_n_bootstrap_draws <- function(bootstrap) {
  if(is.list(bootstrap)) {
    n_draws <- length(bootstrap)
  } else {
    n_draws <- dim(bootstrap)[1]
  }
  n_draws
}

get_bootstrap_draw <- function(bootstrap, draw_i) {
  if(is.list(bootstrap)) {
    bootstrap[[draw_i]]
  } else {
    bootstrap[draw_i,,]
  }
}

is_observation_subset <- function(bootstrap) {
  attr_val <- attr(bootstrap, "observation_subset", exact = TRUE)
  if(is.null(attr_val)) {
    FALSE
  } else {
    attr_val
  }
}
