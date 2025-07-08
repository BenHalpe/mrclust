loglik <- function(m, pi_clust, theta, theta_sd, theta_clust,
                   junk_mixture, df, mu, sig,
                   null_mixture, mu_null, sig_null) {
  i <- 1:m # the theta parameter is a vector of all variants, not a single one
  if (junk_mixture | null_mixture) {
    k <- length(theta_clust) - sum(junk_mixture) - sum(null_mixture)
  } else {
    k <- length(theta_clust)
  }
  tmp_lg_den <- log_norm(i, j = 1:k, m, tmp_lgt2 = k, theta, theta_sd,
                         theta_clust)  # tmp_lg_den[i,j] is the log likelihood of variant i to be in cluster j
  if (null_mixture) {
    tmp_lg_den <- c(tmp_lg_den, null_den(x = theta, mu = mu_null,
                                          sig = sig_null, log = TRUE))
  }
  if (junk_mixture) {
    tmp_lg_den <- c(tmp_lg_den, gen_t(theta, df = df, mu = 0, sig = sig,
                                       log = TRUE)) # same thing as for the substansive clusters, but we use t distribution instead of gaussian
  }

  tmp <- exp(matrix(rep(log(pi_clust), each = m) + tmp_lg_den, nrow = m,
                    ncol = length(pi_clust))) # multiplies the log likelihood for each cluster by its proprtion to get the overall probabilities - done by adding the log proportions to the log likelihood and taking the exponent.
  tmp_rowsum <- rowSums(tmp) # probability per variant
  res <- sum(log(tmp_rowsum)) # total log likelihood of dataset
  return(res)
}
