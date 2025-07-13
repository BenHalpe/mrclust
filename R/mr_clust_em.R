#' MR-Clust mixture model fitting
#'
#' Assessment of clustered heterogeneity in Mendelian randomization analyses
#' using expectation-maximisation (EM) based model fitting of the MR-Clust
#' mixture model. Function output includes both data-tables and a visualisation
#' of the assingment of variants to clusters.
#'
#' @param theta numeric vector of length the number of variants,
#' the i-th element is a ratio-estimate for the i-th genetic variant.
#' @param theta_se numeric vector of length the number of variants,
#' the i-th element is the standard error of the ratio-estimate
#' for the i-th genetic variant.
#' @param bx numeric vector of length the number of variants,
#' the i-th element is the estimated regression coefficient -
#' i.e. beta-x value - relating the i-th genetic variant to the risk-factor.
#' @param by numeric vector of length the number of variants, the i-th element
#' is the estimated regression coefficient - i.e. beta-y value - relating the
#' i-th genetic variant to the outcome.
#' @param bxse numeric vector of length the number of variants, the i-th element
#'  is the standard error of the estimated regression coefficient relating the
#'  i-th genetic variant to the risk-factor.
#' @param byse numeric vector of length the number of variants, the i-th element
#'  is the standard error of the estimated regression coefficient relating the
#'  i-th genetic variant to the outcome.
#' @param obs_names character vector of length the number of variants,
#' the i-th element is the name of the i-th genetic variants - e.g. the rsID.
#' @param max_iter numeric integer denoting the maximum number of iterations to
#' take before stopping the EM-algorithm's search for a maxima in the
#' log-likelihood.
#' @param tol numeric scalar denoting the maximum absolute difference between
#' two computations of the log-likelihood with which we accept that a maxima in
#' the log-likelihood has been computed.
#' @param junk_sd numeric scalar denoting the scale parameter in the generalised
#'  t-distribution
#' @param junk_mean numeric scalar denoting the mean of the generalised
#' t-distribution. By default mean is set to zero.
#' @param min_clust_search numeric integer which denotes the minimum number of
#' clusters searched for in the data - default computes evidence supporting up
#' to K=10 clusters which might explain any clustered heterogeneity in the data.
#' @param stop_bic_iter numeric integer I, for computational efficiency -
#' particularly when analysing large numbers of variants - we can stop the
#' EM-algorithm if the BIC is monotonic increasing over the previous I increases
#'  in the number of clusters K. By default evidence supporting at least 10
#'  clusters in the data is computed and so, for example, if the BIC from models
#'   which assume 6 clusters; 7 clusters; ... or; 10 clusters is monotonic
#'   increasing - in the number of clusters K -then the EM-algorithm is stopped
#'   and the model whose K minimises the BIC is returned.
#' @param results_list character list allowing users to choose whether to return
#'  a table with the variants assigned to: "all" of the clusters; a single
#'  "best" cluster or; both. By default we return both, i.e. results_list =
#'  list("all", "best").
#' @param cluster_membership numeric list which allows users to output a list
#' which, for each cluster, returns the variants assigned to the cluster by
#' stratified by the probability of belonging to the cluster. By default,
#' cluster_membership = list(by_prob = 0.1, bound = 0); so that MRClust
#' returns a list, which for each cluster, outputs the variants assigned
#' to the cluster with probability between (0.9,1); (0.8,0.9);... and finally;
#' (0.1,0), i.e. by probability increments 0.1 from 1 to a lower bound of 0.
#' @param plot_results numeric list which allows users to plot the output of
#' MRClust. By default, plot_results = list("best", min_pr = 0.5); so that the
#' best clustering is plotted with variants assigned to a cluster with
#' probability above 0.5.
#' @param trait_search logical, for each of the non-null and non-junk clusters
#' search phenoscanner for traits associated with the variants.
#' @param trait_pvalue numeric scalar for use with trait_search, representing
#' the maximum p-value with with at least one variant in the cluster must be
#' associated with a trait for it to be returned in the phenoscanner search.
#' Default value is GWA significance, i.e. 5*10^-8.
#' @param proxy_r2 numeric scalar for use with trait search, allowing variants
#' whose r2>=proxy_r2 to be included in the trait search. Default r2=0.8.
#' @param catalogue character, for use with trait search. From Phenoscanner
#' (http://www.phenoscanner.medschl.cam.ac.uk/information/) "the catalogue to be
#'  searched (options: None, GWAS, eQTL, pQTl, mQTL, methQTL)". Default setting
#'  is catalogue = "GWAS".
#' @param proxies character, for use with trait search. From Phenoscanner
#' (http://www.phenoscanner.medschl.cam.ac.uk/information/) "the proxies
#' database to be searched (options: None, AFR, AMR, EAS, EUR, SAS)".
#' Default setting is proxies = "None"
#' @param build integer, for use with trait search. From Phenoscanner
#' (http://www.phenoscanner.medschl.cam.ac.uk/information/) "Human genome
#' build numbers (options: 37, 38; default: 37)". Default setting is build = 37.
#' @return Returned are: estimates of the putative number of clusters in the
#' sample, complete with allocation probabilities and summaries of the
#' association estimates for each variant; plots which visualise the allocation
#' of variants to clusters and; several summaries of the fitting process, i.e.
#' the BIC and likelihood estimates.
#' @export
mr_clust_em <- function(theta, theta_se, bx, by, bxse, byse,
                        obs_names = NULL, max_iter = 5e3, tol = 1e-5,
                        junk_sd = NULL, junk_mean = 0,
                        stop_bic_iter = 5, min_clust_search = 10,
                        results_list = list("all", "best"),
                        cluster_membership = list(by_prob = 0.1, bound = 0),
                        plot_results = list("best", min_pr = 0.5),
                        trait_search = FALSE, trait_pvalue = 1e-5,
                        proxy_r2 = 0.8, catalogue = "GWAS", proxies = "None",
                        build = 37) {

  # previous user defined options now hard-coded
  cluster_sizes <- NULL
  init_clust_means <- NULL
  init_clust_probs <- NULL
  junk_mixture <- TRUE
  df <- 4
  junk_prob <- NULL
  junk_null_aware_bic <- TRUE
  fix_junk_prob <- FALSE
  null_mixture <- TRUE
  null_sd <- NULL
  null_prob <- NULL
  fix_null_prob <- FALSE
  scale_grid_search <- FALSE
  grid_increment <- 0.1
  grid_max <- 2
  clust_size_prior <- FALSE
  bic_prior <- NULL
  rand_num <- 5
  rand_sample <- seq(0.05, 0.4, by = 0.05)

  # initial parameter values

  m <- length(theta)
  cluster_sizes <- 0:m # MISLEADING NAME! this is an array of the possible numbers of substansive clusters. 
  if (is.null(obs_names)) {
    obs_names <- paste0("snp_", 1:m)
  }

  num_clust <- length(cluster_sizes)
  num_clust <- num_clust * rand_num # rand_num is the number of repeats the algorithm will have in order to avoid local optima

  bic_clust <- bic_clust_mx <- vector() # both initialized here. bic_clust stores all bics, bic_clust_mx stores the best bic across random initialization for each number of clusters.
  clust_mn <- vector("list", num_clust)
  clust_pi <- vector("list", num_clust)
  log_like <- vector("list", num_clust)


  count0 <- count_k <- 1

  for (i in cluster_sizes) { # for each selection of number of clusters - this is for the BIC selection
    for (itr in 1:(rand_num + 1)) { #for each random initalization of the k means algorithm
      if (is.null(init_clust_means)) {
        if (i > 0 & i != m) { #substantive clusters
          init_conds <- stats::kmeans(x = theta, centers = i, iter.max = 5e3) #as in section 3.2.1, the initial centers/proportions are calculated using k-means.
          clust_means <- as.numeric(init_conds$centers)
          # init_conds$clusters is the assignment of cluster per variant. table counts the number of appearances of each cluster. dividing by the number of cluster gives us the proprtion pi of each cluster
          clust_probs <- table(init_conds$cluster) / m 
          # printing for debugging, like in matlab
          init_conds
          clust_means
          clust_probs
        } else if (i == m) { # if the number of clusters is the number of variants, set each cluster to include only the variant and uniform proportion
          clust_means <- theta
          clust_probs <- rep(1, m) / m
        } else if (i == 0) { 
          clust_means <- clust_probs <- NULL
        }
      }

      # junk distribution parameters

      if (junk_mixture & is.null(junk_sd)) { #computes the standard deviation of the junk cluster if not directly given by user
        # we want a width that will allow us to capture distant points
        rng_thet <- range(theta) # [min(theta), max(theta)]
        max_disp <- which.max(abs(theta) + 2 * theta_se) # this finds the variant which maximizes how far it is from the center
        sig <- (rng_thet[2] - rng_thet[1] + theta_se[max_disp]) # this is the computed std, it is the range of all thetas with the added safety distance of the max dispersion
      } else if (junk_mixture & !is.null(junk_sd)) {
        sig <- junk_sd
      } else {
        sig <- NULL
      }
      if (junk_mixture & is.null(junk_mean)) {
        mu <- 0
      } else if (junk_mixture & !is.null(junk_mean)) {
        mu <- junk_mean
      } else {
        mu <- NULL
      }
      # null distribution parameters
      # abs(theta / theta_se) is the z score for the null, which measures by how many standard deviation each theta deviates from the mean - the mean in the null cluster is 0 so we just use theta
      # this returns the indices of variants that are marked to belong to the null cluster, since they are statistically close to 0 with high confidence. 
      null_obs <- abs(theta / theta_se) < 1.96
      if (null_mixture & is.null(null_sd)) {
        sig_null <- theta_se
        mu_null <- 0
      } else if (null_mixture & !is.null(null_sd)) {
        sig_null <- null_sd
        mu_null <- 0
      } else {
        sig_null <- mu_null <- NULL
      }

      ####
      if (itr == 1) {
        if (junk_mixture | null_mixture) {
          ### junk_mixture parameters
          sum_pr <- 0
          # next line finds how many variants are junk, by measuring how many are statistically extremely far from the median value.
          # pnorm(y,params) is how likely a sampling of a gaussian with params is to return something smaller than y
          # so we are counting variants where it is less than 5% to be further from the others than
          junk_obs <- sum(2 * (1 - stats::pnorm(theta - stats::median(theta),
                                                0, theta_se)) < 0.05)
          # jk_pr is the proprtion of the junk cluster
          if (junk_mixture & !fix_junk_prob) {
            if (sum(junk_obs) == 0) {
              jk_pr <- 1 / m  # if no extreme outliers were observed, the proprtion of the junk cluster is uniform
              sum_pr <- jk_pr 
            } else {
              jk_pr <- sum(junk_obs) / m # if there are outliers they make up the proprtion of the junk cluster as done with the substantive clusters
              sum_pr <- jk_pr
            }
          } else if (junk_mixture & fix_junk_prob) { # fix as in use user defined value instead of calculating
            jk_pr <- junk_prob
            sum_pr <- jk_pr
          } else {
            jk_pr <- NULL
          }

          #### null_mixture parameters
          if (null_mixture & !fix_null_prob) {
            if (is.null(jk_pr)) {
              tmp_jk_pr <- 0
            } else {
              tmp_jk_pr <- jk_pr
            }
            if (sum(null_obs) == 0) {
              null_pr <- 1 / m
              sum_pr <- sum_pr + null_pr
            } else {
              null_pr <- sum(null_obs) / m # see similar logic fot jk_pr above
              if (null_pr + tmp_jk_pr > (1 - 1 / m)) { # in case the junk and null clusters are proportionally most of the clusters with substantive clusters being very small, fix them to avoid the algorithm treating everything like noise
                tmp_sum <- null_pr + tmp_jk_pr
                # reassign the null and junk proportions by normalizing them with respect to each other and taking extra safety margin of 1/m over all non substantive clusters
                if (is.null(jk_pr)) {
                  null_pr <- null_pr / tmp_sum - 1 / m
                } else {
                  null_pr <- null_pr / tmp_sum - 1 / (2 * m)
                  tmp_jk_pr <- tmp_jk_pr / tmp_sum - 1 / (2 * m)
                  jk_pr <- tmp_jk_pr
                  sum_pr <- tmp_jk_pr
                }
              }
              sum_pr <- sum_pr + null_pr
            }
          } else if (null_mixture & fix_null_prob) {
            null_pr <- null_prob
            sum_pr <- sum_pr + null_pr
          } else {
            null_pr <- NULL
          }
          # c() just concatenates, mu here is the junk mean
          clust_means <- c(clust_means, mu_null, mu)
          clust_probs <- c(clust_probs * (1 - sum_pr), null_pr, jk_pr) # rescale the substantive cluster proprtions to account for null and junk clusters
          k_clust <- cluster_sizes[count_k] + sum(junk_mixture) +
            sum(null_mixture)
        } else {
          k_clust <- cluster_sizes[count_k]
          sig <- mu <- NULL
          sig_null <- mu_null <- NULL
          junk_null_aware_bic <- FALSE
        }
      } else {
        if (null_mixture) {
          null_pr <- sample(rand_sample, 1)
          sum_pr <- null_pr
        } else {
          null_pr <- NULL
          sum_pr <- 0
        }
        if (junk_mixture) {
          jk_pr <- sample(rand_sample, 1)
          sum_pr <- sum_pr + jk_pr
        } else {
          jk_pr <- NULL
        }

        clust_means <- c(clust_means, mu_null, mu)
        clust_probs <- c(clust_probs * (1 - sum_pr), null_pr, jk_pr)
        k_clust <- cluster_sizes[count_k] + sum(junk_mixture) +
          sum(null_mixture) # each sum will be 1 if there are any variants in the junk/null cluster so will flag the existance of the cluster
      }
      #### Expectation Maximization step

      if (scale_grid_search & junk_mixture) { # runs the E-M step multiple times for multiple guesses at the junk cluster standard deviation in order to find the one that leads to the highest log likelihood
        sigs <- sig * seq(1, grid_max, by = grid_increment)
        pi_tmp <- vector("list", length(sigs))
        theta_tmp <- vector("list", length(sigs))
        loglik_tmp <- vector("numeric", length(sigs))
        count2 <- 1

        for (inc in sigs) {
          sig_tmp <- inc
          pi_clust <- clust_probs # initialise cluster class probabilities
          theta_clust <- clust_means # initialise cluster centroids; these can
          #be the same as the std error is assumed to vary between the
          #observations
          count <- 1 # start iteration count
          loglik_diff <- 1
          loglik_iter <- NULL
          while (loglik_diff > tol & count < max_iter) {
            tmp_loglik0 <- loglik(
              m, pi_clust, theta, theta_se, theta_clust, junk_mixture, df, mu,
              sig_tmp, null_mixture, mu_null, sig_null
            )
            tmp_clust <- theta_updt(
              1:k_clust, m, theta, theta_se, theta_clust, pi_clust,
              junk_mixture, df, mu, sig_tmp, null_mixture, mu_null, sig_null
            )
            tmp_pi <- pi_updt(
              1:k_clust, m, pi_clust, theta, theta_se, theta_clust,
              junk_mixture, df, mu, sig_tmp, junk_prob, fix_junk_prob,
              null_mixture, mu_null, sig_null, null_prob, fix_null_prob
            )
            theta_clust <- tmp_clust
            pi_clust <- tmp_pi
            tmp_loglik <- loglik(
              m, pi_clust, theta, theta_se, theta_clust, junk_mixture, df, mu,
              sig_tmp, null_mixture, mu_null, sig_null
            )
            loglik_diff <- abs(tmp_loglik - tmp_loglik0)
            loglik_iter <- c(loglik_iter, tmp_loglik0)
            count <- count + 1
          }
          pi_tmp[[count2]] <- pi_clust
          theta_tmp[[count2]] <- theta_clust
          loglik_tmp[count2] <- tmp_loglik
          count2 <- count2 + 1
        }
        mx <- which.max(loglik_tmp)
        theta_clust <- theta_tmp[[mx]]
        pi_clust <- pi_tmp[[mx]]
        tmp_loglik <- loglik_tmp[mx]
        sig <- sigs[mx]
      } else {
        pi_clust <- clust_probs # initialise cluster class proportions
        theta_clust <- clust_means # initialise cluster centroids; these can be
        #the same as the std error is assumed to vary between the observations
        count <- 1 # start iteration count
        loglik_diff <- 1
        loglik_iter <- NULL # records all log likelihoods calculated during the run

        while (loglik_diff > tol & count < max_iter) {
          tmp_loglik0 <- loglik(m, pi_clust, theta, theta_se, theta_clust,
                                junk_mixture, df, mu, sig, null_mixture,
                                mu_null, sig_null)
          tmp_clust <- theta_updt(
            1:k_clust, m, theta, theta_se, theta_clust, pi_clust, junk_mixture,
            df, mu, sig, null_mixture, mu_null, sig_null
          )
          tmp_pi <- pi_updt(
            1:k_clust, m, pi_clust, theta, theta_se, theta_clust,
            junk_mixture, df, mu, sig, junk_prob, fix_junk_prob,
            null_mixture, mu_null, sig_null, null_prob, fix_null_prob
          )
          theta_clust <- tmp_clust
          pi_clust <- tmp_pi
          tmp_loglik <- loglik(m, pi_clust, theta, theta_se, theta_clust,
                               junk_mixture, df, mu, sig, null_mixture, mu_null,
                               sig_null)
          loglik_diff <- abs(tmp_loglik - tmp_loglik0)
          loglik_iter <- c(loglik_iter, tmp_loglik0)
          count <- count + 1
        }
      }

      clust_mn[[count0]] <- theta_clust
      clust_pi[[count0]] <- pi_clust
      log_like[[count0]] <- loglik_iter
      bic_clust[count0] <- if (!junk_null_aware_bic) {
        -2 * tmp_loglik + (2 * k_clust - 1) * log(m)
      } else {
        -2 * tmp_loglik + (2 * k_clust - 1 - sum(junk_mixture)
                           - sum(null_mixture)) * log(m)
      }
      count0 <- count0 + 1
    }
    bic_clust_mx[count_k] <- max(bic_clust[(count0 -
                                              (rand_num + 2)):(count0 - 1)])

    if ((count_k > stop_bic_iter) & (i >= min_clust_search)) { # as in section 3.3, repeat the algorithm with growing number of clusters until the BIC criterion starts increasing with the number of clusters
      tmp_cond <- TRUE
      for (j in 1:stop_bic_iter) {
        tmp_cond <- (tmp_cond & (bic_clust_mx[count_k - j + 1] >
                                   bic_clust_mx[count_k - j]))
      }
      if (tmp_cond) {
        break
      }
    }
    count_k <- count_k + 1
  }

  results <- clust_prob(bic_clust, clust_mn, clust_pi,
    theta = theta, theta_sd = theta_se, junk_mixture = junk_mixture, df = df,
    junk_mean = mu, junk_sd = sig, null_mixture = null_mixture,
    null_mean = mu_null, null_sd = sig_null, obs_names = obs_names,
    clust_size_prior = clust_size_prior, prior = bic_prior
  )
  results_all <- results[order(results$cluster), ]
  results_best <- best_clust(results)

  if (!is.null(cluster_membership)) {
    variant_clusters <- clust_inc_list(results,
                                       by_prob = cluster_membership$by_prob,
                                       bound = cluster_membership$bound)
  }
  if (!is.null(plot_results)) {
    if (plot_results[[1]] == "best") {
      plot_1 <- junk_clust_plot(results_best, sig = sig, mu = mu,
                                mu_null = mu_null)
      plot_2 <- two_stage_plot(results_best, bx, by, bxse, byse, obs_names)
    } else {
      tmp <- plot_results[[2]]
      tmp_res <- pr_clust(results, prob = tmp)
      plot_1 <- junk_clust_plot(tmp_res, sig = sig, mu = mu, mu_null = mu_null,
                                junk_mixture = junk_mixture,
                                null_mixture = null_mixture)
      plot_2 <- two_stage_plot(tmp_res, bx, by, bxse, byse, obs_names)
    }
    res <- list(
      results = list(all = results_all, best = results_best),
      cluster_membership = variant_clusters, plots = list(two_stage = plot_2,
                                                          split_plot = plot_1),
      log_likelihood = log_like, bic = bic_clust
    )
  } else {
    res <- list(
      results = list(all = results_all, best = results_best),
      cluster_membership = variant_clusters,
      log_likelihood = log_like, bic = bic_clust
    )
  }
  if (trait_search) { # for each cluster, check if its variant match specific biological traits that are less matched by variants outside the cluster. added for interpretability of results.
    trt_search <- pheno_search(variant_clusters = variant_clusters,
                               p_value = trait_pvalue, r2 = proxy_r2,
                               catalogue = catalogue, proxies = proxies,
                               build = build)
    res <- list(
      results = list(all = results_all, best = results_best),
      cluster_membership = variant_clusters, trait_search = trt_search,
      log_likelihood = log_like, bic = bic_clust
    )
  }
  return(res)
}
