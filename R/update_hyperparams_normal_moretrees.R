# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

#'   \code{update_hyperparams_normal} Performs hyperparameter updates and computes 
#'   current value of ELBO in VI algorithm for gaussian outcomes

update_hyperparams_normal_moretrees <- function(X, groups, XtX, W, WtW, y,
                                outcomes_units, ancestors,
                                n, K, p, pL, m, # data
                                prob, mu, Sigma, Sigma_det, tau_t,
                                delta, Omega, Omega_det, # variational params
                                omega, sigma2, tau, rho, # hyperparameters
                                update_hyper = T) { 
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  # Sum of squared residuals
  xi <- mapply(FUN = function(prob, mu) prob * mu,
               prob = prob, mu = mu, SIMPLIFY = F)
  lp <- numeric(n) + 0
  for (v in 1:pL) {
    beta_v <- Reduce(`+`, xi[ancestors[[v]]])
    theta_v <- Reduce(`+`, delta[ancestors[[v]]])
    lp[outcomes_units[[v]]] <- X[outcomes_units[[v]], ] %*% beta_v + 
      W[outcomes_units[[v]], ] %*% theta_v
  }
  ssr <- sum( (y - lp) ^ 2 )
  # Expected sum of squared residuals
  prob_tr <- 0
  ssr_corr <- 0
  for (v in 1:p) {
    prob_tr <- prob_tr + prob[v] * trace_prod(Sigma[[v]], XtX[[v]]) +
      trace_prod(Omega[[v]], WtW[[v]])
    ssr_corr <- ssr_corr + 
      prob[v] * (1 - prob[v]) * crossprod(mu[[v]], XtX[[v]]) %*% mu[[v]]
  }
  expected_ssr <- as.numeric(prob_tr + ssr + ssr_corr)
  
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (v in 1:p) {
    expected_ss_gamma <- expected_ss_gamma + prob[v] * (sum(diag(Sigma[[v]])) + sum(mu[[v]] ^ 2))
  }
  expected_ss_gamma <- as.numeric(expected_ss_gamma + K * sum(tau_t * (1 - prob)))
  
  # Expected sum of squared thetas
  if (m == 0) {
    expected_ss_theta <- 0
  } else {
    expected_ss_theta <- 0
    for (v in 1:p) {
      expected_ss_theta <- expected_ss_theta +
        (sum(diag(Omega[[v]])) + sum(delta[[v]] ^ 2))
    }
  }
  
  # Update hyperparameters ---------------------------------------------------------
  if (update_hyper) {
    if (m != 0) {
      omega <- expected_ss_theta / (m * p)
    }
    tau <- expected_ss_gamma / (K * p)
    rho <- mean(prob)
    sigma2 <- expected_ssr / n
  }
  # Compute ELBO -------------------------------------------------------------------
  # See pg 5 of "variational inference for spike & slab model" document -
  # line numbers correspond to lines in equation
  line1 <- - 1 / (2 * sigma2) * expected_ssr - (n / 2) * log(2 * pi * sigma2)
  line2 <- - expected_ss_gamma / (2 * tau) - K * p * log(2 * pi * tau) / 2 + 
    log(rho ^ sum(prob)) + log((1 - rho) ^ (p - sum(prob)))
  line3 <- - expected_ss_theta / (2 * omega) - (m * p / 2) * log(2 * pi * omega)
  line4 <- (K * sum(prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line5 <- (K / 2) * sum(1 - prob) + 
    (K / 2) * sum(log(2 * pi * tau_t) * (1 - prob))
  line6 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line7 <- (m * p + sum(log(Omega_det)) + m * p * log(2 * pi)) / 2
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7
  # Return -------------------------------------------------------------------------
  return(list(ELBO = ELBO, omega = omega, sigma2 = sigma2, tau = tau, rho = rho))
}