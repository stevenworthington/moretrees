# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

#'   \code{update_hyperparams_normal} Performs hyperparameter updates and computes 
#'   current value of ELBO in VI algorithm for gaussian outcomes

update_hyperparams_normal <- function(X, groups, XtX, W, WtW, y, n, K, G, m, # data
                                prob, mu, Sigma, Sigma_det, tau_t,
                                delta, Omega, Omega_det, # variational params
                                omega, sigma2, tau, rho, # hyperparameters
                                update_hyper = T) { 
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  # Sum of squared residuals
  lp <- W %*% delta
  for (g in 1:G) {
    lp <- lp + prob[g] * X[ , groups[[g]], drop = F] %*% mu[[g]]
  }
  ssr <- sum( (y - lp) ^ 2 )
  # Expected sum of squared residuals
  prob_tr <- 0
  ssr_corr <- 0
  for (g in 1:G) {
    prob_tr <- prob_tr + prob[g] * trace_prod(Sigma[[g]], XtX[[g]])
    ssr_corr <- ssr_corr + 
      prob[g] * (1 - prob[g]) * t(mu[[g]]) %*% XtX[[g]] %*% mu[[g]]
  }
  expected_ssr <- as.numeric(trace_prod(WtW, Omega) + prob_tr + ssr + ssr_corr)
  
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (g in 1:G) {
    expected_ss_gamma <- expected_ss_gamma + prob[g] *
      (sum(diag(Sigma[[g]])) + sum(mu[[g]] ^ 2))
  }
  expected_ss_gamma <- as.numeric(expected_ss_gamma + sum(K * tau_t * (1 - prob)))
  
  # Expected sum of squared thetas
  if (m == 0) {
    expected_ss_theta <- 0
  } else {
    expected_ss_theta <- sum(diag(Omega)) + sum(delta ^ 2)
  }
  
  # Update hyperparameters ---------------------------------------------------------
  if (update_hyper) {
    if (m != 0) {
      omega <- expected_ss_theta / m
    }
    sigma2 <- expected_ssr / n
    tau <- expected_ss_gamma / sum(K)
    rho <- mean(prob)
  }
  # Compute ELBO -------------------------------------------------------------------
  # See pg 5 of "variational inference for spike & slab model" document -
  # line numbers correspond to lines in equation
  line1 <- - 1 / (2 * sigma2) * expected_ssr - (n / 2) * log(2 * pi * sigma2)
  line2 <- - expected_ss_gamma / (2 * tau) - 
    sum(K) * log(2 * pi * tau) / 2 + 
    log(rho ^ sum(prob)) +
    log((1 - rho) ^ (G - sum(prob)))
  line3 <- - expected_ss_theta / (2 * omega) - 
    (m / 2) * log(2 * pi * omega)
  line4 <- (sum(K * prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line5 <- (1 / 2) * sum(K * (1 - prob)) + 
    (1 / 2) * sum(K * log(2 * pi * tau_t) * (1 - prob))
  line6 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line7 <- (m + log(Omega_det) + m * log(2 * pi)) / 2
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7
  # Return -------------------------------------------------------------------------
  return(list(ELBO = ELBO, omega = omega, sigma2 = sigma2, tau = tau, rho = rho))
}