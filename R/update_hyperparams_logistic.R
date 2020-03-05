# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

#'   \code{update_hyperparams_logistic} Performs hyperparameter updates and computes 
#'   current value of ELBO in VI algorithm for bernoulli outcomes.

update_hyperparams_logistic <- function(X, groups, W, y, n, K, G, m, # data
                                prob, mu, Sigma, Sigma_det, tau_t,
                                delta, Omega, Omega_det, 
                                eta, g_eta, # variational params
                                omega, tau, rho, # hyperparameters
                                update_hyper = T) {
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  # Expected linear predictor
  lp <- W %*% delta
  for (g in 1:G) {
    lp <- lp + prob[g] *  X[ , groups[[g]], drop = F] %*% mu[[g]]
  }
  lp <- as.numeric(lp)
  # Expected linear predictor squared
  lp2 <- emulator::quad.tdiag(Omega, W)
  for (g in 1:G) {
    Sigma_g <- Sigma[[g]] + (1 - prob[g]) * tcrossprod(mu[[g]])
    lp2 <- lp2 + prob[g] * quadFormByRow(Sigma_g, X[, groups[[g]], drop = F])
  }
  lp2 <- lp2 + lp ^ 2

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
    tau <- expected_ss_gamma / sum(K)
    rho <- mean(prob)
  }
  # Update eta  --------------------------------------------------------------------
  eta <- sqrt(lp2)
  g_eta <- gfun(eta)
  
  # Compute ELBO -------------------------------------------------------------------
  # See pg 13 of "variational inference for spike & slab model" document -
  # line numbers correspond to lines in equation
  line1 <- (1 / 2) * t(y) %*% lp + 
    sum(logexpit(eta)) - sum(eta) / 2 + g_eta %*% (eta ^ 2)
  line2 <- - g_eta %*% lp2
  line3 <- - expected_ss_gamma / (2 * tau) - 
    sum(K) * log(2 * pi * tau) / 2 + 
    log(rho ^ sum(prob)) +
    log((1 - rho) ^ (G - sum(prob)))
  line4 <- - expected_ss_theta / (2 * omega) - 
    (m / 2) * log(2 * pi * omega)
  line5 <- (sum(K * prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line6 <- (1 / 2) * sum(K * (1 - prob)) + 
    (1 / 2) * sum(K * log(2 * pi * tau_t) * (1 - prob))
  line7 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line8 <- (m + log(Omega_det) + m * log(2 * pi)) / 2
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8
  # Return -------------------------------------------------------------------------
  return(list(ELBO = as.numeric(ELBO), omega = omega, tau = tau, rho = rho,
              eta = eta, g_eta = g_eta))
}