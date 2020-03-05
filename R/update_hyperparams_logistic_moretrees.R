# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

#'   \code{update_hyperparams_logistic_moretrees} Performs hyperparameter updates and computes 
#'   current value of ELBO in VI algorithm for bernoulli outcomes.

update_hyperparams_logistic_moretrees <- function(X, W, y, 
                                        outcomes_units,
                                        ancestors,
                                        n, K, p, m, # data
                                        prob, mu, Sigma, Sigma_det, tau_t,
                                        delta, Omega, Omega_det, 
                                        eta, g_eta,
                                        a_rho, b_rho, # variational params
                                        omega, tau, # hyperparameters
                                        model = "ss",
                                        update_hyper = T) {
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  # Expected linear predictor
  xi <- mapply(FUN = function(prob, mu) prob * mu,
                 prob = prob, mu = mu, SIMPLIFY = F)
  lp <- numeric(n) + 0
  for (v in 1:length(ancestors)) {
    beta_v <- Reduce(`+`, xi[ancestors[[v]]])
    theta_v <- Reduce(`+`, delta[ancestors[[v]]])
    lp[outcomes_units[[v]]] <- X[outcomes_units[[v]], ] %*% beta_v + 
      W[outcomes_units[[v]], ] %*% theta_v
  }
  # Expected linear predictor squared
  Sigma_u <- mapply(FUN = function(prob, Sigma, mu) prob * 
                    (Sigma + (1 - prob) * tcrossprod(mu)),
                    prob = prob, Sigma = Sigma, mu = mu,
                    SIMPLIFY = F)
  lp2 <- numeric(n) + 0
  for (v in 1:length(ancestors)) {
    Sigma_v <- Reduce(`+`, Sigma_u[ancestors[[v]]])
    Omega_v <- Reduce(`+`, Omega[ancestors[[v]]])
    lp2[outcomes_units[[v]]] <- quadFormByRow(Omega_v, W[outcomes_units[[v]], , drop = F]) +
      quadFormByRow(Sigma_v, X[outcomes_units[[v]], , drop = F])
  }
  lp2 <- lp2 + lp ^ 2
  
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (v in 1:p) {
    expected_ss_gamma <- expected_ss_gamma + prob[v] *
      (sum(diag(Sigma[[v]])) + sum(mu[[v]] ^ 2))
  }
  expected_ss_gamma <- as.numeric(expected_ss_gamma + sum(K * tau_t * (1 - prob)))
  
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
  }
  # Update eta  --------------------------------------------------------------------
  eta <- as.numeric(sqrt(lp2))
  g_eta <- as.numeric(gfun(eta))
  
  # Compute ELBO -------------------------------------------------------------------
  # See pg 13 of "variational inference for spike & slab model" document -
  # line numbers correspond to lines in equation
  line1 <- (1 / 2) * crossprod(y, lp) + sum(logexpit(eta)) - sum(eta) / 2 + g_eta %*% (eta ^ 2)
  line2 <- - g_eta %*% lp2
  line3 <- 0
  line4 <- - expected_ss_gamma / (2 * tau) - p * K * log(2 * pi * tau) / 2 
  line5 <- - expected_ss_theta / (2 * omega) - (m * p / 2) * log(2 * pi * omega)
  line6 <- (K * sum(prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line7 <- (K / 2) * sum(1 - prob) + 
    (K / 2) * sum(log(2 * pi * tau_t) * (1 - prob))
  line8 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line9 <- (m * p + sum(log(Omega_det)) + m * p * log(2 * pi)) / 2
  line10 <- lbeta(a_rho, b_rho)
  
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8 + line9 + line10
  # Return -------------------------------------------------------------------------
  return(list(ELBO = as.numeric(ELBO), omega = omega, tau = tau,
              eta = eta, g_eta = g_eta))
}