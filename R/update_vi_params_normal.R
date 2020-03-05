# --------------------------------------------------------------------------------- #
# --------- Performs one step in VI optimization for Gaussian outcome  ------------ #
# --------------------------------------------------------------------------------- #

#'   \code{update_hyperparams_logistic} Performs variational updates for bernoulli outcomes.

update_vi_params_normal <- function(X, groups, XtX, W, WtW, y, n, K, G, m, # data
                                 prob, mu, Sigma, Sigma_inv, Sigma_det, tau_t, 
                                 delta, Omega, Omega_inv, Omega_det, # variational params
                                 omega, sigma2, rho, tau, # hyperparams
                                 update_hyper_last) { 
  # Update sparse coefficients ------------------------------------------------------
  pred_g <- matrix(0, nrow = n, ncol = G)
  for (g in 2:G) {
    pred_g[, g] <- prob[g] *  X[ , groups[[g]], drop = F] %*% mu[[g]]
  }
  Wdelta <- W %*% delta
  for (g in 1:G) {
    # Update Sigma_g and tau_t_g only if hyperparameters were updated last round
    if (update_hyper_last) {
      Sigma_inv[[g]] <- XtX[[g]] / sigma2 + diag(nrow = K[g], x = 1 / tau)
      Sigma[[g]] <- solve(Sigma_inv[[g]])
      Sigma_det[g] <- det(Sigma[[g]])
      tau_t[g] <- tau
    }
    # update mu_g
    mu[[g]] <- (1 / sigma2) * Sigma[[g]] %*%
      crossprod( X[ , groups[[g]], drop = F], y - Wdelta - rowSums(pred_g[, -g, drop = F]))
    # update prob_g (pi_g in manuscript)
    u <- 0.5 * crossprod(mu[[g]], Sigma_inv[[g]]) %*% mu[[g]] +
      0.5 * log(Sigma_det[g]) + log(rho / (1 - rho)) - 0.5 * K[g] * log(tau_t[g])
    prob[g] <- expit(u[1, 1])
    # update pred_g
    pred_g[, g] <- prob[g] *  X[ , groups[[g]], drop = F] %*% mu[[g]]
  }
  # Update non-sparse coefficients ---------------------------------------------------
  # Update Omega only if hyperparameters were updated at last step
  if (update_hyper_last) {
    Omega_inv <- WtW / sigma2 + diag(nrow = m, x = 1 / omega)
    if (m != 0) {
      Omega <- solve(Omega_inv)
    }
    if (m == 1) {
      Omega_det <- Omega[1, 1]
    } else {
      Omega_det <- det(Omega)
    }
  }
  # Update delta
  delta <- (1 / sigma2) * Omega %*% t(W) %*% (y - rowSums(pred_g))
  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
              Omega = Omega, Omega_inv = Omega_inv, Omega_det = Omega_det))
}
