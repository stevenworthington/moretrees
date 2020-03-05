# --------------------------------------------------------------------------------- #
# --------- Performs one step in VI optimization for Gaussian outcome  ------------ #
# --------------------------------------------------------------------------------- #

#'   \code{update_hyperparams_logistic} Performs variational updates for bernoulli outcomes.

update_vi_params_logistic <- function(X, groups, W,  xxT, wwT,
                                    y, n, K, G, m, # data
                                    prob, mu, Sigma, Sigma_inv, Sigma_det, tau_t, 
                                    delta, Omega, Omega_inv, Omega_det, 
                                    eta, g_eta, # variational params
                                    omega, rho, tau) { # hyperparams
  # Update sparse coefficients ------------------------------------------------------
  pred_g <- matrix(0, nrow = n, ncol = G)
  for (g in 2:G) {
    pred_g[, g] <- prob[g] *  X[ , groups[[g]], drop = F] %*% mu[[g]]
  }
  Wdelta <- W %*% delta
  # Update Sigma_g and tau_t_g
  for (g in 1:G) {
    xxT_g_eta <- xxT_g_eta_fun_ss(xxT[[g]], K = K[g], g_eta = g_eta)
    Sigma_inv[[g]] <- 2 * xxT_g_eta + diag(x = 1 / tau, nrow = K[g])
    Sigma[[g]] <- solve(Sigma_inv[[g]])
    Sigma_det[g] <- det(Sigma[[g]])
    tau_t[g] <- tau
    # update mu_g
    mu[[g]] <- Sigma[[g]] %*% crossprod( X[ , groups[[g]], drop = F],
                              y / 2 - 2 * g_eta * (Wdelta + rowSums(pred_g[, -g, drop = F])))
    # update prob_g (pi_g in manuscript)
    u <- 0.5 * crossprod(mu[[g]], Sigma_inv[[g]]) %*% mu[[g]] +
      0.5 * log(Sigma_det[g]) + log(rho / (1 - rho)) - 0.5 * K[g] * log(tau_t[g])
    prob[g] <- expit(u[1, 1])
    # update pred_g
    pred_g[, g] <- prob[g] *  X[ , groups[[g]], drop = F] %*% mu[[g]]
  }
  # Update non-sparse coefficients ---------------------------------------------------
  wwT_g_eta <- xxT_g_eta_fun_ss(wwT, m, g_eta)
  Omega_inv <- 2 * wwT_g_eta + diag(1 / omega, m)
  if (m != 0) {
    Omega <- solve(Omega_inv)
  }
  if (m == 1) {
    Omega_det <- Omega[1, 1]
  } else {
    Omega_det <- det(Omega)
  }
  # Update delta
  delta <- Omega %*% crossprod(W, y / 2 - 2 * g_eta * rowSums(pred_g))
  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
              Omega = Omega, Omega_inv = Omega_inv, Omega_det = Omega_det))
}