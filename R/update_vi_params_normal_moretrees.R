# --------------------------------------------------------------------------------- #
# --------- Performs one step in VI optimization for Gaussian outcome  ------------ #
# --------------------------------------------------------------------------------- #

#'   \code{update_hyperparams_logistic} Performs variational updates for bernoulli outcomes.

update_vi_params_normal_moretrees <- function(X, groups, XtX, W, WtW, y, n, K, p, pL, m, # data
                                 ancestors, outcomes_nodes, outcomes_units,
                                 prob, mu, Sigma, Sigma_inv, Sigma_det, tau_t, 
                                 delta, Omega, Omega_inv, Omega_det, # variational params
                                 omega, sigma2, rho, tau, # hyperparams
                                 update_hyper_last) { 
  # Update sparse coefficients ------------------------------------------------------
  xi <- mapply(`*`, prob, mu, SIMPLIFY = F)
  Wtheta <- numeric(n) + 0
  for (u in 1:pL) {
    theta_u <- Reduce(`+`, delta[ancestors[[u]]])
    Wtheta[outcomes_units[[u]]] <- W[outcomes_units[[u]], ] %*% theta_u
  }
  for (v in 1:p) {
    # Update Sigma_v and tau_t_v only if hyperparameters were updated last round
    if (update_hyper_last) {
      Sigma_inv[[v]] <- XtX[[v]] / sigma2 + diag(nrow = K, x = 1 / tau)
      Sigma[[v]] <- solve(Sigma_inv[[v]])
      Sigma_det[v] <- det(Sigma[[v]])
      tau_t[v] <- tau
    }
    # update mu_g
    leaf_descendants <- outcomes_nodes[[v]]
    mu[[v]] <- matrix(0, nrow = K, ncol = 1)
    for (u in leaf_descendants) {
      anc_u_mv <- setdiff(ancestors[[u]], v)
      beta_u_mv <- Reduce(`+`, xi[anc_u_mv])
      units_u <- outcomes_units[[u]]
      mu[[v]] <- mu[[v]] + crossprod(X[units_u, ],
                y[units_u] - (X[units_u, ] %*% beta_u_mv + Wtheta[units_u]))
    }
    mu[[v]] <- (1 / sigma2) * Sigma[[v]] %*% mu[[v]]
    # update prob_g (pi_g in manuscript)
    u_v <- 0.5 * crossprod(mu[[v]], Sigma_inv[[v]]) %*% mu[[v]] +
      0.5 * log(Sigma_det[v]) + log(rho / (1 - rho)) - 0.5 * K * log(tau_t[v])
    prob[v] <- expit(u_v)
    # Update xi
    xi[[v]] <- prob[v] * mu[[v]]
  }
  # Update non-sparse coefficients ---------------------------------------------------
  # Update Omega only if hyperparameters were updated at last step
  if (m > 0) {
    Xbeta <- numeric(n) + 0
    for (u in 1:pL) {
      beta_u <- Reduce(`+`, xi[ancestors[[u]]])
      Xbeta[outcomes_units[[u]]] <- X[outcomes_units[[u]], ] %*% beta_u
    }
    for (v in 1:p) {
      # Update Omega_v
      if (update_hyper_last) {
        Omega_inv[[v]] <- WtW[[v]] / sigma2 + diag(nrow = m, x = 1 / omega)
        Omega[[v]] <- solve(Omega_inv[[v]])
        Omega_det[v] <- det(Omega[[v]])
      }
      leaf_descendants <- outcomes_nodes[[v]]
      # Update delta_v
      delta[[v]] <- delta[[v]] * 0
      for (u in leaf_descendants) {
        anc_u_mv <- setdiff(ancestors[[u]], v)
        theta_u_mv <- Reduce(`+`, delta[anc_u_mv])
        units_u <- outcomes_units[[u]]
        delta[[v]] <- delta[[v]] + crossprod(W[units_u, ],
             y[units_u] - (W[units_u, ] %*% theta_u_mv + Xbeta[units_u]))
      }
      delta[[v]] <- (1 / sigma2) * Omega[[v]] %*% delta[[v]]
    }
  }
  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
              Omega = Omega, Omega_inv = Omega_inv, Omega_det = Omega_det))
}
