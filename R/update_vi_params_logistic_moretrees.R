# --------------------------------------------------------------------------------- #
# --------- Performs one step in VI optimization for Gaussian outcome  ------------ #
# --------------------------------------------------------------------------------- #

#' \code{update_vi_logistic_moretrees} Performs variational updates for bernoulli outcomes.

update_vi_params_logistic_moretrees <- function(X, W, y, xxT, wwT,
                                      outcomes_nodes, outcomes_units,
                                      ancestors,
                                      n, p, pL, K, m, # data
                                      prob, mu, Sigma, Sigma_inv, Sigma_det, tau_t, 
                                      delta, Omega, Omega_inv, Omega_det, 
                                      a_rho, b_rho, # variational params
                                      eta, g_eta, 
                                      omega, tau) { # hyperparams
  
  # Update sparse coefficients ------------------------------------------------------
  xi <- mapply(`*`, prob, mu, SIMPLIFY = F)
  Wtheta <- numeric(n) + 0
  for (u in 1:pL) {
    theta_u <- Reduce(`+`, delta[ancestors[[u]]])
    Wtheta[outcomes_units[[u]]] <- W[outcomes_units[[u]], ] %*% theta_u
  }
  xxT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                        xxT = xxT, g_eta = g_eta, K = K)
  for (v in 1:p) {
    leaf_descendants <- outcomes_nodes[[v]]
    # Update Sigma_v and tau_t_v
    Sigma_inv[[v]] <- 2 * Reduce(`+`, xxT_g_eta[leaf_descendants]) + 
      diag(1 / tau, nrow = K)
    Sigma[[v]] <- solve(Sigma_inv[[v]])
    Sigma_det[v] <- det(Sigma[[v]])
    tau_t[v] <- tau
    # Update mu_v
    mu[[v]] <- matrix(0, nrow = K, ncol = 1)
    for (u in leaf_descendants) {
      anc_u_mv <- setdiff(ancestors[[u]], v)
      beta_u_mv <- Reduce(`+`, xi[anc_u_mv])
      units_u <- outcomes_units[[u]]
      mu[[v]] <- mu[[v]] + crossprod(X[units_u, , drop = FALSE],
        (y[units_u] / 2 - 2 * g_eta[units_u] * 
           (X[units_u, , drop = FALSE] %*% beta_u_mv + Wtheta[units_u]))
        )
    }
    mu[[v]] <- Sigma[[v]] %*% mu[[v]]
    # Update u_v
    u_v <- 0.5 * crossprod(mu[[v]], Sigma_inv[[v]]) %*% mu[[v]] +
      0.5 * log(Sigma_det[v]) + digamma(a_rho) - digamma(b_rho) -
      0.5 * K * log(tau_t[v])
    prob[v] <- expit(u_v)
    # Update xi
    xi[[v]] <- prob[v] * mu[[v]]
  }
  
  # Update rho ------------------------------------------------------------------------
  a_rho <- 1 + sum(prob) 
  b_rho <- 1 + p - sum(prob)
  
  # Update non-sparse coefficients ----------------------------------------------------
  if (m > 0) {
    Xbeta <- numeric(n) + 0
    for (u in 1:pL) {
      beta_u <- Reduce(`+`, xi[ancestors[[u]]])
      Xbeta[outcomes_units[[u]]] <- X[outcomes_units[[u]], ] %*% beta_u
    }
    wwT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                        xxT = wwT, g_eta = g_eta, K = m)
    for (v in 1:p) {
      # Update Omega_v
      leaf_descendants <- outcomes_nodes[[v]]
      Omega_inv[[v]] <- 2 * Reduce(`+`, wwT_g_eta[leaf_descendants]) +
        diag(1 / omega, nrow = m)
      Omega[[v]] <- solve(Omega_inv[[v]])
      Omega_det[v] <- det(Omega[[v]])
      # Update delta_v
      delta[[v]] <- delta[[v]] * 0
      for (u in leaf_descendants) {
        anc_u_mv <- setdiff(ancestors[[u]], v)
        units_u <- outcomes_units[[u]]
        theta_u_mv <- Reduce(`+`, delta[anc_u_mv])
        delta[[v]] <- delta[[v]] + crossprod(W[units_u, , drop = FALSE],
                 (y[units_u] / 2 - 2 * g_eta[units_u] *
                 (W[units_u, , drop = FALSE] %*% theta_u_mv + Xbeta[units_u])
                 ) )
      }
      delta[[v]] <- Omega[[v]] %*% delta[[v]]
    }
  }

  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
              Omega = Omega, Omega_inv = Omega_inv, Omega_det = Omega_det,
              a_rho = a_rho, b_rho = b_rho))
}