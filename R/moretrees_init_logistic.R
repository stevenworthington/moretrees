# --------------------------------------------------------------------------------- #
# -------------------- moretrees initial values function -------------------------- #
# --------------------------------------------------------------------------------- #

#' Here's a brief description.
#'   \code{moretrees_init_rand} Randomly generates starting values for moretrees
#'   models. Not recommended if the model is converging slowly!
#' 
#' @export
#' @useDynLib moretrees
#' 
#' @section Model Description:
#' Describe MOReTreeS model and all parameters here.
#' 
#' @param dsgn Design list generated by moretrees_design_tree()
#' @param xxT Computed from exposure matrix X
#' @param wwT Computed from covariate matrix W
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param hyper_fixed Fixed values of hyperparameters to use if update_hyper = FALSE.
#' If family = "bernoulli", this should be a list including the following elements:
#' tau (prior variance for sparse node coefficients)
#' rho (prior node selection probability for sparse node coefficients)
#' omega (prior variance for non-sparse node coefficients)
#' If family = "gaussian", in addition to the above, the list should also include:
#' sigma2 (variance of residuals)
#' @param hyper_random_init If update_hyper = TRUE, this is a list containing the 
#' maximum values of the hyperparameters. Each hyperparameter will be initialised
#' uniformly at random between 0 and the maximum values given by the list elements
#' below. If multiple random restarts are being used, it is recommended
#' to use a large range for these initial values so that the parameter space
#' can be more effectively explored. The list contains the following elements:
#' tau_max (maxmimum of prior sparse node variance)
#' omega_max (maximum of prior non-sparse node variance)
#' sigma2_max (maximum of residual error variance--- for gaussian data only)
#' @param vi_random_init A list with parameters that determine the distributions from
#' which the initial VI parameters will be randomly chosen. All parameters will be randomly
#' selected from independent normal distributions with the standard deviations given by
#' the list elements below. If multiple random restarts are being used, it is recommended
#' to use large standard deviations for these initial values so that the parameter space
#' can be more effectively explored. The list contains the following elements:
#' mu_sd (standard deviation for posterior means of sparse node coefficients)
#' delta_sd (standard deviation for posterior means of non-sparse node coefficients)
#' xi_sd (standard deviation for auxilliary parameters xi--- for bernoulli data only)
#' @return A list containing starting values
#' @examples 
#' @family MOReTreeS functions

moretrees_init_logistic <- function(X, W, y, A,
                                outcomes_units,
                                outcomes_nodes,
                                ancestors,
                                xxT, wwT,
                                update_hyper,
                                hyper_fixed) {
  
  n <- length(y)
  m <- ncol(W)
  p <- length(unique(unlist(ancestors)))
  pL <- length(ancestors)
  K <- ncol(X)
  vi_params <- list()
  hyperparams <- list()
  
  # Get coefficient estimates from maximum likelihood ----------------------------------
  beta_ml <- matrix(0, nrow = p, ncol = K)
  theta_ml <- matrix(0, nrow = p, ncol = m)
  for (v in 1:p) {
    u <- outcomes_nodes[[v]]
    units <- unlist(outcomes_units[u])
    suppressWarnings(suppressMessages(
      if (m > 0){
        mod <- glm(y[units] == 1 ~ 0 + X[units,  , drop = F] 
                   + W[units,  , drop = F],
                   family = "binomial")
      } else {
        mod <- glm(y[units] == 1 ~ 0 + X[units,  , drop = F],
                   family = "binomial")
      }
    ))
    beta_ml[v, ] <- mod$coefficients[1:K]
    if (m > 0) {
      theta_ml[v, ] <- mod$coefficients[(K+1):(K + m)]
    }
  }
  # replace any NA vals with zero
  beta_ml[is.na(beta_ml)] <- 0
  theta_ml[is.na(theta_ml)] <- 0
  # transform to get initial values of mu and delta
  A_inv <- solve(A)
  mu <- A_inv %*% beta_ml
  delta <- A_inv %*% theta_ml
  vi_params$mu <- lapply(1:p, function(v, mu) matrix(mu[v, ], ncol = 1),
                                        mu = mu)
  vi_params$delta <- lapply(1:p, function(v, delta) matrix(delta[v, ], ncol = 1),
                                           delta = delta)
  
  # Set initial values for prob to be high (but not one) -------------------------------
  vi_params$prob <- rep(0.95, p)
  vi_params$a_rho <- 1 + sum(vi_params$prob) # need to initialise a_rho and b_rho using VI updates
  vi_params$b_rho <- 1 + p - sum(vi_params$prob) # so that terms cancel in ELBO.
  
  # Get starting values for eta --------------------------------------------------------
  # Use expected linear predictor squared 
  # (this is close to the real update for eta)
  xi <- mapply(FUN = function(prob, mu) prob * mu,
               prob = vi_params$prob, mu = vi_params$mu, SIMPLIFY = F)
  lp <- numeric(n) + 0
  for (v in 1:pL) {
    beta_v <- Reduce(`+`, xi[ancestors[[v]]])
    theta_v <- Reduce(`+`, vi_params$delta[ancestors[[v]]])
    lp[outcomes_units[[v]]] <- X[outcomes_units[[v]], , drop = F] %*% beta_v +
      W[outcomes_units[[v]], , drop = F ] %*% theta_v
  }
  hyperparams$eta <- abs(lp)
  hyperparams$g_eta <- gfun(hyperparams$eta)
  
  # Get hyperparameter starting values -------------------------------------------------
  if (update_hyper) {
    # If hyperparameters will be updated, initialise them
    hyperparams$omega <- var(as.numeric(delta))
    hyperparams$tau <- var(as.numeric(mu))
  } else {
    # Otherwise, use fixed values
    hyperparams$omega <- hyper_fixed$omega
    hyperparams$tau <- hyper_fixed$tau
  }
  if (m == 0) {
    hyperparams$omega <- 1
  }

  # Sigma and Omega initial values ---------------------------------------------------
  xxT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                      xxT = xxT, g_eta = hyperparams$g_eta, K = K)
  vi_params$Sigma_inv <- lapply(X = outcomes_nodes, 
                      FUN = function(outcomes, x, K, tau) 2 * Reduce(`+`, x[outcomes]) + 
                        diag(1 / tau, nrow = K),
                      x = xxT_g_eta,
                      K = K,
                      tau = hyperparams$tau)
  vi_params$Sigma <- lapply(vi_params$Sigma_inv, solve)
  vi_params$Sigma_det <- sapply(vi_params$Sigma, det)
  vi_params$tau_t <- rep(hyperparams$tau, p)
  if (m > 0) {
    wwT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                        xxT = wwT, g_eta = hyperparams$g_eta, K = m)
    vi_params$Omega_inv <- lapply(X = outcomes_nodes, 
                        FUN = function(outcomes, w, m, omega) 2 * Reduce(`+`, w[outcomes]) + 
                          diag(1 / omega, nrow = m),
                        w = wwT_g_eta,
                        m = m,
                        omega = hyperparams$omega)
    vi_params$Omega <- sapply(vi_params$Omega_inv, solve, simplify = F)
    vi_params$Omega_det <- sapply(vi_params$Omega, det, simplify = T)
  } else {
    vi_params$Omega <- rep(list(matrix(nrow = 0, ncol = 0)), p)
    vi_params$Omega_inv <- rep(list(matrix(nrow = 0, ncol = 0)), p)
    vi_params$Omega_det <- rep(1, p)
  }
  # Compute initial ELBO
  hyperparams <-  update_hyperparams_logistic_moretrees(X = X, 
                                                        W = W,
                                                        y = y, 
                                                        outcomes_units = outcomes_units,
                                                        ancestors = ancestors,
                                                        n = n, K = K, p = p, m = m,
                                                        prob = vi_params$prob, mu = vi_params$mu,
                                                        Sigma = vi_params$Sigma, Sigma_det = vi_params$Sigma_det,
                                                        tau_t = vi_params$tau_t, delta = vi_params$delta,
                                                        Omega = vi_params$Omega, Omega_det = vi_params$Omega_det,
                                                        eta = hyperparams$eta, g_eta = hyperparams$g_eta,
                                                        omega = hyperparams$omega, tau = hyperparams$tau,
                                                        a_rho = vi_params$a_rho, b_rho = vi_params$b_rho,
                                                        update_hyper = F)
  
  return(list(vi_params = vi_params, hyperparams = hyperparams))
}