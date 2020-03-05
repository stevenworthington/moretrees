#' Group spike and slab variable selection with Gaussian outcome
#' 
#' Here's a brief description.
#'   \code{spike_and_slab_logistic} performs group variable selection via a spike
#'   and slab prior for binary data.
#'   The posterior is approximated via variational inference.
#'   This function returns the parameters of the variational approximation.
#' 
#' All the details go here!
#' 
#' @section Model Description:
#'   Describe group spike and slab prior and all parameters here.
#' 
#' @param y Integer vector of length n containing outcomes; 1 = success, 0 = failure.
#' @param X matrix of dimension n x sum(K), where n is the number of units, and
#' K[g] is the number of variables in group g.
#' @param groups A list of length G (number of groups), where groups[[g]] is an integer
#' vector specifying the columns of X that belong to group g.
#' @param W matrix of non-sparse regression coefficients of dimension n x m
#' @param tol Convergence tolerance for ELBO.
#' @param maxiter Maximum number of iterations of the VI algorithm.
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param update_hyper_freq How frequently to update hyperparameters. Default = every 10 iterations.
#' @return A list of variational parameters.
#' @examples
#' @family spike and slab functions

spike_and_slab_logistic <- function(dsgn, initial_values,
                                    tol, max_iter,
                                    update_hyper, 
                                    update_hyper_freq,
                                    hyper_fixed,
                                    print_freq,
                                    hyper_random_init,
                                    vi_random_init) {
  if (is.null(dsgn$W)) {
    W <- matrix(nrow = length(dsgn$y), ncol = 0)
  }
  # Prepare for running algorithm ---------------------------------------------------
  G <- length(dsgn$groups)
  n <- length(dsgn$y)
  K <- sapply(dsgn$groups, length)
  m <- ncol(dsgn$W)
  xxT <- lapply(X = dsgn$groups, FUN = xxT_ss_fun, dat = dsgn$X)
  wwT <- rowOuterProds(dsgn$W)
  if (is.null(initial_values)) {
    # Initial hyperparameter values
    eta <- abs(rnorm(n, mean = 0, sd = vi_random_init$eta_sd))
    g_eta <- gfun(eta)
    if (update_hyper) {
      # If hyperparameters will be updated, randomly initialise them
      hyperparams <- list(omega = runif(1, 0, hyper_random_init$omega_max),
                          tau = runif(1, 0, hyper_random_init$tau_max),
                          rho = runif(1, 0, 1))
    } else {
      # Otherwise, use fixed values
      hyperparams <- hyper_fixed
    }
    
    if (m == 0) {
      hyperparams$omega <- 1
    }
    hyperparams$eta <- eta
    hyperparams$g_eta <- g_eta
    # Variational parameter initial values
    xxT_g_eta <- mapply(FUN = xxT_g_eta_fun_ss, xxT = xxT, K = K, MoreArgs = list(g_eta = g_eta),
                        SIMPLIFY = F)
    Sigma_inv <- mapply(FUN = function(x, K, tau) 2 * x + diag(x = 1 / tau, nrow = K),
                        x = xxT_g_eta, K = K, MoreArgs = list(tau = hyperparams$tau),
                        SIMPLIFY = F)
    Sigma <- lapply(Sigma_inv, solve)
    Sigma_det <- sapply(Sigma, det)
    mu <- lapply(K, rnorm, mean = 0 , sd = vi_random_init$mu_sd)
    mu <- lapply(mu, matrix, ncol = 1)
    prob <- runif(G, 0 , 1)
    tau_t <- rep(hyperparams$tau, G)
    delta <- matrix(rnorm(m, sd = vi_random_init$delta_sd), ncol = 1)
    wwT_g_eta <- xxT_g_eta_fun_ss(wwT, m, g_eta)
    Omega_inv <- 2 * wwT_g_eta + diag(1 / hyperparams$omega, m)
    if (m != 0) {
      Omega <- solve(Omega_inv)
      Omega_det <- det(Omega)
    } else {
      Omega <- matrix(nrow = 0, ncol = 0)
      Omega_det <- 1
    }
    # Put VI parameters in list
    vi_params <- list(mu = mu, prob = prob, Sigma = Sigma,
                      Sigma_inv = Sigma_inv, Sigma_det = Sigma_det,
                      tau_t = tau_t, delta = delta,
                      Omega = Omega, Omega_inv = Omega_inv,
                      Omega_det = Omega_det)
    # Compute initial ELBO
    hyperparams <-  update_hyperparams_logistic(X = dsgn$X, groups = dsgn$groups, W = dsgn$W,
                                                y = dsgn$y, n = n,
                                                K = K, G = G, m = m,
                                                prob = prob, mu = mu,
                                                Sigma = Sigma, Sigma_det = Sigma_det,
                                                tau_t = tau_t,
                                                delta = delta,
                                                Omega = Omega, Omega_det = Omega_det,
                                                eta = hyperparams$eta,
                                                g_eta = hyperparams$g_eta,
                                                omega = hyperparams$omega,
                                                tau = hyperparams$tau,
                                                rho = hyperparams$rho,
                                                update_hyper = F)
    initial_ELBO <- hyperparams$ELBO
  } else {
    # Otherwise, if initial values were supplied
    vi_params <- initial_values$vi_params
    hyperparams <- initial_values$hyperparams
    initial_ELBO <- initial_values$ELBO_track[length(initial_values$ELBO_track)]
  }
  ELBO_track <- numeric(max_iter %/% update_hyper_freq + 1)
  ELBO_track[1] <- initial_ELBO
  ELBO_track2 <- numeric(max_iter + 1)
  ELBO_track2[1] <- initial_ELBO
  # Run algorithm -----------------------------------------------------------------
  i <- 0
  repeat {
    if (i >= max_iter) {
      cat(paste("Iteration", i, "complete.\n"))
      cat("\nWarning: Maximum number of iterations reached!\n")
      break
    }
    i <- i + 1
    if (i %% print_freq == 0) cat("Iteration", i, "\n")
    update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
    vi_params <- update_vi_params_logistic(X = dsgn$X, groups = dsgn$groups, W = dsgn$W,
                                           y = dsgn$y, xxT = xxT, wwT = wwT,
                                           n = n, K = K, G = G, m = m,
                                           prob = vi_params$prob, 
                                           mu = vi_params$mu, 
                                           Sigma = vi_params$Sigma, 
                                           Sigma_inv = vi_params$Sigma_inv, 
                                           Sigma_det = vi_params$Sigma_det, 
                                           tau_t = vi_params$tau_t,
                                           delta = vi_params$delta, 
                                           Omega = vi_params$Omega,
                                           Omega_inv = vi_params$Omega_inv, 
                                           Omega_det = vi_params$Omega_det,
                                           eta = hyperparams$eta,
                                           g_eta = hyperparams$g_eta,
                                           omega = hyperparams$omega, 
                                           rho = hyperparams$rho, 
                                           tau = hyperparams$tau)
    if (!update_hyper_i) {
      hyperparams <- update_hyperparams_logistic(X = dsgn$X, groups = dsgn$groups, W = dsgn$W,
                                                 y = dsgn$y, n = n,
                                                 K = K, G = G, m = m,
                                                 prob = vi_params$prob, 
                                                 mu = vi_params$mu,
                                                 Sigma = vi_params$Sigma, 
                                                 Sigma_det = vi_params$Sigma_det,
                                                 tau_t = vi_params$tau_t,
                                                 delta = vi_params$delta,
                                                 Omega = vi_params$Omega, 
                                                 Omega_det = vi_params$Omega_det,
                                                 eta = hyperparams$eta,
                                                 g_eta = hyperparams$g_eta,
                                                 omega = hyperparams$omega,
                                                 tau = hyperparams$tau,
                                                 rho = hyperparams$rho,
                                                 update_hyper = F)
      ELBO_track2[i + 1] <- hyperparams$ELBO
      if (abs(ELBO_track2[i + 1] - ELBO_track2[i]) < tol) {
        # If we are not updating hyperparameters, we have converged
        if (!update_hyper) break
        # Otherwise, fill in results until next hyperparameter update
        i2 <- ceiling(i / update_hyper_freq) * update_hyper_freq
        if (i2 >= max_iter) {
          ELBO_track2[(i + 2):max_iter] <- hyperparams$ELBO
          i <- max_iter
          cat("Iteration", i, "complete.\n")
          cat("\nWarning: Maximum number of iterations reached!\n")
          break
        }
        ELBO_track2[(i + 2):(i2 + 1)] <- hyperparams$ELBO
        i <- i2
        update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
        update_hyper_im1 <- (i %% update_hyper_freq == 1)
      }
    }
    # Update hyperparameters
    if (update_hyper_i) {
      hyperparams <- update_hyperparams_logistic(X = dsgn$X, groups = dsgn$groups, W = dsgn$W,
                                                 y = dsgn$y, n = n,
                                                 K = K, G = G, m = m,
                                                 prob = vi_params$prob, 
                                                 mu = vi_params$mu,
                                                 Sigma = vi_params$Sigma, 
                                                 Sigma_det = vi_params$Sigma_det,
                                                 tau_t = vi_params$tau_t,
                                                 delta = vi_params$delta,
                                                 Omega = vi_params$Omega, 
                                                 Omega_det = vi_params$Omega_det,
                                                 eta = hyperparams$eta,
                                                 g_eta = hyperparams$g_eta,
                                                 omega = hyperparams$omega,
                                                 tau = hyperparams$tau,
                                                 rho = hyperparams$rho,
                                                 update_hyper = T)
      j <- i %/% update_hyper_freq
      ELBO_track[j + 1] <- hyperparams$ELBO
      ELBO_track2[i + 1] <- hyperparams$ELBO
      if (abs(ELBO_track[j + 1] - ELBO_track[j]) < tol) break
    }
  }
  return(list(vi_params = vi_params, hyperparams = hyperparams,
              ELBO_track = ELBO_track2[1:(i + 1)]))
}