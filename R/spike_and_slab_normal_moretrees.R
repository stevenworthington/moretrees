#' Group spike and slab variable selection with Gaussian outcome
#' 
#' Here's a brief description.
#'   \code{spike_and_slab_normal} performs group variable selection via a spike
#'   and slab prior for continuous, normally distributed data.
#'   The posterior is approximated via variational inference.
#'   This function returns the parameters of the variational approximation.
#' 
#' All the details go here!
#' 
#' @section Model Description:
#'   Describe group spike and slab prior and all parameters here.
#'   
#' @param y Numeric vector of length n of outcome data
#' @param X Matrix of dimension n x sum(K), where n is the number of units, and
#' K[g] is the number of variables in group g.
#' @param groups A list of length G (number of groups), where groups[[g]] is an integer
#' vector specifying the columns of X that belong to group g.
#' @param W Matrix of data with non-sparse regression coefficients of dimension n x m
#' @param tol Convergence tolerance for ELBO.
#' @param max_iter Maximum number of iterations of the VI algorithm.
#' @param nrestarts Number of random re-starts of the VI algorithm. The result that 
#' gives the highest ELBO will be returned. It is recommended to choose nrestarts > 1.
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param update_hyper_freq How frequently to update hyperparameters. 
#' Default = every 50 iterations.
#' @param print_freq How often to print out iteration number. 
#' @return A list of variational parameters.
#' @examples
#' @family spike and slab functions

spike_and_slab_normal_moretrees <- function(dsgn,
                                            initial_values,
                                            tol, max_iter,
                                            update_hyper, 
                                            update_hyper_freq,
                                            hyper_fixed,
                                            print_freq,
                                            hyper_random_init,
                                            vi_random_init) {
  if (is.null(dsgn$W)) {
    W <-  matrix(nrow = length(dsgn$y), ncol = 0)
  }
  # Prepare for running algorithm ---------------------------------------------------
  n <- length(dsgn$y)
  m <- ncol(dsgn$W)
  p <- length(unique(unlist(dsgn$ancestors)))
  pL <- length(dsgn$ancestors)
  K <- ncol(dsgn$X)
  # Computing XtX and WtW so we don't have to do this repeatedly
  XtX <- lapply(dsgn$outcomes_units, function(units) crossprod(dsgn$X[units , ]))
  XtX <- lapply(dsgn$outcomes_nodes, function(outcomes) Reduce(`+`, XtX[outcomes]))
  WtW <- lapply(dsgn$outcomes_units, function(units) crossprod(dsgn$W[units , ]))
  WtW <- lapply(dsgn$outcomes_nodes, function(outcomes) Reduce(`+`, WtW[outcomes]))
  if (is.null(initial_values)) {
    # Initial hyperparameter values
    if (update_hyper) {
      # If hyperparameters will be updated, randomly initialise them
      hyperparams <- list(omega = runif(1, 0, hyper_random_init$omega_max),
                          tau = runif(1, 0, hyper_random_init$tau_max),
                          rho = runif(1, 0, 1),
                          sigma2 = runif(1, 0, hyper_random_init$sigma2_max))
    } else {
      # Otherwise, use fixed values
      hyperparams <- hyper_fixed
    }
    if (m == 0) {
      hyperparams$omega <- 1
    }
    # Variational parameter initial values
    Sigma_inv <- lapply(XtX, FUN = function(XtX, tau, sigma2, K) XtX / sigma2 + 
                          diag(x = 1 / tau, nrow = K), tau = hyperparams$tau,
                        sigma2 = hyperparams$sigma2, K = K)  
    Sigma <- lapply(Sigma_inv, solve)
    Sigma_det <- sapply(Sigma, det)
    mu <- lapply(X = 1:p, FUN = function(i) matrix(rnorm(K), ncol = 1))
    prob <- runif(p, 0, 1)
    tau_t <- rep(hyperparams$tau, p) # this should not be changed; tau_t = tau according to algorithm
    delta <- lapply(X = 1:p, FUN = function(i) matrix(rnorm(m), ncol = 1))
    if (m > 0) {
      Omega_inv <- lapply(WtW, FUN = function(WtW, omega, sigma2, m) WtW / sigma2 + 
                            diag(x = 1 / omega, nrow = m), omega = hyperparams$omega,
                          sigma2 = hyperparams$sigma2, m = m)   
      Omega <- lapply(Omega_inv, solve)
      Omega_det <- sapply(Omega, det)
    } else {
      Omega <- rep(list(matrix(nrow = 0, ncol = 0)), p)
      Omega_inv <- rep(list(matrix(nrow = 0, ncol = 0)), p)
      Omega_det <- rep(1, p)
    }
    
    # Put VI parameters in list
    vi_params <- list(mu = mu, prob = prob, Sigma = Sigma,
                      Sigma_inv = Sigma_inv, Sigma_det = Sigma_det,
                      tau_t = tau_t, delta = delta,
                      Omega = Omega, Omega_inv = Omega_inv,
                      Omega_det = Omega_det)
    # Compute initial ELBO
    hyperparams <-  update_hyperparams_normal_moretrees(X = dsgn$X, XtX = XtX, 
                                                        groups = dsgn$groups,
                                                        W = dsgn$W, WtW = WtW,
                                                        y = dsgn$y, 
                                                        outcomes_units = dsgn$outcomes_units,
                                                        ancestors = dsgn$ancestors,
                                                        n = n, K = K, p = p, pL = pL, m = m,
                                                        prob = prob, mu = mu,
                                                        Sigma = Sigma, Sigma_det = Sigma_det,
                                                        tau_t = tau_t,
                                                        delta = delta,
                                                        Omega = Omega, Omega_det = Omega_det,
                                                        omega = hyperparams$omega,
                                                        sigma2 = hyperparams$sigma2,
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
      cat("Iteration", i, "complete.\n")
      cat("\nWarning: Maximum number of iterations reached!\n")
      break
    }
    i <- i + 1
    update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
    update_hyper_im1 <- (i %% update_hyper_freq == 1) & update_hyper
    vi_params <- update_vi_params_normal_moretrees(X = dsgn$X, groups = dsgn$groups,
                                                   XtX = XtX, W = dsgn$W, WtW = WtW,
                                                   ancestors = dsgn$ancestors,
                                                   outcomes_nodes = dsgn$outcomes_nodes,
                                                   outcomes_units = dsgn$outcomes_units,
                                                   y = dsgn$y, 
                                                   n = n, K = K, p = p, pL = pL, m = m,
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
                                                   omega = hyperparams$omega, 
                                                   sigma2 = hyperparams$sigma2,  
                                                   rho = hyperparams$rho, 
                                                   tau = hyperparams$tau,
                                                   update_hyper_last = update_hyper_im1)
    if (!update_hyper_i) {
      hyperparams <- update_hyperparams_normal_moretrees(X = dsgn$X, groups = dsgn$groups,
                                                         XtX = XtX, W = dsgn$W, WtW = WtW,
                                                         y = dsgn$y, 
                                                         outcomes_units = dsgn$outcomes_units,
                                                         ancestors = dsgn$ancestors,
                                                         n = n, K = K, p = p, pL = pL, m = m,
                                                         prob = vi_params$prob, 
                                                         mu = vi_params$mu,
                                                         Sigma = vi_params$Sigma, 
                                                         Sigma_det = vi_params$Sigma_det,
                                                         tau_t = vi_params$tau_t,
                                                         delta = vi_params$delta,
                                                         Omega = vi_params$Omega, 
                                                         Omega_det = vi_params$Omega_det,
                                                         omega = hyperparams$omega,
                                                         sigma2 = hyperparams$sigma2,
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
      hyperparams <- update_hyperparams_normal_moretrees(X = dsgn$X, groups = dsgn$groups,
                                                         XtX = XtX, W = dsgn$W, WtW = WtW,
                                                         y = dsgn$y,
                                                         outcomes_units = dsgn$outcomes_units,
                                                         ancestors = dsgn$ancestors,
                                                         n = n, K = K, p = p, pL = pL, m = m,
                                                         prob = vi_params$prob, 
                                                         mu = vi_params$mu,
                                                         Sigma = vi_params$Sigma, 
                                                         Sigma_det = vi_params$Sigma_det,
                                                         tau_t = vi_params$tau_t,
                                                         delta = vi_params$delta,
                                                         Omega = vi_params$Omega, 
                                                         Omega_det = vi_params$Omega_det,
                                                         omega = hyperparams$omega,
                                                         sigma2 = hyperparams$sigma2,
                                                         tau = hyperparams$tau,
                                                         rho = hyperparams$rho,
                                                         update_hyper = T)
      j <- i %/% update_hyper_freq
      ELBO_track[j + 1] <- hyperparams$ELBO
      ELBO_track2[i + 1] <- hyperparams$ELBO
      if (abs(ELBO_track[j + 1] - ELBO_track[j]) < tol) break
    } 
    if (i %% print_freq == 0) cat("Iteration", i, "complete.\n")
  }
  return(list(vi_params = vi_params, hyperparams = hyperparams,
              ELBO_track = ELBO_track2[1:(i + 1)]))
}