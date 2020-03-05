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

spike_and_slab_normal <- function(dsgn, initial_values,
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
  G <- length(dsgn$groups)
  n <- length(dsgn$y)
  K <- sapply(dsgn$groups, length)
  m <- ncol(dsgn$W)
  # Computing XtX and WtW so we don't have to do this repeatedly
  XtX <- lapply(dsgn$groups, function(cols) crossprod(dsgn$X[ , cols]))
  WtW <- crossprod(dsgn$W)
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
    Sigma_inv <- lapply(XtX, FUN = function(XtX, tau, sigma2) XtX / sigma2 + 
                          diag(x = 1 / tau, nrow = ncol(XtX)),
                        tau = hyperparams$tau, 
                        sigma2 = hyperparams$sigma2)  
    Sigma <- lapply(Sigma_inv, solve)
    Sigma_det <- sapply(Sigma, det)
    mu <- lapply(K, rnorm, mean = 0 , sd = vi_random_init$mu_sd)
    mu <- lapply(mu, matrix, ncol = 1)
    prob <- runif(G, 0, 1)
    tau_t <- rep(hyperparams$tau, G) # this should not be changed; tau_t = tau according to algorithm
    delta <- matrix(rnorm(m, sd = vi_random_init$delta_sd), ncol = 1)
    Omega_inv <- WtW / hyperparams$sigma2 + diag(x = 1 / hyperparams$omega, nrow = m)
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
    hyperparams <-  update_hyperparams_normal(X = dsgn$X, XtX = XtX, 
                                              groups = dsgn$groups,
                                              W = dsgn$W, WtW = WtW,
                                              y = dsgn$y, n = n,
                                              K = K, G = G, m = m,
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
    vi_params <- update_vi_params_normal(X = dsgn$X, groups = dsgn$groups,
                                         XtX = XtX, W = dsgn$W, WtW = WtW,
                                         y = dsgn$y, n = n, K = K, G = G, m = m,
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
      hyperparams <- update_hyperparams_normal(X = dsgn$X, groups = dsgn$groups,
                                               XtX = XtX, W = dsgn$W, WtW = WtW,
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
      hyperparams <- update_hyperparams_normal(X = dsgn$X, groups = dsgn$groups,
                                               XtX = XtX, W = dsgn$W, WtW = WtW,
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