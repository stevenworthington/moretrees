#' Group spike and slab variable selection with Gaussian outcome
#' 
#' Here's a brief description.
#'   \code{spike_and_slab_logistic_moretrees} performs group variable selection via a spike
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
#' @param X Matrix of dimension n x sum(K), where n is the number of units, and
#' K[g] is the number of variables in group g.
#' @param W Matrix of non-sparse regression coefficients of dimension n x m
#' @param model "ss" or "moretrees". If "ss", regular spike and slab variable selection is implemented.
#' If "moretrees", the multi-outcome moretrees model is fitted.
#' @param groups If model = "ss", groups is a list of length G (number of groups), where groups[[g]] 
#' is an integer vector specifying the columns of X that belong to group g. NULL if model = "moretrees".
#' @param outcomes_units: If model = "moretrees", outcomes_units is a list of length equal to the number 
#' of unique outcomes. Each element of the list is an integer vector indicating which units (entries 
#' of y_reord, rows of X_reord) correspond to each outcomes. NULL if model = "ss".
#' @param outcomes_nodes: If model = "moretrees", outcomes_nodes is a list of length equal to the number 
#' of unique nodes. Each element of the list is an integer vector indicating which outcomes/leaves
#' are descendants of each node. NULL if model = "ss".
#' @param ancestors If model = "moretrees", ancestors ia  list of length equal to the number of unique 
#' outcomes. Each element of the list is an integer vector indicating which nodes on the tree (including 
#' leaves) are ancestors of the corresponding outcome. NULL if model = "ss".
#' @param tol Convergence tolerance for ELBO.
#' @param maxiter Maximum number of iterations of the VI algorithm.
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param update_hyper_freq How frequently to update hyperparameters. Default = every 10 iterations.
#' @return A list of variational parameters.
#' @examples
#' @family spike and slab functions

spike_and_slab_logistic_moretrees <- function(dsgn, initial_values,
                                              random_init,
                                              tol, max_iter,
                                              update_hyper, 
                                              update_hyper_freq,
                                              hyper_fixed,
                                              print_freq,
                                              hyper_random_init,
                                              vi_random_init) {
  if (is.null(dsgn$W)) {
    dsgn$W <- matrix(nrow = length(dsgn$y), ncol = 0)
  }
  # Prepare for running algorithm ---------------------------------------------------
  n <- length(dsgn$y)
  m <- ncol(dsgn$W)
  p <- length(unique(unlist(dsgn$ancestors)))
  pL <- length(dsgn$ancestors)
  K <- ncol(dsgn$X)
  if (K == 1) {
    xxT <- dsgn$X ^ 2
  } else {
    xxT <- rowOuterProds(dsgn$X)
  }
  if (m > 0) {
    if (m == 1) {
      wwT <- dsgn$W ^ 2
    } else {
      wwT <- rowOuterProds(dsgn$W)
    }
  } else {
    wwT <- NULL
  }
  # Initial hyperparameter values
  if (is.null(initial_values)) {
    if (random_init) {
      initial_values <- moretrees_init_rand(X = dsgn$X, W = dsgn$W, y = dsgn$y,
                                            outcomes_units = dsgn$outcomes_units,
                                            outcomes_nodes = dsgn$outcomes_nodes,
                                            ancestors = dsgn$ancestors,
                                            xxT = xxT, wwT = wwT,
                                            update_hyper = update_hyper,
                                            hyper_fixed = hyper_fixed,
                                            vi_random_init = vi_random_init,
                                            hyper_random_init = hyper_random_init)
    } else {
      initial_values <- moretrees_init_logistic(X = dsgn$X, W = dsgn$W, 
                                                y = dsgn$y, A = dsgn$A,
                                                outcomes_units = dsgn$outcomes_units,
                                                outcomes_nodes = dsgn$outcomes_nodes,
                                                ancestors = dsgn$ancestors,
                                                xxT = xxT, wwT = wwT,
                                                update_hyper = update_hyper,
                                                hyper_fixed = hyper_fixed)
    }
  } else {
    # In some cases, we may supply starting values for mu but not delta
    if (m > 0 & nrow(initial_values$vi_params$delta[[1]]) != m) {
      initial_values <- moretrees_init_W_logistic(X = dsgn$X, W = dsgn$W, 
                                                  y = dsgn$y, A = dsgn$A,
                                                  initial_values = initial_values,
                                                  outcomes_units = dsgn$outcomes_units,
                                                  outcomes_nodes = dsgn$outcomes_nodes,
                                                  ancestors = dsgn$ancestors,
                                                  xxT = xxT, wwT = wwT,
                                                  update_hyper = update_hyper,
                                                  hyper_fixed = hyper_fixed)
    } # otherwise, all initial values should already be supplied
  }
  vi_params <- initial_values$vi_params
  hyperparams <- initial_values$hyperparams

  # Initialise ELBO
  ELBO_track <- numeric(max_iter %/% update_hyper_freq + 1)
  ELBO_track[1] <- initial_values$hyperparams$ELBO
  ELBO_track2 <- numeric(max_iter + 1)
  ELBO_track2[1] <- initial_values$hyperparams$ELBO
  # Run algorithm -----------------------------------------------------------------
  i <- 0
  repeat {
    if (i >= max_iter) {
      cat(paste("Iteration", i, "complete.\n"))
      cat("\nWarning: Maximum number of iterations reached!\n")
      break
    }
    i <- i + 1
    if (i %% print_freq == 0 & i > 1) cat("Iteration", i, "; epsilon =",
                                  ELBO_track2[i] - ELBO_track2[i - 1] , "\n")
    update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
    vi_params <- update_vi_params_logistic_moretrees(y = dsgn$y, X = dsgn$X,
                                                     W = dsgn$W, xxT = xxT,
                                                     wwT = wwT,
                                                     outcomes_nodes = dsgn$outcomes_nodes,
                                                     outcomes_units = dsgn$outcomes_units,
                                                     ancestors = dsgn$ancestors,
                                                     n = n, p = p, pL = pL, K = K, m = m,
                                                     prob = vi_params$prob, 
                                                     mu = vi_params$mu, 
                                                     Sigma = vi_params$Sigma, 
                                                     Sigma_inv = vi_params$Sigma_inv, 
                                                     Sigma_det = vi_params$Sigma_det, 
                                                     tau_t = vi_params$tau_t,
                                                     a_rho = vi_params$a_rho,
                                                     b_rho = vi_params$b_rho,
                                                     delta = vi_params$delta, 
                                                     Omega = vi_params$Omega,
                                                     Omega_inv = vi_params$Omega_inv, 
                                                     Omega_det = vi_params$Omega_det,
                                                     eta = hyperparams$eta,
                                                     g_eta = hyperparams$g_eta,
                                                     omega = hyperparams$omega,
                                                     tau = hyperparams$tau)
    if (!update_hyper_i) {
      hyperparams <-  update_hyperparams_logistic_moretrees(X = dsgn$X, 
                                                            W = dsgn$W,
                                                            y = dsgn$y, 
                                                            outcomes_units = dsgn$outcomes_units,
                                                            ancestors = dsgn$ancestors,
                                                            n = n, K = K, p = p, m = m,
                                                            prob = vi_params$prob, mu = vi_params$mu,
                                                            Sigma = vi_params$Sigma, Sigma_det = vi_params$Sigma_det,
                                                            tau_t = vi_params$tau_t, delta = vi_params$delta,
                                                            Omega = vi_params$Omega, Omega_det = vi_params$Omega_det,
                                                            a_rho = vi_params$a_rho, b_rho = vi_params$b_rho,
                                                            eta = hyperparams$eta, g_eta = hyperparams$g_eta,
                                                            omega = hyperparams$omega, tau = hyperparams$tau,
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
      hyperparams <-   update_hyperparams_logistic_moretrees(X = dsgn$X, 
                                                             W = dsgn$W,
                                                             y = dsgn$y, 
                                                             outcomes_units = dsgn$outcomes_units,
                                                             ancestors = dsgn$ancestors,
                                                             n = n, K = K, p = p, m = m,
                                                             prob = vi_params$prob, mu = vi_params$mu,
                                                             Sigma = vi_params$Sigma, Sigma_det = vi_params$Sigma_det,
                                                             tau_t = vi_params$tau_t, delta = vi_params$delta,
                                                             Omega = vi_params$Omega, Omega_det = vi_params$Omega_det,
                                                             a_rho = vi_params$a_rho, b_rho = vi_params$b_rho,
                                                             eta = hyperparams$eta, g_eta = hyperparams$g_eta,
                                                             omega = hyperparams$omega, tau = hyperparams$tau,
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