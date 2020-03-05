# --------------------------------------------------------------------------------- #
# ----------- spike & slab variable selection wrapper function -------------------- #
# --------------------------------------------------------------------------------- #

#' Group spike and slab variable selection with bernoulli or gaussian outcome
#' 
#' Here's a brief description.
#'   \code{spike_and_slab_normal} performs group variable selection via a spike
#'   and slab prior. The posterior is approximated via variational inference.
#'   This function returns coefficient estimates and 95% credible intervals.
#' 
#' All the details go here!
#' 
#' @section Model Description:
#'   Describe group spike and slab prior and all parameters here.
#' 
#' @export
#' @param y Vector of length n containing outcomes data.
#' If family = "bernoulli", y must be an integer vector where 1 = success, 0 = failure.
#' If family = "gaussian", y must be a numeric vector containing continuous data.
#' @param X Matrix of dimension n x sum(K), where n is the number of units, and
#' K[g] is the number of variables in group g.
#' @param groups A list of length G (number of groups), where groups[[g]] is an integer
#' vector specifying the columns of X that belong to group g.
#' @param W Matrix of data with non-sparse regression coefficients of dimension n x m
#' @param family A string specifying the distribution of the outcomes: 
#' either "bernoulli" (for classification) or "gaussian" (for regression)
#' @param ci_level A number between 0 and 1 giving the desired credible interval.
#' For example, ci_level = 0.95 (the default) returns a 95\% credible interval.
#' @param tol Convergence tolerance for ELBO. Default = 1E-8.
#' @param max_iter Maximum number of iterations of the VI algorithm.
#' @param nrestarts Number of random re-starts of the VI algorithm. The result that 
#' gives the highest ELBO will be returned. It is recommended to choose nrestarts > 1.
#' The default is 3.
#' @param keep_restarts If TRUE, the results from all random restarts will be returned.
#' If FALSE, only the restart with the highest ELBO is returned.
#' @param parallel If TRUE, the random restarts will be run in parallel. It is recommended
#' to first set the number of cores using doParallel::registerDoParallel(). Otherwise,
#' the default number of cores specified by the doParallel package will be used.
#' @param log_restarts If TRUE, progress of each random restart will be logged to a text
#' file in log_dir.
#' @param log_dir Directory for logging progress of random restarts.
#' Default is the working directory.
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param update_hyper_freq How frequently to update hyperparameters. 
#' Default = every 50 iterations.
#' @param hyper_fixed Fixed values of hyperparameters to use if update_hyper = FALSE.
#' If family = "bernoulli", this should be a list including the following elements:
#' tau (prior variance for sparse coefficients)
#' rho (prior selection probability for sparse coefficients)
#' omega (prior variance for non-sparse coefficients)
#' If family = "gaussian", in addition to the above, the list should also include:
#' sigma2 (variance of residuals)
#' @param hyper_random_init If update_hyper = TRUE, this is a list containing the 
#' maximum values of the hyperparameters. Each hyperparameter will be initialised
#' uniformly at random between 0 and the maximum values given by the list elements
#' below. If multiple random restarts are being used, it is recommended
#' to use a large range for these initial values so that the parameter space
#' can be more effectively explored. The list contains the following elements:
#' tau_max (maxmimum of prior sparse coefficient variance)
#' omega_max (maximum of prior non-sparse coefficient variance)
#' sigma2_max (maximum of residual error variance--- for gaussian data only)
#' @param vi_random_init A list with parameters that determine the distributions from
#' which the initial VI parameters will be randomly chosen. All parameters will be randomly
#' selected from independent normal distributions with the standard deviations given by
#' the list elements below. If multiple random restarts are being used, it is recommended
#' to use large standard deviations for these initial values so that the parameter space
#' can be more effectively explored. The list contains the following elements:
#' mu_sd (standard deviation for posterior means of sparse coefficients)
#' delta_sd (standard deviation for posterior means of non-sparse coefficients)
#' xi_sd (standard deviation for auxilliary parameters xi--- for bernoulli data only)
#' @param print_freq How often to print out iteration number. 
#' @return A list containing the following elements:
#' 1. estimated coefficients and credible intervals; 
#' 2. outputs from variational inference algorithm
#' @examples
#' @family spike and slab functions

spike_and_slab <- function(y, X, groups, W = NULL, 
                           initial_values = NULL,
                           family = "bernoulli",
                           ci_level = 0.95,
                           tol = 1E-4, max_iter = 1E5,
                           nrestarts = 3,
                           keep_restarts = nrestarts > 1,
                           parallel = nrestarts > 1,
                           log_restarts = nrestarts > 1,
                           log_dir = getwd(),
                           update_hyper = T, update_hyper_freq = 10,
                           print_freq = update_hyper_freq,
                           hyper_fixed = NULL, 
                           hyper_random_init = list(omega_max = 100,
                                                    tau_max = 100,
                                                    sigma2_max = 100),
                           vi_random_init = list(eta_sd = 10,
                                                 mu_sd = 10,
                                                 delta_sd = 10)) {
  # Get correct function
  if (!(family %in% c("bernoulli", "gaussian"))) {
    stop("family must be a string (\"bernoulli\" or \"gaussian\")")
  }
  if (family == "bernoulli") {
    ss_fun <- spike_and_slab_logistic
    y[y == 0] <- -1
  }
  if (family == "gaussian") ss_fun <- spike_and_slab_normal
  if (!update_hyper & is.null(hyper_fixed)) {
    stop("Must supply fixed hyperparameter values if update_hyper = FALSE")
  }
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!is.null(W) & !(is.matrix(W))) stop("If W is not NULL, must be a matrix")
  
  # Setting up parallelization
  if (parallel) {
    `%doRestarts%` <- foreach::`%dopar%`
  } else {
    `%doRestarts%` <- foreach::`%do%`
  }
  
  # Run algorithm
  # mod_restarts <- foreach::foreach(i = 1:nrestarts) %doRestarts% {
  #   if (log_restarts) {
  #     sink(file = paste0(log_dir, "restart_", i, "_log.txt"))
  #     cat("Initialising random restart", i, "...\n\n")
  #   }
    mod_restarts <- list(ss_fun(dsgn = list(y = y, X = X, groups = groups, W = W),
             initial_values = initial_values,
             tol = tol, max_iter = max_iter,
             update_hyper = update_hyper, 
             update_hyper_freq = update_hyper_freq,
             hyper_fixed = hyper_fixed,
             print_freq = print_freq,
             hyper_random_init = hyper_random_init,
             vi_random_init = vi_random_init))
  #   if (log_restarts) {
  #     cat("\nRestart", i, "complete.")
  #     sink()
  #   }
  #   mod
  # }
  
  # Select random restart that gave the highest ELBO
  ELBO_restarts <- sapply(mod_restarts, FUN = function(mod) mod$ELBO_track[length(mod$ELBO_track)])
  best_restart <- which.max(ELBO_restarts)
  mod <- mod_restarts[[best_restart]]
  if (keep_restarts) {
    mod_restarts <- mod_restarts[- best_restart]
  } else {
    rm(mod_restarts)
    mod_restarts <- NULL
  }

  # Compute estimates and credible intervals
  beta_est <- compute_betas(mod, ci_level)
  theta_est <- compute_thetas(mod, ci_level)
  
  # Return
  return(list(sparse_est = beta_est, nonsparse_est = theta_est,
              mod = mod, mod_restarts = mod_restarts))
}
