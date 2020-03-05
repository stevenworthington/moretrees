# --------------------------------------------------------------------------------- #
# ------------------------- moretrees wrapper function ---------------------------- #
# --------------------------------------------------------------------------------- #

#' Here's a brief description.
#'   \code{moretrees} Fits Multi-Outcome Regression with Tree-structured Shrinkage
#'   (MOReTreeS) model to normally-distributed or binary outcome data.
#'   The posterior is approximated via variational inference.
#'   Returns coefficient estimates and 95% credible intervals.
#' 
#' All the details go here!
#' 
#' @export
#' @useDynLib moretrees
#' 
#' @section Model Description:
#' Describe MOReTreeS model and all parameters here.
#' 
#' @param y Vector of length n containing outcomes data.
#' If family = "bernoulli", y must be an integer vector where 1 = success, 0 = failure.
#' If family = "gaussian", y must be a numeric vector containing continuous data.
#' @param X An n x K matrix of exposure data, where K is the dimension of the exposure.
#' Grouping of the outcomes will be based on their relationships with the variables in X.
#' @param W Matrix of covariates of dimension n x m.
#' Coefficients for these variables do not affect grouping of the outcomes.
#' @param outcomes Character vector of length n. outcomes[i] is a string indicating the 
#' outcome experienced by unit i.
#' @param tr A directed igraph object. This is a tree representing the relationships
#' among the outcomes. The leaves represent individual outcomes, and internal nodes
#' represent outcome categories consisting of their leaf descendants. All nodes
#' of tr must have unique names as given by names(V(tr)). The names of the leaves must 
#' be equal to the unique elements of outcomes.
#' @param method = "matrix" or "tree". "matrix" uses a transformation of the 
#' design matrix to fit the MOReTreeS model; "tree" uses the information in tr. "matrix" may
#' be more efficient for small trees; "tree" may be more efficient for large trees. (?)
#' @param W_method = "shared" if information about the effect of variables in W wil be shared
#' across the outcomes according to the tree structure. If W_method = "individual", the effect of
#' W will be estimated separately for each outcome (no infromation sharing).
#' @param family A string specifying the distribution of the outcomes: 
#' either "bernoulli" (for classification) or "gaussian" (for regression)
#' @param ci_level A number between 0 and 1 giving the desired credible interval. 
#' For example, ci_level = 0.95 (the default) returns a 95\% credible interval
#' @param get_ml If TRUE, moretrees will also return the maximum likelihood estimates of the
#' coefficients for each outcome group discovered by the model. The default is FALSE.
#' @param tol Convergence tolerance for ELBO. Default = 1E-8.
#' @param maxiter Maximum number of iterations of the VI algorithm.
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
#' @param print_freq How often to print out iteration number. 
#' @return A list containing the following elements:
#' 1. estimated coefficients and credible intervals; 
#' 2. outputs from variational inference algorithm
#' @examples 
#' @family MOReTreeS functions

moretrees <- function(X, W = NULL, y, outcomes, tr,
                      random_init = FALSE,
                      initial_values = NULL,
                      method = "tree",
                      W_method = "shared",
                      family = "bernoulli",
                      ci_level = 0.95,
                      get_ml = FALSE,
                      update_hyper = T, update_hyper_freq = 50,
                      print_freq = update_hyper_freq,
                      hyper_fixed = NULL,
                      tol = 1E-8, max_iter = 5000,
                      nrestarts = 3,
                      keep_restarts = nrestarts > 1,
                      parallel = nrestarts > 1,
                      log_restarts = nrestarts > 1,
                      log_dir = getwd(),
                      hyper_random_init = list(omega_max = 100,
                                               tau_max = 100,
                                               sigma2_max = 100),
                      vi_random_init = list(eta_sd = 10,
                                            mu_sd = 10,
                                            delta_sd = 10)) {
  
  if (!(family %in% c("bernoulli", "gaussian"))) {
    stop("family must be a string (\"bernoulli\" or \"gaussian\")")
  }
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!is.null(W) & !(is.matrix(W))) stop("If W is not NULL, must be a matrix")
  if (family == "bernoulli" & method == "matrix") ss_fun <- spike_and_slab_logistic
  if (family == "gaussian" & method == "matrix") ss_fun <- spike_and_slab_normal
  if (family == "bernoulli" & method == "tree") ss_fun <- spike_and_slab_logistic_moretrees
  if (family == "gaussian" & method == "tree") ss_fun <- spike_and_slab_normal_moretrees
  if (!(length(get_ml) == 1 & is.logical(get_ml))) stop("get_ml must be either TRUE or FALSE")
  if (!update_hyper & is.null(hyper_fixed)) {
    stop("Must supply fixed hyperparameter values if update_hyper = FALSE")
  }
  
  # Get MOReTreeS design elements
  if (method == "matrix") {
    dsgn <- moretrees_design_matrix(X = X, W = W, y = y,
                                    outcomes = outcomes, tr = tr,
                                    W_method = W_method)
  }
  if (method == "tree") {
    dsgn <- moretrees_design_tree(X = X, W = W, y = y,
                                    outcomes = outcomes, tr = tr)
  }
  
  
  # Setting up parallelization
  if (parallel) {
    `%doRestarts%` <- foreach::`%dopar%`
  } else {
    `%doRestarts%` <- foreach::`%do%`
  }
  
  # # Run algorithm
  # mod_restarts <- foreach::foreach(i = 1:nrestarts) %doRestarts% {
  #   if (log_restarts) {
  #     sink(file = paste0(log_dir, "restart_", i, "_log.txt"))
  #     cat("Initialising random restart", i, "...\n\n")
  #   }
  #   mod <- 
  mod_restarts <- list(ss_fun(dsgn = dsgn,
           random_init = random_init,
           initial_values = initial_values,
           update_hyper = update_hyper, 
           update_hyper_freq = update_hyper_freq,
           print_freq = print_freq,
           hyper_fixed = hyper_fixed,
           tol = tol,
           max_iter = max_iter,
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
  
  # Compute MOReTreeS exposure coefficient estimates from model output
  betas <- moretrees_compute_betas(mod = mod, ci_level = ci_level,
                                   outcomes = outcomes,
              A_leaves = dsgn$A[names(igraph::V(tr))[igraph::V(tr)$leaf], ])
  
  # Compute MOReTreeS covariate coefficient estimates from model output
  if (!is.null(W)) {
    theta_est <- moretrees_compute_thetas(mod = mod, ci_level = ci_level, 
                                          m = ncol(W), W_method = W_method, method = method,
                                          A_leaves = dsgn$A[names(igraph::V(tr))[igraph::V(tr)$leaf], ])
  } else {
    theta_est <- NULL
  }
  
  
  # Get maximum likelihood estimates by group for comparison
  if (get_ml) {
    beta_ml <- ml_by_group(X = X, W = W, y = y, outcomes = outcomes,
                           outcome_groups = betas$beta_moretrees$outcomes,
                           ci_level = ci_level,
                           family = family)
  } else {
    beta_ml <- NULL
  }
  
  # Return results
  return(list(beta_est = betas$beta_est,
              beta_moretrees = betas$beta_moretrees,
              beta_ml = beta_ml, 
              theta_est = theta_est,
              mod = mod,
              mod_restarts = mod_restarts))
}