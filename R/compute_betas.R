# --------------------------------------------------------------------------------- #
# --------- computing sparse estimates from spike and slab model output ----------- #
# --------------------------------------------------------------------------------- #

#' Here's a brief description.
#' \code{compute_betas} Computes sparse effect estimates and credible
#'   intervals from spike & slba model output
#' 
#' All the details go here!
#' 
#' @section Details
#' 
#' @export
#' @param mod List containing outputs from spike and slab VI algorithm
#' @param ci_level A number between 0 and 1 giving the desired credible interval.
#' For example, ci_level = 0.95 (the default) returns a 95\% credible interval.
#' @return A matrix containing estimated coefficients and credible intervals.
#' @examples 
#' @family spike and slab functions

compute_betas <- function(mod, ci_level) {
  
  # Compute estimated betas
  node_select <- mod$vi_params$prob >= 0.5
  beta_est <- mapply(function(x, y) as.vector(x * y),
                   mod$vi_params$mu, node_select,
                   SIMPLIFY = F)
  
  # Compute credible intervals
  beta_sd_est <- mapply(function(Sigma, node_select) sqrt(diag(as.matrix(Sigma)) * node_select),
                       Sigma = mod$vi_params$Sigma, 
                       node_select = node_select, 
                       SIMPLIFY = F)
  z <- qnorm(ci_level + (1 - ci_level) / 2)
  beta_ci_l <- mapply(function(beta, sd) beta - z * sd,
                      beta = beta_est, sd = beta_sd_est)
  beta_ci_u <- mapply(function(beta, sd) beta + z * sd,
                      beta = beta_est, sd = beta_sd_est)
  
  # Combine results
  beta_est <- mapply(function(est, cil, ciu) data.frame(est = est, cil = cil, ciu = ciu),
                     est = beta_est, cil = beta_ci_l, ciu = beta_ci_u,
                     SIMPLIFY = F)
  beta_names <- sapply(1:length(beta_est), function(i) paste0("group",i))
  names(beta_est) <- beta_names
  
  # Return results
  return(beta_est)
  
}