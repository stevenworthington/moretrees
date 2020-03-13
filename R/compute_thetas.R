# --------------------------------------------------------------------------------- #
# ------- computing non-sparse estimates from spike and slab model output --------- #
# --------------------------------------------------------------------------------- #

#' Here's a brief description.
#' \code{compute_thetas} Computes non-sparse effect estimates and credible
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

compute_thetas <- function(mod, ci_level, method, dsgn) {
  
  theta_est <- data.frame(est = as.numeric(mod$vi_params$delta))
  theta_sd <- as.matrix(mod$vi_params$Omega) %>%
    diag %>% sqrt
  z <- qnorm(ci_level + (1 - ci_level) / 2)
  theta_est$cil <- theta_est$est - z * theta_sd
  theta_est$ciu <- theta_est$est + z * theta_sd
  
  # Return results
  return(theta_est)
  
}