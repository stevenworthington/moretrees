# --------------------------------------------------------------------------------- #
# ----- computing MOReTreeS estimates for exposure effects from model output ------ #
# --------------------------------------------------------------------------------- #

#' Here's a brief description.
#' \code{moretrees_compute_betas} Computes MOReTreeS estimates and credible
#'   intervals for exposure effects from model output
#' 
#' All the details go here!
#' 
#' @section Details
#' 
#' @export
#' @param mod List containing outputs from spike and slab VI algorithm
#' @param ci_level A number between 0 and 1 giving the desired credible interval.
#' For example, ci_level = 0.95 (the default) returns a 95\% credible interval.
#' @param A_leaves pL x p sparse ancestor Matrix where rows correspond to leaves
#' of tree (outcomes) and columns correspond to nodes on tree. Results in mod
#' must have same ordering as columns of A_leaves.
#' @return A list containing the following elements:
#' 1. beta_moretrees: estimated coefficients and credible intervals by group;
#' 2. beta_est: estimated coefficients and credible intervals by outcome;
#' @examples
#' @family MOReTreeS functions

moretrees_compute_betas <- function(mod, ci_level, A_leaves, outcomes) {
  
  # Get betas from xis
  p <- length(mod$vi_params$prob)
  node_select <- mod$vi_params$prob >= 0.5
  xi_est <- mapply(function(x, y) as.numeric(x * y),
                   mod$vi_params$mu, node_select,
                   SIMPLIFY = T) %>% t
  K <- nrow(mod$vi_params$mu[[1]])
  if (K == 1) xi_est <- t(xi_est)
  beta_est <- as.matrix(A_leaves %*% xi_est) %>%
    as.data.frame
  beta_names <- sapply(1:K, function(i) paste0("est",i))
  colnames(beta_est) <- beta_names
  
  # Compute credible intervals
  xi_var_est <- mapply(function(Sigma, node_select) diag(as.matrix(Sigma)) * node_select,
                       mod$vi_params$Sigma, node_select, SIMPLIFY = T) %>% t
  if (K == 1) xi_var_est <- t(xi_var_est)
  beta_sd_est <- as.matrix(A_leaves %*% xi_var_est) %>% sqrt
  z <- qnorm(ci_level + (1 - ci_level) / 2)
  beta_ci_l <- beta_est - z * beta_sd_est
  colnames(beta_ci_l) <- sapply(1:K, function(i) paste0("cil",i))
  beta_ci_u <- beta_est + z * beta_sd_est
  colnames(beta_ci_u) <- sapply(1:K, function(i) paste0("ciu",i))
  beta_est <- cbind(beta_est, beta_ci_l, beta_ci_u)
  beta_est$group <- as.numeric(beta_est$est1) %>%
    as.factor %>% as.integer
  
  # Get estimated coefficients and CIs by group
  beta_moretrees <- beta_est[!duplicated(beta_est), ]
  beta_moretrees <- beta_moretrees[order(beta_moretrees$group), ]
  row.names(beta_moretrees) <- NULL
  beta_moretrees$outcomes <- vector(mode = "list", length = nrow(beta_moretrees))
  G <- nrow(beta_moretrees)
  for (i in 1:G) {
    beta_moretrees$outcomes[[i]] <- row.names(beta_est)[beta_est$group == i]
  }
  
  # Get counts of observations by group
  beta_moretrees$n_obs <- sapply(beta_moretrees$outcomes,
                       function(out, dat) sum(dat %in% out),
                       dat = outcomes)
  
  # re-order columns for readability
  cols <- c("est", "cil", "ciu")
  cols <- sapply(1:K, function(i) paste0(cols, i), simplify = T) %>%
    as.vector
  beta_moretrees <- beta_moretrees[ , c("group", cols, "n_obs", "outcomes")]
  
  # Return results
  return(list(beta_moretrees = beta_moretrees, beta_est = beta_est))
  
}