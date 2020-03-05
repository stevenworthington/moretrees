# --------------------------------------------------------------------------------- #
# ---------- Useful mathematical functions not implemented in base R  ------------- #
# --------------------------------------------------------------------------------- #

trace_prod <- function(A, B) {
  # computes trace(A %*% B) for square matrices A, B
  sum(Matrix::t(B) * A)
}

log1p_exp <- function(x) {
  # computes log(1 + exp(x))
  if (x > 20) {
    return(x)
  } else {
    return(log1p(exp(x)))
  }
}

log1p.exp.vec <- function(x){
  # computes log(1 + exp(x)) for x a vector
  y <- x
  which.small <- x <= 20
  y[which.small] <- log1p(exp(x[which.small]))
  return(y)
}

logexpit <- function(x) {
  # computes -log(1+exp(-x)) for x a vector
  -log1p.exp.vec(-x)
}

expit <- function(x) {
  # computes 1/(1+exp(-x)) for x a vector
  exp(logexpit(x))
}

gfun <- function(x) {
  # computes function g described in manuscript;
  # needed for normal approx to logistic likelihood
  (expit(x) - 1 / 2) / (2 * x)
}

xxT_g_eta_fun <- function(units, xxT, g_eta, K) {
  vec <- g_eta[units] %*% xxT[units, ]
  mat <- matrix(0, nrow = K, ncol = K)
  mat[lower.tri(mat, diag = T)] <- vec
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)] 
  return(mat)
}

xxT_ss_fun <- function(groups, dat) {
  X_g <- dat[ , groups, drop = F]
  if (ncol(X_g) == 1) {
    xxT <- X_g ^ 2
  } else {
    xxT <- rowOuterProds(X_g)
  }
  return(xxT)
}

xxT_g_eta_fun_ss <- function(xxT, K, g_eta) {
  vec <- g_eta %*% xxT
  mat <- matrix(0, nrow = K, ncol = K)
  mat[lower.tri(mat, diag = T)] <- vec
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)] 
  return(mat)
}

quadFormByRow <- function(Sigma, X) rowSums(tcrossprod(X, Sigma) * X)
