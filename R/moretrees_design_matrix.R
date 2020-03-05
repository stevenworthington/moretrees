# --------------------------------------------------------------------------------- #
# -------- Code for converting design matrix + tree into MOReTreeS ---------------- #
# -------- Design Matrix ---------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

#' Here's a brief description.
#'   \code{moretrees_design_matrix} converts outcome, exposure, and covariate data
#'   into format suitable for analysis using MOReTreeS.
#' 
#' All the details go here!
#' 
#' @export
#' @section Model Description:
#'   Describe group spike and slab prior and all parameters here.
#' 
#' @param y Vector of length n containing outcomes data.
#' If family = "bernoulli", y must be an integer vector where 1 = success, 0 = failure.
#' If family = "gaussian", y must be a numeric vector containing continuous data.
#' @param X An n x K matrix of exposure data, where K is the dimension of the exposure.
#' Grouping of the outcomes will be based on their relationships with the variables in X.
#' @param W Matrix of covariates of dimension n x m.
#' Coefficients for these variables do not affect grouping of the outcomes.
#' Default is NULL (no covariates).
#' @param outcomes is a character vector of length n, where entry i
#  tells us which outcome is represented by unit i
#' @param tr is an igraph tree, where the leaves represent outcomes
#' @param W_method = "shared" if information about the effect of variables in W wil be shared
#' across the outcomes according to the tree structure. If W_method = "individual", the effect of
#' W will be estimated separately for each outcome (no infromation sharing).
#' @return A list containing the following elements:
#' Xstar: a MOReTreeS design matrix of size n x (K * p) for the exposure. 
#' Note that the rows of Xstar[[i]] have been ro-ordered according to the vector ord,
#' and so may have a different ordering from X.
#' groups: a list of length p indicating the columns of X which correspond to the same
#' node on the tree. Used for spike and slab variable selection on the nodes.
#' Wstar: the MOReTreeS design matrices for covariates. If share_W = T, this will 
#' be a sparse Matrix of dimension n x (p * m), where m is the number of columns of W. 
#' If share_W = F, Wstar is just W with re-ordered rows. Wstar is NULL if W is NULL.
#' y_reord: Re-ordered outcome vector.
#' A: A sparse Matrix of dimension p x p, where p is the number of nodes in tr.
#' A_ij = 1 if j = i or node j is an ancestor of node i; A_ij = 0 otherwise.
#' ord = integer vector indicating order of rows in Xstar relative to input matrix X.
#' @examples
#' @family MOReTreeS functions

moretrees_design_matrix <- function(y, X, W = NULL, outcomes, tr, W_method = "shared") {
  # Some checks
  if (!is.character(outcomes)) stop("outcomes is not a character object")
  if (!igraph::is.igraph(tr)) stop("tr is not a graph object")
  if (!igraph::is.directed(tr)) stop
  if (!(W_method %in% c("shared", "individual"))) {
    stop("W_method must be either \"shared\" or \"individual\"")
  } 
  if (is.integer(y) & !(sum(y %in% c(0, 1)) == length(y))) 
    stop("y contains values other than zero or one")
  
  nodes <- names(igraph::V(tr))
  leaves <- names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0])
  if(!setequal(unique(outcomes), leaves)) {
    stop("Not all outcomes are leaves of tree")
  }
  # Re-order nodes to have internal nodes first, then leaves
  nodes <- c(nodes[!(nodes %in% leaves)], leaves)
  
  # Extract relevant parameters
  p <- length(nodes)
  pL <- length(leaves)
  K <- ncol(X)
  n <- nrow(X)
  
  # Get ancestor matrix
  A <- igraph::as_adjacency_matrix(tr, sparse = T)
  A <- A[nodes, nodes] # re-order rows/columns to mirror nodes
  A <- Matrix::expm(Matrix::t(A))
  A[A > 0 ] <- 1 
  A <- Matrix::Matrix(A, sparse = T)
  
  # Sort by outcomes, where order is specified by ordering in tr
  ord <- order(ordered(outcomes, levels = leaves))
  X <- X[ord, , drop = F]
  y <- y[ord]
  outcomes <- outcomes[ord]
  
  # Get list of MOReTreeS exposure design matrices for each node
  Xstar <- Matrix::Matrix(nrow = n, ncol = 0, sparse = T)
  # names(Xstar) <- nodes
  for (k in 1:K) {
    # Get design matrix for variable k
    Xmat_k <- sapply(leaves, function(v) X[outcomes == v, k], simplify = F)
    Xmat_k <- Matrix::bdiag(Xmat_k)
    Xmat_k <- cbind(Matrix::Matrix(0, nrow = n, ncol = p - pL), Xmat_k)
    Xstar <- cbind(Xstar, Xmat_k %*% A)
    rm(Xmat_k)
  }
  # Get list of variable groups for selection
  groups <- sapply(1:p, function(v) (1:K) * p - (p - v),simplify = F)
  # Remove X
  rm(X)
  
  # Get covariate design matrix
  if (!is.null(W)) {
    W <- W[ord, , drop = F]
    if (W_method == "shared") {
      m <- ncol(W)
      Wstar <- Matrix::Matrix(nrow = n, ncol = 0, sparse = T)
      for (j in 1:m) {
        # Get design matrix for variable j
        Wmat_j <- sapply(leaves, function(v) W[outcomes == v, j], simplify = F)
        Wmat_j <- Matrix::bdiag(Wmat_j)
        Wmat_j <- cbind(Matrix::Matrix(0, nrow = n, ncol = p - pL), Wmat_j)
        Wstar <- cbind(Wstar, Wmat_j %*% A)
        rm(Wmat_j)
      }
    } else {
      Wstar <- Matrix::bdiag(sapply(leaves, function(v) W[outcomes == v, ], simplify = F))
    }
    rm(W)
  } else {
    Wstar <- NULL
  }
  
  # Replace y = 0 with y = -1 for compatibility with moretrees algorithm
  if (is.integer(y)) y[y == 0] <- -1
  
  return(list(X = as.matrix(Xstar), groups = groups, W = as.matrix(Wstar), y = y, A = A))
}