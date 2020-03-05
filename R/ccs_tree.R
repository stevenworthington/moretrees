# --------------------------------------------------------------------------------- #
# --------------- Code for creating tree of outcomes for CCS codes ---------------- #
# --------------------------------------------------------------------------------- #

#' \code{ccs_tree} Produces tree of Clinical Classification Software (CCS) codes
#' based on root node.
#' 
#' All the details go here!
#' 
#' @section Details
#' 
#' @export
#' @importFrom magrittr %>%
#' @param group character string specifying that only codes beginning with group 
#' will be included in tree The default is NULL (returns full tree of CCS codes)
#' @return A list containing the following elements:
#' tr = directed igraph object; the tree of outcomes.
#' ccs_icd_mapping = data frame specifying the ICD9 codes corresponding to each 
#' CCS code. Useful for converting from ICD9 to CCS.
#' @examples 
#' @family tree functions

ccs_tree <- function(group = NULL) {
  # Get data.frame showing mapping from ICD9 to multilevel CCS
  ccs_icd <- data.frame(icd = unlist(icd::icd9_map_multi_ccs[[1]]))
  for (i in 1:4) {
    ccs_list <- icd::icd9_map_multi_ccs[[i]]
    ccs_df <- data.frame(icd = unlist(ccs_list))
    ccs_df$ccs <-  ccs_list %>% 
      names %>% # names of the list entries are the CCS codes
      sapply(FUN = function(nm) rep(nm, length(ccs_list[[nm]]))) %>%
      unlist
    names(ccs_df)[names(ccs_df) == "ccs"] <- paste0("l", i)
    ccs_icd <- merge(ccs_icd, ccs_df, by = "icd", all.x = T, all.y = F)
  }
  
  # CCS codes only
  ccs <- ccs_icd
  ccs$icd <- NULL
  ccs <- ccs[!duplicated(ccs), ]
  
  # Order CCS codes appropriately
  ccs_levels <- matrix(nrow = nrow(ccs), ncol = 4)
  for (i in 1:nrow(ccs_levels)) {
    splt <- as.integer(stringr::str_split(ccs$l4[i], "\\.")[[1]])
    if (length(splt) != 4) {
      splt <- as.integer(stringr::str_split(ccs$l3[i], "\\.")[[1]])
      if (length(splt) != 3) {
        splt <- as.integer(stringr::str_split(ccs$l2[i], "\\.")[[1]])
        if (length(splt) != 2) {
          splt <- as.integer(ccs$l1[i])
        }
      }
      splt <- c(splt, rep(0, 4 - length(splt)))
    }
    ccs_levels[i, ] <- splt
  }
  ccs_levels <- as.data.frame(ccs_levels)
  names(ccs_levels) <- c("l1", "l2", "l3", "l4")
  ccs <- ccs[order(ccs_levels$l1, ccs_levels$l2, ccs_levels$l3, ccs_levels$l4), ]
  
  # Keep only diseases starting with group
  if (!is.null(group)) {
    group_length <- length(stringr::str_split(group, "\\.")[[1]])
    group_lvl <- paste0("l", group_length)
    ccs <- ccs[ccs[ , group_lvl] == group, ]
  } else {
    ccs <- cbind(l0 = "0", ccs)
    group_length <- 0
  }
  
  # Add zeros
  ccs_zeros <- ccs
  for (i in 1:nrow(ccs_zeros)) {
    which_zero <- which(ccs_zeros[i, ] == " ")
    for (j in sort(which_zero)) {
      ccs_zeros[i, j] <- paste0(ccs_zeros[i, j - 1], ".0")
    }
  }
  
  # Make tree
  ccs_mat <- ccs_zeros
  if (group_length > 1) {
    ccs_mat <- ccs_mat[ , sapply(group_length:4, function(x) paste0("l", x))]
  }
  ccs_mat <- as.matrix(ccs_mat)
  edges <- ccs_mat[ , c(1, 2)]
  if (ncol(ccs_mat) > 2) {
    for (i in 3:ncol(ccs_mat)) {
      edges <- rbind(edges, ccs_mat[ , c(i - 1, i)])
    }
  }
  edges <- edges[edges[ , 2] != " ", ]
  edges <- edges[!duplicated(edges), ]
  tr <- igraph::graph_from_edgelist(e = edges, directed = T)
  leaves <- igraph::V(tr)[igraph::degree(tr, mode = "out") == 0]
  igraph::V(tr)$leaf <- FALSE
  igraph::V(tr)$leaf[igraph::V(tr) %in% leaves] <- TRUE
  
  # ICD mapping
  ccs$id <- 1:nrow(ccs)
  names(ccs_zeros) <- sapply(names(ccs_zeros), function(nm) paste0(nm, "_0"))
  ccs_zeros$id <- 1:nrow(ccs)
  ccs$keep <- T
  ccs_icd <- merge(ccs_icd, ccs, by = c("l1", "l2", "l3", "l4"),
                   all.x = T, sort = FALSE)
  ccs_icd <- subset(ccs_icd, keep)
  ccs_icd <- merge(ccs_icd, ccs_zeros, by = "id",
                   all.x = T, sort = FALSE)
  ccs_icd$keep <- NULL
  ccs_icd$id <- NULL
  for (i in 2:4) {
    nm <- paste0("l", i)
    nm_p <- paste0("l", i - 1)
    ccs_icd[ccs_icd[ , nm] == " ", nm] <- ccs_icd[ccs_icd[ , nm] == " ", nm_p]
  }
  ccs_icd <- ccs_icd[ , c("icd", "l4", "l4_0")]
  names(ccs_icd)[names(ccs_icd) == "l4"] <- "ccs_original"
  names(ccs_icd)[names(ccs_icd) == "l4_0"] <- "ccs_added_zeros"
  
  return(list(tr = tr, ccs_icd_mapping = ccs_icd))
}