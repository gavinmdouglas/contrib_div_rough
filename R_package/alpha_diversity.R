library(ape)


faiths_pd <- function(tips_in_sample, tree) {
  
  # Simplified version of equivalent functions from picante.
  # Requires tree to be rooted and will only compute distances
  # based on edges in tree after pruning to tips in samples
  # (i.e., distance to overall root is not included).
  
  if (length(tips_in_sample) == 1) { return(0) }
  
  if (is.null(tree$edge.length)) {
    stop("Stopping - no branch lengths in tree.")
  }
  
  if (! ape::is.rooted(tree)) {
    stop("Stopping - tree is unrooted.")
  }
  
  if (length(which(! tips_in_sample %in% tree$tip.label) > 0)) {
    stop("Stopping - some features in sample are not found as tips in the tree.")
  }
  
  if (length(which(! tree$tip.label %in% tips_in_sample) > 0)) {
    tree <- ape::drop.tip(phy = tree,
                          tip = tree$tip.label[which(! tree$tip.label %in% tips_in_sample)],
                          trim.internal = TRUE)
  }
  
  return(sum(tree$edge.length))

}


richness <- function(x) {
  length(which(x > 0))
}


shannon_index <- function(x) {
  x <- x[which(x > 0)]
  x <- x / sum(x)
  return(-1 * sum(x * log2(x)))
}


berger_parker_dominance <- function(x) {
  max(x) / sum(x)
}


ENS_pie <- function(x) {
  x <- x / sum(x)
  return(1 / sum(x**2))
}


fishers_alpha <- function(x) {
  
  # Simplified version of function implemented in vegan.
  # Returns 0 if fit cannot be performed
  # (which is assumed to always be the case when there are
  # fewer than three unique observations).
  
  x <- x[which(x > 0)]
  
  if (length(x) < 3) { return(0) }
  
  unique_obs <- length(x)
  total_obs <- sum(x)
  
  func2scan <- function(x, unique_obs, total_obs) { x * log(1 + total_obs / x) - unique_obs }
  
  options(warn = -1)
  uniroot_out <- try(uniroot(func2scan,
                             c(1, 50),
                             extendInt = "upX",
                             unique_obs = unique_obs, 
                             total_obs = total_obs)$root,
                    silent = TRUE)
  options(warn = 0)
  
  if (class(uniroot_out) == "try-error") {
    return(0)
  } else {
   return(uniroot_out)
  }
}


heaps_evenness <- function(x) {
  if (richness(x) == 1) { return(0) }
  (exp(shannon_index(x)) - 1) / (richness(x) - 1)
}


margalefs_richness <- function(x) {
  (richness(x) - 1) / log(sum(x))
}


mcintoshs_dominance <- function(x) {
  total_obs <- sum(x)
  (total_obs - sqrt(sum(x**2))) / (total_obs - sqrt(total_obs))
}


mcintoshs_evenness <- function(x) {
  x <- x[which(x > 0)]
  total_obs <- sum(x)
  unique_obs <- length(x)
  sqrt(sum(x**2)) / sqrt((total_obs - unique_obs + 1)**2 + unique_obs - 1)
}


menhinicks_richness <- function(x) {
  x <- x[which(x > 0)]
  length(x) / sqrt(sum(x))
}


pielous_evenness <- function(x) {
  x <- x[which(x > 0)]
  if (length(x) == 1) { return(0) }
  shannon_index(x) / log(length(x))
}


gini_simpson_index <- function(x) {
  x <- x / sum(x)
  return(1 - sum(x ** 2))
}


simpsons_evenness <- function(x) {
  x <- x / sum(x)
  return(inverse_simpson_index(x) / richness(x))
}


inverse_simpson_index <- function(x) {
  x <- x / sum(x)
  return(1 / sum(x ** 2))
}


# A list containing different functions.
# This is done for conveniently calling them programatically
calc_diversity <- list()
calc_diversity[["richness"]] <- richness
calc_diversity[["shannon_index"]] <- shannon_index
calc_diversity[["berger_parker_dominance"]] <- berger_parker_dominance
calc_diversity[["ENS_pie"]] <- ENS_pie
calc_diversity[["faiths_pd"]] <- faiths_pd
calc_diversity[["fishers_alpha"]] <- fishers_alpha
calc_diversity[["heaps_evenness"]] <- heaps_evenness
calc_diversity[["margalefs_richness"]] <- margalefs_richness
calc_diversity[["mcintoshs_dominance"]] <- mcintoshs_dominance
calc_diversity[["mcintoshs_evenness"]] <- mcintoshs_evenness
calc_diversity[["menhinicks_richness"]] <- menhinicks_richness
calc_diversity[["pielous_evenness"]] <- pielous_evenness
calc_diversity[["gini_simpson_index"]] <- gini_simpson_index
calc_diversity[["simpsons_evenness"]] <- simpsons_evenness
calc_diversity[["inverse_simpson_index"]] <- inverse_simpson_index
