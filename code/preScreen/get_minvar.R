# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
# get.minvar
#
# A function to adaptively determine a minimum variability cutoff for a given
# sieve analysis.  This function only works for binary or dichotomous features.
#
# Arguments (all default values pertain to the context of the HVTN 702 sieve 
#            analysis):
#   -- "n":  the number of participants in the study 
#   -- "n.vx":  the number of study participants in the vaccine/treatment arm 
#   -- "p":  the total number of features being considered for the analysis
#   -- "fwer":  the target FWER error rate significance level
#
# Example:  get.minvar (n=234, n.vx=130, p=400, fwer=0.05)
#
#   (The above example, which is using the default settings, returns the value 
#   16, indicating that the minimum variability threshold should be 16 in order 
#   to have a chance of achieving a significant result with an FWER-adjusted
#   p-value less than 0.05.)
#
# -- C.A. Magaret
#    cmagaret@fredhutch.org
#    24 October 2023
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


get.minvar <- function (n=234, n.vx=130, p=400, fwer=0.05) {
  #browser()
  # set up a vector of our minvar range
  minvar <- 2:n
  
  # loop through each minvar cutoff
  for (iter in 1:(length (minvar) - 1)) {
    
    # set up our idealized data
    data.tmp <- rep (0, n)
    trt.tmp <- rep (0, n)
    trt.tmp[1:n.vx] <- 1
    
    # set up the treatment differences by our minority class so that we catch
    # cases break along the diagnal; all minority cases in one treatment arm and all majority cases in the other treatment arm
    # the least-significant option
    # Li : Since all cases are assigned to the minority mark group, the VE difference between match vs. mismatch would be the greatest. 
    if (n.vx <= n / 2) {
      data.tmp[n:(n - minvar[iter] + 1)] <- 1
    } else {
      data.tmp[1:minvar[iter]] <- 1
    }
    
    # determine our p-value and return the cutoff if we reach FWER significance
    p.val.tmp <- fisher.test (table (data.tmp, trt.tmp))$p.value
    if (p.val.tmp * p < fwer) {
      return (minvar[iter])
    }
  }
  
  # otherwise it looks like we never reach significance
  return (NA)
}


# Variability filter functions --------------------------------------------

#' \code{varFilterThresFisher} returns the minimum number of events of each type needed for Fisher's exact test to detect a sieve effect (p < \code{pThres}) 
#' if all failures of one type were in one treatment arm
#'
#' A 2x2 contingency table underlies the procedure, i.e., a binary mark and 2 treatment arms are required.
#' 
#' @param eventIndWithMark an indicator of an event with an observed mark (events with a missing mark are coded as 0)
#' @param tx an indicator of one of the two treatment arms. The arguments \code{eventIndWithMark} and \code{tx} are only used to 
#' compute the total number of events with an observed mark in each arm.
#' @param pThres Fisher's exact test p-value threshold used for defining a sieve effect
#'
#' @return The minimum required total number of events of each type.
varFilterThresFisher <- function(eventIndWithMark, tx, pThres=0.05){
  n <- tapply(eventIndWithMark, tx, sum)
  
  if (any(n<5)){
    idxP <- idxV <- numeric(0)
    while (length(c(idxP, idxV))==0){
      n <- n + 1
      
      pFisher <- t(sapply(1:min(n, 50), function(i, n){
        splitP <- matrix(c(n[1]-i, i, n[2], 0), nrow=2)
        pSplitP <- fisher.test(splitP)$p.value
        
        splitV <- matrix(c(n[1], 0, n[2]-i, i), nrow=2)
        pSplitV <- fisher.test(splitV)$p.value
        
        return(c(i, pSplitP, pSplitV))
      }, n=n))
      
      idxP <- which(pFisher[, 2] < pThres)
      idxV <- which(pFisher[, 3] < pThres)
    }
  } else {
    pFisher <- t(sapply(1:min(n, 50), function(i){
      # contingency table:
      #         P        V
      # Type 1  n[1]-i   n[2]
      # Type 2  i        0
      # all events of type 2 are in the P arm
      splitP <- matrix(c(n[1] - i, i, n[2], 0), nrow=2)
      pSplitP <- fisher.test(splitP)$p.value
      
      # contingency table:
      #         P        V
      # Type 1  n[1]     n[2]-i
      # Type 2  0        i
      # all events of type 2 are in the V arm
      splitV <- matrix(c(n[1], 0, n[2]-i, i), nrow=2)
      pSplitV <- fisher.test(splitV)$p.value
      
      return(c(i, pSplitP, pSplitV))
    })) 
  }
  
  return(min(which(pFisher[, 2] < pThres), which(pFisher[, 3] < pThres)))
}

