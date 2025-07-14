# ====================================================================
# SSD Method Implementation (Hagar, Maleyeff et al.)
# 
# This script contains R code to implement the Sample Size Determination (SSD)
# method proposed in the manuscript:
# "An Efficient Approach to Design Bayesian Platform Trials with an Application 
#  to the SSTaRLeT Trial" by Hagar, Maleyeff et al. (2025). 
#
# The code was developed by Luke Hagar and Lara Maleyeff.
# ====================================================================

# Helper functions
expit <- function(x){1/(1 + exp(-x))}
logit <- function(x){log(x) - log(1-x)}

#' Compute posterior weights for a robust mixture prior using marginal likelihoods
#' 
#' This function updates the weights of a robust mixture prior composed of several
#' informative1 Beta components and a weakly informative1 component. The update is
#' based on newly observed binomial data (y successes out of n individuals). 
#'
#' @param w Vector of weights for the informative Beta components (sum to 1)
#' @param w_inf Overall weight on the informative portion (between 0 and 1)
#' @param y Number of successes in the current data
#' @param n Number of total observations in the current data
#' @param a Vector of alpha parameters for the informative Beta components
#' @param b Vector of beta parameters for the informative Beta components
#' @param a0 Alpha parameter for the weakly informative Beta component (default = 1)
#' @param b0 Beta parameter for the weakly informative Beta component (default = 1)
#' 
#' @return A vector of updated posterior weights for the informative components
wpost = function(w, w_inf, y, n, a, b, a0 = 1, b0 = 1) {
  p = NULL
  H = length(w)
  for (h in 1:H){
    p[h] = beta(a[h] + y, b[h] + n - y)/beta(a[h], b[h])
  }
  p0 = beta(a0 + y, b0 + n - y)/beta(a0, b0)
  wp = w_inf*w*p/(sum(w_inf*w*p) + (1-w_inf)*p0)
  return(wp)
}

#' Sample from the posterior distribution under a finite mixture of Beta priors
#' 
#' This function simulates posterior draws from a robust MAP prior consisting of a 
#' weighted mixture of informative1 Beta distributions and a weakly informative component.
#'
#' @param w Vector of updated posterior weights for the informative components
#' @param a Vector of alpha parameters for the informative components
#' @param b Vector of beta parameters for the informative components
#' @param y Number of successes in the current data
#' @param n Number of total observations in the current data
#' @param M Number of posterior samples to draw (default = 1000)
#' 
#' @return A vector of M samples from the robust posterior distribution
betaFMsamp = function(w, a, b, y, n, M = 1000) {
  x = rmultinom(M, 1, p = c(w, 1-sum(w)))
  a = c(a, 1)
  b = c(b, 1)
  as = t(x) %*% a + y
  bs = t(x) %*% b + n - 1
  s = mapply(rbeta, shape1 = as, shape2 = bs, MoreArgs = list(n = 1))
  return(c(s))
}

#' Estimate Posterior Probability of Superiority
#'
#' This function estimates the posterior probability that a treatment effect exceeds 
#' a prespecified non-inferiority or superiority margin, using a log-odds transformation 
#' and kernel density estimation to account for uncertainty in the effect estimate.
#'
#' @param treatment Numeric vector of treatment group event probabilities.
#' @param control Numeric vector of control group event probabilities.
#' @param lmargin Numeric value specifying the threshold on the log-odds scale above which 
#' the treatment is considered superior or non-inferior.
#'
#' @return A scalar representing the estimated posterior probability that the treatment 
#' effect exceeds the specified margin.
estimate_pp <- function(treatment, control, lmargin) {
  ldiff <- log(treatment - control + 1) - log(1 - (treatment - control))
  kd <- density(na.omit(ldiff))
  prob <- mean(pnorm(lmargin, ifelse(is.finite(ldiff), ldiff, max(na.omit(ldiff)) + 1), kd$bw, lower.tail = FALSE))
  if (prob < 0.00004) prob <- pnorm(lmargin, mean(na.omit(ldiff)), sd(na.omit(ldiff)), lower.tail = FALSE)
  return(prob)
}


#' Simulate a three-arm (plus control) Bayesian platform trial with interim analysis and arm adaptation
#'
#' This function simulates a platform trial with one control and two active treatment arms,
#' followed by optional addition of a third treatment. The design includes an interim analysis
#' where one or both arms can be dropped based on inferiority, and posterior probabilities are
#' computed at both interim and final stages using informative or uninformative priors.
#'
#' @param i Integer; replicate index (used for parallel simulation).
#' @param cRates Vector of control arm event rates for each endpoint: AE, non-completion, and non-tolerability.
#' @param eff Vector of treatment effects (added to control rates) applied equally to Trt1, Trt2, and Trt3.
#' @param n.int Integer; number of patients accrued by interim analysis trigger.
#' @param n.final Integer; total number of patients at final analysis.
#' @param M Integer; number of posterior draws per parameter (default = 1000).
#' @param margin Vector of non-inferiority/superiority margins for each endpoint.
#' @param alloc Initial vector of allocation ratios before interim trigger (e.g., c(1,2,2) for 1:2:2 Control:Trt1:Trt2).
#' @param a0, b0 Scalars; hyperparameters for weakly informative Beta prior.
#' @param a1, b1 Vectors; hyperparameters for informative1 (4R10) Beta prior.
#' @param w1 Vector of prior weights for components of informative1 prior.
#' @param a2, b2 Vectors; hyperparameters for informative2 (2R20) Beta prior.
#' @param w2 Vector or scalar; prior weights for informative2 prior.
#' @param w_inf Scalar; total weight on informative components for robust MAP prior.
#' @param informative1 Logical; whether to use informative prior on control (4R10) AE rate.
#' @param informative2 Logical; whether to use informative prior on Trt1 (2R20) AE rate.
#'
#' @return A named numeric vector of length 30:\cr
#' - First 6 values (`tau1`) are posterior probabilities of inferiority at interim (3 endpoints × 2 arms).\cr
#' - Next 24 values (`tau2`) are posterior probabilities at final analysis under 4 possible drop scenarios:
#'     \itemize{
#'       \item 0: No arms dropped (Trt1, Trt2, Trt3 vs. control).
#'       \item 1: Trt2 dropped (Trt1, Trt3 vs. control).
#'       \item 2: Trt1 dropped (Trt2, Trt3 vs. control).
#'       \item 3: Both dropped (Trt3 vs. control only).
#'     }
#' Probabilities are calculated for AE, completion, and non-tolerability endpoints.
oneIAaddArm <- function(i, cRates = c(0.02, 0.25, 0.25), eff = c(0, 0, 0), 
                        n.int = 700, n.final = 1800, M = 1000, 
                        margin = c(0.04, 0.1, 0.1), alloc = c(1, 2, 2), 
                        a0 = 1, b0 = 1, 
                        a1 = c(3, 16, 36, 12), b1 = c(57, 379, 2853, 430), w1 = rep(1/4, 4),
                        a2 = c(9), b2 = c(434), w2 = 1, 
                        w_inf = 0.5, informative1 = TRUE, informative2 = TRUE) {
  
  # Transform margins to log-odds scale
  lmargin <- log(margin + 1) - log(1 - margin)
  
  # Compute outcome rates for Trt1, Trt2, and Trt3
  tRates = cRates + eff
  
  # Allocation before interim
  prand <- alloc / sum(alloc)
  x.int <- round(n.int * prand)

  # Define timing for drop decision and subsequent accrual
  n.decision = n.int + 300
  if (n.decision >= n.final) {
    n.decision = n.int
  }
  x.mid.total <- n.decision - n.int
  x.final.total <- n.final - n.decision
  
  # Simulate pre-interim events
  yc.int = sapply(cRates, function(z) rbinom(1, x.int[1], z))
  yt1.int = sapply(tRates, function(z) rbinom(1, x.int[2], z))
  yt2.int = sapply(tRates, function(z) rbinom(1, x.int[3], z))
  
  # Posterior samples for control pre-interim (AE only)
  thetac.int <- if (informative1) {
    wp.int <- wpost(w1, w_inf, yc.int[1], x.int[1], a1, b1, a0, b0)
    betaFMsamp(wp.int, a1, b1, yc.int[1], x.int[1], M)
  } else {
    rbeta(M, 1 + yc.int[1], 1 + x.int[1] - yc.int[1])
  }
  
  # Posterior samples for Trt1 pre-interim (AE only)
  thetat1.int <- if (informative2) {
    wp.int <- wpost(w2, w_inf, yt1.int[1], x.int[2], a2, b2, a0, b0)
    betaFMsamp(wp.int, a2, b2, yt1.int[1], x.int[2], M)
  } else {
    rbeta(M, 1 + yt1.int[1], 1 + x.int[2] - yt1.int[1])
  }
  
  # Trt2 posterior (uninformative)
  thetat2.int = rbeta(M, 1 + yt2.int[1], 1 + x.int[3] - yt2.int[1])
  
  # Interim inferiority probabilities
  tau1 <- numeric(6)
  tau1[1] <- estimate_pp(thetat1.int, thetac.int, lmargin[1])
  tau1[4] <- estimate_pp(thetat2.int, thetac.int, lmargin[1])
  
  for (k in 2:3) {
    thetac.int.k <- rbeta(M, 1 + yc.int[k], 1 + x.int[1] - yc.int[k])
    thetat1.int.k <- rbeta(M, 1 + yt1.int[k], 1 + x.int[2] - yt1.int[k])
    thetat2.int.k <- rbeta(M, 1 + yt2.int[k], 1 + x.int[3] - yt2.int[k])
    
    tau1[k] <- estimate_pp(thetat1.int.k, thetac.int.k, lmargin[k])
    tau1[k+3] <- estimate_pp(thetat2.int.k, thetac.int.k, lmargin[k])
  }
  
  names(tau1) <- c("tau1_Trt1_AE", "tau1_Trt1_Comp", "tau1_Trt1_NTol",
                   "tau1_Trt2_AE", "tau1_Trt2_Comp", "tau1_Trt2_NTol")
  

  # Final analysis for all drop scenarios
  tau2.list <- list()
  drop.labels <- c("tau2_no_drop", "tau2_trt2_dropped", "tau2_trt1_dropped", "tau2_both_dropped")
  
  for (drop in 0:3) {
    active1 <- drop %in% c(0, 1)
    active2 <- drop %in% c(0, 2)

    # Allocation strategy depends on which arms are active
    alloc_mid <- rep(0.5 / 3, 3); alloc_mid[4] <- 0.5
    alloc_final <- switch(as.character(drop),
                          "0" = c(rep(0.5 / 3, 3),0.5),
                          "1" = c(0.25, 0.25, 0, 0.5),
                          "2" = c(0.25, 0, 0.25, 0.5),
                          "3" = c(0.5, 0, 0, 0.5)
    )
    
    x.final <- round(x.final.total * alloc_final) + round(x.mid.total * alloc_mid)
    yc.final <- sapply(cRates, function(z) rbinom(1, x.final[1], z)) 
    yt1.final <- sapply(tRates, function(z) rbinom(1, x.final[2], z))
    yt2.final <- sapply(tRates, function(z) rbinom(1, x.final[3] , z))
    yt3.final <- sapply(tRates, function(z) rbinom(1, x.final[4], z))
    
    wp.final <- wpost(w1, w_inf, yc.int[1] + yc.final[1], x.int[1] + x.final[1], a1, b1, a0, b0)
    
    tau2_vec <- c()
    tau2_labels <- c()
    
    # Trt1 comparisons
    if (active1) {
      for (k in 1:3) {
        thetac.k <- if (informative1 & k == 1) {
          betaFMsamp(wp.final, a1, b1, yc.int[1] + yc.final[1], x.int[1] + x.final[1], M)
        } else {
          rbeta(M, 1 + yc.int[k] + yc.final[k], 1 + x.int[1] + x.final[1] - yc.int[k] - yc.final[k])
        }
        
        thetat1.k <- if (informative2 & k == 1) {
            wp <- wpost(w2, w_inf, yt1.int[1] + yt1.final[1], x.int[2] + x.final[2], a2, b2, a0, b0)
            betaFMsamp(wp, a2, b2, yt1.int[1] + yt1.final[1], x.int[2] + x.final[2], M)
          } else {
            rbeta(M, 1 + yt1.int[k] + yt1.final[k], 1 + x.int[2] + x.final[2] - yt1.int[k] - yt1.final[k])
          }
        tau2_vec <- c(tau2_vec, estimate_pp(thetat1.k, thetac.k, lmargin[k]))
        tau2_labels <- c(tau2_labels, paste0(drop.labels[drop+1], "_Trt1_", c("AE", "Comp", "NTol")[k]))
      }
    }
    
    # Trt2 comparisons
    if (active2) {
      for (k in 1:3) {
        thetac.k <- if (informative1 & k == 1) {
          betaFMsamp(wp.final, a1, b1, yc.int[1] + yc.final[1], x.int[1] + x.final[1], M)
        } else {
          rbeta(M, 1 + yc.int[k] + yc.final[k], 1 + x.int[1] + x.final[1] - yc.int[k] - yc.final[k])
        }
        thetat2.k <- rbeta(M, 1 + yt2.int[k] + yt2.final[k], 1 + x.int[3] + x.final[3] - yt2.int[k] - yt2.final[k])
        tau2_vec <- c(tau2_vec, estimate_pp(thetat2.k, thetac.k, lmargin[k]))
        tau2_labels <- c(tau2_labels, paste0(drop.labels[drop+1], "_Trt2_", c("AE", "Comp", "NTol")[k]))
      }
    }
    
    # Trt3 comparisons
    for (k in 1:3) {
      thetac.k <- if (informative1 & k == 1) {
        betaFMsamp(wp.final, a1, b1, yc.int[1] + yc.final[1], x.int[1] + x.final[1], M)
      } else {
        rbeta(M, 1 + yc.int[k] + yc.final[k], 1 + x.int[1] + x.final[1] - yc.int[k] - yc.final[k])
      }
      thetat3.k <- rbeta(M, 1 + yt3.final[k], 1 + x.final[4] - yt3.final[k])
      tau2_vec <- c(tau2_vec, estimate_pp(thetat3.k, thetac.k, lmargin[k]))
      tau2_labels <- c(tau2_labels, paste0(drop.labels[drop+1], "_Trt3_", c("AE", "Comp", "NTol")[k]))
    }
    
    names(tau2_vec) <- tau2_labels
    tau2.list[[drop+1]] <- tau2_vec
  }
  
  return(c(tau1,  unlist(tau2.list)))
}

#' Estimate operating characteristics using linear approximations
#' 
#' This function estimates various operating characteristics of the trial as
#' a function of the sample size. Note: we use interim sample size in the examples.
#' 
#' @param ml Matrix of posterior probabilities at a smaller sample size (e.g., n = 600).
#' @param mu Matrix of posterior probabilities at a larger sample size (e.g., n = 1000).
#' @param nl Integer; sample size corresponding to ml.
#' @param nu Integer; sample size corresponding to mu.
#' @param lb Lower bound for the sample size range to evaluate.
#' @param ub Upper bound for the sample size range to evaluate.
#' @param gamma Vector of inferiority thresholds (length 3).
#' @param kappa Vector of non-inferiority thresholds (length 3).
#'
#' 
#' @return A matrix with operating characteristics for each sample size (rows):
#'   - Column 1: interim sample size in control group
#'   - Columns 2-4: probability of dropping Trt1, Trt2, or both at interim
#'   - Columns 5-13: final-stage success probabilities for AE, completion, and tolerability
#'     for Trt1, Trt2, Trt3
#'   - Column 14: probability that any treatment demonstrates success on AE
#' 
getOCsLin <- function(ml, mu, nl, nu, lb, ub, gamma, kappa){
  
  logit <- function(x){log(x) - log(1-x)}
  l1 <- function(x) -logit(x)  # Used for complementary probabilities
  
  # Logit-transform posterior probabilities for interim (1:6) and final (7:30) stages
  ll.gamma <- apply(ml[,1:6], 2, logit)
  ll.kappa <- apply(ml[,7:30], 2, l1)
  lu.gamma <- apply(mu[,1:6], 2, logit)
  lu.kappa <- apply(mu[,7:30], 2, l1)
  
  ll <- cbind(ll.gamma, ll.kappa)
  lu <- cbind(lu.gamma, lu.kappa)
  
  # Replace -Inf and Inf with finite values based on data range
  for (j in 1:ncol(ll)){
    ll[,j] <- ifelse(ll[,j] == -Inf, min(subset(ll, is.finite(ll[,j]))[,j]) - 1, ll[,j])
    ll[,j] <- ifelse(ll[,j] == Inf, max(subset(ll, is.finite(ll[,j]))[,j]) + 1, ll[,j])
  }
  for (j in 1:ncol(lu)){
    lu[,j] <- ifelse(lu[,j] == -Inf, min(subset(lu, is.finite(lu[,j]))[,j]) - 1, lu[,j])
    lu[,j] <- ifelse(lu[,j] == Inf, max(subset(lu, is.finite(lu[,j]))[,j]) + 1, lu[,j])
  }
  
  # Add row index for tracking
  ll <- cbind(ll, seq(1, nrow(ll), 1))
  lu <- cbind(lu, seq(1, nrow(lu), 1))
  
  # Linear interpolation on the logit scale
  slopes <- ints <- NULL
  for (j in 1:ncol(ml)){
    ll_s <- ll[order(ll[,j]),j]
    lu_s <- lu[order(lu[,j]),j]
    l_slope <- (lu_s - ll_s)/(nu-nl)
    l_int <- ll_s - l_slope*nl
    l_slope[ll[order(ll[,j]),ncol(ll)]] <- l_slope
    l_int[ll[order(ll[,j]),ncol(ll)]] <- l_int
    slopes <- cbind(slopes, l_slope)
    ints <- cbind(ints, l_int)
  }

  ns <- seq(lb, ub,1)
  
  for (j in 1:ncol(ml)){
    temp.mat <- matrix(0, nrow = nrow(slopes), ncol=length(ns))
    for (i in 1:length(ns)){
      temp.mat[,i] <- ints[,j] + slopes[,j]*ns[i]
    }
    assign(paste0("l", j), temp.mat)
  }
  
  for (j in 1:6) {
    assign(paste0("y", j), get(paste0("l", j)) > logit(gamma[(j-1) %% 3 + 1]))
  }

  for (j in 7:30) {
    assign(paste0("y", j), get(paste0("l", j)) > logit(kappa[(j-7) %% 3 + 1]))
  }

  # At least one of the three endpoints needs to be inferior to stop first treatment
  stop1 <- colMeans(y1 + y2 + y3 > 0)
  # At least one of the three endpoints needs to be inferior to stop second treatment
  stop2 <- colMeans(y4 + y5 + y6 > 0)
  stopboth <- colMeans((((y1 + y2 + y3) > 0) + ((y4 + y5 + y6 > 0))) > 1)
  
  nrep <- nrow(get("y1"))
  n.col <- ncol(get("y1"))

  # drop0 to drop3 are matrices (one row per simulation repetition
  # and one column per sample size) which denote which arms are 
  # dropped for each sample size and simulation repetition combination
  stop1mat <- (y1 + y2 + y3) > 0
  stop2mat <- (y4 + y5 + y6) > 0
  drop0 <- ((stop1mat == 0) & (stop2mat == 0))
  drop1 <- ((stop1mat == 0) & (stop2mat == 1))
  drop2 <- ((stop1mat == 1) & (stop2mat == 0))
  drop3 <- ((stop1mat == 1) & (stop2mat == 1))
  
  # Estimate operating characteristics
  ae1 <- ifelse(drop0, y7, 
                ifelse(drop1, y16,
                      ifelse(drop2, matrix(0, nrow = nrep, ncol = n.col),
                             matrix(0, nrow = nrep, ncol = n.col))))
  
  comp1 <- colMeans(ifelse(drop0, y8, 
                    ifelse(drop1, y17,
                          ifelse(drop2, matrix(0, nrow = nrep, ncol = n.col),
                                matrix(0, nrow = nrep, ncol = n.col)))))
  
  ntol1 <- colMeans(ifelse(drop0, y9, 
                    ifelse(drop1, y18,
                           ifelse(drop2, matrix(0, nrow = nrep, ncol = n.col),
                                  matrix(0, nrow = nrep, ncol = n.col)))))
  
  ae2 <- ifelse(drop0, y10, 
                ifelse(drop2, y22,
                       ifelse(drop1, matrix(0, nrow = nrep, ncol = n.col),
                          matrix(0, nrow = nrep, ncol = n.col))))
  
  comp2 <- colMeans(ifelse(drop0, y11, 
                    ifelse(drop2, y23,
                           ifelse(drop1, matrix(0, nrow = nrep, ncol = n.col),
                              matrix(0, nrow = nrep, ncol = n.col))))) 
  
  ntol2 <- colMeans(ifelse(drop0, y12, 
                    ifelse(drop2, y24,
                           ifelse(drop1, matrix(0, nrow = nrep, ncol = n.col),
                                matrix(0, nrow = nrep, ncol = n.col))))) 
  
  ae3 <- ifelse(drop0, y13, 
                ifelse(drop1, y19,
                       ifelse(drop2, y25, y28)))
  
  comp3 <- colMeans(ifelse(drop0, y14, 
                    ifelse(drop1, y20,
                           ifelse(drop2, y26, y29))))
  
  ntol3 <- colMeans(ifelse(drop0, y15, 
                    ifelse(drop1, y21,
                           ifelse(drop2, y27, y30))))
  
  any <- ae1 + ae2 + ae3 > 0
  
  ae1 <- colMeans(ae1)
  ae2 <- colMeans(ae2)
  ae3 <- colMeans(ae3)
  any <- colMeans(any)
  
  # Combine these into one matrix; the first column denotes the sample size for that row
  return(cbind(ns, stop1, stop2, stopboth, ae1, comp1, ntol1, ae2, comp2, ntol2,
               ae3, comp3, ntol3,any))
}

#' Apply decision rules to simulated trial outputs
#'
#' This function computes estimated operating characteristics from a single sampling distribution 
#' estimate by applying decision rules for early stopping and final success. These rules are based on
#' user-specified posterior probability thresholds for inferiority (gamma) and success (kappa).
#'
#' @param trial.outputs A numeric matrix with 30 columns corresponding to posterior probabilities from 
#'                      each simulated trial. The first 6 columns (tau1) represent interim probabilities 
#'                      for Trt1 and Trt2 across AE, completion, and non-tolerability endpoints. The 
#'                      remaining 24 columns (tau2) are final probabilities under four drop scenarios.
#' @param gamma Vector of length 3 specifying thresholds for early stopping at interim analysis.
#' @param kappa Vector of length 3 specifying thresholds for final-stage success for each endpoint.
#' @param n Integer; sample size used for labeling the results.
#'
#' @return A named numeric vector with estimated operating characteristics:
#'   - ns: sample size
#'   - stop1, stop2, stopboth: probabilities of stopping Trt1, Trt2, or both
#'   - ae1, comp1, ntol1: final success probabilities for Trt1 on AE, completion, and non-tolerability
#'   - ae2, comp2, ntol2: final success probabilities for Trt2
#'   - ae3, comp3, ntol3: final success probabilities for Trt3 (always included)
#'   - any: probability that any treatment achieves success on AE
decisions = function(trial.outputs, gamma = c(0.2, 0.5, 0.5), kappa = c(0.95, 0.99, 0.99), n){
  
  # Copy of trial output matrix
  z = trial.outputs
  # Flip direction of tau2 probabilities (columns 7–30) to reflect success (not inferiority)
  for (j in 7:ncol(trial.outputs)){
    z[,j] <- 1 - z[,j]
  }
  
  # Apply stopping rules at interim (tau1, columns 1–6)
  for (j in 1:6) {
    assign(paste0("y", j), z[,j] > gamma[(j-1) %% 3 + 1])
  }
  # Apply success rules at final (tau2, columns 7–30)
  for (j in 7:30) {
    assign(paste0("y", j), z[,j] > kappa[(j-7) %% 3 + 1])
  }
  
  # At least one of the three endpoints needs to be inferior to stop first treatment
  stop1 <- mean(y1 + y2 + y3 > 0)
  # At least one of the three endpoints needs to be inferior to stop second treatment
  stop2 <- mean(y4 + y5 + y6 > 0)
  
  # Add a tracker for stopping both arms
  stopboth <- mean((((y1 + y2 + y3) > 0) + ((y4 + y5 + y6 > 0))) > 1)
  
  nrep <- length(get("y1"))

  # Get drop vectors
  stop1vec <- (y1 + y2 + y3) > 0
  stop2vec <- (y4 + y5 + y6) > 0
  
  drop0 <- ((stop1vec == 0) & (stop2vec == 0))  # No arms dropped
  drop1 <- ((stop1vec == 0) & (stop2vec == 1))  # Trt2 dropped
  drop2 <- ((stop1vec == 1) & (stop2vec == 0))  # Trt1 dropped
  drop3 <- ((stop1vec == 1) & (stop2vec == 1))  # Both dropped
  
  # Estimate operating charcteristics
  ae1 <- ifelse(drop0, y7, 
                ifelse(drop1, y16,
                       ifelse(drop2, rep(0, nrep),
                              rep(0, nrep))))
  
  comp1 <- mean(ifelse(drop0, y8, 
                           ifelse(drop1, y17,
                                  ifelse(drop2, rep(0, nrep),
                                         rep(0, nrep)))))
  
  ntol1 <- mean(ifelse(drop0, y9, 
                           ifelse(drop1, y18,
                                  ifelse(drop2, rep(0, nrep),
                                         rep(0, nrep)))))
  
  ae2 <- ifelse(drop0, y10, 
                ifelse(drop2, y22,
                       ifelse(drop1, rep(0, nrep),
                              rep(0, nrep)))) 
  
  comp2 <- mean(ifelse(drop0, y11, 
                           ifelse(drop2, y23,
                                  ifelse(drop1, rep(0, nrep),
                                         rep(0, nrep))))) 
  
  ntol2 <- mean(ifelse(drop0, y12, 
                           ifelse(drop2, y24,
                                  ifelse(drop1, rep(0, nrep),
                                         rep(0, nrep))))) 
  
  ae3 <- ifelse(drop0, y13, 
                ifelse(drop1, y19,
                       ifelse(drop2, y25, y28)))
  
  comp3 <- mean(ifelse(drop0, y14, 
                           ifelse(drop1, y20,
                                  ifelse(drop2, y26, y29))))
  
  ntol3 <- mean(ifelse(drop0, y15, 
                           ifelse(drop1, y21,
                                  ifelse(drop2, y27, y30))))
  
  any <- ae1 + ae2 + ae3 > 0
  
  ae1 <- mean(ae1)
  ae2 <- mean(ae2)
  ae3 <- mean(ae3)
  any <- mean(any)
  
  # Return all results in a named vector
  return(c(ns = n, stop1 = stop1, stop2 = stop2, stopboth = stopboth, ae1 = ae1, comp1 = comp1, ntol1 = ntol1, 
           ae2 = ae2, comp2 = comp2, ntol2 = ntol2, ae3 = ae3, comp3 = comp3, ntol3 = ntol3,any = any))
}

#' Tune AE decision threshold to control FWER under the global null
#'
#' This function finds the optimal decision threshold (U1) for the AE endpoint 
#' at the final analysis to ensure the family-wise error rate (FWER) does not 
#' exceed a prespecified value. This is done by iteratively updating U1 using 
#' a bisection method until the FWER constraint is satisfied.
#'
#' @param trial.outputs A numeric matrix with 30 columns:\cr
#'   - Columns 1–6: Interim posterior probabilities of inferiority (tau1).\cr
#'   - Columns 7–30: Final-stage posterior probabilities under 4 drop scenarios (tau2).
#' @param gamma Vector of thresholds for interim stopping (length 3).
#' @param FWER Desired maximum family-wise type I error rate (e.g., 0.05).
#'
#' @return A named numeric vector with:\cr
#'   - U1: Tuned AE threshold.\cr
#'   - ae1, ae2, ae3: Type I error rates for each treatment.\cr
#'   - any: Overall FWER (probability any AE test is falsely positive).
getU1 <- function(trial.outputs, gamma = c(0.2, 0.5, 0.5), FWER = 0.05) {
  
  # Flip tau2 probabilities to success scale
  z = trial.outputs
  for (j in 7:ncol(trial.outputs)){
    z[,j] <- 1 - z[,j]
  }
  
  # Interim stopping decisions (tau1) based on gamma
  for (j in 1:6) {
    assign(paste0("y", j), z[,j] > gamma[(j-1) %% 3 + 1])
  }
  
  # Compute stopping indicators
  stop1vec <- (y1 + y2 + y3) > 0
  stop2vec <- (y4 + y5 + y6) > 0
  drop0 <- !stop1vec & !stop2vec
  drop1 <- !stop1vec & stop2vec
  drop2 <- stop1vec & !stop2vec
  drop3 <- stop1vec & stop2vec

  # Initialize a lower and upper bound for U1
  lower <- 0.5
  upper <- 1
  
  nrep <- length(y1)
  
  # Binary search for optimal U1 threshold
  while((upper-lower) > 0.0002){
    
    # Current candidate value for U1
    U1 <- round(0.5*(upper+lower), 4)
    
    # Evaluate AE endpoint outcomes at candidate U1
    for (j in seq(7, 28, 3)) {
      assign(paste0("y", j), z[,j] > U1)
    }
    
    # Apply drop rules and compute false positive rates for each treatment
    ae1 <- ifelse(drop0, y7, 
                  ifelse(drop1, y16,
                         ifelse(drop2, rep(0, nrep),
                                rep(0, nrep))))
    
    ae2 <- ifelse(drop0, y10, 
                  ifelse(drop2, y22,
                         ifelse(drop1, rep(0, nrep),
                                rep(0, nrep)))) 
    
    ae3 <- ifelse(drop0, y13, 
                  ifelse(drop1, y19,
                         ifelse(drop2, y25, y28)))
    
    any <- ae1 + ae2 + ae3 > 0
    
    ae1 <- mean(ae1)
    ae2 <- mean(ae2)
    ae3 <- mean(ae3)
    any <- mean(any)
    
    # Update bounds based on FWER constraint
    if (any <= FWER){
      upper <- U1
      ae1.final <- ae1
      ae2.final <- ae2
      ae3.final <- ae3
      any.final <- any
    } else {
      lower <- U1
    }
  }
  
  # Return optimal threshold and error rates
  return(c(U1 = upper, ae1 = ae1.final, ae2 = ae2.final, ae3 = ae3.final, any = any.final))
}

#' Determine optimal interim sample size to achieve power targets
#'
#' This function identifies the smallest sample size (from a matrix of
#' operating characteristics) that achieves the desired power for all three AE comparisons 
#' at the final analysis. It assumes that `getOCsLin()` has been run over a range of 
#' sample sizes and returns power estimates for each.
#'
#' @param lin.mat A numeric matrix returned by `getOCsLin()`, with:\cr
#'   - Column 1: interim control sample sizes.\cr
#'   - Columns `inds`: power estimates for AE endpoints for Trt1, Trt2, Trt3.
#' @param inds Integer vector (length 3) specifying column indices in `lin.mat` that 
#' correspond to the power estimates for each AE comparison (default = c(5, 8, 11)).
#' @param pwr Desired power for each comparison (default = 0.8).
#'
#' @return A named numeric vector with:\cr
#'   - n: optimal interim sample size.\cr
#'   - pwr1, pwr2, pwr3: estimated power for Trt1, Trt2, and Trt3 AE comparisons.
getN <- function(lin.mat, inds = c(5, 8, 11), pwr = 0.8){
  
  # Find smallest row index where power exceeds threshold for each comparison
  row1 <- min(which(lin.mat[,inds[1]] > pwr))
  row2 <- min(which(lin.mat[,inds[2]] > pwr))
  row3 <- min(which(lin.mat[,inds[3]] > pwr))
  
  # Ensure power target was achievable in provided sample size range
  if(!is.finite(row1 + row2 + row3)){
    stop("Range of sample sizes is not large enough")
  }
  
  # Optimal sample size is the one where all comparisons achieve power
  row.all <- max(row1, row2, row3)
  
  # Return selected sample size and power values
  return(c(n = as.numeric(lin.mat[row.all, 1]), pwr1 = as.numeric(lin.mat[row.all, inds[1]]),
           pwr2 = as.numeric(lin.mat[row.all, inds[2]]), pwr3 = as.numeric(lin.mat[row.all, inds[3]])))
}
