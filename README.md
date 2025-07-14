# Sample Size Determination for Bayesian Platform Trials (SSTaRLeT)

This repository contains R code to implement the Sample Size Determination (SSD) method proposed in:

> **"An Efficient Approach to Design Bayesian Platform Trials with an Application to the SSTaRLeT Trial"**\
> *Luke Hagar, Lara Maleyeff et al. (2025)*

The SSD method is designed for Bayesian platform trials with features such as:

- Interim arm dropping based on posterior inferiority probabilities
- Sample size re-allocation if one or more arms are dropped
- Non-inferiority evaluation using informative robust MAP priors

---

## Required Packages

```r
install.packages(c("foreach", "doParallel", "doSNOW"))
```

---

## Files

### `hagar_maleyeff_ssd.R`

This file contains the source code implementing the SSD procedure and operating characteristic functions described in the manuscript.

| Function        | Purpose                                                                 |
|----------------|-------------------------------------------------------------------------|
| `wpost()`       | Updates posterior weights under a robust MAP prior                      |
| `betaFMsamp()`  | Samples from posterior distribution under the robust MAP prior          |
| `estimate_pp()` | Estimates posterior probability that a treatment effect exceeds a margin|
| `oneIAaddArm()` | Simulates a single trial replicate with interim and final decision points|
| `getOCsLin()`   | Approximates operating characteristics across sample sizes              |
| `decisions()`   | Computes trial outcomes for a given set of thresholds                   |
| `getU1()`       | Tunes the final AE decision threshold to control FWER                  |
| `getN()`        | Selects interim sample size to achieve target power for all comparisons |

### `example_execution.R`

This file contains a full example script that:
- Runs the sample size determination (SSD) procedure, including tuning of the final-stage threshold `U1` and selecting an optimal interim sample size
- Performs a simulation study to compare the estimated power from `getOCsLin()` against simulated power (used as a proxy for the true power)

---

## Example: Sample Size Determination Procedure

```r
library(foreach)
library(doParallel)
library(doSNOW)

# Set up parallel processing 
cores = detectCores()
cl <- makeSOCKcluster(cores[1]-1)

# Set the number of simulated trials per setting
m <- 10000
registerDoSNOW(cl)
# Create a progress bar for tracking simulations
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Sequence of interim sample sizes to evaluate -- only need 2 to use SSD
n.ints <- c(600,1000)

# Set final sample sizes to be 2.5x the interim sizes (as in manuscript)
c.final <- 2.5
n.finals <- ceiling(c.final*n.ints)

# Scenario 1: All Treatments are Clearly Acceptable
ii <- 1
for (jj in 1:length(n.ints)){
  # Simulate m = 1000 trials at each sample size
  sim.res <- foreach(k=1:m, .combine=rbind,
                     .options.snow=opts) %dopar% {
                       
                       oneIAaddArm(i, eff = c(0, 0, 0), n.int = n.ints[jj], n.final = n.finals[jj],
                                   M = 1000)
                     }
  # Save results for each sample size to a separate .csv file
  write.csv(sim.res, paste0("probs_n_int_",n.ints[jj], "_scen_",ii,".csv"), row.names = FALSE)
}

# Scenario 0: Null scenario where all treatments are Unacceptable — used to tune U1 threshold
sim.res <- foreach(k=1:m, .combine=rbind,
                   .options.snow=opts) %dopar% {
                     oneIAaddArm(i, eff = c(0.04, 0.1, 0.1), n.int = 600, n.final = 2.5*600,
                                 M = 1000)
                   }
write.csv(sim.res, paste0("probs_n_int_600_scen_0.csv"), row.names = FALSE)

# Tune final-stage AE threshold to control type I error (FWER)
# U1 = 0.975 in manuscript
temp <- read.csv("probs_n_int_600_scen_0.csv")
getU1(temp, FWER = 0.05)

# Sample size determination based on interpolated OC curves
low.mat <- read.csv("probs_n_int_600_scen_1.csv")
up.mat <- read.csv("probs_n_int_1000_scen_1.csv")

# nl, nu are parameterized in terms of the number of total observations at the interim analysis
# lb and ub are the lower and upper bounds for the total observations at the interim analysis
lin.temp <- getOCsLin(ml = low.mat, mu = up.mat, nl = 600, 
                      nu = 1000, 
                      lb = 400, ub = 1200, gamma = c(0.2, 0.5, 0.5), kappa = c(0.975, 0.99, 0.99))

# Return optimal interim sample size meeting power criteria (targeting 95\% power)
# Final sample size (using ratio in manuscript) is 2.5*n 
getN(lin.temp, pwr = 0.95)
```

---

## Plotting Estimated vs Simulated Characteristics

```r
# Simulate across additional sample sizes
n.ints <- setdiff(seq(400, 1200, 100), c(600,1000))
n.finals <- ceiling(c.final*n.ints)
ii <- 1
for (jj in 1:length(n.ints)){
  # Simulate m = 1000 trials at each sample size
  sim.res <- foreach(k=1:m, .combine=rbind,
                     .options.snow=opts) %dopar% {
                       
                       oneIAaddArm(i, eff = c(0, 0, 0), n.int = n.ints[jj], n.final = n.finals[jj],
                                   M = 1000)
                     }
  # Save results for each sample size to a separate .csv file
  write.csv(sim.res, paste0("probs_n_int_",n.ints[jj], "_scen_",ii,".csv"), row.names = FALSE)
}

# Load and evaluate simulated operating characteristics
dat.res <- NULL
for (jj in 1:length(n.ints)){
  trial.res <- read.csv(paste0("probs_n_int_",n.ints[jj], "_scen_",ii,".csv"))
  dat.res <- rbind(dat.res, decisions(trial.res, n = n.ints[jj], 
                                      gamma = c(0.2, 0.5, 0.5), kappa = c(0.975, 0.99, 0.99)))
}

# Plot estimated vs. simulated operating characteristics (main comparisons)
for (j in 2:4){
  pdf(file = paste0("Fig", j, ".pdf"),   # The directory you want to save the file in
      width = 5, # The width of the plot in inches
      height = 5)
  plot(lin.temp[,1], lin.temp[,j], type = "l", ylim = c(0,1),
       ylab = "Stopping Probability", xlab = bquote(italic(n)['int']),
       main = colnames(dat.res)[j])
  lines(n.ints, dat.res[,j], col = "red")
  legend("topright", c("Estimated", "Simulated"), col = c("black", "red"), lty = 1)
  dev.off()
}

# Plot remaining operating characteristics
for (j in 5:14){
  pdf(file = paste0("Fig", j, ".pdf"),   # The directory you want to save the file in
      width = 5, # The width of the plot in inches
      height = 5)
  plot(lin.temp[,1], lin.temp[,j], type = "l", ylim = c(0,1),
       ylab = "Stopping Probability", xlab = bquote(italic(n)['int']),
       main = colnames(dat.res)[j])
  lines(n.ints, dat.res[,j], col = "red")
  legend("bottomright", c("Estimated", "Simulated"), col = c("black", "red"), lty = 1)
  dev.off()
}

```

---

## Citation

If you use this code, please cite:

**Hagar, L., Maleyeff, L., et al. (2025)**\
*"An Efficient Approach to Design Bayesian Platform Trials with an Application to the SSTaRLeT Trial"*\
X*.

---

## Authors

- Luke Hagar, [luke.hagar@mail.mcgill.ca](mailto:luke.hagar@mail.mcgill.ca).
- Lara Maleyeff, [lara.maleyeff@mail.mcgill.ca](mailto:lara.maleyeff@mail.mcgill.ca).
- Shirin Golchi
- Dick Menzies


