#' ICC to Type I Error Inflation
#'
#' This function allows you to calculate the inflated Type-I error rate for a particular sample size and ICC value.
#' @param lvl_2_units The number of level-2 units.
#' @param lvl_1_units The number of level-1 units.
#' @param target_icc The ICC level. 
#' @param replications The number of simulations to be run. Defaults to 1000. 
#' @param icc_acc The accuracy for the simulated ICC. Defaults to +- .01.
#' @keywords icc
#' @keywords type i error
#' @export
#' @examples
#' icc_type_i(lvl_2_units = 20, lvl_1_units = 10, target_icc = .1, replications = 1000, icc_acc = .01)

# create function
icc_type_i <- function(lvl_2_units, lvl_1_units, target_icc, replications = 1000, icc_acc = .01) {
  
  # ensure all arguments are numeric
  if (!is.numeric(lvl_2_units)) {
    stop("Argument level_2_units is not numeric.")
  } else if (!is.numeric(lvl_1_units)) {
    stop("Argument level_1_units is not numeric.")
  } else if (!is.numeric(target_icc)) {
    stop("Argument target_icc is not numeric.")
  } else if (!is.numeric(replications)) {
    stop("Argument replications is not numeric.")
  } else if (!is.numeric(icc_acc)) {
    stop("Argument icc_acc is not numeric.")
  }
  
  # ensure target_icc is between 0 and 1
  if (target_icc > 1 | target_icc < 0) {
    stop("Argument target_icc must be between 0 and 1.")
  }
  
  # calculate total number of units and number of units per condition
  total_units          <- lvl_2_units * lvl_1_units # total number of units accross all levels
  subjects_cond        <- total_units / 2 # the number of subjects per condition; subjects arbitrarly separated into a treatment and a control condition
  
  # set target icc
  icc_inf_bound <- target_icc - icc_acc # set the lowest acceptable simulated ICC level
  icc_sup_bound <- target_icc + icc_acc # set the highest acceptable simulated ICC level
  
  # set up model values
  sigma_a       <- (-66579.2859 * target_icc^10) + 
                   (278665.7676 * target_icc^9)  - 
                   (501224.1950 * target_icc^8)  + 
                   (506823.0658 * target_icc^7)  - 
                   (316468.9923 * target_icc^6)  + 
                   (126275.3718 * target_icc^5)  - 
                   (32250.17988 * target_icc^4)  + 
                   (5162.324206 * target_icc^3)  - 
                   (498.8644784 * target_icc^2)  + 
                   (35.70130764 * target_icc)    + 
                   (.1601755428) # sigma_a for a given icc
  mu_a          <- 0 # mean for intercept
  mu_b          <- 3 # mean for the slope
  sigma_b       <- 4 # get sd for slope
  sigma_y       <- 1 # get sd for outcome variable
  rho           <- 0.56 # between group correlation parameter

  # set counter values
  counter <- 0
  
  # preallocate memory for output
  zr <- data.frame(count           = rep(NA, replications),
                   t_val           = rep(NA, replications),
                   p_val           = rep(NA, replications),
                   icc             = rep(NA, replications),
                   kish_correction = rep(NA, replications),
                   kish_t_val      = rep(NA, replications),
                   kish_p_val      = rep(NA, replications))
  
  # set group information 
  group        <- rep(1:lvl_2_units, each = lvl_1_units) # create a group variable with each level level 2 variable repeated the number of level 1 units 
  cond         <- gl(n = 2,           k = subjects_cond) # create a condition faactor with two levels, each repeated the number of subjects per condition 
  lvl_2_groups <- gl(n = lvl_2_units, k = lvl_1_units) # create a factor with levels equal to the number of level 2 units and with each level repeated equal to the number of level 1 units
  
  # get variance and covariances
  sigma_ab <- array(c(sigma_a^2, # variance for the intercept
                      rho * sigma_a * sigma_b, # covariance between the intercept and the slope
                      rho * sigma_a * sigma_b, # covariance between the intercept and the slope
                      sigma_b^2), # variance for the slope
                    dim = c(2, 2)) # number of dimensions for the matrixmv
  
  # simulate data
  while (counter != replications) { 
    
    # simulate from a multivariate normal distribution
    ab <- MASS::mvrnorm(n     = lvl_2_units,   # draw the number of samples equivalent to level 1 units
                  mu    = c(mu_a, mu_b), # set means for the variables
                  Sigma = sigma_ab)      # set the covariance matrix
    a <- ab[ , 1] # simulation of intercept
    b <- ab[ , 2] # simulation of slope
    x <- rnorm(n = total_units) # sample the number of total units from a normal distribution
    
    # sample from a random normal distribution
    data <- rnorm(n    = total_units, # sample the total number of units
                  mean = a[group] + b * x, # set the mean to expected group mean 
                  sd   = sigma_y) # set the standard deviation to the standard deviation of the outcome variable
    
    # predict simulated data from condition
    lm_fit        <- lm(data ~ cond) # predict simulated data from condition
    mod_sum       <- summary(lm_fit) # summarize model
    cond_estimate <- mod_sum$coefficients[2] # extract condition slope
    cond_se       <- mod_sum$coefficients[4] # extract condition SE
    t_val         <- mod_sum$coefficients[6] # extract condition t
    p_val         <- mod_sum$coefficients[8] # extract condition p
    
    # predict simulated data from groups
    lm_lvl_2_groups_fit <- anova(lm(data ~ lvl_2_groups))
    
    # calculate icc
    mse_within  <- lm_lvl_2_groups_fit[2,2] / lm_lvl_2_groups_fit[2,1] # calculate residual mean squares error
    mse_groups  <- lm_lvl_2_groups_fit[1,2] / lm_lvl_2_groups_fit[1,1] # calculate group mean squares error
    mse_between <- (mse_groups - mse_within) / lvl_1_units # calculate between mean squares error
    icc         <- mse_between / (mse_within + mse_between) # calcualte icc
    
    # skip loop if icc isn't between the values of interest
    if (icc < icc_inf_bound || icc > icc_sup_bound) { 
      next
    }
    
    # calculate t-test using kish's correction
    kish_correction <- sqrt(1 + (lvl_1_units - 1) * icc) # calculate correction value for the given ICC
    kish_t_val      <- t_val / kish_correction # reduce t_val by correction value
    kish_p_val      <- 2 * pt(-abs(kish_t_val), anova(lm_fit)[2,1]) # calculate p-val for kish corrected degrees of freedom
    
    # increment counter by 1
    counter <- counter + 1 
    
    # assign calculated values to zr
    zr[counter,           "count"] <- counter
    zr[counter,           "t_val"] <- t_val
    zr[counter,           "p_val"] <- p_val
    zr[counter,             "icc"] <- icc
    zr[counter, "kish_correction"] <- kish_correction
    zr[counter,      "kish_t_val"] <- kish_t_val
    zr[counter,      "kish_p_val"] <- kish_p_val
    
    # notify of ICC saved
    if (counter %% (100) == 0) {
      message(paste0("Replication ", 
                     counter, 
                     " saved."))
    }
  }
  
  # calculate false positive rate
  type_1_error <- sum(zr$p_val < .05) / length(zr$p_val) * 100
  
  # calculate mean icc
  mean_icc <- mean(zr$icc)
  
  # provide the average ICC for the simulated data
  message(paste0("The average ICC for the simulated data was ", round(mean_icc, 2), "."))
  
  # return type_1_error
  return(type_1_error)
  
}
  