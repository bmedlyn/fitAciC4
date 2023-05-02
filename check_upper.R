# Adding extra function to allow for possibility that Jmax is never limiting
# Check two things: does having a Jmax actually give you a better RMSE? 
# Does the fitted Jmax give Aj < Ac at the upper end of the curve? 
# Only take Jmax if both things are true

check_upper <- function(data, pars, RdRatio = 0.01) {
  
  # First sort data by Ci
  data <- subset(data, Ci > 0)
  data <- data[order(data$Ci),]
  meanTleaf <- mean(data$Tleaf)
  Jmax_upper <- 10000
    
  # fit Vcmax and Vpmax to full curve
  vs <- fitAciC4num_V(data,RdRatio)

  # calculate ssq with only Vcmax and Vpmax fitted
  fitted_U <- with(data,AciC4(Vcmax=as.numeric(vs["Vcmax"]),
                              VPMAX25=as.numeric(vs["Vpmax"]),
                              JMAX25=Jmax_upper,
                              PPFD=PARi,Ci=Ci,Tleaf=meanTleaf,
                              RD0=RdRatio*as.numeric(vs["Vcmax"])))
  ssq_U <- sum((data$Photo-fitted_U$ALEAF)^2)
  
  # calculate the fitted values with best parameters previously
  fitted_L <- with(data,AciC4(Vcmax=pars$Vcmax,
                              VPMAX25=pars$Vpmax,
                              JMAX25=pars$Jmax,
                              PPFD=PARi,Ci=Ci,Tleaf=meanTleaf,
                              RD0=pars$Rd0))
  ssq_L <- sum((data$Photo-fitted_L$ALEAF)^2)
  
  # if either SSQ is better OR upper A is Vcmax limited then take the new fit
  len <- length(fitted_L$ALEAF)
  if ( (ssq_U < ssq_L) | (fitted_L$Aj[len] > fitted_L$Ac[len])) {
  
  # new best-fit parameters
    rmse <- sqrt(ssq_U/length(data$Ci))
    pars <- data.frame(vs["Vcmax"],vs["Vpmax"],NA,RdRatio*vs["Vcmax"],
                     vs["VcmaxSE"],vs["VpmaxSE"],NA,
                     rmse,NA)
    names(pars) <- c("Vcmax","Vpmax","Jmax","Rd0","VcmaxSE","VpmaxSE","JmaxSE",
                   "RMSE","trans_pt")
  }
  return(pars)
  
}