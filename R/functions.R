
### LIST OF FUNCTIONS
# 
# renameLi6800 - give standard column names to LI6800 outputs
# fAC4v - Wrapper to calculate Av from function AciC4
# fAC4j - Wrapper to calculate Aj from function AciC4
# fitAciC4num_J - Fit Jmax to light-limited portion of curve
# fitAciC4num_V - Fit Vcmax and Vpmax to enzyme-limited portion of curve
#   Currently assuming RD0 = 0.01 Vcmax, other parameters defaults
# empirical - Non-rectangular hyperbola benchmark
# fit_empirical - Fits the empirical benchmark
#
# fitAciC4trans - Fit full A-Ci curve. This works by looping over all
#   possible breakpoints of the curve, fitting the two halves of the curve
#   each time, then finding the breakpoint that gives the lowest SSQ
# citrans - Find Ci value at which Aj = Ac
# visfit - make plot
# do_the_lot - call fitting and plotting routines for one curve
###

# wrappers for C4 function, so it only returns A, to which we are fitting
# This allows us to fit numerical solution
fAC4v <- function(Vcmax,Vpmax25,Jmax25,PPFD,Ci,Tleaf,RdRatio) {
  Rd0 <- RdRatio * Vcmax
  AciC4(Vcmax=Vcmax,VPMAX25=Vpmax25,JMAX25=Jmax25,RD0=Rd0,
        PPFD=PPFD,Ci=Ci,Tleaf=Tleaf)$Ac
}

fAC4j <- function(Vcmax,Vpmax25,Jmax25,PPFD,Ci,Tleaf,Rd0) {
  AciC4(Vcmax=Vcmax,VPMAX25=Vpmax25,JMAX25=Jmax25,RD0=Rd0,
        PPFD=PPFD,Ci=Ci,Tleaf=Tleaf)$Aj
}

# fit Jmax to light-limited portion of curve; given Rd0
fitAciC4num_J <- function(data, Vcmax, Vpmax, Rd0) {
  
  # Tleaf is not vectorised - just use mean
  meanTleaf = mean(data$Tleaf)
  # call nls to fit
  fit <- try(nls(Photo ~ fAC4j(Vcmax=Vcmax, Vpmax25=Vpmax, Jmax25=JMAX25, 
                               PPFD=PARi, Ci=Ci, Tleaf=meanTleaf, Rd0=Rd0),
             start=list(JMAX25=100),trace=FALSE,
             control=(nls.control(warnOnly=TRUE)),data = data))
  coefs <- coef(summary(fit))
  ret <- coefs[1:2]
  names(ret) <- c("Jmax","JmaxSE")
  return(ret)
}

# fit Vcmax and Vpmax to light-limited portion of curve
fitAciC4num_V <- function(data,RdRatio) {
  
  # Tleaf is not vectorised - just use mean
  meanTleaf = mean(data$Tleaf)
  # call nls to fit
  fit <- try(nls(Photo ~ fAC4v(Vcmax, Vpmax25, Jmax25=100, 
                               PPFD=PARi, Ci=Ci, Tleaf=meanTleaf, RdRatio=RdRatio),
                 start=list(Vcmax=25, Vpmax25 = 20),trace=FALSE,
                 control=(nls.control(warnOnly=TRUE)),data = data))
  coefs <- coef(summary(fit))
  ret <- coefs[1:4]
  names(ret) <- c("Vcmax","Vpmax","VcmaxSE","VpmaxSE")
  return(ret)
}

# change names from 6800 to what's needed
renameLi6800 <- function(data) {
  
  data <- renameCol(data,src = c("A","gsw","VPDleaf","Qin"),
                    tgt=c("Photo","Cond","VpdL","PARi"))
  return(data)
}

# Fit ACiC4 by fitting the curve in two parts
# Enzyme to low Ci and Light limitation to high Ci
# Possible transition points lie from between points 3,4 up to between points n-2,n-1 
# Find the transition point that minimises the sum of squares overall
# Fixing an RdRatio
fitAciC4trans <- function(data,RdRatio = 0.01) {
  
  # First sort data by Ci
  data <- subset(data, Ci > 0)
  data <- data[order(data$Ci),]
  meanTleaf <- mean(data$Tleaf)
  
  # Find length of curve
  len <- length(data$Ci)
  # Need at least 5 points for this to work
  if (len < 5) return
  
  # Set up ssq array
  ssq <- jmax <- vcmax <- vpmax <- Rd0 <- jmaxSE <- vcmaxSE <- vpmaxSE <- c()
  
  # loop over possible transition points
  for (i in 1:(len-4)) {
    # cut dataset in 2
    d1 <- data[1:(i+2),]
    d2 <- data[(i+3):len,]
    
    # fit Vcmax and Vpmax to lower portion of curve
    vs <- fitAciC4num_V(d1,RdRatio)
    vcmax[i] <- vs["Vcmax"]
    vpmax[i] <- vs["Vpmax"]
    vcmaxSE[i] <- vs["VcmaxSE"]
    vpmaxSE[i] <- vs["VpmaxSE"]
    Rd0[i] <- RdRatio*vcmax[i]
    
    # fit Jmax to upper portion of curve, given Rd0
    js <- fitAciC4num_J(d2,vcmax[i],vpmax[i],Rd0[i])
    jmax[i] <- js["Jmax"]
    jmaxSE[i] <- js["JmaxSE"]
    
    # calculate ssq
    fitted <- with(data,AciC4(Vcmax=vcmax[i],VPMAX25=vpmax[i],
                              JMAX25=jmax[i],
                              PPFD=PARi,Ci=Ci,Tleaf=meanTleaf,RD0=Rd0[i]))
    ssq[i] <- sum((data$Photo-fitted$ALEAF)^2)
    
  }
  
  # find minimum ssq
  i_trans <- which.min(ssq)
  rmse <- sqrt(ssq[i_trans]/length(data$Ci))
  
  # best-fit parameters
  pars <- data.frame(vcmax[i_trans],vpmax[i_trans],jmax[i_trans],Rd0[i_trans],
                     vcmaxSE[i_trans],vpmaxSE[i_trans],jmaxSE[i_trans],
                     rmse,i_trans+2)
  names(pars) <- c("Vcmax","Vpmax","Jmax","Rd0","VcmaxSE","VpmaxSE","JmaxSE",
                   "RMSE","trans_pt")
  return(pars)
  
}

# visualise fit
visfit <- function(data,pars) {
  
  # calculate fitted values
  meanTleaf <- mean(data$Tleaf)
  meanPAR <- mean(data$PARi)
  Ciseq <- seq(20,max(data$Ci),by=20)
  fittedAnum <- with(data,AciC4(Vcmax=pars$Vcmax,VPMAX25=pars$Vpmax,
                                JMAX25=pars$Jmax,
                                PPFD=meanPAR,Ci=Ciseq,Tleaf=meanTleaf,RD0=pars$Rd0))
  
  cit <- ifelse(!is.numeric(pars$ci_trans),NA,round(pars$ci_trans))
  title <- paste0(data$ID[1]," Vc ",round(pars$Vcmax),
                  " Vp ",round(pars$Vpmax),
                  " J ",round(pars$Jmax),
                  " RMSE ",round(pars$RMSE,1),
                  " pt ",pars$trans_pt,
                  " Cit ",cit)
  #make plot
  with(data,plot(Ci, Photo,ylim=c(0,max(Photo)*2), main=title))
  with(fittedAnum,points(Ci, Aj, col="red",type="l"))
  with(fittedAnum,points(Ci, Ac, col="blue",type="l"))
  points(fittedAnum$Ci,fittedAnum$ALEAF,col="black",type="l")
  abline(v=cit,lty=3)
}

# Function to calculate ci transition point
# Find the Ci where Ac = Aj
citrans <- function(data,pars) {
  
  fci <- function(ci) {
      a <- AciC4(Vcmax=Vcmax,VPMAX25=Vpmax25,JMAX25=Jmax25,
                 PPFD=PPFD,Ci=ci,Tleaf=Tleaf,RD0=Rd0)
      return(a$Aj - a$Ac)
  }
  
  Vcmax <- pars$Vcmax
  Vpmax25 <- pars$Vpmax
  Jmax25 <- pars$Jmax
  Rd0 <- pars$Rd0
  PPFD <- mean(data$PARi)
  Tleaf <- mean(data$Tleaf)
  
  ci_trans <- try(uniroot(fci,interval=c(5,1500))$root)
  if (!is.numeric(ci_trans)) ci_trans <- NA
  return(ci_trans)

}

# Empirical benchmark - try non-rectangular hyperbola
# with four parameters
# theta*A^2 - (alpha*Ci + Amax)*A + alpha*Ci*Amax - Rd = 0
non_rect_hyp <- function(Ci,Amax,alpha,theta,Rd) {
  
  a <- theta
  b <- -(alpha*Ci + Amax)
  c <- alpha*Ci*Amax
  disc <- b^2 - 4*a*c
  ret <- (-b - sqrt(disc)) / (2*a) + Rd
  return(ret)
  
}

# Fit benchmark
# may not be possible to fit all four parameters
fit_empirical <- function(data) {
  
  # fit with nls
  fit <- try(nls(Photo ~ non_rect_hyp(Ci=Ci,Amax,alpha,theta=0.7,Rd),
                 start=list(Amax=max(data$Photo), alpha=0.5,
                            Rd=0.5),trace=FALSE,
                 control=(nls.control(warnOnly=TRUE)),data = data))
  coefs <- coef(fit)
  fittedAemp <- with(data,non_rect_hyp(data$Ci,Amax=coefs[1],
                                    alpha=coefs[2],theta=0.7,Rd=coefs[3]))
  
  empssq <- sum((data$Photo-fittedAemp)^2)
  emp_RMSE <- sqrt(empssq/length(fittedAemp))
  title <- paste0("RMSE ",round(emp_RMSE,1), "Amax ",round(coefs[1]),
                  " alpha ",round(coefs[2],2)," Rd ",round(coefs[3],2))
  with(data,plot(Ci,Photo,ylim=c(0,2*max(Photo)),main=title))
  points(data$Ci,fittedAemp,col="red",type="l")

  return(emp_RMSE)
}

# Function to do the lot
# Fit curve, make plot, return params and ci transition point
do_the_lot <- function(data,RdRatio=0.01) {
  
  # fit curve
  pars <- fitAciC4trans(data,RdRatio)
  pars <- check_upper(data,pars,RdRatio)
  if (!is.null(pars)) {
    pars$ci_trans <- citrans(data,pars)
    visfit(data,pars)
  }
  pars$emp_RMSE <- fit_empirical(data)
  pars_a <-fitAciC4_ainsy(data)
  pars$AmaxLA <- pars_a["Amax"]
  pars$VpmaxLA <- pars_a["Vpmax"]
  return(pars)
  
}

