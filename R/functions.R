
### LIST OF FUNCTIONS
# 
# renameLi6800 - give standard column names to LI6800 outputs
# fAC4v - Wrapper to calculate Av from function AciC4
# fAC4j - Wrapper to calculate Aj from function AciC4
# fitAciC4num_J - Fit Jmax to light-limited portion of curve
# fitAciC4num_V - Fit Vcmax and Vpmax to enzyme-limited portion of curve
#   Currently assuming RD0 = 0.01 Vcmax, other parameters defaults
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
fAC4v <- function(Vcmax,Vpmax25,Jmax25,PPFD,Ci,Tleaf,RdRatio=0.01) {
  Rd0 <- RdRatio * Vcmax
  AciC4(Vcmax=Vcmax,VPMAX25=Vpmax25,JMAX25=Jmax25,RD0=Rd0,
        PPFD=PPFD,Ci=Ci,Tleaf=Tleaf)$Ac
}

fAC4j <- function(Vcmax,Vpmax25,Jmax25,PPFD,Ci,Tleaf,RdRatio=0.01) {
  Rd0 <- RdRatio * Vcmax
  AciC4(Vcmax=Vcmax,VPMAX25=Vpmax25,JMAX25=Jmax25,RD0=Rd0,
        PPFD=PPFD,Ci=Ci,Tleaf=Tleaf)$Aj
}

# fit Jmax to light-limited portion of curve
fitAciC4num_J <- function(data) {
  
  # Tleaf is not vectorised - just use mean
  meanTleaf = mean(data$Tleaf)
  # call nls to fit
  ret <- try(nls(Photo ~ fAC4j(Vcmax=20, Vpmax25=20, Jmax25=JMAX25, PPFD=PARi, Ci=Ci, Tleaf=meanTleaf),
             start=list(JMAX25=100),trace=FALSE,
             control=(nls.control(warnOnly=TRUE)),data = data))
  coefs <- coef(ret)
  return(coefs)
}

# fit Vcmax and Vpmax to light-limited portion of curve
fitAciC4num_V <- function(data) {
  
  # Tleaf is not vectorised - just use mean
  meanTleaf = mean(data$Tleaf)
  # call nls to fit
  ret <- try(nls(Photo ~ fAC4v(Vcmax, Vpmax25, Jmax25=100, PPFD=PARi, Ci=Ci, Tleaf=meanTleaf),
                 start=list(Vcmax=25, Vpmax25 = 20),trace=FALSE,
                 control=(nls.control(warnOnly=TRUE)),data = data))
  coefs <- coef(ret)
  return(coefs)
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
fitAciC4trans <- function(data) {
  
  # First sort data by Ci
  data <- subset(data, Ci > 0)
  data <- data[order(data$Ci),]
  meanTleaf <- mean(data$Tleaf)
  
  # Find length of curve
  len <- length(data$Ci)
  # Need at least 5 points for this to work
  if (len < 5) return
  
  # Set up ssq array
  ssq <- jmax <- vcmax <- vpmax <- c()
  
  # loop over possible transition points
  for (i in 1:(len-4)) {
    # cut dataset in 2
    d1 <- data[1:(i+2),]
    d2 <- data[(i+3):len,]
    
    # fit Vcmax and Vpmax to lower portion of curve
    vs <- fitAciC4num_V(d1)
    vcmax[i] <- vs[1]
    vpmax[i] <- vs[2]
    
    # fit Jmax to upper portion of curve
    js <- fitAciC4num_J(d2)
    jmax[i] <- js[1]
    
    # calculate ssq
    fitted <- with(data,AciC4(Vcmax=vcmax[i],VPMAX25=vpmax[i],
                              JMAX25=jmax[i],
                              PPFD=PARi,Ci=Ci,Tleaf=meanTleaf))
    ssq[i] <- sum((data$Photo-fitted$ALEAF)^2)
    
  }
  
  # find minimum ssq
  i_trans <- which.min(ssq)
  rmse <- sqrt(ssq[i_trans]/length(data$Ci))
  
  # best-fit parameters
  pars <- data.frame(vcmax[i_trans],vpmax[i_trans],jmax[i_trans],
                     rmse,i_trans)
  names(pars) <- c("Vcmax","Vpmax","Jmax","RMSE","trans_pt")
  return(pars)
  
}

# visualise fit
visfit <- function(data,pars) {
  
  # calculate fitted values
  meanTleaf <- mean(data$Tleaf)
  fittedAnum <- with(data,AciC4(Vcmax=pars$Vcmax,VPMAX25=pars$Vpmax,
                                JMAX25=pars$Jmax,
                                PPFD=PARi,Ci=Ci,Tleaf=meanTleaf))
  fittedAnum <- fittedAnum[order(fittedAnum$Ci),]
  
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
}

# Function to calculate ci transition point
# Find the Ci where Ac = Aj
citrans <- function(data,pars) {
  
  fci <- function(ci) {
      a <- AciC4(Vcmax=Vcmax,VPMAX25=Vpmax25,JMAX25=Jmax25,PPFD=PPFD,Ci=ci,Tleaf=Tleaf)
      return(a$Aj - a$Ac)
  }
  
  Vcmax <- pars$Vcmax
  Vpmax25 <- pars$Vpmax
  Jmax25 <- pars$Jmax
  PPFD <- mean(data$PARi)
  Tleaf <- mean(data$Tleaf)
  
  ci_trans <- try(uniroot(fci,interval=c(5,1500))$root)
  if (!is.numeric(ci_trans)) ci_trans <- NA
  return(ci_trans)

}

# Function to do the lot
# Fit curve, make plot, return params and ci transition point
do_the_lot <- function(data) {
  
  # fit curve
  pars <- fitAciC4trans(data)
  if (!is.null(pars)) {
    pars$ci_trans <- citrans(data,pars)
    visfit(data,pars)
  }
  return(pars)
  
}

