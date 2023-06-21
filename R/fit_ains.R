
# implement fitting routine according to Ainsworth and Leakey
# this fits Vpmax to initial slope (Ci < 150 ppm) - assumes sufficient points
# Vcmax is fitted to the whole curve

# Vpmax:  f = ((x*Vpmax)/(x+Kp))-Rm

# Vmax: f = R+(a*x+Amax-sqrt((a*x+Amax)^2-4*O*a*x*Amax))/(2*O) 
# (Dynamic Fit Options: 200 Total Number of Fits; 1000 Max Number of Iterations)

# initial slope function
vpf <- function(Ci,Vpmax,Kp,Rm) {
  Ci*Vpmax/(Ci+Kp)-Rm
}

# fit overall function
fitAciC4_ainsy <- function(data,lowCi = 150) {
  
  # fit Vpmax and Rm to initial slope
  data <- data[order(data$Ci),]
  low <- subset(data,Ci < lowCi)
  Kp <- 80
  fit_low <- try(nls(Photo ~ vpf(Ci=Ci,Vpmax,Kp,Rm),
                 start=list(Vpmax=25,Rm=0.5),trace=FALSE,
                 control=(nls.control(warnOnly=TRUE)),data = low))
  
  # fit 4 parameters to full curve
  fit_all <- try(nls(Photo ~ non_rect_hyp(Ci=Ci,Amax,alpha,theta=0.7,Rd),
                 start=list(Amax=max(data$Photo), alpha=0.1,
                            Rd=0.5),trace=FALSE,
                 control=(nls.control(warnOnly=TRUE)),data = data))
  ret1 <- coef(fit_low)
  ret2 <- coef(fit_all)
  ret <- c(ret1,ret2)
  
  #Ciseq <- seq(20,max(data$Ci),by=20)
  Ciseq <- seq(-100,1500,by=20)
  fittedAvp <- with(data, vpf(Ci=Ciseq,Vpmax=ret["Vpmax"],Kp,Rm=ret["Rm"]))
  fittedAvc <- with(data,non_rect_hyp(Ciseq,Amax=ret["Amax"],
                                       alpha=ret["alpha"],theta=0.7,Rd=ret["Rd"]))
  title <- paste0("Ains ",data$ID[1],
                  " Vpmax ",round(ret["Vpmax"],1)," Amax ",round(ret["Amax"],1))
  with(data,plot(Ci,Photo,ylim=c(-10,1.2*max(Photo)),main=title,
                 xlim=c(-100,1500)))
  points(Ciseq,fittedAvc,col="red",type="l")
  points(Ciseq,fittedAvp,col="blue",type="l") 
  abline(v=0,lty=2)
  abline(v=150,lty=2)
  return(ret)
  
}