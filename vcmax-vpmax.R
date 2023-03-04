
# Script exploring effect of parameters Vcmax and Vpmax on A.enzyme
# Suggests that Vpmax is fittable but Vcmax may not be if Jmax is low
# Consistent with SE being low for Vpmax and Jmax but high for Vcmax
# Would a better way be to assume that Jmax and Vcmax co-limit? 

Aenzyme <- function(Ci,Vcmax25,Vpmax25) {

  # DEFAULTS
  PPFD=1500
  Tleaf = 25
  Vpr=80      
  alpha=0.0		  
  gbs=3e-3	 
  O2=210		 
  x=0.4		 
  THETA=0.7
  Q10 = 2.3
  RD0=1
  RTEMP=25
  TBELOW=0
  DAYRESP=1
  Q10F=2
  FRM=0.5

  TK <- Tleaf+273.15
  .Rgas <- function()8.314
  
  # Temperature effects on Vcmax, Vpmax and Jmax (Massad et al. 2007)
  # This function returns value between 0 and 1.
  Arrhenius <- function(TK, Ea, Hd, DS){
    exp( Ea*(TK-298)/(298*.Rgas()*TK) ) * 
      (1+exp( (298*DS-Hd)/(298*.Rgas()) )) / 
      (1+exp( (TK*DS-Hd)/(TK*.Rgas()) )) 
  }
  
  # Half the reciprocal for Rubisco specificity (NOT CO2 compensation point)
  low_gammastar <- 1.93e-4
  
  # Michaelis-Menten coefficients for CO2 (Kc, mu mol mol-1) and 
  # O (Ko, mmol mol-1) and combined (K)
  Kc <- 650*Q10^((Tleaf-25)/10)
  Kp <- 80*Q10^((Tleaf-25)/10)
  Ko <- 450*Q10^((Tleaf-25)/10)
  K <- Kc*(1+O2/Ko)
  
  # T effects according to Massad et al. (2007)
  Vcmax <- Vcmax25*Arrhenius(TK, 67294, 144568, 472)
  Vpmax <- Vpmax25*Arrhenius(TK, 70373, 117910, 376)
  
  # Day leaf respiration, umol m-2 s-1
  Rd <- RD0 * Q10^((Tleaf-RTEMP)/10) * DAYRESP
  Rm <-  FRM*Rd
  
# PEP carboxylation rate
  Vp <- pmin(Ci*Vpmax/(Ci+Kp),Vpr)

# Quadratic solution for enzyme limited C4 assimilation
  a.c <- 1 - (alpha*Kc)/(0.047*Ko)
  b.c <- -( (Vp-Rm+gbs*Ci) + (Vcmax-Rd) + gbs*K + 
            alpha*low_gammastar/0.047*( low_gammastar*Vcmax+Rd*Kc/Ko ) )
  c.c <- (Vcmax-Rd)*(Vp-Rm+gbs*Ci) - (Vcmax*gbs*low_gammastar*O2 + Rd*gbs*K)

  A.enzyme <- (-b.c - sqrt(b.c^2 - 4*a.c*c.c)) / (2*a.c)

  return(A.enzyme)
}

Ci = seq(20,500,by=20)
Vcmax25 <- seq(10,100,by=10)
Vpmax25 <- seq(10,100,by=10)

# Effect of varying Vcmax - sets asymptote 
# will not be fitted when Jmax is strongly limiting
plot(Ci,Aenzyme(Ci,Vcmax25=50,Vpmax25=50),col="red",pch=19)
for (i in 1:length(Vcmax25)) {
  points(Ci,Aenzyme(Ci,Vcmax25[i],Vpmax25=50),type="l")
}
# Effect of varying Vpmax - sets initial slope
for (i in 1:length(Vpmax25)) {
  points(Ci,Aenzyme(Ci,Vcmax25=50,Vpmax25=Vpmax25[i]),type="l")
}
# Varying both at same time
for (i in 1:length(Vcmax25)) {
  points(Ci,Aenzyme(Ci,Vcmax25=Vcmax25[11-i],Vpmax25=Vpmax25[i]),type="l")
}


