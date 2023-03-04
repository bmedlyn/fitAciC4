# Code fits A-Ci function for C4 photosynthesis

# To do: 
# Check the ci transition point
# add SE's
# test if Vpmax and Vcmax are identifiable
# implement an empirical model as a benchmark

source("R/loadPackages.R") # required packages
source("R/AciC4.R")        # implements C4 model
source("R/functions.R")    # all fitting functions

# Example dataset from Vinod Jacob (Hawkesbury Institute for the Environment)
# Note example is for kangaroo grass, photosynthetic rates low
vin <- read_xlsx("Data/Vin_themeda_Aci_curves_aug22.xlsx")
# Set column headings
vin <- renameLi6800(vin)
# Check dataset looks OK. Each curve is identified with individual ID
with(vin,plot(Ci,Photo,col=as.factor(ID)))
# Filter out incorrect points
vin <- subset(vin, Ci > 0)

# split dataset into individual curves
vinlist <- split(vin,vin$ID)
# function do_the_lot fits the curve, finds the Ci transition point,
# and visualises the fit
# apply function to first curve to test
do_the_lot(vinlist[[1]])

# make pdf file with all curves
pdf(file="Output/vin_out.pdf")
vinout <- lapply(vinlist,do_the_lot)
dev.off()

# Put output in df
vin_smry <- bind_rows(vinout)

# investigate fitted parameters
with(vin_smry,plot(Vcmax,Vpmax))
with(vin_smry,plot(Vcmax,Jmax))
with(vin_smry,plot(ci_trans,Vpmax))
with(vin_smry,plot(RMSE,Vcmax))
with(vin_smry,plot(RMSE,Jmax/Vcmax))
with(vin_smry,plot(RMSE,Vpmax/Vcmax))
