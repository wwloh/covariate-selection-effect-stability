options(scipen = 2, digits = 3)

# rename for the covariate selection
var.list <- paste0("L.",1:p)
colnames(mydata) <- c("i",var.list,"treat","Y")
rm(rawData)
summary(mydata)

# load helper functions
source("../selection-stability.R")

# order covariates based on priority to be confounders ====================== 
L.ordered <- ForwardSelect_DS(Y=mydata$Y,
                              X=mydata[,var.list],
                              A=mydata$treat)

# convert back to original names
L.names.ordered <- L.names[as.integer(lapply(sapply(
  L.ordered$ordered, strsplit, split="L."),"[",2))]
L.names.ordered
# [1] "Pretest"  "GradeExp" "AfAm"     "Cauc"     "AikPos"   "Grade"   
# [7] "CRAlg"    "AikNeg"   "CRMeas"   "CRSpace"  "Other"    "AmInd"   
# [13] "PR"       "MexAm"    "CRCD"     "AsPIFi"   "CRNum"    "Sex" 

# calculate standardized differences ==========================================
L.stable <- StdDiffEst_Ordered_DRAIPW(L.ordered$ordered,mydata)
round(L.stable$est,2)
#       diff   se
# [1,] -1.90 1.79
# [2,] -1.97 1.72
# [3,] -1.86 1.71
# [4,] -1.66 1.66
# [5,] -1.64 1.63
# [6,] -1.83 1.50
# [7,] -2.44 1.59
# [8,] -2.06 1.37
# [9,] -1.47 1.21
# [10,] -1.40 1.20
# [11,] -0.42 0.27
# [12,] -0.11 0.27
# [13,] -0.26 0.22
# [14,] -0.31 0.21
# [15,] -0.22 0.17
# [16,] -0.12 0.12
# [17,]  0.12 0.11
# [18,]  0.00 0.00

L.stable$selected_orbit
# [1] 2

L.names.selected <- L.ordered$ordered[1:L.stable$selected_orbit]

# calculate treatment effect estimator given selected covariates
est <- OneDR_AIPW_Est(Ls=L.names.selected,mydata,return.se=TRUE)
round(est[1:3],2)
#   EST    SE    PV 
# -4.32  1.04  0.00
rm(est)

# convert back to original names
L.names[as.integer(lapply(sapply(
  L.names.selected, strsplit, split="L."),"[",2))]
# [1] "Pretest"  "GradeExp"
