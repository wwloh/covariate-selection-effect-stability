options(scipen = 3, digits = 2)

# rename for the covariate selection
var.list <- paste0("L.",1:p)
colnames(mydata) <- c("i",var.list,"treat","Y")
rm(rawData,L.Data)
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
# [1] "NegEmo" "Extra"  "Age"    "Consci" "Female" "Agree"  "Open"

# calculate standardized differences ==========================================
L.stable <- StdDiffEst_Ordered_DRAIPW(L.ordered$ordered,mydata)
L.stable$est
#         diff    se
# [1,] 0.0102 0.0165
# [2,] 0.0129 0.0130
# [3,] 0.0147 0.0114
# [4,] 0.0070 0.0081
# [5,] 0.0038 0.0048
# [6,] 0.0038 0.0043
# [7,] 0.0000 0.0000

L.stable$selected_orbit
# [1] 6

L.names.selected <- L.ordered$ordered[1:L.stable$selected_orbit]

# calculate treatment effect estimator given selected covariates
est <- OneDR_AIPW_Est(Ls=L.names.selected,mydata,return.se=TRUE)
round(est[1:3],2)
#  EST   SE   PV 
# 0.23 0.05 0.00
rm(est)

# convert back to original names
L.names[as.integer(lapply(sapply(
  L.names.selected, strsplit, split="L."),"[",2))]
# [1] "NegEmo" "Extra"  "Age"    "Consci" "Female" "Agree"
