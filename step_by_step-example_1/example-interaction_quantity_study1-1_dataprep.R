rm(list=ls())
library("osfr")
library("data.table")

# download data file from OSF
osf_download(subset(osf_ls_files(osf_retrieve_node("8r7gz")),name=="Study 1"),
             conflicts="skip")

## prep data ##################################################################
rawData <- data.table::fread("Study 1/S1_data.csv")
rawData <- cbind("id"=1:nrow(rawData),rawData)
setkey(rawData)

str(rawData)

# treatment
## compare moderate/high vs low
rawData[, treat := (condition=="moderate" | condition=="high")*1L]

# outcome: life satisfaction
rawData[, life_satisf := mean(c(ls1,ls2,ls3,ls4,ls5), na.rm=TRUE),  by=id]
# rescale to between 0 and 1
rawData[, life_satisf := (life_satisf-1)/7]

# covariates: big five
rawData[, c("bfi_extra", "bfi_agree", "bfi_consci", "bfi_neuro", "bfi_open") := 
          list(
            mean(c(bfi1, bfi6, bfi11, bfi16, bfi21, bfi26, bfi31, bfi36), 
                 na.rm=TRUE),
            mean(c(bfi2, bfi7, bfi12, bfi17, bfi22, bfi27, bfi32, bfi37, bfi42), 
                 na.rm=TRUE),  
            mean(c(bfi3, bfi8, bfi13, bfi18, bfi23, bfi28, bfi33, bfi38, bfi43), 
                 na.rm=TRUE),  
            mean(c(bfi4, bfi9, bfi14, bfi19, bfi24, bfi29, bfi34, bfi39),
                 na.rm=TRUE),
            mean(c(bfi5, bfi10, bfi15, bfi20, bfi25, bfi30, bfi35, bfi40, bfi41, bfi44), 
                 na.rm=TRUE)),  by=id]

# mean-center variables
rawData[, c("bfi_extra", "bfi_agree", "bfi_consci", "bfi_neuro", "bfi_open") := 
          list(bfi_extra - mean(bfi_extra),
               bfi_agree - mean(bfi_agree),
               bfi_consci - mean(bfi_consci),
               bfi_neuro - mean(bfi_neuro),
               bfi_open - mean(bfi_open))] 

summary(rawData)

# other covariates
rawData[, gender := (gender==2)*1L]
L.Data <- rawData[, list(age,gender,
                         bfi_extra,bfi_agree,bfi_consci,bfi_neuro,bfi_open)]
setnames(L.Data, c("Age","Female",
                   "Extra","Agree","Consci","NegEmo","Open"))
setkey(L.Data)
L.Data

# covariates
L.names <- colnames(L.Data)
(p <- length(L.names))
mydata <- data.frame("i"=1:nrow(rawData), 
                     L.Data,
                     "treat"=rawData$treat,
                     "Y"=rawData$life_satisf)
round(head(mydata,4),2)
#   i Age Female Extra Agree Consci NegEmo  Open treat    Y
# 1 1  17      0 -0.89 -0.33  -1.41  -0.01  0.65     1 0.51
# 2 2  18      0 -1.76 -0.44  -1.30   0.86 -0.55     1 0.71
# 3 3  18      0 -1.51 -0.89  -0.52  -0.89  0.15     0 0.11
# 4 4  18      0 -1.01 -0.89  -0.63   1.49  0.55     1 0.00
