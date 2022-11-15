rm(list=ls())
library("data.table")
library("osfr")

# download data file from OSF
osf_download(osf_ls_files(osf_retrieve_node("8t3s9")),conflicts="skip")

## prep data ##################################################################
rawData <- read.csv("Ex1Data.csv")
names(rawData)

# covariates
L.names <- c("Sex", "Grade", "GradeExp", 
             "AfAm", "AmInd", "Cauc", "MexAm", "AsPIFi", "PR", "Other", 
             "AikPos", "AikNeg", 
             "CRNum", "CRAlg", "CRSpace", "CRMeas", "CRCD", 
             "Pretest")
(p <- length(L.names))
mydata <- data.frame("i"=1:nrow(rawData), 
                     rawData[,L.names],
                     "treat"=rawData$Treat,
                     "Y"=rawData$Outcome)
head(mydata,4)
#   i Sex Grade GradeExp AfAm AmInd Cauc MexAm AsPIFi PR Other AikPos
# 1 1   0     9        4    0     0    1     0      0  0     0     36
# 2 2   1     9        4    0     0    1     0      0  0     0     31
# 3 3   0     9        4    0     0    1     0      0  0     0     25
# 4 4   1     9        3    0     0    1     0      0  0     0     30
#   AikNeg CRNum CRAlg CRSpace CRMeas CRCD Pretest treat  Y
# 1     15     1     2       2      2    1      23     0 32
# 2     14     4     1       3      2    1      23     0 32
# 3     27     2     1       2      2    1      24     0 38
# 4     33     1     1       2      2    1      24     0 33