rm(list=ls())
libraries_check <- c("data.table","lslx","regsem","glmnet",
                     "twang","survey","CBPS")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs,repos="http://lib.ugent.be/CRAN/")
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# load helper functions
source("selection-stability.R")
source("selection-other_methods.R")

# simulation settings
simsets <- expand.grid(n=c(400),p=c(20,40),r=c(1,2,3)*0.5,x=c(0,1),z=c(1,1.5))

n_sims <- 2
# initialize for parallel MC jobs
args <- nrow(simsets)
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
  simsets <- simsets[rep(1:nrow(simsets),each=100),]
  nrow(simsets)
  n_sims <- 10
}
(seed <- as.integer(args[1]))
# check data-generating parameters
simsets[seed,]

(n <- simsets[seed,"n"]) # sample size
(p <- simsets[seed,"p"]) # number of covariates
(r <- simsets[seed,"r"]) # association b/w treatment and instrument
(x <- simsets[seed,"x"]) # normal (0) or non-normal (1) covariates
(z <- simsets[seed,"z"]) # association b/w outcome and outcome-only predictors

# indices for (C)onfounders, outcome (P)redictors, and (I)nstruments
lC <- 1:2
lP <- 3:4
lI <- 5:6

# load data generation
source("selection-sims-data_generating.R")

One_sim <- function() {
  mydata <- One_ObsData()
  var.list <- paste0("L.",1:p)
  
  L.selected <- list()
  L.selected[["none"]] <- "1" # ignoring all covariates
  L.selected[["unwanted.only"]] <- paste0("L.",lI) # instruments only
  L.selected[["target"]] <- paste0("L.",c(lC,lP)) # target set of covariates
  
  est.lslx <- NULL
  # variable selection procedures for the outcome model =======================
  # LASSO
  L.selected[["glmnet"]] <- OneData_glmnet(mydata,var.list,double.select=FALSE)
  
  if ("lslx" %in% (.packages())) {
    # Semi-Confirmatory SEM
    res.lslx <- OneData_lslx(mydata,var.list,double.select=FALSE)
    L.selected[["lslx"]] <- res.lslx$selected
    rm(res.lslx)
  }
  # regularized SEM
  if ("regsem" %in% (.packages())) {
    L.selected[["regsem"]] <- OneData_regsem(mydata,var.list,double.select=FALSE)
  }
  
  # order covariates based on priority to be confounders ====================== 
  L.ordered <- ForwardSelect_DS(Y=mydata$Y,X=mydata[,var.list],A=mydata$treat,
                                A.binary.par_model="mle")
  L.stable <- StdDiffEst_Ordered_DRAIPW(L.ordered$ordered,mydata,k=3,
                                        A.binary.par_model="mle")
  L.selected[["stability_mle"]] <- L.ordered$ordered[1:L.stable$selected_orbit]
  rm(L.stable,L.ordered)
  
  # check for empty adjustment sets and sort ==================================
  L.selected <- lapply(L.selected, function(l.select) {
    if (any(!grepl("L",l.select) | is.na(l.select)) | length(l.select)==0) {
      l.select <- "1" # ignoring all covariates
    }
    return(sort(l.select))
  })
  
  # calculate treatment effect estimate for each selected covariate set =======
  res <- lapply(L.selected, function(l.select) {
    trt <- OneDR_AIPW_Est(Ls=l.select,mydata,return.se=TRUE,bal.tab=TRUE,
                          A.binary.par_model="mle")
    return(unlist(c(
      trt,
      "numL.sel"=sum(grepl("L",l.select)), # number of covariates selected 
      "lcselect"=sum(paste0("L.",lC) %in% l.select) # number of confounders 
    )))
  })
  
  if (all(c("twang","CBPS") %in% (.packages()))) {
    res.ipw <- lapply(c("GBM","CBPS"), function(m)
      OneData_IPW(mydata,var.list=var.list,method=m))
    names(res.ipw) <- c("ipw.GBM","ipw.CBPS")
    res <- c(res,res.ipw)
    rm(res.ipw)
  }
  
  return( unlist(c(simsets[seed,],res)) )
}

ptm=proc.time()[3]
sim_res <- replicate(n=n_sims,expr=One_sim(),simplify=FALSE)
proc.time()[3]-ptm
# 15 mins per sim

warnings()

subfolder <- "selection-sim3/"
myfile <- "selection-sim3-"
save(sim_res,file=paste0(subfolder,myfile,seed,".Rdata"))
cat(seed,"\n")
q()
