rm(list=ls())
libraries_check <- c("data.table","lslx","regsem","SignifReg","hdm","glmnet",
                     "twang","survey","SBdecomp","CBPS")
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
simsets <- expand.grid(n=c(400),p=c(20,40),r=c(1,2,3)*0.8,x=c(0,1),
                       z=c(2,4),q=c(0,1),a=c(0,1))
## reduce strength of associations with binary treatment
simsets[simsets$a==1,"r"] <- simsets[simsets$a==1,"r"]/8*5

n_sims <- 2
# initialize for parallel MC jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
  simsets <- simsets[rep(1:nrow(simsets),each=100),]
  nrow(simsets)
  n_sims <- 10
}
(seed <- as.integer(args[1]))

(n <- simsets[seed,"n"]) # sample size
(p <- simsets[seed,"p"]) # number of covariates
(r <- simsets[seed,"r"]) # association b/w treatment and instrument or collider
(x <- simsets[seed,"x"]) # continuous (0) or binary (1) covariates
(z <- simsets[seed,"z"]) # number of instruments or colliders
(q <- simsets[seed,"q"]) # whether there is collider bias
(a <- simsets[seed,"a"]) # continuous (0) or binary (1) treatment

# indices for (C)onfounders, outcome (P)redictors, and (I)nstruments
lC <- 1:2
lP <- 3:4
lI <- 5+(0:(z-1))

# load data generation
source("selection-sims-data_generating.R")

One_sim <- function() {
  mydata <- One_ObsData()
  var.list <- paste0("L.",1:p)
  
  L.selected <- list()
  L.selected[["none"]] <- "1" # ignoring all covariates
  L.selected[["unwanted.only"]] <- paste0("L.",lI) # instruments or colliders
  L.selected[["target"]] <- paste0("L.",c(lC,lP)) # target set of covariates
  
  est.lslx <- NULL
  # variable selection procedures for the outcome model only ==================
  # # or double selection using treatment model as well
  for (ds in c(FALSE,TRUE)) {
    # LASSO
    L.selected[[paste0("glmnet",ifelse(ds,".DS",""))]] <- 
      OneData_glmnet(mydata,var.list,double.select=ds)
    ## significance-based variable selection
    L.selected[[paste0("SignifReg",ifelse(ds,".DS",""))]] <- 
      OneData_SignifReg(mydata,var.list,double.select=ds)
    if ("lslx" %in% (.packages())) {
      # Semi-Confirmatory SEM
      res.lslx <- OneData_lslx(mydata,var.list,double.select=ds)
      L.selected[[paste0("lslx",ifelse(ds,".DS",""))]] <- res.lslx$selected
      # post-selection treatment coefficient
      est.lslx[[ds+1]] <- res.lslx[2]
      rm(res.lslx)
    }
    # regularized SEM
    if ("regsem" %in% (.packages())) {
      L.selected[[paste0("regsem",ifelse(ds,".DS",""))]] <-
        OneData_regsem(mydata,var.list,double.select=ds)
    }
  }
  
  # double selection using rigorous LASSO =====================================
  res.hdm <- OneData_hdm(mydata,var.list)
  
  # order covariates based on priority to be confounders ====================== 
  L.ordered <- ForwardSelect_DS(Y=mydata$Y,X=mydata[,var.list],A=mydata$treat)
  L.stable <- StdDiffEst_Ordered_lm(L.ordered$ordered,mydata,k=3,robust=TRUE)
  L.selected[["stability"]] <- L.ordered$ordered[1:L.stable$selected_orbit]
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
    trt <- OneTrtCoef_Est(Ls=l.select,mydata)
    return(unlist(list(
      trt,
      "numL.sel"=sum(grepl("L",l.select)), # number of covariates selected 
      "lcselect"=sum(paste0("L.",lC) %in% l.select) # number of confounders 
    )))
  })
  
  res <- c(res,
           unlist(est.lslx), # post-selection inference using lslx
           res.hdm[2]) # double-selection point estimates
  
  return( unlist(c(simsets[seed,],res)) )
}

## for profiling code
# library(profvis);library(ggplot2)
# profvis(One_sim())

ptm=proc.time()[3]
sim_res <- replicate(n=n_sims,expr=One_sim(),simplify=FALSE)
proc.time()[3]-ptm
# 8 mins per sim

warnings()

subfolder <- "selection-sim1/"
myfile <- "selection-sim1-"
save(sim_res,file=paste0(subfolder,myfile,seed,".Rdata"))
cat(seed,"\n")
q()
