rm(list=ls())
# setup follows sim 1 =========================================================
# read in preamble
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what=character(), 
                     skip=start-1, nlines=end-start+1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}
source2("selection-sim1.R",2,31)

# check seed: change here if needed
seed <- seed+nrow(simsets)/2 # limitation of max seed number
# check data-generating parameters
simsets[seed,]

source2("selection-sim1.R",32,45)
rm(source2)

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
  # variable selection procedures for the treatment and outcome models ========
  for (ds in c(FALSE,TRUE)) {
    # LASSO
    L.selected[[paste0("glmnet",ifelse(ds,".DS",""))]] <- 
      OneData_glmnet(mydata,var.list,double.select=ds)
  }
  # variable selection procedures for the outcome model only ==================
  ds <- FALSE
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
  
  # double selection using rigorous LASSO =====================================
  res.hdm <- OneData_hdm(mydata,var.list)
  
  # order covariates based on priority to be confounders ====================== 
  ## estimation method for propensity score model
  for (abin in c("mle","CBPS")) {
    ## MLE or CBPS
    L.ordered <- ForwardSelect_DS(Y=mydata$Y,X=mydata[,var.list],A=mydata$treat,
                                  A.binary.par_model=abin)
    L.stable <- StdDiffEst_Ordered_DRAIPW(L.ordered$ordered,mydata,k=3,
                                          A.binary.par_model=abin)
    L.selected[[paste0("stability_",abin)]] <- 
      L.ordered$ordered[1:L.stable$selected_orbit]
    rm(L.stable,L.ordered)
  }
  
  # check for empty adjustment sets and sort ==================================
  L.selected <- lapply(L.selected, function(l.select) {
    if (any(!grepl("L",l.select) | is.na(l.select)) | length(l.select)==0) {
      l.select <- "1" # ignoring all covariates
    }
    return(sort(l.select))
  })
  
  # calculate treatment effect estimate for each selected covariate set =======
  res <- NULL
  for (abin in c("mle","CBPS")) {
    if (abin=="mle") {
      ## MLE for treatment model
      meths.abin <- grep("CBPS",names(L.selected),value=TRUE,invert=TRUE)  
    } else {
      ## CBPS for treatment model
      meths.abin <- grep("mle",names(L.selected),value=TRUE,invert=TRUE)
    }
    res.abin <- lapply(L.selected[meths.abin], function(l.select) {
      trt <- OneDR_AIPW_Est(Ls=l.select,mydata,return.se=TRUE,bal.tab=TRUE,
                            A.binary.par_model=abin)
      return(unlist(c(
        trt,
        "numL.sel"=sum(grepl("L",l.select)), # number of covariates selected 
        "lcselect"=sum(paste0("L.",lC) %in% l.select) # number of confounders 
      )))
    })
    names(res.abin) <- paste0(names(res.abin),".",abin)
    res <- c(res,res.abin)
    rm(meths.abin,res.abin)
  }
  
  if (all(c("twang","CBPS") %in% (.packages()))) {
    res.ipw.ls <- lapply(
      list("all"=var.list, # all covariates
           "target"=L.selected$target), # target covariates
      function(x) {
        res.ipw <- lapply(c("GBM","CBPS"), function(m)
          OneData_IPW(mydata,var.list=x,method=m))
        names(res.ipw) <- c("ipw.GBM","ipw.CBPS")
        return(res.ipw)
      })
    res <- c(res,unlist(res.ipw.ls))
    rm(res.ipw.ls)
  }
  
  res <- c(res,
           unlist(est.lslx), # post-selection inference using lslx
           res.hdm[2]) # double-selection point estimates
  
  return( unlist(c(simsets[seed,],res)) )
}

ptm=proc.time()[3]
sim_res <- replicate(n=n_sims,expr=One_sim(),simplify=FALSE)
proc.time()[3]-ptm
# 25 mins per sim

warnings()

subfolder <- "selection-sim2/"
myfile <- "selection-sim2-"
save(sim_res,file=paste0(subfolder,myfile,seed,".Rdata"))
cat(seed,"\n")
q()
