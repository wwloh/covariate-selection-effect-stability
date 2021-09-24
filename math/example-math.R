rm(list=ls())

# process data ================================================================
rawdata <- read.csv("osfstorage-archive/Ex1Data.csv")
names(rawdata)
n <- nrow(rawdata) # sample size

# covariates
L.names <- c("Sex", "Grade", "GradeExp", 
             "AfAm", "AmInd", "Cauc", "MexAm", "AsPIFi", "PR", "Other", 
             "AikPos", "AikNeg", 
             "CRNum", "CRAlg", "CRSpace", "CRMeas", "CRCD", 
             "Pretest")
(p <- length(L.names))
mydata <- data.frame("i"=1:n, 
                     rawdata[,L.names],
                     "treat"=rawdata$Treat,
                     "Y"=rawdata$Outcome)
var.list <- paste0("L.",1:p)
colnames(mydata) <- c("i",var.list,"treat","Y")
rm(rawdata)
summary(mydata)

# analysis ====================================================================
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
source("../selection-stability.R")
source("../selection-other_methods.R")

L.selected <- list()
L.selected[["none"]] <- "1" # ignoring all covariates
L.selected[["all"]] <- var.list # including all covariates

est.df <- est.sel.idx <- L.names.ordered <- NULL
# order covariates based on priority to be confounders ====================== 
for (abin in c("mle","CBPS")) {
  ## MLE or CBPS
  L.ordered <- ForwardSelect_DS(Y=mydata$Y,X=mydata[,var.list],A=mydata$treat,
                                A.binary.par_model=abin)
  L.stable <- StdDiffEst_Ordered_DRAIPW(L.ordered$ordered,mydata,k=3,
                                        A.binary.par_model=abin)
  L.selected[[paste0("stability_",abin)]] <- 
    L.ordered$ordered[1:L.stable$selected_orbit]

  # ordered covariates in decreasing priority for confounding adjustment
  L.names.ordered[[abin]] <- L.names[as.integer(lapply(sapply(
    L.ordered$ordered, strsplit, split="L."),"[",2))]

  # make plots ==================================================================
  filename <- "example-math-plots-stability-"
  # Standardized difference between effect estimators in different orbits
  pdf(paste0(filename,abin,"-1-std_diff.pdf"),width=6,height=4)
  plot.df <- data.frame(L.stable$est)
  plot.df[,"std.diff"] <- ifelse(plot.df[,"se"]==0,0.0,plot.df[,1]/plot.df[,2])
  plot(plot.df[,"std.diff"], type="b",
       xlab="Adjustment set size",ylab="Std. Diff.")
  sel.idx <- L.stable$selected_orbit
  points(sel.idx,plot.df[sel.idx,"std.diff"],pch=19)
  abline(h=0,lty=3)
  dev.off()
  
  # Q statistic
  pdf(paste0(filename,abin,"-2-Qstat.pdf"),width=6,height=4)
  plot.df <- c(NA,L.stable[[3]],NA)
  plot(log(plot.df), type="b",
       xlab="Adjustment set size",ylab="Q statistic (log)")
  points(sel.idx,log(plot.df[sel.idx]),pch=19)
  dev.off()

  # Effect estimates 
  plot.df <- data.frame(do.call(rbind,lapply(0:p, function(x) {
    if (x==0) {
      l.select <- "1" # now includes empty adjustment set
    } else {
      l.select <- L.ordered$ordered[1:x]
    }
    OneDR_AIPW_Est(Ls=l.select,mydata,return.se=TRUE,bal.tab=TRUE,
                   A.binary.par_model=abin)
  })))
  
  plot.df[,"CI.l"] <- plot.df[,"EST"]-qnorm(.975)*plot.df[,"SE"]
  plot.df[,"CI.u"] <- plot.df[,"EST"]+qnorm(.975)*plot.df[,"SE"]

  est.df[[abin]] <- plot.df
  est.sel.idx[[abin]] <- sel.idx
  rm(L.stable,L.ordered)
}

my_ylim <- range(unlist(lapply(est.df, function(x) 
  range(x[c("CI.l","CI.u")],na.rm=TRUE))))

for (abin in c("mle","CBPS")) { 
  plot.df <- est.df[[abin]]
  sel.idx <- est.sel.idx[[abin]]
  
  pdf(paste0("example-math-plots-stability-",abin,"-3-eff_est.pdf"),
      width=6,height=4)
  plot(0:p,plot.df[,"EST"], 
       ylim=my_ylim, xaxt="n",
       xlab="",ylab="Effect Estimate")
  axis(1, at=0:p, labels=c(NA,L.names.ordered[[abin]]), 
       tck=.01, cex.axis=0.9, srt=45, col.ticks = "grey", las=2)
  points(sel.idx,plot.df[sel.idx+1,"EST"],pch=19)
  abline(h=0,lty=3)
  for (x in 0:p) {
    lines(rep(x,2),plot.df[x+1,c("CI.l","CI.u")])
  }
  dev.off()
}

save.image("example-math.Rdata")

# variable selection procedures ===============================================
est.lslx <- NULL
# variable selection procedures for the outcome model only ==================
ds <- FALSE
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

save.image("example-math.Rdata")

# double selection using rigorous LASSO =====================================
res.hdm <- OneData_hdm(mydata,var.list)

# check for empty adjustment sets and sort ==================================
L.selected <- lapply(L.selected, function(l.select) {
  if (any(!grepl("L",l.select) | is.na(l.select)) | length(l.select)==0) {
    l.select <- "1" # ignoring all covariates
  }
  return(sort(l.select))
})

# calculate treatment effect estimate for each selected covariate set =======
res <- lapply(L.selected, function(l.select) {
  return(unlist(c(
    "CBPS"=OneDR_AIPW_Est(Ls=l.select,mydata,return.se=TRUE,bal.tab=TRUE,
                          A.binary.par_model="CBPS"),
    "MLE"=OneDR_AIPW_Est(Ls=l.select,mydata,return.se=TRUE,bal.tab=TRUE,
                         A.binary.par_model="mle"),
    "numL.sel"=sum(grepl("L",l.select)) # number of covariates selected 
  )))
})

if (all(c("twang","CBPS") %in% (.packages()))) {
  res.ipw <- lapply(c("GBM","CBPS"), function(m)
    OneData_IPW(mydata,var.list,method=m))
  names(res.ipw) <- c("all.ipw.GBM","all.ipw.CBPS")
  res <- c(res,res.ipw)
}

res <- c(res,
         unlist(est.lslx), # post-selection inference using lslx
         res.hdm[1:2]) # double-selection point estimates

res <- unlist(res)

save.image("example-math.Rdata")
q()

load("example-math.Rdata")

# selected using each method
L.selected.indicators <- data.frame(do.call(cbind,lapply(
  L.selected, function(l.select) 
    (1:length(L.names)) %in% 
    sort(as.integer(sapply(strsplit(l.select, split="L."),"[",2))))))
L.selected.indicators <- cbind(L.names,L.selected.indicators)

# sort original names in terms of ordered using CBPS
all.equal(L.names[order(match(L.names,L.names.ordered$CBPS))],
          L.names.ordered$CBPS)

library("xtable")
print.meths <- c("stability_CBPS","glmnet.DS","regsem","lslx","SignifReg",
                 "none","all")
print(xtable(L.selected.indicators[order(match(L.names,L.names.ordered$CBPS)),
                                   c("L.names",print.meths)]), 
      include.rownames=FALSE)

(meths <- grep("EST",names(res),value=TRUE))
res.summary <- NULL
for (mm in meths) {
  res.mm.pt <- data.frame(
    "meth"=strsplit(mm,".EST")[[1]],
    "DS"=grepl("DS",mm,ignore.case=TRUE),
    "EST"=res[mm],
    "SE"=res[gsub("EST","SE",mm)],
    row.names=NULL)
  res.mm.pt[,"CI.u"]=res.mm.pt["EST"]+qnorm(.975)*res.mm.pt["SE"]
  res.mm.pt[,"CI.l"]=res.mm.pt["EST"]-qnorm(.975)*res.mm.pt["SE"]
  
  # covariate balance
  for (cb in paste0("std.eff.sz.",c("median","max"))) {
    if (gsub("EST",cb,mm) %in% names(res)) {
      res.mm.pt[,cb]=as.numeric(res[gsub("EST",cb,mm)])
    } else {
      res.mm.pt[,cb]=NA
    }
  }
  
  res.summary <- c(res.summary, list(res.mm.pt))
  rm(res.mm.pt)
}
res.est <- do.call(rbind,res.summary)

print.meths.est <- sapply(print.meths, function(x) 
  ifelse(x %in% c("none","lslx"), paste0(x,".MLE"),paste0(x,".CBPS")))
res.est[match(print.meths.est,res.est$meth),]

xtable(t(res.est[match(print.meths.est,res.est$meth),-(1:2)]))
