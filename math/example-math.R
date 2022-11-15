rm(list=ls())
library("osfr")

# download data file from OSF
osf_download(osf_ls_files(osf_retrieve_node("8t3s9")),conflicts="skip")

## prep data ##################################################################
rawdata <- read.csv("Ex1Data.csv")
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
libraries_check <- c("data.table","lslx","regsem","glmnet","twang","survey")
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

# order covariates based on priority to be confounders ====================== 
L.ordered <- ForwardSelect_DS(Y=mydata$Y,X=mydata[,var.list],A=mydata$treat,
                              A.binary.par_model="mle")
L.stable <- StdDiffEst_Ordered_DRAIPW(L.ordered$ordered,mydata,
                                      A.binary.par_model="mle")
L.selected[["stability"]] <- L.ordered$ordered[1:L.stable$selected_orbit]

# ordered covariates in decreasing priority for confounding adjustment
L.names.ordered <- L.names[as.integer(lapply(sapply(
  L.ordered$ordered, strsplit, split="L."),"[",2))]

# make plots ==================================================================
filename <- "example-math-plots-stability-"
# Standardized difference between effect estimators in different orbits
pdf(paste0(filename,"1-std_diff.pdf"),width=6,height=4)
plot.df <- data.frame(L.stable$est)
plot.df[,"std.diff"] <- ifelse(plot.df[,"se"]==0,0.0,plot.df[,1]/plot.df[,2])
plot(plot.df[,"std.diff"], type="b",
     xlab="Adjustment set size",ylab="Std. Diff.", 
     main="New Math Curriculum")
sel.idx <- L.stable$selected_orbit
points(sel.idx,plot.df[sel.idx,"std.diff"],pch=19)
abline(h=0,lty=3)
dev.off()
# Q statistic
pdf(paste0(filename,"2-Qstat.pdf"),width=6,height=4)
plot.df <- c(NA,L.stable[[3]],NA)
plot(log(plot.df), type="b",
     xlab="Adjustment set size",ylab="Q statistic (log)", 
     main="New Math Curriculum")
points(sel.idx,log(plot.df[sel.idx]),pch=19)
dev.off()
rm(L.stable)

# Effect estimates 
pdf("example-math-plots-stability-3-eff_est.pdf",width=6,height=4)
plot.df <- data.frame(t(sapply(0:p, function(x) {
  if (x==0) {
    l.select <- "1" # now includes empty adjustment set
  } else {
    l.select <- L.ordered$ordered[1:x]
  }
  OneDR_AIPW_Est(Ls=l.select,mydata,return.se=TRUE)
})))
plot.df[,"l"] <- plot.df$EST-qnorm(.975)*plot.df$SE
plot.df[,"u"] <- plot.df$EST+qnorm(.975)*plot.df$SE
plot(0:p,plot.df[,"EST"], 
     ylim=range(c(0,plot.df[,c("l","u")])),xaxt="n",
     xlab="",ylab="Effect Estimate", main="New Math Curriculum")
axis(1, at=0:p, labels=c(NA,L.names.ordered),
     # tck=.01, cex.axis=0.5, 
     srt=45, col.ticks = "grey", las=2)
points(sel.idx,plot.df[sel.idx+1,"EST"],pch=19)
abline(h=0,lty=3)
for (x in 0:p) {
  lines(rep(x,2),plot.df[x+1,c("l","u")])
}
dev.off()

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

save.image("example-math-mle_only.Rdata")

# check for empty adjustment sets and sort ==================================
L.selected <- lapply(L.selected, function(l.select) {
  if (any(!grepl("L",l.select) | is.na(l.select)) | length(l.select)==0) {
    l.select <- "1" # ignoring all covariates
  }
  return(sort(l.select))
})

# calculate treatment effect estimate for each selected covariate set =======
res <- lapply(L.selected, function(l.select) {
  OneDR_AIPW_Est(Ls=l.select,mydata,return.se=TRUE,bal.tab=TRUE,
                 A.binary.par_model="mle")
})

if (all(c("twang","CBPS") %in% (.packages()))) {
  res.ipw <- lapply(c("GBM","CBPS"), function(m)
    OneData_IPW(mydata,var.list=var.list,method=m))
  names(res.ipw) <- c("ipw.GBM","ipw.CBPS")
  res <- c(res,res.ipw)
  rm(res.ipw)
}

res <- unlist(res)

save.image("example-math-mle_only.Rdata")
q()

load("example-math-mle_only.Rdata")

# selected using each method
L.selected.indicators <- data.frame(do.call(cbind,lapply(
  L.selected, function(l.select) 
    (1:length(L.names)) %in% 
    sort(as.integer(sapply(strsplit(l.select, split="L."),"[",2))))))
L.selected.indicators <- cbind(L.names,L.selected.indicators)

# sort original names in terms of ordered using stability
all.equal(L.names[order(match(L.names,L.names.ordered))],
          L.names.ordered)

library("xtable")
print.meths <- c("stability","glmnet","regsem","lslx","none","all")
print(xtable(L.selected.indicators[order(match(L.names,L.names.ordered)),
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
res.est[match(print.meths,res.est$meth),]

xtable(t(res.est[match(print.meths,res.est$meth),-(1:2)]))
