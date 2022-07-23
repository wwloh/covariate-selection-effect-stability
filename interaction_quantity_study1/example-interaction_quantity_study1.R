rm(list=ls())
library("data.table")

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
rawData[, life_satisf := (life_satisf-min(life_satisf))/
          diff(range(life_satisf))]

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
var.list <- paste0("L.",1:p)
colnames(mydata) <- c("i",var.list,"treat","Y")
rm(rawData,L.Data)

summary(mydata)

# analysis ====================================================================
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
source("../selection-stability.R")
source("../selection-other_methods.R")

L.selected <- list()
L.selected[["none"]] <- "1" # ignoring all covariates
L.selected[["all"]] <- var.list # including all covariates

# order covariates based on priority to be confounders ====================== 
L.ordered <- ForwardSelect_DS(Y=mydata$Y,X=mydata[,var.list],A=mydata$treat,
                              A.binary.par_model="mle")
L.stable <- StdDiffEst_Ordered_DRAIPW(L.ordered$ordered,mydata,k=3,
                                      A.binary.par_model="mle")
L.selected[["stability"]] <- L.ordered$ordered[1:L.stable$selected_orbit]

# ordered covariates in decreasing priority for confounding adjustment
L.names.ordered <- L.names[as.integer(lapply(sapply(
  L.ordered$ordered, strsplit, split="L."),"[",2))]

# make plots ==================================================================
filename <- "example-interaction_quantity_study1-plots-stability-"
# Standardized difference between effect estimators in different orbits
pdf(paste0(filename,"1-std_diff.pdf"),width=6,height=4)
plot.df <- data.frame(L.stable$est)
plot.df[,"std.diff"] <- ifelse(plot.df[,"se"]==0,0.0,plot.df[,1]/plot.df[,2])
plot(plot.df[,"std.diff"], type="b",
     xlab="Adjustment set size",ylab="Std. Diff.", 
     main="Life Satisfaction")
sel.idx <- L.stable$selected_orbit
points(sel.idx,plot.df[sel.idx,"std.diff"],pch=19)
abline(h=0,lty=3)
dev.off()
# Q statistic
pdf(paste0(filename,"2-Qstat.pdf"),width=6,height=4)
plot.df <- c(NA,L.stable[[3]],NA)
plot(log(plot.df), type="b",
     xlab="Adjustment set size",ylab="Q statistic (log)", 
     main="Life Satisfaction")
points(sel.idx,log(plot.df[sel.idx]),pch=19)
dev.off()
rm(L.stable)

# Effect estimates 
pdf("example-interaction_quantity_study1-plots-stability-3-eff_est.pdf",width=6,height=4)
plot.df <- data.frame(t(sapply(0:p, function(x) {
  if (x==0) {
    l.select <- "1" # now includes empty adjustment set
  } else {
    l.select <- L.ordered$ordered[1:x]
  }
  OneTrtCoef_Est(Ls=l.select,mydata)
})))
plot.df[,"l"] <- plot.df$EST-qnorm(.975)*plot.df$SE
plot.df[,"u"] <- plot.df$EST+qnorm(.975)*plot.df$SE
plot(0:p,plot.df[,"EST"], 
     ylim=range(c(0,plot.df[,c("l","u")])),xaxt="n",
     xlab="",ylab="Effect Estimate", main="Life Satisfaction")
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

save.image("example-interaction_quantity_study1.Rdata")

# check for empty adjustment sets and sort ==================================
L.selected <- lapply(L.selected, function(l.select) {
  if (any(!grepl("L",l.select) | is.na(l.select)) | length(l.select)==0) {
    l.select <- "1" # ignoring all covariates
  }
  return(sort(l.select))
})

# calculate treatment effect estimate for each selected covariate set =======
res <- lapply(L.selected, function(l.select) {
  OneTrtCoef_Est(Ls=l.select,mydata)
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

save.image("example-interaction_quantity_study1.Rdata")
q()

load("example-interaction_quantity_study1.Rdata")

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
  
  res.summary <- c(res.summary, list(res.mm.pt))
  rm(res.mm.pt)
}
res.est <- do.call(rbind,res.summary)
res.est[match(print.meths,res.est$meth),]

xtable(t(res.est[match(print.meths,res.est$meth),-(1:2)]))
