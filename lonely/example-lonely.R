rm(list=ls())

# process data ================================================================
df <- read.csv("data_processed_for wen.csv")
names(df)
summary(df)

# possibly affected by treatment
df[,c("ls","bored","diagnose","compliance",
      "contact_offline","contact_online")] <- NULL

# recode missing ethnicity as separate level
df[which(is.na(df),arr.ind=TRUE)[,"row"],]
df[is.na(df$ethnicity),"ethnicity"] <- max(df$ethnicity,na.rm=TRUE)+1L

# continuous variables
L.cont <- c("age","burger","bfi_extra","bfi_agree","bfi_consci",
            "bfi_negemo","bfi_open","pre_lonely","pre_ls","pre_bored",
            "lonely")
for (cn in L.cont) {
  df[,cn] <- as.numeric(df[,cn])  
}

# categorical covariates to be recoded as factors
L.fact <- names(df)[!(names(df) %in% c(L.cont,"isolation_duration","ls"))]
for (cn in L.fact) {
  df[,cn] <- as.factor(df[,cn])  
}
rm(cn,L.cont,L.fact)

summary(df)

# create model matrix with dummy variables for categorical variables
rawdata <- data.frame(model.matrix(lm(lonely~.,data=df))[,-1])
rawdata <- cbind(rawdata,"lonely"=df$lonely)
n <- nrow(rawdata) # sample size

# covariates
(L.names <- colnames(rawdata)[!(colnames(rawdata) %in% 
                                  c("isolation_duration","lonely"))])
(p <- length(L.names))
mydata <- data.frame("i"=1:n, 
                     rawdata[,L.names],
                     "treat"=rawdata$isolation_duration,
                     "Y"=rawdata$lonely)
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

# order covariates based on priority to be confounders ====================== 
L.ordered <- ForwardSelect_DS(Y=mydata$Y,X=mydata[,var.list],A=mydata$treat)
L.stable <- StdDiffEst_Ordered_lm(L.ordered$ordered,mydata,k=3,robust=TRUE)
L.selected[["stability"]] <- L.ordered$ordered[1:L.stable$selected_orbit]

# ordered covariates in decreasing priority for confounding adjustment
L.names.ordered <- L.names[as.integer(lapply(sapply(
  L.ordered$ordered, strsplit, split="L."),"[",2))]

# make plots ==================================================================
filename <- "example-lonely-plots-stability-"
# Standardized difference between effect estimators in different orbits
pdf(paste0(filename,"1-std_diff.pdf"),width=6,height=4)
plot.df <- data.frame(L.stable$est)
plot.df[,"std.diff"] <- ifelse(plot.df[,"se"]==0,0.0,plot.df[,1]/plot.df[,2])
plot(plot.df[,"std.diff"], type="b",
     xlab="Adjustment set size",ylab="Std. Diff.", main="Loneliness")
sel.idx <- L.stable$selected_orbit
points(sel.idx,plot.df[sel.idx,"std.diff"],pch=19)
abline(h=0,lty=3)
dev.off()
# Q statistic
pdf(paste0(filename,"2-Qstat.pdf"),width=6,height=4)
plot.df <- c(NA,L.stable[[3]],NA)
plot(log(plot.df), type="b",
     xlab="Adjustment set size",ylab="Q statistic (log)")
points(sel.idx,log(plot.df[sel.idx]),pch=19)
dev.off()
rm(L.stable)

# Effect estimates 
pdf("example-lonely-plots-stability-3-eff_est.pdf",width=6,height=4)
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
     ylim=range(plot.df[,c("l","u")]),xaxt="n",
     xlab="",ylab="Effect Estimate")
axis(1, at=0:p, labels=c(NA,L.names.ordered),
     tck=.01, cex.axis=0.5, srt=45, col.ticks = "grey", las=2)
points(sel.idx,plot.df[sel.idx+1,"EST"],pch=19)
abline(h=0,lty=3)
for (x in 0:p) {
  lines(rep(x,2),plot.df[x+1,c("l","u")])
}
dev.off()

# variable selection procedures ===============================================
est.lslx <- NULL
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

save.image("example-lonely.Rdata")

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
  OneTrtCoef_Est(Ls=l.select,mydata)
})

res <- c(res,
         unlist(est.lslx), # post-selection inference using lslx
         res.hdm[1:2]) # double-selection point estimates

res <- unlist(res)

save.image("example-lonely.Rdata")
q()

load("example-lonely.Rdata")

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
print.meths <- c("stability","glmnet.DS","regsem","lslx","SignifReg.DS",
                 "none","all")
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
