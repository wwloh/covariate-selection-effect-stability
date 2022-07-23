rm(list=ls())
library("data.table")

subfolder <- "selection-sim3/"
myfiles <- list.files(subfolder)
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
simres <- NULL
for (ll in myfiles) {
  load(paste0(subfolder,ll))
  simres <- c(simres,sim_res) # observed estimate
  rm(sim_res)
}
sim_res <- do.call(rbind,simres); rm(simres)

simres <- data.table(sim_res)
setkey(simres)
(simsets <- simres[,.N,by=list(n,p,r,x,z)]) # unique sims per setting
par_idx <- 1:(ncol(simsets)-1)

# selecting confounders vs. number selected ===================================
(meths <- grep("lcselect",names(simres),value=TRUE))
meths <- meths[-(1:3)]
simres.summary <- NULL
for (mm in meths) {
  
  # both confounders selected
  simres.mm.lcselect <- simres[,as.list(unlist(lapply(.SD, function(i) 
    mean(i==2)))),
    by=list(n,p,r,x,z),.SDcols=mm]
  
  # number of covariates selected
  simres.mm.numL.sel <- simres[,as.list(unlist(lapply(.SD, function(i)
    c(mean(i),quantile(i,probs=c(.25,.75)))/p))),
    by=list(n,p,r,x,z),.SDcols=gsub("lcselect","numL.sel",mm)]
  
  simres.mm <- merge(simres.mm.lcselect,simres.mm.numL.sel)
  setnames(simres.mm,max(par_idx)+(1:4),
           c("confounders",paste0("total",c("",".lo",".up"))))
  simres.mm[, "meth" := strsplit(mm,".lcselect")[[1]]]
  simres.mm[, "DS" := grepl("DS",mm,ignore.case=TRUE) | 
              grepl("stability",mm,ignore.case=TRUE)]
  setcolorder(simres.mm,c(par_idx,ncol(simres.mm)-(1:0)))
  simres.summary[[strsplit(mm,".lcselect")[[1]]]] <- simres.mm
  rm(simres.mm,simres.mm.lcselect,simres.mm.numL.sel)
}

## methods to be plotted
meths <- NULL
meths[[1]] <- c("stability_mle","glmnet","regsem","lslx")
meths.nicenames <- NULL
meths.nicenames[[1]] <- c("stability","LASSO","RegSEM","SC-SEM")

# make plots ==================================================================
for (i in c(1:24)) {
  dat.plot <- rbindlist(lapply(simres.summary, function(dt) 
    dt[simsets[i,..par_idx]]))
  my_xlim <- range(unlist(dat.plot[,list(total.lo,total.up)]))
  my_ylim <- c(0,1)
  my_pch <- c(19,0:2)
  
  filename <- paste(names(simsets)[par_idx],simsets[i,..par_idx],
                    sep="_",collapse="-")
  filename <- gsub(pattern="[.]",replacement="",x=filename) # remove periods
  pdf(paste0("figures/plot-sim3-selected-",filename,".pdf"),
      width=4,height=4,bg=NULL)
  for (mm in 1:length(meths)) {
    plot(my_xlim,my_ylim,type="n",
         xlab="Prop. of all covariates selected",
         ylab="Prob. both confounders selected",
         main=paste0("J=",simsets[i,p], " covariates"))
    for (mm.i in meths[[mm]]) {
      points(dat.plot[meth==mm.i,list(total,confounders)],
             pch=my_pch[which(mm.i==meths[[mm]])],cex=1.25)
      lines(dat.plot[meth==mm.i,list(total.lo,total.up)],
            dat.plot[meth==mm.i,rep(confounders,2)])
    }
    legend("bottomright",legend=meths.nicenames[[mm]], border=FALSE, 
           pch=my_pch, pt.cex=1.25,cex=1, bty="n")
  }
  dev.off()
}

# point estimates and p-values ================================================
simres.summary.selected <- simres.summary
simres.summary <- NULL
meths <- grep("EST",names(simres),value=TRUE)
for (mm in meths) {
  simres.mm.pt <- simres[,as.list(unlist(lapply(.SD, function(i) {
    c("est"=mean(i,na.rm=TRUE),
      "ese"=sd(i,na.rm=TRUE),
      "rmse"=sqrt(mean(i^2,na.rm=TRUE)),
      "pr.NAs"=mean(is.na(i)))
  }))),by=list(n,p,r,x,z),.SDcols=mm]
  
  simres.mm.se <- simres[,as.list(unlist(lapply(.SD, function(i)
    mean(i,na.rm=TRUE)))),
    by=list(n,p,r,x,z),.SDcols=gsub("EST","SE",mm)]
  
  simres.mm.pv <- simres[,as.list(unlist(lapply(.SD, function(i) 
    mean(i<=0.05,na.rm=TRUE)))),
    by=list(n,p,r,x,z),.SDcols=gsub("EST","PV",mm)]
  
  simres.mm <- merge(simres.mm.pt,simres.mm.se)
  simres.mm <- merge(simres.mm,simres.mm.pv)
  setnames(simres.mm,max(par_idx)+(1:6),c("bias","ese","rmse",
                                          "pr.NAs","ase","typeI"))
  simres.mm[, "meth" := strsplit(mm,".EST")[[1]]]
  simres.mm[, "DS" := grepl("DS",mm,ignore.case=TRUE) | 
              grepl("stability",mm,ignore.case=TRUE)]
  setcolorder(simres.mm,c(par_idx,ncol(simres.mm)-(1:0)))
  
  # covariate balance
  if (gsub("EST","std.eff.sz.max",mm) %in% colnames(simres)) {
    simres.mm.covbal <- simres[,as.list(unlist(lapply(.SD, function(i) 
      c(mean(i,na.rm=TRUE), quantile(i,probs=1,na.rm=TRUE))))),
      by=list(n,p,r,x,z),.SDcols=gsub("EST","std.eff.sz.max",mm)]
    simres.mm <- merge(simres.mm, simres.mm.covbal)
  } else {
    simres.mm[,paste0(gsub("EST","std.eff.sz.max",mm),1:2) := NA]
  }
  setnames(simres.mm,ncol(simres.mm)-(1:0),paste0("max_sb.",c("av","max")))

  # merge with covariate selection summaries
  mm.name <- strsplit(mm,".EST")[[1]]
  if (mm.name %in% names(simres.summary.selected)) {
    simres.mm <- merge(simres.mm, simres.summary.selected[[mm.name]],
                       by=c(names(simres.mm)[par_idx],"meth","DS"))
  } else {
    simres.mm[,c("confounders","total","total.lo","total.up") := NA]
  }
  simres.summary[[mm.name]] <- simres.mm
  rm(simres.mm,simres.mm.pt,simres.mm.se,simres.mm.pv)
}
(simres.summary <- rbindlist(simres.summary,fill=TRUE))
setkey(simres.summary)

# more informative column names
setnames(simres.summary,par_idx,c("sample_size_n",
                                  "number_covariates_J",
                                  "coefficient_k",
                                  "nonnormal_covariates_1",
                                  "coefficient_yonly"))

write.csv(simres.summary, 
          file=paste0(strsplit(subfolder,"/")[[1]],"-results.csv"))

# for tables in paper =========================================================
library("xtable")
printcols <- max(par_idx)+c(1,3:5,7:8)
# include standardized bias summaries
printcols <- c(printcols,max(par_idx)+c(9:10))

xtable(simres.summary[simsets[1,..par_idx],..printcols])
