# regularized SEM
OneData_regsem <- function(mydata,var.list,double.select=TRUE,
                           use.method="lasso") {
  model.regsem <- '
      # treatment model
      treat~L
      # outcome model
      Y~L+treat
    '
  if (!double.select) {
    # outcome model only
    model.regsem <- gsub("treat~L\n","",model.regsem)
  }
  # replace with covariate names and coefficients
  coef.trt <- gsub("L.","a",var.list)
  coef.out <- gsub("L.","b",var.list)
  
  model.regsem <- gsub("treat~L\n",paste0(
    "treat~",paste(coef.trt,var.list,sep="*",collapse="+"),"\n"),
    model.regsem,fixed=TRUE)
  model.regsem <- gsub("Y~L+treat\n",paste0(
    "Y~treat+",paste(coef.out,var.list,sep="*",collapse="+"),"\n"),
    model.regsem,fixed=TRUE)
  
  coef_to_penalize <- coef.out
  if (grepl("treat~",model.regsem)) {
    coef_to_penalize <- c(coef.trt,coef.out)
  }
  
  # fit full model using lavaan
  model.full <- sem(model.regsem, data = mydata, estimator = "ML")
  
  # initialize result
  if (use.method=="ridge") {
    select.regsem <- NA  
  } else {
    select.regsem <- var.list
  }
  
  # fit the regularized SEM
  if (use.method %in% c("lasso","ridge")) {
    fit.OK <- tryCatch(system.time(
      fit.regsem <- regsem::cv_regsem(
        model=model.full,
        pars_pen=coef_to_penalize,
        n.lambda = 100, # default in glmnet
        # mult.start = TRUE, # suggested in warning message but results in error
        alpha=ifelse(test=use.method=="lasso",yes=0,no=1), 
        type=use.method, # '1 = ridge, 0 = lasso' (opposite of glmnet) 
        verbose=FALSE))[3], error=function(cond) return(NA))
  } else if (use.method=="stabsel") {
    # stability-based LASSO for variable selection: takes too long
    fit.OK <- tryCatch(system.time(
      fit.regsem <- regsem::stabsel(
        model=model.full, data=mydata,
        n.lambda=100, from=1e-6,to=1,
        pars_pen=coef_to_penalize))[3], error=function(cond) return(NA))
  }
  # check if successfully converged
  if (all(!is.na(fit.OK))) {
    if (use.method %in% c("lasso","stabsel")) {
      # all covariates that are selected  
      coef.est <- fit.regsem$final_par
      parest.table <- data.frame(do.call(rbind,strsplit(names(coef.est)," ")))
      colnames(parest.table) <- c("lhs","op","rhs")
      parest.table <- cbind(parest.table,"est"=coef.est)
      row.names(parest.table) <- NULL
      
      select.regsem <- var.list[var.list %in% parest.table[
        parest.table$op=="->" & parest.table$est!=0,"lhs"]]
      select.regsem <- unique(select.regsem)
      if (length(select.regsem)==0) {
        select.regsem <- "1" # intercept only
      }
    } else {
      # ridge regression coefficient estimates
      select.regsem <- as.numeric(fit.regsem$final_par["treat -> Y"])
    }
  }
  return(select.regsem)
}

# Semi-Confirmatory SEM
OneData_lslx <- function(mydata,var.list,double.select=TRUE) {
  model.lslx <- '
      # treatment model
      treat<~L
      # outcome model
      Y<~L+treat
    '
  if (!double.select) {
    # outcome model only
    model.lslx <- gsub("treat<~L\n","",model.lslx)
  }
  # replace with covariate names and coefficients
  model.lslx <- gsub("treat<~L\n",paste0(
    "treat<~",paste(var.list,collapse="+"),"\n"),
    model.lslx,fixed=TRUE)
  model.lslx <- gsub("Y<~L+treat\n",paste0(
    "Y<~free()*treat+",paste(var.list,collapse="+"),"\n"),
    model.lslx,fixed=TRUE)
  
  coef_to_penalize <- paste0("Y<-",var.list)
  if (grepl("treat<~",model.lslx)) {
    coef_to_penalize <- c(paste0("treat<-",var.list),coef_to_penalize)
  }
  
  lslx_sem <- lslx$new(model = model.lslx, data = mydata, verbose = FALSE)
  lslx_sem$penalize_coefficient(name = coef_to_penalize, verbose = FALSE)
  lslx_sem$fit_lasso(start_method = "heuristic", verbose = FALSE)
  fit.OK <- tryCatch(system.time(
    coef.est <- lslx_sem$extract_coefficient(selector="abic")
  )[3], error=function(cond) return(NA))
  
  # initialize results
  select.lslx <- est.lslx <- NULL
  # check if successfully converged
  if (all(!is.na(fit.OK))) {
    # all covariates that are selected
    parest.table <- data.frame(do.call(
      rbind,lapply(strsplit(names(coef.est),"/g"),function(cname) {
        if (grepl("<->",cname)) {
          lrhs <- strsplit(cname,"<->")[[1]]
          c(lrhs[1],"~~",lrhs[2])
        } else {
          lrhs <- strsplit(cname,"<-")[[1]]
          c(lrhs[1],"<-",lrhs[2])
        }
      })
    ))
    colnames(parest.table) <- c("lhs","op","rhs")
    parest.table <- cbind(parest.table,"est"=coef.est)
    row.names(parest.table) <- NULL
  
    select.lslx <- var.list[var.list %in% parest.table[
      parest.table$op=="<-" & parest.table$rhs%in%var.list & 
        parest.table$est!=0,"rhs"]]
    select.lslx <- unique(select.lslx)
    
    # treatment coefficient estimate and SE
    posi.OK <- tryCatch(system.time(
      est.lslx <- lslx_sem$test_coefficient(selector="abic",
                                            standard_error="sandwich",
                                            debias="one_step",
                                            inference="scheffe")
    )[3], error=function(cond) return(NA))
    # Error: No non-zero parameters are selected and hence 
    # post-selection inference cannot be implemented.
    if (all(!is.na(posi.OK))) {
      est.lslx <- est.lslx["Y<-treat/g",c("estimate","standard_error","p_value")]
      est.lslx <- as.numeric(est.lslx)
    }
  }
  if (length(select.lslx)==0 || all(is.null(est.lslx))) {
    select.lslx <- "1" # intercept only
    # simple linear regression
    est.lslx <- coef(summary(lm(Y~treat,data=mydata)))["treat",c(1:2,4)]
  }
  names(est.lslx) <- c("EST","SE","PV")
  
  res.lslx <- list(select.lslx,est.lslx)
  names(res.lslx) <- c("selected",
                       paste0("post_lslx",ifelse(double.select,".DS","")))
  return(res.lslx)
}

OneData_glmnet <- function(mydata,var.list,double.select=TRUE) {
  select.glmnet <- NULL
  A.binary <- length(unique(mydata$treat))==2
  if (double.select) {
    # treatment model
    fit.OK <- tryCatch(system.time(
      fitA.cv <- cv.glmnet(
        x=as.matrix(mydata[,var.list]),
        y=mydata$treat,
        family=ifelse(test=A.binary,yes="binomial",no="gaussian"),
        type.measure=ifelse(test=A.binary,yes="class",no="mse"),
        standardize=TRUE,trace.it=FALSE,
        alpha=1 # LASSO
        ))[3],error=function(cond) return(NA))
    # check if successfully converged
    if (all(!is.na(fit.OK))) {
      # all covariates that are selected
      select.glmnet[["A"]] <- coef(fitA.cv, s="lambda.min")[,1]
      rm(fitA.cv)
    }
    rm(fit.OK)
  }
  # outcome model with treatment and all covariates
  fit.OK <- tryCatch(system.time(
    fitY.cv <- cv.glmnet(x=as.matrix(mydata[,c("treat",var.list)]),
                         y=mydata$Y,
                         family="gaussian",type.measure="mse",
                         standardize=TRUE,trace.it=FALSE,
                         alpha=1, # LASSO
                         penalty.factor=c(0,rep(1,length(var.list)))))[3],
    error=function(cond) return(NA))
  # check if successfully converged
  if (all(!is.na(fit.OK))) {
    # all covariates that are selected
    select.glmnet[["Y"]] <- coef(fitY.cv, s="lambda.min")[,1]
    rm(fitY.cv)
  }
  rm(fit.OK)
  
  select.glmnet <- unique(unlist(lapply(select.glmnet, function(fit_) 
    names(fit_)[abs(fit_)>.Machine$double.eps])))
  select.glmnet <- var.list[var.list %in% select.glmnet]
  if (length(select.glmnet)==0) {
    select.glmnet <- "1" # intercept only
  }
  return(select.glmnet)
}

OneData_SignifReg <- function(mydata,var.list,double.select=TRUE) {
  data.SignifReg <<- mydata[,c("Y","treat",var.list)]
  
  # forward selection can only start with intercept-only model
  ## outcome model
  fitY.SignifReg <<- lm(Y~1,data=data.SignifReg)
  lmY.SignifReg <- SignifReg(fitY.SignifReg, scope=as.formula(
    paste("~.+",paste(c("treat",var.list),collapse="+"))))
  selectAY.SignifReg <- names(lmY.SignifReg$coefficients)
  
  # treatment model
  if (double.select) {
    fitA.SignifReg <<- lm(treat~1,data=data.SignifReg)
    lmA.SignifReg <- SignifReg(fitA.SignifReg, scope=as.formula(
      paste("~.+",paste(var.list,collapse="+"))))
    selectA.SignifReg <- names(lmA.SignifReg$coefficients)
    selectAY.SignifReg <- c(selectAY.SignifReg,selectA.SignifReg)
  }
  
  select.SignifReg <- var.list[var.list %in% selectAY.SignifReg]
  
  if (length(select.SignifReg)==0) {
    select.SignifReg <- "1" # intercept only
  }
  return(select.SignifReg)
}

OneData_hdm <- function(mydata,var.list) {
  # Following example in \S4.4 of 
  # https://cran.r-project.org/web/packages/hdm/vignettes/hdm.pdf
  y = mydata$Y
  d = mydata$treat
  X = as.matrix(mydata[,var.list])
  
  lasso.effect = rlassoEffect(x = X, y = y, d = d, method = "partialling out")
  doublesel.effect = rlassoEffect(x = X, y = y, d = d, method = "double selection")
  
  hdm_res <- list(
    c("EST"=lasso.effect$alpha,
      "SE"=lasso.effect$se,
      "PV"=lasso.effect$pval),
    c("EST"=as.numeric(doublesel.effect$alpha),
      "SE"=as.numeric(doublesel.effect$se),
      "PV"=as.numeric(doublesel.effect$pval)))
  names(hdm_res) <- c("rlassoEffect.po","rlassoEffect.ds")
  
  return(hdm_res)
}

OneData_IPW <- function(mydata,var.list,method) {
  # initialize results
  ate <- rep(NA,6)
  names(ate) <- c("EST","SE","PV","ipw.sd","std.eff.sz.median","std.eff.sz.max")
  # propensity score model with all covariates
  ps.form <- as.formula(paste0("treat~",paste(var.list,collapse="+")))
  if (method=="GBM") {
    fit_OK <- tryCatch(
      ps.fit <- twang::ps(ps.form, estimand="ATE", verbose=FALSE, data=mydata,
                          stop.method="es.max"),
      error=function(cond) return(NA))
  } else if (method=="CBPS") {
    fit_OK <- tryCatch(
      ps.fit <- CBPS::CBPS(formula=ps.form, data=mydata, ATT=0, 
                           iterations=10000),
      error=function(cond) return(NA))  
  }
  if (all(!is.na(fit_OK))) {
    # inverse probability of treatment weights
    if (method=="GBM") {
      ipw <- twang::get.weights(ps.fit, stop.method="es.max")
    } else if (method=="CBPS") {
      ipw <- ps.fit$weights
    }
    # IPW estimator without any outcome model
    mydata.ipw <- mydata
    mydata.ipw[,"w"] <- ipw
    design.ps <- survey::svydesign(ids=~1, weights=~w, data=mydata.ipw)
    rm(mydata.ipw)
    glm.fit_OK <- tryCatch(
      glm.fit <- svyglm(formula=Y ~ treat, design=design.ps),
      error=function(cond) return(NA))
    if (all(!is.na(glm.fit_OK))) {
      ate[c("EST","SE","PV")] <- as.numeric(
        summary(glm.fit)$coef["treat",c(1:2,4)])
      # SD of weights
      ate["ipw.sd"] <- sd(ipw)
      # median abs. std. bias across all covariates (Belitser et al., 2011)
      if (method=="GBM") {
        std.eff.sz <- twang::bal.table(ps.fit)$es.max.ATE$std.eff.sz
      } else if (method=="CBPS") {
        std.eff.sz <- apply(balance(ps.fit)$balanced,1,function(x) 
          diff(x[c("0.std.mean","1.std.mean")]))
      }
      ate[5:6] <- quantile(abs(std.eff.sz),probs=c(0.5,1),na.rm=TRUE)
    }
  }
  names(ate) <- c("EST","SE","PV","ipw.sd","std.eff.sz.median","std.eff.sz.max")
  return(ate)
}
