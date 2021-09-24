# covariate prioritization using double selection
ForwardSelect_DS <- function(Y,X,A,A.binary.par_model="mle") {
  Y.binary <- length(sort(unique(Y)))==2
  A.binary <- length(sort(unique(A)))==2
  X.names <- colnames(X)
  p <- length(X.names)
  n <- length(Y)
  X.curr <- NULL
  X.crit <- NULL
  while(length(X.curr) < p) {
    # find covariate among those not already in the model with minimin p-value
    X.cands <- X.names[!(X.names %in% X.curr)]
    Xj.crit <- sapply(X.cands, function(X.cand) {
      # create dataset with current and candidate covariates
      data.X.cand <- cbind("i"=1:n,"Y"=Y,"treat"=A,
                           X[,c(X.curr,X.cand),drop=FALSE])
      
      # treatment and outcomes models with the included covariates
      ps.X.cand <- as.formula(paste0("treat~",paste(
        c(X.curr,X.cand),collapse="+")))
      if (A.binary) {
        if (A.binary.par_model=="mle") {
          ps.X.cand.lm <- glm(ps.X.cand,family=binomial(link = "logit"),
                               data=data.X.cand)
          ps.X.cand.pv <- coef(summary(ps.X.cand.lm))[X.cand,"Pr(>|z|)"]
        } else if (A.binary.par_model=="CBPS") {
          ps.X.cand.lm <- CBPS::CBPS(formula=ps.X.cand,data=data.X.cand,ATT=0,
                                     iterations=10000)
          ## avoid using summary that is automatically printed
          ps.X.cand.est_std <- ( ps.X.cand.lm$coefficients/
                                   sqrt(diag(ps.X.cand.lm$var)) )[,1]
          ps.X.cand.pv <- pnorm(abs(ps.X.cand.est_std),lower.tail=FALSE)*2
          if (X.cand %in% names(ps.X.cand.pv)) {
            ps.X.cand.pv <- as.numeric(ps.X.cand.pv[X.cand])
          } else {
            ps.X.cand.pv <- as.numeric(ps.X.cand.pv[length(ps.X.cand.pv)])
          }
        }
      } else {
        ps.X.cand.lm <- lm(ps.X.cand,data=data.X.cand)
        ps.X.cand.pv <- coef(summary(ps.X.cand.lm))[X.cand,"Pr(>|t|)"]
      }
      
      out.X.cand <- as.formula(paste0("Y~",paste(
        c("treat",X.curr,X.cand),collapse="+")))
      
      if (Y.binary) {
        out.X.cand.lm <- glm(out.X.cand,family=binomial(link = "logit"),
                             data=data.X.cand)
        out.X.cand.pv <- coef(summary(out.X.cand.lm))[X.cand,"Pr(>|z|)"]
      } else {
        out.X.cand.lm <- lm(out.X.cand,data=data.X.cand)
        out.X.cand.pv <- coef(summary(out.X.cand.lm))[X.cand,"Pr(>|t|)"]
      }
      
      ## check for errors in treatment model coefficient estimates
      ps.OK <- all(!is.na(c(coef(ps.X.cand.lm),vcov(ps.X.cand.lm))))
      if (!ps.OK) {
        ps.X.cand.pv <- Inf
      }
      ## check for errors in outcome model coefficient estimates
      out.OK <- all(!is.na(c(coef(out.X.cand.lm),vcov(out.X.cand.lm))))
      if (!out.OK) {
        out.X.cand.pv <- Inf
      }
      return( min(ps.X.cand.pv,out.X.cand.pv,na.rm=TRUE) )
    })
    if (min(Xj.crit)<=1) {
      X.add <- X.cands[which.min(Xj.crit)] # valid p-values
    } else {
      X.add <- sample(x=X.cands,size=1) # randomly select a covariate
    }
    X.curr <- c(X.curr,X.add)
    X.crit <- c(X.crit,min(Xj.crit))
  }
  return(list("ordered"=X.curr,"pvalues"=X.crit))
}

# calculate treatment effect estimator and SE given a set of covariates
OneTrtCoef_Est <- function(Ls,mydata,return.lm=FALSE) {
  # outcome model
  out.Ls <- as.formula(paste0("Y~",paste(c("treat",Ls),collapse="+")))
  # linear regression model
  out.lm <- lm(out.Ls, data=mydata)
  if (return.lm==FALSE) {
    # fit the treatment model given covariates from the outcome model
    treat.Ls <- as.formula(paste0("treat~",paste(Ls,collapse="+")))
    treat.lm <- lm(treat.Ls, data=mydata)
    fitA.resids <- residuals(treat.lm)
    fitY.resids <- residuals(out.lm)
    # estimated influence functions using OLS estimators
    infn.est <- fitA.resids*fitY.resids/mean(mydata$treat*fitA.resids)
    est.var <- var(infn.est)/length(infn.est)
    if (abs(est.var)<.Machine$double.eps) {
      est.var <- 0.0 # avoid warning message if sqrt of small negative value
    }
    res <- c("EST"=coef(summary(out.lm))["treat","Estimate"],
             "SE"=sqrt(est.var))
    res["PV"] <- pnorm(abs(res["EST"])/res["SE"],lower.tail=FALSE)*2 # 2-sided
  } else {
    res <- out.lm
  }
  return( res ) 
}

TreatmentHeterogeneity <- function(x,s,k=NA) {
  if (is.na(k)) k <- length(x)
  # indices for subsets of size k
  windows <- sapply(0:(length(x)-k), function(i) (1:k)+i, simplify="matrix")
  # calculate Cochran's Q for each subset
  window.summ <- apply(windows, 2, function(i) {
    xi <- x[i]
    wi <- s[i]^(-2) # weights
    wi[is.infinite(wi)] <- 0 # zero weights if zero std err
    xbarw <- weighted.mean(x=xi,w=wi)
    CochranQ <- sum(wi*((xi-xbarw)^2))
    return( CochranQ )
  })
  # index for center of each window
  names(window.summ) <- windows[median(1:k),]
  return(window.summ)
}

StdDiffEst_Ordered_lm <- function(L.ordered,mydata,k,robust=FALSE) {
  # L.ordered: sequence of ordered covariate indices
  p_ <- length(L.ordered)
  diags.list <- lapply(1:p_, function(x) {
    OneTrtCoef_Est(Ls=L.ordered[1:x],mydata=mydata,return.lm=TRUE)
  })
  # differences with benchmark and their SE
  diags <- do.call(rbind,lapply(1:p_, function(x) {
    delta.est <- coef(diags.list[[x]])["treat"]-coef(diags.list[[p_]])["treat"]
    if (!robust) {
      # using CPH (1995) Eq.(18)
      ## requires homoscedasticity assumption
      delta.var <- vcov(diags.list[[p_]])["treat","treat"]-
        (vcov(diags.list[[x]])["treat","treat"]*
           (summary(diags.list[[p_]])$sigma^2)/(summary(diags.list[[x]])$sigma^2))
    } else {
      # using influence function from SIM paper 
      n_ <- nrow(mydata)
      infn.treat <- sapply(c(x,p_), function(x_) {
        # fit the treatment model given covariates from the outcome model
        fitA <- lm(
          as.formula(paste0("treat~",paste(L.ordered[1:x_],collapse="+"))),
          data=mydata)
        fitA.resids <- residuals(fitA)
        fitY.resids <- residuals(diags.list[[x_]])
        # estimated influence functions using OLS estimators
        return( fitA.resids*fitY.resids/mean(mydata$treat*fitA.resids) )
      })
      delta.var <- var(infn.treat[,1]-infn.treat[,2])/n_
    }
    if (abs(delta.var)<.Machine$double.eps) {
      delta.var <- 0.0 # avoid warning message if sqrt of small negative value
    }
    c("diff"=delta.est,"se"=sqrt(delta.var))
  }))
  
  # Q statistics for given window width
  Qk <- TreatmentHeterogeneity(x=diags[,"diff.treat"],
                               s=diags[,"se"],k=k)
  crits <- as.integer(names(which.min(Qk)))
  return(list("est"=diags,"selected_orbit"=crits,Qk))
}


# calculate DR-AIPW effect estimator
OneDR_AIPW_Est <- function(Ls, mydata, return.se=FALSE, bal.tab=FALSE,
                           A.binary.par_model="mle") {
  # initialize results
  if (return.se) {
    tau <- rep(NA,6)
    tau.names <- c("EST","SE","PV","ipw.sd",
                   "std.eff.sz.median","std.eff.sz.max")
    names(tau) <- tau.names
  } else {
    tau <- list("ate"=NA,"infn"=rep(NA,nrow(mydata)))
  }
  
  # propensity score model ====================================================
  ps.Ls <- as.formula(paste0("treat~",paste(Ls,collapse="+")))
  if (A.binary.par_model=="mle") {
    stm <- tryCatch(
      system.time(
        ps.fit <- glm(ps.Ls,family=binomial("logit"), data=mydata)),
      error=function(cond) return(NA))
  } else if (A.binary.par_model=="CBPS") {
    stm <- tryCatch(
      system.time(
        ps.fit <- CBPS::CBPS(formula=ps.Ls, data=mydata, ATT=0, 
                             iterations=10000)),
      error=function(cond) return(NA))
  }
  if (any(is.na(stm))) {
    return(tau)
  } else {
    ps.hat <- predict(ps.fit,type="response")
    # inverse probability weights
    ipw <- (1-mydata$treat)/(1-ps.hat) + mydata$treat/ps.hat
  }
  # error if there is an extreme propensity score estimate (up to machine error)
  if(any(pmin(ps.hat,1-ps.hat) < (.Machine$double.eps*1e1))) {
    return(tau)
  }
  
  # parametric outcome model ==================================================
  out.Ls <- as.formula(paste0("Y~",paste(c("treat",Ls),collapse="+")))
  Y.binary <- length(sort(unique(mydata$Y)))==2
  if (Y.binary) {
    # logistic regression model
    stm.out <- tryCatch(
      system.time(
        out.fit <- glm(out.Ls, family=binomial("logit"), data=mydata)),
      error=function(cond) return(NA))
  } else {
    # linear regression model
    stm.out <- tryCatch(
      system.time(
        out.fit <- glm(out.Ls, family=gaussian("identity"), data=mydata)),
      error=function(cond) return(NA))
  }
  if (any(is.na(stm.out))) {
    return(tau)
  } else {
    if(!out.fit$converged | out.fit$boundary | 
       any(is.na(coef(out.fit,na.rm=FALSE)))) {
      return(tau)
    }
  }
  
  # for predicting counterfactuals ============================================
  if (all(Ls %in% colnames(mydata))) {
    mydata.A1 <- mydata.A0 <- mydata[,c("treat",Ls)]
  } else {
    mydata.A1 <- mydata.A0 <- mydata[,"treat",drop=FALSE]
  }
  mydata.A1[,"treat"] <- 1
  mydata.A0[,"treat"] <- 0
  Y.A1.hat <- predict(out.fit,newdata=mydata.A1,type="response")
  Y.A0.hat <- predict(out.fit,newdata=mydata.A0,type="response")
  
  infn.est <- (mydata$treat/ps.hat)*(mydata$Y-Y.A1.hat) - 
    ((1-mydata$treat)/(1-ps.hat))*(mydata$Y-Y.A0.hat) +
    (Y.A1.hat-Y.A0.hat)
  
  # one-step plug-in estimator
  ate.hat <- mean(infn.est)
  # individual influence functions
  infn.est <- infn.est - ate.hat

  if (return.se) {
    tau["EST"] <- ate.hat
    tau["SE"] <- sqrt(var(infn.est)/length(infn.est))
    tau["PV"] <- pnorm(abs(tau["EST"])/tau["SE"],lower.tail=FALSE)*2 # 2-sided
    # SD of weights
    tau["ipw.sd"] <- sd(ipw)
    if (bal.tab==TRUE) {
      # median abs. std. bias across all covariates (Belitser et al., 2011)
      if (all(Ls %in% colnames(mydata))) {
        std.eff.sz <- cobalt::bal.tab(
          mydata[,Ls,drop=FALSE], treat=mydata$treat, weights=ipw,
          continuous="std",binary="std",s.d.denom="pooled")$Balance[,"Diff.Adj"]
        tau[5:6] <- quantile(abs(std.eff.sz),probs=c(0.5,1),na.rm=TRUE)
      }
    }
  } else {
    tau <- list("ate"=ate.hat,"infn"=infn.est)
  }
  return(tau)
}

StdDiffEst_Ordered_DRAIPW <- function(L.ordered,mydata,k,...) {
  # L.ordered: sequence of ordered covariate indices
  p_ <- length(L.ordered)
  diags.list <- lapply(1:p_, function(x) {
    OneDR_AIPW_Est(Ls=L.ordered[1:x],mydata=mydata,...)
  })
  diags <- do.call(rbind,lapply(diags.list, "[[", "ate"))
  score.list <- do.call(rbind,lapply(diags.list, "[[", "infn"))
  # max(abs(rowMeans(score.list)),na.rm=TRUE)<.Machine$double.eps*1e1
  
  # differences with benchmark and their SE
  diags <- do.call(rbind,lapply(1:p_, function(x) {
    delta.est <- diags[x,1]-diags[p_,1]
    delta.var <- var(score.list[x,]-score.list[p_,])/ncol(score.list)
    c("diff"=delta.est,"se"=sqrt(delta.var))
  }))
  
  # Q statistics for given window width
  Qk <- TreatmentHeterogeneity(x=diags[,"diff"],
                               s=diags[,"se"],k=k)
  crits <- as.integer(names(which.min(Qk)))
  return(list("est"=diags,"selected_orbit"=crits,Qk))
}
