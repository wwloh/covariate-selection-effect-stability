One_ObsData <- function() {
  # covariates
  l <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
  if (x==1) {
    l.lognormal <- matrix(exp(rnorm(n=n*(p/2))),nrow=n,ncol=p/2)
    l.lognormal <- scale(l.lognormal, center=TRUE, scale=TRUE)
    l[,2*(1:(p/2))-1] <- l.lognormal
    rm(l.lognormal)
  }
  
  # treatment
  g_raw <- rep(0,p)
  g_raw[lC] <- 1 # true confounders
  g_raw[lI] <- r # instruments
  linear.ps <- (l %*% g_raw)[,1]
  # binary treatment
  X <- rbinom(n,1,exp(linear.ps)/(1+exp(linear.ps)))
  rm(linear.ps,g_raw)

  # outcome unaffected by treatment
  b_raw <- rep(0,p)
  b_raw[lC] <- 1 # true confounders
  b_raw[lP] <- z # predictors of outcome only
  linear.out <- (l %*% b_raw)[,1]
  y0 <- linear.out + rnorm(n=n, sd=max(abs(linear.out)))
  
  mydata <- data.frame("i"=1:n,"L"=l,"treat"=X,"Y"=y0)
  return(mydata)
}
