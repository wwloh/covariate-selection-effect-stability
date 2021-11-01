One_ObsData <- function() {
  # covariates
  if (x==0) {
    l <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
  } else {
    l <- matrix(rbinom(n=n*p,size=1,prob=0.1),nrow=n,ncol=p)
    # skewness = (1-2*.01)/sqrt(.01*(1-.01))
  }
  u <- matrix(0,nrow=n,ncol=2)
  if (q==1) {
    # unmeasured confounders to possibly induce collider bias
    u <- matrix(rnorm(n=n*2,sd=1/4),nrow=n,ncol=2)
    # covariates that induce collider bias and should not be adjusted for
    for (k in lI) {
      l[,k] <- rnorm(n=n,mean=2*rowSums(u),sd=sqrt(1/2)) # variance=1
    }
  }
  
  # treatment
  g_raw <- rep(0,p)
  g_raw[lC] <- 1 # true confounders
  g_raw[lI] <- r*(1-q) # instruments when q==0
  nu <- r*q # to induce collider bias when q==1
  linear.ps <- (l %*% g_raw)[,1] + nu*u[,1]
  if (a==0) {
    # continuous treatment
    X <- linear.ps + rnorm(n=n,sd=sqrt(abs(linear.ps)/max(abs(linear.ps))))
    # rescale X
    X <- scale(X)[,1]
  } else {
    # binary treatment
    X <- rbinom(n,1,exp(linear.ps)/(1+exp(linear.ps)))
  }
  rm(linear.ps,g_raw)

  # outcome unaffected by treatment
  b_raw <- rep(0,p)
  b_raw[lC] <- 1 # true confounders
  b_raw[lP] <- 1 # predictors of outcome only
  linear.out <- (l %*% b_raw)[,1] + 2*nu*u[,2]
  y0 <- linear.out + rnorm(n=n, sd=max(abs(linear.out)))
  
  mydata <- data.frame("i"=1:n,"L"=l,"treat"=X,"Y"=y0)
  return(mydata)
}
