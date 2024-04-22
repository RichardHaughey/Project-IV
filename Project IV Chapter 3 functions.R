# Functions used in Chapter 3: Bayesian Inference

# required packages
library(mvtnorm)
library(MASS)

# 3.2 Euler-Maruyama
# 3.2.1  Example: Geometric Brownian motion

# Priors for log parameters of geometric Brownian motion
# theta1 = log(mu)
# b1 = standard deviation of normal distribution
GeoBM_theta1_prior <- function(theta1,b1){
  return(dnorm(theta1,0,b1,log = TRUE))
}

# theta1 = log(sigma)
# b2 = standard deviation of normal distribution
GeoBM_theta2_prior <- function(theta2,b2){
  return(dnorm(theta2,0,b2,log = TRUE))
}

# Geometric Brownian motion analytic log likelihood
# theta = vector of log-parameters
# st = process observations
# dt = time step
GeoBM_loglike <- function(theta, st, dt){
  n <- length(st)
  diff <- log(st[-1]) - log(st[-n])
  return(sum(dnorm(diff ,
                   (exp(theta[1])-exp(2*theta[2])/2)*dt ,
                   exp(theta[2])*sqrt(dt) , 
                   log = TRUE )))
}

# log posterior for parameters of geometric Brownian motion model 
# theta = vector of log-parameters
# st = process observations
# b = vector of hyper parameters (b1,b2)
# dt = time step
GeoBM_lpost=function(theta, st, b, dt=1)
{
  lprior1 = GeoBM_theta1_prior(theta[1],b[1])
  lprior2 = GeoBM_theta2_prior(theta[2],b[2])
  llike = GeoBM_loglike(theta,st,dt)
  return(lprior1+lprior2+llike)
}

# Metropolis-Hastings algorithm for geometric Brownian motion
# N = number of runs
# st = process observations
# dt = time step
# b = vector of hyper parameters (b1,b2)
# psi = initial value of log-parameters
# Sig = innovation variance matrix
GeoBM_MH=function(N, st, dt, b, psi, Sig) 
{  
  mat = matrix(rep(0,2*N),ncol = 2)
  mat[1,] = psi
  for (i in 2:N) 
  {
    can = rmvnorm(1,psi,Sig)
    laprob = GeoBM_lpost(can,st,b,dt)-GeoBM_lpost(psi,st,b,dt)
    u = runif(1)
    if (log(u) < laprob){ 
      psi = can
    } 
    mat[i,] = psi
  }
  return(mat)
}


# Geometric Brownian motion log likelihood under Euler-Maruyama approximation
# theta = vector of log-parameters
# st = process observations
# dt = time step
GeoBM_EM_llike <- function(theta, st, dt)
{
  n <- length(st)
  mu <- exp(theta[1])
  sig <- exp(theta[2])
  sum <- 0
  for(i in 2:n){
    st_dt <- st[i-1]
    sum <- sum + dnorm(st[i], st_dt + st_dt*mu*dt, sig*st_dt*sqrt(dt), log = TRUE)
  }
  return(sum)
}

# log posterior for parameters of geometric Brownian motion model under Euler-Maruyama approximation
# theta = vector of log-parameters
# st = process observations
# b = vector of hyper parameters (b1,b2)
# dt = time step
GeoBM_EM_lpost=function(theta, st, b, dt)
{
  lprior1 = GeoBM_theta1_prior(theta[1],b[1])
  lprior2 = GeoBM_theta2_prior(theta[2],b[2])
  llike = GeoBM_EM_llike(theta,st,dt)
  return(lprior1+lprior2+llike)
}

# Metropolis-Hastings algorithm for geometric Brownian motion under  Euler-Maruyama approximation
# N = number of runs
# st = process observations
# dt = time step
# b = vector of hyper parameters (b1,b2)
# psi = initial value of log-parameters
# Sig = innovation variance matrix
GeoBM_EM_MH=function(N, st, dt, b, psi, Sig) 
{  
  mat = matrix(rep(0,2*N),ncol = 2)
  mat[1,] = psi
  for (i in 2:N) 
  {
    can = rmvnorm(1,psi,Sig)
    laprob = GeoBM_EM_lpost(can,st,b,dt)-GeoBM_EM_lpost(psi,st,b,dt)
    u = runif(1)
    if (log(u) < laprob){ 
      psi = can
    } 
    mat[i,] = psi
  }
  return(mat)
}



# 3.3 Data augmentation
# 3.3.1 Single site update

# Single site update for geometric Brownian motion
# nruns = number of iterations 
# x_obs = process observations
# dt = time step
# m = number of unobserved data points to be simulated +1
# b = vector of hyper parameters (b1,b2)
# psi = initial value of log-parameters
# Sig = innovation variance matrix
ssu_MwG <- function(nruns, x_obs, dt, m, b, psi, Sig){
  
  # initialise
  dtau <- dt/m
  theta.mat = matrix(rep(0,2*nruns),ncol = 2)
  theta.mat[1,] <-  psi
  n <- length(x_obs)
  N <- m*(n-1)+1
  x <- rep(NA,N)
  for(i in 1:(n-1)){
    for(j in 1:m){
      x[(i-1)*m + j] <- (x_obs[i]*(m-(j-1))+x_obs[i+1]*(j-1))/m
    }
  }
  x[N] <- x_obs[n]
  x.mat <- matrix(NA, ncol = N, nrow = nruns)
  x.mat[1,] <- x
  
  accept.theta <- 1
  accept.x <- rep(1,(n-1)*(m-1))
  
  # updating
  
  for (k in 2:nruns) 
  {
    # Updating theta
    theta.can = rmvnorm(1,psi,Sig)
    theta.laprob = EM_lpost(theta.can,x,b,dtau)-EM_lpost(psi,x,b,dtau)
    u = runif(1)
    if (log(u) < theta.laprob){ 
      psi = theta.can
      accept.theta <- accept.theta+1
    } 
    theta.mat[k,] = psi
    
    # Updating x
    for(i in 1:(n-1)){
      for(j in 2:m){
        mu <- (x[(i-1)*m+j-1]+x[(i-1)*m+j+1])/2
        sd <- exp(psi[2])*x[(i-1)*m+j-1]*sqrt(dtau/2)
        x.can <- rnorm(1,mu,sd)
        
        sd.prev <- exp(psi[2])*x[(i-1)*m+j-1]*sqrt(dtau)
        sd.curr <- exp(psi[2])*x[(i-1)*m+j]*sqrt(dtau)
        sd.can <- exp(psi[2])*x.can*sqrt(dtau)
        
        alpha.prev <- x[(i-1)*m+j-1] + exp(psi[1])*x[(i-1)*m+j-1]*dtau
        alpha.curr <- x[(i-1)*m+j] + exp(psi[1])*x[(i-1)*m+j]*dtau
        alpha.can <- x.can + exp(psi[1])*x.can*dtau
        
        x.laprob1 <- (dnorm(x.can,alpha.prev,sd.prev,log=TRUE) +
                        dnorm(x[(i-1)*m+j+1],alpha.can,sd.can,log=TRUE) -
                        dnorm(x[(i-1)*m+j],alpha.prev,sd.prev,log=TRUE) -
                        dnorm(x[(i-1)*m+j+1],alpha.curr,sd.curr,log=TRUE) )
        x.laprob2 <- dnorm(x[(i-1)*m+j],mu,sd,log=TRUE)-dnorm(x.can,mu,sd,log=TRUE)
        x.laprob <- x.laprob1 + x.laprob2
        u <- runif(1)
        if (log(u) < x.laprob){ 
          x[(i-1)*m+j] <- x.can
          accept.x[(i-1)*(m-1)+j] <- accept.x+1
        } 
      }
    }
    x.mat[k,] <- x
    
  }
  print(paste("theta acceptance rate =",accept.theta/nruns))
  print(paste("x acceptance rate =",accept.x/nruns))
  return(list(x = x.mat, theta = theta.mat))
}


# 3.3.2 Modified diffusion bridge

# Modified diffusion bridge proposal function between two observations of geometric Brownian motion
# T0,T1 = observations tow simulated augmented data between
# m = number of unobserved data points to be simulated +1
# dt = time step
# sigma = volatility parameter
bridge_prop <- function(T0, T1, m, dt, sigma){
  dtau <- dt/m
  bridge <- rep(NA,m+1)
  bridge[1] <- T0
  bridge[m+1] <- T1
  for (i in 2:m) {
    prop <- rnorm(1,bridge[i-1] + (T1-bridge[i-1])/(m+2-i), sigma*bridge[i-1]*sqrt(dtau*(m+1-i)/(m+2-i)))
    if(prop<=0.001){
      prop <- 0.001
    }
    bridge[i] <- prop
  }
  return(bridge)
}

# MDB likelihood under Euler approximation of geometric Brownian motion
# theta = log-parameter vector
# st = observations of process
# dt = time step
bridge_llike <- function(theta, st, dt)
{
  n <- length(st)
  m <- n-1
  T1 <- st[n]
  mu <- exp(theta[1])
  sig <- exp(theta[2])
  sum <- 0
  for(i in 2:(n-1))
  {
    sum <- sum + dnorm(st[i],st[i-1] + (T1-st[i-1])/(m+2-i), sig*st[i-1]*sqrt(dt*(m+1-i)/(m+2-i)), log = TRUE)
  }
  return(sum)
}


# Acceptance probability of MDB proposal for geometric Brownian motion model
# theta = log-parameter vector
# can = (candidate) bridge proposal
# curr = current bridge value
# dt = time step
# m = number of unobserved data points to be simulated +1
# T1 = observation at end of interval
A <- function(theta, can, curr, dt, m, T1)
{
  dtau <- dt/m
  laprob1 <- EM_llike(theta,can,dtau)-EM_llike(theta,curr,dtau)
  laprob2 <- bridge_llike(theta,curr,dtau)-bridge_llike(theta,can,dtau)
  return(laprob1+laprob2)
}

# Metropolis-Hastings with componentwise transitions for Geometric BM with 
# MDB proposal for unobserved data
# nruns = number of iterations 
# st = process observations
# dt = time step
# m = number of unobserved data points to be simulated +1
# b = vector of hyper parameters (b1,b2)
# psi = initial value of log-parameters
# tuneMat = innovation variance matrix
MDB_MwG <- function(nruns, st, dt, m, b, psi,tuneMat)
{
  
  # initialise
  dtau <- dt/m
  n <- length(st)
  N <- m*(n-1)+1
  stau <- rep(NA,N)  # holds entire path incl. between obs
  for(i in 1:(n-1)){
    for(j in 1:m){
      stau[(i-1)*m + j] <- (st[i]*(m-(j-1))+st[i+1]*(j-1))/m
    }
  }
  stau[N] <- st[n]
  s.mat <- matrix(NA, ncol = N, nrow = nruns)
  s.mat[1,] <- stau
  p<-length(psi)
  psi.mat <- matrix(0,ncol=p,nrow=nruns)
  psi.mat[1,] <- psi
  count <- 1
  count.s <- rep(1,(n-1))
  
  # updating
  
  for (k in 2:nruns)
  {
    #parameter updating (conditional on current path)
    can <- mvrnorm(1,psi,tuneMat)
    laprob <- EM_lpost(can, stau, b, dtau)-EM_lpost(psi, stau, b, dtau)
    u <- runif(1)
    if(log(u)<laprob)
    {
      psi <- can
      count <- count+1
    }
    psi.mat[k,] <- psi
    
    #segment updating (conditional on current params)
    for(i in 1:(n-1))
    {
      s.can <- bridge_prop(st[i],st[i+1],m,dt,exp(psi[2]))
      s.curr <- stau[seq((i-1)*m + 1,i*m+1,length.out=(m+1))]
      s.laprob <- A(psi,s.can,s.curr,dt,m,st[i+1])
      u <- runif(1)
      if (log(u) < s.laprob)
      {
        stau[seq((i-1)*m + 1,i*m+1,length.out=(m+1))] <- s.can
        count.s[i] <- count.s[i]+1
      }
    }
    s.mat[k,] <- stau
  }
  
  print(count/nruns) 
  print(count.s/nruns) 
  return(list(psi.mat,s.mat))
}

