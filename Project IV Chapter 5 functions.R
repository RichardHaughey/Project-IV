# Functions used in Chapter 5: Bayesian inference for stochastic kinetic models

# 5.1 Inference under the LNA
# 5.1.1 Exact and fully observed system 

# Data simulation of exact Lokta-Volterra MJP via Gillespie's direct method
# Tmax = simulation end time
# cvec = rate constant vector
# init = initial state
LV_gillespie <- function(Tmax,cvec,init){
  S <- matrix(c(1,0,
                -1,1,
                0,-1),
              ncol = 2 , byrow = TRUE
  )
  x <- matrix(c(0,init), nrow=1, ncol = 3)
  colnames(x) <- c("time","X_1", "X_2")
  c1 <- cvec[1]
  c2 <- cvec[2]
  c3 <- cvec[3]
  xt <- init
  x1 <- init[1]
  x2 <- init[2]
  t <- 0
  while (t < Tmax) {
    h1 <- c1 * x1
    h2 <- c2 * x1 * x2
    h3 <- c3 * x2
    h0 <- h1 + h2 + h3
    if(h0 > 0){
      dt <- rexp(1, h0)
      j <- sample(c(1,2,3),1, prob = c(h1,h2,h3)/h0)
      xt <- xt + S[j,]
      x <- rbind(x,c(t+dt,xt))
      x1 <- xt[1]
      x2 <- xt[2]
      t <- t + dt
    }
    else{
      t <- Tmax+1
      x <- rbind(x,c(t,0,0))
    }
  }
  return(x[-length(x[,1]),])
}

# Collecting reactions at end of time interval
X_t <- LV_gillespie(Tmax,theta,init)
t <- ceiling(X_t[,1]/dt)+1
x <- matrix(NA, nrow=N, ncol = 2)
for (j in 1:length(t)) {
  x[t[j],] <- X_t[j,2:3] 
}
s <- sum(is.na(x[,1]))
while (s>0) {
  k <- which(is.na(x[,1]))
  x[k,1] <- x[k-1,1]
  x[k,2] <- x[k-1,2]
  s <- sum(is.na(x[,1]))
}



# LNA deterministic path ODE for Lokta-Volterra model
# time = time interval to integrate system over
# eta = initial value
# parameters = vector of rate parameters
eta_ode_system <- function(t, eta, parameters) {
  c1 <- parameters[1]
  c2 <- parameters[2]
  c3 <- parameters[3]
  
  detadt1 <- c1 * eta[1] - c2 * eta[1] * eta[2]
  detadt2 <- c2 * eta[1] * eta[2] - c3 * eta[2]
  
  return(list(c(detadt1, detadt2)))
}

# LNA residual mean ODE system for Lokta-Volterra model
# time = time interval to integrate system over
# m = initial value
# parameters = vector of rate parameters and initial state
m_ode_system <- function(t, m, parameters) {
  c1 <- parameters[1]
  c2 <- parameters[2]
  c3 <- parameters[3]
  x1 <- parameters[4]
  x2 <- parameters[5]
  
  m <- matrix(m, nrow = 2, ncol = 1)
  
  Ft <- matrix(c(c1-c2*x2,-c2*x1,c2*x2,c2*x1-c3), nrow = 2, byrow = TRUE)  # Jacobian matrix
  
  dmdt <- Ft %*% m 
  
  return(list(as.vector(dmdt)))
}

# LNA residual variance ODE system for Lokta-Volterra model
# time = time interval to integrate system over
# V = initial value
# parameters = vector of rate parameters and initial state
V_ode_system <- function(t, V, parameters) {
  c1 <- parameters[1]
  c2 <- parameters[2]
  c3 <- parameters[3]
  x1 <- parameters[4]
  x2 <- parameters[5]
  
  V <- matrix(V, nrow = 2, ncol = 2)
  
  Ft <- matrix(c(c1-c2*x2,-c2*x1,c2*x2,c2*x1-c3), nrow = 2, byrow = TRUE)  # Jacobian matrix
  beta <- matrix(c(c1*x1+c2*x1*x2, -c2*x1*x2,
                   -c2*x1*x2, c2*x1*x2+c3*x2),
                 nrow = 2, byrow = TRUE)
  
  dVdt <- Ft %*% V + V %*% t(Ft) + beta
  
  return(list(as.vector(dVdt)))
}

# Priors for log parameters of Lokta-Volterra model
LNA_theta1_prior <- function(theta1,b1){
  return(dnorm(theta1,0,b1,log = TRUE))
}

LNA_theta2_prior <- function(theta2,b2){
  return(dnorm(theta2,0,b2,log = TRUE))
}

LNA_theta3_prior <- function(theta3,b3){
  return(dnorm(theta3,0,b3,log = TRUE))
}

# LNA log likelihood function for Lokta-Volterra model
# theta. = vector of log(rate constants)
# x. = observations
# dt. = time step between observations
LNA_llike <- function(theta., x., dt.)
{
  n <- length(x.[,1])
  C <- exp(theta.)
  sum <- 0
  for(i in 2:n){
    t_eval <- c((i-2)*dt.,(i-1)*dt.)
    
    parameters <- c(C, x.[i-1,])

    eta_step <- ode(y = x.[i-1,], times = t_eval, func = eta_ode_system, parms = parameters)
    eta_t <- eta_step[2,-1]

    V_initial_state <- matrix(c(0.0, 0.0, 0.0, 0.0), nrow = 2, byrow = TRUE)

    V_step <- ode(y = as.vector(V_initial_state), times = t_eval, func = V_ode_system, parms = parameters)
    Vt <- matrix(V_step[2, -1], ncol = 2, byrow = TRUE)
    
    sum <- sum + dmvnorm(x.[i,], mean = eta_t, sigma = Vt, log = TRUE)
  }
  return(sum)
}

# log posterior for parameters of Lokta-Volterra model 
# theta = vector of log-parameters
# xt = observations
# b = vector of hyper parameters (b1,b2,b3)
# dt = time step
LNA_lpost=function(theta, xt, b, dt)
{
  lprior1 = LNA_theta1_prior(theta[1],b[1])
  lprior2 = LNA_theta2_prior(theta[2],b[2])
  lprior3 = LNA_theta3_prior(theta[3],b[3])
  llike = LNA_llike(theta,xt,dt)
  return(lprior1+lprior2+lprior3+llike)
}

# Random walk metropols algorithm for Lokta-volterra model with likelihood estimated under the LNA
# nruns = number of iterations 
# xt = observations
# dt = time step
# b = vector of hyper parameters (b1,b2,b3)
# psi = initial value of log-parameters
# tuneMat = innovation variance matrix
LNA_RWM <- function(nruns=400, xt=x, dt=1, b=c(1,1,1), psi=log(c(0.5,0.0025,0.3)),tuneMat=diag(c(0.0001,0.0001,0.0001)))
{
  n <- length(xt[,1])
  p <- length(psi)
  psi.mat <- matrix(0,ncol=p,nrow=nruns)
  psi.mat[1,] <- psi
  count <- 1
  
  for (k in 2:nruns)
  {
    can <- mvrnorm(1,psi,tuneMat)
    laprob <- LNA_lpost(can, xt, b, dt)-LNA_lpost(psi, xt, b, dt)
    u <- runif(1)
    if(log(u)<laprob)
    {
      psi <- can
      count <- count+1
    }
    psi.mat[k,] <- psi
  }
  print(count/nruns)
  return(exp(psi.mat))
}


# 5.1.2 Partially observed system 

# Simulate noisy observations of Lokta-Volterra model on regular time grid
X_t <- LV_gillespie(Tmax,theta,init)
t <- ceiling(X_t[,1]/dt)+1
x <- matrix(NA, nrow=N, ncol = 2)
for (j in 1:length(t)) {
  x[t[j],] <- X_t[j,2:3] 
}
s <- sum(is.na(x[,1]))
while (s>0) {
  k <- which(is.na(x[,1]))
  x[k,1] <- x[k-1,1]
  x[k,2] <- x[k-1,2]
  s <- sum(is.na(x[,1]))
}

P <- diag(1,2) 
sig <- 5

y <- t(P%*%t(x)) + round(matrix(rnorm(2*(Tmax+1),0,sig), ncol = 2))
y[which(y < 0)] <- 0


# forward filter log-likelihood estimation for parameters of Lokta-Volterra model
# yt. = noisy observations
# dt. = time step
# tehta. = vector of log-parameters
# a. = initial state prior mean
# B. = initial state prior variance
# P. = observation matrix 
# Sigma. = observation error variance matrix
fwd_filter <- function(yt., dt., theta., a., B., P., Sigma.){
  n <- length(yt.[,1])
  p <- length(yt.[1,])
  sum <- dmvnorm(yt.[1,], t(P.)%*%a., t(P.)%*%B.%*%P. + Sigma., log = TRUE)
  
  a <- a. + B.%*%P. %*% solve(t(P.)%*%B.%*%P. + Sigma.) %*% (yt.[1,]- t(P.)%*%a.)
  B <- B. - B.%*%P. %*% solve(t(P.)%*%B.%*%P. + Sigma.) %*% t(P.)%*%B.
  
  C <- exp(theta.)
  
  for (i in 2:n) {
    t_eval <- c((i-2)*dt.,(i-1)*dt.)
    
    parameters <- c(C, a)
    
    eta_step <- ode(y = a, times = t_eval, func = eta_ode_system, parms = parameters)
    eta_t <- eta_step[2,-1]
    
    V_step <- ode(y = as.vector(B), times = t_eval, func = V_ode_system, parms = parameters)
    Vt <- matrix(V_step[2, -1], ncol = 2, byrow = TRUE)
    
    a <- eta_t + Vt%*%P %*% solve(t(P.)%*%Vt%*%P. + Sigma.) %*% (yt.[i,]- t(P.)%*%eta_t)
    B <- Vt - Vt%*%P %*% solve(t(P.)%*%Vt%*%P. + Sigma.) %*% t(P.)%*%Vt
    
    sum <- sum + dmvnorm(yt.[i,], t(P.)%*%eta_t, t(P.)%*%Vt%*%P. + Sigma., log = TRUE)
  }
  return(sum)
}

# log posterior for parameters of Lokta-Volterra model with 
# likelihood estimate coming from the forward filter
# theta = vector of log-parameters
# st = process observations
# b = vector of hyper parameters (b1,b2)
# dt = time step
ff_lpost=function(yt_, dt_, theta_, a_, B_, P_, Sigma_, b_)
{
  lprior1 = theta1_prior(theta_[1],b_[1])
  lprior2 = theta2_prior(theta_[2],b_[2])
  lprior3 = theta3_prior(theta_[3],b_[3])
  llike = fwd_filter(yt_,dt_,theta_,a_,B_,P_,Sigma_)
  return(lprior1+lprior2+lprior3+llike)
}


# Random walk metropols algorithm for Lokta-volterra model with likelihood estimated via forward filter
# nruns = number of iterations 
# yt = observations
# dt = time step
# b = vector of hyper parameters (b1,b2,b3)
# psi = initial value of log-parameters
# tuneMat = innovation variance matrix
# Sigma = observation error variance matrix
# P = observation matrix
# a = initial state prior mean
# B = initial state prior variance
RWM_LNA.partial.obs <- function(nruns=400, yt=y, dt=1, b=c(10,10,10), psi=log(c(0.5,0.0025,0.3)),
                                tuneMat=diag(c(0.0001,0.0001,0.0001)), 
                                Sigma=diag(sig,2),P=diag(1,2),a=c(71, 79),B=diag(0,2))
{
  n <- length(yt[,1])
  p <- length(psi)
  psi.mat <- matrix(0,ncol=p,nrow=nruns)
  psi.mat[1,] <- psi
  count <- 1
  
  for (k in 2:nruns)
  {
    can <- mvrnorm(1,psi,tuneMat)
    laprob <- ff_lpost(yt,dt,can,a,B,P,Sigma,b)-ff_lpost(yt,dt,psi,a,B,P,Sigma,b)
    u <- runif(1)
    if(log(u)<laprob)
    {
      psi <- can
      count <- count+1
    }
    psi.mat[k,] <- psi
  }
  print(count/nruns)
  return(exp(psi.mat))
}


# 5.3 Comparison
# 5.3.1 Death model

# Simulation of death model data
# x0 = initial value
# theta = rate constant
# dt = time step
# Tmax = simulation end time
death_sim <- function(x0,theta,dt,Tmax){
  N <- Tmax/dt+1
  x <- rep(NA,N)
  x[1] <- x0
  for (i in 2:N) {
    x[i] <- rbinom(1,x[i-1],exp(-theta*dt))
  }
  return(x)
}


## Death model analysis with analytic likelihood
death_An_llike <- function(theta, x, dt)
{
  n <- length(x)
  C <- exp(theta)
  sum <- 0
  for(i in 2:n){
    sum <- sum + dbinom(x[i], x[i-1],exp(-C*dt), log = TRUE)
  }
  return(sum)
}

# analytic posterior
death_An_lpost=function(theta, xt, b, dt)
{
  lprior = death_theta_prior(theta,b)
  llike = death_An_llike(theta,xt,dt)
  return(lprior+llike)
}

# Metropolis with analytic likelihood
death_RWM_An <- function(nruns=400, xt=x1, dt=0.1, b=10, psi=log(0.5),tuneMat=0.3)
{
  psi.mat <- rep(NA,nruns)
  psi.mat[1] <- psi
  count <- 1
  
  for (k in 2:nruns)
  {
    can <- rnorm(1,psi,tuneMat)
    laprob <- death_An_lpost(can, xt, b, dt)-death_An_lpost(psi, xt, b, dt)
    u <- runif(1)
    if(log(u)<laprob)
    {
      psi <- can
      count <- count+1
    }
    psi.mat[k] <- psi
  }
  print(count/nruns)
  return(exp(psi.mat))
}


## Death model analysis with LNA

# theta prior for death model
death_theta_prior <- function(theta,b){
return(dnorm(theta,0,b,log = TRUE))
}

# log-likelihood under LNA for death model
death_LNA_llike <- function(theta, x, dt)
{
  n <- length(x)
  C <- exp(theta)
  sum <- 0
  for(i in 2:n){
    sum <- sum + dnorm(x[i], x[i-1]*exp(-C*dt),
                       sqrt(x[i-1]*exp(-C*dt)*(1-exp(-C*dt))),
                       log = TRUE)
  }
  return(sum)
}

# log-parameter posterior for death model
death_LNA_lpost=function(theta, xt, b, dt)
{
  lprior = death_theta_prior(theta,b)
  llike = death_LNA_llike(theta,xt,dt)
  return(lprior+llike)
}

# Random walk metropolis with LNA likelihood for death model
death_RWM_LNA <- function(nruns=400, xt=x1, dt=0.1, b=10, psi=log(0.5),tuneMat=0.3)
{
  # initialise
  psi.mat <- rep(NA,nruns)
  psi.mat[1] <- psi
  psi.lpost <- death_LNA_lpost(psi, xt, b, dt)
  count <- 1
  
  for (k in 2:nruns)
  {
    can <- rnorm(1,psi,tuneMat)
    can.lpost <- death_LNA_lpost(can, xt, b, dt)
    laprob <- can.lpost - psi.lpost
    u <- runif(1)
    if(log(u)<laprob)
    {
      psi <- can
      psi.lpost <- can.lpost
      count <- count+1
    }
    psi.mat[k] <- psi
  }
  print(count/nruns)
  return(exp(psi.mat))
}


## Death model analysis with PMMH
# state forward propagation
death_xtSim <- function(x,theta,dt){
  x_sim <- rbinom(1,x,exp(-theta*dt))
  return(x_sim)
}

# Function to implement particle filter to estimate log-likelihood for Death model
# N = number of particles
death_BPF=function(N,ltheta,ydata,dt)
{
  endT=length(ydata)
  theta <- exp(ltheta)
  x=rep(NA,N) # Samples of current state x_t
  p_y_given_theta <- 0
  
  # Loop for times t=1,...,T
  for(i in 2:endT)
  {
    for(j in 1:N)
    {
      dx <- death_xtSim(ydata[i-1],theta,dt) # Propagate
      x[j]= dx 
    }
    hit.prop <- sum(x==ydata[i])
    p_y_given_theta <- p_y_given_theta + log(hit.prop) - log(N)
    if(hit.prop == 0){break}
  }
  return(p_y_given_theta)
}

# log-posterior for rate constants for death model with likelihood estimate via PMMH
death_BPF_lpost=function(theta, xt, b, dt, N)
{
  lprior = death_theta_prior(theta,b)
  llike = death_BPF(N,theta,xt,dt)
  return(lprior+llike)
}

# Random walk Metropolis algorithm for death model with likelihood estimated via PMMH
death_RWM_PMMH <- function(nruns=400, xt=x2, dt=0.1, b=10, psi=log(0.5),tuneMat=0.1, N=200)
{
  # initialise
  psi.mat <- rep(NA,nruns)
  psi.mat[1] <- psi
  psi.lpost <- -Inf
  while (psi.lpost == -Inf) {
    psi.lpost <- death_BPF_lpost(psi, xt, b, dt, N)
  }
  count <- 1
  
  for (k in 1:nruns)
  {
    can <- rnorm(1,psi,tuneMat)
    can.lpost <- death_BPF_lpost(can, xt, b, dt, N)
    laprob <- can.lpost - psi.lpost
    u <- runif(1)
    if(log(u)<laprob)
    {
      psi <- can
      psi.lpost <- can.lpost
      count <- count+1
    }
    psi.mat[k] <- psi
  }
  print(count/nruns)
  return(exp(psi.mat))
}


# 5.3.1 Epidemic model

# Epidemic model data simulation via Gillespie's direct method
sir_gillespie <- function(Tmax, cvec, init) {
  S <- matrix(c(-1, 1,
                0, -1),
              ncol = 2, byrow = TRUE)
  
  x <- matrix(0, nrow = 100*Tmax, ncol = 3)
  x[1, ] <- c(0, init)
  colnames(x) <- c("time", "S", "I")
  
  c1 <- cvec[1]
  c2 <- cvec[2]
  
  x1 <- init[1]
  x2 <- init[2]
  
  t <- 0
  idx <- 1
  
  while (t < Tmax) {
    h1 <- c1 * x1 * x2
    h2 <- c2 * x2
    h0 <- h1 + h2
    
    if (h0 > 0) {
      dt <- rexp(1, h0)
      j <- sample(c(1, 2), 1, prob = c(h1, h2) / h0)
      xt <- c(x1, x2) + S[j, ]
      idx <- idx + 1
      
      if (idx > nrow(x)) {
        x <- rbind(x, matrix(0, nrow = 100, ncol = 3))
      }
      
      x[idx, ] <- c(t + dt, xt)
      x1 <- xt[1]
      x2 <- xt[2]
      t <- t + dt
    } else {
      t <- Tmax + 1
      idx <- idx + 1
      x[idx, ] <- c(t, 0, 0)
    }
  }
  
  return(x[1:(idx-1), ])
}

X_t <- sir_gillespie(Tmax,cvec,init)
t <- ceiling(X_t[,1]/dt)+1
x1 <- matrix(NA, nrow=steps+1, ncol = 2)
for (j in 1:length(t)) {
  x1[t[j],] <- X_t[j,2:3] 
}
s <- sum(is.na(x[,1]))
while (s>0) {
  k <- which(is.na(x[,1]))
  x1[k,1] <- x1[k-1,1]
  x1[k,2] <- x1[k-1,2]
  s <- sum(is.na(x1[,1]))
}

x1 <- cbind(x1,120-x1[,1]-x1[,2]) # removed population
r1 <- x1[2:(steps+1),3]-x1[1:steps,3] # change in removed population

y1 <- rbinom(length(r1),r1,0.8) # partial observations


# rate constant priors for Epidemic model
sir_sir_theta1_prior <- function(theta1,b1){
  return(dnorm(theta1,0,b1,log = TRUE))
}

sir_sir_theta2_prior <- function(theta2,b2){
  return(dnorm(theta2,0,b2,log = TRUE))
}

# LNA deterministic path ODE for SIR model
sir_eta_ode_system <- function(t, eta, parameters) {
  c1 <- parameters[1]
  c2 <- parameters[2]
  
  detadt1 <- - c1 * eta[1] * eta[2]
  detadt2 <- c1 * eta[1] * eta[2] - c2 * eta[2]
  
  return(list(c(detadt1, detadt2)))
}

# LNA fundemental matrix ODE system for SIR model
sir_G_ode_system <- function(t, G, parameters) {
  c1 <- parameters[1]
  c2 <- parameters[2]
  x1 <- parameters[3]
  x2 <- parameters[4]
  
  G <- matrix(G, nrow = 2, ncol = 2)
  
  Ft <- matrix(c(-c1*x2,-c1*x1,c1*x2,c1*x1-c2), nrow = 2, byrow = TRUE)  # Jacobian matrix
  
  dGdt <- Ft %*% G 
  
  return(list(as.vector(dGdt)))
}

# LNA residual variance ODE system for SIR model
sir_V_ode_system <- function(t, V, parameters) {
  c1 <- parameters[1]
  c2 <- parameters[2]
  x1 <- parameters[3]
  x2 <- parameters[4]
  
  V <- matrix(V, nrow = 2, ncol = 2)
  
  Ft <- matrix(c(-c1*x2,-c1*x1,c1*x2,c1*x1-c2), nrow = 2, byrow = TRUE)  # Jacobian matrix
  beta <- matrix(c(c1*x1*x2, -c1*x1*x2,
                   -c1*x1*x2, c1*x1*x2+c2*x2),
                 nrow = 2, byrow = TRUE)
  
  dVdt <- Ft %*% V + V %*% t(Ft) + beta
  
  return(list(as.vector(dVdt)))
}

# forward filter likelihood for SIR model with binomial observation model
# lambda = binomial observation model success rate parameter
fwd_filter_sir <- function(yt, dt, ltheta, a, B, P, lambda){
  suppressWarnings({
    n <- length(yt)
    sum <- 0
    
    C <- exp(ltheta)
    
    for (i in 1:n) {
      t_eval <- c((i-1)*dt,(i)*dt)
      
      parameters <- c(C, a)
      
      eta_step <- ode(y = a, times = t_eval, func = sir_eta_ode_system, parms = parameters)
      eta_t <- eta_step[2,-1];
      
      G_step <- ode(y = c(1,0,0,1), times = t_eval, func = sir_G_ode_system, parms = parameters)
      Gt <- matrix(G_step[2, -1], ncol = 2)
      
      V_step <- ode(y = as.vector(B), times = t_eval, func = sir_V_ode_system, parms = parameters)
      Vt <- matrix(V_step[2, -1], ncol = 2)
      
      E_dXt <-  eta_t - a
      
      V_dXt <- Vt + B - B%*%t(Gt) - Gt%*%B
      
      V_obs <- lambda*(1-lambda)*t(P)%*%E_dXt
      
      if(sum(is.na(lambda^2*t(P)%*%V_dXt%*%P + V_obs)) > 0){
        sum <- -Inf
        break
      }
      if(lambda^2*t(P)%*%V_dXt%*%P + V_obs < 0){
        sum <- -Inf
        break
      }
      
      sum <- sum + dnorm(yt[i],lambda*t(P)%*%E_dXt,sqrt(lambda^2*t(P)%*%V_dXt%*%P + V_obs),log = TRUE)
      
      C_XY <- lambda*(Vt-Gt%*%B)%*%P
      
      a <- eta_t + C_XY %*% solve(lambda^2*t(P)%*%V_dXt%*%P + V_obs) %*% (yt[i]- lambda*t(P)%*%E_dXt)
      B <- Vt - C_XY %*% solve(lambda^2*t(P)%*%V_dXt%*%P + V_obs) %*% t(C_XY)
    }
    return(sum)
  })
}

# log-posterior
sir_ff_lpost=function(yt, dt, theta, a, B, P, lambda, b)
{
  lprior1 = sir_theta1_prior(theta[1],b[1])
  lprior2 = sir_theta2_prior(theta[2],b[2])
  llike = fwd_filter_sir(yt,dt,theta,a,B,P,lambda)
  return(lprior1+lprior2+llike)
}

# Random walk Metropolis for SIR modlel with binomial observations
sir_RWM_FF <- function(nruns=400, yt=y1, dt=10, b=c(100,100), psi=log(c(0.00091, 0.082)),
                   tuneMat=diag(c(0.05,0.05)),
                   lambda=0.8,P=c(-1,-1),a=c(119,1),B=diag(0,2))
{
  p <- length(psi)
  psi.mat <- matrix(0,ncol=p,nrow=nruns)
  psi.mat[1,] <- psi
  psi.lpost <- -Inf
  while (psi.lpost == -Inf) {
    psi.lpost <- sir_ff_lpost(yt,dt,psi,a,B,P,lambda,b)
  }
  count <- 1
  
  for (k in 2:nruns)
  {
    can <- mvrnorm(1,psi,tuneMat)
    can.lpost <- sir_ff_lpost(yt,dt,can,a,B,P,lambda,b)
    laprob <- can.lpost-psi.lpost
    u <- runif(1)
    if(log(u)<laprob)
    {
      psi <- can
      psi.lpost <- can.lpost
      count <- count+1
    }
    psi.mat[k,] <- psi
  }
  print(count/nruns)
  return(exp(psi.mat))
}


# 5.4 Real world example: Eyam plague data

# data 
St <- c(254, 235, 201, 153, 121, 110, 97, 83)
It <- c(7, 14, 22, 29, 20, 8, 8, 0)
ydata <- matrix(c(St,It),ncol=2,byrow = FALSE)
t_obs <- c(0,0.5,1,1.5,2,2.5,3,4)

# priors
sir_theta1_prior <- function(theta1,b1){
  return(dnorm(theta1,0,b1,log = TRUE))
}
sir_theta2_prior <- function(theta2,b2){
  return(dnorm(theta2,0,b2,log = TRUE))
}
# prior for observation noise variance
sir_lsigma_prior <- function(lsigma,b3){
  return(dgamma(exp(lsigma),1,2,log = TRUE))
}

# Exact observation likelihood estimate for SIR model
LNA_llike <- function(theta, yt, times)
{
  n <- length(yt[,1])
  C <- exp(theta)
  sum <- 0
  for(i in 2:n){
    t_eval <- c(times[i-1],times[i])
    
    parameters <- c(C, yt[i-1,])
    
    eta_step <- ode(y = yt[i-1,], times = t_eval, func = eta_ode_system, parms = parameters)
    eta_t <- eta_step[2,-1]
    
    V_step <- ode(y = rep(0,4), times = t_eval, func = V_ode_system, parms = parameters)
    Vt <- matrix(V_step[2, -1], ncol = 2, byrow = TRUE)
    
    sum <- sum + dmvnorm(yt[i,], eta_t, Vt, log = TRUE)
  }
  return(sum)
}

# forward filter likelihood estimate for noisy SIR observations
fwd_filter <- function(yt, times, theta, a, B, P, Sigma){
  
  if(Sigma[1,1] < 10^(-300)){return(LNA_llike(theta, yt, times))}
  else{
    n <- length(times)
    sum <- dmvnorm(yt[1,], t(P)%*%a, t(P)%*%B%*%P + Sigma, log = TRUE)
    
    a <- a + B%*%P %*% solve(t(P)%*%B%*%P + Sigma) %*% (yt[1,]- t(P)%*%a)
    B <- B - B%*%P %*% solve(t(P)%*%B%*%P + Sigma) %*% t(P)%*%B
    
    C <- exp(theta)
    
    for (i in 2:n) {
      t_eval <- c(times[i-1],times[i])
      
      parameters <- c(C, a)
      
      eta_step <- ode(y = a, times = t_eval, func = eta_ode_system, parms = parameters)
      eta_t <- eta_step[2,-1]
      
      V_step <- ode(y = as.vector(B), times = t_eval, func = V_ode_system, parms = parameters)
      Vt <- matrix(V_step[2, -1], ncol = 2, byrow = TRUE)
      
      a <- eta_t + Vt%*%P %*% solve(t(P)%*%Vt%*%P + Sigma) %*% (yt[i,]- t(P)%*%eta_t)
      B <- Vt - Vt%*%P %*% solve(t(P)%*%Vt%*%P + Sigma) %*% t(P)%*%Vt
      
      sum <- sum + dmvnorm(yt[i,], t(P)%*%eta_t, t(P)%*%Vt%*%P + Sigma, log = TRUE)
    }
    return(sum)
  }
}

# log-posterior with forward filter log-likelihood
ff_lpost=function(yt, t, theta, a, B, P, b)
{
  lprior1 = sir_theta1_prior(theta[1],b[1])
  lprior2 = sir_theta2_prior(theta[2],b[2])
  lprior3 = sir_lsigma_prior(theta[3],b[3])
  llike = fwd_filter(yt,t,theta[1:2],a,B,P,diag(exp(theta[3]),2))
  return(lprior1+lprior2+lprior3+llike)
}

# log-posterior with exact observation log-likelihood estimate under LNA
LNA_lpost=function(yt, t, theta, b)
{
  lprior1 = sir_theta1_prior(theta[1],b[1])
  lprior2 = sir_theta2_prior(theta[2],b[2])
  llike = LNA_llike(theta,yt,t)
  return(lprior1+lprior2+llike)
}

# Random walk Metropolis for SIR model with likelihood estimate via forward filter
RWM_FF <- function(nruns=400, yt=ydata, times=t_obs, b=c(100,100,5), psi=log(c(0.0196, 3.204, 1)),
                   tuneMat=diag(c(0.05,0.05,0.05)),
                   P=diag(1,2),a=c(254,7),B=diag(0,2))
{
  p <- length(psi)
  psi.mat <- matrix(0,ncol=p,nrow=nruns)
  psi.mat[1,] <- psi
  psi.lpost <- ff_lpost(yt,times,psi,a,B,P,b)
  count <- 1
  
  for (k in 2:nruns)
  {
    can <- mvrnorm(1,psi,tuneMat)
    can.lpost <- ff_lpost(yt,times,can,a,B,P,b)
    laprob <- can.lpost-psi.lpost
    u <- runif(1)
    if(log(u)<laprob)
    {
      psi <- can
      psi.lpost <- can.lpost
      count <- count+1
    }
    psi.mat[k,] <- psi
  }
  print(count/nruns)
  return(exp(psi.mat))
}

# Random walk Metropolis for SIR model with likelihood estimate via LNA
MwG_LNA <- function(nruns=400, yt=ydata, times=t_obs, b=c(100,100), psi=log(c(0.0196, 3.204)),
                    tuneMat=diag(c(0.05,0.05)))
{
  p <- length(psi)
  psi.mat <- matrix(0,ncol=p,nrow=nruns)
  psi.mat[1,] <- psi
  psi.lpost <- LNA_lpost(yt,times,psi,b)
  count <- 1
  
  for (k in 2:nruns)
  {
    can <- mvrnorm(1,psi,tuneMat)
    can.lpost <- LNA_lpost(yt,times,can,b)
    laprob <- can.lpost-psi.lpost
    u <- runif(1)
    if(log(u)<laprob)
    {
      psi <- can
      psi.lpost <- can.lpost
      count <- count+1
    }
    psi.mat[k,] <- psi
  }
  print(count/nruns)
  return(exp(psi.mat))
}