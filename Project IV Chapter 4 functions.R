# Functions used in Chapter 4: Multivariate SDEs

# required packages
library(mvtnorm)
library(MASS)

# 4 Multivariate SDEs
# 4.5 Simulation study: Lokta Volterra

# Numerical integration for Lokta-Volterra ODE model
# init = initial state
# theta = rate constants
# Tmax = simulation end time
# dt = time step
lv_ode <- function(init,theta, Tmax, dt){
  n <- Tmax/dt + 1
  lv <- matrix(NA, nrow = n, ncol = 2)
  lv[1,] <- init
  
  c1 <- theta[1]
  c2 <- theta[2]
  c3 <- theta[3]
  
  pop <- sum(init)
  for(j in 2:n){
    x <- lv[j-1,1]
    y <- lv[j-1,2]
    lv[j,1] <- x + (c1*x  - c2*x*y)*dt
    lv[j,2] <- y + (c2*x*y - c3*y)*dt
  }
  return(lv)
}

# Simulation of Lokta-Volterra MJP via Gillespie's direct method
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

# Lokta-Volterra simulation via LNA without restart
# Tmax = simulation end time
# dt = time step
# cvec = rate constant vector
# init = initial state
LV_LNA_no_restart <- function(Tmax,dt,cvec,init){
  N <- Tmax/dt + 1
  x <- matrix(NA, nrow=N, ncol = 2)
  colnames(x) <- c("X_1", "X_2")
  xt <- init
  x[1,] <- xt
  eta_t <- xt
  m_t <- xt-eta_t
  V_t <- c(0,0,0,0)
  for (i in 2:N) {
    t_eval <- c((i-2)*dt,(i-1)*dt)
    parameters <- c(c1 = cvec[1],c2 = cvec[2],c3 = cvec[3],x1 = eta_t[1],x2 = eta_t[2])

    eta_step <- ode(y = eta_t, times = t_eval, func = eta_ode_system, parms = parameters)
    eta_dt <- eta_step[2,-1]

    m_step <- ode(y = m_t, times = t_eval, func = m_ode_system, parms = parameters)
    m_dt <- matrix(m_step[2, -1], ncol = 1, byrow = TRUE)

    V_step <- ode(y = V_t, times = t_eval, func = V_ode_system, parms = parameters)
    V_dt <- matrix(V_step[2, -1], ncol = 2, byrow = TRUE)
    
    xt <- mvrnorm(1,eta_dt+m_dt, V_dt)
    xt[which(xt < 0)] <- 0
    x[i,] <- xt
    eta_t <- eta_dt
    m_t <- xt-eta_t
    V_t <- c(0,0,0,0)
  }
  return(x)
}

# Lokta-Volterra simulation via LNA with restart
# Tmax = simulation end time
# dt = time step
# cvec = rate constant vector
# init = initial state
LV_LNA <- function(Tmax,dt,cvec,init){
  N <- Tmax/dt + 1
  x <- matrix(NA, nrow=N, ncol = 2)
  colnames(x) <- c("X_1", "X_2")
  xt <- init
  x[1,] <- xt
  for (i in 2:N) {
    t_eval <- c((i-2)*dt,(i-1)*dt)
    
    # Set up initial conditions, parameters, and time grid
    # eta_initial_state <- c(x1 = xt[1], x2 = xt[2])
    parameters <- c(c1 = cvec[1],c2 = cvec[2],c3 = cvec[3])
    
    # Call ode function with parameters
    eta_step <- ode(y = xt, times = t_eval, func = eta_ode_system, parms = parameters)
    eta_t <- eta_step[2,-1]
    
    
    # Set up initial conditions, parameters, and time grid
    V_initial_state <- matrix(c(0.0, 0.0, 0.0, 0.0), nrow = 2, byrow = TRUE)
    parameters2 <- c(parameters,x1 = xt[1],x2 = xt[2])
    
    # Call ode function with initial state as a matrix
    V_step <- ode(y = as.vector(V_initial_state), times = t_eval, func = V_ode_system, parms = parameters2)
    Vt <- matrix(V_step[2, -1], ncol = 2, byrow = TRUE)
    
    xt <- mvrnorm(1,eta_t, Vt)
    xt[which(xt < 0)] <- 0
    x[i,] <- xt
  }
  return(x)
}

# Lokta-Volterra simulation via Euler-Maruyama approximation of the CLE
# Tmax = simulation end time
# dt = time step
# cvec = rate constant vector
# init = initial state
LV_CLE <- function(Tmax,dt,cvec,init){
  N <- Tmax/dt + 1
  x <- matrix(NA, nrow=N, ncol = 2)
  colnames(x) <- c("X_1", "X_2")
  xt <- init
  x1 <- init[1]
  x2 <- init[2]
  x[1,] <- xt
  c1 <- cvec[1]
  c2 <- cvec[2]
  c3 <- cvec[3]
  S <- matrix(c(1,-1,0,
                0,1,-1),
              ncol = 3,byrow = TRUE)
  for (i in 2:N) {
    h1 <- c1 * x1
    h2 <- c2 * x1 * x2
    h3 <- c3 * x2
    h <- c(h1,h2,h3)
    alpha <- S%*%h
    beta <- matrix(c(c1*x1+c2*x1*x2, -c2*x1*x2,
                     -c2*x1*x2, c2*x1*x2+c3*x2),
                   nrow = 2, byrow = TRUE)
    xt <- mvrnorm(1,xt+alpha*dt, beta*dt)
    xt[which(xt < 0)] <- 0
    x[i,] <- xt
    x1 <- xt[1]
    x2 <- xt[2]
  }
  return(x)
}