# Functions used in Chapter 2: Stochastic differential equations

# 2.1 Ordinary differential equations
# 2.1.1 Application: SIR epidemic

# SIR numerical integration function
# init = initial state
# beta = beta parameter
# gamma = gamma parameter
# Tmax = simulation end time
# dt = time step
ep_ode <- function(init, beta, gamma, Tmax, dt){
  n <- Tmax/dt + 1
  sir <- matrix(NA, nrow = n, ncol = 3)
  sir[1,] <- init
  pop <- sum(init)
  for(j in 2:n){
    s <- sir[j-1,1]
    i <- sir[j-1,2]
    r <- sir[j-1,3]
    sir[j,1] <- s + (-beta*i*s/pop)*dt
    sir[j,2] <- i + (beta*i*s/pop - gamma*i)*dt
    sir[j,3] <- r + (gamma*i)*dt
  }
  return(sir)
}


# 2.1.2 Brownian motions

# Standard Brownian motion skeleton path simulation function
# Tmax = simulation end time
# dt = time step
BM <- function(Tmax=100,dt=1){
  n <- Tmax+1
  wt <- rep(0,n)
  for (i in 2:n) {
    wt[i] <- rnorm(1,wt[i-1],sqrt(dt))
  }
  return(wt)
}

# Geonetric Brownian motion skeleton path simulation function
# Tmax = simulation end time
# dt = time step
# a = drift
# b = diffusion coefficient
GenBM <- function(Tmax=100,dt=1,a=0,b=1){
  n <- Tmax+1
  xt <- rep(0,n)
  for (i in 2:n) {
    xt[i] <- rnorm(1,xt[i-1]+a*dt,b*sqrt(dt))
  }
  return(xt)
}

# Geometric Brownian motion skeleton path simulation function
# Tmax = simulation end time
# dt = time step
# m = mu (i.e mean rate of return in Black-Scholes model)
# s = sigma (i.e. volatility in Black-Scholes model)
GeoBM <- function(Tmax=100,dt=1,m=0.05,s=0.1){
  n <- Tmax+1
  wt <- rep(0,n)
  st <- rep(1,n)
  for (i in 2:n) {
    delta_t <- t[i]-t[i-1]
    wt[i] <- rnorm(1,wt[i-1],sqrt(delta_t))
    st[i] <- exp((m-s^2/2)*t[i]+s*wt[i])
  }
  return(st)
}



# 2.2 Stochastic calculus
# 2.2.5  Solving SDEs numerically

# Geometric Brownian motion skeleton path simulation function with analytic transition density
# mu = mean rate of return
# sigma = volatility
# dt = time step
# Tmax = simulation end time
# init = initial state
An_GeoBM <- function(mu = 0.05, sigma = 0.1, dt = 1, Tmax = 100, init = 1){
  n <- Tmax/dt + 1
  st <- rep(NA, n)
  st[1] <- init
  for(i in 2:n){
    st[i] <- st[i-1]*exp((mu-sigma^2/2)*dt + sigma*rnorm(1,0,sqrt(dt)))
  }
  return(st)
} 

# Geometric Brownian motion skeleton path simulation function with Euler-Maruyama transition density
# mu = mean rate of return
# sigma = volatility
# dt = time step
# Tmax = simulation end time
# init = initial state
EM_GeoBM <- function(mu = 0.05, sigma = 0.1, dt = 1, Tmax = 100, init = 1){
  n <- Tmax/dt + 1
  st <- rep(NA, n)
  st[1] <- init
  for(i in 2:n){
    st_dt <- st[i-1]
    st[i] <- st_dt + mu*st_dt*dt + st_dt*sigma*rnorm(1,0,sqrt(dt))
  }
  return(st)
}


