library(tidyverse)
library('fitR')
library(reshape2)
library('coda')


#Data Clean -----------------------------

load("./data-raw/mers_times.RData")

mers_times2<- melt(mers_times[,2:4])

epi <- dcast(mers_times2,  value ~ variable, length) %>%
  rename("times" = value, "exp" = exp_date, "onset" = onset_date, "conf" = conf_date) %>%
  mutate(times = as.numeric(times))

times = data.frame(times = 1:56)
epi = full_join(times, epi, by = "times")

epi[is.na(epi)] <- 0

rm(mers_times, mers_times2, times)



# Coding Model Expanded  -----------------------------
SIR$name <- c("SEIC model with w, D0, D1 param")
SIR$state.names <- c("S", "E", "I", "C", "Exp", "Inc", "Con")
SIR$theta.names <- c("beta", "L", "D0", "D1", "w")


SIR$simulate <- function (theta, init.state, times) 
{
  SIR_ode <- function(time, state, parameters) {
    
    beta <- parameters[["beta"]]
    L <- parameters[["L"]]
    D0 <- parameters[["D0"]]
    D1<- parameters[["D1"]]
    w <- parameters[["w"]]
    
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C <- state[["C"]]
    Exp <- state[["Exp"]]
    Inc <- state[["Inc"]]
    Con <- state[["Con"]]
    
    N <- S + E + I + C
    
    if (time < 18) {
      dS <- -1*beta * S * I/N
      S1 = max(S,0) #prevent from going below zero
      dE <- 1*beta * (S1 * I)/N - (E/L)
      
      dI <- (E/L) - (I/D0)
      dC <- I/D0
      
      dExp <- beta * S1 * I/N
      dInc <- (E/L)
      dCon <- (I/D0)
      
    }
    else {
      dS <- -w*beta * S * I/N
      S1 = max(S,0) #prevent from going below zero
      dE <- (w*beta * (S1 * I)/N) - (E/L)
      
      dI <- (E/L) - (I/D1)
      dC <- I/D1
      
      dExp <- w*beta * S1 * I/N
      dInc <- (E/L)
      dCon <- (I/D1)
    }
    
    return(list(c(dS, dE, dI, dC, dExp, dInc, dCon)))
  }
  
  # put incidence at 0 in init.state
  init.state["Exp"] <- 0
  init.state["Inc"] <- 0
  init.state["Con"] <- 0
  
  traj <- data.frame(ode(y = init.state, times = times, method = "ode45",
                         func = SIR_ode, parms = theta)) #events = list(data = eventdat))) #ode45
  
  # compute incidence of each time interval
  traj$Exp <- c(0, diff(traj$Exp))
  traj$Inc <- c(0, diff(traj$Inc))
  traj$Con <- c(0, diff(traj$Con))
  
  return(traj)
}


#theta<- c(beta = .99,  L = 8.19 , D0 = 9.26, D1 = 4.05, w = .09) #8.07 #3.77
#init.state <- c(S = 10000, E = 0, I = 1, C = 0) #S = 51413925
#times <- 1:57
#traj <- SIR$simulate(theta, init.state, times)


# Priors -----------------------------
SIR$dprior <- function(theta, log = FALSE) {
  
  ## gamma prior on beta
  log.prior.beta <- dgamma(theta[["beta"]], shape =  1.5, rate = 2.0,  log = TRUE)
  
  # gamma prior on L
  log.prior.L <- dgamma(theta[["L"]], shape =  4.44, rate = .55,  log = TRUE)
  
  ## gamma prior on D0
  log.prior.D0 <-dgamma(theta[["D0"]], shape =  3.28, rate = .48,  log = TRUE)
  
  ## gamma prior on D1
  log.prior.D1 <-dgamma(theta[["D1"]], shape =  3.28, rate = .48,  log = TRUE)
  
  ## gamma prior on w
  log.prior.w <-dgamma(theta[["w"]], shape = 2, rate = 2,  log = TRUE)
  
  #sum
  log.sum <- log.prior.beta + log.prior.L + log.prior.D0 + log.prior.D1 + log.prior.w
  
  return(log.sum)
  #  return(ifelse(log, log.sum, exp(log.sum)))
}

# Log Likelihood of data  -----------------------------

#* Exposed ----------------------


dTrajObs_E <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 5:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    dens <- dens + 
      dpois(x = data.point[["exp"]], #dispersion param
            lambda = model.point[["Exp"]] + .0001,
            log = log)
    
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}

#* Infection  ----------------------
dTrajObs_I <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    dens <- dens + 
      dpois(x = data.point[["onset"]], #dispersion param
            lambda = model.point[["Inc"]] + .0001,
            log = log)
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}


#* Confirmed ----------------------
dTrajObs_C <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 10:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    dens <- dens + 
      dpois(x = data.point[["conf"]], #dispersion param
            lambda = model.point[["Con"]] + .0001,
            log = log)
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}

# Calculating Posterior -----------------------------

# dlogposterior function
my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  # calculate the log-likelihood of `theta`
  log.likelihood <- dTrajObs_E(fitmodel, theta, init.state, data, log = TRUE) +
    dTrajObs_I(fitmodel, theta, init.state, data, log = TRUE) +
    dTrajObs_C(fitmodel, theta, init.state, data, log = TRUE)
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  return(log.posterior)
  
}


#wrapper function for log posterior
logPosterior_trunc <- function(theta) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = theta,
                          init.state = c(S = 51413925, E = 0, I = 1, C = 0),
                          data = epi))
  
}

logPosterior_trunc(c(beta = .99,  L = 8.19 , D0 = 9.26, D1 = 4.05, w = .09))

# Run MCMC  -----------------------------

# This runs pretty well except there is a flat part and also beta isn't great. might want to tweak beta sd
mcmc.epi2 <- mcmcMH(target = logPosterior_trunc,
                   init.theta = c(beta = .6, L =6, D0 = 7, D1 = 5, w = .06),
                   proposal.sd = c(.03, .8, .8, .8, .02),
                   n.iterations = 15000,
                   adapt.size.start = 1000,
                   adapt.size.cooling=0.999,
                   limits = list(lower = c(beta = 0, L = 0, D0 = 0, D1 = 0, w = 0)))


trace <- mcmc.epi2$trace
mcmc.trace <- mcmc(trace)
summary(mcmc.trace)

acceptanceRate <- 1 - rejectionRate(mcmc.trace)
effectiveSize(mcmc.trace)
plot(mcmc.trace)

mcmc.trace.burned <- burnAndThin(mcmc.trace, burn = 1000, thin = 5)
plot(mcmc.trace.burned)

autocorr.plot(mcmc.trace.burned)
plotESSBurn(mcmc.trace)


correlationPlot(data.frame(mcmc.trace))

