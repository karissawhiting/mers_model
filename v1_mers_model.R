library(tidyverse)
library('fitR')
library(reshape2)
library('coda')

# Data Clean -----------------------------

load("./data-raw/mers_times.RData")

mers_times2<- melt(mers_times[,2:4])

epi <- dcast(mers_times2,  value ~ variable, length) %>%
  rename("times" = value, "exp" = exp_date, "onset" = onset_date, "conf" = conf_date) %>%
  mutate(times = as.numeric(times))

times = data.frame(times = 1:56)
epi = full_join(times, epi, by = "times")

epi[is.na(epi)] <- 0

epi_melt <- melt(epi, id = "times")

ggplot(epi_melt, aes(x = times, y = value, color = variable)) + geom_line()

# Coding Model -----------------------------
SIR$name
SIR$state.names <- c("S", "E", "I", "C", "Exp", "Inc", "Con")
SIR$theta.names <- c("beta", "L", "D")

SIR$simulate <- function (theta, init.state, times) 
{
  SIR_ode <- function(time, state, parameters) {
    beta <- parameters[["beta"]]
    L <- parameters[["L"]]
    D <- parameters[["D"]]
    
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C <- state[["C"]]
    Exp <- state[["Exp"]]
    Inc <- state[["Inc"]]
    Con <- state[["Con"]]
    
    N <- S + E + I + C
    
    dS <- -beta * S * I/N
  
    S1 = max(S,0) #prevent from going below zero
    
    dE <- beta * S1 * I/N - (E/L)
    dI <- (E/L) - (I/D)
    dC <- I/D
    
    dExp <- beta * S1 * I/N
    dInc <- (E/L)
    dCon <- (I/D)
    
    return(list(c(dS, dE, dI, dC, dExp, dInc, dCon)))
  }
  
  # put incidence at 0 in init.state
  init.state["Exp"] <- 0
  init.state["Inc"] <- 0
  init.state["Con"] <- 0
  
  traj <- data.frame(ode(y = init.state, times = times, 
                               func = SIR_ode, parms = theta, method = "ode45"))
  
  # compute incidence of each time interval
  traj$Exp <- c(0, diff(traj$Exp))
  traj$Inc <- c(0, diff(traj$Inc))
  traj$Con <- c(0, diff(traj$Con))
  
  return(traj)
}



theta<- c(beta = 100,  L = .1 , D = 100) #8.07 #3.77
init.state <- c(S = 51413925, E = 0, I = 1, C = 0)
times <- 1:57
traj <- SIR$simulate(theta, init.state, times)

#theta<- c(beta = .01, L = .01 , D = .01) #8.07 #3.77
#traj <- SIR$simulate(theta, init.state, times)
#traj[traj < 0]

#plotTraj(traj)

# Priors -----------------------------
SIR$dprior <- function(theta, log = FALSE) {
  
  ## gamma prior on beta
  log.prior.beta <- dgamma(theta[["beta"]], shape =  1.5, rate = 2.0,  log = TRUE)
  
  # gamma prior on D
  log.prior.L <- dgamma(theta[["L"]], shape =  4.44, rate = .55,  log = TRUE)
  
  ## gamma prior on D
  log.prior.D <-dgamma(theta[["D"]], shape =  3.28, rate = .48,  log = TRUE)
  
  log.sum <- log.prior.beta + log.prior.L + log.prior.D
  
#  return(ifelse(log, log.sum, exp(log.sum)))
}

#SIR$dprior(theta)

# Log Likelihood of data  -----------------------------

# dPointsObs - calculates likelihood of an observed prevalance. e
# evaluates with respect to a Poisson dist centered around I member of model point.
# what is the probability of observing 18 cases when the prevalance is 31.
# assumes observation follows Poisson dist centered around current prevalance 


#* Exposure ----------------------
SIR$dPointObs_E <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a neg binomial process
  dnbinom(x = data.point[["exp"]], size = 10, #dispersion param
          prob = 10/(10+model.point[["Exp"]]),
          log = log)
  
}

SIR$dPointObs_E(epi[1,], traj[1,], theta, log = FALSE) 
SIR$dPointObs_E(epi[2,], traj[2,], theta, log = FALSE) 

 # calculates total log likelihood of a set of parameters and initial state.
dTrajObs_E <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
    if(traj[i,][["Exp"]] != 0) {
      data.point <- unlist(data[i, ])
      model.point <- unlist(traj[i, ])
      dens <- dens + fitmodel$dPointObs_E(data.point = data.point, 
                                          model.point = model.point, theta = theta, log = TRUE)
    }
    else {
      dens <- dens + 0
    }
  }
  return(ifelse(log, dens, exp(dens)))
}

ll_E<- dTrajObs_E(SIR, theta, init.state, epi, log = TRUE)


#* Infection  ----------------------
SIR$dPointObs_I <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a neg binomial process
  dnbinom(x = data.point[["onset"]], size = 10, #dispersion param
          prob = 10/(10+model.point[["Inc"]]),
          log = log)
  
}

SIR$dPointObs_I(epi[1,], traj[1,], theta, log = FALSE) 
SIR$dPointObs_I(epi[2,], traj[2,], theta, log = FALSE) 

# calculates total log likelihood of a set of parameters and initial state.
dTrajObs_I <- function (fitmodel, theta, init.state, data, log = TRUE) 
{
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
    if(traj[i,][["Inc"]] != 0) {
      data.point <- unlist(data[i, ])
      model.point <- unlist(traj[i, ])
      dens <- dens + fitmodel$dPointObs_I(data.point = data.point, 
                                          model.point = model.point, theta = theta, log = TRUE)
    }
    else {
      dens <- dens + 0
    }
  }
  return(ifelse(log, dens, exp(dens)))
}

LL_I <- dTrajObs_I(SIR, theta, init.state, epi, log = TRUE)

#* Confirmed ----------------------
SIR$dPointObs_C <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a neg binomial process
  dnbinom(x = data.point[["conf"]], size = 10, #dispersion param
          prob = 10/(10+model.point[["Con"]]),
          log = log)
}


# calculates total log likelihood of a set of parameters and initial state.
dTrajObs_C <- function (fitmodel, theta, init.state, data, log = TRUE) 
{
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
    if(traj[i,][["Con"]] != 0) {
      data.point <- unlist(data[i, ])
      model.point <- unlist(traj[i, ])
      dens <- dens + fitmodel$dPointObs_C(data.point = data.point, 
                                          model.point = model.point, theta = theta, log = TRUE)
    }
    else {
      dens <- dens + 0
    }
  }
  return(ifelse(log, dens, exp(dens)))
}

LL_C <- dTrajObs_C(SIR, theta, init.state, epi, log = TRUE)

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

my_dLogPosterior(SIR, theta, init.state, epi)


#wrapper function for log posterior
logPosterior_trunc <- function(theta) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = theta,
                          init.state = c(S = 51413925, E = 0, I = 1, C = 0),
                          data = epi))
  
}

my_dLogPosterior_epi3(theta)

mcmc.epi <- my_mcmcMH(target = logPosterior_trunc,
                       init.theta = c(beta = 10, L =2, D = 3),
                       proposal.sd = c(0.1, 0.1, .1),
                       n.iterations = 1000)




# Diagnostics  ----------------------

trace <- mcmc.epi$trace

