library(tidyverse)
library('fitR')
library(reshape2)
library('coda')
library(lattice)
library(MCMCvis)
library(wesanderson)
library(deSolve)

#Data Clean -----------------------------

load("./data-raw/mers_times.RData")

mers_times2<- melt(mers_times[,2:4])

epi <- dcast(mers_times2,  value ~ variable, length) %>%
  rename("times" = value, "exp" = exp_date, "onset" = onset_date, "conf" = conf_date) %>%
  mutate(times = as.numeric(times))

times = data.frame(times = 1:56)
epi = full_join(times, epi, by = "times")

epi[is.na(epi)] <- 0
epi<- epi[1:56,]

rm(mers_times, mers_times2, times)

#save(epi, file = "./data/epi.RData")

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

SIR$theta.names <- c("beta", "L", "D0", "D1", "w", "n")

# Priors -----------------------------
SIR$dprior <- function(theta, log) {
  
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
  
  ## gamma prior on n
  log.prior.n <-dgamma(theta[["n"]], shape = 3.125, rate = .3125,  log = TRUE)
  
  #sum
  log.sum <- log.prior.beta + log.prior.L + log.prior.D0 + log.prior.D1 + log.prior.w + log.prior.n
  
  return(log.sum)
  #  return(ifelse(log, log.sum, exp(log.sum)))
}

#SIR$dprior(c(beta = .99,  L = 8.19 , D0 = 9.26, D1 = 4.05, w = .09, n = 10))
# Log Likelihood of data  -----------------------------

#* Exposed ----------------------

#size = m/(n-1), prob = 1/n
dTrajObs_E <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 5:49) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    if(model.point[["Exp"]] > 0) {
      dens <- sum(dens,  
        (dnbinom(x = data.point[["exp"]],
              mu = (model.point[["Exp"]]),
              size = ((model.point[["Exp"]])/(theta[["n"]]-1)),
              log = log)), na.rm = TRUE)
    }
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}
#dTrajObs_E(SIR, theta, init.state, epi, log = TRUE)

#* Infection  ----------------------
dTrajObs_I <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:53) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    if(model.point[["Inc"]] > 0) {
    dens <- sum(dens,
      (dnbinom(x = data.point[["onset"]],
              mu = model.point[["Inc"]],
              size = ((model.point[["Inc"]])/(theta[["n"]]-1)),
              log = log)), na.rm = TRUE)
    }
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}
#dTrajObs_I(SIR, theta, init.state, epi, log = TRUE)

#theta<- c(beta = .99,  L = 8.19 , D0 = 9.26, D1 = 4.05, w = .09, n = 2)
#dTrajObs_I(SIR, theta = theta, init.state = init.state, data = epi, log = TRUE)



#* Confirmed ----------------------
dTrajObs_C <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 10:55) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    if(model.point[["Con"]] > 0) {
    dens <- sum(dens,  
      (dnbinom(x = data.point[["conf"]], #dispersion param
              mu = (model.point[["Con"]]),
              size = ((model.point[["Con"]])/(theta[["n"]]-1)),
              log = log)), na.rm = TRUE)
    }
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}

#dTrajObs_C(SIR, theta, init.state, epi, log = TRUE)
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
                          init.state = c(S = 10000, E = 0, I = 1, C = 0),
                          data = epi))
  
}

#logPosterior_trunc(c(beta = .99,  L = 8.19 , D0 = 9.26, D1 = 4.05, w = .09, n = 10))

# Run MCMC  -----------------------------
# look at beta L and D simultaneous. 3-D plot. 
# try different intial  vlaues

# parameters
adapt.size.start <- 1000
adapt.size.cooling <- 0.999
adapt.shape.start <- 3000
adapt.shape.stop <- 17000
n.iterations <- 50000
proposal.sd <- c(beta = .03, L =.5, D0 = .5, D1 = .2, w = .03, n = .05)


mcmc.epi <- mcmcMH(target = logPosterior_trunc,
                    init.theta = c(beta = .6, L =6, D0 = 11, D1 = 3, w = .06, n = 6),
                    proposal.sd = proposal.sd,
                    n.iterations = n.iterations,
                    adapt.size.start = adapt.size.start,
                    adapt.size.cooling = adapt.size.cooling,
                   adapt.shape.start = adapt.shape.start, 
                   adapt.shape.stop = adapt.shape.stop,
                   limits = list(lower = c(beta = 0, L = 0, D0 = 0, D1 = 0, w = 0, n = 1.000000001)))

save(mcmc.epi,  file= "./data/mcmc.epi.RData")
#mcmc.epi$covmat.empirical


trace <- mcmc.epi$trace
mcmc.trace <- mcmc(trace)
summary(mcmc.trace)
save(mcmc.trace,  file= "./data/trace.RData")



mcmc.epi$covmat.empirical

acceptanceRate <- 1 - rejectionRate(mcmc.trace)
#effectiveSize(mcmc.trace)
#plot(mcmc.trace, denseplot = TRUE)

mcmc.trace.burned <- burnAndThin(mcmc.trace, burn = 5000, thin = 3)
acceptanceRate <- 1 - rejectionRate(mcmc.trace.burned)
plot(mcmc.trace.burned)

autocorr.plot(mcmc.trace.burned)
plotESSBurn(mcmc.trace)

#xyplot(x = mcmc.trace.burned)
#densplot(mcmc.trace['beta',])
#MCMCtrace(mcmc.trace)
#MCMCchains(mcmc.trace)

ggplot(data,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.4)

b<- data.frame(mcmcdf[,'beta'][10000:50000])
L <- data.frame(mcmcdf[,'L'][10000:50000])
D0 <- data.frame(mcmcdf[,'D0'][10000:50000])
D1 <- data.frame(mcmcdf[,'D1'][10000:50000])
w <- data.frame(mcmcdf[,'w'][10000:50000])
n <- data.frame(mcmcdf[,'n'][10000:50000])


ggplot(b,aes(x=b)) + 
  geom_density(fill = "pink") + 
  geom_density(aes(L), fill = "blue", alpha = .5) + 
  geom_density(aes(D0), fill = "red", alpha = .5) + 
  geom_density(aes(D1), fill = "green", alpha = .5) +
  geom_density(aes(w), fill = "orange", alpha = .5) + 
  geom_density(aes(n), fill = "yellow", alpha = .5) 

ggplot(w,aes(x=w)) + 
  geom_density()
