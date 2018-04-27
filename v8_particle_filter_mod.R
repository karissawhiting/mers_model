#Second attempt

library(tidyverse)
library('fitR')
library(reshape2)
library('coda')
library(lattice)
library(MCMCvis)
library(wesanderson)
library(deSolve)
library(adaptivetau)

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
SIR$theta.names <- c("beta", "L", "D1", "w")

SIR$simulate <- function (theta, init.state, times) 
{
  SIR_transitions <- list(c(S = -1, E = 1, Exp = 1), c(E = -1, I = 1, Inc = 1), 
                          c(I = -1, C = 1, Con = 1))
  
  SIR_rateFunc <- function(state, theta, t) {
    beta <- theta[["beta"]]
    L <- theta[["L"]]
 #   D0<- theta[["D0"]]
    D1<- theta[["D1"]]
    w <- theta[["w"]]
    
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C <- state[["C"]]
    Exp <- state[["Exp"]]
    Inc <- state[["Inc"]]
    Con <- state[["Con"]]
    
    N <- S + E + I + C
    
    return(c(w*beta * S * I/N, E/L, I/D1))
}
  # put incidence at 0 in init.state
  init.state["Exp"] <- 1
  init.state["Inc"] <- 1
  init.state["Con"] <- 1
  
  traj <- simulateModelStochastic(theta, init.state, times, 
                                  SIR_transitions, SIR_rateFunc)
  traj$Exp <- c(0, diff(traj$Exp))
  traj$Inc <- c(0, diff(traj$Inc))
  traj$Con <- c(0, diff(traj$Con))
  return(traj)
}

theta<- c(beta = .99,  L = 8.19 , D0 = 9.26, D1 = 4.05, w = .09)
init.state <- c(S = 10000, E = 0, I = 1, C = 0) 
times <- 1:100

SIR$simulate(theta, init.state, times)

my_particleFilter(SIR, theta, init.state, epi, n.particles = 10)




SIR$rPointObs <- function (model.point, theta) 
{
  obs.point <- rpois(n = 1, lambda = model.point[["Inc"]])
  return(c(obs = obs.point))
}


SIR$dprior <- function (theta, log) {
  log.prior.beta <- dgamma(theta[["beta"]], shape = 1.5, rate = 2, 
                           log = TRUE)
  log.prior.L <- dgamma(theta[["L"]], shape = 4.44, rate = 0.55, 
                        log = TRUE)
  log.prior.D0 <- dgamma(theta[["D0"]], shape = 3.28, rate = 0.48, 
                         log = TRUE)
  log.prior.D1 <- dgamma(theta[["D1"]], shape = 3.28, rate = 0.48, 
                         log = TRUE)
  log.prior.w <- dgamma(theta[["w"]], shape = 2, rate = 2, 
                        log = TRUE)
  log.prior.n <- dgamma(theta[["n"]], shape = 3.125, rate = 0.3125, 
                        log = TRUE)

  log.sum <- log.prior.beta + log.prior.L + log.prior.D0 + 
    log.prior.D1 + log.prior.w + log.prior.n 
  return(log.sum)
}

SIR$dPointObs <- function (data.point, model.point, theta, log = FALSE) 
{
  return(dpois(x = data.point[["onset"]], lambda = model.point[["Inc"]], 
               log = log))
}



