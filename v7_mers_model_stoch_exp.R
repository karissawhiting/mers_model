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

save(epi, file = "./data/epi.RData")

# Coding Model Expanded  -----------------------------
SIR$name <- c("SEIC model")
SIR$state.names <- c("S", "E", "I", "C", "Exp", "Inc", "Con")
SIR$theta.names <- c("beta", "L", "D1")

SIR$simulate <- function (theta, init.state, times) 
{
  SIR_transitions <- list(c(S = -1, E = 1, Exp = 1), c(E = -1, I = 1, Inc = 1), 
                            c(I = -1, C = 1, Con = 1))
  
  SIR_rateFunc <- function(state, theta, t) {
    beta <- theta[["beta"]]
    L <- theta[["L"]]
    D1<- theta[["D1"]]
    
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C <- state[["C"]]
    Exp <- state[["Exp"]]
    Inc <- state[["Inc"]]
    Con <- state[["Con"]]
    
    N <- S + E + I + C
    
    return(c(beta * S * I/N, E/L, I/D1))
  }
  # put incidence at 0 in init.state
  init.state["Exp"] <- 1
  init.state["Inc"] <- 1
  init.state["Con"] <- 1
  
  traj <- simulateModelStochastic(theta, init.state, times, 
                                  SIR_transitions, SIR_rateFunc)
  traj$Inc <- c(0, diff(traj$Inc))
  return(traj)
}

theta<- c(beta = .99,  L = 8.19 , D1 = 4.05, w = .09)
init.state <- c(S = 10000, E = 0, I = 1, C = 0) 
times <- 1:100

sim<- SIR$simulate(theta, init.state, times)












$rPointObs
function (model.point, theta) 
{
  obs.point <- rpois(n = 1, lambda = theta[["rho"]] * model.point[["Inc"]])
  return(c(obs = obs.point))
}

$dprior
function (theta, log = FALSE) 
{
  log.prior.R0 <- dunif(theta[["R0"]], min = 1, max = 50, log = TRUE)
  log.prior.latent.period <- dunif(theta[["D_lat"]], min = 0, 
                                   max = 10, log = TRUE)
  log.prior.infectious.period <- dunif(theta[["D_inf"]], min = 0, 
                                       max = 15, log = TRUE)
  log.prior.temporary.immune.period <- dunif(theta[["D_imm"]], 
                                             min = 0, max = 50, log = TRUE)
  log.prior.probability.long.term.immunity <- dunif(theta[["alpha"]], 
                                                    min = 0, max = 1, log = TRUE)
  log.prior.reporting.rate <- dunif(theta[["rho"]], min = 0, 
                                    max = 1, log = TRUE)
  log.sum = log.prior.R0 + log.prior.latent.period + log.prior.infectious.period + 
    log.prior.temporary.immune.period + log.prior.probability.long.term.immunity + 
    log.prior.reporting.rate
  return(ifelse(log, log.sum, exp(log.sum)))
}

$dPointObs
function (data.point, model.point, theta, log = FALSE) 
{
  return(dpois(x = data.point[["obs"]], lambda = theta[["rho"]] * 
                 model.point[["Inc"]], log = log))
}

attr(,"class")
[1] "fitmodel"