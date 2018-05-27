
library('devtools')
install_github('sbfnk/fitR')
library('fitR')

#Overview of Functions  ------------------
data(SIR)
SIR$name
SIR$state.names
SIR$theta.names

#dprior - calculates prior density
#dPointObs = calculates the log likelihood of a given data point with respect to the model
#rPointObs - generate observations from a model run

#simulate (part of SIR fitmode object) - runs the model
# provide parameters, initial states, and times
simuate <- function(theta,init.state,times) {
  SIR_ode <- function(time, state, parameters) {
    
    ## parameters
    beta <- parameters[["R0"]] / parameters[["D_inf"]]
    nu <- 1 / parameters[["D_inf"]]
    
    ## states
    S <- state[["S"]]
    I <- state[["I"]]
    R <- state[["R"]]
    
    N <- S + I + R
    
    dS <- -beta * S * I/N
    dI <- beta * S * I/N - nu * I
    dR <- nu * I
    
    return(list(c(dS, dI, dR)))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = SIR_ode,
                               parms = theta,
                               method = "ode45"))
  
  return(trajectory)
}

theta<- c(R0 = 3, D_inf = 2)
init.state <- c(S = 999, I = 1, R = 0)
times <- 1:50
traj <- SIR$simulate(theta, init.state, times)

head(traj)
plotTraj(traj)

#dprior- look at priors or parameters. calculates priors of parameters using uniform distributions
SIR$dprior(theta)

# dPointsObs - calculates likelihood of an observed prevalance. e
# evaluates with respect to a Poisson dist centered around I member of model point.
# what is the probability of observing 18 cases when the prevalance is 31.
# assumes observation follows Poisson dist centered around current prevalance 

SIR$dPointObs(data.point = c(obs= 18), model.point = c(I = 31), theta, log = TRUE)
dPointObs <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a Poisson process
  return(dpois(x = data.point[["obs"]],
               lambda = model.point[["I"]],
               log = log))
}


#Example ------------------

data(epi)
head(epi1)

plotTraj(epi1)

# dTrajObs - calculates total log likelihood of a set of parameters and initial state.
# simulates SIR model with given param and states, evaluates log likelihood at every point in given observed data
#returns sum of all log likelihoods

dTrajObs <- function (fitmodel, theta, init.state, data, log = FALSE) 
{
  times <- c(0, data$time)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i + 1, ])
    dens <- dens + fitmodel$dPointObs(data.point = data.point, 
                                      model.point = model.point, theta = theta, log = TRUE)
  }
  return(ifelse(log, dens, exp(dens)))
}

dTrajObs(SIR, theta, init.state, epi1, log = TRUE)

# rPointsObs generated randomly sampled data point given the model trajectory
#generates from Poisson process. Differs from simulate because encodes uncertainty involved
#in obervations 

SIR$rPointObs(model.point = c(I = 31), theta)

#rTrajObs simulated model using parameters and states, applys rPointObs at every point to generate 
#observations, returns simulated trjectory with added column for observations

rTrajObs <- function (fitmodel, theta, init.state, times) 
{
  traj <- fitmodel$simulate(theta, init.state, times)
  obs <- ddply(traj, "time", fitmodel$rPointObs, theta = theta)
  traj_obs <- join(traj, obs, by = "time")
  return(traj_obs)
}
obs.traj <- rTrajObs(SIR, theta, init.state, epi1$time)
head(obs.traj)

#Calculate Posterior ---------------------

# calculates value of posterior theta given data 
# This is a function that takes 4 arguments:
# - fitmodel, a fitmodel object that defines the model dynamics,
#   prior and likelihoods.
# - theta, a named vector of parameters
# - init.state,  a named vector of initial state
# - data, the data set we are fitting the model to
# It returns the posterior for the given model, parameters, initial
# state and data.
my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  
  # calculate the fitmodel prior for parameter vector theta using
  # fitmodel$dprior, and assign to variable log.prior
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  
  # calculate the log-likelihood of `theta`
  # and `init.state` with respect to the data using `dTrajObs`
  # and assign to a variable `log.likelihood`    
  log.likelihood <- dTrajObs(fitmodel, theta, init.state, data, log = TRUE)
  
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

my_dLogPosterior(SIR, theta, init.state, epi1)

plotFit(SIR, theta, init.state, epi1)
theta<- c(R0 = 10, D_inf = .5)
plotFit(SIR, theta, init.state, epi1)


# Explore Posterior for parameter est  ---------------------

#want to maximum log posterior

theta<- c(R0 = 10, D_inf = 2)
my_dLogPosterior(SIR, theta, init.state, epi1)

r0 = seq(0,6, by = .2)
post <- c()
for (i in 1:length(r0)){
  theta = c(R0 = r0[i], D_inf = 2)
  post[i] = my_dLogPosterior(SIR, theta, init.state, epi1)
}


plot(x = r0, y = post)

# try changing dpriors







