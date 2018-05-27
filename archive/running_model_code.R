library(ggplot2)

theta<- c(beta = 0.75,  L = 8.16870 , D = 6.83) #8.07 #3.77
init.state <- c(S = 51413925, E = 0, I = 1, C = 0)
times <- 1:57
traj <- SIR$simulate(theta, init.state, times)


ggplot(data = traj, aes(x = time)) + 
  geom_line(aes(y = S), color = "blue") + 
  geom_line(aes(y = E), color = "purple") +
  geom_line(aes(y = I), color = "green") +
  geom_line(aes(y = C), color = "orange") +
  geom_line(aes(y = Exp), color = "green") +
  geom_line(aes(y = Inc), color = "purple" ) +
  geom_line(aes(y = Con), color = "orange") +
  labs(y = "") #





SIR$dprior(theta)

SIR$dPointObs_E(epi[1,], traj[1,], theta, log = FALSE) 
SIR$dPointObs_E(epi[2,], traj[2,], theta, log = FALSE) 
dTrajObs_E(SIR, theta, init.state, epi, log = TRUE)

SIR$dPointObs_I(epi[2,], traj[2,], theta, log = FALSE) 
dTrajObs_E(SIR, theta, init.state, epi, log = TRUE)

SIR$dPointObs_C(epi[2,], traj[2,], theta, log = FALSE) 
dTrajObs_C(SIR, theta, init.state, epi, log = TRUE)

my_dLogPosterior(SIR, theta, init.state, data)

logPosterior_trunc(theta)


# Troubleshoot


SIR$dPointObs_E <- function(data.point, model.point, theta, log = FALSE){
    
    ## the prevalence is observed through a Poisson process
    return(dpois(x = data.point[["exp"]],
                 lambda = model.point[["Exp"]],
                 log = log))
  }
}




for (i in 1:57){
  x <- SIR$dPointObs_E(epi[i,], traj[i,], theta, log = TRUE)
  print(i)
  print(x)
  
}

SIR$dPointObs_E(epi[27,], traj[27,], theta, log = TRUE)

dnbinom(x = 9, size = 10, #dispersion param
        prob = 10/(10 + 0),
        log = log)

#when model point is zero there is problem. prob = 1, and throws error on dbinom 
dnbinom(x = 9, size = 10, #dispersion param
        prob = traj[27,][["Exp"]]/(traj[27,][["Exp"]] + 10),
        log = log)

dnbinom(x = 9, size = 10, #dispersion param
        mu = 1,
        log = log)




# Debug ---------------

SIR$dPointObs_E <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a neg binomial process
  return(return(dpois(x = data.point[["exp"]],
                      lambda = model.point[["Exp"]],
                      log = log)))
  
}

function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a Poisson process
  return(dpois(x = data.point[["obs"]],
               lambda = model.point[["I"]],
               log = log))
}

SIR$dPointObs_E(epi[2,], traj[2,], theta, log = FALSE) 

# calculates total log likelihood of a set of parameters and initial state.
dTrajObs_E <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    
    dens <- dens + fitmodel$dPointObs_E(data.point = data.point, 
                                        model.point = model.point, theta = theta, log = TRUE)
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}

dTrajObs_E (SIR, theta, init.state, epi, log = TRUE)



my_dLogPosterior(SIR, theta, init.state, epi)


dpois(x = .4, #dispersion param
      lambda = 0,
      log = log)

dpois(x = 0, #dispersion param
      lambda = 0, log = log)






my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  # calculate the log-likelihood of `theta`
  log.likelihood <- dTrajObs_E(fitmodel, theta, init.state, data, log = TRUE) 
  dTrajObs_I(fitmodel, theta, init.state, data, log = TRUE) +
    dTrajObs_C(fitmodel, theta, init.state, data, log = TRUE)
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  return(log.posterior)
  
}


theta=  c(beta = 3, L = 2.67397118966884, D = 0.460165930030681)
log.likelihood <- dTrajObs_E(SIR, theta, init.state, epi, log = TRUE) 
dTrajObs_I(SIR, theta, init.state, epi, log = TRUE) 
dTrajObs_C(SIR, theta, init.state, epi, log = TRUE)
# calulate the log-posterior using the log-prior and log-likelihood
log.posterior <- log.prior + log.likelihood


traj_test<- SIR$simulate(theta, init.state, times)


dTrajObs_I <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
    print(dens)
    data.point <- unlist(data[i, ])
    print(data.point)
    model.point <- unlist(traj[i, ])
    print(model.point)
    dens <- dens + fitmodel$dPointObs_I(data.point = data.point, 
                                        model.point = model.point, theta = theta, log = TRUE)
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}

SIR$dPointObs_I <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a neg binomial process
  dpois(x = data.point[["onset"]], #dispersion param
        lambda = model.point[["Inc"]],
        log = log)
  
}

SIR$dPointObs_I(epi[1,], traj[1,], theta, log = FALSE)

epi[1,][["onset"]]
traj[1,][["Inc"]]

dpois(x = 0, #dispersion param
      lambda = 0,
      log = FALSE)


dTrajObs_I(SIR, theta, init.state, epi, log = TRUE)


############
#neg binomial

theta<- c(beta = .99,  L = 8.19 , D0 = 9.26, D1 = 4.05, w = .09) #8.07 #3.77
init.state <- c(S = 10000, E = 0, I = 1, C = 0) #S = 51413925
times <- 1:57
traj <- SIR$simulate(theta, init.state, times)


