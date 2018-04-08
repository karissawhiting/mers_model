#################################
#try skipping when trajectory is zero
dTrajObs_I <- function (fitmodel, theta, init.state, data = epi, log = TRUE) {
  times <- c(0, data$times)
  traj <<- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
 #   print(dens)
    data.point <- unlist(data[i, ])
#    print(data.point[["onset"]])
    model.point <- unlist(traj[i, ])
#    print(model.point[["Inc"]])
    if (model.point[["Inc"]] < 0){
      next
    } else {
      dens <- dens + 
        dpois(x = data.point[["onset"]], #dispersion param
            lambda = model.point[["Inc"]] + .001,
            log = log)
    
  } 
  # return(ifelse(log, dens, exp(dens)))
  }
  return(dens)
}

dTrajObs_I(fitmodel = SIR, theta = theta, init.state = init.state, data = epi, log = TRUE)

#################################
# try with adding small value to trajectory
dTrajObs_I2 <- function (fitmodel, theta, init.state, data = epi, log = TRUE) {
  times <- c(0, data$times)
  traj <<- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
#    print(data.point[["onset"]])
    model.point <- unlist(traj[i, ])
 #   print(model.point[["Inc"]])
    dens <- dens + 
      dpois(x = data.point[["onset"]], #dispersion param
            lambda = model.point[["Inc"]] + .001,
            log = log)
      
    # return(ifelse(log, dens, exp(dens)))
  }
  return(dens)
}

theta <- c(beta = 3.31048030200639,  L = .8066304953947 , D = 0.345418814011432)

dTrajObs_C(fitmodel = SIR, theta = theta, init.state = init.state, data = epi, log = TRUE)



#################################
# debug full log posterior

theta<- c(beta = 3.31048030200639,  L = .8066304953947 , D = 0.345418814011432) #8.07 #3.77

theta<- c(beta = 3.31048030200639,  L = .8066304953947 , D = 0.3) #

logPosterior_trunc <- function(theta) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = theta,
                          init.state = c(S = 51413925, E = 0, I = 1, C = 0),
                          data = epi))
  
}



logPosterior_trunc(theta)



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


theta<- c(beta = 0.327696,  L = 8.16870 , D = 6.967140) #8.07 #3.77
init.state <- c(S = 51413925, E = 0, I = 1, C = 0)
times <- 1:57
#traj <- SIR$simulate(theta, init.state, times)

logPosterior_trunc(theta)


#######

#Try a grid of theta
my_dLogPosterior

beta <- seq(.4, 3, by = .1)
L <- seq(6, 9, by = .1)
D<- seq(.5, 6, by = .1)

cbin



##########################################
#found example of poisson error 

# this causes errors
mcmcMH(target = logPosterior_trunc, 
                       init.theta = c(beta = .1, L =2, D = 2),
                        proposal.sd = c(.01, .1, .1),
                        n.iterations = 1000,
                        limits = list(lower = c(beta = 0, L = 0, D = 0)))



theta<- c(beta = .1, L =2, D = 2)
logPosterior_trunc(theta)

my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  print(log.prior)
  # calculate the log-likelihood of `theta`
  log.likelihood <- dTrajObs_E(fitmodel, theta, init.state, data, log = TRUE) +
    dTrajObs_I(fitmodel, theta, init.state, data, log = TRUE) +
    dTrajObs_C(fitmodel, theta, init.state, data, log = TRUE)
  
  print(log.likelihood)
  
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  
  print(log.posterior)
  return(log.posterior)
  
}


my_dLogPosterior(SIR, theta, init.state, epi)



##############



dTrajObs_E <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <<- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 5:56) {
    print(i)
    data.point <- unlist(data[i, ])
    print(data.point[["exp"]])
    model.point <- unlist(traj[i+1, ])
    print(model.point[["Exp"]])
    dens <- dens + 
      dpois(x = data.point[["exp"]], #dispersion param
            lambda = (model.point[["Exp"]] + .0000000001),
            log = log)
    
  }
  print(dens)
  # return(ifelse(log, dens, exp(dens)))
}


dens <- 0
for (i in 5:56) {
  print(i)
  print(epi[["exp"]][[i]])
  print(traj[["Exp"]][[i+1]] + .000001)
  dens[i] <- dpois(x = epi[["exp"]][[i]],
          lambda = (traj[["Exp"]][[i+1]] + .0000000001),
          log = log)

}

for (i in 5:56) {
  l <- (traj[["Exp"]][[i+1]] + .0001)
  print(i)
  print(epi[["exp"]][[i]])
  print(l)
den <- dpois(x = epi[["exp"]][[i]], #dispersion param
      lambda = l,
      log = log)
print(den)
}


dpois(x = 0, #dispersion param
      lambda = -7.589595e-08,
      log = log)



dens <- 0
for (i in 5:56) {
  dens <- dens + 
    dpois(x = data.point[["exp"]], #dispersion param
          lambda = (model.point[["Exp"]] + .0000000001),
          log = log)
  
}





dTrajObs_E(SIR, theta = c(beta = .1, L =2, D = 2), 
           init.state = c(S = 51413925, E = 0, I = 1, C = 0), epi, log = TRUE)

dpois(x = 1, #dispersion param
      lambda = 0.04517858,
      log = log)

#* Infection  ----------------------
dTrajObs_I <- function (fitmodel = SIR, theta, init.state = init.state, data = epi, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    dens <- dens + 
      dpois(x = data.point[["onset"]], #dispersion param
            lambda = model.point[["Inc"]] + .0000000001,
            log = log)
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}


#* Confirmed ----------------------
dTrajObs_C <- function (fitmodel = SIR, theta, init.state = init.state, data = epi, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 10:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    dens <- dens + 
      dpois(x = data.point[["conf"]], #dispersion param
            lambda = model.point[["Con"]] + .0000000001,
            log = log)
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}

