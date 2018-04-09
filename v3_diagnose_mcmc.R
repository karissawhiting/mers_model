library(BayesianTools)

setup = createBayesianSetup(logPosterior_trunc,
                            lower = c(0,0,0,0))
out = runMCMC(setup, sampler = "DEzs", settings = NULL)


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


mcmc.epi <- mcmcMH(target = logPosterior_trunc,
                   init.theta = c(beta = .6, L =6, D0 = 7, D1 = 5, w = .01),
                   proposal.sd = c(.1, .5, .5, .5, .03),
                   n.iterations = 10000,
                   adapt.size.start = 100,
                   adapt.size.cooling=0.999,
                   limits = list(lower = c(beta = 0, L = 0, D0 = 0, D1 = 0, w = 0)))



out = runMCMC(setup, sampler = "DEzs", settings = NULL)