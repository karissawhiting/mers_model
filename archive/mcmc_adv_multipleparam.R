library('coda')

# More than one parameter -------------------------

#metropolis hasting but with multivariate normal as proposal function
plotTraj(data = epi3)

# dlogposterior function
my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  # calculate the log-likelihood of `theta`
  log.likelihood <- dTrajObs(fitmodel, theta, init.state, data, log = TRUE)
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  return(log.posterior)
  
}

#wrapper function for dlog_posterio function for 2 paramters
my_dLogPosterior_epi3 <- function(theta) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = theta,
                          init.state = c(S = 999, I = 1, R = 0),
                          data = epi3))
  
}

mcmc.epi3 <- my_mcmcMH(target = my_dLogPosterior_epi3,
                       init.theta = c(R0 = 1, D_inf = 2),
                       proposal.sd = c(0.01, 0.1),
                       n.iterations = 1000)


plot(trace, type = 'l')


# Diagnostics -----------------------------

head(trace)

# convert trace object for use in code package - makes it into a chain object
mcmc.trace <- mcmc(trace)
summary(mcmc.trace) #show mean, sd empircal of each param. naive se (of mean) is adjuted for sample size
# time series SE corrects for autocorrelation

#acceptance rate
acceptanceRate <- 1 - rejectionRate(mcmc.trace)
acceptanceRate

#effective sample size (for autocorrel of samp) - estimate for number of indepdennt samples taking into account autocorelations that were generated
effectiveSize(mcmc.trace)

#* Mixing -----------------------------

#can see burnout and also is there are too many
#sticks (steps) or it moves too much in one direction
# can look at range of param
plot(mcmc.trace)

#burn and thin
mcmc.trace.burned <- burnAndThin(mcmc.trace, burn = 1000)
plot(mcmc.trace.burned)


#autocorrelation

#high dregree of correlaiton between samples means slow mixing
# should get smaller as lag k increases
autocorr.plot(mcmc.trace.burned)

#if autocorrelation gets higher we can thin
mcmc.trace.burned.thinned <- burnAndThin(mcmc.trace.burned, thin = 5)
autocorr.plot(mcmc.trace.burned.thinned)


#truncate proposal function for limited param supports
my_dLogPosterior_epi3 <- function(theta) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = theta,
                          init.state = c(S = 999, I = 1, R = 0),
                          data = epi3))
  
}
trace <- mcmcMH(target = my_dLogPosterior_epi3,
                init.theta = c(R0 = 1, D_inf = 2),
                proposal.sd = c(0.1, 0.01),
                n.iterations = 1000,
                limits = list(lower = c(R0 = 0, D_inf = 0)))

#* Adaptive Proposal Functions  -----------------------------

#The best proposal distribution is the one that best matches the target distribution.

trace$covmat.empirical

#can also adapt start size, shape and cooling

trace <- mcmcMH(target = my_dLogPosterior_epi3,
                init.theta = c(R0 = 1, D_inf = 2),
                proposal.sd = c(1, 0.5),
                n.iterations = 5000,
                adapt.size.start = 100,
                adapt.shape.start = 500,
                adapt.size.cooling=0.999,
                limits = list(lower = c(R0 = 0, D_inf = 0)))
