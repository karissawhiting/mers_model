library('fitR')

# MCMC alg ---------------------------------------
my_mcmcMH <- function(target, init.theta, proposal.sd, n.iterations) {
  target.theta.current<- target(init.theta)
  theta.current <- init.theta
  samples <- theta.current
  accepted <- 0
  
  for (i.iteration in seq_len(n.iterations)) {
    
    #draw a new theta from proposal guassian distribution
    theta.proposed <- rnorm(n = length(theta.current),
                            mean = theta.current,
                            sd = proposal.sd)
    
    #set names of proposed (needed for fitmodel)
    names(theta.proposed) <- names(theta.current)
    
    #evaluate target at the proposed theta 
    target.theta.proposed <- target(theta.proposed)
    
    # compute metropolis-hastings ratio (acceptance probability)
    log.acceptance <- target.theta.proposed - target.theta.current
    
    #draw a random number from uniform
    r<- runif(1)
    
    #must exponentiate because we calculated log
    
    if (r < exp(log.acceptance)) {
      
      #if accepted
      theta.current <- theta.proposed
      
      #update current value of target
      target.theta.current <- target.theta.proposed
      
      #update number of accepted proposals
      accepted <- accepted + 1
      
    }
    
    # add the current theta to the vector of samples
    samples <- rbind(samples, theta.current, deparse.level = 0)
    
    # print current state of the chain and acceptance rate
    message("iteration: ", i.iteration, ", chain:", paste(theta.current, collapse=" "),
            ", acceptance rate:", accepted / i.iteration)
  }
  return(samples)
}
# MCMC alg ---------------------------------------

#make it log norm
dnorm.log <- function(theta) {
    return(dnorm(x = theta, mean = 0, sd = 1, log = TRUE))
  }

#dlogposterior function
my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  
  # calculate the log-likelihood of `theta`
  log.likelihood <- dTrajObs(fitmodel, theta, init.state, data, log = TRUE)
  
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}
#adjust log posterio function to jsut take R0 as a value
my_dLogPosterior_R0_epi1 <- function(R0) {
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = c(R0 = R0, D_inf = 2),
                          init.state = c(S = 999, I = 1, R = 0),
                          data = epi1))
}

my_dLogPosterior_R0_epi1(R0 = 3)

#initiate MCMC
#starting.value <- 100 # starting value for MCMC
#sigma <- 1 # standard deviation of MCMC
#iter <- 10000
trace <- my_mcmcMH(target = my_dLogPosterior_R0_epi1,
                   init.theta = 1.3,
                   proposal.sd = .04, n.iterations = 1000)

plot(trace, type = 'l')
hist(trace, freq = FALSE, breaks = 1000)






