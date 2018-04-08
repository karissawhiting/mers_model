my_mcmcMH <- function(target, init.theta, proposal.sd, n.iterations) {
  
  # evaluate the function "target" at "init.theta", and assign to
  # a variable called target.theta.current.
  target.theta.current <- target(init.theta)
  
  # initialise variables to store the current value of theta, the
  # vector of samples, and the number of accepted runs
  theta.current <- init.theta
  samples <- theta.current
  accepted <- 0
  
  # run MCMC for n.iteration interations
  for (i.iteration in seq_len(n.iterations)) {
    
    # draw a new theta from the (Gaussian) proposal distribution
    # and assign to a variable called theta.proposed.  
    # See "?rnorm for more information
    # Note that this step is vectorized for any arbitratry theta 
    # which will be useful when we will sample from a multivariate
    # target distribution
    theta.proposed <- rnorm(n = length(theta.current),
                            mean = theta.current,
                            sd = proposal.sd)
    print(paste("theta.proposed: ", theta.proposed))
    
    # Note that 'rnorm' returns an unnamed vector, but the functions of
    # 'fitmodel' need a named parameter vector. We therefore set
    # the names of theta.proposed to be the same as the names of
    # theta.current
    names(theta.proposed) <- names(theta.current)
    
    # evaluate the function target at the proposed theta and
    # assign to a variable called target.theta.proposed
    target.theta.proposed <- target(theta.proposed)
    print(paste("target.theta.proposed: ", target.theta.proposed))
    
    # compute Metropolis-Hastings ratio (acceptance probability). Since
    # the multivariate Gaussian is symmetric, we don't need to consider
    # the proposal distribution here
    log.acceptance <- target.theta.proposed - target.theta.current
    
    # draw random number number between 0 and 1 using "runif" and assign to
    # a variable called r.
    r <- runif(1)
    
    # test acceptance by comparing the random number to the
    # Metropolis-Hastings ratio (acceptance probability) (using
    # "exp" because we calculated the logarithm of the
    # Metropolis-Hastings ratio before)
    
    print(exp(log.acceptance))
    print(paste("log.accep: : ",log.acceptance))
    
    if (r < exp(log.acceptance)) {
      
      # if accepted:
      # change the current value of theta to the proposed theta
      theta.current <- theta.proposed
      
      # updated the current value of the target
      target.theta.current <- target.theta.proposed
      
      # update number of accepted proposals
      accepted <- accepted + 1
    }
    
    # add the current theta to the vector of samples
    # Note that we use `rbind` in order to deal with multivariate 
    # target. So if `theta` is a vector then `samples` is a matrix.
    samples <- rbind(samples, theta.current, deparse.level=0)
    
    # print current state of chain and acceptance rate
    # use paste() to deal with the case where `theta` is a vector
    message("iteration: ", i.iteration, ", chain:", paste(theta.current, collapse=" "),
            ", acceptance rate:", accepted / i.iteration, ", error:", exp(log.acceptance))
    
  }
  
  # return the trace of the chain (i.e., the vector of samples)
  return(samples)
}
