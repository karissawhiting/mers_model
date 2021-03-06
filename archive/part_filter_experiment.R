# This is a function that takes four parameters:
# - fitmodel: a fitmodel object
# - theta: named numeric vector. Values of the parameters for which the marginal log-likelihood is desired.
# - init.state: named numeric vector. Initial values of the state variables.
# - data: data frame. Observation times and observed data.
# The function returns the value of the marginal log-likelihood
my_particleFilter2 <- function(fitmodel, theta, init.state, data, n.particles) {
  
  ## Initialisation of the algorithm
  
  # Marginal log-likelihood is set to 0 and will be updated during the filtering steps
  margLogLike <- 0
  
  # Particle states can be stored in a list
  state.particles  <- rep(list(init.state), n.particles)
  
  # Weight: initially equal for all the particles 
  # particle weight can be stored in a vector
  weight.particles <- rep(1/n.particles, length = n.particles)
  
  # Initialise time variable
  current.time <- 0

  ## Loop over observation times: resample, propagate, weight
  for(i in seq_len(nrow(data))){
    
    # Extract next data point (must be a vector)
    data.point <- unlist(data[i, ])
    next.time <- data.point["times"]
    
    # Resample particles according to their weights. 
    # You can use the `sample` function of R
    # (normalisation of the weights is done in the function)
    index.resampled <- sample(x = n.particles,
                              size = n.particles,
                              replace = TRUE,
                              prob = weight.particles)
    state.particles <- state.particles[index.resampled]
    print(state.particles)
    ## Loop over particles: propagate and weight

    for(p in 1:n.particles){
      
      # Extract current state of the particle 
      current.state.particle <- state.particles[[p]]
      
      print(list(current.state.particle, theta, init.state, current.time, next.time))

      # Propagate the particle from current observation time 
      # to the next one using the function `fitmodel$simulate`
      traj <- fitmodel$simulate(theta = theta,
                                init.state = current.state.particle,
                                times = c(current.time,next.time))
      

      # Extract state of the model at next observation time
      # Also make sure that model.point is a vector
      model.point <- unlist(traj[2,fitmodel$state.names])
      
      print(paste("modpoint", model.point))
      print(paste("datpoint", data.point))
      # Weight the particle with the likelihood of the observed 
      # data point using the function `fitmodel$dPointObs`
      weight.particles[p] <-
        fitmodel$dPointObs(data.point = data.point +1,
                           model.point = model.point+1,
                           theta = theta)
      print(weight.particles)
      # Update state of the p particle
      state.particles[[p]] <- model.point
      
    }
    
    # Increment time
    current.time <- next.time
    
    ## Increment the marginal log-likelihood
    # Add the log of the mean of the particles weights
    margLogLike <- margLogLike + log(mean(weight.particles))
  }
  
  ## Return marginal log-likelihood
  return(margLogLike)
  
}


my_particleFilter2(SIR, theta, init.state, epi, n.particles = 20)





data(SEIT4L_stoch)

# load data
data(FluTdC1971)

# theta close to the mean posterior estimate of the deterministic SEIT4L
# model
theta <- c(R0 = 7, D_lat = 1, D_inf = 4, alpha = 0.5, D_imm = 10, rho = 0.65)

# init state as before
init.state <- c(S = 279, E = 0, I = 2, T1 = 3, T2 = 0, T3 = 0, T4 = 0, L = 0, 
                Inc = 0)

# run the particle filter with 20 particles
my_particleFilter(SEIT4L_stoch, theta, init.state, data = FluTdC1971, n.particles = 20)
## [1] -124.5306
#data(SIR)
