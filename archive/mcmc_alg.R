trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 31

x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
y <- trueA * x + trueB + rnorm(n = sampleSize, mean = 0, sd = trueSd)
plot(x,y, main = "Test Data")

likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  
  pred = a*x + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

slopevalues <- function(x){return(likelihood(c(x, trueB, trueSd)))}
slopelikelihoods <- lapply(seq(3, 7, by=.05), slopevalues )
plot (seq(3, 7, by=.05), slopelikelihoods , type="l", xlab = "values of slope parameter a", ylab = "Log likelihood")


prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min=0, max=10, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+sdprior)
}

#The product of the prior and likelihood is the quantity MCMC will be working on

posterior <- function(param){
  return (likelihood(param) + prior(param))
}


######## Metropolis algorithm ################
proposalfunction <- function(param){
  return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue = c(4,0,10)
chain = run_metropolis_MCMC(startvalue, 10000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))



######## Another Version ################
#The wolf pups dataset
x <- c(5, 8, 7, 5, 3, 4, 3, 9, 5, 8, 5, 6, 5, 6, 4, 7)

# The log posterior density of the Laplace distribution model, when assuming 
# uniorm/flat priors. The Laplace distribution is not part of base R but is
# available in the VGAM package.
model <- function(pars) {
  sum(VGAM::dlaplace(x, pars[1], exp(pars[2]), log = TRUE))
}

# The Metropolis-Hastings algorithm using a Uniform(-0.5, 0.5) proposal distribution 
metrop <- function(n_samples, model, inits) {
  samples <- matrix(NA, nrow = n_samples, ncol = length(inits))
  samples[1,] <- inits
  for(i in 2:n_samples) {
    curr_log_dens <- model(samples[i - 1, ])
    proposal <- samples[i - 1, ] + runif(length(inits), -0.5, 0.5)
    proposal_log_dens <- model(proposal)
    if(runif(1) < exp(proposal_log_dens - curr_log_dens)) {
      samples[i, ] <- proposal
    } else {
      samples[i, ] <- samples[i - 1, ]
    }
  }
  samples
}

samples <- metrop(n_samples = 1000, model, inits = c(0,0))
# Plotting a traceplot
plot(samples[,1], type = "l", ylab = expression(Location ~ mu), col = "blue")
# Calculating median posterior and 95% CI discarding the first 250 draws as "burnin".
quantile(samples[250:1000,1], c(0.025, 0.5, 0.975))

