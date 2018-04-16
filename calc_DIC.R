#DIC for neg binomial


load(file= "./mcmc.epi.RData")
#mcmc.epi$covmat.empirical

mcmc.epi <- mcmc.epi.pois
trace <- mcmc.epi$trace
mcmc.trace <- mcmc(trace)
summary(mcmc.trace)
tracedf<- data.frame(mcmc.trace)

log.likelihood <- function(fitmodel, theta, init.state, data, log = TRUE) {
  dTrajObs_E(fitmodel, theta, init.state, data, log = TRUE) +
    dTrajObs_I(fitmodel, theta, init.state, data, log = TRUE) +
    dTrajObs_C(fitmodel, theta, init.state, data, log = TRUE)
}

lik<- apply(tracedf[,1:6], 1, function(x) log.likelihood(fitmodel = SIR, theta = x, 
                                                         init.state = c(S = 10000, E = 0, I = 1, C = 0),
                                                         data = epi, log = TRUE))


dic <- function(lik, tracedf) {
  d.bar <- -2*mean(lik)
  theta.bar = apply(tracedf[,1:6], 2, mean)
  d.hat <- -2*log.likelihood(fitmodel = SIR, theta = theta.bar, 
                        init.state = c(S = 10000, E = 0, I = 1, C = 0),
                        data = epi, log = TRUE)
  pD <- d.bar - d.hat
  pV <- var(-2*lik)/2
  DIC=pD+d.bar
  return(DIC)
}

dic(lik, tracedf)


