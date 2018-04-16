
#model code that runs
theta<- c(beta = .01,  L = 100 , D0 = .1, D1 = .1, w = .09, n = 10)
init.state <- c(S = 10000, E = 0, I = 1, C = 0)
times <- 1:56
traj <- SIR$simulate(theta, init.state, times)


data = epi
dens= 0
n = 10
for (i in 5:nrow(data)) {
  data.point <- unlist(data[i, ])
  model.point <- unlist(traj[i, ])
  dens <- dens + 
    dnbinom(x = data.point[["exp"]], #dispersion param
            mu = model.point[["Exp"]] + .0001,
            size = (model.point[["Exp"]]/(n-1)),
            log = log)
  print(dens)
}
dpois(x =1, lambda = 0, log = TRUE)
dpois(x =1, lambda = 0, log = FALSE)



dTrajObs_E <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 5:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    dens <- dens + 
      dpois(x = data.point[["exp"]], #dispersion param
            lambda = model.point[["Exp"]] + .0001,
            log = log)
    
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}

dTrajObs_E(SIR, theta, init.state, epi, log = TRUE) 

dTrajObs_log <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 5:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    dens <- dens +
      dpois(x = log(data.point[["exp"]]), #dispersion param
            lambda = log(model.point[["Exp"]] + .0001),
            log = log)
    
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}


dTrajObs_log(SIR, theta, init.state, epi, log = TRUE) 



theta<- c(beta = .99,  L = 8.19 , D0 = 9.26, D1 = 4.05, w = .09, n = 2)
dTrajObs_E(SIR, theta = theta, init.state = init.state, data = epi, log = TRUE)



dnbinom(x = 1,
        mu = (0 + .000001),
        size = ((0+ .000001)/9),
        log = log)
  
dnbinom(x = 1,
        mu = (0 + .00000000000000000001),
        size = ((0+ .000000000000000000001)/9),
        log = log)

dnbinom(x = 1,
        mu = (0),
        size = ((0)/9),
        log = log)


for (i in 1:nrow(epi)) {
  data.point <- unlist(epi[i, ])
  print(data.point[["onset"]])
}

# ISSUE WITH negbinomial errors- errors when zero exp
theta<- c(beta = .8,  L = 7 , D0 = 12, D1 = 4, w = .09, n = 3)
dTrajObs_E <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 5:49) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    print(data.point[["exp"]])
    print(model.point[["Exp"]])
    if(model.point[["Exp"]] > 0){
    dens <- dens + 
      dnbinom(x = data.point[["exp"]],
              mu = (model.point[["Exp"]]),
              size = ((model.point[["Exp"]])/(theta[["n"]]-1)),
              log = log)
    } 
    print(dens)
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}
dTrajObs_E(SIR, theta, init.state, epi, log = TRUE)

traj <- SIR$simulate(theta, init.state, times)

theta < cbind(beta = seq(.6, 1.25, by = .01, L =  )

theta<- c(beta = .8,  L = 7 , D0 = 12, D1 = 4, w = .09, n = 3)
dTrajObs_E <- function (fitmodel, theta, init.state, data, log = TRUE) {
  times <- c(0, data$times)
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  for (i in 5:49) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    if(model.point[["Exp"]] > 0) {
      dens <- sum(dens, 
        (dnbinom(x = data.point[["exp"]],
                mu = (model.point[["Exp"]]),
                size = ((model.point[["Exp"]])/(theta[["n"]]-1)),
                log = log)), na.rm = TRUE)
    }
  }
  return(dens)
  # return(ifelse(log, dens, exp(dens)))
}
dTrajObs_E(SIR, theta, init.state, epi, log = TRUE)



bug<- debug(target = logPosterior_trunc,
            init.theta = c(beta = .6, L =6, D0 = 7, D1 = 5, w = .06, n = 6),
            proposal.sd = proposal.sd,
            n.iterations = 5000,
            adapt.size.start = adapt.size.start ,
            adapt.size.cooling= adapt.size.cooling,
            adapt.shape.start = 100, 
            adapt.shape.stop = 1000,
            limits = list(lower = c(beta = 0, L = 0, D0 = 0, D1 = 0, w = 0, n = 1.000000001)))

#mcmc.epi2 <- mcmcMH(target = logPosterior_trunc,
#                   init.theta = c(beta = .6, L =6, D0 = 7, D1 = 5, w = .06, n = 7),
#                   proposal.sd = c(.03, .8, .8, .8, .02, .8),
#                   n.iterations = 10000,
#                   adapt.size.start = 1000,
#                   adapt.size.cooling=0.999,
#                   limits = list(lower = c(beta = 0, L = 0, D0 = 0, D1 = 0, w = 0, n = )))


#Some trace plot code

#xyplot(x = mcmc.trace.burned)
#densplot(mcmc.trace['beta',])
#MCMCtrace(mcmc.trace)
#MCMCchains(mcmc.trace)

ggplot(data,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.4)

b<- data.frame(mcmcdf[,'beta'][10000:50000])
L <- data.frame(mcmcdf[,'L'][10000:50000])
D0 <- data.frame(mcmcdf[,'D0'][10000:50000])
D1 <- data.frame(mcmcdf[,'D1'][10000:50000])
w <- data.frame(mcmcdf[,'w'][10000:50000])
n <- data.frame(mcmcdf[,'n'][10000:50000])


ggplot(b,aes(x=b)) + 
  geom_density(fill = "pink") + 
  geom_density(aes(L), fill = "blue", alpha = .5) + 
  geom_density(aes(D0), fill = "red", alpha = .5) + 
  geom_density(aes(D1), fill = "green", alpha = .5) +
  geom_density(aes(w), fill = "orange", alpha = .5) + 
  geom_density(aes(n), fill = "yellow", alpha = .5) 

ggplot(w,aes(x=w)) + 
  geom_density()








