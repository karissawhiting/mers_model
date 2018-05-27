
SIR$dPointObs_E <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a neg binomial process
  return(dnbinom(x = data.point[["exp"]], size = 10, #dispersion param
                 prob = model.point[["Exp"]],
                 log = log))
}

SIR$dPointObs_E(unlist(epi[2, ]), (unlist(traj[2, ])/10))


dens = 0
for (i in 2:nrow(epi)) {
if(traj[i,][["Exp"]] != 0) {
  data.point <- unlist(epi[i, ])
  model.point <- unlist(traj[i, ])
  print(data.point)
  print(model.point)
  dens <- dens + fitmodel$dPointObs_E(data.point = data.point, 
                                      model.point = model.point, theta = theta, log = TRUE)
  }
}






data.point <- unlist(epi[3, ])
model.point <- unlist(traj[3, ])

dnbinom(x = data.point[["conf"]], size = 10, #dispersion param
        prob = 10/ (10+ model.point[["Con"]]),
        log = log)

dnbinom(x = data.point[["conf"]], size = 10, #dispersion param
        prob = 10/ (10+ model.point[["Con"]]))


#what if zerp 
dnbinom(x = 0, size = 10, #dispersion param
        prob = 0,
        log = log)


dnbinom(x = 0, size = 10, #dispersion param
        prob = .5,
        log = FALSE)






SIR$dPointObs_E <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a neg binomial process
  dpois(x = data.point[["exp"]], #dispersion param
        lambda = model.point[["Exp"]],
        log = log)
  
}


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




