# Data Clean -----------------------------

load("./data-raw/mers_times.RData")

mers_times2<- melt(mers_times[,2:4])

epi <- dcast(mers_times2,  value ~ variable, length) %>%
  rename("times" = value, "exp" = exp_date, "onset" = onset_date, "conf" = conf_date) %>%
  mutate(times = as.numeric(times))

times = data.frame(times = 1:56)
epi = full_join(times, epi, by = "times")

epi[is.na(epi)] <- 0

rm(mers_times, mers_times2, times)



# Coding Model Expanded  -----------------------------
SIR$name <- c("SEIC model with w, D0, D1 param")
SIR$state.names <- c("S", "E", "I", "C", "Exp", "Inc", "Con")
SIR$theta.names <- c("beta", "L", "D0", "D1", "w")


SIR$simulate <- function (theta, init.state, times) 
{
  SIR_ode <- function(time, state, parameters) {
    
    beta <- parameters[["beta"]]
    L <- parameters[["L"]]
    D0 <- parameters[["D0"]]
    D1<- parameters[["D1"]]
    w <- parameters[["w"]]
    
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C <- state[["C"]]
    Exp <- state[["Exp"]]
    Inc <- state[["Inc"]]
    Con <- state[["Con"]]
    
    N <- S + E + I + C
    
    if (time < 18) {
      dS <- -1*beta * S * I/N
      S1 = max(S,0) #prevent from going below zero
      dE <- 1*beta * (S1 * I)/N - (E/L)
      
      dI <- (E/L) - (I/D0)
      dC <- I/D0
      
      dExp <- beta * S1 * I/N
      dInc <- (E/L)
      dCon <- (I/D0)
      
    }
    else {
      dS <- -w*beta * S * I/N
      S1 = max(S,0) #prevent from going below zero
      dE <- (w*beta * (S1 * I)/N) - (E/L)
      
      dI <- (E/L) - (I/D1)
      dC <- I/D1
      
      dExp <- w*beta * S1 * I/N
      dInc <- (E/L)
      dCon <- (I/D1)
    }
    
    return(list(c(dS, dE, dI, dC, dExp, dInc, dCon)))
  }
  
  # put incidence at 0 in init.state
  init.state["Exp"] <- 0
  init.state["Inc"] <- 0
  init.state["Con"] <- 0
  
  traj <- data.frame(ode(y = init.state, times = times, method = "ode45",
                         func = SIR_ode, parms = theta)) #events = list(data = eventdat))) #ode45
  
  # compute incidence of each time interval
  traj$Exp <- c(0, diff(traj$Exp))
  traj$Inc <- c(0, diff(traj$Inc))
  traj$Con <- c(0, diff(traj$Con))
  
  return(traj)
}



theta<- c(beta = .99,  L = 8.19 , D0 = 9.26, D1 = 4.05, w = .09) #8.07 #3.77
init.state <- c(S = 10000, E = 0, I = 1, C = 0) #S = 51413925
times <- 1:57
traj <- SIR$simulate(theta, init.state, times)

#plot conf data against model
ggplot(data = traj, aes(x = time)) + 
  geom_line(aes(y = Con), color = "green") +
  #  geom_line(aes(y = Inc), color = "purple" ) +
  #  geom_line(aes(y = Con), color = "orange") +
  geom_point(data = epi, aes(y = conf, x = times), color = "red") +
  labs(y = "") + xlim(10,57) +
  ggtitle("Confirmed Cases, Start May 20th")

#plot Inc  data against model
ggplot(data = traj, aes(x = time)) + 
  geom_line(aes(y = Inc), color = "green") +
  #  geom_line(aes(y = Inc), color = "purple" ) +
  #  geom_line(aes(y = Con), color = "orange") +
  geom_point(data = epi, aes(y = onset, x = times), color = "red") +
  labs(y = "") + xlim(1,67) + 
  ggtitle("Incidence, Start May 11th")


#plot exposed data against model 
ggplot(data = traj, aes(x = time)) + 
  geom_line(aes(y = Exp), color = "green") +
  #  geom_line(aes(y = Inc), color = "purple" ) +
  #  geom_line(aes(y = Con), color = "orange") +
  geom_point(data = epi, aes(y = exp, x = times), color = "red") +
  labs(y = "") + xlim(4,57) +
  ggtitle("Exposed Cases, Start May 13th")


