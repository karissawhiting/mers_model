# Plot Cases -------------------
load("./data-raw/mers_times.RData")

mers_times2<- melt(mers_times[,2:4])

epi <- dcast(mers_times2,  value ~ variable, length) %>%
  rename("times" = value, "exp" = exp_date, "onset" = onset_date, "conf" = conf_date) %>%
  mutate(times = as.numeric(times))

times = data.frame(times = 1:56)
epi = full_join(times, epi, by = "times")

epi[is.na(epi)] <- 0

epi_melt <- melt(epi, id = "times")

ggplot(epi_melt, aes(x = times, y = value, color = variable)) + geom_line() + 
  geom_line(aes(I, traj[1:56,]))


##########
ggplot(data = traj, aes(x = time)) + 
  geom_line(aes(y = S), color = "blue") + 
  geom_line(aes(y = E), color = "purple") +
  geom_line(aes(y = I), color = "green") +
  geom_line(aes(y = C), color = "orange") +
  geom_line(aes(y = Exp), color = "green") +
#  geom_line(aes(y = Inc), color = "purple" ) +
#  geom_line(aes(y = Con), color = "orange") +
  geom_point(data = epi, aes(y = onset, x = times), color = "red") +
  labs(y = "") 

#plot I against I
ggplot(data = traj, aes(x = time)) + 
  geom_line(aes(y = Exp), color = "green") +
  #  geom_line(aes(y = Inc), color = "purple" ) +
  #  geom_line(aes(y = Con), color = "orange") +
  geom_point(data = epi, aes(y = onset, x = times), color = "red") +
  labs(y = "") 
