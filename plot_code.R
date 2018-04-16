library(reshape2)
library(tidyverse)
library(fields)
library(scatterplot3d)
library(ggplot2)

# Plot Cases -------------------
load("./data-raw/mers_times.RData")

mers_times2<- melt(mers_times[,2:4])

epi <- dcast(mers_times2,  value ~ variable, length) %>%
  rename("times" = value, "exp" = exp_date, "onset" = onset_date, "conf" = conf_date) %>%
  mutate(times = as.numeric(times))

times = data.frame(times = 1:56)
epi = full_join(times, epi, by = "times")

epi[is.na(epi)] <- 0

epi_melt <- melt(epi, id = "times") %>%
  filter(times != -10)

jpeg("./plots/epi_data", width = 8, height = 5, units = "in", res = 300)
ggplot(epi_melt, aes(x = times, y = value, color = variable)) + geom_line(size = 1.5) + 
  geom_point(size = 3) +
  ggtitle("2015 South Korea Outbreak MERS Incidence") + 
  scale_color_manual(name=" ", breaks=c("onset", "conf", "exp"),
                     labels=c("Infected", "Confirmed", "Exposed"),values=c("#C64743",
                             "#7355A3",
                             "#2589A4")) +
  xlab("Days Since Outbreak Start (May 11th)") +
  ylab("Cases")+
  theme_minimal(base_size = 21) +
  theme(legend.justification = c(1.1, 1.1), legend.position = c(1, 1))

dev.off()

# Compare Trajectory with Final Param -------------------

theta<- c(beta = .992,  L = 7.335 , D0 = 7.615, D1 = 4.339, w = .0814, n = 3.889) #8.07 #3.77
init.state <- c(S = 10000, E = 0, I = 1, C = 0) #S = 51413925
times <- 1:57
traj <- SIR$simulate(theta, init.state, times)

jpeg("./plots/exp_data", width = 7.5, height = 3, units = "in", res = 300)
#plot exp
ggplot(data = traj, aes(x = time)) + 
  geom_line(aes(y = Exp), color = "#2589A4", size = 1) +
  #  geom_line(aes(y = Inc), color = "purple" ) +
  #  geom_line(aes(y = Con), color = "orange") +
  geom_point(data = epi, aes(y = exp, x = times), color = "black") +
  labs(y = "") +
  theme_minimal(base_size = 18) + xlab("Days of Outbreak") + xlim(c(0,57)) + 
  ggtitle("Exposed Cases")
dev.off()
#plot I 

jpeg("./plots/inc_data", width = 7.5, height = 3, units = "in", res = 300)
ggplot(data = traj, aes(x = time)) + 
  geom_line(aes(y = Inc), color = "#C64743", size = 1) +
  #  geom_line(aes(y = Inc), color = "purple" ) +
  #  geom_line(aes(y = Con), color = "orange") +
  geom_point(data = epi, aes(y = onset, x = times), color = "black") +
  labs(y = "") +
  theme_minimal(base_size = 18) + xlab("Days of Outbreak") + xlim(c(0,57)) + 
  ggtitle("Infected Cases")
dev.off()

#plot Con
jpeg("./plots/con_data", width = 7.5, height = 3, units = "in", res = 300)
ggplot(data = traj, aes(x = time)) + 
  geom_line(aes(y = Con), color = "#7355A3", size = 1) +
  #  geom_line(aes(y = Inc), color = "purple" ) +
  #  geom_line(aes(y = Con), color = "orange") +
  geom_point(data = epi, aes(y = conf, x = times), color = "black") +
  labs(y = "") + 
  theme_minimal(base_size = 18) + xlab("Days of Outbreak") + xlim(c(0,57)) + 
  ggtitle("Confirmed Cases")
dev.off()

# Map of Cases -------------------
# Read in data
setwd("~/Repositories/mers_model")
mers_raw <- read_csv("./data-raw/KCDC_mers_cl.csv")

# clean data and turn to dates
mers <- mers_raw %>% 
  mutate(onset_date = as.Date(onset_date, format="%m/%d/%y"),
         conf_date = as.Date(conf_date, format="%m/%d/%y"), 
         f_exp_date = as.Date(f_exp_date, format="%m/%d/%y"),
         l_exp_date = as.Date(l_exp_date, format="%m/%d/%y"),
         dis_death_date = as.Date(dis_death_date, format="%m/%d/%y"),
         delay_on_conf =  conf_date - onset_date)

# check for na's
map(mers, ~sum(is.na(.x)))


# Prior Densities -------------------
library(ggplot2);library(reshape2)

x <- data.frame(Beta=rgamma(1000, 1.5, 2.0),
                L=rgamma(1000, 4.44,.55),
                D=rgamma(1000,3.28, .48),
                w=rgamma(1000,2,2), 
                n=rgamma(1000, 3.126, .3125)
)

data<- melt(x)
ggplot(data,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.4) +
  facet_grid(variable ~ ., scale = "free_y") +
  ggtitle("Prior Distributions of Parameters") +
  theme_bw() + 
  scale_fill_discrete(name = "Parameters") +
  theme(legend.justification = c(1.05, 1.05), legend.position = c(1, 1)) + 
  scale_fill_discrete(name=" ",
                       breaks=c("Beta", "L", "D", "w", "n"),
                       labels=c("expression(Beta) ~ G(1.5, 2.0)", 
                                "L ~ G(4.44, .55)", 
                                "D ~ G(3.28, .48)", 
                                "w ~ G(2, 2)", 
                                "n ~  G(3.126, .3125)"))

#gsave("./plots/prior_dist", plot = last_plot(), device = "jpeg", path = NULL, width = 6, height = 3.5,
#       scale = 1)

#png(), insert plot dev.off()


# Posterior Densities -------------------
mcmcdf<- MCMCchains(mcmc.trace)

mcmcdf<-as.data.frame(mcmcdf)[,1:6]

mcmcdfm<- melt(mcmcdf)


ggplot(mcmcdfm ,aes(x=value, fill = variable, alpha = .8)) + 
  geom_density() + 
  facet_grid(variable ~ ., scale = "free_y") + 
  scale_fill_manual(values=c("#E37621",
                             "#C64743",
                             "#7355A3",
                             "#253EA4",
                             "#3766D4", 
                             "#2589A4")) +
  theme_minimal()




# Both Densities -------------------

x <- data.frame(beta=rgamma(1000, 1.5, 2.0),
                L=rgamma(1000, 4.44,.55),
                D0=rgamma(1000,3.28, .48),
                D1=rgamma(1000,3.28, .48),
                w=rgamma(1000,2,2), 
                n=rgamma(1000, 3.126, .3125), 
                dist = rep("Prior", 1000)
)
load(file= "./trace.RData")

mcmcdf<- data.frame(trace)
mcmcdf<-as.data.frame(mcmcdf)[20000:50000,1:6]

mcmcdf<- cbind(mcmcdf, dist = rep("Posterior", nrow(mcmcdf)))


y <- rbind(x, mcmcdf)
mel<- melt(y, by = "dist")

jpeg("./plots/prior_n_post_dist.jpeg", width = 7, height = 8, units = "in", res = 300)
ggplot(mel,aes(x=value, fill = variable, alpha = .9)) + 
  geom_density(size = .7) + 
  xlim(0,20) + 
  facet_grid(variable ~ dist, scale = "free_y") +
  scale_fill_manual(values=c("#E37621",
                             "#C64743",
                             "#2589A4",
                             "#7355A3",
                             "#253EA4",
                             "#3766D4")) +
  theme_minimal(base_size = 21) + theme(legend.position="none")

#ggsave("./plots/prior_n_post_dist", device = "jpeg", path = NULL, width = 4.5, height = 3.5,
#       scale = 1)
dev.off()

### #3D Plots - on burned and thinned
mcmc.trace.burned <- burnAndThin(mcmc.trace, burn = 5000, thin = 5)
mcmcdf<-as.data.frame(mcmc.trace.burned)[,1:6]

mcmcdf2<- unique(mcmcdf[,1:6])
#persp(x = 1:2501, mcmcdf2[,'D0'], mcmcdf2[,'L'], col="lightblue",main="Perspective Plot")
persp(den3d1, box=FALSE)

library(MASS)
library(plotly)

#DO and L
den3d1 <- kde2d(mcmcdf2[,'D0'], mcmcdf2[,'L'])
cols1 <- c("#2589A4",
           "#7355A3",
           "#C64743"
           )
plot_ly(x=den3d1$x, y=den3d1$y, z=den3d1$z) %>% add_surface(colors = cols1) %>% 
  layout(
    scene = list(
      xaxis = list(title = "D0"),
      yaxis = list(title = "L"),
      zaxis = list(title = "Dens")
    ))



# D0 and Beta
den3d2 <- kde2d(mcmcdf2[,'D0'], mcmcdf2[,'beta'])
#persp(den3d2, box=FALSE)


cols1 <- c("#2589A4",
           "#7355A3",
           "#C64743"
)
plot_ly(x=den3d2$x, y=den3d2$y, z=den3d1$z) %>% add_surface(colors = cols1) %>% 
  layout(
    scene = 
    list(
      xaxis = list(title = "D0"),
      yaxis = list(title = "Beta"),
      zaxis = list(title = "Dens")
      ))


# L and Beta
den3d3 <- kde2d(mcmcdf2[,'L'], mcmcdf2[,'beta'])
#persp(den3d2, box=FALSE)

cols1 <- c("#2589A4",
           "#7355A3",
           "#C64743"
)
plot_ly(x=den3d3$x, y=den3d3$y, z=den3d3$z) %>% 
  add_surface(colors = cols1) %>% 
  layout(
    scene = 
      list(
        xaxis = list(title = "L"),
        yaxis = list(title = "Beta"),
        zaxis = list(title = "Dens")
      ))




