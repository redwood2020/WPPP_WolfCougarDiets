# Lauren Satterfield adapted from Briana Abrahms' Code

#load packages
library(tidyverse)
library(graphics)
library(gridExtra)
library(quantreg)

# load relevant data from bootstrapping run (Oct19.22 was a bootstrap of 2000 instances)
load(file="Oct21.22_bs10000.RData")

# prepare Shannon's values for plotting ----
# Shannon's H index - diversity ------------------------------

# Shannon's H value - dietary diversity
# means
cougars_H_mean <- shannon_H_plot[c(1,3,5,7),2] #cougar means
wolves_H_mean <- shannon_H_plot[c(2,4,6,8),2] #wolf means
# lower 95%
cougars_H_low <- shannon_H_plot[c(1,3,5,7),3]
wolves_H_low <- shannon_H_plot[c(2,4,6,8),3]
# upper 95%
cougars_H_high <- shannon_H_plot[c(2,4,6,8),4]
wolves_H_high <- shannon_H_plot[c(2,4,6,8),4]

# Shannon's J index - evenness ------------------------------

# Shannon's J value - dietary evennness
# means
cougars_J_mean <- shannon_J_plot[c(1,3,5,7),2] #cougar means
wolves_J_mean <- shannon_J_plot[c(2,4,6,8),2] #wolf means
# lower 95%
cougars_J_low <- shannon_J_plot[c(1,3,5,7),3]
wolves_J_low <- shannon_J_plot[c(2,4,6,8),3]
# upper 95%
cougars_J_high <- shannon_J_plot[c(1,3,5,7),4]
wolves_J_high <- shannon_J_plot[c(2,4,6,8),4]


# pianka's means
overlap_mean <- pianka_plot[,2] #pianka means
overlap_low <- pianka_plot[,3] #pianka lower
overlap_high <- pianka_plot[,4] #pianka upper

wc_shannons_H = data.frame(Species = c(rep("Cougar",4), rep("Wolf",4)),
                  diversity = c(cougars_H_mean, wolves_H_mean),
                  low_d = c(cougars_H_low, wolves_H_low),
                  high_d = c(cougars_H_high, wolves_H_high),
                  overlap = rep(overlap_mean,2))

wc_shannons_J = data.frame(Species = c(rep("Cougar",4), rep("Wolf",4)),
                 evenness = c(cougars_J_mean, wolves_J_mean),
                 low_e = c(cougars_J_low, wolves_J_low),
                 high_e = c(cougars_J_high, wolves_J_high),
                 overlap = rep(overlap_mean,2))

# Shannon's H - diversity
plot(cougars_H_mean, overlap_mean, xlim=c(0,1))
plot(wolves_H_mean, overlap_mean, xlim=c(0,1))

# Shannon's J - evennness
plot(cougars_J_mean, overlap_mean, xlim=c(0,1))
plot(wolves_J_mean, overlap_mean, xlim=c(0,1))
# ---------------------------------------------------

# plot Shannon's H 
# linear model
m <- lm(overlap~diversity, data = wc_shannons_H)
# 95% quantile, 2 tailed
Sd_low <- rq(overlap ~ diversity, data = wc_shannons_H, tau = 0.025) #lower quantile
Sd_high <- rq(overlap ~ diversity, data = wc_shannons_H, tau = 0.975) #upper quantile

# plot Shannon's H - diversity ----
wc_shannons_H <- wc_shannons_H %>% 
  mutate(low_line_d = Sd_low$coefficients[1] + diversity * Sd_low$coefficients[2],
         high_line_d = Sd_high$coefficients[1] + diversity * Sd_high$coefficients[2]) 

shan_H_plot <- ggplot(wc_shannons_H, aes(x=diversity, y=overlap)) +
  geom_point(aes(shape=Species), size=6, color='blue')+
  geom_smooth(aes(diversity, overlap),method = 'lm', se=FALSE, color='blue') +
  geom_line(aes(diversity, low_line_d), linetype = "dashed", color='blue') +
  geom_line(aes(diversity, high_line_d), linetype = "dashed", color='blue') +
  theme_classic() +
  theme(text = element_text(size=15)) +
  scale_x_continuous(limits = c(0,1)) +
  xlab("\n Shannon's Diversity Index (H)") +
  ylab("Pianka's Index (mean) \n") +
  annotate(geom="text", x=0.8, y=1.0, label="r^2 = 0.5623; p = 0.0322 \n y = -0.2276*x + 0.9742", color="blue", size = 5)

shan_H_plot

summary(lm(overlap~diversity, data=wc_shannons_H))

# plot Shannon's J index - evenness ----------------------------

# linear model
m <- lm(overlap ~ evenness, data = wc_shannons_J)
# 95% quantile, 2 tailed
Se_low <- rq(overlap ~ evenness, data = wc_shannons_J, tau = 0.025) #lower quantile
Se_high <- rq(overlap ~ evenness, data = wc_shannons_J, tau = 0.975) #upper quantile

# Shannon's J - evennness
wc_shannons_J <- wc_shannons_J %>% 
  mutate(low_line_e = Se_low$coefficients[1] + evenness * Se_low$coefficients[2],
         high_line_e = Se_high$coefficients[1] + evenness * Se_high$coefficients[2]) 

shan_J_plot <- ggplot(wc_shannons_J, aes(x=evenness, y=overlap)) +
  geom_point(aes(shape=Species), size=6, color='darkred')+
  geom_smooth(aes(evenness, overlap),method = 'lm', se=FALSE, color='darkred') +
  geom_line(aes(evenness, low_line_e), linetype = "dashed", color='darkred') +
  geom_line(aes(evenness, high_line_e), linetype = "dashed", color='darkred') +
  theme_classic() +
  theme(text = element_text(size=15)) +
  scale_x_continuous(limits = c(0,1)) +
  xlab("\n Shannon's Evenness Index (J)") +
  ylab("Pianka's Index (mean) \n")+
  annotate(geom="text", x=0.8, y=1.0, label="r^2 = 0.5650; p = 0.0315 \n y = -0.2376*x + 0.9746", color="darkred", size = 5)

summary(lm(overlap~evenness, data=wc_shannons_J))

# plot Shannon's H and Shannon's J in two panels of same plot ----
grid.arrange(shan_H_plot, shan_J_plot, nrow = 2)
