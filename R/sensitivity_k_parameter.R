
#   Larval tuna piscivory

# Author: Daniel Ottmann
# Created: June 2020
# Last update: May 2021


###########################################
#       Readme

# This is a sensibility analysis to evaluate how the parameter k of predator capture success affects the predator model


#######################################################################################################################

##################################################################
# Clear environment:
rm(list = ls())


#############################
# Load packages:
library(tidyverse)
library(cowplot)


#################################################################################

#################################################################################
# Define parameters:

# Hours of light:
# Obtained from https://www.timeanddate.com/sun/spain/palma?month=7
# Daylight in July 1st: 6:25-21:20 -> 15 hours
light <- 15

# Proportion of fish shape to a Circle:
prey_shape_proportion <- 0.5 

# Behavioral anatomical ratio:
behav_anatomical_ratio <- 0.5  # Job & Bellwood 1996

# Larvae only use ~0.5 of their visual area because they mostly look upward; Hilder et al 2019
used_visual_area <- 0.5

# Set prey and predator cut-offs:
prey_cutoff <- 4.1 
predator_cutoff <- 7.5 

# Coeficient for the probability of capture success:
k <- 10


#######################################################################################

#######################################################################################


############################
# Sensitivity analysis for k:

# First evaluate the probability of capture success as a funtion of predator:prey length ratio:
# Set predator size:
SLbft <- predator_cutoff - .1

# Set prey sizes:
SLalb <- c(seq(2.5, predator_cutoff - .1, by = .1))

# Set k values:
k <- c(0, 1, 3, 6, 10, 15, 20)

# Put all combinations in a table:
df_k <- expand.grid(SLalb = SLalb, SLbft = SLbft, k = k)

# Edit the table to get probability of capture success as a function of relative size proportion:
df_k <- df_k %>%
  mutate(real_length_proportion = SLalb/SLbft,
         length_proportion = (SLalb - (prey_cutoff)) / (SLbft - ((prey_cutoff))),
         P = (1 - length_proportion)^k,
         k = as.factor(k)) %>%
  filter(length_proportion >= 0)


# Plot it out:
p <- ggplot() +
  geom_line(data = df_k, aes(x = length_proportion, y = P, color = k, linetype = k)) +
  labs(y = "Probability of capture success", x = "Prey:Predator length ratio", linetype = "k", color = "k") +
  theme_bw() +
  theme(panel.grid = element_blank())


#########################################################################################
# Now evaluate how each k affects survival of a simulated prey population:

# Make a list and a simulation id:
plot_list <- list()
id <- 0

# Loop the simulation for each k value:
for (i in k) {
  
  # Ad 1 to the simulation id:
  id <- id + 1
  
  # Simulate the prey population:
  sim_prey <- data.frame(length_fresh = seq(3.4, predator_cutoff - .1, by = .1),
                         density = c(0.01, 0.035, 0.07, 0.115, 0.15, 0.17, 0.1775, .181, 0.18, 0.1785, 0.175, 0.17,
                                     0.1625, 0.1525, 0.14, 0.125, 0.113, 0.101, 0.09, 0.08, 0.071, 0.063, 0.056, 
                                     0.048, 0.0415, 0.036, 0.031, 0.027, 0.023, 0.02, 0.017, 0.014, 0.01, 0.008,
                                     0.007, 0.0065, 0.006, 0.0055, 0.00527,0.00515, 0.005))
  
  # Simulate predators:
  sim_predator <- data.frame(length_fresh = c(8, 10),
                             density = 0.005)
  
  # Combine both data frames:
  simulation <- rbind(sim_prey, sim_predator)
  
    # Make some shenanigans to calculate predation of each predator over each prey cohort:
  simulation <- simulation %>%
    filter(length_fresh <= predator_cutoff) %>% 
    slice(rep(row_number(), nrow(sim_predator))) %>%
    arrange(length_fresh) %>%
    mutate(prey_length = length_fresh / 1000,
           prey_density = density,
           predator_length = rep(sim_predator$length_fresh / 1000, n() / nrow(sim_predator)),
           predator_density = rep(sim_predator$density, n() / nrow(sim_predator)),
           predator_swimming_velocity = 3 * predator_length * 3600,
           prey_swimming_velocity = 3 * prey_length * 3600,
           msa = 4.699 * ((predator_length * 1000.)^ -1.129),  # Hilder et al (2019) 
           visual_r = prey_shape_proportion * 0.5 * prey_length / tan(msa * 0.5 * pi/180), # Job & Bellwood 1996 JFB 48:952-963, converted from degrees to radians (1 degree =  pi/180 radians)
           visual_r = behav_anatomical_ratio * visual_r,
           length_proportion = (prey_length - (prey_cutoff / 1000)) / (predator_length - ((prey_cutoff) / 1000)),
           probability_capture_success = case_when(prey_length <= (prey_cutoff / 1000) ~ 1,
                                                   T ~ (1 - length_proportion)^i),
           encounter_rate_h = pi * (visual_r^2) * used_visual_area * predator_density * 
             (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen & Strickler 1977
           survival_h = exp(-encounter_rate_h * probability_capture_success),
           survival_d = exp(-encounter_rate_h * probability_capture_success * sum(light)),
           ingestion_rate_h = encounter_rate_h * probability_capture_success * prey_density) %>% 
    as_data_frame() %>%
    dplyr::select(-length_fresh, -density)
  
  # Make another table for the surviving prey:
  simulation2 <- simulation %>%
    dplyr::select(prey_length, predator_length, prey_density, survival_d) %>%
    group_by(prey_length, prey_density) %>%
    summarise(total_survival_d = prod(survival_d)) %>%
    mutate(prey_density_24 = prey_density * total_survival_d)
  
    # Make a plot:
  p2 <- ggplot(data = simulation) +
    geom_col(aes(x = prey_length * 1000, y = prey_density / nrow(sim_predator)), width = .085) +
    lims(x = c(2.5, 7.5)) + # 10.5
    geom_line(aes(x = prey_length * 1000, y = survival_d * .4/1, colour = as.factor(predator_length * 1000)), alpha = .85) +
    scale_color_discrete(name = paste0("k = ",i,"
Predator (mm)")) +
    scale_y_continuous(sec.axis = sec_axis(~ .* (1/.4), name = "Probability of surviving one day"), limits = c(0, .4)) +
    labs(x = "Standard length (mm)", y = expression("Larvae m"^-3*"")) +
    geom_col(data = simulation2, aes(x = prey_length * 1000, y = prey_density_24), fill = "black", width = .085) +
    geom_vline(xintercept = c(prey_cutoff), linetype = "dotted", alpha = .5) + # , predator_cutoff
    theme_bw() +
    theme(panel.grid = element_blank())
  
  # Put each plot in a list:
  plot_list[[id]] <- p2

}

# Plot it out:
png(filename = "plots/sensitivity_k.png", width = 21, height = 30, units = "cm", res = 300)
plot_grid(p, plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]], ncol = 2,labels = "auto")
dev.off()


#########################################################################################
# With  k = 10, how does predator:prey size ratio affect probability of capture success?

# Set k value and predator and prey lengths:
k <- 10
prey_length <- seq(3.4, predator_cutoff - .1, by = .1) / 1000
predator_length <- seq(predator_cutoff, 13, by = .1)  / 1000

# Put all predator-prey size combinations in a table:
df <- expand.grid(prey_length = prey_length, predator_length = predator_length, k = k)

# Edit table:
df <- df %>%
  mutate(predator_prey_ratio = predator_length / prey_length,
         length_proportion = (prey_length - (prey_cutoff / 1000)) / (predator_length - ((prey_cutoff) / 1000)),
         probability_capture_success = case_when(prey_length <= (prey_cutoff / 1000) ~ 1,
                                                 T ~ (1 - length_proportion)^k)) 

# Set a range of predator:prey size ratios:
predator_prey_ratio <- seq(0, 10, by = .01)

# Calculate probability of capture success by Miller et al (1988):
probability_capture_success <- (100-((predator_prey_ratio + 3.37)/44.7)^-2.28) / 100

# Make a data frame with Miller's equation:
miller_df <- data.frame(predator_prey_ratio = predator_prey_ratio, probability_capture_success = probability_capture_success)
miller_df <- miller_df %>%
  filter(probability_capture_success >= 0)

# Plot it all otgether:
p <- ggplot() +
  geom_point(data = df, aes(x = predator_prey_ratio, y = probability_capture_success, color = prey_length * 1000), alpha = .8) +
  geom_smooth(data = df, aes(x = predator_prey_ratio, y = probability_capture_success), color = "black", se = F) +
  geom_line(data = miller_df, aes(x = predator_prey_ratio, y = probability_capture_success), color = "red") +
  xlim(c(0, 10)) +
  ylim(c(0, 1)) +
  ylab("Probability of capture success") +
  xlab("Predator:Prey size ratio") +
  labs(color = "k = 10
Albacore
SL (mm)") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Outfile plot:
png(filename = "plots/sensitivity_capture_success_vs_predator_prey_ratio_k10.png", width = 10, height = 7, units = "cm", res = 300)
p
dev.off()


#                                          END OF SCRIPT
#################################################################################################