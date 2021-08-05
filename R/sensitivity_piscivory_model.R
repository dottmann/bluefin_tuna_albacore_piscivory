
#   sensitivity analysis on visual radius and clearence rate

# Author: Daniel Ottmann
# Created: December 2020
# Last update: May 2021


###########################################
#       Readme

# This script evaluates how sensible is the visual radius and clearance rate of our piscivory model to some of the equation parameters


#######################################################################################################################

##################################################################
# Clear environment:
rm(list = ls())

#############################
# Load packages:
library(tidyverse)

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)


###########################################################################

###########################################################################
# Define parameters:

# Parameter modification values:
parameter_mod <- seq(0.7, 1.3, by = .01)

# Predator and prey SL (m):
predator_length <- 8 / 1000

prey_length <- 4.1 / 1000

# Behavioral anatomical ratio:
behav_anatomical_ratio <- 0.5  # Job & Bellwood 1996

# Larvae only use ~0.5 of their visual area because they mostly look upward; Hilder et al 2019
used_visual_area <- 0.5

# Swimming velocities (m/h):
predator_swimming_velocity = 3 * predator_length * 3600   # Reglero et al 2015
prey_swimming_velocity = 3 * prey_length * 3600   # Reglero et al 2015

# MSA:
msa = 4.699 * ((predator_length * 1000.)^ -1.129)   # Hilder et al 2019

# Proportion of fish shape to a Circle:
prey_shape_proportion <- 0.5 

# Visual radius:
visual_r = prey_shape_proportion * behav_anatomical_ratio * 0.5 * prey_length / tan(msa * 0.5 * pi/180) # Job & Bellwood 1996 JFB 48:952-963, converted from degrees to radians (1 degree =  pi/180 radians)

# Clearance rate:
clearance_rate <- pi * (visual_r^2) * used_visual_area * 
  (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity) # Gertisen &Strickler 1977



###########################################################################
# Prepare sensitivity tables for each parameter:

# Let's start with behav/anatomical ratio:
behav_anatomical_ratio_sensitivity <- expand.grid(parameter_mod = parameter_mod, behav_anatomical_ratio = behav_anatomical_ratio)

behav_anatomical_ratio_sensitivity <- behav_anatomical_ratio_sensitivity %>%
  mutate(behav_anatomical_ratio_s = behav_anatomical_ratio * parameter_mod,
         visual_r = prey_shape_proportion * behav_anatomical_ratio_s * 0.5 * prey_length / tan(msa * 0.5 * pi/180),
         clearance_rate = pi * (visual_r^2) * used_visual_area * 
           (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen & Strickler 1977
         parameter = "behav_anatomical_ratio") %>%
  dplyr::select(parameter_mod, visual_r, clearance_rate, parameter)


# Next comes prey shape proportion:
prey_shape_proportion_sensitivity <- expand.grid(parameter_mod = parameter_mod, prey_shape_proportion = prey_shape_proportion)

prey_shape_proportion_sensitivity <- prey_shape_proportion_sensitivity %>%
  mutate(prey_shape_proportion_s = prey_shape_proportion * parameter_mod,
         behav_anatomical_ratio = .5,
         visual_r = prey_shape_proportion_s * behav_anatomical_ratio * 0.5 * prey_length / tan(msa * 0.5 * pi/180),
         clearance_rate = pi * (visual_r^2) * used_visual_area * 
           (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen & Strickler 1977
         parameter = "prey_shape_proportion_sensitivity") %>%
  dplyr::select(parameter_mod, visual_r, clearance_rate, parameter)


# Next comes used visual area:
used_visual_area_sensitivity <- expand.grid(parameter_mod = parameter_mod, used_visual_area = used_visual_area)

used_visual_area_sensitivity <- used_visual_area_sensitivity %>%
  mutate(used_visual_area_s = used_visual_area * parameter_mod,
         visual_r = prey_shape_proportion * behav_anatomical_ratio * 0.5 * prey_length / tan(msa * 0.5 * pi/180),
         clearance_rate = pi * (visual_r^2) * used_visual_area_s * 
           (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen &Strickler 1977
         parameter = "used_visual_area") %>%
  dplyr::select(parameter_mod, visual_r, clearance_rate, parameter)


# Next comes predator length:
predator_length_sensitivity <- expand.grid(parameter_mod = parameter_mod, predator_length = predator_length)

predator_length_sensitivity <- predator_length_sensitivity %>%
  mutate(predator_length_s = predator_length * parameter_mod,
         predator_swimming_velocity = 3 * predator_length_s * 3600,
         msa = 4.699 * ((predator_length_s * 1000.)^ -1.129),
         visual_r = prey_shape_proportion * behav_anatomical_ratio * 0.5 * prey_length / tan(msa * 0.5 * pi/180),
         clearance_rate = pi * (visual_r^2) * used_visual_area * 
           (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen &Strickler 1977
         parameter = "predator_length") %>%
  dplyr::select(parameter_mod, visual_r, clearance_rate, parameter)


# Next comes prey length:
prey_length_sensitivity <- expand.grid(parameter_mod = parameter_mod, prey_length = prey_length)

prey_length_sensitivity <- prey_length_sensitivity %>%
  mutate(prey_length_s = prey_length * parameter_mod,
         prey_swimming_velocity = 3 * prey_length_s * 3600,
         visual_r = prey_shape_proportion * behav_anatomical_ratio * 0.5 * prey_length_s / tan(msa * 0.5 * pi/180),
         clearance_rate = pi * (visual_r^2) * used_visual_area * 
           (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen &Strickler 1977
         parameter = "prey_length") %>%
  dplyr::select(parameter_mod, visual_r, clearance_rate, parameter)


# Next comes predator swimming velocity:
predator_swimming_velocity_sensitivity <- expand.grid(parameter_mod = parameter_mod, predator_swimming_velocity = predator_swimming_velocity)

predator_swimming_velocity_sensitivity <- predator_swimming_velocity_sensitivity %>%
  mutate(predator_swimming_velocity_s = predator_swimming_velocity * parameter_mod,
         visual_r = prey_shape_proportion * behav_anatomical_ratio * 0.5 * prey_length / tan(msa * 0.5 * pi/180),
         clearance_rate = pi * (visual_r^2) * used_visual_area * 
           (prey_swimming_velocity^2 + 3 * predator_swimming_velocity_s^2) / (3 * predator_swimming_velocity_s), # Gertisen &Strickler 1977
         parameter = "predator_swimming_velocity") %>%
  dplyr::select(parameter_mod, visual_r, clearance_rate, parameter)


# Next comes prey swimming velocity:
prey_swimming_velocity_sensitivity <- expand.grid(parameter_mod = parameter_mod, prey_swimming_velocity = prey_swimming_velocity)

prey_swimming_velocity_sensitivity <- prey_swimming_velocity_sensitivity %>%
  mutate(prey_swimming_velocity_s = prey_swimming_velocity * parameter_mod,
         visual_r = prey_shape_proportion * behav_anatomical_ratio * 0.5 * prey_length / tan(msa * 0.5 * pi/180),
         clearance_rate = pi * (visual_r^2) * used_visual_area * 
           (prey_swimming_velocity_s^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen &Strickler 1977
         parameter = "prey_swimming_velocity") %>%
  dplyr::select(parameter_mod, visual_r, clearance_rate, parameter)


# Next comes msa:
msa_sensitivity <- expand.grid(parameter_mod = parameter_mod, msa = msa)

msa_sensitivity <- msa_sensitivity %>%
  mutate(msa_s = msa * parameter_mod,
         visual_r = prey_shape_proportion * behav_anatomical_ratio * 0.5 * prey_length / tan(msa_s * 0.5 * pi/180),
         clearance_rate = pi * (visual_r^2) * used_visual_area * 
           (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen &Strickler 1977
         parameter = "msa") %>%
  dplyr::select(parameter_mod, visual_r, clearance_rate, parameter)


# Finally comes visual radius:
visual_r_sensitivity <- expand.grid(parameter_mod = parameter_mod, visual_r = visual_r)

visual_r_sensitivity <- visual_r_sensitivity %>%
  mutate(visual_r = visual_r * parameter_mod,
         clearance_rate = pi * (visual_r^2) * used_visual_area * 
           (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen &Strickler 1977
         parameter = "visual_r") %>%
  dplyr::select(parameter_mod, visual_r, clearance_rate, parameter)



# Combine it all in a data frame:
df_sensitivity <- rbind(behav_anatomical_ratio_sensitivity, used_visual_area_sensitivity, predator_length_sensitivity,
                        predator_length_sensitivity, prey_length_sensitivity, predator_swimming_velocity_sensitivity,
                        prey_swimming_velocity_sensitivity, prey_shape_proportion_sensitivity, msa_sensitivity, visual_r_sensitivity)


# Plot it out:

# Some parameters have the same effect, so we group them together.

# Visual Radius:
p <- df_sensitivity %>%
  filter(parameter %in% c("msa", "predator_length", "prey_length")) %>% 
  mutate(parameter = case_when(parameter == "msa" ~ "Minimum separable angle",
                              parameter == "predator_length" ~ "Standard length BFT",
                              parameter == "prey_length" ~ "Standard length ALB,
Behavioral:anatomical ratio 
& larval shape:cicle",
                              T ~ parameter)) %>%
  ggplot() +
  geom_line(aes(x = parameter_mod, y = visual_r, color = parameter, linetype = parameter)) +
  labs(x = "Parameter modification", y = "Visual radius (m)") +
  theme_bw() +
  theme(panel.grid = element_blank())

png(filename = "plots/sensitivity_visual_radius.png", width = 12, height = 5, units = "cm", res = 700)
p
dev.off()


# Clearance rate
p <- df_sensitivity %>%
  filter(parameter %in% c("msa", "predator_length", "predator_swimming_velocity", "visual_r", "prey_swimming_velocity", "used_visual_area", "prey_length")) %>% 
  mutate(parameter = case_when(parameter == "msa" ~ " Minimum separable angle",
                              parameter == "predator_length" ~ "Standard length BFT",
                              parameter == "prey_length" ~ "Standard length ALB",
                              parameter == "visual_r" ~ "Behavioral:anatomical ratio
Larval shape:cicle
& Visual radius",
                              parameter == "predator_swimming_velocity" ~ "Swimming velocity BFT",
                              parameter == "prey_swimming_velocity" ~ "Swimming velocity ALB",
                              parameter == "used_visual_area" ~ "Effective visual area",
                              T ~ parameter)) %>%
  ggplot() +
  geom_line(aes(x = parameter_mod, y = clearance_rate, color = parameter, linetype = parameter)) +
  labs(x = "Parameter modification", y = expression("Clearance rate " (m^3*h^-1*""))) +
  theme_bw() +
  theme(panel.grid = element_blank())

png(filename = "plots/sensitivity_clearance_rate.png", width = 12, height = 5, units = "cm", res = 700)
p
dev.off()


#                        END OF SCRIPT
##################################################################