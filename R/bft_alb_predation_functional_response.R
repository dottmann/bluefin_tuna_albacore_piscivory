
#   Size-dependent functional response of larval bluefin tuna preying on albacore

# Author: Daniel Ottmann
# Created: May 2020
# Last update: May 2021


###########################################
#       Readme

# This script evaluates functional response of predation of larval bluefin tuna on larval albacore
# We simulate range of densities of post-flexion predators (using the real distribution) and evaluate their 
# effect on observed populations and on a virtual population of prey


##################################################################

##################################################################
# Clear environment:
rm(list = ls())


#############################
# Load packages:
library(tidyverse)
library(grid)
library(png)

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)


##################################################################
# Load data
load("data/data.RData")
load("data/id_operations.RData")


# Load images:
img_ysl <- readPNG("media/ysl.png", native = T)
img_preflex <- readPNG("media/preflex.png", native = T)


##################################################################

##################################################################


#########################
# Define model constants:


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


###############################################################################

###############################################################################

# Simulate different predatory pressures:

# Filter the size-densities of predators:
predator_size_distribution <-  data %>%
  filter(species == "Thunnus thynnus", length_fresh >= predator_cutoff)


# What is the range of predator density across size-classes?
range(predator_size_distribution$mean_density)
# Cool, so we want to set 0.032 as the greatest density and 0.001 as the lowest density

# Bin predator size classes by 0.1 mm:
range_lengths <- seq(predator_cutoff, 13, by = 0.1)

# Make an output vector for the for-loop:
density <- c()

# Simulate densities following an exponential decline from 0.032 to 0.001:
for (i in 1:length(range_lengths) - 1) {
  density <- rbind(density, print(round(.0005 + 0.0315 * .926^i, 5)))
}


# Make a table with the density of each simulated size class:
sim_predator <- data.frame(length_fresh = range_lengths,
                       density = density)


# Now plot observations on top of simulated predators:
p <- ggplot() +
  geom_col(data = sim_predator, aes(x = length_fresh, y = density), show.legend = FALSE, width = .09, fill = "darkorange1") +
  labs(x = "Standard length (mm)", y = expression("BFT density (larvae m"^-3*")")) +
  geom_point(data = predator_size_distribution, aes(x = length_fresh, y = mean_density), size = 1, alpha = .5) +
  theme_bw() +
  theme(panel.grid = element_blank())

png(filename = "plots/functional_response_predator_abundances.png", width = 12, height = 8, units = "cm", res = 400)
p
dev.off()


# The greatest density of predators found in a single station (including all class sizes) was 0.167985 BFT / m^3
0.167985 / sum(density)

# Set values of maximum predator pressure:
# Extreme = 0.375; High = 0.18; Mid =  0.045; low = 0.005
x <- c(0.005, 0.045, 0.18, 0.375)

# Run a for loop for each maximum predator pressure:
for(j in x) {
  
  # First, simulate the predator population:
  pressure <- seq(0, j, by = .005)  
  
  sim_predator <- data.frame(length_fresh = range_lengths,
                             density = density)
  
  sim_predator <- sim_predator %>% 
    slice(rep(row_number(), length(pressure))) %>%
    mutate(pressure = rep(pressure, each = length(range_lengths)),
           mean_density = density * pressure) %>%
    dplyr::select(-density)

  # Now simulate the prey population:
  sim_prey <- data.frame(length_fresh = seq(3.4, predator_cutoff - .1, by = .1),
                         mean_density = c(0.01, 0.035, 0.07, 0.115, 0.15, 0.17, 0.1775, .181, 0.18, 0.1785, 0.175, 0.17,
                                          0.1625, 0.1525, 0.14, 0.125, 0.113, 0.101, 0.09, 0.08, 0.071, 0.063, 0.056, 
                                          0.048, 0.0415, 0.036, 0.031, 0.027, 0.023, 0.02, 0.017, 0.014, 0.01, 0.008,
                                          0.007, 0.0065, 0.006, 0.0055, 0.00527,0.00515, 0.005))
  
  # Prepare an ouput table:
  df1 <- data.frame()
  

  # Run another forloop:
  for (i in unique(sim_predator$pressure)) {
    
    sim_predator_i <- sim_predator %>%
      filter(pressure == i) %>%
      dplyr::select(-pressure)
    
    
    
    # Combine both data frames of predators and prey:
    simulation <- rbind(sim_prey, sim_predator_i)
    
    
    # Make some shenanigans to calculate predation of each predator over each prey cohort:
    simulation_i <- simulation %>%
      filter(length_fresh < predator_cutoff) %>% 
      slice(rep(row_number(), nrow(sim_predator_i))) %>%
      arrange(length_fresh) %>%
      
      mutate(prey_length = length_fresh / 1000,
             prey_density = mean_density,
             predator_length = rep(sim_predator_i$length_fresh / 1000, n() / nrow(sim_predator_i)),
             predator_density = rep(sim_predator_i$mean_density, n() / nrow(sim_predator_i)),
             prey_dry_weight = 0.0008 * exp(prey_length * 1000 *  0.9052),   # Reglero et al 2018
             predator_dry_weight = 0.0008 * exp(predator_length * 1000 *  0.9052),   # Reglero et al 2018
             
             predator_swimming_velocity = 3 * predator_length * 3600,
             prey_swimming_velocity = 3 * prey_length * 3600,
             msa = 4.699 * ((predator_length * 1000.)^ -1.129),  # Hilder et al (2019) 
             visual_r = prey_shape_proportion * 0.5 * prey_length / tan(0.5 * msa * pi/180), # Job & Bellwood 1996 JFB 48:952-963, converted from degrees to radians (1 degree =  pi/180 radians)
             visual_r = behav_anatomical_ratio * visual_r,
             length_proportion = (prey_length - (prey_cutoff / 1000)) / (predator_length - ((prey_cutoff) / 1000)),
             probability_capture_success = case_when(prey_length <= (prey_cutoff / 1000) ~ 1,
                                                     T ~ (1 - length_proportion)^k),
             encounter_rate_h = pi * (visual_r^2) * used_visual_area * predator_density * 
               (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen & Strickler 1977
             survival_h = exp(-encounter_rate_h * probability_capture_success),
             survival_d = exp(-encounter_rate_h * probability_capture_success * light)) %>%
      
      as_data_frame() %>%
      dplyr::select(-length_fresh, -mean_density)
    
    simulation2 <- simulation_i %>%
      dplyr::select(prey_length, predator_length, prey_density, survival_d) %>%
      group_by(prey_length, prey_density) %>%
      summarise(total_survival_d = prod(survival_d)) %>%
      mutate(prey_density_24 = prey_density * total_survival_d)
    
    # What was the tota survival rate?
    total_predator_density <- sum(unique(simulation_i$predator_density))
    pressure_group <- i
    prey_t0 <- sum(sim_prey$mean_density)
    prey_t24 <- sum(simulation2$prey_density_24)
    probability_survival <- prey_t24 / prey_t0
    df1 <- rbind(df1, data.frame(pressure_group, total_predator_density, prey_t0, prey_t24, probability_survival))
    
    }
  
  # Make a plot:
  p <- ggplot(data = simulation_i)  +
    
    annotation_custom(rasterGrob(img_ysl, interpolate = TRUE), xmin = 2, xmax = 3.9, ymin = .18, ymax = .22) +
    annotation_custom(rasterGrob(img_preflex, interpolate = TRUE), xmin = 4.8, xmax = 7.4, ymin = .17, ymax = .23) +
    
    geom_vline(xintercept = c(prey_cutoff), linetype = "dotted", alpha = .5) +
    geom_col(aes(x = prey_length * 1000, y = prey_density / nrow(sim_predator_i)), width = .09) +
    lims(x = c(2.2, 7.5)) + # 13.1
    labs(x = "Standard length (mm)", y = expression("#")) +         # ALB density (larvae m"^-3*")
    geom_col(data = simulation2, aes(x = prey_length * 1000, y = prey_density_24), fill = "black", width = .09) +
    geom_line(data = simulation2, aes(x = prey_length * 1000, y = total_survival_d * .19/1), show.legend = FALSE, alpha = .9, color = "red") +
    scale_y_continuous(sec.axis = sec_axis(~ .* (1/.19), name = "#"), limits = c(0, .2)) +      # "ALB probability to survie 1 day"
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title.y.right = element_text(color = "red"),
          axis.text.y.right = element_text(color = "red"))
  
  
  png(filename = paste0("plots/functional_response_simulated_station_pressure_", j, ".png"), width = 5, height = 5, units = "cm", res = 700)
  print(p)
  dev.off()
  
  # Save the last data frame to prevent overwriting:
  df1_sim <- df1
  
}


##############################################################


##################################################
# Run it for the entire data set:

# Create output dataframes:
df3 <- data.frame()

# Run the forloop: (This took some 20 minuts to run wirh an i7-8750H CPU @ 2.20GHz - go for coffee!!)
time0 <- Sys.time()
for(j in id_operations) {
  sim_prey <- data %>%
    as.data.frame() %>%
    filter(id == j, species == "Thunnus alalunga", length_fresh < predator_cutoff) %>%
    dplyr::select(length_fresh, mean_density) 
  
  
  # Prepare an ouput table:
  df2 <- data.frame()
  
  
  for (i in unique(sim_predator$pressure)) {
    
    sim_predator_i <- sim_predator %>%
      filter(pressure == i) %>%
      dplyr::select(-pressure)
    
    
    
    # Combine both data frames:
    simulation <- rbind(sim_prey, sim_predator_i)
    
    
    # Make some shenanigans to calculate predation of each predator over each prey cohort:
    simulation_i <- simulation %>%
      filter(length_fresh < predator_cutoff) %>% 
      slice(rep(row_number(), nrow(sim_predator_i))) %>%
      arrange(length_fresh) %>%
      
      mutate(prey_length = length_fresh / 1000,
             prey_density = mean_density,
             predator_length = rep(sim_predator_i$length_fresh / 1000, n() / nrow(sim_predator_i)),
             predator_density = rep(sim_predator_i$mean_density, n() / nrow(sim_predator_i)),
             prey_dry_weight = 0.0008 * exp(prey_length * 1000 *  0.9052),   # Reglero et al (2018)
             predator_dry_weight = 0.0008 * exp(predator_length * 1000 *  0.9052),   # Reglero et al (2018)
             
             predator_swimming_velocity = 3 * predator_length * 3600,
             prey_swimming_velocity = 3 * prey_length * 3600,
             msa = 4.699 * ((predator_length * 1000.)^ -1.129),  # Hilder et al (2019) 
             visual_r = prey_shape_proportion * 0.5 * prey_length / tan(0.5 * msa * pi/180),# Job & Bellwood 1996 JFB 48:952-963, converted from degrees to radians (1 degree =  pi/180 radians)
             visual_r = behav_anatomical_ratio * visual_r,
             length_proportion = (prey_length - (prey_cutoff / 1000)) / (predator_length - ((prey_cutoff) / 1000)),
             probability_capture_success = case_when(prey_length <= (prey_cutoff / 1000) ~ 1,
                                                     T ~ (1 - length_proportion)^k),
             encounter_rate_h = pi * (visual_r^2) * used_visual_area * predator_density * 
               (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen & Strickler 1977
             survival_h = exp(-encounter_rate_h * probability_capture_success),
             survival_d = exp(-encounter_rate_h * probability_capture_success * light)) %>% 
      
      as_data_frame() %>%
      dplyr::select(-length_fresh, -mean_density)
    
    simulation2 <- simulation_i %>%
      dplyr::select(prey_length, predator_length, prey_density, survival_d) %>%
      group_by(prey_length, prey_density) %>%
      summarise(total_survival_d = prod(survival_d)) %>%
      mutate(prey_density_24 = prey_density * total_survival_d)
    
    # What was the tota survival rate?
    total_predator_density <- sum(unique(simulation_i$predator_density))
    pressure_group <- i
    id <- j
    prey_t0 <- sum(sim_prey$mean_density)
    prey_t24 <- sum(simulation2$prey_density_24)
    probability_survival <- prey_t24 / prey_t0
    df2 <- rbind(df2, data.frame(id, pressure_group, total_predator_density, prey_t0, prey_t24, probability_survival))
    
  }
  
  df3 <- rbind(df3, df2)
  
}
Sys.time() -time0


# Calculate expected mean survival of all observed stations at increasing piscivorous pressure:
# Make id a factor:
df3$id <- as.factor(df3$id)

# Now calculate mean survival +- SE:
df3_mean <- df3 %>%
  mutate(id = as.factor(id)) %>%
  filter(!is.na(probability_survival)) %>%
  group_by(total_predator_density) %>%
  summarise(mean_survival = mean(probability_survival),
            sd_survival = sd(probability_survival),
            sd_up_surrvival = mean_survival + sd_survival,
            sd_low_surrvival = mean_survival - sd_survival,
            error_survival = sd_survival / sqrt(n()),
            error_up_surrvival = mean_survival + error_survival,
            error_low_surrvival = mean_survival - error_survival)

# What is the mean survival at the greatest predator pressure (0.168 BFT/m^3):
min(df3_mean$mean_survival)


# Pick the 4 rows of the simulated data representing predator densities = 0.002, 0.02, 0.080 and 0.168 BFT/m^3:
df12_sim <- df1_sim[c(2, 10, 37, 76),]  


# Before plotting, let's arrange the predator data set to include it in the plot as a % observations:
predator_size_distribution <- predator_size_distribution %>%
  group_by(id) %>%
  summarise(total_mean_density = sum(mean_density)) %>% 
  group_by(density_range = cut(total_mean_density, breaks = seq(0.000001, 0.170001, by = 0.005))) %>% 
  summarise(n = n()) %>%
  mutate(percent = 100 * n / length(id_operations))

# Make a vector with the mid value of each density range:
my_density <- c(.0025, .0075, .0125, .0175, .0225, .0325, .0375, .0425, .0475, .0525, .0625, .0775, .1025, .1325, .1675)

# Put it in a data frame:
predator_size_distribution <- data.frame(total_mean_density = my_density, percent_observations = predator_size_distribution$percent)
  

# Make the plot:
p <- ggplot() + 
  
  geom_col(data = predator_size_distribution, aes(x = total_mean_density, y = percent_observations),  fill = "darkorange1", color = "white")  + 
  
  geom_line(data = df3, aes(x = total_predator_density, y = 5.5 * probability_survival, 
                            group = as.factor(id)), colour = "gray28", show.legend = FALSE, alpha = .15, size = .1) +
  geom_line(data = df3_mean, aes(x = total_predator_density, y =  5.5 * mean_survival))  +
  geom_line(data = df1_sim, aes(x = total_predator_density, y =  5.5 * probability_survival),
            color = "black", linetype = "dashed", show.legend = FALSE, alpha = 1, size = .3) +
  geom_point(data = df12_sim, aes(x = total_predator_density, y =  5.5 * probability_survival), colour = "black", shape = 22, fill = "white", size = 2) +
  
  scale_y_continuous(sec.axis = sec_axis(~ (.) * (100/5.5), name = "% ALB surviving one day")) +
  
  labs(x = expression("BFT density (larvae m"^-3*")"), y = expression("% observed stations at given BFT density")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y.left = element_text(colour = "darkorange1"),  # #00A5FF darkorange1 black
        axis.text.y.left = element_text(color = "darkorange1"),
        axis.title.x.bottom = element_text(colour = "darkorange1"),
        axis.text.x.bottom  = element_text(color = "darkorange1"))

png(filename = "plots/functional_response_probability_surviving_all_stations.png", width = 6, height = 10, units = "cm", res = 400)
p
dev.off()


#                 END OF SCRIPT
############################################################