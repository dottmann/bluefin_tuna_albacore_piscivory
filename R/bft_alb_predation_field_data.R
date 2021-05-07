
#   Temporal trends of predatory pressure in the field

# Author: Daniel Ottmann
# Created: May 2020
# Last update: May 2021


###########################################
#       Readme

# This script evaluates predation of larval bluefin tuna on larval albacore in the western Mediterranean Sea
# We compare this with changes in the abundances of albacore and bluefin tuna


##################################################################

##################################################################


#############################
# Load packages:
library(tidyverse)


##################################################################
# LOAD  DATA

load("data/data.RData")
load("data/df_postflexion.RData")
load("data/data2.RData")


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


################################################################################

################################################################################


# Now let's get values for the entire data set:

# Prepare an ouput table:
df1 <- data.frame(row.names = c("id_operation", "id_netcolector", "prey_t0", "prey_t24", "prey_removed"))


for(i in unique(data$id_operation)) {
  
  # First prepare the denisty and length data of predaotrs and prey:
  observations_df <- data %>%
    filter(id_operation %in% i)
  
  # Make a data frame with prey and predator info:
  prey_df <- observations_df %>%
    filter(species == "Thunnus alalunga",length_fresh < predator_cutoff)
  
  predator_df <- observations_df %>%
    filter(species ==  "Thunnus thynnus", length_fresh >= predator_cutoff)
  
  # Combine both data frames:
  observations_df <- rbind(prey_df, predator_df)
  
  # Make some shenanigans to calculate predation of each predator cohort over each prey cohort, leaving intermediate sizes away:
  observations_df <- observations_df %>%
    filter(length_fresh < predator_cutoff) %>%                # Here we pick only the prey
    slice(rep(row_number(), nrow(predator_df))) %>%      # And replicate the prey vector as many times ar predators we have to allow each pairwise interaction
    mutate(prey_length = length_fresh / 1000,                 # Convert length to meters
           prey_density = mean_density,                       # Define is prey density values
           predator_length = rep(predator_df$length_fresh / 1000, n() / nrow(predator_df)),  # Replicate predator lengths (in meters) for each interaction
           predator_density = rep(predator_df$mean_density, n() / nrow(predator_df)),        # And assign their density
           
           # Now run the model and get the output:
           predator_swimming_velocity = 3 * predator_length * 3600 ,
           prey_swimming_velocity = 3 * prey_length * 3600,
           msa = 4.699 * ((predator_length * 1000.)^ -1.129),
           visual_r = prey_shape_proportion * 0.5 * prey_length / tan(0.5 * msa * pi/180),
           visual_r = behav_anatomical_ratio * visual_r,
           length_proportion = (prey_length - (prey_cutoff / 1000)) / (predator_length - ((prey_cutoff) / 1000)),
           probability_capture_success = case_when(prey_length <= (prey_cutoff / 1000) ~ 1,
                                                   T ~ (1 - length_proportion)^k),
           encounter_rate_h = pi * (visual_r^2) * used_visual_area * predator_density * 
             (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gertisen &Strickler 1977
           survival_h = exp(-encounter_rate_h * probability_capture_success),
           survival_d = exp(-encounter_rate_h * probability_capture_success * sum(light)),
           ingestion_rate_h = encounter_rate_h * probability_capture_success * prey_density)  %>% 
    
    # Format tibble to data frame and remove junk:
    as_data_frame() %>%
    dplyr::select(-id_operation, -length_fresh, -mean_density)
  
  # Get a table with predictied prey abundances after 1 day:
  observations_df2 <- observations_df %>%
    dplyr::select(prey_length, predator_length, prey_density, survival_d) %>%
    group_by(prey_length, prey_density) %>%
    summarise(total_survival_d = prod(survival_d)) %>%
    mutate(prey_density_24 = prey_density * total_survival_d)
  
  
  # What was the tota survival rate?
  prey_t0 <- sum(prey_df$mean_density)
  prey_t24 <- sum(observations_df2$prey_density_24)
  id_operation <- i
  df1 <- rbind(df1, data.frame(id_operation, prey_t0, prey_t24))
}

# The 0s result from stations where there were either no prey or no predators

# Fix that and get some values of survival and mortality:
df1 <- df1 %>%
  mutate(prey_t24 = case_when(prey_t24 == 0 ~ prey_t0,
                              T ~ prey_t24),
         prey_removed = prey_t0 - prey_t24,
         proportion_survivors = prey_t24 / prey_t0,
         proportion_killed = 1 - proportion_survivors)

head(df1)

# Some NA's belong to stations that had predators but no prey.
# To remove them keeping those samples that had predators but no prey, we do the follwoing:

# Keep track of stations with predator but no prey:
id_noprey_stations <- subset(df1, is.na(proportion_survivors))$id_operation

temp <- data2 %>%
  dplyr::select(id_operation, year, talalunga_n)

# Combine temp with df1:
df1 <- left_join(temp, df1, by = "id_operation")


# Filter real NA from prey_t0 and positive prey_t0:
temp1 <- subset(df1, is.na(prey_t0) & talalunga_n == 0 | prey_t0 > 0 & talalunga_n > 0)

# Filter  NA's that belong to sations that had predaotrs but not prey:
temp2 <- subset(df1, id_operation %in% id_noprey_stations)

# Combine the tables:
df1 <- rbind(temp1, temp2)


# What is the range of prey survival?
range(subset(df1, !is.na(proportion_survivors))$proportion_survivors)

# Fill NAs in prey_t0, prey_t24 and prey_removed with 0s:
df1 <- df1 %>%
  mutate(prey_t0 = case_when(is.na(prey_t0) ~ 0,
                             T ~ prey_t0),
         prey_t24 = case_when(is.na(prey_t24) ~ 0,
                             T ~ prey_t24),
         prey_removed = case_when(is.na(prey_removed) ~ 0,
                                  T ~ prey_removed))


# Clean environment:
rm(list = c("temp", "temp1", "temp2", "prey_df", "predator_df", "observations_df", "observations_df2", "i"))


######################################################################

######################################################################

# Get a data frame with the risk of dying for a 4.1 mm albacore:
data2_temp <- data2 %>%
  dplyr::select(id_operation, year)

df_risc <- data %>%
  filter(id_operation %in% data2_temp$id_operation, 
         species ==  "Thunnus thynnus", length_fresh >= predator_cutoff)

# Make some shenanigans to calculate predation of each predator cohort over each prey cohort, leaving intermediate sizes away:
df_risc <- df_risc %>%
  arrange(id_operation) %>%
  mutate(predator_length = length_fresh / 1000,
         predator_density = mean_density,
         prey_length = prey_cutoff / 1000,
         predator_swimming_velocity = 3 * predator_length * 3600,
         prey_swimming_velocity = 3 * prey_length * 3600,
         msa = 4.699 * ((predator_length * 1000.)^ -1.129),  # Hilder PE, Battaglene SC, Hart NS, Collin SP, Cobcroft JM (2019) 
         visual_r = prey_shape_proportion * 0.5 * prey_length / tan(0.5 * msa * pi/180),# Job & Bellwood 1996 JFB 48:952-963, converted from degrees to radians (1 degree =  pi/180 radians)
         visual_r = behav_anatomical_ratio * visual_r,
         length_proportion = (prey_length - (prey_cutoff / 1000)) / (predator_length - ((prey_cutoff) / 1000)),
         probability_capture_success = case_when(prey_length <= (prey_cutoff / 1000) ~ 1,
                                                 T ~ (1 - length_proportion)^k),
         encounter_rate_h = pi * (visual_r^2) * used_visual_area * predator_density * 
           (prey_swimming_velocity^2 + 3 * predator_swimming_velocity^2) / (3 * predator_swimming_velocity), # Gerritsen & Strickler 1977
         survival_h = exp(-encounter_rate_h * probability_capture_success),
         survival_d = exp(-encounter_rate_h * probability_capture_success * sum(light))) %>% 
  as_data_frame() %>%
  dplyr::select(-length_fresh, -mean_density) %>%
  group_by(id_operation) %>%
  summarise(total_survival_d = prod(survival_d))


# Next, we have to join both tables and replace survival == 1 where needed:
df_risc <- data2_temp %>% left_join(df_risc, by = "id_operation") %>%
  mutate(total_survival_d = case_when(is.na(total_survival_d) ~ 1,
                                      T ~ total_survival_d),
         probability_mortality_d = 1 - total_survival_d)


# Change in frequency of stations with 53% chance of dying in one day:

total_stations <- df_risc %>%
  mutate(period = case_when(year < 2012 ~ "Before BFT recovery",
                            T  ~ "After BFT recovery")) %>%
  group_by(year, period) %>%
  summarise(total_stations = n())

stations_50 <- df_risc %>%
  filter(probability_mortality_d > .53) %>%
  mutate(period = case_when(year < 2012 ~ "Before BFT recovery",
                            T  ~ "After BFT recovery")) %>%
  group_by(year, period) %>%
  summarise(stations_50 = n())

stations <- left_join(total_stations, stations_50, by = c("year", "period"))

stations %>%
  mutate(stations_50 = replace_na(stations_50, 0),
         perc_stations_50 = stations_50 / total_stations) %>%
  group_by(period) %>%
  summarise(mean_perc_stations_50 = mean(perc_stations_50))


# Clean environment:
rm(list = c("stations", "stations_50", "total_stations"))


#####################################################################

#####################################################################

#######################
# Plot temporal trends:

# Select the variables of interest of each table and put them in a single table:
data2_temp <- data2 %>%
  dplyr::select(id_operation, talalunga_dens)

df_risc <- df_risc %>%
  dplyr::select(-year)

df2 <- left_join(df1, df_risc, by = "id_operation")

df2 <- left_join(df2, data2_temp, by = "id_operation")

# Save this as df3 and aff period:
df3 <- df2 %>%
  mutate(period = case_when(year < 2012 ~ "Before tuna recovery",
                            T ~ "After tuna recovery"))

# Percent stations with presence of postflexion:
df3_presence <- df3 %>%
  mutate(postflexion_presence = case_when(probability_mortality_d == 0 ~ 0,
                                          T ~ 1)) %>%
  group_by(year, period) %>%
  summarise(percent_presence_postflex = 100 * sum(postflexion_presence) / n()) 


df2 <- df2 %>%
  group_by(year) %>%
  summarise(mean_probability_mortality_d = mean(probability_mortality_d),
         mean_prey_t0 = mean(prey_t0),
         sd_prey_t0 = sd(prey_t0),
         se_prey_t0 = sd_prey_t0 / sqrt(n()),
         mean_probability_surviving_d = mean(total_survival_d),
         mean_talalunga_dens = mean(talalunga_dens)) %>%
  mutate(period = case_when(year < 2012 ~ "Before tuna recovery",
                            T ~ "After tuna recovery"))

# With larval densities as bars:

p <- ggplot() +
  
  geom_point(data = df2, aes(x = year, y = mean_prey_t0), color = "#00A5FF", size = .5) +
  geom_line(data = df2, aes(x = year, y = mean_prey_t0), color = "#00A5FF", alpha = .1) +
  geom_line(data = df2, aes(x = year, y = mean_prey_t0, group = period), color = "#00A5FF") +
  geom_errorbar(data = df2, aes(x = year, y = mean_prey_t0, 
                                ymin = mean_prey_t0 - se_prey_t0, 
                                ymax = mean_prey_t0 + se_prey_t0), 
                width = .3, size = .6, color = "#00A5FF", alpha = .8) +
  
  geom_jitter(data = df3, aes(x = year, y = -0 + total_survival_d * 1/20, group = year), fill = "#FF6666", color = "#FF6666", alpha = .9, size = 1) +
  scale_y_continuous(sec.axis = sec_axis(~(. + 0) * (20), name = "Probability of surviving 1 day")) +
  
  geom_point(data = df_postflexion , aes(x = year, y = mean_predator_dens * 4), color = "darkorange1", size = .5) +
  geom_line(data = df_postflexion , aes(x = year, y = mean_predator_dens * 4, group = period), color = "darkorange1") +
  geom_line(data = df_postflexion , aes(x = year, y = mean_predator_dens * 4), color = "darkorange1", alpha = .1) +
  geom_errorbar(data = df_postflexion, aes(x = year, y = mean_predator_dens,
                                           ymin = ((mean_predator_dens - se_predator_dens) * 4),
                                           ymax = ((mean_predator_dens + se_predator_dens) * 4)),
                width = .3, size = .6, color = "darkorange1", alpha = .8) +
#  scale_fill_manual(values = c("#00A5FF", "darkorange1")) +
#  scale_y_continuous(sec.axis = sec_axis(~ (. -0) * (1/4), name = expression("Piscivorous BFT (larvae m"^-3*")"))) +
  
  geom_point(data = df3_presence, aes(x = year, y = percent_presence_postflex * 1/800), size = .6) +
  geom_line(data = df3_presence, aes(x = year, y = percent_presence_postflex * 1/800), alpha = .1) +
  geom_line(data = df3_presence, aes(x = year, y = percent_presence_postflex * 1/800, group = period)) +
#   scale_y_continuous(sec.axis = sec_axis(~(. + 0) * (800), name = "% stations with
# presence of piscivorous BFT")) +
  
  labs(x = "Year", y = expression("Mean albacore density (larvae m"^-3*")")) +
  
  
  theme_bw() +
  
  theme(panel.grid = element_blank(),
        axis.title.y.right = element_text(color = "#FF3333"), 
        axis.text.y.right = element_text(color = "#FF3333"),  
        axis.title.y.left = element_text(color = "#00A5FF"),
        axis.text.y.left = element_text(color = "#00A5FF"))

png(filename = "plots/probability_surviving_temporal_trends.png", width = 10, height = 8, units = "cm", res = 300)
p
dev.off()


# Do a Welch Two Sample t-test to copare if predator pressure before and after the BFT stock recovery differs:
t.test(subset(df3, period == "Before tuna recovery")$total_survival, subset(df3, period == "After tuna recovery")$total_survival)


#                 END OF PRESENTATION
###########################################################