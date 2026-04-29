# This script processes runs a number of sensitivity analyses
# on the results presented in the main text.
# Where relevant, comments direct to the core analyses for which sensitivity is
# being assessed.

# Written by Michael Noonan


#Load in any necessary packages
library(mgcv)

#Import the MP datasets
source("scripts/04_microplastics_data_import.R")

#-------------------------------------------------------------
# Sensitivity analyses on the species imbalance
#-------------------------------------------------------------

#The following steps assess the sensitivity of the core findings to the
#species imbalance by sequentially dropping a species and refitting the models

#---------------------
# Maximum HFI
#---------------------

#Test for a correlation between maximum HFI exposure and blood MP concentration
fit <- gam(mp_ml ~ max_HFI + s(species, bs = "re") + s(max_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)

#Excluding giant armadillos
fit <- gam(mp_ml ~ max_HFI + s(species, bs = "re") + s(max_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$species != "Priodontes_maximus",],
           method = "REML")

summary(fit)

#Tapir only
fit <- gam(mp_ml ~ max_HFI,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Tapirus_terrestris",],
           method = "REML")

summary(fit)

#---------------------
# Distance to development
#---------------------

#Test for a correlation between mean distance to development and blood MP concentration
fit <- gam(mp_ml ~ mean_dist_development + s(species, bs = "re") + s(mean_dist_development, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)

#Excluding giant armadillos
fit <- gam(mp_ml ~ mean_dist_development + s(species, bs = "re") + s(mean_dist_development, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$species != "Priodontes_maximus",],
           method = "REML")

summary(fit)

#Tapir only
fit <- gam(mp_ml ~ mean_dist_development,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Tapirus_terrestris",],
           method = "REML")

summary(fit)


#---------------------
# Distance to nearest waste treatment facility
#---------------------

#Test for a correlation between mean distance to waste treatment facility and blood MP concentration
fit <- gam(mp_ml ~ mean_dist_waste + s(species, bs = "re") + s(mean_dist_waste, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)

#Excluding giant armadillos
fit <- gam(mp_ml ~ mean_dist_waste + s(species, bs = "re") + s(mean_dist_waste, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$species != "Priodontes_maximus",],
           method = "REML")

summary(fit)

#Tapir only
fit <- gam(mp_ml ~ mean_dist_waste,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Tapirus_terrestris",],
           method = "REML")

summary(fit)



#-------------------------------------------------------------
# Sensitivity analyses on the disturbance gradient imbalance
#-------------------------------------------------------------

#---------------------
# Maximum HFI
#---------------------

#Test for a correlation between max HFI and blood MP concentration
#only for data points below the 95th quantile of the max HFI exposure
fit <- gam(mp_ml ~ max_HFI + s(species, bs = "re") + s(max_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$max_HFI < quantile(mp_data$max_HFI, .95),],
           method = "REML")

summary(fit)

#Next, repeat for the 90th quantile of the mean HFI exposure

#Test for a correlation between max HFI and blood MP concentration
#only for data points below the 90th quantile of the max HFI exposure
fit <- gam(mp_ml ~ max_HFI + s(species, bs = "re") + s(max_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$max_HFI < quantile(mp_data$max_HFI, .90),],
           method = "REML")

summary(fit)


#---------------------
# Distance to development
#---------------------

#Test for a correlation between distance to development and blood MP concentration
#only for data points below the 95th quantile of the distance to development
fit <- gam(mp_ml ~ mean_dist_development + s(species, bs = "re") + s(mean_dist_development, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$mean_dist_development < quantile(mp_data$mean_dist_development, .95),],
           method = "REML")

summary(fit)

#Next, repeat for the 90th quantile

#Test for a correlation between distance to development and blood MP concentration
#only for data points below the 90th quantile of the distance to development
fit <- gam(mp_ml ~ mean_dist_development + s(species, bs = "re") + s(mean_dist_development, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$mean_dist_development < quantile(mp_data$mean_dist_development, .90),],
           method = "REML")

summary(fit)


#---------------------
# Distance to waste treamt
#---------------------

#Test for a correlation between distance to development and blood MP concentration
#only for data points below the 95th quantile of the distance to development
fit <- gam(mp_ml ~ mean_dist_development + s(species, bs = "re") + s(mean_dist_development, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$mean_dist_development < quantile(mp_data$mean_dist_development, .95),],
           method = "REML")

summary(fit)

#Next, repeat for the 90th quantile

#Test for a correlation between distance to development and blood MP concentration
#only for data points below the 90th quantile of the distance to development
fit <- gam(mp_ml ~ mean_dist_development + s(species, bs = "re") + s(mean_dist_development, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$mean_dist_development < quantile(mp_data$mean_dist_development, .90),],
           method = "REML")

summary(fit)