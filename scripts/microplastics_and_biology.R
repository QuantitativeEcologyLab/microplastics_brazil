# This script processes runs the analyses investigating the relationship
# between the interspecific variation in microplastic abundances,
# concentrations, and polymer types. 

# Written by Michael Noonan

#Load in any necessary packages
library(mgcv)

#Import the MP datasets
source("scripts/data_import.R")


#-------------------------------------------------------------
# Relationship with sex in tapirs and anteaters
#-------------------------------------------------------------

#Test for a relationship between sex and blood MP concentration in tapirs
fit <- gam(mp_ml ~ sex,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Tapirus_terrestris",],
           method = "REML")

summary(fit)


#Test for a relationship between age and blood MP concentration in anteaters
fit <- gam(mp_ml ~ sex,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Myrmecophaga_tridactyla",],
           method = "REML")

summary(fit)

#-------------------------------------------------------------
# Plastic concentrations vs. age in tapirs
#-------------------------------------------------------------


#Test for a correlation between age and blood MP concentration
fit <- gam(mp_ml ~ age,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Tapirus_terrestris",],
           method = "REML")

summary(fit)



#--------------------------------------------------------------------
# Plastic concentrations vs. body weight
#--------------------------------------------------------------------

#MPs vs. body weight in all three species
fit <- gam(mp_ml ~ weight + s(species, bs = 're'),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)

#Tapirs body weight in males vs. females
fit <- gam(mp_ml ~ weight + s(sex, weight, bs = 're'),
           family = tw(link = "log"),
           data = mp_data[which(mp_data$species == "Tapirus_terrestris"),],
           method = "REML")

summary(fit)

#Giant anteater body weight in males vs. females
fit <- gam(mp_ml ~ weight + s(sex, weight, bs = 're'),
           family = tw(link = "log"),
           data = mp_data[which(mp_data$species == "Myrmecophaga_tridactyla"),],
           method = "REML")

summary(fit)



#Giant armadillo body weight (only have data on females)
fit <- gam(mp_ml ~ weight,
           family = tw(link = "log"),
           data = mp_data[which(mp_data$species == "Priodontes_maximus"),],
           method = "REML")

summary(fit)


#-------------------------------------------------------------
# Plastic concentrations vs. body mass
#-------------------------------------------------------------

#Test for a correlation between body mass and blood MP concentration across all species
fit <- gam(mp_ml ~ weight + s(species, bs = "re") + s(weight, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)

#Test for a correlation between body mass and blood MP concentration in tapirs
fit <- gam(mp_ml ~ weight,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Tapirus_terrestris",],
           method = "REML")

summary(fit)


#Test for a correlation between body mass and blood MP concentration in anteaters
fit <- gam(mp_ml ~ weight,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Myrmecophaga_tridactyla",],
           method = "REML")

summary(fit)


#Test for a correlation between body mass and blood MP concentration in armadillos
fit <- gam(mp_ml ~ weight,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Priodontes_maximus",],
           method = "REML")

summary(fit)
