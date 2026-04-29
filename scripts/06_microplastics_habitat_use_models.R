# This script processes runs the analyses investigating the relationship
# between the observed concentrations of microplastics in blood
# and the animals' patterns of movement and habitat use.

# Written by Michael Noonan


#Load in any necessary packages
library(mgcv)

#Import the MP datasets
source("scripts/04_microplastics_data_import.R")

#Drop any individuals without movement data
mp_data <- mp_data[!is.na(mp_data$hr),]

#-------------------------------------------------------------
# Relationship with home-range size
#-------------------------------------------------------------

#Test for a correlation between home range size and blood MP concentration
fit <- gam(mp_ml ~ hr + s(species, bs = "re") + s(hr, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


#-------------------------------------------------------------
# Relationship with diffusion rate
#-------------------------------------------------------------

#Test for a correlation between home range size and blood MP concentration
fit <- gam(mp_ml ~ diffusion + s(species, bs = "re") + s(diffusion, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


#-------------------------------------------------------------
# Correlation with mean distance to development
#-------------------------------------------------------------

#Test for a correlation between mean distance to development and blood MP concentration
fit <- gam(mp_ml ~ mean_dist_development + s(species, bs = "re") + s(mean_dist_development, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)

#Test for a correlation between mean distance to development and blood MP concentration
fit <- gam(mp_ml ~ mean_dist_development + s(species, bs = "re") + s(mean_dist_development, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data[mp_data$species != "Priodontes_maximus",],
           method = "REML")

summary(fit)

#Test for a correlation between mean distance to development and blood MP concentration
fit <- gam(mp_ml ~ mean_dist_development,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Tapirus_terrestris",],
           method = "REML")

summary(fit)

#-------------------------------------------------------------
# Correlation with mean distance to waste treatment facilities
#-------------------------------------------------------------


#Test for a correlation between mean distance to development and blood MP concentration
fit <- gam(mp_ml ~ mean_dist_waste + s(species, bs = "re") + s(mean_dist_waste, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


#-------------------------------------------------------------
# Correlation with max human footprint index
#-------------------------------------------------------------

#Test for a correlation between max HFI and blood MP concentration
fit <- gam(mp_ml ~ max_HFI + s(species, bs = "re") + s(max_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


#Predict MP concentrations at max HFI
PREDS <- predict(fit,
                 newdata = (data.frame(max_HFI = 0, species = "null")),
                 type = "link",
                 exclude = "s(species)",
                 se.fit = TRUE)

#mean
exp(PREDS$fit)

#min
exp(PREDS$fit - 1.96*PREDS$se.fit)

#max
exp(PREDS$fit + 1.96*PREDS$se.fit)


#Predict MP concentrations at max HFI
PREDS <- predict(fit,
                 newdata = (data.frame(mean_HFI = max(mp_data$max_HFI), species = "null")),
                 type = "link",
                 exclude = "s(species)",
                 se.fit = TRUE)

#mean
exp(PREDS$fit)

#min
exp(PREDS$fit - 1.96*PREDS$se.fit)

#max
exp(PREDS$fit + 1.96*PREDS$se.fit)


#-------------------------------------------------------------
# Correlation with mean human footprint index
#-------------------------------------------------------------

#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ mean_HFI + s(species, bs = "re") + s(mean_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


#-------------------------------------------------------------
# Correlation with native forests
#-------------------------------------------------------------


#Test for a correlation between native forests and blood MP concentration
fit <- gam(mp_ml ~ Native_forest + s(species, bs = "re") + s(Native_forest, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


#-------------------------------------------------------------
# Correlation with agricultural land
#-------------------------------------------------------------


#Test for a correlation between the amount of agricultural land in the HR and blood MP concentration
fit <- gam(mp_ml ~ Agriculture + s(species, bs = "re") + s(Agriculture, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


#-------------------------------------------------------------
# Correlation with water and wetlands
#-------------------------------------------------------------


#Test for a correlation between water and wetlands and blood MP concentration
fit <- gam(mp_ml ~ Water + s(species, bs = "re") + s(Water, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)





