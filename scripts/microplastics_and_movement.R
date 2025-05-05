# This script processes runs the analyses investigating the relationship
# between the observed concentrations of microplastics in blood
# and the animals' patterns of movement and habitat use.

# Written by Michael Noonan


#Load in any necessary packages
library(mgcv)

#Import the MP datasets
source("scripts/data_import.R")

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
# Correlation with mean human footprint index
#-------------------------------------------------------------

#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ mean_HFI + s(species, bs = "re") + s(mean_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)

#Predict MP concentrations at max HFI
PREDS <- predict(fit,
                 newdata = (data.frame(mean_HFI = 0, species = "null")),
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

#Final steps are to run a series of sensitivity analysis on the HFI effect
# due to the unblanced distribution of HFI values

#First step is to remove any data points above the 95th quantile of the mean HFI exposure
mp_data_subset <- mp_data[mp_data$mean_HFI < quantile(mp_data$mean_HFI, .95),]

#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ mean_HFI + s(species, bs = "re") + s(mean_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data_subset,
           method = "REML")

summary(fit)

#Next, remove any data points above the 90th quantile of the mean HFI exposure
mp_data_subset <- mp_data[mp_data$mean_HFI < quantile(mp_data$mean_HFI, .90),]

#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ mean_HFI + s(species, bs = "re") + s(mean_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data_subset,
           method = "REML")

summary(fit)

#Note how in both instances, the effect size appears stronger,
# though the direction and significance are unaffected.


#-------------------------------------------------------------
# Correlation with max human footprint index
#-------------------------------------------------------------

#Test for a correlation between max HFI and blood MP concentration
fit <- gam(mp_ml ~ max_HFI + s(species, bs = "re") + s(max_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


#Final steps are to run a series of sensitivity analysis on the HFI effect

#First step is to remove any data points above the 95th quantile of the max HFI exposure
mp_data_subset <- mp_data[mp_data$max_HFI < quantile(mp_data$max_HFI, .95),]

#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ max_HFI + s(species, bs = "re") + s(max_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data_subset,
           method = "REML")

summary(fit)

#Next, remove any data points above the 95th quantile of the max HFI exposure
mp_data_subset <- mp_data[mp_data$max_HFI < quantile(mp_data$max_HFI, .90),]

#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ max_HFI + s(species, bs = "re") + s(max_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data_subset,
           method = "REML")

summary(fit)


#Note how again, the direction and significance are unaffected
# by removing the upper quantile of maximum HFI exposure.


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




