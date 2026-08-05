# This script processes runs the analyses investigating the relationship
# between the interspecific variation in microplastic abundances,
# concentrations, and polymer types. 

# Written by Michael Noonan

#Load in any necessary packages
library(mgcv)

#Import the MP datasets
source("scripts/04_microplastics_data_import.R")

#Exclude the giant armadillos for the main text analyses
mp_data <- mp_data[mp_data$species != "Priodontes_maximus",]

#-------------------------------------------------------------
# Differences in microplastic concentrations between species
#-------------------------------------------------------------

#Identify individuals without MPs in their blood
mp_data[which(mp_data$mp_ml == 0),]

#Average concentration
mean(mp_data$mp_ml)

#Range
range(mp_data$mp_ml)

#Estimate the mean concentration and 95%CIs in lowland tapirs
tapir_mean <- gam(mp_ml ~ 1,
                  family = tw(link = "log"),
                  data = mp_data[mp_data$species == "Tapirus_terrestris",])
preds <- predict(tapir_mean, se.fit = TRUE)
exp(preds$fit[1])
exp(preds$fit[1] - 1.96*preds$se.fit[1])
exp(preds$fit[1] + 1.96*preds$se.fit[1])


#Estimate the mean concentration and 95%CIs in giant anteaters
anteater_mean <- gam(mp_ml ~ 1,
                     family = tw(link = "log"),
                     data = mp_data[mp_data$species == "Myrmecophaga_tridactyla",])
preds <- predict(anteater_mean, se.fit = TRUE)
exp(preds$fit[1])
exp(preds$fit[1] - 1.96*preds$se.fit[1])
exp(preds$fit[1] + 1.96*preds$se.fit[1])


#Estimate the mean concentration and 95%CIs in giant armadillos
armadillo_mean <- gam(mp_ml ~ 1,
                      family = tw(link = "log"),
                      data = mp_data[mp_data$species == "Priodontes_maximus",])
preds <- predict(armadillo_mean, se.fit = TRUE)
exp(preds$fit[1])
exp(preds$fit[1] - 1.96*preds$se.fit[1])
exp(preds$fit[1] + 1.96*preds$se.fit[1])


#Test for any differences between giant anteaters and tapirs
#Sample size for the giant armadillos is too small for inter-specific comparison
fit <- gam(mp_ml ~ species,
           family = tw(link = "log"),
           data = mp_data[mp_data$species != "Priodontes_maximus",])

summary(fit)


#---------------------------------------------------------------------
# Differences in polymer concentrations between species
#---------------------------------------------------------------------

#Generate the pairwise combinations of species excluding giant armadillos
all_pairs <- combn(unique(mp_data[mp_data$species != "Priodontes_maximus","species"]), 2, simplify = F)

mp_data$biome <- as.factor(mp_data$biome)

#Run the analyses
results <- list()
for(i in 1:length(all_pairs)){
  
  #Subset the dataset to the ith pair of species
  data_sub <-  mp_data[mp_data$species %in% all_pairs[[i]],]
  
  #Empty list for storing results
  res <- list()
  
  #Loop over the vector of polymer names
  for(j in 1:length(polymers)){
    
    
    #Fit the GLM for the jth polymer
    fit <- gam(formula(paste(polymers[j]," ~ species + s(biome, bs = 're')")),
               family = tw(link = "log"),
               data = data_sub,
               method = "REML")
    
    #Assemble the results into a data frame
    res[[j]] <- data.frame(pair = paste(all_pairs[[i]][1], all_pairs[[i]][2]),
                           species1 = all_pairs[[i]][1],
                           species2 = all_pairs[[i]][2],
                           polymer = polymers[j],
                           coef_spp = sub("species","",names(summary(fit)$p.coeff[2])),
                           beta = unname(summary(fit)$p.coeff[2]),
                           t = unname(summary(fit)$p.t[2]),
                           p = unname(summary(fit)$p.pv[2]))
    
  } #closes the loop over the polymer names
  
  results[[i]] <- do.call(rbind,res)
  rm(res)
  
}#closes the loop over the species pairs

#Convert from a list to a data frame
results <- do.call(rbind,results)

#Adjust the p valuesfor multiple comparisons using Benjamini & Hochberg's correction
results$p.adjusted <-  p.adjust(results$p, method = "BH")

#Which were significant
results[results$p.adjusted < 0.05,]




#---------------------------------------------------------------------
# Size distribution across species
#---------------------------------------------------------------------

#Average size
mean(sizes$Length)

#Standard deviation
sd(sizes$Length)

#Range
range(sizes$Length)
