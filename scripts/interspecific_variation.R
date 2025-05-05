# This script processes runs the analyses investigating the relationship
# between the interspecific variation in microplastic abundances,
# concentrations, and polymer types. 

# Written by Michael Noonan

#Load in any necessary packages
library(mgcv)

#Import the MP datasets
source("scripts/data_import.R")

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


#Test for any differences across species using a likelihood ratio test
fit <- gam(mp_ml ~ species,
           family = tw(link = "log"),
           data = mp_data)

fit_null <- gam(mp_ml ~ 1,
                family = tw(link = "log"),
                data = mp_data)

anova(fit_null,fit, test = "Chisq")


#---------------------------------------------------------------------
# Differences in polymer concentrations between species
#---------------------------------------------------------------------

#Generate the pairwise combinations of species
all_pairs <- combn(unique(mp_data$species), 2, simplify = F)

#Run the analyses
results <- list()
for(i in 1:length(all_pairs)){
  
  #Subset the dataset to the ith pair of species
  data_sub <-  mp_data[mp_data$species %in% all_pairs[[i]],]
  
  #Empty list for storing results
  res <- list()
  
  #Loop over the vector of polymer names
  for(j in 1:length(polymer_names)){
    
    
    #Fit the GLM for the jth polymer
    fit <- gam(formula(paste(polymer_names[j]," ~ species")),
               family = tw(link = "log"),
               data = data_sub,
               method = "REML")
    
    #Assemble the results into a data frame
    res[[j]] <- data.frame(pair = paste(all_pairs[[i]][1], all_pairs[[i]][2]),
                              species1 = all_pairs[[i]][1],
                              species2 = all_pairs[[i]][2],
                              polymer = polymer_names[j],
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
results[which(p.adjust(results$p, method = "BH") < 0.05),]




#---------------------------------------------------------------------
# Size distribution across species
#---------------------------------------------------------------------

#Average size
mean(sizes$Length)

#Standard deviation
sd(sizes$Length)

#Range
range(sizes$Length)
