# This script processes runs the analyses investigating the differences
# in microplastic concentrations between the different biomes

# Written by Michael Noonan

#Load in any necessary packages
library(mgcv)

#Import the MP datasets
source("scripts/data_import.R")


#--------------------------------------------------------------------
# Plastic concentrations in the different biomes
#--------------------------------------------------------------------

#MPs vs. body weight in all three species
fit <- gam(mp_ml ~ biome + s(species, bs = 're'),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)



#---------------------------------------------------------------------
# Differences in polymer concentrations between biomes
#---------------------------------------------------------------------

#Generate the pairwise combinations of biomes
all_pairs <- combn(unique(mp_data$biome), 2, simplify = F)

#Run the analyses
results <- list()
for(i in 1:length(all_pairs)){
  
  #Subset the dataset to the ith pair of biomes
  data_sub <-  mp_data[mp_data$biome %in% all_pairs[[i]],]
  
  #Empty list for storing results
  res <- list()
  
  #Loop over the vector of polymer names
  for(j in 1:length(polymer_names)){
    
    
    #Fit the GLM for the jth polymer
    fit <- gam(formula(paste(polymer_names[j]," ~ biome")),
               family = tw(link = "log"),
               data = data_sub,
               method = "REML")
    
    #Assemble the results into a data frame
    res[[j]] <- data.frame(pair = paste(all_pairs[[i]][1], all_pairs[[i]][2]),
                           biome1 = all_pairs[[i]][1],
                           biome2 = all_pairs[[i]][2],
                           polymer = polymer_names[j],
                           coef_biome = sub("species","",names(summary(fit)$p.coeff[2])),
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

