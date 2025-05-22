#Load in any necessary packages
library(FactoMineR)
library(ggplot2)
library(ellipse)
library(randomForest)
library(caret)
library(ggridges)
library(gridExtra)
library(mgcv)
library(vegan)
library(fmsb)

#Import the MP datasets
source("scripts/data_import.R")

#---------------------------------------------------------------------
# Differences in polymer concentrations
#---------------------------------------------------------------------

#Generate the pairwise combinations
PAIRS <- combn(unique(mp_data$species),2, simplify = F)

#Run the analyses
results <- list()
for(i in 1:length(PAIRS)){
  
  #Subset the dataset
  data_sub <-  mp_data[mp_data$species %in% PAIRS[[i]],]
  
  res <- list()
  
  #Loop over the polymers
  for(j in 44:56){
    
    formula(paste(names(data_sub)[j]," ~ species"))
    
    fit <- gam(formula(paste(names(data_sub)[j]," ~ species")),
               family = tw(link = "log"),
               data = data_sub,
               method = "REML")
    
    res[[j-12]] <- data.frame(pair = paste(PAIRS[[i]][1], PAIRS[[i]][2]),
                              species1 = PAIRS[[i]][1],
                              species2 = PAIRS[[i]][2],
                              polymer = names(data_sub)[j],
                              coef_spp = sub("species","",names(summary(fit)$p.coeff[2])),
                              beta = unname(summary(fit)$p.coeff[2]),
                              t = unname(summary(fit)$p.t[2]),
                              p = unname(summary(fit)$p.pv[2]))
    
  }
  
  results[[i]] <- do.call(rbind,res)
  
}

results <- do.call(rbind,results)

#Adjust for multiple comparisons using Benjamini & Hochberg's correction
results$p.adjusted <-  p.adjust(results$p, method = "BH")
results[which(p.adjust(results$p, method = "BH") < 0.05),]
