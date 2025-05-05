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

#Load in the data
metadata <- read.csv("data/mp_data/All_general.csv")
polymers <- read.csv("data/mp_data/All_PP.csv")
data_long <- merge(x = metadata, y = polymers, by.x = "sample", by.y = "Sample")

#Some data carpentry to get the format correct for analysis
polymers <- reshape(polymers, idvar = "Sample", timevar = "Microplastic", direction = "wide")
names(polymers)[-1] <- gsub('Adjusted.', '', names(polymers)[-1])
names(polymers)[-1] <- gsub(' ', '_', names(polymers)[-1])
names(polymers)[5] <- "POLYPROPYLENE"

data <- merge(x = metadata, y = polymers, by.x = "sample", by.y = "Sample")


#---------------------------------------------------------------------
# Differences in polymer concentrations
#---------------------------------------------------------------------

#Generate the pairwise combinations
PAIRS <- combn(unique(data$species),2, simplify = F)

#Run the analyses
results <- list()
for(i in 1:length(PAIRS)){
  
  #Subset the dataset
  data_sub <-  data[data$species %in% PAIRS[[i]],]
  
  res <- list()
  for(j in 13:25){
    
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




#---------------------------------------------------------------------
# Figure 2A - heatmap of polymer abundances
#---------------------------------------------------------------------


# First, reorder the 'Sample' column by its numeric component
ORDER <- metadata[order(metadata$species), "sample"]
data_long$sample <- factor(data_long$sample, levels = ORDER, ordered = TRUE)


#Adjust the polymer names
data_long$Microplastic <- stringr::str_to_title(data_long$Microplastic)
data_long$Microplastic[data_long$Microplastic == "Abs"] <- "ABS"
data_long$Microplastic[data_long$Microplastic == "Pet"] <- "PET"
data_long$Microplastic[data_long$Microplastic == "Polyvinyl Chloride"] <- "PVC"

# Heatmap
a <- 
ggplot(data_long, aes(Microplastic, sample, fill= log(Adjusted+1))) + 
  ggtitle("A") +
  geom_tile(alpha = 0.95) +
  scico::scale_fill_scico(palette = "lipari",
                          name = "Particles/mL",
                          breaks = c(0,log(10),log(50),log(400)),
                          labels = c(0,10,50,400)) +
  geom_hline(yintercept = 21.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = 26.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_text(size=4,
                                    family = "sans",
                                    angle = 90,
                                    face = "bold",
                                    color = "black",
                                    hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.position = "top",
        legend.title = element_text(size=5, family = "sans", face = "bold", vjust = -2, hjust = 0.5),
        legend.text = element_text(size=4, family = "sans", face = "bold", vjust = 4),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 10,
                               barheight = 0.3, direction = "horizontal"))


ggsave(a,
       width = 2.375, height = 5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_2_left.png")


#---------------------------------------------------------------------
# Figure 1C - E - Radar plots of polymer abundances
#---------------------------------------------------------------------


mean_polymers <- aggregate(Adjusted ~ species + Microplastic, FUN = "mean", data = data_long)
mean_polymers <- reshape(data = mean_polymers,
                         idvar = "species",
                         timevar = "Microplastic",
                         v.names = "Adjusted",
                         direction = "wide")
#Adjust the polymer names
names(mean_polymers)[-1] <- sub("Adjusted.","",names(mean_polymers)[-1])


#Define the min and max (specific requirement of the fmsb package)
mean_polymers[4,] <- 40 #Max
mean_polymers[5,] <- 0 #Min

row.names(mean_polymers) <- c("Myrmecophaga_tridactyla",
                              "Priodontes_maximus",
                              "Tapirus_terrestris",
                              "Max",
                              "Min")

mean_polymers <- mean_polymers[c(4:5,1:3),]

# Generate and save the figures
png(filename = "figures/figure_2_right.png",
    width = 2.375, height = 5, units = "in",
    bg = "transparent",
    res = 600)

par(mfrow = c(3,1),
    mar = c(0.2,0.1,0.2,0.1),
    font.axis = 2,
    font = 2)
# Create the radar charts for the anteaters
radarchart(mean_polymers[c(1:2,3),-1],
           axistype = 1,
           pty = NA,
           pcol = "#619b8a",
           pfcol = adjustcolor("#619b8a", alpha.f = 0.3),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.5,
           cglwd = 0.4,
           vlcex = 0.48)
title(main = "B",
      cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

# Create the radar chart for the armadillos
radarchart(mean_polymers[c(1:2,4),-1],
           axistype = 1,
           pty = NA,
           pcol = "#bb3e03",
           pfcol = adjustcolor("#bb3e03", alpha.f = 0.3),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.5,
           cglwd = 0.4,
           vlcex = 0.48)
title(main = "C", cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

# Create the radar chart for the tapirs
radarchart(mean_polymers[c(1:2,5),-1],
           axistype = 1,
           pty = NA,
           pcol = "#005f73",
           pfcol = adjustcolor("#005f73", alpha.f = 0.3),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.5,
           cglwd = 0.4,
           vlcex = 0.48)
title(main = "D", cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

dev.off()




dist_matrix <- vegdist(data[,13:25], method = "mahalanobis")
adonis2(dist_matrix ~ species, data = data)

#----------------------------------------------------------------------
# Principal component analysis on the abundance of different polymer types
#----------------------------------------------------------------------


res.pca <- PCA(data[,13:25], graph = FALSE) # Conduct a PCA on the data
PC1 <- res.pca$ind$coord[,1] #Store individual coordinates of PC1 as a vector
PC2 <- res.pca$ind$coord[,2] #Store individual coordinates of PC2 as a vector
PCs.ID <- data.frame(cbind(PC1,PC2)) #Bind the coordinates together as a dataframe
PCs.ID$species <- (data$species) #Add in species info to the dataset
PCs.ID$biome <- (data$biome) #Add in biome info to the dataset


#Define axis labels based on % of data explained across each dimension of the PCA
DIM_1 <- paste("PCA Dimension 1 (", round(res.pca$eig[1,2], 1), "%)")
DIM_2 <- paste("PCA Dimension 2  (", round(res.pca$eig[2,2], 1), "%)")


#Then make the figure
ggplot(PCs.ID, aes(x=PC1, y=PC2, color = species, fill = species), guide = "none") +
  stat_ellipse(geom = "polygon", alpha = 0.2)+
  geom_hline(aes(yintercept=0), linetype="dashed", lwd = 0.1) +
  geom_vline(aes(xintercept=0), linetype="dashed", lwd = 0.1) +
  theme_bw() +
  geom_point(size=0.3, aes(color = species)) +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73")) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73")) +
  ylab(DIM_2) +
  xlab(DIM_1) + 
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(size=8, family = "sans", face = "bold"),
        axis.title.y  = element_text(size=8, family = "sans", face = "bold"),
        plot.title = element_text(size=8, hjust = 0, family = "sans"),
        axis.text.y  = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans"),
        legend.position=c(0.75,0.9),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.3, "cm"),
        legend.key = element_blank())

ggsave(file="Figures/Species_PCA.png",
       width = 3.23,
       height=3,
       units = "in",
       dpi = 600)




#----------------------------------------------------------------------
# Random forest classifying species
#----------------------------------------------------------------------

#Run the random forest model identifying species
species.mod <- randomForest(y = as.factor(mp_data$species),
                            x = mp_data[,44:56],
                            mtry =  2,
                            ntree= 20000,
                            importance=TRUE,
                            proximity = TRUE,
                            keep.forest=TRUE,
                            replace = TRUE)

species.mod

#Check prediction accuracy
pred <- predict(species.mod, newdata = mp_data[,44:56])
confusionMatrix(as.factor(mp_data$species), as.factor(pred))

varImpPlot(species.mod, type=1, scale = FALSE)

#Get top 3 polymers
POLYMERS <- row.names(species.mod$importance)[order(species.mod$importance[,"MeanDecreaseAccuracy"],
                                                    decreasing = TRUE)][1:4]
POLYMERS

res.pca <- PCA(species.mod$proximity, graph = FALSE) # Conduct a PCA on the proximity matrix
PC1 <- res.pca$ind$coord[,1] #Store individual coordinates of PC1 as a vector
PC2 <- res.pca$ind$coord[,2] #Store individual coordinates of PC2 as a vector
PCs.ID <- data.frame(cbind(PC1,PC2)) #Bind the coordinates together as a dataframe
PCs.ID$species <- mp_data$species #Add in species to the df

#Define axis labels based on % of data explained across each dimension of the PCA
DIM_1 <- paste("Dim 1 (", round(res.pca$eig[1,2], 1), "%)")
DIM_2 <- paste("Dim 2 (", round(res.pca$eig[2,2], 1), "%)")


#Then make the figure
PCA_FIG <- 
  ggplot(PCs.ID, aes(x=PC1, y=PC2, color = species, fill = species)) +
  ggtitle("A")+
  geom_hline(aes(yintercept=0), linetype="dashed", lwd = 0.1, col = "grey70", alpha = 0.8) +
  geom_vline(aes(xintercept=0), linetype="dashed", lwd = 0.1, col = "grey70", alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.2, segments = 200, size = 0.1, show.legend = FALSE) +
  theme_bw() +
  ylab(DIM_2) +
  xlab(DIM_1) + 
  geom_point(size=0.1, aes(color = species)) +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73"),
                    guide = "none") +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73")) +
  ylab(DIM_2) +
  xlab(DIM_1) + 
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0, size = 10, family = "sans", face = "bold"),
        axis.title.x  = element_text(size=8, family = "sans", face = "bold"),
        axis.title.y  = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y  = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans"),
        legend.position=c(0.8,0.9),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.3, "cm"),
        legend.key = element_blank())

ggsave(PCA_FIG,
       file="Figures/Species_RF_PCA.png",
       width = 3.23,
       height=3,
       units = "in",
       dpi = 600)

#Figure of the primary polymers
X_LAB <- paste(POLYMERS[1], "Concentration (particles/mL)")

a <- 
  ggplot(mp_data, aes(x=mp_data[,POLYMERS[1]], y=mp_data$species, fill = mp_data$species)) +
  geom_density_ridges(scale = 5, alpha=0.6, linewidth = 0.2) +
  theme_ridges() +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x=X_LAB, y="Species")+
  ggtitle("B")+
  theme(plot.title = element_text(hjust = 0, size = 10, family = "sans", face = "bold"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size=5, family = "sans"),
        axis.title.y  = element_blank(),
        axis.text.y  = element_text(size=5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size=4, family = "sans"),
        legend.position=c("none"))

X_LAB <- paste(POLYMERS[2], "Concentration (particles/mL)")

b <- 
  ggplot(mp_data, aes(x=mp_data[,POLYMERS[2]], y=mp_data$species, fill = mp_data$species)) +
  geom_density_ridges(scale = 5, alpha=0.6, linewidth = 0.2) +
  theme_ridges() +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x=X_LAB, y="Species")+
  ggtitle("C")+
  theme(plot.title = element_text(hjust = 0, size = 10, family = "sans", face = "bold"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size=5, family = "sans"),
        axis.title.y  = element_blank(),
        axis.text.y  = element_text(size=5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size=4, family = "sans"),
        legend.position=c("none"))

X_LAB <- paste(POLYMERS[4], "Concentration (particles/mL)")

c <- 
  ggplot(mp_data, aes(x=mp_data[,POLYMERS[4]], y=mp_data$species, fill = mp_data$species)) +
  geom_density_ridges(scale = 5, alpha=0.6, linewidth = 0.2) +
  theme_ridges() +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x=X_LAB, y="Species")+
  ggtitle("D")+
  theme(plot.title = element_text(hjust = 0, size = 10, family = "sans", face = "bold"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size=5, family = "sans"),
        axis.title.y  = element_blank(),
        axis.text.y  = element_text(size=5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size=4, family = "sans"),
        legend.position=c("none"))

RIGHT <- arrangeGrob(a, b, c, ncol = 1)

FIG <- grid.arrange(PCA_FIG, RIGHT, ncol = 2, widths = c(0.6, 0.4))

ggsave(FIG,
       file="Figures/Concentration_Density_Plot.png",
       width = 6.86,
       height=4,
       units = "in",
       dpi = 600)






#----------------------------------------------------------------------
# Random forest classifying biome
#----------------------------------------------------------------------

#The data are too unbalanced for this model to run

BIOMES <- data$biome
BIOMES <- as.factor(BIOMES)

#Run the random forest model identifying species
biome.mod <- randomForest(y = BIOMES,
                          x = data[,14:25],
                          mtry =  3,
                          ntree= 20000,
                          importance=TRUE,
                          proximity = TRUE,
                          keep.forest=TRUE,
                          replace = TRUE)

biome.mod
varImpPlot(biome.mod, type=1, scale = FALSE)

#Check prediction accuracy
pred <- predict(biome.mod, newdata = data[,14:25])
confusionMatrix(BIOMES, as.factor(pred))

#Get top 3 polymers
POLYMERS <- row.names(biome.mod$importance)[order(biome.mod$importance[,"MeanDecreaseAccuracy"],
                                                  decreasing = TRUE)][1:4]
POLYMERS

res.pca <- PCA(biome.mod$proximity, graph = FALSE) # Conduct a PCA on the proximity matrix
PC1 <- res.pca$ind$coord[,1] #Store individual coordinates of PC1 as a vector
PC2 <- res.pca$ind$coord[,2] #Store individual coordinates of PC2 as a vector
PCs.ID <- data.frame(cbind(PC1,PC2)) #Bind the coordinates together as a dataframe
PCs.ID$Bioma <- BIOMES #Add in species to the df

#Define axis labels based on % of data explained across each dimension of the PCA
DIM_1 <- paste("Dim 1 (", round(res.pca$eig[1,2], 1), "%)")
DIM_2 <- paste("Dim 2 (", round(res.pca$eig[2,2], 1), "%)")



#Then make the figure
PCA_FIG <- 
  ggplot(PCs.ID, aes(x=PC1, y=PC2, color = Bioma, fill = Bioma)) +
  ggtitle("A") +
  geom_hline(aes(yintercept=0), linetype="dashed", lwd = 0.1, col = "grey70", alpha = 0.8) +
  geom_vline(aes(xintercept=0), linetype="dashed", lwd = 0.1, col = "grey70", alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.2, segments = 200, size = 0.1, show.legend = FALSE) +
  theme_bw() +
  ylab(DIM_2) +
  xlab(DIM_1) + 
  geom_point(size=0.1, aes(color = Bioma)) +
  scale_fill_manual(values = c("#57cc99", "#fca311", "#22577a"),
                    guide = "none") +
  scale_colour_manual(values = c("#57cc99", "#fca311", "#22577a")) +
  ylab(DIM_2) +
  xlab(DIM_1) + 
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0, size = 10, family = "sans", face = "bold"),
        axis.title.x  = element_text(size=8, family = "sans", face = "bold"),
        axis.title.y  = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y  = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans"),
        legend.position=c(0.1,0.9),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.3, "cm"),
        legend.key = element_blank())

ggsave(PCA_FIG,
       file="Figures/Biome_RF_PCA.png",
       width = 3.23,
       height=3,
       units = "in",
       dpi = 600)

#Figure of the primary polymers
X_LAB <- paste(POLYMERS[1], "Concentration (particles/mL)")

a <- 
  ggplot(data, aes(x=data[,POLYMERS[1]], y=data$Bioma, fill = data$Bioma)) +
  geom_density_ridges(scale = 5, alpha=0.6, linewidth = 0.2) +
  theme_ridges() +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x=X_LAB, y="Species")+
  ggtitle("A")+
  theme(plot.title = element_text(hjust = 0, size = 10, family = "sans", face = "bold"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size=5, family = "sans"),
        axis.title.y  = element_blank(),
        axis.text.y  = element_text(size=5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size=4, family = "sans"),
        legend.position=c("none"))

X_LAB <- paste(POLYMERS[2], "Concentration (particles/mL)")

b <- 
  ggplot(data, aes(x=data[,POLYMERS[2]], y=data$Bioma, fill = data$Bioma)) +
  geom_density_ridges(scale = 5, alpha=0.6, linewidth = 0.2) +
  theme_ridges() +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x=X_LAB, y="Species")+
  ggtitle("B")+
  theme(plot.title = element_text(hjust = 0, size = 10, family = "sans", face = "bold"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size=5, family = "sans"),
        axis.title.y  = element_blank(),
        axis.text.y  = element_text(size=5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size=4, family = "sans"),
        legend.position=c("none"))

X_LAB <- paste(POLYMERS[4], "Concentration (particles/mL)")

c <- 
  ggplot(data, aes(x=data[,POLYMERS[4]], y=data$Bioma, fill = data$Bioma)) +
  geom_density_ridges(scale = 5, alpha=0.6, linewidth = 0.2) +
  theme_ridges() +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x=X_LAB, y="Species")+
  ggtitle("C")+
  theme(plot.title = element_text(hjust = 0, size = 10, family = "sans", face = "bold"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size=5, family = "sans"),
        axis.title.y  = element_blank(),
        axis.text.y  = element_text(size=5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size=4, family = "sans"),
        legend.position=c("none"))

RIGHT <- arrangeGrob(a, b, c, ncol = 1)

FIG <- grid.arrange(PCA_FIG, RIGHT, ncol = 2, widths = c(0.6, 0.4))

ggsave(FIG,
       file="Figures/Biome_Differences.png",
       width = 6.86,
       height=4,
       units = "in",
       dpi = 600)

