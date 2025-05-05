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
# Correlation with age in tapirs
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



a <- 
  ggplot(data = mp_data[mp_data$species == "Tapirus_terrestris",], aes(x = age, y = mp_ml)) +
  ggtitle("A)") +
  geom_smooth(method = "gam",
              formula = y ~ x,
              method.args = list(family = tw(link = "log")),
              col = "black",
              fill = "grey80",
              linewidth = 0.2,
              linetype = "solid") +
  geom_point(col = "#005f73", size = 0.4) +
  ylab("MP concentration (particles/mL)") +
  xlab(expression(bold(Age~(Years))))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        legend.text  = element_text(size=5, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "inside",
        legend.position.inside  = c(0.7,0.95),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  coord_cartesian(ylim = c(50,1250))









#-------------------------------------------------------------
# Correlation with body mass
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

#-------------------------------------------------------------
# Differences between biomes
#-------------------------------------------------------------

#Test for differences between biomes in tapirs
fit <- gam(mp_ml ~ biome,
           family = tw(link = "log"),
           data = mp_data[mp_data$species == "Tapirus_terrestris",],
           method = "REML")

summary(fit)

d <- 
  ggplot(data = mp_data, aes(x = max_HFI, y = mp_ml)) +
  ggtitle("D)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "solid") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("MP concentration (particles/mL)") +
  xlab("Maximum human footprint index in home range") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  coord_cartesian(ylim = c(5,170))



#-------------------------------------------------------------
# Correlation with agricultural land
#-------------------------------------------------------------


#Test for a correlation between the ammount of agg land in the HR and blood MP concentration
fit <- gam(mp_ml ~ Agriculture + s(species, bs = "re") + s(Agriculture, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)

e <- 
  ggplot(data = mp_data, aes(x = Agriculture, y = mp_ml)) +
  ggtitle("E)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("MP concentration (particles/mL)") +
  xlab("Proportion of home range in agricultural land") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  coord_cartesian(ylim = c(5,170))




#-------------------------------------------------------------
# Correlation with water and wetlands
#-------------------------------------------------------------


#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ Water + s(species, bs = "re") + s(Water, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


f <- 
  ggplot(data = mp_data, aes(x = Water, y = mp_ml)) +
  ggtitle("F)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "solid") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("MP concentration (particles/mL)") +
  xlab("Proportion of water and wetlands in home range") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  coord_cartesian(ylim = c(5,170))



#-------------------------------------------------------------
# Correlation with water and wetlands
#-------------------------------------------------------------


#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ Native_forest + s(species, bs = "re") + s(Native_forest, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


g <- 
  ggplot(data = mp_data, aes(x = Native_forest, y = mp_ml)) +
  ggtitle("G)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("MP concentration (particles/mL)") +
  xlab("Proportion of natural forest in home range") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  coord_cartesian(ylim = c(5,170))


FIG <-
  grid.arrange(a,b,c,d,e,f,
               ncol=2,
               nrow=3)


#Save the figures
ggsave(FIG,
       width = 4.75, height = 5.1, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/Movement_Trends.png")

