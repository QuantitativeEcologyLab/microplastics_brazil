

#Import the MP datasets
source("scripts/data_import.R")




#---------------------------------------------------------------------
# Figure 3A; correlation with age in tapirs
#---------------------------------------------------------------------


a <-
ggplot(data = mp_data[which(mp_data$species == "Tapirus_terrestris"),],
       aes(x = age,
           y = mp_ml+1,
           col = sex,
           fill = sex)) +
  ggtitle("A") +
  geom_smooth(method = "gam",
              formula = y ~ x,
              method.args = list(family = tw(link = "log")),
              linewidth = 0.2,
              linetype = "solid",
              show.legend = FALSE) +
  geom_point(size = 0.4) +
  scale_colour_manual(values = c("#fca311", "#14213d"), labels = c("Female", "Male"), name = "") +
  scale_fill_manual(values = c("#fca311", "#14213d"), labels = c("Female", "Male"), name = "") +
  xlab("Age (years)") +
  ylab("Microplastic concentration (particles/mL)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_log10(breaks = c(1,10+1,100+1,1000+1,10000,100000,1000000,10000000,100000000,1000000000,10000000000,100000000000,1000000000000,10000000000000),
                labels = c(0,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000,10000000000,100000000000,1000000000000,10000000000000)) +
  annotation_logticks(#outside = TRUE,
    sides = "l",
    linewidth = 0.1,
    short = unit(0.05, "cm"),
    mid = unit(0.05, "cm"),
    long = unit(0.08, "cm"))



#---------------------------------------------------------------------
# Figure 3B-D; differences between sexes
#---------------------------------------------------------------------

#Some processing to get the names in the right order/format
mp_data$species_common <- as.character(mp_data$species)
mp_data[mp_data$species_common == "Tapirus_terrestris","species_common"] <- "Lowland tapir"
mp_data[mp_data$species_common == "Myrmecophaga_tridactyla","species_common"] <- "Giant anteater"
mp_data[mp_data$species_common == "Priodontes_maximus","species_common"] <- "Giant armadillo"
mp_data$species_common <- factor(mp_data$species_common, levels = c("Lowland tapir", "Giant armadillo", "Giant anteater"), ordered = TRUE)

#Generate the figure
mid <-
ggplot(data = mp_data, aes(x = sex,
                           y = mp_ml+1,
                           col = sex,
                           fill = sex,
                           alpha = 0.5)) +
  ggtitle("B") +
  geom_boxplot(size = 0.1, outlier.size = 0.2, outlier.shape = 16, outlier.alpha = 0) +
  geom_jitter(size = 0.5, shape = 16, position=position_jitter(height=0, width=0.2)) +
  scale_fill_manual(values = c("#fca311", "#14213d")) +
  scale_colour_manual(values = c("#fca311", "#14213d")) +
  facet_wrap(vars(species_common),  ncol = 3) +
  ylab("Microplastics (particles/mL)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=4, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        #strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_x_discrete(breaks = c("female", "male"), labels = c("Female","Male")) +
  scale_y_log10(breaks = c(1,10+1,100+1,1000+1),
                labels = c(0,10,100,1000),
                limits = c(1,1000)) +
  annotation_logticks(#outside = TRUE,
    sides = "l",
    linewidth = 0.1,
    short = unit(0.05, "cm"),
    mid = unit(0.05, "cm"),
    long = unit(0.08, "cm"))





#---------------------------------------------------------------------
# Figure 3E-G; correlations with body weight
#---------------------------------------------------------------------


#Generate the figure
bot <-
  ggplot(data = mp_data, aes(x = weight,
                             y = mp_ml+1,
                             col = sex,
                             fill = sex,
                             alpha = 0.5)) +
  ggtitle("E") +
  geom_smooth(method = "gam",
              formula = y ~ x,
              method.args = list(family = tw(link = "log")),
              linewidth = 0.2,
              linetype = "solid",
              show.legend = FALSE) +
  geom_point(size = 0.5, shape = 16) +
  scale_fill_manual(values = c("#fca311", "#14213d")) +
  scale_colour_manual(values = c("#fca311", "#14213d")) +
  facet_wrap(vars(species_common),  ncol = 3, scales = "free_x") +
  ylab("Microplastics (particles/mL)") +
  xlab("Body weight (kg)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=4, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        #strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_log10(breaks = c(1,10+1,100+1,1000+1,10000+1),
                labels = c(0,10,100,1000,10000)) +
  annotation_logticks(#outside = TRUE,
    sides = "l",
    linewidth = 0.1,
    short = unit(0.05, "cm"),
    mid = unit(0.05, "cm"),
    long = unit(0.08, "cm"))




#---------------------------------------------------------------------
# Assemble and save
#---------------------------------------------------------------------


figure_3 <-
  grid.arrange(a,mid,bot,
               ncol=1,
               nrow=3,
               heights=c(2,1,1))


#Save the figures
ggsave(figure_3,
       width = 4.75, height = 5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_3.png")

