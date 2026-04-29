# This script imports the various microplastics datasets and runs some
# data carpentry steps to get them into the correct format for any subsequent analyses

# Written by Michael Noonan




#-------------------------------------------------------------
# Averages across all animals
#-------------------------------------------------------------

#Import the MP concentration dataset
mp_data <- read.csv("data/mp_data/blood_microplastics_individual_summary.csv")
mp_data$sex <- as.factor(mp_data$sex)
mp_data$species <- as.factor(mp_data$species)
mp_data[which(mp_data$name == "Alvinho/Alvaro"),"name"] <- "Alvinho_Alvaro"

#Import the estimated movement metrics for the 3 species
tapir_data <- read.csv("data/movement_summaries/tapir_movement_summary.csv")
anteater_data <- read.csv("data/movement_summaries/anteater_movement_summary.csv")
armadillo_data <- read.csv("data/movement_summaries/armadillo_movement_summary.csv")
move_data <- dplyr::bind_rows(tapir_data, anteater_data, armadillo_data)

#Import the estimated land use for the 3 species
tapir_data <- read.csv("data/movement_summaries/tapir_land_use_summary.csv")
anteater_data <- read.csv("data/movement_summaries/anteater_land_use_summary.csv")
armadillo_data <- read.csv("data/movement_summaries/armadillo_land_use_summary.csv")
land_use <- dplyr::bind_rows(tapir_data, anteater_data, armadillo_data)
land_use[is.na(land_use)] <- 0

#Merge with the microplastics information
mp_data <- merge(x = mp_data, y = move_data, by.x = c("name", "species"), by.y = c("ID", "binomial"), all.x = TRUE)
mp_data <- merge(x = mp_data, y = land_use, by.x = c("name", "species"), by.y = c("ID", "binomial"), all.x = TRUE)

#Convert hr to km^2
mp_data$hr <- mp_data$hr*1e-6

#Rescale distances from meters to km
mp_data$mean_dist_development <- mp_data$mean_dist_development/1000
mp_data$min_dist_development <- mp_data$min_dist_development/1000
mp_data$mean_dist_waste <- mp_data$mean_dist_waste/1000

#Import the polymer concentration data
polymers <- read.csv("data/mp_data/blood_microplastics_polymer_profiles.csv")

#Some data carpentry to get the format correct for analysis
polymers <- reshape(polymers, idvar = "Sample", timevar = "Microplastic", direction = "wide")

#Adjust the polymer names
names(polymers)[-1] <- gsub('Adjusted.', '', names(polymers)[-1])
names(polymers)[-1] <- stringr::str_to_title(names(polymers)[-1])
names(polymers)[2] <- "ABS"
names(polymers)[8] <- "Polypropylene"
names(polymers)[10] <- "PVC"
names(polymers)[13] <- "PET"

#Create a vector with the polymer names for future use
polymer_names <- names(polymers)[-1]

#Merge the polymer concentrations with the other information
mp_data <- merge(x = mp_data, y = polymers, by.x = "sample", by.y = "Sample")

#Pool the two giant armadillos with repeat sampling
stacy <- mp_data[mp_data$name == "Stacy",][1,]
stacy[,10:ncol(stacy)] <- colMeans(mp_data[mp_data$name == "Stacy",-c(1:9)], na.rm = TRUE)
mp_data[mp_data$name == "Stacy",] <- stacy

laura <- mp_data[mp_data$name == "Laura",][1,]
laura[,10:ncol(laura)] <- colMeans(mp_data[mp_data$name == "Laura",-c(1:9)], na.rm = TRUE)
mp_data[mp_data$name == "Laura",] <- laura

mp_data <- mp_data[!duplicated(mp_data),]

#Remove the intermediate objects from memory
rm(move_data)
rm(land_use)
rm(tapir_data)
rm(anteater_data)
rm(armadillo_data)
rm(polymers)
rm(stacy)
rm(laura)

#-------------------------------------------------------------
# Averages across all animals in long format
#-------------------------------------------------------------

#Load in the polymer data
polymers <- read.csv("data/mp_data/blood_microplastics_polymer_profiles.csv")
mp_data_long <- merge(x = mp_data[,c("sample","species", "biome")],
                      y = polymers,
                      by.x = "sample",
                      by.y = "Sample")

#Adjust the polymer names
mp_data_long$Microplastic <- stringr::str_to_title(mp_data_long$Microplastic)
mp_data_long$Microplastic[mp_data_long$Microplastic == "Abs"] <- "ABS"
mp_data_long$Microplastic[mp_data_long$Microplastic == "Polypropylene "] <- "Polypropylene"
mp_data_long$Microplastic[mp_data_long$Microplastic == "Pet"] <- "PET"
mp_data_long$Microplastic[mp_data_long$Microplastic == "Polyvinyl Chloride"] <- "PVC"

#Adjust the column names
names(mp_data_long)[4] <- "polymer"
names(mp_data_long)[5] <- "concentration"

#Remove the intermediate objects from memory
rm(polymers)



#-------------------------------------------------------------
# Particle Size data
#-------------------------------------------------------------


#Load in the particle size data
sizes <- read.csv("data/mp_data/blood_microplastics_particle_sizes.csv")

#Load in the polymer id key
key <- read.csv("data/mp_data/polymer_classification_key.csv")

sizes <- merge(x = key, y = sizes, by.x = "polymer_full", by.y = "Material")
sizes <- merge(x = mp_data[,c("sample","species")], y = sizes, by.x = "sample", by.y = "Sample")

#Remove the intermediate objects from memory
rm(key)



#-------------------------------------------------------------
# Individual Polymer Counts
#-------------------------------------------------------------

#Note: Builds off of the sizes data frame generated above

#Build a dataframe of the counts of polymers in each sample
polymers <- data.frame(table(sizes[,c("sample","polymer")]))

#Reshape from long to wide
polymers <- reshape(polymers, idvar = "sample", timevar = "polymer", direction = "wide")

#Rename colums
names(polymers)[2:14] <- colnames(table(sizes[,c("sample","polymer")]))

#Add in meta data on biomes and species
polymers <- merge(x = mp_data[,c("sample","species","biome")], y = polymers, by.x = "sample", by.y = "sample")


#-------------------------------------------------------------
# Blank Control data
#-------------------------------------------------------------


#Load in the particle size data in the controls
<<<<<<< HEAD
control <- read.csv("data/mp_data/blank_controls_particle_sizes.csv")[,1:4]
=======
control <- read.csv("data/mp_data/blank_controls_particle_sizes.csv")[,1:3]
>>>>>>> 10c3184d928c37a1697cb86521eefd8d627b0ffd

#Clean up the polymer names
polymer_names <- c(
  "POLYVINYL CHLORIDE" = "PVC",
  "POLYESTER" = "Polyester",
  "Silicon" = "Silicon",
  "POLYPROPYLENE " = "Polypropylene", # trailing space preserved intentionally
  "Other" = "Other",
  "PA" = "Polyamide",
  "TEFLON" = "Teflon",
  "ABS" = "ABS",
  "POLYSTYRENE" = "Polystyrene",
  "POLYAMIDE" = "Polyamide"
)

control$polymer <- polymer_names[control$Polymer]

#Rename columns
names(control)[1] <- "length"
names(control)[2] <- "width"
names(control)[3] <- "sample"


#Remove the intermediate objects from memory
control$Polymer <- NULL
rm(polymer_names)



polymers <- c("ABS", "Cellulose", "PET", "Polyamide", "Polyester",
              "Polyethylene", "Polypropylene", "Polystyrene",
              "PVC", "Rubber", "Silicon", "Teflon", "Other")

#---------------------
#Get polymer counts per sample
#---------------------

control_polymers <- aggregate(width ~ polymer + sample, data = control, FUN = length)
names(control_polymers)[3] <- "count"

# Build complete grid of all polymer x sample combinations
complete <- expand.grid(polymer = polymers, sample = unique(control$sample),
                        stringsAsFactors = FALSE)

# Merge and fill zeros
control_polymers <- merge(complete, control_polymers, by = c("polymer", "sample"), all.x = TRUE)
control_polymers$count[is.na(control_polymers$count)] <- 0

#Remove the intermediate objects from memory
rm(complete)