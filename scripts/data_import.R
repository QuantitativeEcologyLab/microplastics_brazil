# This script imports the various microplastics datasets and runs some
# data carpentry steps to get them into the correct format for any subsequent analyses

# Written by Michael Noonan




#-------------------------------------------------------------
# Averages across all animals
#-------------------------------------------------------------

#Import the MP concentration dataset
mp_data <- read.csv("data/mp_data/All_general.csv")
mp_data$sex <- as.factor(mp_data$sex)
mp_data$species <- as.factor(mp_data$species)
mp_data[which(mp_data$name == "Alvinho/Alvaro"),"name"] <- "Alvinho_Alvaro"

#Import the estimated movement metrics for the 3 species
tapir_data <- read.csv("data/movement_data/tapir_movement.csv")
anteater_data <- read.csv("data/movement_data/anteater_movement.csv")
armadillo_data <- read.csv("data/movement_data/armadillo_movement.csv")
move_data <- dplyr::bind_rows(tapir_data, anteater_data, armadillo_data)

#Import the estimated land use for the 3 species
tapir_data <- read.csv("data/movement_data/tapir_land_use.csv")
anteater_data <- read.csv("data/movement_data/anteater_land_use.csv")
armadillo_data <- read.csv("data/movement_data/armadillo_land_use.csv")
land_use <- dplyr::bind_rows(tapir_data, anteater_data, armadillo_data)
land_use[is.na(land_use)] <- 0

#Merge with the microplastics information
mp_data <- merge(x = mp_data, y = move_data, by.x = c("name", "species"), by.y = c("ID", "binomial"), all.x = TRUE)
mp_data <- merge(x = mp_data, y = land_use, by.x = c("name", "species"), by.y = c("ID", "binomial"), all.x = TRUE)

#Convert hr to km^2
mp_data$hr <- mp_data$hr*1e-6


#Import the polymer concentration data
polymers <- read.csv("data/mp_data/All_PP.csv")

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


#Remove the intermediate objects from memory
rm(move_data)
rm(land_use)
rm(tapir_data)
rm(anteater_data)
rm(armadillo_data)
rm(polymers)


#-------------------------------------------------------------
# Averages across all animals in long format
#-------------------------------------------------------------

#Load in the polymer data
polymers <- read.csv("data/mp_data/All_PP.csv")
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
sizes <- read.csv("data/mp_data/All_size.csv")

#Load in the polymer id key
key <- read.csv("data/mp_data/polymer_key.csv")

sizes <- merge(x = key, y = sizes, by.x = "polymer_full", by.y = "Material")
sizes <- merge(x =  mp_data[,c("sample","species")], y = sizes, by.x = "sample", by.y = "Sample")

rm(key)




#-------------------------------------------------------------
# Blank Control Size data
#-------------------------------------------------------------


#Load in the particle size data in the controls
control <- read.csv("data/mp_data/Control_particles_size.csv")[,1:3]

