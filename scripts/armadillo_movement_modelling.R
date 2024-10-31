# This script processes the movement data for the giant armadillos
# The main steps carried out by this script are:
# 1. import the gps location data
# 2. set a prior for the measurement error
# 3. fit error informed continuous time movement models
# 4. extract land use
# 5. save all outputs for subsequent analyses


# Written by Michael Noonan

# Last updated: October 30 2024

# Load the requisite packages
library(ctmm) # for fitting the movement models
library(terra) # For handling the habitat rasters
library(weights) # For calculating the weighted proportions of the land class types

# Source any custom function
source("scripts/fit_mods.R") # simplifies analysing gps data for multiple individuals

#-------------------------------------------------------------
# Fit movement models and save information 
#-------------------------------------------------------------

# Import the tapir gps data and convert to a telemetry object
# Note: these data are not on GitHub
DATA <- read.csv("data/gps_data/armadillo_gps_data.csv")
DATA <- as.telemetry(DATA, mark.rm = TRUE)


#Define a prior for the error structure of the GPS location data
# first supply the point estimates
# 20-meter 2D error at HDOP=1
# 10-meter 3D error at HDOP=1
uere(DATA) <- c(20,10)

# extract the UERE object
PRIOR <- uere(DATA)
# set the DOF as 2 for wide credible intervals on the UERE
PRIOR$DOF[] <- 2
# assign prior to data
uere(DATA) <- PRIOR


#Plot to confirm the import worked as expected
plot(DATA, col = viridis::viridis(length(DATA)))

#Confirm which animals from the microplastics dataset have movement data
mp_data <- read.csv("data/mp_data/All_general.csv")
mp_IDs <- mp_data[which(mp_data$species == "Priodontes_maximus"),"name"]
mp_IDs[mp_IDs %in% names(DATA)]
mp_IDs[!mp_IDs %in% names(DATA)]

#2 animals have gps data, 1 is missing gps data (this is correct)


#Define the file paths for saving the movement model outputs
Model_path <- paste("~/Dropbox/UBC/Projects/microplastics_brazil/results/armadillo_movement/movement_models", sep = "")
HR_path <- paste("~/Dropbox/UBC/Projects/microplastics_brazil/results/armadillo_movement/akdes", sep = "")
Fig_Path <- paste("~/Dropbox/UBC/Projects/microplastics_brazil/results/armadillo_movement/figures", sep = "")


#Then walk through each individual sequentially
for(j in 1:length(DATA)){
  
  cat("Working on individual ", j, " of ", length(DATA), "\n")
  
  #Extract the current individual
  
  if(class(DATA) == "list") {cilla <- DATA[[j]]} else {cilla <- DATA}
  
  RESULTS <- tryCatch(
    {
      RESULTS <- CTMM_FIT(cilla,
                          Model_path = Model_path,
                          Fig_Path = Fig_Path,
                          HR_path = HR_path,
                          error = TRUE,
                          binomial = "Priodontes_maximus")
    }, error=function(err) {
      message(cat("Model fitting failed, returning NAs \n"))
      
      RESULTS <- as.data.frame(t(unlist(c("Priodontes_maximus",
                                          cilla@info[1],
                                          rep(NA, 23)))))
      
      return(RESULTS)
    }
  )
  
  if(j == 1){
    
    write.table(RESULTS,
                file = "~/Dropbox/UBC/Projects/microplastics_brazil/data/movement_data/armadillo_movement.csv",
                row.names=FALSE,
                col.names=TRUE,
                sep=",",
                append=TRUE)
    
  } else {
    
    write.table(RESULTS,
                file = "~/Dropbox/UBC/Projects/microplastics_brazil/data/movement_data/armadillo_movement.csv",
                row.names=FALSE,
                col.names=FALSE,
                sep=",",
                append=TRUE)
  }
  
}#Closes the loop that runs over the telemetry object



#-------------------------------------------------------------
# Extract land use
#-------------------------------------------------------------


# Import the human footprint index raster
# Note: these rasters are not on GitHub
HFI <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/HFI.tif")
land_types <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/brasil_coverage_2022.tif")


RES <- list()

# Then walk through each individuals for that species
for(i in 1:length(DATA)){
  
  # Generate a brief message to keep track of progress
  cat("Working on individual ", i, " of ", length(DATA), "\n")
  
  
  #Import the HR estimate for the ith animal
  PATH <- file.path("~/Dropbox/UBC/Projects/microplastics_brazil/results/armadillo_movement/akdes",
                    paste("AKDE_",
                          DATA[[i]]@info[1],
                          ".rda",
                          sep = ""))
  load(PATH)
  
  
  #Data carpentry to get the home range PDF into the correct format for extracting values
  HR <- rast(raster(AKDE, DF = "PMF"))
  HR <- project(HR, crs(land_types))
  HR.df <- terra::as.data.frame(HR, xy = TRUE, na.rm = TRUE)
  
  
  #Extract habitat values
  HR.df$HFI <- extract(HFI, HR.df[,1:2])[,2]
  HR.df$land_class <- extract(land_types, HR.df[,1:2])[,2]
  HR.df$land_class[HR.df$land_class %in% c("1","3", "4","5","6","49","29")] <- "Native_forest"
  #HR.df$land_class[HR.df$land_class %in% c("11")] <- "Wetland"
  HR.df$land_class[HR.df$land_class %in% c("9")] <- "Forestry"
  HR.df$land_class[HR.df$land_class %in% c("12","15")] <- "Pasture"
  HR.df$land_class[HR.df$land_class %in% c("18","19","20","21","39","40","41","62","36","46","47","35","48")] <- "Agriculture"
  HR.df$land_class[HR.df$land_class %in% c("24","25","30")] <- "Development"
  HR.df$land_class[HR.df$land_class %in% c("11","26","33")] <- "Water"
  
  # Use the home range PDF to calculate the weighted proportions of time spent the different land class types
  PROPS <- round(wpct(HR.df$land_class, HR.df$layer)*100,2)
  PROPS2 <- data.frame(class = names(PROPS),
                       proportion = as.numeric(PROPS))
  PROPS <- data.frame(t(PROPS2))[2,]
  names(PROPS) <- PROPS2$class
  
  
  res <- data.frame(binomial = "Priodontes_maximus")
  res$ID <- AKDE@info$identity
  res$mean_HFI <- sum(HR.df$layer*HR.df$HFI)
  res$max_HFI <- max(HR.df$HFI)
  res <- cbind(res,PROPS)
  
  RES[[i]] <- res
  
} # Closes the loop that runs over the telemetry object (i.e., i)


res <- do.call(dplyr::bind_rows, RES)
res[is.na(res)] <- 0

# Save the land use data as a csv
write.table(res,
            file = "~/Dropbox/UBC/Projects/microplastics_brazil/data/movement_data/armadillo_land_use.csv",
            row.names=FALSE,
            col.names=TRUE,
            sep=",")



