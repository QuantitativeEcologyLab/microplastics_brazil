# This script processes runs a number of sensitivity analyses
# on the results presented in the main text.
# Where relevant, comments direct to the core analyses for which sensitivity is
# being assessed.

# Written by Michael Noonan


#Load in any necessary packages
library(mgcv)

#Import the MP datasets
source("scripts/04_microplastics_data_import.R")
