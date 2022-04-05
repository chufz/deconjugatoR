###############################
mzmlFilesDir <- "" # Change here for foldername with mzML files in each folder
###############################
PrecMassesFile <- paste0(mzmlFilesDir, "/Features.txt") # File in each folder, with a inclusion list of "mz@RT" values
###############################
#Load in Functions
source("1_extractMSMS.R")
source("2_makeMatrix.R")
source("3_plotBreakdown.R")
###############################
#Extract related MS2 spectra
extractMSMS(mzmlFilesDir, PrecMassesFile)
#list Folders with Features
features <- list.dirs(mzmlFilesDir, full.names = T)
features <- features[-1]
# For each feature, align Fragments and Plot Breakdown Curve
for(feature in features){
    matrix <- makeMatrix(feature, ppm_value=25)
    plotBreakdown(matrix, outputdir=feature)
}

