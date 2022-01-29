###############################################
### Script for loading MS2 data from one file #
### truncating by a specific neutral loss ###
### performing in-silico deconjugation ########
### performing spectral library search ########
### using MassBank as MSQL query #############
###############################################
### Set variables => Change for your system ###
file_mzml <- "testfile.mzML" # change for your mzML File
output_csv <- "librarysearch_result.csv" #Directory where output will be stored
output_pdf <- "librarysearch_result.pdf"
NL_screening <- TRUE # If TRUE, spectra without a NL will be removed
insilico_deconjugation <- TRUE # Shall in-silico deconjugation be performed?, For FALSE, normal library search will be preformed for all spectra containing the neutral loss
plot_Spectra <- TRUE # Shall spectra be plotted?
NL <- 176.0321 # NL for Glucuronide
mz_ppm <- 10 # Mass accuracy for NL search
int_tresh <- 5 # Clean spectra, removing fragments below x% of the base peak intensity after Nl search
dp_tresh <- 0.2 #  Dot product score which is considered as spectral match 
###############################################
### Load or install packages ##################
###############################################
if (!require("RMariaDB")) install.packages("RMariaDB")
library(RMariaDB)

if (!require("pander")) install.packages("pander")
library(pander)

if (!require("devtools")) install.packages("devtools")

if (!require("Spectra")) devtools::install_github("rformassspectrometry/Spectra")
library(Spectra)

if (!require("MetaboAnnotation")) devtools::install_github("rformassspectrometry/MetaboAnnotation")
library(MetaboAnnotation)

if (!require("MsBackendMassbank")) devtools::install_github("rformassspectrometry/MsBackendMassbank")
library(MsBackendMassbank)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("mzR")) BiocManager::install("mzR")

###############################################
### Workflow ##################################
###############################################
# Load the data
ms2data <- Spectra(file_mzml)
# restrict to MS2 data
ms2data <- filterMsLevel(ms2data, msLevel = 2L)
# neutral loss screening
if(NL_screening){
    # remove MS2 with Precursor smaller than NL +70, otherwise negative Precursor masses exist in data
    ms2data <- filterPrecursorMz(ms2data, c((NL + 70), 800))
    # restrict to MS2 with abundant Neutral Loss
    ms2data <- ms2data[which(containsNeutralLoss(ms2data, neutralLoss = NL, ppm = mz_ppm ))]
    if(length(ms2data) < 1){
        "No spectra with associated neutral loss found." 
    }else{ message("Number of spectra containing a neutral loss: ", length(ms2data))}
}
# change to write-able MsBackend
ms2data <- setBackend(ms2data, backend = MsBackendDataFrame()) 
#  change precursor mass, remove fragments above precursor mass
if(insilico_deconjugation){
    # Change precursorMass to deconjugated mass
    ms2data@backend@spectraData@listData$precursorMz <- precursorMz(ms2data) - NL
    # define modification function for truncating spectra
    removeAbovePrecursor <- function() {
        function(x, precursorMz, ...) {
            x[!(x[,1] >= precursorMz),,drop=FALSE]
        }
    }
    # apply modification
    ms2data <- addProcessing(ms2data, removeAbovePrecursor(), spectraVariables =c("precursorMz"))
    ms2data <- applyProcessing(ms2data, parallel = FALSE)
}
# remove fragments below x% of base peak intensity
low_int <- function(x, ...) {
    x > max(x, na.rm = TRUE) * (int_tresh / 100 )
}
ms2data <- filterIntensity(ms2data, intensity = low_int)
# normalize intensities
norm_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}
ms2data <- addProcessing(ms2data, norm_int)
# Library search, load MassBank, Connect to the MassBank MySQL 
con <- dbConnect(MariaDB(), user = "massbank", dbname = "MassBank",
                 host = "localhost", pass = "massbank")
mbank <- Spectra(con, source = MsBackendMassbankSql())
mbank <- addProcessing(mbank, norm_int)
# Perform library search in MassBank
prm <- MatchForwardReverseParam(ppm = mz_ppm, #requirePrecursor = TRUE,
                           THRESHFUN = function(x) which(x >= dp_tresh))
mtch <- matchSpectra(ms2data, mbank, param = prm)
mtch_sub <- mtch[whichQuery(mtch)]
# Write results in console
pandoc.table(style = "rmarkdown",
             as.data.frame(spectraData(mtch_sub, c("rtime", "target_compound_name",
                                             "score"))))
# Save output in csv file
write.csv2(spectraData(mtch_sub, c("rtime", "target_compound_name",
                                  "score", "reverse_score",  "collisionEnergy", "precursorMz", "precursorIntensity", "peaksCount", "dataOrigin", "target_splash", "target_smiles")), 
          output_csv)

# plot Head Tail plots in pdf with the corresponding Molecular Structure
if(plot_Spectra){
    pdf(output_pdf)
    splash <- spectraData(target(mtch_sub), "splash") 
    for(i in 1:length(mtch_sub)){
        highest_score  <- which.max(spectraData(mtch_sub[i], "score")[,1])
        targetid <- spectraData(mtch_sub[i], "target_splash")[highest_score,1]
        plotSpectraMirror(query(mtch_sub)[i], target(mtch_sub)[which(splash$splash == targetid)],
            title(spectraData(mtch_sub[i], "target_spectrum_name")[highest_score,1]))
    }
    dev.off()
}

