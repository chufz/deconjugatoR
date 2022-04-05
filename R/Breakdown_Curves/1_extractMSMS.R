extractMSMS <- function(sampledir, PrecMassesFile, mzdiff=0.001, rtdiff=30, puritytresh=0.5){
  ########################################
  if(!requireNamespace("msPurity", quietly = TRUE)) BiocManager::install("msPurity") # will install msPurity in case it is not installed yet
  library(stringr)
  library(msPurity)
  library(mzR)
  ########################################
  message("[DDextract] Read in data")
  peakids <- read.table(file=PrecMassesFile, header= F)
  mz <- as.numeric(unlist(lapply(strsplit(peakids[,1], "@"), function(x){x[1]})))
  rt <- as.numeric(sapply(strsplit(peakids[,1], "@"),"[[",2))
  message("[DDextract] Get mzML files")
  mzMLpths <- list.files(path= sampledir, pattern=".mzML", full.names = T)
  ########################################
  message("[DDextract] Run MSPurity, please wait...")
  pa <- purityA(mzMLpths)
  ms_result <- pa@puritydf
  ########################################
  message("[DDextract] Search in result for Peakid- related entries")
  res <- list()
  for(i in 1:length(mz)){
    mzwindow <- c(mz[i]-mzdiff,mz[i]+mzdiff)
    rtwindow <- c(rt[i]*60-rtdiff,rt[i]*60+rtdiff)
    scans <- which(ms_result$precursorMZ > mzwindow[1] & ms_result$precursorMZ < mzwindow[2] & ms_result$precursorRT > rtwindow[1] & ms_result$precursorRT < rtwindow[2])
    res[[i]] <- ms_result[scans,]
  }
  ########################################
  message("[DDextract] Create folder for each Peakid where MSMS was found, store a resultfile")
  message("[DDextract] Extract MS2 scans and save as txt")
  dirnames <- str_replace(peakids[,1], "@", "_")
  no_scans <- unlist(lapply(res, nrow))
  if(sum(no_scans)> 0){
    for(i in 1:length(res)){
      if(no_scans[i]>0){
        #create directory
        if(!dir.exists(paste0(sampledir,"/", dirnames[i]))){dir.create(paste0(sampledir,"/", dirnames[i]))}
        #save summary file
        write.csv(res[[i]], paste0(sampledir,"/", dirnames[i], "/Result_mspurity.csv"))
        #extract MSMS
        for(j in 1:nrow(res[[i]])){
          if(res[[i]]$aPurity[j]> puritytresh){
            aq_scan <- res[[i]]$acquisitionNum[j]
            file <- paste0(sampledir,res[[i]]$filename[j])
            file2 <- tools::file_path_sans_ext(basename(file))
            #open file
            mzml <- openMSfile(file)
            pl <- mzR::peaks(mzml, aq_scan)
            ce <- mzR::header(mzml)$collisionEnergy[aq_scan]
            write.table(pl, file=paste0(sampledir,"/", dirnames[i],"/",file2,"_",aq_scan,"_", ce, ".txt"), row.names=F, col.names=F)
            rm(mzml)
          }else{message("[DDextract] MSMS purity score is lower than applied threshold")}
        }
      }
    }
  }else{message("[DDextract] No MSMS has been found to the associated Peakids")}
  ########################################
}
