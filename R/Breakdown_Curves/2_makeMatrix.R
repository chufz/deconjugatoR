makeMatrix <- function(feature, ppm_value){
  library(stringr)
  library(dplyr)
  files <- list.files(feature, pattern=".txt", full.names = T)
  files2 <- tools::file_path_sans_ext(list.files(feature, pattern=".txt"))
  list_spectra <- list()
  list_spectra_rounded <- list()
  # get collision enegy from filename
  ce <- str_sub(files2, -2, -1)
  #get first file of each collision energy
  ce_unique <- unique(ce)
  ce_unique_rows <- vector()
  for(f in 1:length(ce_unique)){
    ce_unique_rows[f] <- which(ce== ce_unique[f])[1]
  }
  # load in files
  for(f in 1:length(files)){
  list_spectra[[f]] <- read.table(files[f], header=F)
  colnames(list_spectra[[f]]) <- c("mz", "int")
  }
  # only keep one of the similar CE files
  list_spectra_unique <- list_spectra[ce_unique_rows]
  # bin all mz values together in each spectra
  mz_all <- unlist(lapply(list_spectra_unique, "[", 1))
  # generate mz that belong together
  ppm <- function(mass, ppm){return(mass*ppm*1e-6)}
  # generate unique mz values that can be aligned
  mz_all_min <- mz_all - ppm(mz_all, ppm_value)
  mz_all_max <-  mz_all + ppm(mz_all, ppm_value)
  mzcluster <- vector(length=length(mz_all))
  mzcluster[1] <- mz_all[1]
  for(i in 2:length(mz_all)){
      # check if there is a preexisting same mz value before
      mz_all_min_before <- mz_all_min[1:(i-1)]
      mz_all_max_before <- mz_all_max[1:(i-1)]
      
      x <- which(mz_all_min_before < mz_all[i])
      y <- which(mz_all_max_before[x] > mz_all[i])
      chose <- x[y]
      if(!length(y)==0){
          mzcluster[i] <- mz_all[chose[1]]
      }else{mzcluster[i] <- mz_all[i]}
  }
  # replace mz by unique mz
  a <- 1
  b <- 0
  for(i in 1:length(list_spectra_unique)){
   b <- b + nrow(list_spectra_unique[[i]])
   list_spectra_unique[[i]]$mz <- mzcluster[a:b]
   a <- a + nrow(list_spectra_unique[[i]])
  }
  # align same fragments
  matrix <- list_spectra_unique[[1]]
  for(i in 2:length(list_spectra_unique)){
    matrix <- full_join(matrix, list_spectra_unique[[i]], by="mz")
  }
  colnames(matrix) <- c("mz", ce_unique)
  return(matrix)
}