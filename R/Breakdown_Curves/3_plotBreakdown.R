plotBreakdown <- function(matrix, outputdir, frag=8 ){
  library(ggplot2)
  library(tidyr)
  # replace NA with zero
  matrix[is.na(matrix)] <- 0
  # trunctuate matrix by fragments which have the 8 highest max values
  max_mz <- apply(matrix, 1, max)
  max_mz2 <- sort(max_mz, decreasing = T)[1:frag]
  matrix <- matrix[which(max_mz %in% max_mz2),]
  # generate ggplotstyle matrix
  gg_matrix <- gather(matrix, key="CE" , value="Intensity", -mz )
  gg_matrix$mz <- as.factor(round(gg_matrix$mz,4))
  ggplot(gg_matrix,aes(CE, Intensity, group=mz, color=mz)) +
    geom_point() +
    theme_classic() +
    ggtitle(basename(outputdir))
  ggsave(paste0(outputdir,"/Plot.png"))
}