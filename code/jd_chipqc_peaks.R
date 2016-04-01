# load libraries
library(ChIPQC)

# function to generate chipqc object and peak table with read counts in rpm 
# exports ChIPQCsample object and peak table
# ARGS: 

peak_rpm <- function(reads, peaks) {
  # create chromosome vector
  chr <- paste("chr",1:19,sep="")
  chr <- c(chr, "chrX","chrY")
  
  # create chipQCsample object
  qcs <- ChIPQCsample(reads, peaks, annotation = "mm9", chromosomes=chr)
  
  # extract total reads
  tr <- QCmetrics(qcs)[[1]]
  
  # extract peak table with counts
  qcs_peaks <- peaks(qcs)
  
  # convert to dataframe
  qcs_peaks <- as.data.frame(qcs_peaks)
  
  # append total reads
  qcs_peaks$Total_Rds <- tr
  
  # calculate and append rpm
  qcs_peaks$rpm <- (qcs_peaks$Counts / qcs_peaks$Total_Rds) * 1000000
    
  # extract df name 
  Sample <- strsplit(peaks, split=".peaks.atac.bed.txt")
  Sample <- unlist(Sample)
  
  # append name
  qcs_peaks <- cbind(Sample, qcs_peaks)
  
  # write out csv file
  write.csv(qcs_peaks, file=paste(Sample,".csv", sep=""), row.names=F)
  
  # write ChIPQCsample
  save(qcs, file=paste(Sample,".RDa",sep=""))
}

#### Function implementation ####
# set working directory
setwd("/groups/cbdm-db/jrd26/ATAC_crossbatch/temp/")

# get names of peak files
peaks <- list.files(".", pattern = "*.bed.txt")

# get names of bam files
reads <- list.files(".", pattern = "*.uniq.bam$")

# loop over each peak/bam pair in folder
for(i in 1:length(reads)) {
  peak_rpm(reads[i], peaks[i])
}



