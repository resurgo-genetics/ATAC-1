## ------- a5.db.counts.R
## last modified 2/25/16

library(DiffBind)

##----- create sample sheet --------------------------------------------------------------------------------
path <- "/groups/cbdm-db/jrd26/ATAC_crossbatch/R/Analysis5/"

Peaks <- list.files(paste0(path, "bed/2rpm-filtered/"))

# create columns for sample sheet
SampleID <- unlist(strsplit(Peaks, split = "-2rpm.bed"))
Tissue <- c(rep('Colon',4),rep('Muscle',2),rep('Spleen',8),rep("VAT",2))
Factor <- unlist(strsplit(SampleID, split = "_rep1.2"))
Factor <- unlist(strsplit(Factor, split = "_rep1"))
Factor <- unlist(strsplit(Factor, split = "_rep2"))
Factor <- unlist(strsplit(Factor, split = "_rep3"))
Condition <- "Treg"
Treatment <- c('Batch1.2','Batch2','Batch1.2','Batch2',rep('Batch2',2),'Batch1','Batch2',
               'Batch1','Batch2',rep('Batch2',2),rep('Batch1',4))
Replicate <- rep(c('1','2'),8)
PeakCaller <- "homer"
PeakFormat <- "bed"
PeakPath <- list.files(paste0(path, "bed/2rpm-filtered/"), full.names = T)
bamPath <- list.files(paste0(path,"reads/"), full.names=T)

samples <- cbind(SampleID = SampleID, 
                 Tissue = Tissue, 
                 Factor = Factor, 
                 Condition = Condition,
                 Treatment = Treatment,
                 Replicate = Replicate,
                 bamReads = bamPath,
                 Peaks = PeakPath,
                 PeakCaller = PeakCaller,
                 PeakFormat = PeakFormat)
write.csv(samples, "samples.csv", row.names=F)


##-----loading data-----------------------------------------------------------------------------------
# create DBA object
ttreg <- dba(sampleSheet = "samples.csv", bRemoveRandom = T, minOverlap = 2)

# make consensus peakset from replicates
ttreg <- dba.peakset(ttreg, consensus = DBA_FACTOR)

# make consensus peakset from samples and count reads
ttreg <- dba.count(ttreg, peaks = ttreg$masks$Consensus, score = DBA_SCORE_RPKM)

# save dba to file
n <- 2  # save rpm filter value
dba.save(ttreg, file=paste0("a5-ttreg",n,'rpm-filtered'), dir=path, pre="dba_", ext="RData", bMinimize=F)

# export consensus peakset with read counts
# set filter
n <- 2
suffix <- paste('a5-ttreg-consensus-',n,'rpm-filtered-rpkm',sep="")
ttreg_con <- dba.peakset(ttreg, ttreg$masks$Consensus, bRetrieve = T, 
                         writeFile=paste0(path,suffix,'.txt'), DataType=DBA_DATA_FRAME)

write.csv(ttreg_con, file=paste0(path, suffix, '.csv'), row.names=F)

# save read count correlation heatmap
suffix <- paste('a5-ttreg-',n,'rpm-reads-heat',sep="")
pdf(paste0(path,suffix,'.pdf'), width=10, height=10, pagecentre = T)
par(oma = c(3,2,2,3))
dba.plotHeatmap(ttreg)
dev.off()

# save read count PCA
suffix <- paste('a5-ttreg-',n,'rpm-pca',sep="")
pdf(paste0(path,suffix,'.pdf'), width=10, height=10, pagecentre = T)
par(oma = c(3,2,2,3))
dba.plotPCA(ttreg, attributes=DBA_TISSUE, label=DBA_FACTOR, score=DBA_SCORE_RPKM)
dev.off()
