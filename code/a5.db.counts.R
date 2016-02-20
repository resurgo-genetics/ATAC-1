## ------- a5.db.counts.R
## last modified 2/16/16

library(DiffBind)

# create sample sheet
path <- "/groups/cbdm-db/jrd26/ATAC_crossbatch/R/Analysis5/"

Peaks <- list.files(paste0(path, "bed/rpm-filtered/"))

# create columns for sample sheet
SampleID <- unlist(strsplit(Peaks, split = ".2p5rpm.bed"))
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
PeakPath <- list.files(paste0(path, "bed/rpm-filtered/"), full.names = T)
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
ttreg <- dba.peakset(ttreg, consensus = -DBA_REPLICATE)

# make consensus peakset from samples and count reads
ttreg <- dba.count(ttreg, peaks = ttreg$masks$Consensus, score = DBA_SCORE_RPKM)

# save dba to file
dba.save(ttreg, file="ttreg", dir=path, pre="dba_", ext="RData", bMinimize=F)

# export consensus peakset with read counts
ttreg_con <- dba.peakset(ttreg, ttreg$masks$Consensus, bRetrieve = T, 
                         writeFile=paste0(path,"a5-ttreg-consensus.rpkm.txt"), DataType=DBA_DATA_FRAME)
write.csv(ttreg_con, file=paste0(path, "a5-ttreg-consensus.rpkm.csv"), row.names=F)

# save read count correlation heatmap
pdf(paste0(path,"a5.ttreg.reads.heat.pdf"), width=10, height=10, pagecentre = T)
par(oma = c(3,2,2,3))
dba.plotHeatmap(ttreg)
dev.off()