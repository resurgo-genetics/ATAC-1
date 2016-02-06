#!/bin/bash
#
# atac_mapping_PE_2.sh
# AUTHOR: JRD/DPZ
# LAST MODIFIED: 12/19/15
#
# run where the fastq.gz files are
# may have to make executable with chmod +x atac_mapping_PE.sh
# usage: ./atac_mapping_PE.sh  dir_where_fastq.gz_files_are prefix


#### Load libraries and modules
source ~/.bash_profile
export R_LIBS="/groups/cbdm_lab/dp133/R_libraries"

module load seq/fastqc/0.10.1
module load seq/cutadapt/1.8.3 
module load seq/fastx/0.0.13
module load seq/bowtie/2.2.4
module load seq/samtools/0.1.19
module load seq/BEDTools/2.19.0
module load stats/R/3.1.2
module load seq/sickle/1.2


#### Set working directory and prefix for output files (args)
echo "working directory : " $1
echo "prefix : " $2
prefix=$2

cd $1


#### Gunzip files
echo "Decompressing files..."
gunzip -fv *R1.fastq.gz
gunzip -fv *R2.fastq.gz
mv *R1.fastq $prefix.R1.fastq
mv *R2.fastq $prefix.R2.fastq


#### Generate fastqc output
echo "Running fastqc on Read1..."
mkdir $prefix.R1.fastqc/
fastqc -o $prefix.R1.fastqc/ $prefix.R1.fastq

echo "Running fastqc on Read2..."
mkdir $prefix.R2.fastqc/
fastqc -o $prefix.R2.fastqc/ $prefix.R2.fastq


#### Filter reads on quality
echo "Filtering reads on quality..."
sickle pe -f $prefix.R1.fastq -r $prefix.R2.fastq -o $prefix.filtered_R1.fastq -p $prefix.filtered_R2.fastq -t sanger -s $prefix.filtered_singles.fastq


#### Clip adapters from sequences
# e = maximum error rate, default = 0.1
# m = minimum length, throw away reads shorter than N bases
echo "Clipping adapter from 5' side..."
cutadapt -g AGATGTGTATAAGAGACAG -G CTGTCTCTTATACACATCT -e 0.1 -m 20 -o $prefix.trim1_R1.fastq -p $prefix.trim1_R2.fastq $prefix.filtered_R1.fastq $prefix.filtered_R2.fastq

echo "Clipping adapter from 3' side..."
cutadapt -a AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT -e 0.1 -m 20 -o $prefix.trim2_R1.fastq -p $prefix.trim2_R2.fastq $prefix.trim1_R1.fastq $prefix.trim1_R2.fastq


#### Align reads and keep only reads mapping to one location
# -p N		number of threads to use for multi-thread software
# -x		location of indexed genome
# -S <>		name of .sam output file
echo "Mapping reads to mm9..."
bowtie2 -p 8 -x /groups/shared_databases/bowtie2_indexes/mm9 -X 1000 --fr -1 $prefix.trim2_R1.fastq -2 $prefix.trim2_R2.fastq -S $prefix.btout2.sam


#### Keep only reads mapping to one location
samtools view -hS -F 4 $prefix.btout2.sam > $prefix.mapped.sam # keep only mapped reads
sed '/XS:/d' $prefix.mapped.sam > $prefix.mapped_1alignmentonly.sam # remove multiple aligned reads
samtools view -bS $prefix.mapped_1alignmentonly.sam > $prefix.mapped_1align.bam # make a bam file


### Remove duplicates 
echo "Removing duplicates..."
samtools sort $prefix.mapped_1align.bam $prefix.mapped_1align.sorted
java -Xms1024m -jar /opt/picard-1.130/picard.jar MarkDuplicatesWithMateCigar INPUT=$prefix.mapped_1align.sorted.bam OUTPUT=$prefix.sorted.uniq.bam METRICS_FILE=$prefix.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true SKIP_PAIRS_WITH_NO_MATE_CIGAR=true
# removes PCR duplicates
samtools index $prefix.sorted.uniq.bam


### Count reads per chromosome
echo "Counting reads per chromosome..."
samtools view $prefix.sorted.uniq.bam | cut -f 3 | sort | uniq -c | sed -e 's/^[ \t]*//' > $prefix.chromMap.txt


### Write the mapping statistics file
echo "Mapping stats..."
TOTAL_RD=$(expr $(cut -f 1 -d' ' *.R1.fastq | wc -l) / 4) # total reads
FILTERED_RD=$(expr $(cut -f 1 -d' ' *filtered_R1.fastq | wc -l) / 4) # after filtering
TRIMMED_RD=$(expr $(cut -f 1 -d' ' *trim2_R1.fastq | wc -l) / 4) # after trimming
UNMAPPED_RD=$(samtools view -S -f 4 *btout2.sam | cut -f 1 | uniq | wc -l) # no good alignments after mapping
MAPPED_RD=$(samtools view -S -F 4 *btout2.sam | cut -f 1 | uniq | wc -l) # a good alignment was found, only count one alignment per read
TOTAL_ALIGNS=$(samtools view -S -F 4 *btout2.sam | cut -f 1 | wc -l) # total number of alignments, including reads with multiple good alignments
UNIQ_ALIGNS=$(samtools view -S *mapped_1alignmentonly.sam | cut -f 1 | wc -l) # reads with only a single good alignment
UNIQ_ALIGNS_ONE_COPY=$(samtools view -S *mapped_1alignmentonly.sam | cut -f 10 | sort | uniq | wc -l) # number of reads with single alignment and their duplicate reads removed
printf "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\n" TOTAL_RD	FILTERED_RD	TRIMMED_RD	UNMAPPED_RD	MAPPED_RD	TOTAL_ALIGNS	UNIQ_ALIGNS	UNIQ_ALIGNS_ONE_COPY > $prefix.map_stats.txt
printf "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\n" $TOTAL_RD $FILTERED_RD $TRIMMED_RD $UNMAPPED_RD $MAPPED_RD $TOTAL_ALIGNS $UNIQ_ALIGNS $UNIQ_ALIGNS_ONE_COPY >> $prefix.map_stats.txt


### gzip output files
echo "Gzipping output files..."
gzip *.sam
gzip *.fastq
gzip *.bam


### TO DO: Create a QC report
# echo "QC report..."
# cp /groups/cbdm_lab/dp133/scripts/ATAC/atac* .
# Rscript atac_knit.R -w . -r $prefix.rmd -o $prefix.html
