#!/bin/bash
# run where the fastq.gz files are
# may have to make executable with chmod +x atac_mapping_SR.sh
# usage: ./atac_mapping_PE_to_SR.sh  dir_where_fastq.gz_files_are prefix

#### Load libraries and modules
source ~/.bash_profile
# export R_LIBS="/groups/cbdm_lab/dp133/R_libraries"

module load seq/fastqc/0.10.1
module load seq/cutadapt/1.8.3 
module load seq/fastx/0.0.13
module load seq/bowtie/2.2.4
module load seq/samtools/0.1.19
module load seq/BEDTools/2.19.0
module load seq/homer/4.6
module load stats/R/3.1.2
module load seq/sickle/1.2

#### Set working directory and prefix for output files (args)
echo "working directory : " $1
echo "prefix : " $2
prefix=$2

cd $1

#### Gunzip files
echo "Decompressing files..."
gunzip -fv *.fastq.gz
mv *.fastq $prefix.fastq

 
#### Generate fastqc output
echo "Running fastqc..."
mkdir $prefix.fastqc/
fastqc -o $prefix.fastqc/ $prefix.fastq


#### Filter reads on quality
echo "Filtering and trimming reads on quality..."
sickle se -t sanger -f $prefix.fastq -o $prefix.filtered.fastq


#### Clip adapters from reads
# e = maximum error rate, default = 0.1
# m = minimum length, throw away reads shorter than N bases
echo "Clipping adapter from 5' side..."
cutadapt -g AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT -e 0.1 -m 20 -o $prefix.trim1.fastq $prefix.filtered.fastq

echo "Clipping adapter from 3' side..."
cutadapt -a AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -e 0.1 -m 20 -o $prefix.trim2.fastq $prefix.trim1.fastq


#### Align reads 
# -p N		number of threads to use for multi-thread software
# -x		location of indexed genome
# -S <>		name of .sam output file
echo "Mapping reads to mm9..."
bowtie2 -p 8 -x /groups/shared_databases/bowtie2_indexes/mm9 -X 1000 $prefix.trim2.fastq -S $prefix.btout2.sam


#### Keep only reads mapping to one location
samtools view -hS -F 4 $prefix.btout2.sam > $prefix.mapped.sam # keep only mapped reads
sed '/XS:/d' $prefix.mapped.sam > $prefix.mapped_1alignmentonly.sam # remove multiple aligned reads
samtools view -bS $prefix.mapped_1alignmentonly.sam > $prefix.mapped_1align.bam # make a bam file


### Remove duplicates 
echo "Removing duplicates..."
# sort sam file
samtools sort $prefix.mapped_1align.bam $prefix.mapped_1align.sorted
# remove duplicates
# REMOVE_DUPLICATES = if true do not write duplicates to the output file,instead of writing them with the appropriate flags; default = false
# ASSUME_SORTED = if true, assume the input file is coordinate sorted even if header says otherwise; default = false
java -Xms1024m -jar /opt/picard-1.130/picard.jar MarkDuplicates INPUT=$prefix.mapped_1align.sorted.bam OUTPUT=$prefix.sorted.dedup.bam METRICS_FILE=$prefix.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
# create bam index
samtools index $prefix.sorted.dedup.bam


### Count reads per chromosome
echo "Counting reads per chromosome..."
samtools view $prefix.sorted.dedup.bam | cut -f 3 | sort | uniq -c | sed -e 's/^[\t]*//' > $prefix.chromMap.txt


### Write the mapping statistics file
echo "Mapping stats..."
TOTAL_RD=$(expr $(cut -f 1 -d' ' *.fastq | wc -l) / 4) # total reads
FILTERED_RD=$(expr $(cut -f 1 -d' ' *filtered.fastq | wc -l) / 4) # after filtering
TRIMMED_RD=$(expr $(cut -f 1 -d' ' *trim2.fastq | wc -l) / 4) # after trimming
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


### Create a QC report
# To Do: create .Rmd





























