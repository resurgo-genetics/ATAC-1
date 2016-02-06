#!/bin/bash
# run where the .sorted.uniq.bam.gz files are
# usage: ./atac_peaks_PE.sh  dir_where_files_are prefix
# 
#
# author: Joanna R. DiSpirito
# last modified: 2/4/16


#### Load modules
source ~/.bash_profile

module load seq/samtools/0.1.19
module load seq/BEDtools/2.19.0
module load seq/homer/4.6
module load seq/UCSC-tools

#### Set working directory and prefix for output files (args)
echo "working directory : " $1
echo "prefix : " $2
prefix=$2

cd $1

#### Gunzip files
echo "Decompressing files..."
gunzip -fv *sorted.uniq.bam.gz


#### makeTagDirectory and generate QC files
echo "Making Tag Directory..."
makeTagDirectory $prefix.homer.tagdir/ -genome /groups/shared_databases/genomes/mm9.fa -tbp 1 -checkGC -illuminaPE *sorted.uniq.bam &> $prefix.tagdir.stdout.txt


#### make bedGraph file
echo "Making bedGraph file..."
makeUCSCfile $prefix.homer.tagdir/ -fsize 1e50 -norm 1e6 -name $prefix -o $prefix.bedGraph


#### make bigWig file
echo "Making biwig file..."

# move chrom sizes from jrd/resources directory to current directory
cp /groups/cbdm_lab/jrd26/JRD_Resources/mm9.chrom.sizes .

# unzip bedGraph
gunzip -fv $prefix.bedGraph.gz

# sort bedGraph and delete last line
sort -k1,1 -k2,2n $prefix.bedGraph >  $prefix_sorted1.bedGraph
sed '$d'  $prefix_sorted1.bedGraph  > $prefix_sorted.bedGraph

# convert to bw
bedGraphToBigWig $prefix_sorted.bedGraph mm9.chrom.sizes $prefix.bw


#### Find peaks (altered original script to add different settings)
echo "Finding peaks in factor mode..."
findPeaks $prefix.homer.tagdir/ -style factor -norm 1000000 -o $prefix.peaks.factor.txt &> $prefix.stdout.peaks.factor.txt

echo "Finding peaks in histone mode..."
findPeaks $prefix.homer.tagdir/ -style histone -norm 1000000 -o $prefix.peaks.histone.txt &> $prefix.stdout.peaks.histone.txt

echo "Finding peaks in dnase mode..."
findPeaks $prefix.homer.tagdir/ -style dnase -norm 1000000 -o $prefix.peaks.dnase.txt &> $prefix.stdout.peaks.dnase.txt

echo "Merging results of dnase and factor modes..."
mergePeaks $prefix.peaks.dnase.txt $prefix.peaks.factor.txt -o $prefix.peaks.atac.txt &> $prefix.stdout.merge.txt
 

#### Compress files
echo "zipping bedGraph file..."
gzip $prefix.bedGraph

echo "zipping bigWig file..."
gzip $prefix.bw

echo "zipping bam file..."
gzip *sorted.uniq.bam


########################### DEPRECATED ######################################################
# deprecated Amit lab settings
# findPeaks homer_tagdir/ -L 0 -C 3 -size 1000 -minDist 1000 -tbp 1 -o peaks.txt

# deprecated DPZ settings
# echo "Finding peaks using user settings..."
# findPeaks $prefix.homer.tagdir/ -C 3 -norm 1000000 -o $prefix.peaks.user.txt &> $prefix.stdout.peaks.user.txt