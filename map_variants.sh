#!/bin/bash

#extract contigs from a given pool
# ids=$1
# origin=$2

# perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $ids $origin

#filter the reads fastq file
# awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 21 && length(seq) <= 25) {print header, seq, qheader, qseq}}' < your.fastq > filtered.fastq


echo "This program maps reads onto a reference.\n"
echo "./map_variants Reference reads k s\n"
#time to start
STARTTIME=$(date +%s)

#read the reference and read files from command line
REF=$1
READS=$2
k=$3
s=$4

#index the reference using SMALT.
echo "indexing $REF as reference file"
smalt index -k $k -s $s index.$k.$s $1 #this will index and rename the index based on the values supplied at the commandline

#run SMALT to map the reads against the conntigs
echo "running SMALT with k=$k and s=$s parameters"
smalt map -f sam -o $REF.sam index.$k.$s $READS

#convert the sam file to bam file
echo "converting $REF.sam to corresponding BAM "
samtools view -bS $REF.sam >$REF.bam

#sort the bam file
echo "sorting $REF.bam"
samtools sort $REF.bam $REF.sorted

#index the bam file
echo "indexing $REF.sorted.bam"
samtools index $REF.sorted.bam $REF.sorted.bai

#filter reads out?
echo "filtering reads with lengths >= 300 and a mapping quality of >=30"
sambamba view -f bam --filter='sequence_length >= 300 and mapping_quality >=30' $REF.sorted.bam >$REF.l300.q30.bam

#sort  and index the filtered bam file
samtools sort $REF.l300.q30.bam $REF.l300.q30.sorted
samtools index $REF.l300.q30.sorted.bam $REF.l300.q30.sorted.bai

#get a mpileup bcf file
echo "Calling variants with samtools mpileup"
samtools mpileup -u -f $REF $REF.l300.q30.sorted.bam >$REF.l300.q30.sorted.bcf

#convert bcf to human readable vcf
echo "converting BCF to VCF format "
bcftools view -vcg $REF.l300.q30.sorted.bcf >$REF.l300.q30.vcf

#create a bed file
echo "creating a bed file for variants"
sed -e 's/chr//' $REF.l300.q30.vcf | awk '{OFS="\t"; if (!/^#/) {print $1,$2-1,$2,$4"/"$5,$6}}' >$REF.l300.q30.bed

echo "pipeline finished"

ENDTIME=$(date +%s)

echo "It took $(($ENDTIME - $STARTTIME)) seconds to run this job!"
