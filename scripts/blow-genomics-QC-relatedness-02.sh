#### repeat mask

#!/bin/bash

module purge 

module load samtools/1.21
module load bedtools/2.31.0

for BAM in $(cat SRR_accession_keys.txt) # list all of these Bph00161_nuDNA_aligned_sorted_unique.bam with absolute path
do

INDIR=directory/path/here
OUTDIR=directory/path/here
BED=/humpback_NCBI_ASM4183430v1/GCA_041834305.1_ASM4183430v1_repeats.bed
#shortid=`echo $BAM | cut -f1 -d "_"`

bedtools intersect -abam ${INDIR}/${BAM}_ASM4183430v1_nuDNA_aligned_sorted_unique.bam -b ${BED} -v > ${OUTDIR}/${BAM}_fastp_nuDNA_aligned_sorted_unique_RepeatMasker.bam # -v flag REMOVES rather than keeps the regions that are highlighted in the mask .gff file

samtools index ${OUTDIR}/${BAM}_fastp_nuDNA_aligned_sorted_unique_RepeatMasker.bam

done

module purge


#### replatedness

#!/bin/bash


module purge
module load ngsrelate/2.0 

input=ngsrelate.refined.noMinDepth.prep.glf.gz
samples=sorted_bam_file_paths_refined_filteredbydepth.txt
num_inds=58

# count the # of sites:
zcat ngsrelate.refined.noMinDepth.prep.mafs.gz | tail -n +2 | wc -l > sites_count_refined_noMinDepth.txt

# extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat ngsrelate.refined.noMinDepth.prep.mafs.gz | cut -f5 |sed 1d > freq_noMinDepth 

ngsRelate -g ${input} -n ${num_inds} -f freq_noMinDepth -z ${samples} -p 10 -O ngsrelate.refined.run.noMinDepth.out

module purge