##################################################################################################
##################################################################################################
#### index reference genome

#!/bin/bash
module load samtools/1.17
module load bwa/0.7.17

for file in *.fna # this might need to be changed to .fna for ref genomes from NCBI
do 

ref=${file}

samtools faidx ${ref}

bwa index ${ref}

done




##################################################################################################
##################################################################################################
#### md5sum checks

#!/bin/bash
# Bash script to calculate MD5 checksums for all .fastq files in a directory

# Directory containing .fastq files
dir=directory/path/here

# Log file to store MD5 checksum results
log="md5.checksums.log"

# Navigate to the directory (optional if running script from within the directory)
cd "$dir" || exit

# Ensure the log file is empty or create it if it doesn't exist
> "$log"

# Loop through each .fastq file
for file in *.fastq; do
    # Calculate MD5 checksum
    md5sum "$file" >> "$log"
done



##################################################################################################
##################################################################################################
#### count raw reads

#!/bin/bash

for dir in $(cat dirs.txt)
do
cd $dir
READS_fwd=$(zcat ${dir}_1.fq.gz | awk 'END {print NR/4}')
READS_rev=$(zcat ${dir}_2.fq.gz | awk 'END {print NR/4}')
total=$(awk -v f=$READS_fwd -v r=$READS_rev 'BEGIN {print f + r}')
echo -e "${dir}\t${READS_fwd}\t${READS_rev}\t${total}" >> ../raw.reads.count.txt
echo "finished $dir"
cd ..
done


##################################################################################################
##################################################################################################
#### fastp clean up of raw data

#!/bin/bash
module purge
module load fastp/0.23.4

for file in $(cat sorted_file_names.txt) # make this new
do
longID=`echo ${file} | cut -f1,2,3 -d "_"`
#shortID=`echo ${file} | cut -f1,2,3 -d "_" | cut -f4 -d "."`
path=path/to/files # update this

fastp --in1 ${path}/${longID}_1.fq.gz --in2 ${path}/${longID}_2.fq.gz --out1 ${longID}_fastp_1.fq.gz --out2 ${longID}_fastp_2.fq.gz --thread 15

#echo "${path}/${longID}_1.fq.gz"

done

module purge




##################################################################################################
##################################################################################################
#### mapping script

#!/bin/bash
module purge
module load bwa/0.7.17
module load samtools/1.17

#set genome ref file
nuref=GCA_041834305.1_ASM4183430v1_genomic.fna


for file in *1.fq.gz

do

id=`echo ${file} | cut -f1,2,3 -d "_"`

#align to the reference genome
bwa mem -t 20 ${nuref} ${id}_1.fq.gz ${id}_2.fq.gz > ${id}_nuDNA_aligned.sam

#convert to bam file
samtools view -@ 20 -T ${nuref} -b ${id}_nuDNA_aligned.sam | samtools sort -@ 20 -o ${id}_nuDNA_aligned_sorted.bam

#remove duplicates
samtools rmdup -S ${id}_nuDNA_aligned_sorted.bam ${id}_nuDNA_aligned_sorted_unique.bam

samtools index ${id}_nuDNA_aligned_sorted_unique.bam
STATS=$(samtools depth -a ${id}_nuDNA_aligned_sorted_unique.bam -@ 5 | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}')
echo -e "${id}\t${STATS}" >> mapped_blow_coverage.txt

done

module purge