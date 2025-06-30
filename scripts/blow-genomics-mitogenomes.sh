#### DOWNLOADING MITOGENOME DATA FROM NCBI

# 1. NC_006927.1 Megaptera novaeangliae mitochondrion, complete genome:
# https://www.ncbi.nlm.nih.gov/nuccore/NC_006927.1?report=fasta

megaptera.novaeangliae.mitogenome.NC_006927.1.AntarcticOcean.fna


# 2. PP475430.1 Megaptera novaeangliae voucher NIST KW2013002 mitochondrion, complete genome
# https://www.ncbi.nlm.nih.gov/nuccore/PP475430.1?report=fasta

megaptera.novaeangliae.mitogenome.PP475430.1.Hawaii.fna

# first 70 bases
GTTAATGTAGCTTAAACACTCACAAAGCAAGACACTGAAAATGTCTAGATGGGTCTAATCAACCCCATTG

# last 70 bases
TTATAAATCAATACTAAATCTGACACAAGCCCAATAATGAAAATACATGAACCCTATCCCTATCCAATAC

# 3. AP006467.1 Megaptera novaeangliae mitochondrial DNA, complete genome
# https://www.ncbi.nlm.nih.gov/nuccore/AP006467.1

megaptera.novaeangliae.mitogenome.AP006467.1.AntarcticOcean.fna

# 4. MF409246.1 Megaptera novaeangliae mitochondrion, complete genome
# https://www.ncbi.nlm.nih.gov/nuccore/MF409246.1

megaptera.novaeangliae.mitogenome.MF409246.1.Arnason.fna

# 5.FJ590425.1 Megaptera novaeangliae mitochondrion, partial genome (Carraher et al. "Efficient Mitogenomic Sequencing of Cetaceans") - New Zealand - 12,571 bp
# https://www.ncbi.nlm.nih.gov/nuccore/FJ590425.1 
megaptera.novaeangliae.mitogenome.partial.FJ590425.1.NewZealand.fna


# indexing Hawai'i mitogenome

#!/bin/bash

module load samtools/1.17
module load bwa/0.7.17

for file in *.fna # this might need to be changed to .fna for ref genomes from NCBI
do 

ref=${file}

samtools faidx ${ref}

bwa index ${ref}

done

### mapping samples to mitogenome from Hawaii

#!/bin/bash


module purge

module load bwa/0.7.17
module load samtools/1.17

#set genome ref file

nuref=megaptera.novaeangliae.mitogenome.PP475430.1.Hawaii.fna

for dir in $(cat sample.directories.txt)

do

cd ${dir}

for file in $(cat fwd_fileNames.txt) # insert .txt file here instead of: *clean_artifacts_trimmed_1.fq.gz

do

id=`echo ${file} | cut -f1,2,3,4 -d "_"`

#align to the reference genome
bwa mem -t 20 ${nuref} ${id}_1.fq.gz ${id}_2.fq.gz > ${id}_mtDNA_aligned.sam

#convert to bam file
samtools view -@ 20 -T ${nuref} -b ${id}_mtDNA_aligned.sam | samtools sort -@ 20 -o ${id}_mtDNA_aligned_sorted.bam

#remove duplicates
samtools rmdup -S ${id}_mtDNA_aligned_sorted.bam ${id}_mtDNA_aligned_sorted_unique.bam

samtools index ${id}_mtDNA_aligned_sorted_unique.bam
STATS=$(samtools depth -a ${id}_mtDNA_aligned_sorted_unique.bam -@ 5 | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}')
echo -e "${id}\t${STATS}" >> mapped_mtDNA_coverage.txt

done

cd ..

done

module purge


## Merge mtDNA bam files across sequencing runs:

#!/bin/bash


module purge
module load samtools/1.17

for dir in $(cat sample.directories.txt)

do 

cd ${dir} 

ls *_fastp_mtDNA_aligned_sorted_unique.bam > bam.files.to.merge.mtDNA.txt

samtools merge -o ${dir}_merged_fastp_mtDNA_aligned_sorted_unique.bam -b bam.files.to.merge.mtDNA.txt -@ 5

samtools index ${dir}_merged_fastp_mtDNA_aligned_sorted_unique.bam

cd ..

done

### calc depth on merged runs:

vim depth.mtDNA.sh

#!/bin/bash

module load samtools/1.17

for dir in $(cat sample.directories.txt)

do 

cd ${dir} 

# calculate mean and standard deviation of bam coverage
samtools index ${dir}_merged_fastp_mtDNA_aligned_sorted_unique.bam
STATS=$(samtools depth -a ${dir}_merged_fastp_mtDNA_aligned_sorted_unique.bam -@ 5 | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}')
echo -e "${dir}\t${STATS}" >> mapped_mtDNA_merged_coverage.txt
echo "finished ${dir}"

cd ..

done



### map stats for mtDNA: 

#!/bin/bash

module purge
module load samtools/1.17

for dir in $(cat sample.directories.txt)
do

cd ${dir}

ls *_fastp_mtDNA_aligned_sorted.bam | cut -f1,2,3 -d "_" > map.stats.files.mtDNA.txt

# counting only mapped (primary aligned) reads
for file in $(cat map.stats.files.mtDNA.txt)
do
echo ${file}_fastp_mtDNA_aligned_sorted_unique.bam >> mapped.reads.mtDNA.txt
samtools view -c -F 260 ${file}_fastp_mtDNA_aligned_sorted_unique.bam >> mapped.reads.mtDNA.txt 
done

cd .. 

done

module purge





### consensus mitogenomes

#!/bin/bash


module purge
module load angsd/0.940-2

for file in *_fastp_mtDNA_aligned_sorted_unique.bam
do
nuref=megaptera.novaeangliae.mitogenome.PP475430.1.Hawaii.fna
bams=global_mitogenome_bams.txt
#shortID=`basename ${bams} | cut -d "/" -f8 | cut -d "_" -f1`
shortID=`basename ${file} | cut -d "_" -f1`

angsd -i ${file} -ref ${nuref} -doFasta 2 -doCounts 1 -out /consensus_mitogenomes/${shortID}_mitogenome_consensus
done

module purge


### be careful because the header of all of the consensus mitogenomes will be the name of the reference (in this case PP475430.1)
zcat MN-b-001_mitogenome_consensus.fa.gz | head

### so need to append the header with the sample name, keeping what its been mapped to for later reference: such as 
MN-b-001_mappedto_PP475430.1

sed -i 's/^>PP475430.1/>MN-b-001_mappedto_PP475430.1/' MN-b-001_mitogenome_consensus.fa


### now loop through all files to do this to all of them:
gunzip *.fa.gz

for file in *.fa
do
  # Extract the base sample name from the filename (e.g., MN-b-001)
  sample=$(basename "$file" _mitogenome_consensus.fa)

  # Use sed to modify the header inside each .fa file
  sed -i "s/^>/>${sample}_mappedto_/" "$file"

done


sed -i 's/^>PP475430.1/>MN-b-007_mappedto_PP475430.1/' MN-b-007_mitogenome_consensus.fa



###### merge bam files for replicate blow samples


#!/bin/bash

module purge
module load samtools/1.17

# samtools merge -o output.bam input1.bam input2.bam ... -@ threads

id=BCX0121
s1=MN-b-010
s2=MN-b-109
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCX0171
s1=MN-b-035
s2=MN-b-043
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCX0427
s1=MN-b-002
s2=MN-b-090
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam


id=BCX0694
s1=MN-b-013 # low coverage
s2=MN-b-081
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCX1225
s1=MN-b-045
s2=MN-b-103
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCX1393
s1=MN-b-018
s2=MN-b-096
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCX1484
s1=MN-b-023
s2=MN-b-079
s3=MN-b-095
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam ${s3}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCX0427 # intentionally leaving out MN-b-004 because it has low KING values for nuDNA with all other samples
s1=MN-b-002
s2=MN-b-090
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam


id=BCX1552
s1=MN-b-020
s2=MN-b-127
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCX1596
s1=MN-b-066
s2=MN-b-087
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCX1712
s1=MN-b-007
s2=MN-b-078
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCX2052
s1=MN-b-015
s2=MN-b-082
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCY0117
s1=MN-b-009
s2=MN-b-093
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam


id=BCY0474
s1=MN-b-084
s2=MN-b-085
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCY0863
s1=MN-b-048
s2=MN-b-104
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCY0910
s1=MN-b-055
s2=MN-b-088
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCY0918
s1=MN-b-036
s2=MN-b-092
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCY0965
s1=MN-b-003
s2=MN-b-083
s3=MN-b-099
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam ${s3}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCZ0195
s1=MN-b-008
s2=MN-b-091
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCZ0254
s1=MN-b-040
s2=MN-b-069
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCZ0333
s1=MN-b-111
s2=MN-b-115
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCZ0442
s1=MN-b-063
s2=MN-b-106
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=BCZ0482
s1=MN-b-064
s2=MN-b-126
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=CSX0213
s1=MN-b-028
s2=MN-b-080
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

id=CSY0149
s1=MN-b-073
s2=MN-b-074
S3=MN-b-117
samtools merge -o merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam ${s1}_fastp_mtDNA_aligned_sorted_unique.bam ${s2}_fastp_mtDNA_aligned_sorted_unique.bam ${s3}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5
samtools index merged_mitogenomes/${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam

module purge

### calc depth on merged runs:
MN-b-001-K_merged_fastp_mtDNA_aligned_sorted_unique.bam

vim depth.mtDNA.sh

#!/bin/bash

module load samtools/1.17

for id in $(cat merged_mtDNA_bam_list_IDs.txt)

do 

# calculate mean and standard deviation of bam coverage
samtools index ${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam
STATS=$(samtools depth -a ${id}_merged_fastp_mtDNA_aligned_sorted_unique.bam -@ 5 | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}')
echo -e "${id}\t${STATS}" >> mapped_mtDNA_merged_coverage.txt
echo "finished ${id}"

done


### depth for unmerged data

ls *.bam | cut -f1 -d "_" > unmerged_mtDNA_bam_list_IDs.tx

#!/bin/bash

module load samtools/1.17

for id in $(cat unmerged_mtDNA_bam_list_IDs.txt)

do 

# calculate mean and standard deviation of bam coverage
samtools index ${id}_fastp_mtDNA_aligned_sorted_unique.bam
STATS=$(samtools depth -a ${id}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5 | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}')
echo -e "${id}\t${STATS}" >> mapped_mtDNA_unmerged_coverage.txt
echo "finished ${id}"

done


### flagstat for endogenous content


#!/bin/bash

module purge
module load samtools/1.17

for dir in $(cat sample.directories.endogenous.txt)

do 

cd ${dir} 

ls *_fastp_mtDNA_aligned_sorted.bam > bam.files.to.merge.mtDNA.endogenous.txt

samtools merge -o ${dir}_merged_fastp_mtDNA_aligned_sorted.bam -b bam.files.to.merge.mtDNA.endogenous.txt -@ 5

samtools index ${dir}_merged_fastp_mtDNA_aligned_sorted.bam

samtools flagstat ${dir}_merged_fastp_mtDNA_aligned_sorted.bam -@ 5 > ${dir}.mtDNA.flagstat_output.txt

cd ..

done



##### pull the total reads and mapped reads out for these

echo -e "Filename\tTotal_Reads\tMapped_Reads" > mtDNA_mapped_reads_summary.txt

for file in *.mtDNA.flagstat_output.txt; do
    total=$(awk 'NR==1 {print $1}' "$file")
    mapped=$(awk '/ mapped / {print $1; exit}' "$file")
    echo -e "$file\t$total\t$mapped" >> mtDNA_mapped_reads_summary.txt
done





########### suarez mitogenomes:


#!/bin/bash


module purge

module load bwa/0.7.17
module load samtools/1.17

#set genome ref file

nuref=megaptera.novaeangliae.mitogenome.PP475430.1.Hawaii.fna

for file in $(cat SRR_accession_keys3_menendez.txt) # insert .txt file here instead of: *clean_artifacts_trimmed_1.fq.gz

do

# SRR25114440_1.fastq.gz
#id=`echo ${file} | cut -f1 -d "_"`
id=$file

#align to the reference genome
bwa mem -t 20 ${nuref} ${id}_fastp_R1.fastq.gz ${id}_fastp_R2.fastq.gz > ${id}_fastp_mtDNA_aligned.sam

#convert to bam file
samtools view -@ 20 -T ${nuref} -b ${id}_fastp_mtDNA_aligned.sam | samtools sort -@ 20 -o ${id}_fastp_mtDNA_aligned_sorted.bam

#remove duplicates
samtools rmdup -S ${id}_fastp_mtDNA_aligned_sorted.bam ${id}_fastp_mtDNA_aligned_sorted_unique.bam

samtools index ${id}_fastp_mtDNA_aligned_sorted_unique.bam
STATS=$(samtools depth -a ${id}_fastp_mtDNA_aligned_sorted_unique.bam -@ 5 | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}')
echo -e "${id}\t${STATS}" >> mapped_mtDNA_fastp_coverage.txt

done

module purge



#### call consensus on these:


#!/bin/bash


module purge
module load angsd/0.940

for file in *_mtDNA_aligned_sorted_unique.bam
do
nuref=megaptera.novaeangliae.mitogenome.PP475430.1.Hawaii.fna

#shortID=`basename ${bams} | cut -d "/" -f8 | cut -d "_" -f1`
shortID=`basename ${file} | cut -d "_" -f1`

angsd -i ${file} -ref ${nuref} -doFasta 2 -doCounts 1 -out ${shortID}_mitogenome_consensus

done

module purge


### be careful because the header of all of the consensus mitogenomes will be the name of the reference (in this case PP475430.1)
zcat SRR25114445_mitogenome_consensus.fa.gz | head

### so need to append the header with the sample name, maybe keeping what its been mapped to for later reference: such as 
SRR25114445_mappedto_PP475430.1

sed -i 's/^>PP475430.1/>SRR25114445_mappedto_PP475430.1/' SRR25114445_mitogenome_consensus.fa


### now loop through all files to do this to all of them:
gunzip *.fa.gz

for file in *.fa
do
  # Extract the base sample name from the filename (e.g., MN-b-001)
  sample=$(basename "$file" _mitogenome_consensus.fa)
  echo $sample
  # Use sed to modify the header inside each .fa file
  sed -i "s/^>/>${sample}_mappedto_/" "$file"

done