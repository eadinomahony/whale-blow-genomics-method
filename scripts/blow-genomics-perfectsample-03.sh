### blow sample perfect sample approach (error estimation)
https://www.popgen.dk/angsd/index.php/Error_estimation

# 1. map fin whale DNA Zoo reference genome to fin whale DNA Zoo reference genome
projects/ref_genomes/Bph_DNAzoo_mappedtoMn_reference_ASM4183430v1
# fin whale is the outgroup

# map the raw reads to the reference genome; raw reads downloaded from ENA
https://www.ebi.ac.uk/ena/browser/view/SRR16970345?show=reads

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR169/045/SRR16970345/SRR16970345_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR169/045/SRR16970345/SRR16970345_2.fastq.gz


# fastp

###########################################################################
###########################################################################
#!/bin/bash

module purge
module load fastp/0.23.4

for file in SRR16970345_1.fastq.gz
do
longID=`echo ${file} | cut -f1 -d "_"`
shortID=`echo ${file} | cut -f1 -d "_"`
path=directory/path/here

fastp --in1 ${path}/${longID}_1.fastq.gz --in2 ${path}/${longID}_2.fastq.gz --out1 ${shortID}_fastp_1.fastq.gz --out2 ${shortID}_fastp_2.fastq.gz --thread 15


done

module purge



### map and calculate depth:

###########################################################################
###########################################################################
#!/bin/bash

module load bwa/0.7.17
module load samtools/1.17

#set genome ref file
nuref=GCA_041834305.1_ASM4183430v1_genomic.fna

for file in SRR16970345_fastp_1.fastq.gz

do

id=`echo ${file} | cut -f1 -d "_"`

#align to the reference genome
bwa mem -t 20 ${nuref} ${id}_fastp_1.fastq.gz ${id}_fastp_2.fastq.gz > ${id}_Bph_DNAzoo_mappedtoMn_ASM4183430v1_nuDNA_aligned.sam

#convert to bam file
samtools view -@ 20 -T ${nuref} -b ${id}_Bph_DNAzoo_mappedtoMn_ASM4183430v1_nuDNA_aligned.sam | samtools sort -@ 20 -o ${id}_Bph_DNAzoo_mappedtoMn_ASM4183430v1_nuDNA_aligned_sorted.bam

#remove duplicates
samtools rmdup -S ${id}_Bph_DNAzoo_mappedtoMn_ASM4183430v1_nuDNA_aligned_sorted.bam ${id}_Bph_DNAzoo_mappedtoMn_ASM4183430v1_nuDNA_aligned_sorted_unique.bam

#calculate depth
samtools index ${id}_Bph_DNAzoo_mappedtoMn_ASM4183430v1_nuDNA_aligned_sorted_unique.bam
STATS=$(samtools depth -a ${id}_Bph_DNAzoo_mappedtoMn_ASM4183430v1_nuDNA_aligned_sorted_unique.bam -@ 5 | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}')
echo -e "${id}\t${STATS}" >> mapped_${id}_Bph_DNAzoo_coverage.txt
echo "finished ${id}"

done


#### repeat mask outgroup

###########################################################################
###########################################################################
#!/bin/bash

module purge 

module load samtools/1.21
module load bedtools/2.31.0


INDIR=directory/path/here
OUTDIR=directory/path/here
BED=GCA_041834305.1_ASM4183430v1_repeats.bed
BAM=SRR16970345_Bph_DNAzoo_mappedtoMn_ASM4183430v1_nuDNA_aligned_sorted_unique.bam

shortid=`echo $BAM | cut -f1,2,3,4 -d "_"`

bedtools intersect -abam ${INDIR}/${BAM} -b ${BED} -v > ${OUTDIR}/${shortid}_fastp_nuDNA_aligned_sorted_unique_RepeatMasker.bam # -v flag REMOVES rather than keeps the regions that are highlighted in the mask .gff file

samtools index ${OUTDIR}/${shortid}_fastp_nuDNA_aligned_sorted_unique_RepeatMasker.bam

module purge




# 2. consensus of the perfect sample in fasta format (-ref below)

# simply choose the highest coverage sample:
# Humpback DNA Zoo reference genome:
SRR17854487


###########################################################################
###########################################################################
#!/bin/bash

module purge
module load angsd/0.940-2


output=directory/path/here
autosomes=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt
perfectbam=SRR17854487_fastp_nuDNA_aligned_sorted_unique_RepeatMasker.bam

# generate consensus fasta from your perfect individual
angsd -i $perfectbam -doCounts 1 -nThreads 20 -doFasta 2 -minMapQ 30 -minQ 30 -out ${output}/SRR17854487_humpback_perfectindiv_consensus

module purge


###########################################################################
###########################################################################
#!/bin/bash

module purge
module load angsd/0.940

output=directory/path/here
autosomes=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt
perfectbam=SRR17854487_fastp_nuDNA_aligned_sorted_unique_RepeatMasker.bam
bamlist=perfectindiv_bamlist.txt # list of all bam files for error rate calculation
outgroupbam=SRR16970345_Bph_DNAzoo_mappedtoMn_fastp_nuDNA_aligned_sorted_unique_RepeatMasker.bam


# generate consensus fasta from the outgroup (DNA Zoo Bph mapped to Mn ref)
# -doFasta 2 uses the most common base
angsd -i $outgroupbam -doFasta 2 -doCounts 1 -nThreads 20 -minMapQ 30 -minQ 30 -out $output/SRR16970345_Bph_DNAzoo_mappedtoMn_outgroup_consensus
# will need to gunzip the output before the next step


gunzip SRR17854487_humpback_perfectindiv_consensus.fa.gz # perfect indiv
gunzip SRR16970345_Bph_DNAzoo_mappedtoMn_outgroup_consensus.fa.gz # outgroup


perf=SRR17854487_humpback_perfectindiv_consensus.fa # fasta file of consensus of perfect indiv
outgroup=SRR16970345_Bph_DNAzoo_mappedtoMn_outgroup_consensus.fa # containing outgroup fasta consensus

# error estimation using outgroup + perfect individual
angsd -doAncError 1 -anc ${outgroup} -ref ${perf} -nThreads 20 -out ${output}/methodspaper_blow_MnMap_BphOutgroup_error -bam $bamlist -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -rf ${autosomes}

# -ref = should not be your reference, but a fasta file containing the consensus of a 'perfect' sample 

module purge



###########################################################################
###########################################################################


wget https://raw.githubusercontent.com/ANGSD/angsd/refs/heads/master/R/estError.R

### id file can be any sample IDs but has to be same order as provided bam list above

awk -F'/' '{split($NF,a,"_fastp"); print a[1]}' input.txt > output.txt
awk -F'/' '{split($NF,a,"_fastp"); print a[1]}' perfectindiv_bamlist.txt > error.ids.txt
sed -i 's/\r$//' error.ids.txt
cat -v error.ids.txt | head


#!/bin/bash

module load gcc/11.2.0 R/4.3.1

Rscript estError.R file=methodspaper_blow_MnMap_BphOutgroup_error.ancError indNames=error.ids.txt

