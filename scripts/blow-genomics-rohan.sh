### ROHan is a Bayesian framework to estimate local rates of heterozygosity, infer runs of homozygosity (ROH) and compute global rates of heterozygosity.
# https://github.com/grenaud/rohan

# Due to the lack of heterozygosous sites, ROHs can cause an underestimate of global estimates of heterozygosity. 
# Such global estimates of rates of heterozygous sites is useful to infer population genetics parameters such as effective population size (Ne).

# Do not apply any filters as mapping quality and base quality are all informative in the model. 
# Duplicate removal is very recommended. Simply use the fail QC flag to remove reads/fragments that have failed basic quality control (e.g. for duplicates). 
# Simply sort and index and provide ROHan the same reference used for mapping.

### using bam files from BEFORE repeat masking is advised -- those ending in
nuDNA_aligned_sorted_unique.bam



###########################################################################
###########################################################################

#### ROHan on merged data:

#!/bin/bash

module purge
module load rohan/20230903
module load samtools/1.17

ref=GCA_041834305.1_ASM4183430v1_genomic.fna
chr=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt
outpath=/mergedData_analyses/rohan
bamfolder=/methods_paper_bams/merged_bams

for file in $(cat merged_bam_list.txt)

do

id=`echo ${file} | cut -f4 -d "_" | cut -f2 -d "/"`

samtools view -hb -F 3844 -@ 20 ${file} -o ${outpath}/${id}_mergedblow_filtered_rohan.bam
samtools index ${outpath}/${id}_mergedblow_filtered_rohan.bam

rohan --rohmu 2e-5 -o ${id}_mergedblow_rohan -t 40 ${ref} ${id}_mergedblow_filtered_rohan.bam

done

module purge


###########################################################################
###########################################################################

### ROHan public data:
ls -d -1 "$PWD/"*_fastp_nuDNA_aligned_sorted_unique.bam > public_bam_list.txt

#!/bin/bash

module purge
module load rohan/20230903
module load samtools/1.17

ref=GCA_041834305.1_ASM4183430v1_genomic.fna
chr=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt
outpath=/mergedData_analyses/rohan

for file in $(cat public_bam_list.txt)

do

id=`echo ${file} | cut -f2 -d "_" | cut -f2 -d "/"`

samtools view -hb -F 3844 -@ 20 ${file} -o ${outpath}/${id}_public_filtered_rohan.bam
samtools index ${outpath}/${id}_public_filtered_rohan.bam

rohan --rohmu 2e-5 -o ${id}_public_rohan -t 40 ${ref} ${id}_public_filtered_rohan.bam

done

module purge



###########################################################################
###########################################################################

#### ROHan DFO biopsies


ls -d -1 "$PWD/"*_fastp_nuDNA_aligned_sorted_unique.bam > dfo_biopsy_bam_list.txt

#!/bin/bash

module purge
module load rohan/20230903
module load samtools/1.17

ref=GCA_041834305.1_ASM4183430v1_genomic.fna
chr=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt
outpath=/mergedData_analyses/rohan


for file in $(cat dfo_biopsy_bam_list.txt)

do

id=`echo ${file} | cut -f4 -d "_" | cut -f2 -d "/"`

samtools view -hb -F 3844 -@ 20 ${file} -o ${outpath}/${id}_dfo_biopsy_filtered_rohan.bam
samtools index ${outpath}/${id}_dfo_biopsy_filtered_rohan.bam

rohan --rohmu 2e-5 -o ${id}_dfo_biopsy_rohan -t 40 ${ref} ${id}_dfo_biopsy_filtered_rohan.bam

done

module purge

###########################################################################
###########################################################################

##### ROHan on unmerged data and subsampled biopsies:

ls -d -1 "$PWD/"*_fastp_nuDNA_aligned_sorted_unique.bam > unmerged_bam_list.txt

#!/bin/bash

module purge
module load rohan # rohan/20230903
module load samtools/1.21

ref=GCA_041834305.1_ASM4183430v1_genomic.fna
chr=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt
outpath=/methods_paper/rohan

for file in $(cat unmerged_bam_list.txt)

do

id=`echo ${file} | cut -f3 -d "_" | cut -f2 -d "/"`

samtools view -hb -F 3844 -@ 20 ${file} -o ${outpath}/${id}_unmergedblow_filtered_rohan.bam
samtools index ${outpath}/${id}_unmergedblow_filtered_rohan.bam

rohan --rohmu 2e-5 -o ${id}_unmergedblow_rohan -t 40 ${ref} ${id}_unmergedblow_filtered_rohan.bam

done

module purge


#### subsampled biopsies

#!/bin/bash

module purge
module load rohan/20230903
module load samtools/1.17

ref=GCA_041834305.1_ASM4183430v1_genomic.fna
chr=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt
outpath=/methods_paper/rohan

for file in $(cat subsampled_biopsies_bam_list.txt)

do

id=`echo ${file} | cut -f9 -d "/" | cut -f1,2,3 -d "_"`

samtools view -hb -F 3844 -@ 20 ${file} -o ${outpath}/${id}_biop_filtered_rohan.bam
samtools index ${outpath}/${id}_biop_filtered_rohan.bam

rohan --rohmu 2e-5 -o ${id}_biop_rohan -t 40 ${ref} ${id}_biop_filtered_rohan.bam

done

module purge




#Blow sample plot 
#Use the files looking like this: SRR25114440_public_rohan.mid.hmmrohl

#Skips the first line, grabs the length of the ROH and prints the sample ID in the second column

for file in *mid.hmmrohl
do
shortID=$(echo "${file}" | cut -f1,2,3 -d "_")
awk -v ID="${shortID}" 'NR > 1 {print $5, $6, ID}' "${file}" >> roh.summary.txt
done

#Then added this header afterwards
ROH_LENGTH VALIDATED_SITES ID


#### extract info for unmerged data to summarise in ms:

output="roh_lines8_11_extracted.txt"
echo -e "Sample\tLine8\tLine11" > $output

for file in *_unmergedblow_rohan.summary.txt; do
  sample=$(basename "$file" | cut -d'_' -f1)
  
  # Extract lines 8 and 11
  line8=$(sed -n '8p' "$file" | tr -d '\n' | sed 's/\s\+/ /g')
  line11=$(sed -n '11p' "$file" | tr -d '\n' | sed 's/\s\+/ /g')
  
  echo -e "${sample}\t${line8}\t${line11}" >> $output
done

