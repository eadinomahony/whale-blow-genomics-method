#### genome-wide heterozygosity 

### min depth 3x, submit script per sample:

#!/bin/bash

module purge
module load angsd/0.940-2

bams=bam_MN-b-006.txt
ref=GCA_041834305.1_ASM4183430v1_genomic.fna
chr=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt

for sample in $(cat bam_MN-b-006.txt)
do

id=`echo ${sample} | cut -f1 -d "." | cut -f9 -d "/"`

angsd -b ${bams} -ref ${ref} -anc ${ref} -rf ${chr} -out ${id}_mindepth3x -nThreads 8 -GL 2 -minMapQ 30 -minQ 30 -doSaf 1 -doMajorMinor 1 -setMinDepth 3 -doCounts 1 -skipTriallelic 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1

realSFS -fold 1 ${id}_mindepth3x.saf.idx > ${id}_mindepth3x_est.ml

done

module purge



### min depth 7x, submit script per sample:

#!/bin/bash

module purge
module load angsd/0.940-2

bams=bam_MN-b-006.txt
ref=GCA_041834305.1_ASM4183430v1_genomic.fna
chr=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt

for sample in $(cat bam_MN-b-006.txt)
do

id=`echo ${sample} | cut -f1 -d "." | cut -f9 -d "/"`

angsd -b ${bams} -ref ${ref} -anc ${ref} -rf ${chr} -out ${id}_mindepth7x -nThreads 8 -GL 2 -minMapQ 30 -minQ 30 -doSaf 1 -doMajorMinor 1 -setMinDepth 7 -doCounts 1 -skipTriallelic 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1

realSFS -fold 1 ${id}_mindepth7x.saf.idx > ${id}_mindepth7x_est.ml

done

module purge


### min depth 10x, submit script per sample:

#!/bin/bash

module purge
module load angsd/0.940-2

bams=bam_MN-b-006.txt
ref=GCA_041834305.1_ASM4183430v1_genomic.fna
chr=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt

for sample in $(cat bam_MN-b-006.txt)
do

id=`echo ${sample} | cut -f1 -d "." | cut -f9 -d "/"`

angsd -b ${bams} -ref ${ref} -anc ${ref} -rf ${chr} -out ${id}_mindepth10x -nThreads 8 -GL 2 -minMapQ 30 -minQ 30 -doSaf 1 -doMajorMinor 1 -setMinDepth 10 -doCounts 1 -skipTriallelic 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1

realSFS -fold 1 ${id}_mindepth10x.saf.idx > ${id}_mindepth10x_est.ml

done

module purge


#### use min depth 3x for Figure 3; and use this script for calculation of heterozygosity on merged data, biopsy data and publicly available data