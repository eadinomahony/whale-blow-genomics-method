### unmerged and public data

#!/bin/bash

module purge
module load angsd/0.940
module load pcangsd/1.36.1  

ref=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt
bams=postQC_unmerged_public_bam_paths_norelatives_perfindiv.txt
outpath=/pcangsd/unmerged_biopsies_public_mindepth3x

angsd -GL 2 \
    -doGlf 2 \
    -doMajorMinor 1 \
    -minMapQ 30 -minQ 30 -SNP_pval 1e-6 \
    -doMaf 1 \
    -minMaf 0.05 \
    -rf ${ref}  \
    -nThreads 10 \
    -bam ${bams} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -skipTriallelic 1  -setMinDepthInd 3 -out ${outpath}/pcangsd.prep.unmerged.public.norelatives.perfindiv

### count the # of sites:
zcat ${outpath}/pcangsd.prep.unmerged.public.norelatives.perfindiv.mafs.gz | tail -n +2 | wc -l > sites_count_unmerged_public_perfindiv.txt

###

beagle=pcangsd.prep.unmerged.public.norelatives.perfindiv.beagle.gz # re-enter the beagle.gz file name for p1

pcangsd --beagle ${beagle} -o unmerged.blow.public.pca.output --threads 20

module purge



### merged and public data

#!/bin/bash

module purge
module load angsd/0.940
module load pcangsd/1.36.1 

ref=GCA_041834305.1_ASM4183430v1_genomic_autosomes.txt
bams=postQC_merged_public_bam_paths_norelatives_perfindiv.txt
outpath=/mergedData_analyses/pcangsd


angsd -GL 2 \
    -doGlf 2 \
    -doMajorMinor 1 \
    -minMapQ 30 -minQ 30 -SNP_pval 1e-6 \
    -doMaf 1 \
    -minMaf 0.05 \
    -rf ${ref}  \
    -nThreads 10 \
    -bam ${bams} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -skipTriallelic 1  -setMinDepthInd 3 -out ${outpath}/pcangsd.prep.merged.public.norelatives.perfindiv

### count the # of sites:
zcat ${outpath}/pcangsd.prep.merged.public.norelatives.perfindiv.mafs.gz | tail -n +2 | wc -l > sites_count_merged_public_perfindiv.txt

###

beagle=pcangsd.prep.merged.public.norelatives.perfindiv.beagle.gz # re-enter the beagle.gz file name for p1

pcangsd --beagle ${beagle} -o merged.public.pca.output --threads 20

module purge