### subsample the biopsy genomes down to 0.5x, 1x, 2x, 3x ...

# This typically happens at the read level before aligning to the reference genome, ensuring that the downstream processes (alignment, variant calling, etc) reflect the desired coverage levels.
# Subsampling is done by randomly selecting a fraction of the reads from the FASTQ files to achieve the desired coverage.

# Required reads = Total reads x (Desired coverage/Current coverage)

# BCX0860
174511598 X (0.5 / 15.0644)
#fraction for 0.5x = 0.0331908

1 / 15.0644 = 0.0663817

2 / 15.0644 = 0.1327633

3 / 15.0644 = 0.1991450

### set seed (-s) for reproducibility and consistency across subsamples -- If you use the same seed, the reads selected for lower coverage levels (e.g., 1x) will be a subset of the reads selected for higher coverage levels (e.g., 2x, 3x). This ensures that the subsamples are nested and consistent.
### Facilitates Comparisons --- In analyses where you compare results across different coverage levels, having the same set of reads included in all subsamples (up to the desired fraction) eliminates randomness as a confounding factor.

#!/bin/bash

module purge 
module load seqtk/1.4

# subsampling to different coverage levels:
# 0.5x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.0331908 | gzip > subsampled_0.5x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.0331908 | gzip > subsampled_0.5x_BCX0860_R2.fastq.gz

# 1x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.0663817 | gzip > subsampled_1x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.0663817 | gzip > subsampled_1x_BCX0860_R2.fastq.gz

# 2x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.1327633 | gzip > subsampled_2x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.1327633 | gzip > subsampled_2x_BCX0860_R2.fastq.gz

# 3x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.199145 | gzip > subsampled_3x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.199145 | gzip > subsampled_3x_BCX0860_R2.fastq.gz

module purge 



### Now do everything to the other biopsy sample:

# Required reads = Total reads x (Desired coverage/Current coverage)

# BCX1596
184604714 X (0.5 / 15.8526)
#fraction for 0.5x = 0.0315406

1 / 15.8526 = 0.0630811

2 / 15.8526 = 0.1261623

3 / 15.8526 = 0.1892434

#!/bin/bash

module purge 
module load seqtk/1.4

# subsampling to different coverage levels:
# 0.5x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_10---IDT_i5_10.BCX1596_R1.fastq.gz 0.0315406 | gzip > subsampled_0.5x_BCX1596_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_10---IDT_i5_10.BCX1596_R2.fastq.gz 0.0315406 | gzip > subsampled_0.5x_BCX1596_R2.fastq.gz

# 1x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_10---IDT_i5_10.BCX1596_R1.fastq.gz 0.0630811 | gzip > subsampled_1x_BCX1596_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_10---IDT_i5_10.BCX1596_R2.fastq.gz 0.0630811 | gzip > subsampled_1x_BCX1596_R2.fastq.gz

# 2x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_10---IDT_i5_10.BCX1596_R1.fastq.gz 0.1261623 | gzip > subsampled_2x_BCX1596_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_10---IDT_i5_10.BCX1596_R2.fastq.gz 0.1261623 | gzip > subsampled_2x_BCX1596_R2.fastq.gz

# 3x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_10---IDT_i5_10.BCX1596_R1.fastq.gz 0.1892434 | gzip > subsampled_3x_BCX1596_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_10---IDT_i5_10.BCX1596_R2.fastq.gz 0.1892434 | gzip > subsampled_3x_BCX1596_R2.fastq.gz

module purge 


############################################################################################
############################################################################################

# March 27th 2025
# Serially downsample the biopsies to then estimate heterozygosity -- checking if 2x is reliable enough estimate of heterozygosity with the blow samples
# 15x, 14x, 13x, 12, 11, 10, 9, ... 5x, 4.5x, 4x, 3.5x, 3x, 2.5x, 2x, 1.5x, 1x, 0.5x, 0.25x



# BCX0860
174511598 X (0.5 / 15.0644)
#fraction for 0.5x = 0.0331908

0.25 / 15.0644 = 0.01659542
#0.5 / 15.0644 = 0.0331908
#1 / 15.0644 = 0.0663817
1.5 / 15.0644 = 0.0995725
#2 / 15.0644 = 0.1327633
2.5 / 15.0644 = 0.1659542
#3 / 15.0644 = 0.1991450
3.5 / 15.0644 = 0.2323358
4 / 15.0644 = 0.2655267
4.5 / 15.0644 = 0.2987175
5 / 15.0644 = 0.3319083
6 / 15.0644 = 0.39829
7 / 15.0644 = 0.4646717
8 / 15.0644 = 0.5310533
9 / 15.0644 = 0.597435
10 / 15.0644 = 0.6638167
11 / 15.0644 = 0.7301983
12 / 15.0644 = 0.79658
13 / 15.0644 = 0.8629617
14 / 15.0644 = 0.9293434


#!/bin/bash

module purge 
module load seqtk

# subsampling to different coverage levels:
# 0.25x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.01659542 | gzip > subsampled_0.25x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.01659542 | gzip > subsampled_0.25x_BCX0860_R2.fastq.gz

# 1.5x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.0995725 | gzip > subsampled_1.5x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.0995725 | gzip > subsampled_1.5x_BCX0860_R2.fastq.gz

# 2.5x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.1659542 | gzip > subsampled_2.5x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.1659542 | gzip > subsampled_2.5x_BCX0860_R2.fastq.gz

# 3.5x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.2323358 | gzip > subsampled_3.5x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.2323358 | gzip > subsampled_3.5x_BCX0860_R2.fastq.gz

# 4x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.2655267 | gzip > subsampled_4x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.2655267 | gzip > subsampled_4x_BCX0860_R2.fastq.gz

# 4.5x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.2987175 | gzip > subsampled_4.5x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.2987175 | gzip > subsampled_4.5x_BCX0860_R2.fastq.gz

# 5x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.3319083 | gzip > subsampled_5x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.3319083 | gzip > subsampled_5x_BCX0860_R2.fastq.gz

# 6x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.39829 | gzip > subsampled_6x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.39829 | gzip > subsampled_6x_BCX0860_R2.fastq.gz

# 7x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.4646717 | gzip > subsampled_7x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.4646717 | gzip > subsampled_7x_BCX0860_R2.fastq.gz

# 8x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.5310533 | gzip > subsampled_8x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.5310533 | gzip > subsampled_8x_BCX0860_R2.fastq.gz

# 9x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.597435 | gzip > subsampled_9x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.597435 | gzip > subsampled_9x_BCX0860_R2.fastq.gz

# 10x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.6638167 | gzip > subsampled_10x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.6638167 | gzip > subsampled_10x_BCX0860_R2.fastq.gz

# 11x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.7301983 | gzip > subsampled_11x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.7301983 | gzip > subsampled_11x_BCX0860_R2.fastq.gz

# 12x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.79658 | gzip > subsampled_12x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.79658 | gzip > subsampled_12x_BCX0860_R2.fastq.gz

# 13x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.8629617 | gzip > subsampled_13x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.8629617 | gzip > subsampled_13x_BCX0860_R2.fastq.gz

# 14x
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R1.fastq.gz 0.9293434 | gzip > subsampled_14x_BCX0860_R1.fastq.gz
seqtk sample -s 42 NS.LH00487_0027.007.IDT_i7_58---IDT_i5_58.BCX0860_R2.fastq.gz 0.9293434 | gzip > subsampled_14x_BCX0860_R2.fastq.gz


module purge 

### now for other biopsy


# BCX1596
184604714 X (0.5 / 15.8526)
#fraction for 0.5x = 0.0315406

0.25 / 15.8526 = 0.01577028
#0.5 / 15.8526 = 0.0315406
#1 / 15.8526 = 0.0630811
1.5 / 15.8526 = 0.0946217
#2 / 15.8526 = 0.1261623
2.5 / 15.8526 = 0.1577028
#3 / 15.8526 = 0.1892434
3.5 / 15.8526 = 0.2207840
4 / 15.8526 = 0.2523245
4.5 / 15.8526 = 0.2838651
5 / 15.8526 = 0.3154057
6 / 15.8526 = 0.3784868
7 / 15.8526 = 0.4415679
8 / 15.8526 = 0.5046491
9 / 15.8526 = 0.5677302
10 / 15.8526 = 0.6308113
11 / 15.8526 = 0.6938925
12 / 15.8526 = 0.7569736
13 / 15.8526 = 0.8200548
14 / 15.8526 = 0.8831359

# now just change the fractions using sed:
#0.25
sed -i 's/0.01659542/0.01577028/' subsample.reads.seqtk.additional.BCX1596.sh
#1.5
sed -i 's/0.0995725/0.0946217/' subsample.reads.seqtk.additional.BCX1596.sh
#2.5
sed -i 's/0.1659542/0.1577028/' subsample.reads.seqtk.additional.BCX1596.sh
#3.5
sed -i 's/0.2323358/0.2207840/' subsample.reads.seqtk.additional.BCX1596.sh
#4
sed -i 's/0.2655267/0.2523245/' subsample.reads.seqtk.additional.BCX1596.sh
#4.5
sed -i 's/0.2987175/0.2838651/' subsample.reads.seqtk.additional.BCX1596.sh
#5
sed -i 's/0.3319083/0.3154057/' subsample.reads.seqtk.additional.BCX1596.sh
#6
sed -i 's/0.39829/0.3784868/' subsample.reads.seqtk.additional.BCX1596.sh
#7
sed -i 's/0.4646717/0.4415679/' subsample.reads.seqtk.additional.BCX1596.sh
#8
sed -i 's/0.5310533/0.5046491/' subsample.reads.seqtk.additional.BCX1596.sh
#9
sed -i 's/0.597435/0.5677302/' subsample.reads.seqtk.additional.BCX1596.sh
#10
sed -i 's/0.6638167/0.6308113/' subsample.reads.seqtk.additional.BCX1596.sh
#11
sed -i 's/0.7301983/0.6938925/' subsample.reads.seqtk.additional.BCX1596.sh
#12
sed -i 's/0.79658/0.7569736/' subsample.reads.seqtk.additional.BCX1596.sh
#13
sed -i 's/0.8629617/0.8200548/' subsample.reads.seqtk.additional.BCX1596.sh
#14
sed -i 's/0.9293434/0.8831359/' subsample.reads.seqtk.additional.BCX1596.sh