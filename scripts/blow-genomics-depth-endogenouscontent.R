################################################################################
################################################################################

#### Éadin O'Mahony et al. (submitted June 2025) Nature Methods
#### Endogenous content and depth of blow samples

################################################################################
################################################################################

# Set working directory to folder that this R file is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## dependencies 
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("viridis")
library(dplyr)
library(tidyr)
library(viridis)
library(scales)
library(ggplot2)
library(ggpubr)

################################################################################
################################################################################

# mean raw reads
reads <- read.table("../data/raw.reads.blow.methodspaper.txt", header = F)
mean(reads$V4) # total reads
mean(reads$V2) # forward reads

endogenous <- read.csv("/Users/eadinomahony/Desktop/Desktop/WORK/PhD/BlowProject/Humpback-blow-genomics/data/endogenous_content_flagstat.csv", header = TRUE)
head(endogenous)

#bring in list of methods paper samples to filter down
sample_list <- read.table("../data/methods_paper_samples_sorted.txt", header = F)

colnames(sample_list) <- c("sample"); head(sample_list)

#remove the KingFisher samples:
endog_methodspaper <- endogenous[-1,]
endog_methodspaper <- endog_methodspaper[-2,]
endog_methodspaper <- endog_methodspaper[-3,]
head(endog_methodspaper)

# filter down the samples to methods paper samples:
filtered_endogenous <- endog_methodspaper %>%
  filter(sample %in% sample_list$sample)

head(filtered_endogenous)
nrow(filtered_endogenous) # 57 :)

mean(filtered_endogenous$percentage)
# mean endogenous content = 83.59497 %

sd(filtered_endogenous$percentage)
# std dev = 12.79 %

min(filtered_endogenous$percentage)
max(filtered_endogenous$percentage)

### unmapped reads
filtered_endogenous <- filtered_endogenous %>% 
  mutate(unmapped = total_reads - reads_mapped)

filtered_endogenous <- filtered_endogenous %>% 
  mutate(unmapped_percentage = 100 - percentage)

head(filtered_endogenous)

mean(filtered_endogenous$unmapped)
# 18,090,163

################################################################################
################################################################################

# depth unmerged samples

depth <- read.csv("../data/depth-methodspaper.csv", header = T)
head(depth)
nrow(depth)

# need to remove MN-b-075 because it was re-IDed to a whale without a duplicate:
depth <- depth %>%
  filter(sample != "MN-b-075")

nrow(depth)
# 58 

mean(depth$depth)
# 2.336773x

mean(depth$var)
# 2.24031x

min(depth$depth)
# 0.19x

max(depth$depth)
# 3.45x

# Think about removing the samples < 1x depth (MN-b-117, MN-b-013, MN-b-094)

depth_filtered <- depth %>%
  filter(depth > 1.0)

nrow(depth_filtered)

mean(depth_filtered$depth)
# 2.45x

mean(depth_filtered$var)
# 2.32x

min(depth_filtered$depth)



################################################################################
################################################################################

# correlation with lower depth and lower endogenous content?

# test for normality first:
shapiro.test(depth$depth)
shapiro.test(filtered_endogenous$endogenous_content)
# Both p ≤ 0.05, the data is not normally distributed → use Spearman's correlation.

#make sure they line up with each other:
total_df <- inner_join(depth, filtered_endogenous, by = "sample")
cor.test(total_df$depth, total_df$endogenous_content, method = "spearman")

# visualise
endog_depth_corr_plot <- ggplot(total_df, aes(x = depth, y = percentage)) +
  geom_point() +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm", color = "#009988") + # Linear regression line
  labs(x = "Depth of Coverage",
       y = "Endogenous Content (%)")+
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    size = 4,
    label.x = 2,
    label.y = 52,
    color = "black"
  ) 

#ggsave("../figures/endog_depth_corr_plot.png", plot = endog_depth_corr_plot, width = 12, height = 8, units = "cm", dpi = 300)

################################################################################
################################################################################

# depth merged samples, unique whales

depth_merged <- read.csv("../data/merged_depth.csv", header = T)
head(depth_merged)

new_merged_depth <- readr::read_table("../data/mapped_merged_duplicates_coverage.txt", col_names = FALSE)

depth_merged <- new_merged_depth

# remove BCX0860, because it's merged with a 3x biopsy:
#depth_merged <- depth_merged %>%
#  filter(X1 != "BCX0860")

# remove BCX0147, because of false field ID:
depth_merged <- depth_merged %>%
  filter(X1 != "BCX0147")

mean(depth_merged$X2)

mean(depth_merged$X3) # standard deviation

min(depth_merged$X2)

max(depth_merged$X2)

nrow(depth_merged)

#sort by increasing depth for Table S2:
sorted_merged <- depth_merged[order(depth_merged$merged_ID),]
#write.csv(x = sorted_merged, file = "../data/sorted_merged_depth_methodspaper.csv")


################################################################################
################################################################################

#### corr plot between nuDNA and mtDNA depths
"North Pacific (unmerged blow)" = "#7CAE00", "North Pacific (biopsy)" = "#FFB000", "North Pacific (merged blow)" = "#648FFF")
cols <- c("unmerged" = "#009988", "merged" = "#648FFF")

depth_merged <- read.csv("../data/nuDNA_vs_mtDNA_depth.csv", header = T)
head(depth_merged)

depth_corr <- ggplot(depth_merged, aes(x = nuDNA_depth, y = mtDNA_depth, color = Type)) +
  geom_point(aes(shape = Type), size = 3, show.legend = TRUE) +
  geom_smooth(aes(linetype = Type), method = "lm", se = TRUE, size = 0.8) +
  
  # stat_cor for "merged"
  stat_cor(
    data = filter(depth_merged, Type == "merged"),
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    size = 5.5,
    label.x = 4.9,    
    label.y = 150,  
    color = "black" 
  ) +
  
  # stat_cor for "unmerged"
  stat_cor(
    data = filter(depth_merged, Type == "unmerged"),
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    size = 5.5,
    label.x = 2.2,
    label.y = 28,
    color = "black"
  ) +
  
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c("merged" = 8, "unmerged" = 6)) +
  scale_linetype_manual(values = c("merged" = "solid", "unmerged" = "solid")) +
  theme_bw(base_size = 12) +
  xlab("Nuclear genome depth") +
  ylab("Mitogenome depth") +
  theme(
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 18),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )
depth_corr

#ggsave(filename = "../figures/corr_nuDNA_mtDNA_depth.pdf", plot = depth_corr, width = 8, height = 5.5, dpi = 300)

################################################################################
################################################################################

### boxplot with endogenous content for nuclear and mtDNA

endog_both <- read.csv("../data/endogenous_nuDNA_mtDNA.csv", header = T)
head(endog_both)

#bring in list of methods paper samples to filter down
sample_list <- read.table("../data/methods_paper_samples_sorted.txt", header = F)

colnames(sample_list) <- c("sample"); head(sample_list)

#remove the KingFisher samples:
endog_both <- endog_both[-1,]
endog_both <- endog_both[-2,]
endog_both <- endog_both[-3,]
head(endog_both)

# filter down the samples to methods paper samples:
filtered_both <- endog_both %>%
  filter(sample %in% sample_list$sample)

head(filtered_both)
nrow(filtered_both) # 57 :)

mean(filtered_both$percentage)
# mean endogenous content = 83.59497 %

mean(filtered_both$mtDNA_percent)
# 0.0378971

sd(filtered_both$mtDNA_percent)
# 0.01216623

# shift into long format
library(tidyverse)

# Reshape your data from wide to long format
endog_long <- filtered_both %>%
  select(sample, endogenous_content, mtDNA_endogenous) %>%
  pivot_longer(cols = c(endogenous_content, mtDNA_endogenous),
               names_to = "Metric",
               values_to = "Value")
head(endog_long)


#try a log scale:


"#F5793A", "#A95AA1", "#85C0F9", "#0F2080"
"#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000"
"#0077BB", "#EE7733", "#33BBEE", "#CC3311", "#009988", "#EE3377", "#BBBBBB"

ggplot(endog_long, aes(x = Metric, y = Value)) +
  geom_boxplot() +
  scale_y_log10(labels = scales::percent_format(accuracy = 0.01)) +
  scale_x_discrete(labels = c(
    "endogenous_content" = "Total DNA",
    "mtDNA_endogenous" = "mtDNA"
  )) +
  theme_bw(base_size = 12) +
  labs(x = "", y = "Endogenous content (log scale)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18, margin = margin(t = 16, r = 16, b = 16, l = 16)),
    axis.text = element_text(size = 16)
  )



### just show nuDNA endogenous:

endog_plot <- ggplot(filtered_both, aes(x = extraction, y = percentage)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(
    "percentage" = "nuDNA"
  )) +
  scale_y_continuous(limits = c(0, 100))+
  theme_bw(base_size = 12) +
  labs(x = "nuDNA", y = "Endogenous content (%)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18), #, margin = margin(t = 16, r = 16, b = 16, l = 16)
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16)
  )

endog_plot

#ggsave("../figures/endogenous_content_boxplot2.pdf", plot = endog_plot, width = 7, height = 5.5, dpi = 300)
