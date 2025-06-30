##### perfect sample approach blow samples methods papers

# Set working directory to folder that this R file is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(purrr)
library(ggbreak)

errorest <- readr::read_table("../data/errorEst.txt")
head(errorest)
tail(errorest)

errorest <- errorest %>%
  slice(-1:-102) %>% 
  select(1:2) %>%      
  rename(sample = 1, errorRate = 2)
head(errorest)

# clean some more
errorest <- errorest %>%
  mutate(sample = str_replace_all(sample, '[\\"/]', ''))
head(errorest)


# colour according to depth of coverage:
depth_merged <- rename(depth_merged, sample = X1)

errorest_with_depth <- errorest %>%
  left_join(depth_merged, by = "sample")

errorest_with_depth <- errorest_with_depth %>%
  mutate(nuDNA_depth = if_else(sample == "BCX0860_subsampled_0.5x", 0.5, nuDNA_depth)) %>%
  mutate(nuDNA_depth = if_else(sample == "BCX1596_subsampled_0.5x", 0.5, nuDNA_depth)) %>%
  mutate(nuDNA_depth = if_else(sample == "BCX0860_subsampled_1x", 1, nuDNA_depth)) %>%
  mutate(nuDNA_depth = if_else(sample == "BCX1596_subsampled_1x", 1, nuDNA_depth)) %>%
  mutate(nuDNA_depth = if_else(sample == "BCX0860_subsampled_2x", 2, nuDNA_depth)) %>%
  mutate(nuDNA_depth = if_else(sample == "BCX1596_subsampled_2x", 2, nuDNA_depth)) %>%
  mutate(nuDNA_depth = if_else(sample == "BCX0860_subsampled_3x", 3, nuDNA_depth)) %>%
  mutate(nuDNA_depth = if_else(sample == "BCX1596_subsampled_3x", 3, nuDNA_depth))


errorest_with_depth <- errorest_with_depth %>%
  mutate(color_group = case_when(
    nuDNA_depth < 1 ~ "red",
    nuDNA_depth < 1.5 ~ "orange",
    nuDNA_depth < 2 ~ "yellow",
    TRUE           ~ "black"
  ))

# Plot with colors mapped

ggplot(errorest_with_depth, aes(x = sample, y = errorRate, color = color_group)) +
  geom_point(size = 3) +
  scale_y_break(c(0.7, 4.4), scales = 0.5) +
  coord_cartesian(ylim = c(0, 5.5)) +  
  scale_color_identity() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red")


### filter the plot to remove any biopsy subsamples:

errorest_with_depth_filt <- errorest_with_depth %>%
  slice(-2:-5) %>%
  slice(-3:-6)

error.plot <- ggplot(errorest_with_depth_filt, aes(x = sample, y = errorRate, color = color_group)) +
  geom_point(size = 3) +
  scale_y_break(c(0.7, 4.4), scales = 0.5) +
  coord_cartesian(ylim = c(0, 5.5)) +  
  scale_color_identity() +
  theme_bw()+
  theme(    
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red")+
  xlab("Sample") +
  ylab("Error rates (%)"); error.plot

#ggsave(filename = "../figures/errorRatesplot.png", error.plot, width = 18, height = 9, units = "in", dpi = 300)

### filter to just the unmerged samples:

unmerged_samples <- errorest_with_depth_filt %>%
  filter(grepl("^MN-b-", sample))

ggplot(unmerged_samples, aes(x = sample, y = errorRate, color = color_group)) +
  geom_point(size = 3) +
  scale_y_break(c(0.7, 4.4), scales = 0.5) +
  coord_cartesian(ylim = c(0, 5.5)) +  
  scale_color_identity() +
  theme_bw()+
  theme(    
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red")+
  xlab("Sample") +
  ylab("Error rates (%)")

### filter to just the refined list (<1x and lower qual replicates removed)

final_list <- c("MN-b-006", "MN-b-007","MN-b-008", "MN-b-010", "MN-b-015", 
                "MN-b-023", "MN-b-027", "MN-b-036", "MN-b-040", "MN-b-063", "MN-b-064",
                "MN-b-073", "MN-b-080", "MN-b-081", "MN-b-083", "MN-b-085", 
                "MN-b-087", "MN-b-088", "MN-b-093","MN-b-096","MN-b-103","MN-b-104",
                "MN-b-115","MN-b-127")

refined_unmerged <- unmerged_samples %>%
  filter(sample %in% final_list)


ggplot(refined_unmerged, aes(x = sample, y = errorRate, color = color_group)) +
  geom_point(size = 3) +
  #  coord_cartesian(ylim = c(0, 5.5)) +  
  scale_color_identity() +
  theme_bw()+
  theme(    
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red")+
  xlab("Sample") +
  ylab("Error rates (%)")

# Calculate the 95th percentile of errorRate
threshold <- quantile(refined_unmerged$errorRate, 0.95, na.rm = TRUE)

# replot with this new threshold
error.unmerged.plot <- ggplot(refined_unmerged, aes(x = sample, y = errorRate)) +
  geom_point(size = 3) +
  scale_color_identity() +
  theme_bw()+
  theme(    
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size=1)+
  xlab("Unmerged blow samples") +
  ylab("Error rates (%)")

# Filter the dataframe to keep only samples BELOW that threshold
refined_unmerged_threshold <- refined_unmerged %>%
  filter(errorRate <= threshold)
nrow(refined_unmerged_threshold)
#22 samples

################################################################################
################################################################################

##### repeat with merged data:
merged_samples <- errorest_with_depth_filt %>%
  filter(grepl("^BC|^CS", sample))

#remove related individuals (lower depth of pair with KING >0.125)
merged_samples <- merged_samples %>%
  filter(sample != "BCX0860") %>%
  filter(sample != "BCX1596") %>%
  filter(sample != "BCX0171") %>%
  filter(sample != "BCX2052") %>%
  filter(sample != "BCY0965")


# from ngsrelate section:
filtered_data <- ngsrelate.merged.clean %>%
  +   filter(KING >= 0.125)
filtered_data
# A tibble: 5 × 6
#ida     idb             nSites     R0     R1  KING
#<chr>   <chr>            <dbl>  <dbl>  <dbl> <dbl>
#1 BCX0171 BCY0117        4389718 0.0345  0.415 0.213
#2 BCX0427 BCX2052        4378315 0.0299  0.414 0.214
#3 BCX0860 BCX0860_biopsy 4436291 0       5.92  0.461
#4 BCX1484 BCY0965        4478611 0.185   0.349 0.135
#5 BCX1596 BCX1596_biopsy 4404481 0      13.6   0.482

# Calculate the 95th percentile of errorRate
threshold_merged <- quantile(merged_samples$errorRate, 0.95, na.rm = TRUE)

#relabel the biopsies and colour them to show they are the biopsies
library(stringr)

merged_samples <- merged_samples %>%
  mutate(
    biopsy_flag = if_else(str_detect(sample, "_biopsy"), 1, 0),
    sample = str_remove(sample, "_biopsy")  # remove the suffix
  )

merged_samples <- merged_samples %>%
  mutate(
    color_group = if_else(biopsy_flag == 1, "blue", "black")
  )

# replot with this new threshold
error.merged.plot <- ggplot(
  merged_samples,
  aes(x = sample, y = errorRate, color = color_group, shape = factor(biopsy_flag))
) +
  geom_point(size = 3) +
  scale_color_identity() +  # use the actual color values from color_group
  scale_shape_manual(values = c(`0` = 16, `1` = 15)) +  # 16 = circle, 8 = star
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  geom_hline(
    yintercept = threshold_merged,
    linetype = "dashed",
    color = "red",
    size = 1
  ) +
  xlab("Merged samples") +
  ylab("Error rates (%)")

error.merged.plot

# Filter the dataframe to keep only samples BELOW that threshold
merged_threshold <- merged_samples %>%
  filter(errorRate <= threshold_merged)
nrow(merged_threshold)
# 20 individuals


### put the plot together:

error.unmerged.plot
error.merged.plot

library(cowplot)


error_plot <- plot_grid(
  error.unmerged.plot,
  error.merged.plot,
  ncol = 1,
  align = "v",           
  rel_heights = c(1, 1),
  labels = c("A", "B"),
  label_size = 16          
)
error_plot

#ggsave("../figures/errorplot_panel.pdf", error_plot, width = 18, height = 9, units = "in", dpi = 300)
