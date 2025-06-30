################################################################################
################################################################################

##### O'Mahony et al. (submitted June 2025) Nature Methods
#### ROHan plot for unmerged blow samples:

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(dplyr)
library(ggeasy)
library(tidyverse)

################################################################################
################################################################################
#### unmerged ROHan run

roh.summary <- read.csv("../data/roh.summary.txt", sep="")

# remove flagged individuals from perfect sample approach: 
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-103_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-127_unmergedblow_rohan.mid.hmmrohl")

# only retain biopsies subsampled to 2x:
roh.summary <- roh.summary %>%
  filter(ID != "BCX0860_subsampled_0.5x")
roh.summary <- roh.summary %>%
  filter(ID != "BCX0860_subsampled_1x")
roh.summary <- roh.summary %>%
  filter(ID != "BCX0860_subsampled_3x")
roh.summary <- roh.summary %>%
  filter(ID != "BCX1596_subsampled_0.5x")
roh.summary <- roh.summary %>%
  filter(ID != "BCX1596_subsampled_1x")
roh.summary <- roh.summary %>%
  filter(ID != "BCX1596_subsampled_3x")

### remove samples filtered out by QC
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-001_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-002_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-003_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-004_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-009_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-011_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-012_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-013_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-014_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-016_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-017_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-018_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-019_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-020_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-021_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-022_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-024_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-025_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-026_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-028_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-029_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-030_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-031_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-032_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-033_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-034_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-035_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-037_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-038_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-039_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-041_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-042_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-043_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-044_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-045_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-046_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-047_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-048_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-049_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-050_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-051_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-052_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-053_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-054_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-055_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-056_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-057_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-058_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-059_unmergedblow_rohan.mid.hmmrohl")

roh.summary <- roh.summary %>%
  filter(ID != "MN-b-060_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-061_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-062_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-065_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-066_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-067_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-068_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-069_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-070_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-071_unmergedblow_rohan.mid.hmmrohl")


roh.summary <- roh.summary %>%
  filter(ID != "MN-b-072_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-074_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-075_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-076_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-077_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-078_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-079_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-082_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-084_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-086_unmergedblow_rohan.mid.hmmrohl")

roh.summary <- roh.summary %>%
  filter(ID != "MN-b-089_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-090_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-091_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-092_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-094_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-095_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-097_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-098_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-099_unmergedblow_rohan.mid.hmmrohl")

roh.summary <- roh.summary %>%
  filter(ID != "MN-b-100_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-101_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-102_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-103_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-105_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-106_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-107_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-108_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-109_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-110_unmergedblow_rohan.mid.hmmrohl")

roh.summary <- roh.summary %>%
  filter(ID != "MN-b-111_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-116_unmergedblow_rohan.mid.hmmrohl")
roh.summary <- roh.summary %>%
  filter(ID != "MN-b-126_unmergedblow_rohan.mid.hmmrohl")

### pull in merged dataset to get public genome ROH
merged.roh <- read.csv("../data/roh.summary.merged.txt", sep="")

# filter out merged blow samples
biopsy.roh <- merged.roh %>%
  filter(!grepl("_mergedblow$", ID))

roh.summaryB <- rbind(biopsy.roh, roh.summary)

roh.summary2 <- roh.summaryB %>%
  separate(ID, into = c("ID2"), sep = "_", remove = FALSE) %>%
  mutate(VALIDATED_SITES=VALIDATED_SITES/1000000) %>%
  #mutate(bin = cut(MB, breaks = seq(0, 10, by = 3), include.lowest = TRUE, right = FALSE))
  mutate(bin = cut(VALIDATED_SITES, breaks = c(0.007, 1, 5, 10000), labels = c("<1 Mb", "1-5 Mb", ">5 Mb"), include.lowest = TRUE, right = FALSE))

head(roh.summary2)

roh.summary2 <- roh.summary2 %>%
  mutate(Type = case_when(
    grepl("^MN", ID) ~ "Blow",
    grepl("^BC", ID) ~ "Pacific",
    grepl("^SRR", ID) ~ "Atlantic",
    TRUE ~ ID  # Keep the species names unchanged
  ))

roh.summary2 <- roh.summary2 %>% 
  mutate(Type = ifelse(ID == "SRR17854487_public", "Pacific", Type))

roh.summary2 <- roh.summary2 %>% 
  mutate(ID2 = ifelse(ID == "BCX0860_dfo", "BCX0860_15x", ID2))
roh.summary2 <- roh.summary2 %>% 
  mutate(ID2 = ifelse(ID == "BCX1596_dfo", "BCX1596_16x", ID2))
roh.summary2 <- roh.summary2 %>% 
  mutate(ID2 = ifelse(ID == "BCX0860_subsampled_2x", "BCX0860_2x", ID2))
roh.summary2 <- roh.summary2 %>% 
  mutate(ID2 = ifelse(ID == "BCX1596_subsampled_2x", "BCX1596_2x", ID2))

head(roh.summary2)

roh.summary2 <- roh.summary2 %>%
  mutate(ID = factor(ID, levels = roh.summary2 %>%
                       arrange(Type) %>%  # Order by Type: Blow → Biopsy
                       pull(ID) %>% 
                       unique()))

### order within each 'Type' in decreasing counts:

# count number of ROH per ID
roh_order <- roh.summary2 %>%
  group_by(Type, ID) %>%
  summarise(total_roh = n(), .groups = "drop") %>%
  arrange(Type, desc(total_roh)) %>%
  pull(ID)

# apply new order to the ID factor
roh.summary2 <- roh.summary2 %>%
  mutate(ID = factor(ID, levels = roh_order))


### replot 
rohan.plot2 <- roh.summary2 %>%
  ggplot(aes(x = ID, fill = bin)) + # ,color=Type
  geom_bar(position = position_stack(reverse = TRUE)) +
  labs(
    x = "ID",
    y = "Count (ROH)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #   axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, hjust = 1, color = "black"),
    # axis.text.y = element_text(size = 14, color = "black")
  ) +
  easy_remove_x_axis(what = "title") +
  scale_fill_manual(
    values = c("plum1", "#B1E7F9", "#FA1934"),
    name = "ROH length bins"
  )+
  labs(fill = "ROH length bins") +
  scale_x_discrete(
    labels = setNames(
      roh.summary2 %>% pull(ID2),
      roh.summary2$ID
    )
  )+ 
  guides(color = "none") +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = "black")
  )

rohan.plot2

#ggsave(filename = "../figures/rohanplot.pdf", plot = rohan.plot2, width = 10, height = 5.5, dpi = 300, device = "pdf")
