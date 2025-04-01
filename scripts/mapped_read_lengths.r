library(tidyverse)
library(Rsamtools)

# prepare variables for mapping of disease group, batch and effect on to data
group_mapping <- c("S1" = "high_activity",
                   "S2" = "high_activity",
                   "S3" = "high_activity",
                   "S4" = "high_activity",
                   "S5" = "high_activity",
                   "S6" = "high_activity",
                   "S7" = "high_activity",
                   "S8" = "high_activity",
                   "S9" = "high_activity",
                   "S10" = "low_activity",
                   "S11" = "low_activity",
                   "S12" = "low_activity",
                   "S13" = "low_activity",
                   "S14" = "low_activity",
                   "S15" = "low_activity",
                   "S16" = "low_activity",
                   "S17" = "low_activity",
                   "S18" = "low_activity",
                   "S19" = "low_activity",
                   "S20" = "symptomatic_control",
                   "S21" = "symptomatic_control",
                   "S22" = "symptomatic_control",
                   "S23" = "symptomatic_control",
                   "S24" = "symptomatic_control",
                   "S25" = "symptomatic_control",
                   "S26" = "symptomatic_control",
                   "S27" = "symptomatic_control",
                   "S28" = "symptomatic_control",
                   "S29" = "symptomatic_control",
                   "S30" = "symptomatic_control",
                   "S31" = "low_activity",
                   "S32" = "low_activity",
                   "S33" = "symptomatic_control",
                   "S34" = "symptomatic_control",
                   "S35" = "symptomatic_control",
                   "S36" = "symptomatic_control")

sex_mapping <- c("S1" = "female",
                 "S2" = "female",
                 "S3" = "male",
                 "S4" = "female",
                 "S5" = "female",
                 "S6" = "female",
                 "S7" = "male",
                 "S8" = "male",
                 "S9" = "female",
                 "S10" = "male",
                 "S11" = "female",
                 "S12" = "female",
                 "S13" = "female",
                 "S14" = "male",
                 "S15" = "female",
                 "S16" = "female",
                 "S17" = "female",
                 "S18" = "female",
                 "S19" = "female",
                 "S20" = "female",
                 "S21" = "male",
                 "S22" = "female",
                 "S23" = "female",
                 "S24" = "female",
                 "S25" = "female",
                 "S26" = "female",
                 "S27" = "female",
                 "S28" = "female",
                 "S29" = "male",
                 "S30" = "female",
                 "S31" = "male",
                 "S32" = "male",
                 "S33" = "male",
                 "S34" = "female",
                 "S35" = "female",
                 "S36" = "female")

batch_mapping <- c("S1" = "RUN1",
                   "S2" = "RUN2",
                   "S3" = "RUN2",
                   "S4" = "RUN3",
                   "S5" = "RUN3",
                   "S6" = "RUN3",
                   "S7" = "RUN4",
                   "S8" = "RUN4",
                   "S9" = "RUN4",
                   "S10" = "RUN1",
                   "S11" = "RUN1",
                   "S12" = "RUN2",
                   "S13" = "RUN2",
                   "S14" = "RUN2",
                   "S15" = "RUN3",
                   "S16" = "RUN3",
                   "S17" = "RUN4",
                   "S18" = "RUN4",
                   "S19" = "RUN4",
                   "S20" = "RUN1",
                   "S21" = "RUN1",
                   "S22" = "RUN1",
                   "S23" = "RUN2",
                   "S24" = "RUN2",
                   "S25" = "RUN2",
                   "S26" = "RUN3",
                   "S27" = "RUN3",
                   "S28" = "RUN3",
                   "S29" = "RUN4",
                   "S30" = "RUN4",
                   "S31" = "RUN5",
                   "S32" = "RUN5",
                   "S33" = "RUN5",
                   "S34" = "RUN5",
                   "S35" = "RUN5",
                   "S36" = "RUN5")

custom_colors <- c("#9fca91", "#26788c", "#2c3573")

## code for loading in data
base_path <- "L:/LovbeskyttetMapper/Gridion/August/thesis_data/runs"

run_folders <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)

target_folders <- run_folders[grepl("RUN[0-9]+", basename(run_folders))]

# get full paths of all sample input files
input_files <- c()
for (i in 1:length(target_folders)) {
  subfolders <- list.dirs(target_folders[i], full.names = TRUE, recursive = TRUE)
  target_subfolder <- subfolders[grepl("length_90_to_250/mapped_reads", subfolders)]
  input_files <- c(input_files, list.files(target_subfolder, pattern = ".bam$", full.names = TRUE))
}

# create a named list with filename as key and full path as value
data_dict <- setNames(input_files, str_extract(basename(input_files), "^[^_.]+"))

# for each file bath load in the bam file with Rsamtools and save to a list
bam_list <- map(data_dict, function(file_path){
  bam <- scanBam(file_path)
  return(bam)
})

# extract read lengths from bam object and plot as histogram
tibble(read_length = bam_list$S1[[1]]$qwidth) %>%
  ggplot(aes(x = read_length)) +
  geom_histogram(fill = "aquamarine4", 
                 alpha = 0.5, 
                 binwidth = 1,
                 color = "black") +
  labs(title = "Distribution of read lengths", 
       x = "Read length", 
       y = "Read count") +
  theme_minimal()

# extract the relative amount of reads at each of the two specified read lengths
length_fractions <- map(bam_list, function(bam){
  x <- na.omit(bam[[1]]$qwidth)
  length_fraction <- sum(x == 155) / sum(x == 137)
  return(length_fraction)
})


read_lengths <- map(bam_list, function(bam){
  x <- na.omit(bam[[1]]$qwidth)
  length_tables <- table(unlist(x))
  return(length_tables)
})

# create tibble of the extracted relative distribution of reads
data <- tibble(sample_name = names(length_fractions), 
       frac = unlist(length_fractions)) %>%
  # add column for disease group, sex and batch
  mutate(group = group_mapping[sample_name],
         sex = sex_mapping[sample_name],
         batch = batch_mapping[sample_name]) %>%
  # add a column for disease status
  mutate(disease_status = ifelse(group == "symptomatic_control", "symptomatic_control", "multiple_sclerosis")) %>%
  # convert group, disease_status, sex and batch values to factors
  mutate(group = factor(group, levels = c("symptomatic_control", "low_activity", "high_activity")),
         disease_status = factor(disease_status, levels = c("symptomatic_control", "multiple_sclerosis")),
         sex = factor(sex, levels = c("female", "male")),
         batch = factor(batch, levels = unique(batch)))

# plot boxplots of the relative number of reads at each position for each MS groups
data %>%
  ggplot(aes(x = group, y = frac, color = group)) +
  geom_boxplot(aes(fill = group), alpha = 0.3, color = "black") +   # Boxplot
  geom_point(size = 2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Number of reads at long peak divided by number of reads at short peak", 
       y = "Number of reads at long peak for each read at short peak", 
       x = NULL) +
  theme(legend.position = "none") +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors)

# plot boxplots of the relative number of reads at each position for each MS disease status
data %>%
  ggplot(aes(x = disease_status, y = frac, color = disease_status)) +
  geom_boxplot(aes(fill = disease_status), alpha = 0.3, color = "black") +   # Boxplot
  geom_point(size = 2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Number of reads at long peak divided by number of reads at short peak", 
       y = "Number of reads at long peak for each read at short peak", 
       x = NULL) +
  theme(legend.position = "none") +
  scale_fill_manual(values = custom_colors[c(1,3)]) +
  scale_color_manual(values = custom_colors[c(1,3)])

# perform t-test comparing the mean of the relative number of reads at each 
# of two peaks between samples of different MS groups or MS disease status
data %>%
  wilcox.test(frac ~ disease_status, data = .)

data %>%
  filter(group != "low_activity") %>%
  wilcox.test(frac ~ group, data = .)

data %>%
  filter(group != "high_activity") %>%
  wilcox.test(frac ~ group, data = .)

data %>%
  kruskal.test(frac ~ group, data = .)