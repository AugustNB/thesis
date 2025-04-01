library(tidyverse)
library(brglm2)
library(ggrepel)

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

# set custom colors to be used for three disease groups
custom_colors <- c("#9fca91", "#26788c", "#2c3573")

# get filepaths for methylation data 
base_path <- "L:/LovbeskyttetMapper/Gridion/August/thesis_data/methylation_data"
input_files <- list.files(base_path, pattern = "*.bed.gz$", full.names = TRUE)

# create a named list with filename as key and full path as value
data_dict <- setNames(input_files, basename(input_files))

# read in data and add column with sample name
df_list <- lapply(names(data_dict), function(key) {
  df <- read_delim(data_dict[[key]], delim = "\t", col_names = FALSE,
                   show_col_types = FALSE)
  df$sample_name <- str_extract(key, "^[^_.]+")
  return(df)
})

# merge dataframes
all_data <- bind_rows(df_list)

subset_data <- all_data %>% 
  # only use positions with depth of 10 or greater
  filter(X10 >= 10) %>%
  # concatenate position information
  mutate(position = str_c(X1, X2, X3, sep = "_"),
         # add column for disease group, sex and batch
         group = group_mapping[sample_name],
         sex = sex_mapping[sample_name],
         batch = batch_mapping[sample_name],
         # select and rename columns with methylation data
         n_mod = X12,
         n_unmod = X10-X12) %>%
  # add a column for disease status
  mutate(disease_status = ifelse(group == "symptomatic_control", "symptomatic_control", "multiple_sclerosis")) %>%
  # convert group, disease_status, sex and batch values to factors for correct levels for tests
  mutate(group = factor(group, levels = c("symptomatic_control", "low_activity", "high_activity")),
         disease_status = factor(disease_status, levels = c("symptomatic_control", "multiple_sclerosis")),
         sex = factor(sex, levels = c("female", "male")),
         batch = factor(batch, levels = c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")))

# create subset with only 5hmc methylation data
h_subset <- subset_data %>%
  filter(X4 == "h") %>%
  select(c(position, sample_name, group, sex, batch, disease_status, n_mod, n_unmod)) 

# create subset with only 5mc methylation data
m_subset <- subset_data %>%
  filter(X4 == "m" & n_mod >= 1) %>%
  select(c(position, sample_name, group, sex, batch, disease_status, n_mod, n_unmod))

dim(h_subset)
dim(m_subset)

## H_SUBSET PLOTS

# plot methylation fraction for disease_status
plotly::ggplotly(
  ggplot(h_subset, aes(x = disease_status, y = n_mod/(n_mod+n_unmod), 
                       color = disease_status, fill = disease_status)) +
    geom_violin(alpha = 0.2, color = "black") + 
    geom_jitter(width = 0.1, size = 2, alpha = 0.3) + 
    theme_minimal() + 
    theme(legend.position = "none") +
    labs(title = "5hmC methylated fraction at genomic positions with detected methylation",
         x = NULL,
         y = "Fraction methylated",
         color = "MS disease status") +
    scale_color_manual(values = custom_colors[c(1,3)]) +
    scale_fill_manual(values = custom_colors[c(1,3)])
)

# plot methylation fraction for groups
plotly::ggplotly(
  ggplot(h_subset, aes(x = group, y = n_mod/(n_mod+n_unmod), 
                       color = group, fill = group)) +
    geom_violin(alpha = 0.2, color = "black") +
    geom_jitter(width = 0.1, size = 2, alpha = 0.3) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "5hmC methylated fraction at genomic positions with detected methylation",
         x = NULL,
         y = "Fraction methylated",
         color = "MS disease status") +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors)
)

# perform binomial logistic regression correlating 5hmC methylation at each site to
# MS status
h_subset %>%
  glm(cbind(n_mod,n_unmod) ~ disease_status + batch + sex, family = "binomial", data = .) %>%
  summary(.)

h_subset %>%
  glm(cbind(n_mod,n_unmod) ~ group + batch + sex, family = "binomial", data = .) %>%
  summary(.)

h_subset %>%
  filter(group != "high_activity") %>%
  glm(cbind(n_mod,n_unmod) ~ group + batch + sex, family = "binomial", data = .) %>%
  summary(.)

## M_SUBSET PLOTS

# plot methylation fraction for disease_status
plotly::ggplotly(
  ggplot(m_subset, aes(x = disease_status, y = n_mod/(n_mod+n_unmod), 
                       color = disease_status, fill = disease_status)) +
    geom_violin(alpha = 0.2, color = "black") +  # Violin plot with transparency
    geom_jitter(width = 0.1, size = 2, alpha = 0.1) +  # Jittered points
    theme_minimal() +  # Clean theme
    theme(legend.position = "none") +
    labs(title = "5mC methylated fraction at genomic positions with detected methylation",
         x = NULL,
         y = "Fraction methylated",
         color = "MS disease status") +
    scale_color_manual(values = custom_colors[c(1,3)]) +
    scale_fill_manual(values = custom_colors[c(1,3)])
)

# plot methylation fraction for groups
plotly::ggplotly(
  ggplot(m_subset, aes(x = group, y = n_mod/(n_mod+n_unmod), 
                       color = group, fill = group)) +
    geom_violin(alpha = 0.2, color = "black") +  # Violin plot with transparency
    geom_jitter(width = 0.1, size = 2, alpha = 0.1) +  # Jittered points
    theme_minimal() +  # Clean theme
    theme(legend.position = "none") +
    labs(title = "5mC methylated fraction at genomic positions with detected methylation",
         x = NULL,
         y = "Fraction methylated",
         color = "MS disease status") +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors)
)

# perform binomial logistic regression correlating 5mC methylation at each site to
# MS status
m_subset %>%
  glm(cbind(n_mod,n_unmod) ~ disease_status + batch + sex, family = "binomial", data = .) %>%
  summary(.)

m_subset %>%
  glm(cbind(n_mod,n_unmod) ~ group + batch + sex, family = "binomial", data = .) %>%
  summary(.)

df_list <- list(h_subset, m_subset)
names <- c("5hmc", "5mc")

# iterate over 5mc and 5hmc subsets
for (i in seq_along(df_list)){
  tests <- list()
  tests_grouped <- list()
  # iterate over each differentially methylated position
  for (genom_position in unique(df_list[[i]]$position)) {
    df <- filter(df_list[[i]], position == genom_position)
    # check if each group is represented at the given positions
    if (length(unique(df$group)) == 3 & length(unique(df$disease_status)) == 2){
      # check if more than one sex and batch is included at the given position
      if (length(unique(df$sex)) > 1 & length(unique(df$batch)) > 1){
        # perform binomial regression with MS group, with sex and batch as covariables
        glm <- glm(cbind(n_mod, n_unmod) ~ group + sex + batch, family = binomial(logit),
                   data = df,
                   method = "brglmFit", type = "AS_mean")
        tests[[genom_position]] <- glm
        # perform binomial regression with disease_status, with sex and batch as covariables
        glm_grouped <- glm(cbind(n_mod, n_unmod) ~ disease_status + sex + batch, family = binomial(logit),
                           data = df,
                           method = "brglmFit", type = "AS_mean")
        tests_grouped[[genom_position]] <- glm_grouped
      } else{
        # perform binomial regression with MS group, with no covariables
        glm <- glm(cbind(n_mod, n_unmod) ~ group, family = binomial(logit),
                   data = df,
                   method = "brglmFit", type = "AS_mean")
        # perform binomial regression with disease_status, with no covariables
        glm_grouped <- glm(cbind(n_mod, n_unmod) ~ disease_status, family = binomial(logit), 
                           data = df,
                           method = "brglmFit", type = "AS_mean")
        
        tests[[genom_position]] <- glm
        tests_grouped[[genom_position]] <- glm_grouped
      } 
    }
  }
  # extract position, gene, coefficients and p-values from tests
  tests_df_pre <- map_dfr(names(tests), function(x) {tibble(position = x,
                                                            gene = unlist(str_split(unlist(str_split(x, "\\|"))[4], "_"))[1],
                                                            low_activity_coefficient = coef(summary(tests[[x]]))[2],
                                                            low_activity_pval = coef(summary(tests[[x]]))[length(coef(summary(tests[[x]])))-((length(coef(summary(tests[[x]])))/4)-2)],
                                                            high_activity_coefficient = coef(summary(tests[[x]]))[3],
                                                            high_activity_pval = coef(summary(tests[[x]]))[length(coef(summary(tests[[x]])))-((length(coef(summary(tests[[x]])))/4)-3)])
                      })
  # perform fdr correction of p-values
  tests_df <- tests_df_pre %>%
    mutate(adj_low_activity_pval = p.adjust(low_activity_pval, method = "fdr"),
           adj_high_activity_pval = p.adjust(high_activity_pval, method = "fdr"))
  # save tests-objects and dataframe
  new_df <- paste0("result_", names[i])
  new_tests <- paste0("tests_", names[i])
  assign(new_df, tests_df)
  assign(new_tests, tests)
  # extract position, gene, coefficients and p-values from tests
  tests_grouped_df <- map_dfr(names(tests_grouped), function(x) {tibble(position = x, 
                                                                        gene = unlist(str_split(unlist(str_split(x, "\\|"))[4], "_"))[1],
                                                                        coefficient = coef(summary(tests_grouped[[x]]))[2], 
                                                                        pval = coef(summary(tests_grouped[[x]]))[length(coef(summary(tests_grouped[[x]])))-((length(coef(summary(tests_grouped[[x]])))/4)-2)])
                              })
  # perform fdr correction of p-values
  tests_grouped_df <- tests_grouped_df %>%
    mutate(adj_pval = p.adjust(pval, method = "fdr"))
  # save tests-objects and dataframe
  new_df <- paste0("result_grouped_", names[i])
  new_tests <- paste0("tests_grouped_", names[i])
  assign(new_df, tests_grouped_df)
  assign(new_tests, tests_grouped)
}

# plot p-value distribution of tests
plotly::ggplotly(
ggplot(result_grouped_5hmc, aes(x=pval)) +
  geom_histogram(binwidth = 0.05,
                 boundary = 0,
                 fill="orange",
                 color="black") +
  labs(title = "P-value distribution of binomial regression tests for 5hmC",
       y = "Count",
       x = "P-values") +
  theme_bw()
)

# create volcano plot of test results
result_grouped_5mc
ggplot(result_grouped_5hmc, aes(x = coefficient, y = -log10(pval), color = ifelse(pval < 0.05, "Significant", "Not Significant"))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "orange")) +
  theme_minimal() +
  labs(
    title = "Test results for low activity group",
    x = "Coefficient Size",
    y = "-Log10(p-value)",
    color = "Significance 0.05"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  theme_bw() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = gene))
  #geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")

# investigate tests
result_5mc %>%
  arrange(low_activity_pval)

result_5hmc %>%
  arrange(high_activity_pval)

result_grouped_5hmc %>%
  arrange(pval)

result_grouped_5mc %>%
  arrange(pval) %>%
  print(n = 50)


## repeat above procedure and generate methylation plots for q-filtered data

base_path <- "L:/LovbeskyttetMapper/Gridion/August/thesis_data/runs"


run_folders <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)


target_folders <- run_folders[grepl("RUN[0-9]+", basename(run_folders))]

input_files <- c()
for (i in 1:length(target_folders)) {
  subfolders <- list.dirs(target_folders[i], full.names = TRUE, recursive = TRUE)
  target_subfolder <- subfolders[grepl("q_filtered_q7/methylation_data", subfolders)]
  input_files <- c(input_files, list.files(target_subfolder, pattern = "*.bed.gz$", full.names = TRUE))
}

data_dict <- setNames(input_files, basename(input_files))

# read data and add sample_name column
df_list <- lapply(names(data_dict), function(key) {
  df <- read_delim(data_dict[[key]], delim = "\t", col_names = FALSE,
                   show_col_types = FALSE)
  df$sample_name <- str_extract(key, "^[^_.]+")
  return(df)
})

# merge dataframes
all_data <- bind_rows(df_list)

subset_data <- all_data %>% 
  # only use positions with depth of 8 or greater
  filter(X10 >= 6) %>%
  # concatenate position information
  mutate(position = str_c(X1, X2, X3, sep = "_"),
         group = group_mapping[sample_name],
         sex = sex_mapping[sample_name],
         batch = batch_mapping[sample_name],
         n_mod = X12,
         n_unmod = X10-X12) %>%
  mutate(disease_status = ifelse(group == "symptomatic_control", "symptomatic_control", "multiple_sclerosis")) %>%
  mutate(group = factor(group, levels = c("symptomatic_control", "low_activity", "high_activity")),
         disease_status = factor(disease_status, levels = c("symptomatic_control", "multiple_sclerosis")),
         sex = factor(sex, levels = c("female", "male")),
         batch = factor(batch, levels = c("RUN2", "RUN3")))

h_subset <- subset_data %>%
  filter(X4 == "h") %>%
  select(c(position, sample_name, group, sex, batch, disease_status, n_mod, n_unmod)) 

m_subset <- subset_data %>%
  filter(X4 == "m" & n_mod >= 1) %>%
  select(c(position, sample_name, group, sex, batch, disease_status, n_mod, n_unmod))

plotly::ggplotly(
  ggplot(m_subset, aes(x = disease_status, y = n_mod/(n_mod+n_unmod), 
                       color = disease_status, fill = disease_status)) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +  # Jittered points
    theme_minimal() +  # Clean theme
    theme(legend.position = "none") +
    labs(title = "5mC methylated fraction at genomic positions with detected methylation",
         x = NULL,
         y = "Fraction methylated",
         color = "MS disease status") +
    scale_color_manual(values = custom_colors[c(1,3)]) +
    scale_fill_manual(values = custom_colors[c(1,3)])
)

plotly::ggplotly(
  ggplot(h_subset, aes(x = disease_status, y = n_mod/(n_mod+n_unmod), 
                       color = disease_status, fill = disease_status)) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +  # Jittered points
    theme_minimal() +  # Clean theme
    theme(legend.position = "none") +
    labs(title = "5hmC methylated fraction at genomic positions with detected methylation",
         x = NULL,
         y = "Fraction methylated",
         color = "MS disease status") +
    scale_color_manual(values = custom_colors[c(1,3)]) +
    scale_fill_manual(values = custom_colors[c(1,3)])
)
