
# Differences in demographic variables between subjects used for LEiDA analysis and those who were discarded"
# Code adapted by Isabella L.C. Mariani Wigley and Aurora Berto  (06 / 2025).

# The aim of this code is to see whether there is any difference between distributions of demographic 
# variables used in the main article between subjects involved in LEiDA analysis, and in the subsequent statistical analysis, 
# and those who were discarded in the previous steps.

library(dplyr)
library(ggplot2)

################################################################################
# Subjects passed to LEiDA analysis
df_path <- "/path/to/dataset"
df_file <- file.path(df_path,"main_dataset_rsfMRI_FCH_LEiDA_demo.csv")
df <- read.csv(df_file, header = T, stringsAsFactors = F)

# Create the dataset with the variables involved in the correlation analysis for subjects passed to LEiDA
df_leida_analysis <- df %>%
  dplyr::select(src_subject_id, demo_sex_v2, race_ethnicity, interview_age, 
                pubertal_developmental_scale, rel_family_id, 
                rel_group_id, triponderal_mass_index, mri_info_visitid)

# 777 = refuse to answer, 999 = don't know set to NA
df_leida_analysis[] <- lapply(df_leida_analysis, function(x) {
  if (is.numeric(x)) x[x %in% c(777, 999)] <- NA
  x
})

# Drop duplicates based on rel_family_id and rel_group_id, keeping the first occurrence
df_leida_analysis <- df_leida_analysis %>%
  distinct(rel_family_id, rel_group_id, .keep_all = TRUE)

# Convert demo_sex_v2 to factor variable
df_leida_analysis$demo_sex_v2 <- factor(df_leida_analysis$demo_sex_v2, levels=c(1,2), labels=c("M","F"))

# Convert race_ethnicity to factor variable
df_leida_analysis$race_ethnicity <- factor(df_leida_analysis$race_ethnicity, levels = c(1,2,3,4,5),
                                           labels = c("W","B","H","A","other"))

# Extract site variable
df_leida_analysis$site <- sub("_.*", "", df_leida_analysis$mri_info_visitid)

# Remove columns from the dataset
df_leida_analysis <- df_leida_analysis %>%
  dplyr::select(-rel_family_id, -rel_group_id, -mri_info_visitid)


# Make a list with all subjects passed to LEiDA
leida_subjs <- df$src_subject_id

################################################################################
# Subjects discarded from LEiDA analysis
# Load all covariates and merge in the same dataset based on subjects ID, retaining only 
# values for subjects that did not go to LEiDA analysis and that were acquired at baseline.


### Demographic variables
# Define path
file_path <- "/path/to/abcd-data-release-5.0/core/abcd-general/abcd_p_demo.csv"

# Read file
raw_file <- read.csv(file_path)

# Define variables to keep
vars_to_keep <- c("demo_sex_v2","race_ethnicity","src_subject_id")

# Retain only columns to keep -> acquisitions at baseline and subjects that did not do LEiDA
df_demo <- raw_file[raw_file$eventname == "baseline_year_1_arm_1",vars_to_keep, drop=F]
df_demo <- df_demo[!df_demo$src_subject_id %in% df_leida_analysis$src_subject_id, ]

# Convert to categorical variables
df_demo$demo_sex_v2 <- factor(df_demo$demo_sex_v2, levels=c(1,2), labels=c("M","F"))
df_demo$race_ethnicity <- factor(df_demo$race_ethnicity, levels = c(1,2,3,4,5),
                                 labels = c("W","B","H","A","other"))
# Remove variables
rm(file_path, raw_file, vars_to_keep)


### Interview age
# Define path
file_path <- "/path/to/abcd-data-release-5.0/core/abcd-general/abcd_y_lt.csv"

# Read file
raw_file <- read.csv(file_path)

# Define variables to keep
vars_to_keep <- c("interview_age","src_subject_id")

# Retain only columns to keep
df_age <- raw_file[raw_file$eventname == "baseline_year_1_arm_1",vars_to_keep, drop=F]
df_age <- df_age[!df_age$src_subject_id %in% df_leida_analysis$src_subject_id, ]

# Remove variables
rm(file_path, raw_file, vars_to_keep)

# Merge with the previous one based on subjects ID
df_NONleida_analysis <- merge(df_demo, df_age, by = "src_subject_id", keep.all = TRUE)


### Puberty - parents reported
# Define path
file_path_parents <- "/path/to/abcd-data-release-5.0/core/physical-health/ph_p_pds.csv"

# Read file
file_parents <- read.csv(file_path_parents)

# Define variables to keep
vars_to_keep <- c("src_subject_id", "pds_p_ss_female_category", "pds_p_ss_male_category")

# Retain only columns to keep
cols_to_keep_parents <- file_parents[file_parents$eventname == "baseline_year_1_arm_1",vars_to_keep, drop=F]
cols_to_keep_parents <- cols_to_keep_parents[!cols_to_keep_parents$src_subject_id %in%
                                               df_leida_analysis$src_subject_id, ]

### Puberty - young reported
# Define path
file_path_youngs <- "/path/to/abcd-data-release-5.0/core/physical-health/ph_y_pds.csv"

# Read file
file_youngs <- read.csv(file_path_youngs)

# Define variables to keep
vars_to_keep <- c("src_subject_id", "pds_y_ss_female_category", "pds_y_ss_male_category")

# Retain only columns to keep
cols_to_keep_youngs <- file_youngs[file_youngs$eventname == "baseline_year_1_arm_1",vars_to_keep, drop=F]
cols_to_keep_youngs <- cols_to_keep_youngs[!cols_to_keep_youngs$src_subject_id %in%
                                             df_leida_analysis$src_subject_id, ]


# Merge the two datasets on the common key (subjects' ID)
df_merged <- merge(cols_to_keep_parents, cols_to_keep_youngs, by = "src_subject_id", keep.all = TRUE)

# Calculate means dividing by sex
df_merged$pubertal_dev_score <- rowMeans(
  cbind(
    ifelse(!is.na(df_merged$pds_p_ss_female_category), df_merged$pds_p_ss_female_category, df_merged$pds_p_ss_male_category),
    ifelse(!is.na(df_merged$pds_y_ss_female_category), df_merged$pds_y_ss_female_category, df_merged$pds_y_ss_male_category)
  ),
  na.rm = TRUE
)

# Final dataset
df_pds <- data.frame(pubertal_dev_score = df_merged$pubertal_dev_score, 
                     src_subject_id = df_merged$src_subject_id)

# Remove variables
rm(file_path_parents, file_path_youngs, file_parents, file_youngs, vars_to_keep, cols_to_keep_parents, cols_to_keep_youngs)

# Merge with the previous one based on subjects ID
df_NONleida_analysis <- merge(df_NONleida_analysis, df_pds, by = "src_subject_id", keep.all = TRUE)


### MRI device information
# Define path
file_path <- "/path/to/abcd-data-release-5.0/core/imaging/mri_y_adm_info.csv"

# Read file
raw_file <- read.csv(file_path)

# Define variables to keep
vars_to_keep <- c("mri_info_visitid","src_subject_id")

# Retain only columns to keep
df_site <- raw_file[raw_file$eventname == "baseline_year_1_arm_1",vars_to_keep, drop=F]
df_site <- df_site[!df_site$src_subject_id %in% df_leida_analysis$src_subject_id, ]

# Retain only site information and insert it in the dataset
df_site$site <- sub("_.*", "", df_site$mri_info_visitid)
df_site <- df_site %>%
  dplyr::select(-mri_info_visitid)

# Remove variables
rm(file_path, raw_file, vars_to_keep)

# Merge with the previous one based on subjects ID
df_NONleida_analysis <- merge(df_NONleida_analysis, df_site, by = "src_subject_id", keep.all = TRUE)


## Triponderal Mass Index (TMI)
# Define path
file_path <- "/path/to/abcd-data-release-5.0/core/physical-health/ph_y_anthro.csv"

# Read file
raw_file <- read.csv(file_path)

# Define variables to keep
vars_to_keep <- c("anthro_1_height_in","anthroweight1lb","src_subject_id")

# Retain only columns to keep
df_tmi <- raw_file[raw_file$eventname == "baseline_year_1_arm_1",vars_to_keep, drop=F]
df_tmi <- df_tmi[!df_tmi$src_subject_id %in% df_leida_analysis$src_subject_id, ]

# Define function to calculate tmi
convert_lbs_to_kg <- function(weight_lbs) { # Convert lbs to kg
  lbs_to_kg_ratio <- 0.45359237
  weight_kg <- weight_lbs * lbs_to_kg_ratio
  return(weight_kg)
}

convert_in_to_m <- function(height_in) { # Convert in to m
  in_to_m_ratio <- 0.0254
  height_m <- height_in * in_to_m_ratio
  return(height_m)
}

compute_tmi <- function(weight_lbs, height_in) { # Compute tmi
  weight_kg <- convert_lbs_to_kg(weight_lbs)
  height_m <- convert_in_to_m(height_in)
  
  # If height is missing or is zero return NA
  if (is.na(height_m) || height_m == 0) {
    return(NA)
  }
  
  tmi <- weight_kg / (height_m^3)
  return(tmi)
}

# Calculate TMI
df_tmi$triponderal_mass_index <- mapply(compute_tmi, df_tmi$anthroweight1lb, df_tmi$anthro_1_height_in)
df_tmi <- df_tmi %>%
  dplyr::select(-anthro_1_height_in, -anthroweight1lb)

# Keep only subjects with 1.5*1iqr < TMI < 3*3iqr
quantile(df_tmi$triponderal_mass_index, na.rm = TRUE)
Q1 <- 11.434366
Q3 <- 14.755722
IQR <- Q3 - Q1
df_tmi <- df_tmi[df_tmi$triponderal_mass_index >= (Q1 - 1.5*IQR) & df_tmi$triponderal_mass_index <= (Q3 + 3*IQR), ]

# Remove variables
rm(file_path, raw_file, vars_to_keep)

# Merge with the previous one based on subjects ID
df_NONleida_analysis <- merge(df_NONleida_analysis, df_tmi, by = "src_subject_id", keep.all = TRUE)

# Lose the pilot site
df_NONleida_analysis <- df_NONleida_analysis[df_NONleida_analysis$site != "G054",]

# Remove all NA rows
df_NONleida_analysis <- na.omit(df_NONleida_analysis) # 7058

# Remove unnecessary datasets
remove(df_age, df_demo, df_merged, df_pds, df_site, df_tmi, df_file, df_path)

################################################################################
# t-tests
## Continuous variables

### Interview age
# Combine the two dataset based on the variable to control
df_plot <- rbind(
  data.frame(value = df_leida_analysis$interview_age, group = "Included sample"),
  data.frame(value = df_NONleida_analysis$interview_age, group = "Excluded sample")
)
df_plot <- df_plot %>% filter(is.finite(value))

# Density plot 
ggplot(df_plot, aes(x = value, fill = group)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  scale_fill_manual(values = c("Included sample" = "#1F3A5A", "Excluded sample" = "#89B3D9")) +
  labs(title = "Interview age distribution")

# Statistic test: t-test
t.test(df_leida_analysis$interview_age, df_NONleida_analysis$interview_age)

# Mean
mean(df_leida_analysis$interview_age)
mean(df_NONleida_analysis$interview_age)

# Standard deviation
sd(df_leida_analysis$interview_age)
sd(df_NONleida_analysis$interview_age)


### PDS
# Combine the two dataset based on the variable to control
df_plot <- rbind(
  data.frame(value = df_leida_analysis$pubertal_developmental_scale, group = "Included sample"),
  data.frame(value = df_NONleida_analysis$pubertal_dev_score, group = "Excluded sample")
)
df_plot <- df_plot %>% filter(is.finite(value))

# Density plot 
ggplot(df_plot, aes(x = value, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), position = "dodge", bins = 20, color = "black") +
  theme_minimal() +
  scale_fill_manual(values = c("Included sample" = "#66304A", "Excluded sample" = "#E7A3B8")) +
  labs(title = "PDS distribution", y = "density")

# Statistic test: t-test
t.test(df_leida_analysis$pubertal_developmental_scale, df_NONleida_analysis$pubertal_dev_score)

# Mean
median(df_leida_analysis$pubertal_developmental_scale)
median(df_NONleida_analysis$pubertal_dev_score)

# Standard deviation
sd(df_leida_analysis$pubertal_developmental_scale)
sd(df_NONleida_analysis$pubertal_dev_score)


### TMI
# Combine the two dataset based on the variable to control
df_plot <- rbind(
  data.frame(value = df_leida_analysis$triponderal_mass_index, group = "Included sample"),
  data.frame(value = df_NONleida_analysis$triponderal_mass_index, group = "Excluded sample")
)

df_plot <- df_plot %>% filter(is.finite(value))

# Density plot 
ggplot(df_plot, aes(x = value, fill = group)) +
  geom_density(alpha = 0.4) +
  #scale_x_continuous(limits = c(1,50))
  scale_fill_manual(values = c("Included sample" = "#2D3B20", "Excluded sample" = "#A2C199")) +
  theme_minimal() +
  labs(title = "TMI distribution")

# Statistic test: t-test
t.test(df_leida_analysis$triponderal_mass_index, df_NONleida_analysis$triponderal_mass_index)

# Mean
mean(df_leida_analysis$triponderal_mass_index)
mean(df_NONleida_analysis$triponderal_mass_index)

# Standard deviation
sd(df_leida_analysis$triponderal_mass_index)
sd(df_NONleida_analysis$triponderal_mass_index)


## Categorical variables
### Sex
# Combine the two dataset based on the variable to control
df_plot <- rbind(
  data.frame(value = df_leida_analysis$demo_sex_v2, group = "Included sample"),
  data.frame(value = df_NONleida_analysis$demo_sex_v2, group = "Excluded sample")
)

# Calculate proportions per group
df_prop <- df_plot %>%
  group_by(group, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(prop = n / sum(n))

# Plot
ggplot(df_prop, aes(x = value, y = prop, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Included sample" = "#4B2466", "Excluded sample" = "#C3A0D9")) +
  theme_minimal() +
  labs(title = "Sex distribution (percentage)",
       x = "Sex", y = "Percentage")


# Build contingency table
sex_table <- table(
  c(df_leida_analysis$demo_sex_v2, df_NONleida_analysis$demo_sex_v2),
  c(rep("LEiDA", nrow(df_leida_analysis)), rep("non-LEiDA", nrow(df_NONleida_analysis)))
)

# Chi-squared test
print(sex_table)
chisq.test(sex_table)
prop.table(sex_table, margin = 2) 

### Ethnicity
# Combine the two dataset based on the variable to control
df_plot <- rbind(
  data.frame(value = df_leida_analysis$race_ethnicity, group = "Included sample"),
  data.frame(value = df_NONleida_analysis$race_ethnicity, group = "Excluded sample")
)

# Calculate proportions per group
df_prop <- df_plot %>%
  group_by(group, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(prop = n / sum(n))

ggplot(df_prop, aes(x = value, y = prop, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Included sample" = "#7A3F17", "Excluded sample" = "#F3A974")) +
  theme_minimal() +
  labs(title = "Ethnicity distribution (percentage)",
       x = "Ethnicity", y = "Percentage")

# Build contingency table
ethnicity_table <- table(
  c(df_leida_analysis$race_ethnicity, df_NONleida_analysis$race_ethnicity),
  c(rep("LEiDA", nrow(df_leida_analysis)), rep("non-LEiDA", nrow(df_NONleida_analysis)))
)

# Chi-squared test
print(ethnicity_table)
chisq.test(ethnicity_table)
prop.table(ethnicity_table, margin = 2) 

### Site
# Combine the two dataset based on the variable to control
df_plot <- rbind(
  data.frame(value = df_leida_analysis$site, group = "Included sample"),
  data.frame(value = df_NONleida_analysis$site, group = "Excluded sample")
)
df_plot$value <- as.factor(df_plot$value)

# Calculate proportions per group
df_prop <- df_plot %>%
  group_by(group, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(prop = n / sum(n))

ggplot(df_prop, aes(x = value, y = prop, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Included sample" = "#67322D", "Excluded sample" = "#E89A8B")) +
  theme_minimal() +
  labs(title = "Site distribution (percentage)",
       x = "Site", y = "Percentage")

# Build contingency table
site_table <- table(
  c(df_leida_analysis$site, df_NONleida_analysis$site),
  c(rep("LEiDA", nrow(df_leida_analysis)), rep("non-LEiDA", nrow(df_NONleida_analysis)))
)

# Chi-squared test
print(site_table)
chisq.test(site_table)
prop.table(site_table, margin = 2) 
