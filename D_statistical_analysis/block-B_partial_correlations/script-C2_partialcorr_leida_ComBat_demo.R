# Partial correlation analysis: LEiDA harmonized and demographics

# Code to run partial correlation analysis between LEiDA harmonized outcomes and demographics.
# Adapted by Isabella L.C. Mariani Wigley (ilmawi@utu.fi) and Aurora Berto (aurber@utu.fi)
# 05 / 2025

library(dplyr)
library(stringr)
library(ppcor)
library(ggplot2)
library(ggrepel)
library(colorspace)
library(ggdist)
library(patchwork)
library(tidyr)

################################################################################
## Load the dataset
# To change accordingly to your path
df_path <- "/path/to/dataset/folder"
df_file <- file.path(df_path,"main_dataset_rsfMRI_FCH_LEiDA_demo.csv")
df <- read.csv(df_file, header = T, stringsAsFactors = T)

## Load the dataset
# To change accordingly to your path
df_leida_path <- "/path/to/neuroCombat_harmonization/results"
df_leida_file <- file.path(df_leida_path,"dataset_leida_harmonized.csv")
df_leida <- read.csv(df_leida_file, header = T, stringsAsFactors = T)

################################################################################
## Prepare dataset for the analysis

## Interview age, sex, pubertal developmental score, weight at birth and triponderal mass index
df_demo <- df %>%
  dplyr::select(src_subject_id, demo_sex_v2, race_ethnicity, interview_age, 
                pubertal_developmental_scale, rel_family_id, 
                rel_group_id, triponderal_mass_index)

# 777 = refuse to answer, 999 = don't know set to NA
df_demo[] <- lapply(df_demo, function(x) {
  if (is.numeric(x)) x[x %in% c(777, 999)] <- NA
  x
})

# Drop duplicates based on rel_family_id and rel_group_id, keeping the first occurrence
df_demo <- df_demo %>%
  distinct(rel_family_id, rel_group_id, .keep_all = TRUE)

# Remove columns from the dataset
df_demo <- df_demo %>%
  dplyr::select(-rel_family_id, -rel_group_id)


### Discard LEiDA metrics with no correspondence to Yeo et al. rs-networks.
# LEiDA columns
leida_cols <- colnames(df_leida)

# Probabilities and lifetimes to remove
leida_to_remove <- c(
  "P_k5c5","P_k10c5", "P_k12c5", "P_k15c7", "P_k16c9", "P_k18c11", "P_k19c11", "P_k15c12",
  "P_k20c12", "P_k14c13", "P_k17c13", "P_k20c14", "P_k15c15", "P_k20c17", "P_k19c19", 
  "P_k20c20",
  
  "LT_k5c5", "LT_k10c5", "LT_k12c5", "LT_k15c7", "LT_k16c9", "LT_k18c11", "LT_k19c11",
  "LT_k15c12", "LT_k20c12", "LT_k14c13", "LT_k17c13", "LT_k20c14", "LT_k15c15", "LT_k20c17",
  "LT_k19c19", "LT_k20c20"
)

# Extract couples (k, c) from probabilities and lifetimes to remove
k_c_pairs <- str_match(leida_to_remove, "k(\\d+)c(\\d+)")
k_c_pairs <- k_c_pairs[!is.na(k_c_pairs[,1]), 2:3]
k_c_pairs <- as.data.frame(k_c_pairs, stringsAsFactors = FALSE)
k_c_pairs <- mutate(k_c_pairs, V1 = as.integer(V1), V2 = as.integer(V2))

# Find transitions to remove
transitions_to_remove <- sapply(leida_cols, function(col) {
  m <- str_match(col, "^TR_K(\\d+)_C(\\d+)x(\\d+)")
  if (any(!is.na(m))) {
    k <- as.integer(m[2])
    c1 <- as.integer(m[3])
    c2 <- as.integer(m[4])
    any(k_c_pairs$V1 == k & (k_c_pairs$V2 == c1 | k_c_pairs$V2 == c2))
  } else {
    FALSE
  }
})
transitions_to_remove <- leida_cols[transitions_to_remove]

# Merge everything to remove
all_to_remove <- union(leida_to_remove, transitions_to_remove)

# Filter only columns to keep
leida_cols <- setdiff(leida_cols, all_to_remove)

# Final LEiDA dataset with only columns to keep
df_leida <- df_leida[,leida_cols]


## Merge all datasets in a single one
data <- merge(df_demo, df_leida, by = "src_subject_id")

################################################################################
# LEiDA probabilities

### Partial correlation analysis

## Divide in two datasets
X_data <- data %>%
  dplyr::select(demo_sex_v2, race_ethnicity, interview_age, 
                pubertal_developmental_scale, triponderal_mass_index)
names_covariates <- colnames(X_data)

## Extract LEiDA probabilities
y_data <- df_leida %>%
  dplyr::select(matches("P_k"))

# Initialize list to hold full correlation matrices
full_results <- list()

for (K in 2:20){
  
  leida_k <- y_data %>%
    dplyr::select(matches(paste0("P_k",K,"c"))) %>%
    dplyr::select(any_of(leida_cols)) # select only those present in leida_cols
  
  K_results <- list()
  
  for (s in 1:K){
    
    leida_state <- paste0("P_k",K,"c",s)
    
    # Skip if current state has been excluded
    if (!(leida_state %in% leida_cols)) next
    
    X_tmp <- X_data
    X_tmp$tmp_state <- leida_k[[leida_state]]
    X_tmp <- na.omit(X_tmp)
    
    # Perform partial correlation
    test <- pcor(X_tmp[, c("tmp_state", names_covariates)], method = "pearson")
    
    # Extract correlations and p-values between tmp_state and all covariates
    cors <- test$estimate["tmp_state", ]
    pvals <- test$p.value["tmp_state", ]
    
    K_results[[s]] <- data.frame(
      K = K,
      state = s,
      variable = names(cors),
      correlation = cors,
      pvalue = pvals,
      stringsAsFactors = FALSE
    )
  }
  
  full_results[[paste0("K",K)]] <- do.call(rbind, K_results)
}

# Merge all results in a unique dataset
df_all_corrs <- bind_rows(full_results)

# Add fdr correction and variable's name
df_all_corrs <- df_all_corrs %>%
  group_by(K, variable) %>%
  mutate(pvalue_fdr = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()

# ==============================================================================
### Visualization
## Create a label with correlations with p < 0.05
df_all_corrs <- df_all_corrs %>%
  filter(variable != "tmp_state") %>%
  mutate(
    label = ifelse(pvalue_fdr < 0.05,
                   state,
                   NA),
    K = factor(K)
  )

## Rename variables for plot visualization
pretty_names <- c(
  "demo_sex_v2" = "Sex",
  "interview_age" = "Age",
  "pubertal_developmental_scale" = "Pubertal Development",
  "race_ethnicity" = "Ethnicity",
  "triponderal_mass_index" = "TMI",
  "site" = "Site"
)

# Add column with readable names
df_all_corrs <- df_all_corrs %>%
  mutate(variable_pretty = recode(variable, !!!pretty_names))

## Define colors
df_all_corrs <- df_all_corrs %>%
  mutate(variable_colored = ifelse(pvalue_fdr < 0.05, variable_pretty, "nonsignificant"))

# Create categories variable+sign
df_all_corrs <- df_all_corrs %>%
  mutate(
    variable_colored = case_when(
      pvalue_fdr < 0.05 & correlation >= 0 ~ paste0(variable_pretty, "_pos"),
      pvalue_fdr < 0.05 & correlation < 0 ~ paste0(variable_pretty, "_neg"),
      TRUE ~ "nonsignificant"
    )
  )

# Define colors dark/light for each variable
pretty_colors <- c(
  "Sex_pos" = "#C79FDB",
  "Sex_neg" = "#7C3E8B",
  "TMI_pos" = "#9EC39C",
  "TMI_neg" = "#4F7450",
  "Ethnicity_pos" = "#E29C5F",
  "Ethnicity_neg" = "#A85F2B",
  "Age_pos" = "#6F99C2",
  "Age_neg" = "#2C5172",
  "Pubertal Development_pos" = "#D7839C",
  "Pubertal Development_neg" = "#8C3D56",
  "Site_pos" = "#D98B78",
  "Site_neg" = "#8B4A39",
  "nonsignificant" = "black"
)

ggplot(df_all_corrs, aes(x = K, y = pvalue)) +
  
  geom_jitter(data = df_all_corrs[is.na(df_all_corrs$label), ],
              aes(color = variable_colored, size = abs(correlation)),
              width = 0.25, alpha = 0.8) +
  
  geom_text_repel(data = df_all_corrs[!is.na(df_all_corrs$label), ],
                  aes(label = label, color = variable_colored),
                  size = 3.5, fontface = "bold",
                  direction = "y",
                  max.overlaps = Inf, min.segment.length = 0,
                  box.padding = 0.2, point.padding = 0.2, segment.color = NA) +
  
  scale_y_log10(limits = c(1e-6, 1),
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)) +
  
  scale_color_manual(values = pretty_colors, guide = "none") +
  scale_size_continuous(range = c(0.5, 2), guide = "none") +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
  
  labs(x = "Partition model K",
       y = "p-values",
       color = "Variable",
       title = "Partial Correlation between LEiDA probabilities and covariates",
       subtitle = "Significance determined by FDR-corrected p-values < 0.05") +
  
  facet_wrap(~ variable_pretty, scales = "free_y") +
  
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.position = "bottom"
  )

# ==============================================================================
### Half-violin plot
## Race-ethnicity
# Filter only the desired variable
df_results_ethn <- df_all_corrs[df_all_corrs$variable == "race_ethnicity",]
df_results_ethn <- df_results_ethn[!is.na(df_results_ethn$label),]

# List to save all graphs
plot_list <- list()

for (i in 1:nrow(df_results_ethn)) {
  
  leida_state <- paste0("P_k", df_results_ethn$K[i], "c", df_results_ethn$state[i])
  tmp_state <- y_data[[leida_state]]
  
  # Temporary dataset
  X_tmp <- data.frame(
    ethn = factor(X_data$race_ethnicity, labels = c("W", "B", "H", "A", "other")),
    probability = tmp_state
  )
  
  # Half violin plot
  p <- ggplot(X_tmp, aes(x = ethn, y = probability, fill = ethn)) +
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4) +
    # geom_jitter(width = 0.05, size = 0.5, alpha = 0.5) +
    scale_fill_manual(values = c("#F3A974", "#FACDA8", "#F7BA8D", "#FDDDBC", "#E18A5E")) +
    #coord_cartesian(ylim = c(0, 0.7)) +
    labs(title = paste("LEiDA:", leida_state), x = "Ethnicity", y = "Probability") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Combine all graphs in the same figure
# Specify the desired number of columns
n_col <- 4
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
print(combined_plot)

# ------------------------------------------------------------------------------
## Sex 
# Filter only the desired variable
df_results_sex <- df_all_corrs[df_all_corrs$variable == "demo_sex_v2",]
df_results_sex <- df_results_sex[!is.na(df_results_sex$label),]

# List to save all graphs
plot_list <- list()

for (i in 1:nrow(df_results_sex)) {
  
  leida_state <- paste0("P_k", df_results_sex$K[i], "c", df_results_sex$state[i])
  tmp_state <- y_data[[leida_state]]
  
  # Temporary dataset
  X_tmp <- data.frame(
    sex = factor(X_data$demo_sex_v2, labels = c("M", "F")),
    probability = tmp_state
  )
  
  # Half violin plot
  p <- ggplot(X_tmp, aes(x = sex, y = probability, fill = sex)) +
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4) +
    # geom_jitter(width = 0.05, size = 0.5, alpha = 0.5) +
    scale_fill_manual(values = c("#D7B8EB", "#E3C9F2")) +
    coord_cartesian(ylim = c(0, 0.6)) +
    labs(title = paste("LEiDA:", leida_state), x = "Sex", y = "Probability") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Combine all graphs in the same figure
# Specify the desired number of columns
n_col <- 8  
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
print(combined_plot)

# ==============================================================================
### Scatter plot
## Interview age
# Filter only the desired variable
df_results_age <- df_all_corrs[df_all_corrs$variable == "interview_age",]
df_results_age <- df_results_age[!is.na(df_results_age$label),]

# List for plots
plot_list <- list()

for (i in 1:nrow(df_results_age)) {
  
  leida_state <- paste0("P_k", df_results_age$K[i], "c", df_results_age$state[i])
  tmp_state <- y_data[[leida_state]]
  
  # Temporary dataset
  X_tmp <- data.frame(
    age = as.numeric(X_data$interview_age), 
    probability = tmp_state
  )
  
  # Scatterplot
  p <- ggplot(X_tmp, aes(x = age, y = probability, color = age)) +
    geom_point(alpha = 0.6, size = 1.2) +
    # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    scale_color_gradient(low = "#BDDDF1", high = "#89B3D9") +
    coord_cartesian(ylim = c(0, 0.6)) +
    labs(title = paste("LEiDA:", leida_state),
         x = "Interview age",
         y = "Probability") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 14),  
          axis.title = element_text(size = 12),
          axis.text  = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Final combination
n_col <- 4
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
suppressMessages(print(combined_plot))

# ------------------------------------------------------------------------------
## PDS
# Filter only the desired variable
df_results_pds <- df_all_corrs[df_all_corrs$variable == "pubertal_developmental_scale",]
df_results_pds <- df_results_pds[!is.na(df_results_pds$label),]

# List for plots
plot_list <- list()

for (i in 1:nrow(df_results_pds)) {
  
  leida_state <- paste0("P_k", df_results_pds$K[i], "c", df_results_pds$state[i])
  tmp_state <- y_data[[leida_state]]
  
  # Temporary dataset
  X_tmp <- data.frame(
    pds = as.numeric(X_data$pubertal_developmental_scale), 
    probability = tmp_state
  )
  
  # Scatterplot
  p <- ggplot(X_tmp, aes(x = pds, y = probability, color = pds)) +
    geom_point(alpha = 0.6, size = 1.2) +
    # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    scale_color_gradient(low = "#E7A3B8", high = "pink") +
    coord_cartesian(ylim = c(0, 0.6)) +
    labs(title = paste("LEiDA:", leida_state),
         x = "Pubertal developmental scale",
         y = "Probability") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Final combination
n_col <- 2
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
suppressMessages(print(combined_plot))

# ==============================================================================
### Heatmap
## Interview age
df_results_all <- df_all_corrs[!is.na(df_all_corrs$label),]
df_results_var <- df_results_all[df_results_all$variable == "interview_age",]

# Avoid log(0)
df_results_var <- df_results_var %>%
  mutate(log_p = -log10(pvalue),
         log_p = ifelse(is.infinite(log_p), NA, log_p),
         correlation_label = sprintf("%.2f", correlation))

# Calculate dynamic range for p-value logarithmic scale
log_p_min <- min(df_results_var$log_p, na.rm = TRUE)
log_p_max <- max(df_results_var$log_p, na.rm = TRUE)

# Heatmap
ggplot(df_results_var, aes(x = factor(K), y = factor(label), fill = log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = correlation_label), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#F5F5F5", high = "#89B3D9",
    name = expression(-log[10](p)),
    limits = c(log_p_min, log_p_max),
    oob = scales::squish
  ) +
  labs(
    title = "Heatmap | Correlation between LEiDA probabilities and age",
    x = "K",
    y = "Connection (label)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

# ------------------------------------------------------------------------------
## PDS
df_results_all <- df_all_corrs[!is.na(df_all_corrs$label),]
df_results_var <- df_results_all[df_results_all$variable == "pubertal_developmental_scale",]

# Avoid log(0)
df_results_var <- df_results_var %>%
  mutate(log_p = -log10(pvalue),
         log_p = ifelse(is.infinite(log_p), NA, log_p),
         correlation_label = sprintf("%.2f", correlation))

# Calculate dynamic range for p-value logarithmic scale
log_p_min <- min(df_results_var$log_p, na.rm = TRUE)
log_p_max <- max(df_results_var$log_p, na.rm = TRUE)

# Heatmap
ggplot(df_results_var, aes(x = factor(K), y = factor(label), fill = log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = correlation_label), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#F5F5F5", high = "#E7A3B8",
    name = expression(-log[10](p)),
    limits = c(log_p_min, log_p_max),
    oob = scales::squish
  ) +
  labs(
    title = "Heatmap | Correlation between LEiDA probabilities and PDS",
    x = "K",
    y = "Connection (label)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

################################################################################
# LEiDA lifetimes

### Partial correlation analysis

## Divide in two datasets
X_data <- data %>%
  dplyr::select(demo_sex_v2, race_ethnicity, interview_age, 
                pubertal_developmental_scale, triponderal_mass_index)
names_covariates <- colnames(X_data)

## Extract LEiDA probabilities
y_data <- df_leida %>%
  dplyr::select(matches("LT_k"))

# Initialize list to hold full correlation matrices
full_results <- list()

for (K in 2:20){
  
  leida_k <- y_data %>%
    dplyr::select(matches(paste0("LT_k",K,"c"))) %>%
    dplyr::select(any_of(leida_cols)) # select only those present in leida_cols
  
  K_results <- list()
  
  for (s in 1:K){
    
    leida_state <- paste0("LT_k",K,"c",s)
    
    # Skip if current state has been excluded
    if (!(leida_state %in% leida_cols)) next
    
    X_tmp <- X_data
    X_tmp$tmp_state <- leida_k[[leida_state]]
    X_tmp <- na.omit(X_tmp)
    
    # Perform partial correlation
    test <- pcor(X_tmp[, c("tmp_state", names_covariates)], method = "pearson")
    
    # Extract correlations and p-values between tmp_state and all covariates
    cors <- test$estimate["tmp_state", ]
    pvals <- test$p.value["tmp_state", ]
    
    K_results[[s]] <- data.frame(
      K = K,
      state = s,
      variable = names(cors),
      correlation = cors,
      pvalue = pvals,
      stringsAsFactors = FALSE
    )
  }
  full_results[[paste0("K",K)]] <- do.call(rbind, K_results)
}
# Merge all results in a unique dataset
df_all_corrs <- bind_rows(full_results)

# Add fdr correction and variable's name
df_all_corrs <- df_all_corrs %>%
  group_by(K, variable) %>%
  mutate(pvalue_fdr = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()

# ==============================================================================
### Visualization
## Create a label with correlations with p < 0.05
df_all_corrs <- df_all_corrs %>%
  filter(variable != "tmp_state") %>%
  mutate(
    label = ifelse(pvalue_fdr < 0.05,
                   state,
                   NA),
    K = factor(K)
  )

## Rename variables for plot visualization
pretty_names <- c(
  "demo_sex_v2" = "Sex",
  "interview_age" = "Age",
  "pubertal_developmental_scale" = "Pubertal Development",
  "race_ethnicity" = "Ethnicity",
  "site" = "Site",
  "triponderal_mass_index" = "TMI"
)

# Add column with readable names
df_all_corrs <- df_all_corrs %>%
  mutate(variable_pretty = recode(variable, !!!pretty_names))

## Define colors
df_all_corrs <- df_all_corrs %>%
  mutate(variable_colored = ifelse(pvalue_fdr < 0.05, variable_pretty, "nonsignificant"))

# Create categories variable+sign
df_all_corrs <- df_all_corrs %>%
  mutate(
    variable_colored = case_when(
      pvalue_fdr < 0.05 & correlation >= 0 ~ paste0(variable_pretty, "_pos"),
      pvalue_fdr < 0.05 & correlation < 0 ~ paste0(variable_pretty, "_neg"),
      TRUE ~ "nonsignificant"
    )
  )

# Define colors dark/light for each variable
pretty_colors <- c(
  "Sex_pos" = "#D0A5E0",
  "Sex_neg" = "#844191",
  "TMI_pos" = "#A7CEA6",
  "TMI_neg" = "#4B6F4B",
  "Ethnicity_pos" = "#F0B27B",
  "Ethnicity_neg" = "#B56533",
  "Age_pos" = "#81B0D8",
  "Age_neg" = "#2E5C8A",
  "Pubertal Development_pos" = "#E898AC",
  "Pubertal Development_neg" = "#893A4F",
  "Site_pos" = "#E79B80",
  "Site_neg" = "#8C4B3D",
  "nonsignificant" = "black"
)

ggplot(df_all_corrs, aes(x = K, y = pvalue)) +
  
  geom_jitter(data = df_all_corrs[is.na(df_all_corrs$label), ],
              aes(color = variable_colored, size = abs(correlation)),
              width = 0.25, alpha = 0.8) +
  
  geom_text_repel(data = df_all_corrs[!is.na(df_all_corrs$label), ],
                  aes(label = label, color = variable_colored),
                  size = 5, fontface = "bold",
                  direction = "y",
                  max.overlaps = Inf, min.segment.length = 0,
                  box.padding = 0.2, point.padding = 0.2, segment.color = NA) +
  
  scale_y_log10(limits = c(1e-6, 1),
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)) +
  
  scale_color_manual(values = pretty_colors, guide = "none") +
  scale_size_continuous(range = c(0.5, 2), guide = "none") +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
  
  labs(x = "Partition model K",
       y = "p-values",
       color = "Variable",
       title = "Partial Correlation between LEiDA lifetimes and covariates",
       subtitle = "Significance determined by FDR-corrected p-values < 0.05") +
  
  facet_wrap(~ variable_pretty, scales = "free_y") +
  
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.position = "bottom"
  )

# ==============================================================================
### Half-violin plot
## Race-ethnicity
# Filter only the desired variable
df_results_ethn <- df_all_corrs[df_all_corrs$variable == "race_ethnicity",]
df_results_ethn <- df_results_ethn[!is.na(df_results_ethn$label),]

# List to save all graphs
plot_list <- list()

for (i in 1:nrow(df_results_ethn)) {
  
  leida_state <- paste0("LT_k", df_results_ethn$K[i], "c", df_results_ethn$state[i])
  tmp_state <- y_data[[leida_state]]
  
  # Temporary dataset
  X_tmp <- data.frame(
    ethn = factor(X_data$race_ethnicity, labels = c("W", "B", "H", "A", "other")),
    lifetime = tmp_state
  )
  
  # Half violin plot
  p <- ggplot(X_tmp, aes(x = ethn, y = lifetime, fill = ethn)) +
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4) +
    # geom_jitter(width = 0.05, size = 0.5, alpha = 0.5) +
    scale_fill_manual(values = c("#D07B3E", "#D9A77C", "#D99262", "#DBB796", "#B56445" )) +
    coord_cartesian(ylim = c(0, 20)) +
    labs(title = paste("LEiDA:", leida_state), x = "Ethnicity", y = "Lifetime") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Combine all graphs in the same figure
# Specify the desired number of columns
n_col <- 3
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
print(combined_plot)

# ------------------------------------------------------------------------------
## Sex 
# Filter only the desired variable
df_results_sex <- df_all_corrs[df_all_corrs$variable == "demo_sex_v2",]
df_results_sex <- df_results_sex[!is.na(df_results_sex$label),]

# List to save all graphs
plot_list <- list()

for (i in 1:nrow(df_results_sex)) {
  
  leida_state <- paste0("LT_k", df_results_sex$K[i], "c", df_results_sex$state[i])
  tmp_state <- y_data[[leida_state]]
  
  # Temporary dataset
  X_tmp <- data.frame(
    sex = factor(X_data$demo_sex_v2, labels = c("M", "F")),
    lifetime = tmp_state
  )
  
  # Half violin plot
  p <- ggplot(X_tmp, aes(x = sex, y = lifetime, fill = sex)) +
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4) +
    # geom_jitter(width = 0.05, size = 0.5, alpha = 0.5) +
    scale_fill_manual(values = c( "#A87FC3","#B89ED9")) +
    coord_cartesian(ylim = c(0, 20)) +
    labs(title = paste("LEiDA:", leida_state), x = "TMI", y = "Lifetime") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Combine all graphs in the same figure
# Specify the desired number of columns
n_col <- 5
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
print(combined_plot)

# ==============================================================================
### Scatter plot
## Interview age
# Filter only the desired variable
df_results_age <- df_all_corrs[df_all_corrs$variable == "interview_age",]
df_results_age <- df_results_age[!is.na(df_results_age$label),]

# List for plots
plot_list <- list()

for (i in 1:nrow(df_results_age)) {
  
  leida_state <- paste0("LT_k", df_results_age$K[i], "c", df_results_age$state[i])
  tmp_state <- y_data[[leida_state]]
  
  # Temporary dataset
  X_tmp <- data.frame(
    age = as.numeric(X_data$interview_age), 
    lifetime = tmp_state
  )
  
  # Scatterplot
  p <- ggplot(X_tmp, aes(x = age, y = lifetime, color = age)) +
    geom_point(alpha = 0.6, size = 1.2) +
    # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    scale_color_gradient(low = "#6BA6CE", high = "#4A84B2") +
    coord_cartesian(ylim = c(0, 20)) +
    labs(title = paste("LEiDA:", leida_state),
         x = "Interview age",
         y = "Lifetime") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 14),  
          axis.title = element_text(size = 12),  
          axis.text  = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Final combination
n_col <- 3
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
suppressMessages(print(combined_plot))

# ------------------------------------------------------------------------------
## PDS
# Filter only the desired variable
df_results_pds <- df_all_corrs[df_all_corrs$variable == "pubertal_developmental_scale",]
df_results_pds <- df_results_pds[!is.na(df_results_pds$label),]

# List for plots
plot_list <- list()

for (i in 1:nrow(df_results_pds)) {
  
  leida_state <- paste0("LT_k", df_results_pds$K[i], "c", df_results_pds$state[i])
  tmp_state <- y_data[[leida_state]]
  
  # Temporary dataset
  X_tmp <- data.frame(
    pds = as.numeric(X_data$pubertal_developmental_scale), 
    lifetime = tmp_state
  )
  
  # Scatterplot
  p <- ggplot(X_tmp, aes(x = pds, y = lifetime, color = pds)) +
    geom_point(alpha = 0.6, size = 1.2) +
    # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    scale_color_gradient(low = "#B05C77", high = "#D37898") +
    coord_cartesian(ylim = c(0, 20)) +
    labs(title = paste("LEiDA:", leida_state),
         x = "Pubertal developmental scale",
         y = "Lifetime") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Final combination
n_col <- 1
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
suppressMessages(print(combined_plot))

# ------------------------------------------------------------------------------
## TMI
# Filter only the desired variable
df_results_tmi <- df_all_corrs[df_all_corrs$variable == "triponderal_mass_index",]
df_results_tmi <- df_results_tmi[!is.na(df_results_tmi$label),]

# List for plots
plot_list <- list()

for (i in 1:nrow(df_results_tmi)) {
  
  leida_state <- paste0("LT_k", df_results_tmi$K[i], "c", df_results_tmi$state[i])
  tmp_state <- y_data[[leida_state]]
  
  # Temporary dataset
  X_tmp <- data.frame(
    tmi = as.numeric(X_data$triponderal_mass_index), 
    lifetime = tmp_state
  )
  
  # Scatterplot
  p <- ggplot(X_tmp, aes(x = tmi, y = lifetime, color = tmi)) +
    geom_point(alpha = 0.6, size = 1.2) +
    # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    scale_color_gradient(low = "#7F9F78", high = "#5F805C") +
    coord_cartesian(ylim = c(0, 40)) +
    labs(title = paste("LEiDA:", leida_state),
         x = "Triponderal mass index",
         y = "Lifetime") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 14),  
          axis.title = element_text(size = 12),
          axis.text  = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Final combination
n_col <- 3
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
suppressMessages(print(combined_plot))

# ==============================================================================
### Heatmap
## TMI
df_results_all <- df_all_corrs[!is.na(df_all_corrs$label),]
df_results_var <- df_results_all[df_results_all$variable == "triponderal_mass_index",]

# Avoid log(0)
df_results_var <- df_results_var %>%
  mutate(log_p = -log10(pvalue),
         log_p = ifelse(is.infinite(log_p), NA, log_p),
         correlation_label = sprintf("%.2f", correlation))

# Calculate dynamic range for p-value logarithmic scale
log_p_min <- min(df_results_var$log_p, na.rm = TRUE)
log_p_max <- max(df_results_var$log_p, na.rm = TRUE)

# Heatmap
ggplot(df_results_var, aes(x = factor(K), y = factor(label), fill = log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = correlation_label), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#F5F5F5", high = "#8CB887",
    name = expression(-log[10](p)),
    limits = c(log_p_min, log_p_max),
    oob = scales::squish
  ) +
  labs(
    title = "Heatmap | Correlation between LEiDA lifetime and TMI",
    x = "K",
    y = "Connection (label)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

# ------------------------------------------------------------------------------
## Interview age
df_results_all <- df_all_corrs[!is.na(df_all_corrs$label),]
df_results_var <- df_results_all[df_results_all$variable == "interview_age",]

# Avoid log(0)
df_results_var <- df_results_var %>%
  mutate(log_p = -log10(pvalue),
         log_p = ifelse(is.infinite(log_p), NA, log_p),
         correlation_label = sprintf("%.2f", correlation))

# Calculate dynamic range for p-value logarithmic scale
log_p_min <- min(df_results_var$log_p, na.rm = TRUE)
log_p_max <- max(df_results_var$log_p, na.rm = TRUE)

# Heatmap
ggplot(df_results_var, aes(x = factor(K), y = factor(label), fill = log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = correlation_label), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#F5F5F5", high = "#6899C6",
    name = expression(-log[10](p)),
    limits = c(log_p_min, log_p_max),
    oob = scales::squish
  ) +
  labs(
    title = "Heatmap | Correlation between LEiDA lifetime and age",
    x = "K",
    y = "Connection (label)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

# ------------------------------------------------------------------------------
## PDS
df_results_all <- df_all_corrs[!is.na(df_all_corrs$label),]
df_results_var <- df_results_all[df_results_all$variable == "pubertal_developmental_scale",]

# Avoid log(0)
df_results_var <- df_results_var %>%
  mutate(log_p = -log10(pvalue),
         log_p = ifelse(is.infinite(log_p), NA, log_p),
         correlation_label = sprintf("%.2f", correlation))

# Calculate dynamic range for p-value logarithmic scale
log_p_min <- min(df_results_var$log_p, na.rm = TRUE)
log_p_max <- max(df_results_var$log_p, na.rm = TRUE)

# Heatmap
ggplot(df_results_var, aes(x = factor(K), y = factor(label), fill = log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = correlation_label), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#F5F5F5", high = "#D1839F",
    name = expression(-log[10](p)),
    limits = c(log_p_min, log_p_max),
    oob = scales::squish
  ) +
  labs(
    title = "Heatmap | Correlation between LEiDA lifetime and PDS",
    x = "K",
    y = "Connection (label)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

################################################################################
# LEiDA transitions
### Partial correlation analysis
## Divide in two datasets
X_data <- data %>%
  dplyr::select(demo_sex_v2, race_ethnicity, interview_age, 
                pubertal_developmental_scale, triponderal_mass_index)
names_covariates <- colnames(X_data)

## Extract LEiDA transitions
y_data <- df_leida %>%
  dplyr::select(matches("TR_K"))

full_results_trans <- list()

for (K in 2:20) {
  
  leida_k <- y_data %>%
    dplyr::select(matches(paste0("TR_K", K, "_"))) %>%
    dplyr::select(any_of(leida_cols)) # select only those present in leida_cols
  
  K_results <- list()
  
  for (s1 in 1:K) {
    for (s2 in 1:K) {
      
      colname <- paste0("TR_K", K, "_C", s1, "x", s2)
      
      # Skip if current state has been excluded
      if (!(colname %in% leida_cols)) next
      
      X_tmp <- X_data
      X_tmp$tmp_trans <- leida_k[[colname]]
      X_tmp <- na.omit(X_tmp)
      
      test <- pcor(X_tmp[, c("tmp_trans", names_covariates)], method = "pearson")
      
      cors <- test$estimate["tmp_trans", ]
      pvals <- test$p.value["tmp_trans", ]
      
      K_results[[colname]] <- data.frame(
        K = K,
        from = s1,
        to = s2,
        variable = names(cors),
        correlation = cors,
        pvalue = pvals,
        stringsAsFactors = FALSE
      )
    }
  }
  full_results_trans[[paste0("K", K)]] <- do.call(rbind, K_results)
}

df_all_corrs <- bind_rows(full_results_trans)
df_all_corrs <- df_all_corrs %>%
  group_by(K, variable) %>%
  mutate(pvalue_fdr = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()

# ==============================================================================
### Visualization
df_all_corrs <- df_all_corrs %>%
  filter(variable != "tmp_trans") %>%
  mutate(
    label = ifelse(pvalue_fdr < 0.05,
                   paste0(from, "→", to),
                   NA),
    K = factor(K)
  )

## Rename variables for plot visualization
pretty_names <- c(
  "demo_sex_v2" = "Sex",
  "interview_age" = "Age",
  "pubertal_developmental_scale" = "Pubertal Development",
  "race_ethnicity" = "Ethnicity",
  "triponderal_mass_index" = "TMI",
  "site" = "Site"
)

# Add column with readable names
df_all_corrs <- df_all_corrs %>%
  mutate(variable_pretty = recode(variable, !!!pretty_names))

## Define colors
df_all_corrs <- df_all_corrs %>%
  mutate(variable_colored = ifelse(pvalue_fdr < 0.05, variable_pretty, "nonsignificant"))

# Create categories variable+sign
df_all_corrs <- df_all_corrs %>%
  mutate(
    variable_colored = case_when(
      pvalue_fdr < 0.05 & correlation >= 0 ~ paste0(variable_pretty, "_pos"),
      pvalue_fdr < 0.05 & correlation < 0 ~ paste0(variable_pretty, "_neg"),
      TRUE ~ "nonsignificant"
    )
  )

pretty_colors <- c(
  "Sex_pos" = "#B88DC9",
  "Sex_neg" = "#6F3371",
  "TMI_pos" = "#8FB78F",
  "TMI_neg" = "#3F593F",
  "Ethnicity_pos" = "#D98C57",
  "Ethnicity_neg" = "#87421E",
  "Age_pos" = "#5F88B4",
  "Age_neg" = "#1F3C63",
  "Pubertal Development_pos" = "#C47088",
  "Pubertal Development_neg" = "#742940",
  "Site_pos" = "#C77A63",
  "Site_neg" = "#6A3A2E",
  "nonsignificant" = "black"
)

ggplot(df_all_corrs, aes(x = K, y = pvalue)) +
  
  geom_jitter(data = df_all_corrs[is.na(df_all_corrs$label), ],
              aes(color = variable_colored, size = abs(correlation)),
              width = 0.25, alpha = 0.8) +
  
  geom_text_repel(data = df_all_corrs[!is.na(df_all_corrs$label), ],
                  aes(label = label, color = variable_colored),
                  size = 2, fontface = "bold",
                  direction = "y",
                  max.overlaps = Inf, min.segment.length = 0,
                  box.padding = 0.2, point.padding = 0.2, segment.color = NA) +
  
  scale_y_log10(limits = c(1e-19, 1),
                breaks = 10^seq(-19, 0, by = 1)) +
  
  scale_color_manual(values = pretty_colors, guide = "none") +
  scale_size_continuous(range = c(0.5, 2), guide = "none") +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
  
  labs(x = "Partition model K",
       y = "p-values",
       color = "Variable",
       title = "Partial Correlation between LEiDA transitions and covariates",
       subtitle = "Significance determined by FDR-corrected p-values < 0.05") +
  
  facet_wrap(~ variable_pretty, scales = "free_y") +
  
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.position = "bottom"
  )

# ==============================================================================
### Half-violin plot
## Race-ethnicity
# Filter only the desired variable
df_results_ethn <- df_all_corrs[df_all_corrs$variable == "race_ethnicity",]
df_results_ethn <- df_results_ethn[!is.na(df_results_ethn$label),]

# List to save all graphs
plot_list <- list()

for (i in 1:nrow(df_results_ethn)) {
  
  leida_state <- paste0("TR_K", df_results_ethn$K[i], 
                        "_C", df_results_ethn$from[i],
                        "x", df_results_ethn$to[i])
  tmp_state_raw <- y_data[[leida_state]]
  
  # Remove outliers (keep only within mean ± 3*sd)
  mu <- mean(tmp_state_raw, na.rm = TRUE)
  sigma <- sd(tmp_state_raw, na.rm = TRUE)
  tmp_state <- ifelse(
    tmp_state_raw >= (mu - 3 * sigma) & tmp_state_raw <= (mu + 3 * sigma),
    tmp_state_raw,
    NA  # set outliers as NA
  )
  
  # Temporary dataset
  X_tmp <- data.frame(
    ethn = factor(X_data$race_ethnicity, labels = c("W", "B", "H", "A", "other")),
    transition = tmp_state
  )
  
  # Calculate y_min and y_max for each transition
  y_min <- min(tmp_state, na.rm = TRUE)
  y_max <- max(tmp_state, na.rm = TRUE)
  buffer <- (y_max - y_min) * 0.1
  coord_limits <- coord_cartesian(ylim = c(y_min - buffer, y_max + buffer))
  
  # Half violin plot
  p <- ggplot(X_tmp, aes(x = ethn, y = transition, fill = ethn)) +
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4) +
    # geom_jitter(width = 0.05, size = 0.5, alpha = 0.5) +
    scale_fill_manual(values = c("#9C4F2E", "#A67653", "#A4623D", "#A8876B", "#8A3F2C")) +
    coord_limits +
    labs(title = paste("LEiDA:", leida_state), x = "Ethnicity", y = "Transition") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 8))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Combine all graphs in the same figure
# Specify the desired number of columns
n_col <- 9
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
suppressWarnings(print(combined_plot))

# ------------------------------------------------------------------------------
## Sex
# Filter only the desired variable
df_results_sex <- df_all_corrs[df_all_corrs$variable == "demo_sex_v2",]
df_results_sex <- df_results_sex[!is.na(df_results_sex$label),]

# List to save all graphs
plot_list <- list()

for (i in 1:nrow(df_results_sex)) {
  
  leida_state <- paste0("TR_K", df_results_sex$K[i], 
                        "_C", df_results_sex$from[i],
                        "x", df_results_sex$to[i])
  tmp_state_raw <- y_data[[leida_state]]
  
  # Remove outliers (keep only within mean ± 3*sd)
  mu <- mean(tmp_state_raw, na.rm = TRUE)
  sigma <- sd(tmp_state_raw, na.rm = TRUE)
  tmp_state <- ifelse(
    tmp_state_raw >= (mu - 3 * sigma) & tmp_state_raw <= (mu + 3 * sigma),
    tmp_state_raw,
    NA  # set outliers as NA
  )
  
  # Temporary dataset
  X_tmp <- data.frame(
    sex = factor(X_data$demo_sex_v2, labels = c("M", "F")),
    transition = tmp_state
  )
  
  # Calculate y_min and y_max for each transition
  y_min <- min(tmp_state, na.rm = TRUE)
  y_max <- max(tmp_state, na.rm = TRUE)
  buffer <- (y_max - y_min) * 0.1
  coord_limits <- coord_cartesian(ylim = c(y_min - buffer, y_max + buffer))
  
  
  # Half violin plot
  p <- ggplot(X_tmp, aes(x = sex, y = transition, fill = sex)) +
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, 
                 point_colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4) +
    # geom_jitter(width = 0.05, size = 0.5, alpha = 0.5) +
    scale_fill_manual(values = c( "#8E66AD","#9E85C2")) +
    coord_limits +
    labs(title = paste(leida_state), x = "Sex", y = "Transition") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10),
          axis.text  = element_text(size = 7))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Combine all graphs in the same figure
# Specify the desired number of columns
n_col <- 11
combined_plot <- wrap_plots(plot_list, ncol = n_col) & scale_y_continuous()

# Visualize
suppressWarnings(print(combined_plot))

# ==============================================================================
### Scatter plot
## Interview age
# Filter only the desired variable
df_results_age <- df_all_corrs[df_all_corrs$variable == "interview_age",]
df_results_age <- df_results_age[!is.na(df_results_age$label),]

# List for plots
plot_list <- list()

for (i in 1:nrow(df_results_age)) {
  
  leida_state <- paste0("TR_K", df_results_age$K[i], 
                        "_C", df_results_age$from[i],
                        "x", df_results_age$to[i])
  tmp_state_raw <- y_data[[leida_state]]
  
  # Remove outliers (keep only within mean ± 3*sd)
  mu <- mean(tmp_state_raw, na.rm = TRUE)
  sigma <- sd(tmp_state_raw, na.rm = TRUE)
  tmp_state <- ifelse(
    tmp_state_raw >= (mu - 3 * sigma) & tmp_state_raw <= (mu + 3 * sigma),
    tmp_state_raw,
    NA  # set outliers as NA
  )
  
  # Temporary dataset
  X_tmp <- data.frame(
    age = as.numeric(X_data$interview_age), 
    transition = tmp_state
  )
  
  # Calculate y_min and y_max for each transition
  y_min <- min(tmp_state, na.rm = TRUE)
  y_max <- max(tmp_state, na.rm = TRUE)
  buffer <- (y_max - y_min) * 0.1
  coord_limits <- coord_cartesian(ylim = c(y_min - buffer, y_max + buffer))
  
  
  # Scatterplot
  p <- ggplot(X_tmp, aes(x = age, y = transition, color = age)) +
    geom_point(alpha = 0.6, size = 1.2) +
    # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    scale_color_gradient(low = "#3C6D9A", high = "#2A4F70") +
    coord_limits +
    labs(title = paste("LEiDA:", leida_state),
         x = "Interview age",
         y = "Transition") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Final combination
n_col <- 5
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
suppressWarnings(print(combined_plot))

# ------------------------------------------------------------------------------
## PDS
# Filter only the desired variable
df_results_pds <- df_all_corrs[df_all_corrs$variable == "pubertal_developmental_scale",]
df_results_pds <- df_results_pds[!is.na(df_results_pds$label),]

# List for plots
plot_list <- list()

for (i in 1:nrow(df_results_pds)) {
  
  leida_state <- paste0("TR_K", df_results_pds$K[i], 
                        "_C", df_results_pds$from[i],
                        "x", df_results_pds$to[i])
  tmp_state_raw <- y_data[[leida_state]]
  
  # Remove outliers (keep only within mean ± 3*sd)
  mu <- mean(tmp_state_raw, na.rm = TRUE)
  sigma <- sd(tmp_state_raw, na.rm = TRUE)
  tmp_state <- ifelse(
    tmp_state_raw >= (mu - 3 * sigma) & tmp_state_raw <= (mu + 3 * sigma),
    tmp_state_raw,
    NA  # set outliers as NA
  )
  
  # Temporary dataset
  X_tmp <- data.frame(
    pds = as.numeric(X_data$pubertal_developmental_scale), 
    transition = tmp_state
  )
  
  # Calculate y_min and y_max for each transition
  y_min <- min(tmp_state, na.rm = TRUE)
  y_max <- max(tmp_state, na.rm = TRUE)
  buffer <- (y_max - y_min) * 0.1
  coord_limits <- coord_cartesian(ylim = c(y_min - buffer, y_max + buffer))
  
  # Scatterplot
  p <- ggplot(X_tmp, aes(x = pds, y = transition, color = pds)) +
    geom_point(alpha = 0.6, size = 1.2) +
    # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    scale_color_gradient(low = "#7B3A50", high = "#9F4D6B") +
    coord_limits +
    labs(title = paste("LEiDA:", leida_state),
         x = "Pubertal developmental scale",
         y = "Transition") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Final combination
n_col <- 5
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
suppressWarnings(print(combined_plot))

# ------------------------------------------------------------------------------
## TMI
# Filter only the desired variable
df_results_tmi <- df_all_corrs[df_all_corrs$variable == "triponderal_mass_index",]
df_results_tmi <- df_results_tmi[!is.na(df_results_tmi$label),]

# List for plots
plot_list <- list()

for (i in 1:nrow(df_results_tmi)) {
  
  leida_state <- paste0("TR_K", df_results_tmi$K[i], 
                        "_C", df_results_tmi$from[i],
                        "x", df_results_tmi$to[i])
  tmp_state_raw <- y_data[[leida_state]]
  
  # Remove outliers (keep only within mean ± 3*sd)
  mu <- mean(tmp_state_raw, na.rm = TRUE)
  sigma <- sd(tmp_state_raw, na.rm = TRUE)
  tmp_state <- ifelse(
    tmp_state_raw >= (mu - 3 * sigma) & tmp_state_raw <= (mu + 3 * sigma),
    tmp_state_raw,
    NA  # set outliers as NA
  )
  
  # Temporary dataset
  X_tmp <- data.frame(
    tmi = as.numeric(X_data$triponderal_mass_index), 
    transition = tmp_state
  )
  
  # Calculate y_min and y_max for each transition
  y_min <- min(tmp_state, na.rm = TRUE)
  y_max <- max(tmp_state, na.rm = TRUE)
  buffer <- (y_max - y_min) * 0.1
  coord_limits <- coord_cartesian(ylim = c(y_min - buffer, y_max + buffer))
  
  # Scatterplot
  p <- ggplot(X_tmp, aes(x = tmi, y = transition, color = tmi)) +
    geom_point(alpha = 0.6, size = 1.2) +
    # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    scale_color_gradient(low = "#4E664B", high = "#364533") +
    coord_limits +
    labs(title = paste("LEiDA:", leida_state),
         x = "Triponderal mass index",
         y = "Transition") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))
  
  # Add to the list
  plot_list[[i]] <- p
}

## Final combination
n_col <- 7
combined_plot <- wrap_plots(plot_list, ncol = n_col)

# Visualize
suppressWarnings(print(combined_plot))

# ==============================================================================
### Heatmap
## TMI
df_results_all <- df_all_corrs[!is.na(df_all_corrs$label),]
df_results_var <- df_results_all[df_results_all$variable == "triponderal_mass_index",]

# Avoid log(0)
df_results_var <- df_results_var %>%
  mutate(log_p = -log10(pvalue),
         log_p = ifelse(is.infinite(log_p), NA, log_p),
         correlation_label = sprintf("%.2f", correlation))

# Calculate dynamic range for p-value logarithmic scale
log_p_min <- min(df_results_var$log_p, na.rm = TRUE)
log_p_max <- max(df_results_var$log_p, na.rm = TRUE)

# Heatmap
ggplot(df_results_var, aes(x = factor(K), y = factor(label), fill = log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = correlation_label), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#F5F5F5", high = "#4E664B",
    name = expression(-log[10](p)),
    limits = c(log_p_min, log_p_max),
    oob = scales::squish
  ) +
  labs(
    title = "Heatmap | Correlation between LEiDA transitions and TMI",
    x = "K",
    y = "Connection (label)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

# ------------------------------------------------------------------------------
## Interview age
df_results_all <- df_all_corrs[!is.na(df_all_corrs$label),]
df_results_var <- df_results_all[df_results_all$variable == "interview_age",]

# Avoid log(0)
df_results_var <- df_results_var %>%
  mutate(log_p = -log10(pvalue),
         log_p = ifelse(is.infinite(log_p), NA, log_p),
         correlation_label = sprintf("%.2f", correlation))

# Calculate dynamic range for p-value logarithmic scale
log_p_min <- min(df_results_var$log_p, na.rm = TRUE)
log_p_max <- max(df_results_var$log_p, na.rm = TRUE)

# Heatmap
ggplot(df_results_var, aes(x = factor(K), y = factor(label), fill = log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = correlation_label), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#F5F5F5", high = "#4E7CB3",
    name = expression(-log[10](p)),
    limits = c(log_p_min, log_p_max),
    oob = scales::squish
  ) +
  labs(
    title = "Heatmap | Correlation between LEiDA transitions and age",
    x = "K",
    y = "Connection (label)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

# ------------------------------------------------------------------------------
## PDS
df_results_all <- df_all_corrs[!is.na(df_all_corrs$label),]
df_results_var <- df_results_all[df_results_all$variable == "pubertal_developmental_scale",]

# Avoid log(0)
df_results_var <- df_results_var %>%
  mutate(log_p = -log10(pvalue),
         log_p = ifelse(is.infinite(log_p), NA, log_p),
         correlation_label = sprintf("%.2f", correlation))

# Calculate dynamic range for p-value logarithmic scale
log_p_min <- min(df_results_var$log_p, na.rm = TRUE)
log_p_max <- max(df_results_var$log_p, na.rm = TRUE)

# Heatmap
ggplot(df_results_var, aes(x = factor(K), y = factor(label), fill = log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = correlation_label), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#F5F5F5", high = "#7B3A50",
    name = expression(-log[10](p)),
    limits = c(log_p_min, log_p_max),
    oob = scales::squish
  ) +
  labs(
    title = "Heatmap | Correlation between LEiDA transitions and PDS",
    x = "K",
    y = "Connection (label)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))
