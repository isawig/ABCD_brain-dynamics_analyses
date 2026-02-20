## Circos plot to visualize correlations between demographics and first 20 harmonics 

# Isabella L.C. Mariani Wigley, 06 / 2025; ilmawi@utu.fi
# Aurora Berto, 06 / 2025; aurber@utu.fi

###################################################################################

#0 load packages
library("tidyverse")
library("circlize")
library("dplyr")


#1 Prepare the data for Circos plot ---------
#Your dataset should be a matrix, with, for example, correlation coefficients in all cells.

#####
## Load the dataset with all the variables
# To change accordingly to your path
df_path <- "/Users/auroraberto/Desktop/24-25/MSc_thesis/main_article/createDataset"
df_file <- file.path(df_path,"main_dataset_rsfMRI_FCH_LEiDA_demo.csv")
df <- read.csv(df_file, header = T, stringsAsFactors = T)

## Estract subjects ID
library(dplyr)
sites_ID <- df %>%
  dplyr::select(matches("site|src_subject_id"))
sites_names <- unique(levels((sites_ID$site)))

## Create 21 variables with 0/1 levels to indicate subjects belonging to each site
for (site in sites_names) {
  var_name <- make.names(site)  
  col_name <- paste0("site_", var_name)
  sites_ID[[col_name]] <- ifelse(sites_ID$site == site, 1, 0)
}


## Interview age, sex and pubertal developmental score
df_demo <- df %>%
  dplyr::select(demo_sex_v2, race_ethnicity, interview_age, 
                pubertal_developmental_scale,
                rel_family_id, rel_group_id,
                triponderal_mass_index)

# 777 = refuse to answer, 999 = don't know set to NA
df_demo[] <- lapply(df_demo, function(x) {
  if (is.numeric(x)) x[x %in% c(777, 999)] <- NA
  x
})

#####
# POWER
df_power <- df %>%
  dplyr::select(matches("power"))
# Retain the first 20 harmonics 
df_power <- df_power[, c(1:20)]

# ENERGY
df_energy <- df %>%
  dplyr::select(matches("energy"))
# Retain the first 20 harmonics
df_energy <- df_energy[, c(1:20)]


#####
## Merge all datasets in a single one
df_analysis <- cbind(sites_ID, df_demo, df_power, df_energy)

# Drop duplicates based on rel_family_id and rel_group_id, keeping the first occurrence
df_analysis <- df_analysis %>%
  distinct(rel_family_id, rel_group_id, .keep_all = TRUE)

# Remove columns from the dataset
df_analysis <- df_analysis %>%
  dplyr::select(-rel_family_id, -rel_group_id)

# Delete empty columns
df_analysis <- df_analysis[, colSums(!is.na(df_analysis)) > 0]

# Define variables to keep
vars_to_keep <- colnames(df_analysis)
print(vars_to_keep)

# Remove subjects' ID
df_sub <- df_analysis %>%
  dplyr::select(-src_subject_id)

# Delete rows with NA
df_sub <- na.omit(df_sub) # 4453 x 46

## Remove unnecessary datasets
remove(df, df_demo, df_energy, df_power, sites_ID)

##### 
# CIRCOS PLOT 
# Circos plot considering harmonic variables, covariates and site

df_circos <- df_sub %>%
  dplyr::select(-matches(sites_names))
df_circos$site <- as.numeric(df_circos$site)

#####
# Calculate correlation matrix
library(corpcor)
cor_mat <- cor(df_circos, use = "pairwise.complete.obs", method = "pearson")
pcor_matrix <- cor2pcor(cor_mat)
colnames(pcor_matrix) <- rownames(pcor_matrix) <- colnames(df_circos)

# library(ppcor)
# cor_matrix <- ppcor::pcor(df_circos, method = "pearson")

#1.1. Replace all intra-dataset associations with zeros

#There would be other ways to do this.

# cor_matrix_projection <- as.data.frame(cor_matrix$estimate) 
cor_matrix_projection <- as.data.frame(pcor_matrix)
cor_matrix_projection[1:1, 1:1] <- 0 
cor_matrix_projection[2:6, 2:6] <- 0 
cor_matrix_projection[7:46, 7:46] <- 0     

#1.2 Replace all correlations deemed as spurious (based on NRR threshold) with zeros. In this example,
#there is a separate dataframe called NRR, which has TRUE or FALSE in every cell based on whether it passed the threshold.
#Basically, my approach is to pivot both datasets longer, cbind them, create a new variable, remove the old variables, and 
#pivot back wider. 

# Set the threshold and create a boolean matrix based on that
threshold <- 0.03
NRR <- abs(cor_matrix_projection) > threshold
NRR <- as.data.frame(NRR)


#There would be other ways to do this.

cor_matrix_projection_long <- cor_matrix_projection %>% 
  rownames_to_column() %>% 
  pivot_longer(2:47)

NRR_long <- rownames_to_column(NRR) %>% 
  pivot_longer(2:47)

cor_matrix_projection_NRR <- cbind(cor_matrix_projection_long, NRR_long) 

colnames(cor_matrix_projection_NRR) <- c("rowname", "name", "value", "NRR_rowname", "NRR_name", "NRR_value")

cor_matrix_projection_NRR <- cor_matrix_projection_NRR %>% 
  dplyr::mutate(new_data = if_else(NRR_value == F, 0, value)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(rowname, name, new_data) %>% 
  pivot_wider(values_from = "new_data", names_from = "name") %>% 
  column_to_rownames()

#1.3 Lastly, remove the duplicates caused by every variable being both a row and a column, and therefore every
#association being shown twice. The method I've used may not seem like the absolute best, so might require some careful quality
#control, to check and ensure that it doesn't get rid of any real associations that you want to plot.
# circos_matrix_NRR <- cor_matrix_projection_NRR %>%
#   rownames_to_column() %>%
#   pivot_longer(2:26) %>%
#   filter(value == 0 | !duplicated(value)) %>%
#   pivot_wider(values_from = "value", names_from = "name") %>%
#   column_to_rownames() %>%  as.matrix()

circos_matrix_NRR <- cor_matrix_projection_NRR %>%
  rownames_to_column(var = "rowname") %>%    
  pivot_longer(-rowname, names_to = "colname", values_to = "value") %>% 
  # Keep only the upper triangular matrix
  filter(match(colname, rowname) >= match(rowname, rowname)) %>% 
  pivot_wider(names_from = colname, values_from = value) %>%  
  column_to_rownames(var = "rowname") %>% 
  as.matrix()



#2. Prepare the formatting for Circos plot -----------------------------

#2.1 Create a colour list vector for all, so each type of variable is coloured the same.
# Unfortunately I don't know any way to do this less manually or less laboriously. 
grid_col <- c(site= "#EB984E", 
              
              triponderal_mass_index= "#58D68D", demo_sex_v2= "#58D68D", race_ethnicity= "#58D68D",
              interview_age= "#58D68D", pubertal_developmental_scale= "#58D68D", #birth_weight_g= "#58D68D",
              
              Harmonics_power1= "lightblue", Harmonics_power2= "lightblue", Harmonics_power3= "lightblue", 
              Harmonics_power4= "lightblue", Harmonics_power5= "lightblue", Harmonics_power6= "lightblue", 
              Harmonics_power7= "lightblue", Harmonics_power8= "lightblue", Harmonics_power9= "lightblue", 
              Harmonics_power10= "lightblue", Harmonics_power11= "lightblue", Harmonics_power12= "lightblue", 
              Harmonics_power13= "lightblue", Harmonics_power14= "lightblue",  Harmonics_power15= "lightblue", 
              Harmonics_power16= "lightblue",  Harmonics_power17= "lightblue",  Harmonics_power18= "lightblue", 
              Harmonics_power19= "lightblue",  Harmonics_power20= "lightblue",
              
              Harmonics_energy1= "lightgrey", Harmonics_energy2= "lightgrey", Harmonics_energy3= "lightgrey", 
              Harmonics_energy4= "lightgrey", Harmonics_energy5= "lightgrey", Harmonics_energy6= "lightgrey",
              Harmonics_energy7= "lightgrey", Harmonics_energy8= "lightgrey", Harmonics_energy9= "lightgrey", 
              Harmonics_energy10= "lightgrey", Harmonics_energy11= "lightgrey", Harmonics_energy12= "lightgrey",
              Harmonics_energy13= "lightgrey", Harmonics_energy14= "lightgrey",  Harmonics_energy15= "lightgrey", 
              Harmonics_energy16= "lightgrey",  Harmonics_energy17= "lightgrey",  Harmonics_energy18= "lightgrey",
              Harmonics_energy19= "lightgrey",  Harmonics_energy20= "lightgrey"
              
)

#2.2 Create a legend using the ComplexHeatmap package, which we will then add manually to the circos plot.
library("ComplexHeatmap")

variable_types <- as.data.frame(x = c("Site", "Demographics", "FCH power", "FCH energy"))

variable_colours <- as.data.frame(x = c("#EB984E", "#58D68D","lightblue", "lightgrey"))

variable_types <- cbind(variable_types, variable_colours)

colnames(variable_types) <- c("type", "colour")

legend_variables = Legend(at = c(variable_types$type), type = "points", 
                          legend_gp = gpar(col = variable_types$colour), title_position = "topleft", 
                          title = "Variable Type", size = unit(12.5, "points"), grid_height = unit(12.5, "points"),
                          grid_width = unit(12.5, "points"),
                          labels_gp = gpar(fontsize = 12.5),
                          title_gp = gpar(fontsize = 12.5, fontface = "bold"), nrow = 1)

legend_links = Legend(at = c("Positive", "Negative"), type = "lines", 
                      legend_gp = gpar(col = c("#2980B9", "#EB984E"), lwd = 7.5), title_position = "topleft", 
                      title = "Association", grid_height = unit(12.5, "points"),
                      grid_width = unit(12.5, "points"),
                      labels_gp = gpar(fontsize = 12.5),
                      title_gp = gpar(fontsize = 12.5, fontface = "bold"), nrow = 1)

#Create a matrix for the colours of the links - blue for positive values and orange for negative values.
link_colours_NRR <- cor_matrix_projection_NRR %>% 
  rownames_to_column() %>% 
  pivot_longer(2:47) %>% 
  filter(value == 0 | !duplicated(value)) %>% 
  mutate(colour = case_when(value > 0 ~ "#2980B9", value < 0 ~ "#EB984E",)) %>% 
  mutate(colour = adjustcolor(colour, alpha = 0.85)) %>% 
  dplyr::select(-value) %>% 
  pivot_wider(values_from = "colour", names_from = "name") %>% 
  column_to_rownames() %>%  as.matrix()

legend_full = packLegend(legend_links, legend_variables, direction = "horizontal")

#3 Construct and save the circos plot -----------------------------------------------------------------

#Run this line each time I want to edit the diagram, or create a new one, as opposed to build upon it
circos.clear()

#Rotate the plot as you wish. This needs to be run after circos.clear, but before making and calling the plot. 
#So potentially lots of back and forth. 210 seems arbitrary, it was just the most suitable to arrange my own plot
#based on the label names.
circos.par(start.degree = 210, 
           gap.degree = 2,        # riduce lo spazio tra i settori
           track.margin = c(0.01, 0.01),  
           track.height = 0.3)   # altezza della traccia più piccola = cerchio più piccolo

#Create the chord diagram
chordDiagram(circos_matrix_NRR, #input data
             grid.col = grid_col, #colours of the grid
             annotationTrack = "grid", #what to display
             preAllocateTracks = 1, #This creates an empty "track", which we will later populate with labels/names
             col = link_colours_NRR, #colours of the links
             link.border = "white", #link border colour
             order = colnames(circos_matrix_NRR) #order of the edges
             # link.lwd = 0.1
)


#Create sector labels
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  new_label <- switch(sector.name,
                      "triponderal_mass_index" = "TMI",
                      "demo_sex_v2" = "sex",
                      "race_ethnicity" = "ethnicity",
                      "interview_age" = "interview age",
                      "pubertal_developmental_scale" = "PDS",
                      # "birth_weight_g" = "birth weight (g)",
                      
                      "Harmonics_power1" = "H1 power",
                      "Harmonics_power2" = "H2 power",
                      "Harmonics_power3" = "H3 power",
                      "Harmonics_power4" = "H4 power",
                      "Harmonics_power5" = "H5 power",
                      "Harmonics_power6" = "H6 power",
                      "Harmonics_power7" = "H7 power",
                      "Harmonics_power8" = "H8 power",
                      "Harmonics_power9" = "H9 power",
                      "Harmonics_power10" = "H10 power",
                      "Harmonics_power11" = "H11 power",
                      "Harmonics_power12" = "H12 power",
                      "Harmonics_power13" = "H13 power",
                      "Harmonics_power14" = "H14 power",
                      "Harmonics_power15" = "H15 power",
                      "Harmonics_power16" = "H16 power",
                      "Harmonics_power17" = "H17 power",
                      "Harmonics_power18" = "H18 power",
                      "Harmonics_power19" = "H19 power",
                      "Harmonics_power20" = "H20 power",
                      
                      "Harmonics_energy1" = "H1 energy",
                      "Harmonics_energy2" = "H2 energy",
                      "Harmonics_energy3" = "H3 energy",
                      "Harmonics_energy4" = "H4 energy",
                      "Harmonics_energy5" = "H5 energy",
                      "Harmonics_energy6" = "H6 energy",
                      "Harmonics_energy7" = "H7 energy",
                      "Harmonics_energy8" = "H8 energy",
                      "Harmonics_energy9" = "H9 energy",
                      "Harmonics_energy10" = "H10 energy",
                      "Harmonics_energy11" = "H11 energy",
                      "Harmonics_energy12" = "H12 energy",
                      "Harmonics_energy13" = "H13 energy",
                      "Harmonics_energy14" = "H14 energy",
                      "Harmonics_energy15" = "H15 energy",
                      "Harmonics_energy16" = "H16 energy",
                      "Harmonics_energy17" = "H17 energy",
                      "Harmonics_energy18" = "H18 energy",
                      "Harmonics_energy19" = "H19 energy",
                      "Harmonics_energy20" = "H20 energy",
                      sector.name) # default
  
  #print labels 
  circos.text(mean(xlim), ylim[1] + 2, new_label, 
              facing = "clockwise", 
              niceFacing = TRUE, #adjust for human eyes (T/F)
              adj = c(0.1, 0.1), #label rotation
              cex=1.3) #fontsize
  
  
}, bg.border = NA)

# Add the legend to the plot
draw(legend_full, x = unit(0.02, "npc"), y = unit(0.02, "npc"), just = c("left", "bottom"))

#Save
# dev.copy(png,"plots/chord_diagram_NRR.png", units = "px", width=2500, height=2500, bg = "transparent", res = 200)
# dev.off()
circos.clear()




#####
# CORRELATION MATRIX WITH SITES
df_corr <- df_analysis %>%
  dplyr::select(-src_subject_id,-site, -triponderal_mass_index, -demo_sex_v2, -race_ethnicity, -interview_age, -pubertal_developmental_scale)

library(ppcor)
cor_matrix <- pcor(df_corr, method = "pearson")$estimate
rownames(cor_matrix) <- colnames(df_corr)
colnames(cor_matrix) <- colnames(df_corr)

# 2. Estrai i nomi delle variabili
vars <- colnames(df_corr)
site_vars <- vars[1:21]
harmonics_vars <- vars[22:length(vars)]

# 3. Estrai solo le correlazioni tra siti e harmonics
cor_submatrix <- cor_matrix[site_vars, harmonics_vars]

# 4. Visualizza con ggcorrplot (modificato per matrice rettangolare)
library(ggcorrplot)

ggcorrplot(t(cor_submatrix),
           method = "square",
           lab = TRUE,
           lab_size = 2.5,
           colors = c("#EB984E", "white", "#2980B9"),
           outline.col = "gray",
           tl.cex = 8,
           tl.srt = 45)


#####
# HALF-VIOLIN PLOT - ETHNICITY
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(patchwork)

# Filter only the desired variable
df_results <- df_circos %>%
  dplyr::select(race_ethnicity,
                Harmonics_power14,  Harmonics_power18,
                Harmonics_energy14, Harmonics_energy18)

df_results_long <- df_results %>%
  mutate(
    race_ethnicity = factor(
      race_ethnicity,
      labels = c("W", "B", "H", "A", "other")
    )
  ) %>%
  pivot_longer(
    cols = c(
      Harmonics_power14,
      Harmonics_power18,
      Harmonics_energy14,
      Harmonics_energy18
    ),
    names_to = "harmonic_var",
    values_to = "value"
  )


plot_list <- list()

for (v in vars) {
  
  X_tmp <- df_results_long %>%
    filter(harmonic_var == v)
  
  # scegli ylim in base alla variabile
  ylim_use <- if (grepl("power", v)) c(0,1.8) else c(0,1000)
  
  p <- ggplot(X_tmp, aes(x = race_ethnicity, y = value, fill = race_ethnicity)) +
    
    stat_halfeye(
      adjust = 0.5,
      width = 0.6,
      .width = 0,
      justification = -0.3,
      point_colour = NA
    ) +
    
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      alpha = 0.4
    ) +
    
    scale_fill_manual(
      values = c("#F3A974", "#FACDA8", "#F7BA8D", "#FDDDBC", "#E18A5E")
    ) +
    
    coord_cartesian(ylim = ylim_use) +
    
    labs(
      title = v,
      x = "Ethnicity",
      y = "Harmonics value"
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold")
    )
  
  plot_list[[v]] <- p
}


combined_plot <- wrap_plots(plot_list, ncol = 2)
print(combined_plot)

#####
# HALF-VIOLIN PLOT - SEX
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(patchwork)

# Filter only the desired variable
df_results_sex <- df_circos %>%
  dplyr::select(demo_sex_v2, Harmonics_power7) %>%
  mutate(
    sex = factor(demo_sex_v2,
                 levels = c(1, 2),
                 labels = c("M", "F"))
  )

p_sex <- ggplot(
  df_results_sex,
  aes(x = sex, y = Harmonics_power7, fill = sex)
) +
  
  stat_halfeye(
    adjust = 0.5,
    width = 0.6,
    .width = 0,
    justification = -0.3,
    point_colour = NA
  ) +
  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.4
  ) +
  
  scale_fill_manual(
    values = c("M" = "#D7B8EB", "F" = "#E3C9F2")
  ) +
  
  labs(
    title = "Harmonics_power7",
    x = "Sex",
    y = "Harmonics value"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold")
  )
print(p_sex)
