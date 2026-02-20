#!/bin/bash

# Code that runs the first three steps to prepare data for connectome analysis:
# 1- for each subject, convert file from .tsv to .txt and store in a separate folder xcpd_txt
# 2- for each subject, retain only the selected 114 ROIs among the given 156, and store in group_subcort_for_leida
# 3- for each subject, transpose the matrix obtained in previous step, and store in ready_for_leida

# Isabella L.C. Mariani Wigley, 03 / 2025; ilmawi@utu.fi
# Aurora Berto, 03 / 2025; aurber@utu.fi

##### Select the site
site="G031" # USER.adapt!

##### 1st step: get .tsv files in .txt format for each subject
# Base directory where subject folders are located
first_dir="/path/to/data/${site}/xcpd_noncifti" # USER.adapt!

# Output directory for converted .txt files
second_dir="/path/to/data/${site}/temp/xcpd_txt" # USER.adapt!
mkdir -p "$second_dir"

# Loop over each subject folder in the base directory
for subject_dir in "$first_dir"/sub-*; do
    # Check if the subject directory exists
    if [[ -d "$subject_dir" ]]; then
        func_dir="$subject_dir/ses-baselineYear1Arm1/func" # USER.adapt!

        # Check if the func subdirectory exists
        if [[ -d "$func_dir" ]]; then
            # Find files ending with 4S156Parcels_stat-mean_timeseries.tsv and convert to .txt
            for file in "$func_dir"/*4S156Parcels_stat-mean_timeseries.tsv; do # USER.adapt!
                # Check if file exists
                if [[ -f "$file" ]]; then
                    # Extract base filename and replace .tsv with .txt
                    base_filename=$(basename "$file" .tsv)
                    output_file="$second_dir/$base_filename.txt"

                    # Copy and rename the file to .txt format
                    cp "$file" "$output_file"
                    echo "Converted $file to $output_file"
                else
                    echo "No matching files found in $func_dir for subject $(basename "$subject_dir")"
                fi
            done
        else
            echo "Directory $func_dir does not exist for subject $(basename "$subject_dir")"
        fi
    else
        echo "Skipping non-directory item $subject_dir"
    fi
done
echo "All matching files have been converted to .txt and saved to $second_dir"


######################################################################################
##### 2nd step: obtain 114 ROIs for ENIGMA from the 156 selected ones for parcellation 
# Import python
module load python-data

# Directories for input and output
# Input directory is equal to the previous output one (second_dir)
third_dir="/path/to/data/${site}/temp/group_subcort_for_leida" # USER.adapt!
mkdir -p "$third_dir"

# Loop through each .txt file and process all files 
for file in "$second_dir"/*.txt; do
    # Define the output file path
    output_file="$third_dir/$(basename "$file" .txt)_grouped.txt"

    # Python code to calculate group means and keep non-grouped columns
    python3 - <<END
import pandas as pd

# File paths
input_file = "$file"
output_file = "$output_file"

# Define brain region groups with left and right hemisphere labels
BRAIN_REGION_GROUPS = {
    "Left_Accumbens": ["LH-NAC"],
    "Left_Amygdala": ["LH_Amygdala"],
    "Left_Caudate": ["LH-Ca"],
    "Left_Hippocampus": ["LH_Hippocampus"],
    "Left_Pallidum": ["LH-GPe", "LH-GPi"],
    "Left_Putamen": ["LH-Pu"],
    "Left_Thalamus": [
        "LH-Pulvinar", "LH-Anterior", "LH-Medio_Dorsal", "LH-Ventral_Latero_Dorsal",
        "LH-Central_Lateral-Lateral_Posterior-Medial_Pulvinar", "LH-Ventral_Anterior", "LH-Ventral_Latero_Ventral"
    ],
    "Right_Accumbens": ["RH-NAC"],
    "Right_Amygdala": ["RH_Amygdala"],
    "Right_Caudate": ["RH-Ca"],
    "Right_Hippocampus": ["RH_Hippocampus"],
    "Right_Pallidum": ["RH-GPe", "RH-GPi"],
    "Right_Putamen": ["RH-Pu"],
    "Right_Thalamus": [
        "RH-Pulvinar", "RH-Anterior", "RH-Medio_Dorsal", "RH-Ventral_Latero_Dorsal",
        "RH-Central_Lateral-Lateral_Posterior-Medial_Pulvinar", "RH-Ventral_Anterior", "RH-Ventral_Latero_Ventral"
    ]
}

# Load the original data
df = pd.read_csv(input_file, sep='\t')

# Cut cortical regions from original data
df_cortical = df.iloc[:,:100]
# print(f"Number of selected columns: {df_cortical.shape[1]}")

# Initialize a dictionary to store the grouped mean columns
group_means = {}

# Calculate mean time series for each group and store in group_means
for group, columns in BRAIN_REGION_GROUPS.items():
    # Ensure only columns present in the dataframe are used
    valid_columns = [col for col in columns if col in df.columns]
    if valid_columns:
        group_means[group] = df[valid_columns].mean(axis=1)
    else:
        print(f"Warning: Missing columns for group {group}")

# Convert group means to DataFrame
df_mean_subcortical = pd.DataFrame(group_means)
# print(f"Number of subcortical columns generated: {df_mean_subcortical.shape[1]}")

# Check that subcortical means are actually 14 (7 right + 7 left)
if df_mean_subcortical.shape[1]!=14:
    print(f"Expected 14 subcortical groups, but got {df_mean_subcortical.shape[1]}")

# Combine non-grouped columns with the grouped means
output_df = pd.concat([df_cortical, df_mean_subcortical], axis=1)

# Check that ROIs are actually 114
if output_df.shape[1]!=114:
    print(f"Expected 114 ROIs, but got {output_df.shape[1]}")

# Save the combined DataFrame (non-grouped columns + grouped means) to the output file
output_df.to_csv(output_file, sep='\t', index=False)
print(f"Processed and saved: {output_file}")

END
done
echo "Processing complete. Processed all files, saved in '$third_dir'."


######################################################################################
##### 3rd step: transpose matrices and delete column labels
# Define the input and output directories
# Input directory is equal to the previous output one (third_dir)
fourth_dir="/path/to/data/${site}/temp/ready_for_leida" # USER.adapt!

# Create the output directory if it does not exist
mkdir -p "$fourth_dir"

# Loop through all .txt files in the input directory
for input_file in "$third_dir"/*.txt; do
    echo "Processing $input_file..."
    
    # Define the output file name
    base_name=$(basename "$input_file" .txt)
    output_file="$fourth_dir/${base_name}_transposed.txt" 
    
    # Remove the header row and then transpose the matrix
    tail -n +2 "$input_file" | awk '
        {
            for (i=1; i<=NF; i++) {
                a[NR,i] = $i
            }
        }
        NF>max { max = NF }
        END {
            for (i=1; i<=max; i++) {
                for (j=1; j<=NR; j++) {
                    printf "%s%s", a[j,i], (j==NR ? "\n" : "\t")
                }
            }
        }
    ' > "$output_file"
    
    echo "Transpose complete for $input_file. Output saved to $output_file."
done
