#!/bin/bash
# Script 1 - Quality control
# Isabella L.C. Mariani Wigley, 02 / 2025; ilmawi@utu.fi
# Aurora Berto, 02 / 2025; aurber@utu.fi

# Code to list all data folders of subjects that have passed the quality check (framewise displacement, FD < 0.4).
# It takes in input the .txt/.csv file with the list of all subjects that passed QC, and compares with subjects folder with available data.
# The output is the list of all subjects that completed QC and whose data are available for the download.

# Directories
# Watch out: change the current directories according to your data storage path.
base_dir="/path/to/preprocessed/xcpd_noncifti"  # USER.adapt!
output_dir="path/to/destination/folder"         # USER.adapt!

# Create folders to store results if not exist
mkdir -p "$output_dir/complete_fd_below_0.4"
mkdir -p "$output_dir/complete_fd_above_0.4"
mkdir -p "$output_dir/subjs_missing_from_txt"

# Allocate variables
sites=("G010" "G031" "G032" "G075" "G087" "S012" "S013" "S014" "S020" "S021" "S022" "S042" "S053" "S065" "S076" "S086" "S090" "S011" "P023" "P043" "P064")  # USER.adapt!
subj_file="_non_faulty_xcpd_noncifti.txt" # ending of files' names, USER.adapt!

# Load file with subjects to keep: convert in csv format to make analysis easier
if [[ ! -f "1_fd_below_0.4.csv" ]]; then
    if [[ -f "1_fd_below_0.4.xlsx" ]]; then
        echo "The file .csv doesn't exist: converting from .xlsx"
        xlsx2csv 1_fd_below_0.4.xlsx 1_fd_below_0.4.csv
    else
        echo "Error: the required file doesn't exist."
    fi
fi
file_csv="1_fd_below_0.4.csv"

# Operation to do for each site
for site in "${sites[@]}"; do
    echo " "
    echo "Analyzing site: $site"

    # Go to the actual site's raw file
    file_tmp="$base_dir/${site}$subj_file"
    
    # Check if the file exists
    if [[ ! -f "$file_tmp" ]]; then
        echo "The file doesn't exist for site: ${site}$subj_file"
        continue # go to the next site
    fi

    # Length of raw file (number of subjects to analyze)
    lines=$(wc -l < "$file_tmp")

    # Extract subjects from the CSV file from CSV file (keep only 2nd column)
    # Add also a check to delete empty lines if present
    xlsx_tmp=$(awk -F',' -v value="${site}" '$1 == value {print $2}' "$file_csv" | grep -v '^$'| tr '\n' ' ')
    
    # Convert CSV file into array
    IFS=' ' read -r -a subjs_to_keep <<< "$xlsx_tmp"  

    # Number of subjects in CSV file for the given site
    lines_toKeep=${#subjs_to_keep[@]}

    # Create empty output files (to avoid overwritting on the same file with multiple runs)
    output_file_below="${output_dir}/complete_fd_below_0.4/${site}_fd_below.txt"
    output_file_above="${output_dir}/complete_fd_above_0.4/${site}_fd_above.txt"
    output_file_missing="${output_dir}/subjs_missing_from_txt/${site}_missing.txt"
    > "$output_file_below"
    > "$output_file_above"
    > "$output_file_missing"

    # Initialize variables to track classified subject 
    below_tmp=0 # number of subjects who passed the QC
    above_tmp=0 # number of subjects who failed the QC

    # Initialize variables to store classified subjects
    subj_in_file_tmp=() # subjects considered in each site
    missing_subjects=() # array of missing subjects in raw file


    # For each subject in site's raw file
    while IFS= read -r subj; do

        # Store them in a separate array (for further check)
        subj_in_file_tmp+=("$subj") 

        # Check if subject has passed or failed QC and store it in the correspondent array
        if [[ " ${subjs_to_keep[@]} " =~ " $subj " ]]; then
            echo "$subj" >> "$output_file_below"
            ((below_tmp++))
        else
            echo "$subj" >> "$output_file_above"
            ((above_tmp++))
        fi
    done < "$file_tmp"

    # Look for subjects present only in CSV file but not in raw .txt file 
    for subj in "${subjs_to_keep[@]}"; do
        if [[ ! " ${subj_in_file_tmp[@]} " =~ " $subj " ]]; then
            missing_subjects+=("$subj")
            echo "$subj" >> "$output_file_missing"
        fi
    done

    # Show quality check summary
    echo "--- Quality control summary ---"
    echo "Number of subjects evaluated: $lines"
    echo "Number of subjects to compare with: $lines_toKeep"
    echo "Number of subjects who passed QC: $below_tmp"
    echo "Number of subjects who failed QC: $above_tmp"
    if [[ ${#missing_subjects[@]} -gt 0 ]]; then
        echo "Missing subjects from raw file:"
        echo "${missing_subjects[@]}"
    else
        echo "No missing subjects."
    fi

    if ((below_tmp + above_tmp != lines)); then
        echo "Error: mismatch in subject count for site $site."
        echo "Passed: $below_tmp, Failed: $above_tmp, Expected: $lines"
    else
        echo "All subjects for site $site have been correctly classified."
    fi
done
