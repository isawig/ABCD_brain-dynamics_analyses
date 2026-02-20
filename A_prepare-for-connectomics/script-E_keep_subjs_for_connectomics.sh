#!/bin/bash
# Code that runs last four steps to prepare data for connectome analysis:
# 4- For each subject, count the number of ROIs (rows) and display it
# 5- For each subject, count the number of samples (columns) and display it
# 6- Keep only subjects with 114 ROIs and store others in a separate file
# 7- Keep only subjects with at least 380 samples and store others in a separate file

# Code to run on SLURM using sbatch command.

# Isabella L.C. Mariani Wigley, 03 / 2025; ilmawi@utu.fi
# Aurora Berto, 03 / 2025; aurber@utu.fi

#SBATCH --job-name=keep_subjs    # Job name
#SBATCH --output=keep_subjs.out  # Standard output and error log
#SBATCH --error=keep_subjs.err   # Error log
#SBATCH --ntasks=1               # Run a single task
#SBATCH --time=00:10:00          # Time limit hrs:min:sec
#SBATCH --mem=1G                 # Memory limit

# Load any required modules here (if necessary)
# module load <module_name>


###### Site #######
site="S090" # USER.adapt!


### Directories ###
# Input directory 
input_dir="/path/to/data/${site}/temp/ready_for_leida" # USER.adapt!

# Output directory for subjects to keep for LEiDA
target_dir="/path/to/data/${site}/temp/114_rows_380_columns" # USER.adapt!
mkdir -p "$target_dir"

# File to store subjects ID
output_leida="/path/to/data/preparing_for_LEiDA/subjs_passed_to_LEiDA" # USER.adapt!
mkdir -p "$output_leida"
file_leida="$output_leida/${site}_subjs_passed_to_LEiDA.txt"
> "$file_leida"

# Output directories and files for subjects with missing ROIs or frames
output_rows="/path/to/data/preparing_for_LEiDA/subjs_missing_rows" # USER.adapt!
mkdir -p "$output_rows"
file_rows="$output_rows/${site}_subjs_missing_rows.txt"
> "$file_rows"

output_columns="/path/to/data/preparing_for_LEiDA/subjs_missing_columns" # USER.adapt!
mkdir -p "$output_columns"
file_columns="$output_columns/${site}_subjs_missing_columns.txt"
> "$file_columns"


### Iterate for each subject ###
# Initialize a dictionary to store row and column counts
row_counts=()
column_counts=()

# Initialize variables to store discarded and confirmed subjects
sbjs_no_rows=()
sbjs_no_cols=()
sbjs_to_leida=()

# Loop over each file in the directory
for file in "$input_dir"/*.txt; do
  subj_name=$(basename "$file" .txt)

  ### Count ROIs ###
  # Count the number of rows in the current file
  num_rows=$(wc -l < "$file")
  # Increment the count for this specific number of rows
  ((row_counts[$num_rows]++))


  ### Count frames ###
  # Count the number of columns in the first line of the current file (assuming space-delimited)
  num_columns=$(awk '{print NF; exit}' "$file")
  # Increment the count for this specific number of columns
  ((column_counts[$num_columns]++))


  ### Select subjects to keep ###
  # Check the number of rows
  if [ ! "$num_rows" -eq 114 ]; then # USER.adapt!
    echo "$subj_name" | cut -d'_' -f1 >> "$file_rows" # list subject's name
    ((sbjs_no_rows++))
  else
    # Check the number of columns
    if [ "$num_columns" -lt 380 ]; then # USER.adapt!
      echo "$subj_name" | cut -d'_' -f1 >> "$file_columns" # list subject's name
      ((sbjs_no_cols++))
    else 
      # If both are correct, copy in the target directory
      cut -f1-380 "$file" > "$target_dir/${subj_name}_ses-baselineYear1Arm1_task-rest_space-MNI152NLin6Asym_seg-4S156Parcels_stat-mean_timeseries_grouped_transposed.txt" # USER.adapt!
      echo "$subj_name" | cut -d'_' -f1 >> "$file_leida"
      ((sbjs_to_leida++))
    fi
  fi
done

### Display the results ###
# ROIs
echo "Row count summary for all subjects:"
for row_count in "${!row_counts[@]}"; do
  echo "Files with $row_count rows: ${row_counts[$row_count]}"
done

# Frames
echo "Column count summary for all subjects:"
for column_count in "${!column_counts[@]}"; do
  echo "Files with $column_count columns: ${column_counts[$column_count]}"
done

# Counts
echo "Subjects passed to LEiDA: ${sbjs_to_leida}"
echo "Subjects discarded due to wrong number of rows: ${sbjs_no_rows}"
echo "Subjects discarded due to wrong number of columns: ${sbjs_no_cols}"


# Final message
echo "Files with 114 rows and 380 columns have been copied to $target_dir."
echo "Files with missing rows have been saved in $output_rows."
echo "Files with missing columns have been saved in $output_columns."


