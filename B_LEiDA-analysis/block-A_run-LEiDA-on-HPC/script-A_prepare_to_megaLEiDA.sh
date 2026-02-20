#!/bin/bash
## Script to prepare data for LEiDA run on the pooled dataset.
# Save data from all subjects that passed LEiDA preprocessing in the same folder in PUHTI.

# Isabella L.C. Mariani Wigley, 04 / 2025; ilmawi@utu.fi
# Aurora Berto, 04 / 2025; aurber@utu.fi

# Sites
sites=("G010" "G031" "G032" "G075" "G087" "S011" "P023" "P043" "P064" "S012" "S013" "S014" "S020" "S021" "S022" "S042" "S053" "S065" "S076" "S086" "S090") # USER.adapt!

# Output directory, where files will be stored
output_path="/path/to/data/pooledLEiDA/All_subjs_data" # USER.adapt!
mkdir -p "$output_path" # create it if not existing


## Loop over each site
for site in ${sites[@]}; do
    # Set the path of the actual site
    site_path="/path/to/data/${site}/temp/114_rows_380_columns" # USER.adapt!
    passed_subjs="/path/to/data/preparing_for_LEiDA/subjs_passed_to_LEiDA/${site}_subjs_passed_to_LEiDA.txt" # USER.adapt!

    ## Loop over each subject
    for subj in "$site_path"/sub-*; do

        # Get fullname
        filename=$(basename "$subj")
        # Extract subject name
        subj_name=${filename%%_ses-baseline*}

        # Define the path where to store the file and add site in front of the site's name
        output_file="$output_path/${site}_${subj_name}.txt"

        # If the subject has already been downloaded
        if [ -f "$output_file" ]; then
            echo "Subject $subj_name already copied. Skipping."
        
        else
            # Double check: if the subject is in the list of passed subjects
            if grep -q "$subj_name" "$passed_subjs"; then
                # Copy the file in the output folder
                cp "$subj" "$output_file"
                echo "Copied $subj_name to $output_file"
            else
                echo "Subject $subj_name not in list of passed subjects for site $site"
            fi
        fi
    done
done
