#!bin/bash
# Code to check that all subjects have been succesfully downloaded.

# Isabella L.C. Mariani Wigley, 03 / 2025; ilmawi@utu.fi
# Aurora Berto, 03 / 2025; aurber@utu.fi

site="S012" # USER.adapt!

# Load the list of subjects that passed the QC -> subjects to download
subjs_QC_file="/path/to/subjs_passed_QC/${site}_fd_below.txt" # USER.adapt!

# Output directory to store non-downloaded subjects
output_file="/path/to/store/data/${site}/xcpd_noncifti/${site}_non_downloaded_subjs.txt" # USER.adapt!
> "$output_file"

# Use the list in a for loop
while read -r subj; do

 # Extract the subject name from the list
  subject_name="sub-NDAR$subj"
  
  # Define the local path where the subject will be downloaded
  local_path="/path/to/store/data/${site}/xcpd_noncifti/" # USER.adapt!
  

  # Check if the subject has already been downloaded
  if [ -e "$local_path/$subject_name" ]; then
    echo "Subject $subject_name already downloaded. Skipping."
  else
    echo "Non downloaded subject. Saving."
    echo "$subj" >> "$output_file"
  fi

done < "$subjs_QC_file"



