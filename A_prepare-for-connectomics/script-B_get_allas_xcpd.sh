#!bin/bash
# Code to download data from Allas 

# This code allows the download of a specific file from Allas, without the need to download the full set of files.
# It takes as input the list of subjects to download. 
# It downloads mean time-series for each subject, obtained after fMRIPrep and XCP-D preprocessing.

# Isabella L.C. Mariani Wigley, 03 / 2025; ilmawi@utu.fi
# Aurora Berto, 03 / 2025; aurber@utu.fi

# Specify site ID
site="S053" # USER.adapt!

# Load the list of subjects that passed the QC -> subjects to download
subjs_QC_file="/path/to/subjs_passed_QC/${site}_fd_below.txt" # USER.adapt!

# Use the list in a for loop
while read -r subj; do

 # Extract the subject name from the list
  subject_name="sub-NDAR$subj" # USER.adapt!
  
  # Define the local path where the subject will be downloaded
  local_path="/path/to/store/data/${site}/xcpd_noncifti/" # USER.adapt!
  

  # Check if the subject has already been downloaded
  if [ -e "$local_path/$subject_name" ]; then
    echo "Subject $subject_name already downloaded. Skipping."
  else
    echo "Downloading $subject_name..."
    
    # download the .tar file to local computer without extracting the .tar packing
    swift download 2001640-puhti-SCRATCH ABCD_fMRI/${site}/xcpd_36P_noncifti/$subject_name.tar -o $subject_name.tar # USER.adapt!
    # a-get --asis 2001640-puhti-SCRATCH/ABCD_fMRI/${site}/xcpd_36P_noncifti/$subject_name.tar
    # list the content of tar file
    # tar --list -f $subject_name.tar 
    
    # extract only the specific file
    file_needed=($(tar --list -f $subject_name.tar | grep '4S156Parcels_stat-mean_timeseries.tsv')) # USER.adapt!
    tar -x $file_needed -f $subject_name.tar # --target_dir $local_path
    rm $subject_name.tar
  fi

done < "$subjs_QC_file"