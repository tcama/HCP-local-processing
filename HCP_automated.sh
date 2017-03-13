#!/bin/bash
# This script is used to execute all necessary command for preprocessing and 
# processing of a subject from HCP. The script should be run as follows:
#
# ./HCP_automated ******
#
# Where the stars are replaced by the HCP subject ID code: 
#
# ./HCP_automated 100206
#
# This script should be executed from the directory containing all of the 
# necessary processing scripts, or should be altered to properly find the
# location of the needed scripts. There are four resting-state fmri (rsfmri) 
# scans, two for each phase encoding used: right-left (RL) and left-right (LR). 
# This code concatenates the two scans for each phase encoding. 

# set the number of TRs and ROIs included
TRs=2400 # both phase encodings, each 1200 TRs long
ROIs=627 # determined by atlas

# Make directories for data
homedir=$PWD
mkdir "$homedir"/HCP/
mkdir "$homedir"/HCP/data/
mkdir "$homedir"/HCP/data/"$1"
mkdir "$homedir"/HCP/data/"$1"/rfMRI_REST1_RL/
mkdir "$homedir"/HCP/data/"$1"/rfMRI_REST2_RL/
mkdir "$homedir"/HCP/data/"$1"/rfMRI_REST1_LR/
mkdir "$homedir"/HCP/data/"$1"/rfMRI_REST2_LR/

# make directory for holding results
mkdir "$homedir"/HCP/results/
mkdir "$homedir"/HCP/results/"$1"

# Download data via AWS and gunzip (1_RL)
cd "$homedir"/HCP/data/"$1"/rfMRI_REST1_RL/
aws s3 cp s3://hcp-openaccess/HCP_900/"$1"/MNINonLinear/Results/rfMRI_REST1_RL/RibbonVolumeToSurfaceMapping/goodvoxels.nii.gz goodvoxels.nii.gz
aws s3 cp s3://hcp-openaccess/HCP_900/"$1"/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_hp2000_clean.nii.gz rfMRI_REST1_RL_hp2000_clean.nii.gz
gunzip rfMRI_REST1_RL_hp2000_clean.nii.gz
gunzip goodvoxels.nii.gz

# copy needed functions to data directory and begin processing (1_RL)
cd "$homedir"/HCP/data/"$1"/rfMRI_REST1_RL/
cp "$homedir"/AAL627_10202016.nii "$homedir"/HCP/data/"$1"/rfMRI_REST1_RL/
"$homedir"/roi_preprocess rfMRI_REST1_RL_hp2000_clean.nii AAL627_10202016.nii goodvoxels.nii
mv "$homedir"/HCP/data/$1/rfMRI_REST1_RL/goodvoxels_new.csv "$homedir"/HCP/results/"$1"/goodvoxels_new1_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_RL/mean_data.csv "$homedir"/HCP/results/"$1"/mean_data1_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_RL/snr_data.csv "$homedir"/HCP/results/"$1"/snr_data1_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_RL/roi_data.csv "$homedir"/HCP/results/"$1"/roi_data1_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_RL/mc.csv "$homedir"/HCP/results/"$1"/meancenter1_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_RL/voxels_summary.csv "$homedir"/HCP/results/"$1"/voxels_summary1_RL.csv

# Download data via AWS and gunzip (2_RL)
cd "$homedir"/HCP/data/"$1"/rfMRI_REST2_RL/
aws s3 cp s3://hcp-openaccess/HCP_900/"$1"/MNINonLinear/Results/rfMRI_REST2_RL/RibbonVolumeToSurfaceMapping/goodvoxels.nii.gz goodvoxels.nii.gz
aws s3 cp s3://hcp-openaccess/HCP_900/"$1"/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_RL_hp2000_clean.nii.gz rfMRI_REST2_RL_hp2000_clean.nii.gz
gunzip rfMRI_REST2_RL_hp2000_clean.nii.gz
gunzip goodvoxels.nii.gz

# copy needed functions to data directory and begin processing (2_RL)
cd "$homedir"/HCP/data/"$1"/rfMRI_REST2_RL/
cp "$homedir"/AAL627_10202016.nii "$homedir"/HCP/data/"$1"/rfMRI_REST2_RL/
"$homedir"/roi_preprocess rfMRI_REST2_RL_hp2000_clean.nii AAL627_10202016.nii goodvoxels.nii
mv "$homedir"/HCP/data/$1/rfMRI_REST2_RL/goodvoxels_new.csv "$homedir"/HCP/results/"$1"/goodvoxels_new2_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_RL/mean_data.csv "$homedir"/HCP/results/"$1"/mean_data2_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_RL/snr_data.csv "$homedir"/HCP/results/"$1"/snr_data2_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_RL/roi_data.csv "$homedir"/HCP/results/"$1"/roi_data2_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_RL/mc.csv "$homedir"/HCP/results/"$1"/meancenter2_RL.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_RL/voxels_summary.csv "$homedir"/HCP/results/"$1"/voxels_summary2_RL.csv

# change to resutls and concatenate processed results
cd "$homedir"/HCP/results/"$1"/
cat meancenter1_RL.csv meancenter2_RL.csv > meancenter_RL.csv

# Download data via AWS and gunzip (1_LR)
cd "$homedir"/HCP/data/"$1"/rfMRI_REST1_LR/
aws s3 cp s3://hcp-openaccess/HCP_900/"$1"/MNINonLinear/Results/rfMRI_REST1_LR/RibbonVolumeToSurfaceMapping/goodvoxels.nii.gz goodvoxels.nii.gz
aws s3 cp s3://hcp-openaccess/HCP_900/"$1"/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz rfMRI_REST1_LR_hp2000_clean.nii.gz
gunzip rfMRI_REST1_LR_hp2000_clean.nii.gz
gunzip goodvoxels.nii.gz

# copy needed functions to data directory and begin processing (1_LR)
cd "$homedir"/HCP/data/"$1"/rfMRI_REST1_LR/
cp "$homedir"/AAL627_10202016.nii "$homedir"/HCP/data/"$1"/rfMRI_REST1_LR/
"$homedir"/roi_preprocess rfMRI_REST1_LR_hp2000_clean.nii AAL627_10202016.nii goodvoxels.nii
mv "$homedir"/HCP/data/$1/rfMRI_REST1_LR/goodvoxels_new.csv "$homedir"/HCP/results/"$1"/goodvoxels_new1_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_LR/mean_data.csv "$homedir"/HCP/results/"$1"/mean_data1_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_LR/snr_data.csv "$homedir"/HCP/results/"$1"/snr_data1_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_LR/roi_data.csv "$homedir"/HCP/results/"$1"/roi_data1_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_LR/mc.csv "$homedir"/HCP/results/"$1"/meancenter1_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST1_LR/voxels_summary.csv "$homedir"/HCP/results/"$1"/voxels_summary1_LR.csv

# Download data via AWS and gunzip (2_LR)
cd "$homedir"/HCP/data/"$1"/rfMRI_REST2_LR/
aws s3 cp s3://hcp-openaccess/HCP_900/"$1"/MNINonLinear/Results/rfMRI_REST2_LR/RibbonVolumeToSurfaceMapping/goodvoxels.nii.gz goodvoxels.nii.gz
aws s3 cp s3://hcp-openaccess/HCP_900/"$1"/MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_LR_hp2000_clean.nii.gz rfMRI_REST2_LR_hp2000_clean.nii.gz
gunzip rfMRI_REST2_LR_hp2000_clean.nii.gz
gunzip goodvoxels.nii.gz

# copy needed functions to data directory and begin processing (2_LR)
cd "$homedir"/HCP/data/"$1"/rfMRI_REST2_LR/
cp "$homedir"/AAL627_10202016.nii "$homedir"/HCP/data/"$1"/rfMRI_REST2_LR/
"$homedir"/roi_preprocess rfMRI_REST2_LR_hp2000_clean.nii AAL627_10202016.nii goodvoxels.nii
mv "$homedir"/HCP/data/$1/rfMRI_REST2_LR/goodvoxels_new.csv "$homedir"/HCP/results/"$1"/goodvoxels_new2_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_LR/mean_data.csv "$homedir"/HCP/results/"$1"/mean_data2_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_LR/snr_data.csv "$homedir"/HCP/results/"$1"/snr_data2_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_LR/roi_data.csv "$homedir"/HCP/results/"$1"/roi_data2_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_LR/mc.csv "$homedir"/HCP/results/"$1"/meancenter2_LR.csv
mv "$homedir"/HCP/data/$1/rfMRI_REST2_LR/voxels_summary.csv "$homedir"/HCP/results/"$1"/voxels_summary2_LR.csv

# concatenate meancentered data
cd "$homedir"/HCP/results/"$1"/
cat meancenter1_LR.csv meancenter2_LR.csv > meancenter_LR.csv

# calculate correlation matrices
"$homedir"/fcMat_generation meancenter_RL.csv $ROIs $TRs
mv fcMap.csv fcMap_RL.csv
"$homedir"/fcMat_generation meancenter_LR.csv $ROIs $TRs
mv fcMap.csv fcMap_LR.csv

# remove extraneous data files
rm goodvoxels* meancenter*

# move data files to HCP_raw_data and delete zipped downloads
rm -R "$homedir"/HCP/data/"$1"
