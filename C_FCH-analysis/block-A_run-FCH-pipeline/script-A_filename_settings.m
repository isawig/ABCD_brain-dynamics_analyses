%% Settings file for functional connectome harmonics estimation
% -------------------------------------------------------------------------
% Script adapted by Isabella M. Wigley and Aurora Berto for the PONS project
% 05 / 2025; ilmawi@utu.fi, aurber@utu.fi
%
% This code runs cortical harmonic decomposition on ABCD data preprocessed 
% with fMRIPrep and XCP-D.
%
% Initial acquisitions were sampled using the 4S156 reference atlas, 
% comprising 100 cortical ROIs from the Schaefer100 atlas (Schaefer et al., 
% 2018) and 56 subcortical ROIs from the COI168 (Pauli et al., 2018), 
% ThalamusHCP (Najdenovska et al., 2018), and SubcorticalHCP 
% (Glasser et al., 2013) atlases.
%
% In a subsequent step, ROIs were relabeled using a reduced atlas consisting 
% of the same 100 cortical ROIs and 14 subcortical ROIs obtained by averaging 
% preselected subcortical regions to match ENIGMA Toolbox requirements 
% (Larivière et al., 2021). These subcortical regions include left and right: 
% accumbens, amygdala, caudate, hippocampus, pallidum, putamen, and thalamus.
%
% Only subjects with a full set of 114 ROIs and a minimum of 380 frames 
% were retained for analysis.
% -------------------------------------------------------------------------

clear; close all; clc


%% Add the source files and scripts for the pipeline

addpath(genpath('/path/to/FCH/FCH_script_ABCDdata/src-fch')) % USER.adapt!

% Define data path
dataPath = '/path/to/FCH/ABCD_data'; % USER.adapt!
addpath(genpath(dataPath)) 

%% Define the data you will be using
% USER: adapt to your data, this is just giving filenames and ID codes that you 
% will use for the fconn matrices

% Define the site
site = "ALL"; % USER.adapt!

% Define the folder where data are stored for the selected site
% folder = sprintf("%s_subjects_for_FCH",site);
folder = "All_subjs_data"; % USER.adapt!
folderPath = fullfile(dataPath, folder);

% Define the .txt file with subjects' IDs
sbjs_ids = "all_subjs_passed_to_LEiDA.txt"; % USER.adapt!


%% For all subjects, extract .txt file's name and save in a cell

% Get a list of all .txt files in the folder
fileStruct = dir(fullfile(folderPath,"*.txt"));

% Extract the names into a cell array
fileNames = {fileStruct.name};

%% For all subjects, extract their ID and save in a cell

% Read each line into a string array
lines = readlines(sbjs_ids);

% Convert the string array to cell array of character vectors
subjectsNames = cellstr(lines);


%% Clear unnecessary variables
clear fileStruct lines 
