%% SCRIPT to run FCH pipeline on Schaefer100 + 14 subcortical atlas  
% Jetro J. Tuulari, 01 / 2022

% -------------------------------------------------------------------------
% Script adapted by Isabella M. Wigley and Aurora Berto for the PONS project
% 05 / 2025; aurber@utu.fi, ilmawi@utu.fi
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

%% Clean up in case there are QC runs before the modelling

clear; close all; clc

%% Read A_filename_settings.m to define the subjects / files

A_filename_settings

%% Define output path for the current site

outputPath = sprintf('/path/to/FCH/FCH_results/D1_fc/%s',site); % USER.adapt!
if ~exist(outputPath,"dir")
    mkdir(outputPath)
end

%% Import data and create a cell for time courses files "TC"

% Clear TC variable if present, to redefine it
clear TC

% Loop through each subject
for iSubs = 1:length(subjectsNames)-1 

    % Import current subject's data
    tc_tmp = importdata(fullfile(folderPath, fileNames{iSubs}));

    % Store data in TC cell
    TC{iSubs} = tc_tmp;
end

%% Create standardised time courses, combined for all participants

% Initialize matrices
tc_tmp = []; T = [];

% Loop through each subject
for iSubs = 1:length(TC)

    % Select subject's data
    tc_tmp = TC{iSubs};

    % Standardization (remove mean and divide by SD)
    tc_tmp = tc_tmp - repmat(mean(tc_tmp,2),1,size(tc_tmp,2));
    tc_tmp = tc_tmp ./ repmat(std(tc_tmp,[],2),1,size(tc_tmp,2));
    
    % Store standardized results in TC_stand
    TC_stand{iSubs} = tc_tmp;

    % Store frame length in T
    T(iSubs) = size(tc_tmp,2);
end

% Convert TC_stand from cell to matrix 
% (all subjects are being concatenated on the 2nd dimension)
tc_matrix = cell2mat(TC_stand);

% Assing the mean value to zero voxels
tc_matrix(tc_matrix == 0) = mean(mean(tc_matrix));

%% Estimate the mean FC matrix

% Calculate mean FC
FC = pdist(tc_matrix, 'correlation');
FCmatrix = squareform(FC);

% Visualize and save
imagesc(FCmatrix); % for QC
colorbar
clim([0 1])
figname = char('GroupAverageFC');
savefig(fullfile(outputPath, figname))

%% Save the group average FC

save(fullfile(outputPath,'FC.mat'), 'FC', '-v7.3'); 

%% END