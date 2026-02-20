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

clear; close all; clc

%% Read A_filename_settings.m to define the subjects / files
A_filename_settings

%% Define required paths
addpath(genpath('/path/to/FCH/FCH_script_ABCDdata/src-fch')) % USER.adapt!

outputPath = '/path/to/FCH/FCH_results'; % USER.adapt!

%% Load the group average FC

% Define path where the average FC is stored
fcPath = sprintf('/path/to/FCH/FCH_results/D1_fc/%s',site); % USER.adapt!

% Load the group average FC
load(fullfile(fcPath,'FC.mat'))

%% ADJACENCY 

%% Prepare clean save of the outputs

% Define the number of nearest neighbours to create the adjacency
nn = 20;

noH = 113; % number of Harmonics, this is number of ROI's minus 1

% Define adjacency output folder
adjacencyFolder = ['FCH','_nn',num2str(nn),'_H',num2str(noH)];
adjacencyOutput = fullfile(outputPath, adjacencyFolder, site);

if ~exist(adjacencyOutput,"dir")
    mkdir(adjacencyOutput)
end

%% Compute the adjacency

% Compute adjacency (function in src-fch path)
[A] = computeAdjacency(squareform(FC), 'nn', nn, 0);

% Visualization
imagesc(A); % for QC
colorbar
clim([0 1])

% Save
figname = char('GroupAverageAdjacency');
savefig(fullfile(adjacencyOutput,figname))
saveas(gcf,fullfile(adjacencyOutput,'GroupAverageAdjacency.pdf'));
    
%% Save adjacency

save(fullfile(adjacencyOutput, ['Adjacency_basisGroupAverage_knn', num2str(nn), '.mat']), 'A', '-v7.3');
% save('Adjacency.mat', 'A', '-v7.3');

%% HARMONICS 

%% Compute the harmonics

% Define number of vertices (equal to number of ROIs)
nr_vertices = 114; 

% Initialise matrix to store harmonic results
Y = zeros(nr_vertices, size(A,2)-1);

% Harmonic estimation
[Y, V] = computeLaplacianEigenmaps(A, size(A,2)-1, 'symmetric');

%% Save the harmonics file

% Store in the same folder as adjacency results 
harmonicsOutput = adjacencyOutput;

save(fullfile(harmonicsOutput, 'Harmonics.mat'), 'Y', 'V', '-v7.3');
% save('Harmonics.mat', 'Y', 'V', '-v7.3');
% save Harmonics.txt Y -ascii

% Visualize
imagesc(Y); % for QC
colorbar
clim([-0.1 0.3])

% Save
figname = char ('GroupAverageHarmonics');
savefig(fullfile(harmonicsOutput,figname))
saveas(gcf,fullfile(harmonicsOutput,'GroupAverageHarmonics.pdf'));

%% PROJECTIONS

%% Compute Projections prep

% Store in the same folder as other results 
projectionsOutput = adjacencyOutput;

in.method = 'dot'; %'dot';

% read the harmonics
CH = Y;

%% Estimate the projections, i,e. projections of harmonics from the fMRI timecourses

% Loop through each subject
for subs = 1:length(subjectsNames)-1
    
    % Import data for the current subject
    tc_X = importdata(fullfile(folderPath, fileNames{subs}));

    % Normalization (remove mean and divide by SD)
    tc_X = zscore(tc_X'); 

    % Transpose again to make zscore output in required format
    tc_X = tc_X';

    % Estimate projections
    Projections.Data_norm(subs,:,:) = surfFMRI_projectCH(tc_X, CH, [], [], in.method);
end

%% save the results
save(fullfile(projectionsOutput, ['Projections_basisGroupAverage_knn', num2str(nn), '.mat']), 'Projections','-v7.3');
    
%% Print figures of the output, just for checking that it worked

figure, imagesc(squeeze(abs(Projections.Data_norm(1,:,:))))
colorbar
clim([0 1])

figure, imagesc(squeeze(abs(Projections.Data_norm(2,:,:))))
colorbar
clim([0 1])

%% END