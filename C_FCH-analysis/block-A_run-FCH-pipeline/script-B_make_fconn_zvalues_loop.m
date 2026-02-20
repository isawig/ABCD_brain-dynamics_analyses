%% Script to estimate functional connectome matrices 
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

site = "ALL"; % USER.adapt!

% Define output to store fc_X results and images
resultsPath = '/path/to/FCH/FCH_results/B_fconn_z'; % USER.adapt!
if ~exist(resultsPath,"dir")
    mkdir(resultsPath)
end

% Define the folder for the selected site
sitePath = sprintf("%s/%s", resultsPath, site);
if ~exist(sitePath,"dir")
    mkdir(sitePath)
end

% Define path to store figures
figuresPath = fullfile(sitePath, "fconn_z_figures"); 
if ~exist(figuresPath,"dir")
    mkdir(figuresPath)
end

%% loop for reading in the timecourses and making the fc

% Iterate over each subject
for subs = 1:length(subjectsNames)-1 % indices for the fnames and IDcodes
    
    % Import data
    tc_X = importdata(fullfile(folderPath, fileNames{subs}));
    
    % Select only the desired ROIs (commented since all are selected)
    % tc_X_114 = tc_X(1:114,:);

    % Calculate functional connectome matrix
    fc_X = script_make_fcz_114(tc_X);
    
    % Define output file's name
    outputFile = fullfile(sitePath, [subjectsNames{subs} '.txt']);

    % Save .txt file
    writematrix(fc_X, outputFile, 'Delimiter','tab')
    
    % Visualize FC matrix
    imagesc(fc_X)
    colorbar
    clim([-0.4 0.8])

    % Save figure
    figname = char(subjectsNames{subs});
    savefig(fullfile(figuresPath,figname))
end

%% Clear unnecessary variables
% Variables
clear fc_X figname folder outputFile sbjs_ids subs tc_X

% Paths
clear dataPath figuresPath resultsPath sitePath

%% END
