%% This script is for a single group of scans (single SITE)
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

% input: individual fconn matrices with z values
% output: mean, SD, CV matrices and group coefficient of variation

close all; clc

%% Define paths 

% Define input path for the selected site
inputPath = sprintf('/path/to/FCH/FCH_results/B_fconn_z/%s',site); % USER.adapt!

% Define output path to store results
resultsPath = '/path/to/FCH/FCH_results/C_QC'; % USER.adapt!
if ~exist(resultsPath,"dir")
    mkdir(resultsPath)
end

%% Define specific paths for the outcome metrics 

% Group mean
groupMeanPath = fullfile(resultsPath, "GroupMean1");
if ~exist(groupMeanPath,"dir")
    mkdir(groupMeanPath)
end

% Standard deviation
stdPath = fullfile(resultsPath, "GroupSD1");
if ~exist(stdPath,"dir")
    mkdir(stdPath)
end

% Coefficient of variation
cvPath = fullfile(resultsPath, "GroupCV1");
if ~exist(cvPath,"dir")
    mkdir(cvPath)
end

% Group mean difference
groupMeanDifferencePath = fullfile(resultsPath, "case2group_MeanDifference");
if ~exist(groupMeanDifferencePath,"dir")
    mkdir(groupMeanDifferencePath)
end

%% Summary matrices

% Create a matrix of all zeros to work as a basis
MatricesGrouped1=zeros(114, 114);

% Loop for reading the matrix data and concatenate 
for subs = 1:length(subjectsNames)-1 % indices for the fnames and IDcodes

    % Import data returned from B_make_fconn_zvalues_loop.m code
    temp_matrix = importdata(fullfile(inputPath, [subjectsNames{subs} '.txt']));

    % Save each subject's FC matrix on the 3rd dimension
    MatricesGrouped1(:,:,subs+1)=temp_matrix;
end 

%remove the first "basis" matrix that is only zeros
MatricesGrouped1 = MatricesGrouped1(:,:,2:end);

% estimate group mean
GroupMean1 = mean(MatricesGrouped1,3);

% estimate standard deviation
GroupSD1 = std(MatricesGrouped1,[],3);

% estimate group coefficient of variation
GroupCV1 = GroupSD1./GroupMean1;

%% Plot and save the summary statistics

%% Visualize group mean
figure
imagesc(GroupMean1)
colorbar
clim([-0.4 0.8])
xlabel('Brain area label'); ylabel ('Brain area label'); axis square; axis tight
title ('Group mean z-scores')

% Save
sitePath = sprintf("%s/%s", groupMeanPath, site);
if ~exist(sitePath,"dir")
    mkdir(sitePath)
end
save(fullfile(sitePath, 'QC1_GroupMean.mat'), 'GroupMean1', '-v7.3');
saveas(gcf, fullfile(sitePath, 'QC1_GroupMean.pdf'));

%% Visualize group standard deviation
figure 
imagesc(GroupSD1)
colorbar
clim([0 0.5])
xlabel('Brain area label'); ylabel ('Brain area label'); axis square; axis tight
title ('Group SD of z scores')

% Save
sitePath = sprintf("%s/%s", stdPath, site);
if ~exist(sitePath,"dir")
    mkdir(sitePath)
end
save(fullfile(sitePath, 'QC2_GroupSD.mat'), 'GroupSD1', '-v7.3');
saveas(gcf, fullfile(sitePath, 'QC2_GroupSD.pdf'));

%% Visualize coefficient of variation
figure
imagesc(GroupCV1)
colorbar
clim([0 2])
xlabel('Brain area label'); ylabel ('Brain area label'); axis square; axis tight
title ('Group coefficient of variation of z scores')

% Save
sitePath = sprintf("%s/%s", cvPath, site);
if ~exist(sitePath,"dir")
    mkdir(sitePath)
end
save(fullfile(sitePath, 'QC3_GroupCV1.mat'), 'GroupCV1', '-v7.3');
saveas(gcf, fullfile(sitePath, 'QC3_GroupCV1.pdf'));

%% Case-to-group covariance 

% Initialize a vector to store the average covariance for each subject
case2group_MeanDifference = zeros(1, size(MatricesGrouped1,3));

% Loop through each subject
for w=1:size(MatricesGrouped1,3)

    % Initialize vectors to store correlation values and p-values
    % between the current subject and all other subjects
    case2group_Covariances_rho = zeros(1, size(MatricesGrouped1,3)-1);
    case2group_Covariances_pval = zeros(1, size(MatricesGrouped1,3)-1);

    % Compare the current subject with all others
    for w2=1:size(MatricesGrouped1,3)
        if w==w2
           continue % skip self-comparison
        end 

        % Compute correlation and p-value between the two matrices
        [RHO, PVAL] = corrcoef(MatricesGrouped1(:,:,w), MatricesGrouped1(:,:,w2));

        % Store the correlation coefficient and p-value
        case2group_Covariances_rho(1,w2) = RHO(1,2);
        case2group_Covariances_pval(1,w2) = PVAL(1,2);
    end
    
    % Save the mean correlation between the current subject and all others
    case2group_MeanDifference(1, w) = mean(case2group_Covariances_rho);
end

%% plot the covariance as a histogram
figure 
histogram(case2group_MeanDifference)
xlabel('correlation values'); ylabel ('no. of participants'); axis square; axis tight
title ('case-to-group correlations')

% Save
sitePath = sprintf("%s/%s", groupMeanDifferencePath, site);
if ~exist(sitePath,"dir")
    mkdir(sitePath)
end
save(fullfile(sitePath, 'QC4_case2group_MeanDifference.mat'), 'case2group_MeanDifference', '-v7.3');
saveas(gcf, fullfile(sitePath, 'QC4_case2group_MeanDifference.pdf'));

%% END


