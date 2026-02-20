%% Estimate harmonics power
clear; close all; clc

% Isabella L.C. Mariani Wigley, 09 / 2025; ilmawi@utu.fi
% Aurora Berto, 09 / 2025; aurber@utu.fi

%% Set paths and parameters
% load the harmonics projections
input_path = '/path/to/FCH/FCH_results/FCH_nn10_H113/ALL'; % USER.adapt!
load(fullfile(input_path,'Projections_basisGroupAverage_knn10.mat')) % USER.adapt!

% output path to store results
output_path = '/path/to/FCH/FCH_derive_measures/Results'; % USER.adapt!

% define the sample size
no_of_subjects = 6624; % USER.adapt!
    
% initialise the output file    
Harmonics_power = zeros(no_of_subjects,1);

%% Loop over 113 harmonics and calculate power 
% loop and save the values of interest    
for HAR = 1:113 
    
    % harmonics variable X is the dependent variable
    xprep = squeeze(Projections.Data_norm(:,:,HAR+1)); % discard the first harmonic (null eigenvalue)
    xabs = abs(xprep); % absolute values to get the power 
    x = mean(xabs,2); % mean across time is the dependedent variable
    
    % add the values to the output file
    Harmonics_power = [Harmonics_power x];
    
end

% get rid of the zero column used to initialise
Harmonics_power = Harmonics_power(:,2:end);

% save the results as .mat file for later use
save(fullfile(output_path,'Extracted_Harmonics_power'), 'Harmonics_power', '-v7.3');
% prepare a table to save
output = array2table(Harmonics_power);
% save the table in csv format
writetable(output,fullfile(output_path,'Extracted_Harmonics_power.csv'))

%% visualise and save the output

R_power = corrcoef(Harmonics_power);
figure
imagesc(R_power)
set(gca,'CLim',[-0.2 1]);
colorbar
title('Extracted_Harmonics_power_correlation')

figdir = output_path;
    baseFileName =  sprintf('Extracted_Harmonics_power_correlation.png');
    fullFileName = fullfile(figdir, baseFileName);
    saveas(gcf, fullFileName);

%% NED