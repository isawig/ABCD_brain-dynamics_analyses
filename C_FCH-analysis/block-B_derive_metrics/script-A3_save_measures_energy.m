%% Code to calculate energy from harmonics projections and eigenvalues

% Isabella L.C. Mariani Wigley, 09 / 2025; ilmawi@utu.fi
% Aurora Berto, 09 / 2025; aurber@utu.fi

clear; close all; clc

%% Set paths and load files

% load the harmonics projections
input_path = '/path/to/FCH/FCH_results/FCH_nn10_H113/ALL'; % USER.adapt!
load(fullfile(input_path,'Projections_basisGroupAverage_knn10.mat')) 
load(fullfile(input_path,'Harmonics.mat'))

% output path to store results
output_path = '/path/to/FCH/FCH_derive_measures/Results'; % USER.adapt!

%% Define parameters 
% define the sample size
no_of_subjects = 6624;
    
% initialise the output file    
Harmonics_energy = zeros(no_of_subjects,1);

% From Selen    
names = fieldnames(Projections); % subject x timepoints x harmonics

%% Estimate energy for all subjects and harmonics
% loop and save the values of interest    
for HAR = 1:113 

    % harmonics variable X is the dependent variable
    xprep = squeeze(abs(Projections.Data_norm(:,:,HAR+1))); % exclude first harmonic (null eigenvalue)

    % initialise the matrix
        xe = [];
        % for loop to estimate energy across time
        for i = 1:no_of_subjects
            M = ((xprep(i,:)) .* V(HAR+1,1)).^2; 
            ET = sum(M);
            xe = [xe,ET];
        end
     % transpose to column vector to enable statistics
     xe = xe';     
     % add the values to the output file
     Harmonics_energy = [Harmonics_energy xe];
end

% get rid of the zero column used to initialise
Harmonics_energy = Harmonics_energy(:,2:end);

%% Save resulting files

% save the results as .mat file for later use
save(fullfile(output_path,'Extracted_Harmonics_energy'), 'Harmonics_energy', '-v7.3');

% prepare a table to save
output = array2table(Harmonics_energy);
% save the table in csv format
writetable(output,fullfile(output_path,'Extracted_Harmonics_energy.csv'))