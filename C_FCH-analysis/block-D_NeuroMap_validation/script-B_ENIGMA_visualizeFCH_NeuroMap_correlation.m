%% ENIGMA script to plot FCH validation with NeuroMap

% Visualizes a table with first row corresponding to seven reference
% gradients by NeuroMap, and the first column representing the first six
% low-frequency harmonics. Then each cell is filled with correlation value
% estimated between each pair of harmonic and gradient (bolded if p <
% 0.05), calculated in script code_A_compare_FCH_to_gradients.m

% Isabella L.C. Mariani Wigley, 09 / 2025; ilmawi@utu.fi
% Aurora Berto, 09 / 2025; aurber@utu.fi

clear; close all; clc

%% Number nn
nn = 10;

%% Load data and set path
addpath(genpath('/path/to/ENIGMA/matlab')) % USER.adapt!

input_path = sprintf('/path/to/FCH/FCH_results/FCH_nn%d_H113/ALL',nn); % USER.adapt!
load(fullfile(input_path,"Harmonics.mat"))

%% Output path
% Output folders: where figures will be saved
output_dir = sprintf('/path/to/FCH/FCH_results/NeuroMap_correlation'); % USER.adapt!
output_cortical = fullfile(output_dir, sprintf('/Cortical/tmp/H%d',nn)); % USER.adapt!
output_subcortical = fullfile(output_dir, sprintf('/Subcortical/tmp/H%d',nn)); % USER.adapt!

% Create output folders to store figures
if ~exist(output_cortical,'dir')
    mkdir(output_cortical)
end

if ~exist(output_subcortical,'dir')
    mkdir(output_subcortical)
end

%% Set parameters

noh = 10; % number of harmonics to visualize
noh_names = strcat('H',string(noh));

N_cortical_areas = 100;
% order_cort = [1:2:N_cortical_areas N_cortical_areas:-2:2];

%% 1. Plot brain states

%% Save cortical figures
% for each harmonics eigenvector
for h = 1:noh

    % Set the harmonic name
    H_tmp = sprintf("H%d",h);

    % Select the column 
    V = Y(1:100, h+1);
    % V_ordered = V(:,order_cort);
    % V(V<0) = 0;

    % Map parcelled data to the cortical surface
    Schaefer100 = parcel_to_surface(V,'schaefer_100_conte69');

    % Project the results on the surface brain
    f = figure('Visible', 'off'); % not displayed
    [a, cb] = plot_cortical(Schaefer100, 'surface_name', 'conte69','color_range',[-.2 .2]);
    delete(findall(gcf,'Type','ColorBar'))

    % Retain only the two external views 
    delete(a([2, 3]));
    external_views = a([1, 4]);

    % Put them closer
    set(external_views(1), 'Position', [0.12, 0, 0.35, 0.8]); 
    set(external_views(2), 'Position', [0.53, 0, 0.35, 0.8]); 

    % Delete left and inner cortical visualizations and colorbar
    delete(a([2, 3]));
    delete(cb);

    % Set the right hemisphere in the center and zoom it
    set(a(1, 4)); % 'Position', [0.2, 0.1, 0.6, 0.8])

    % Save the figure
    file_name = sprintf('ext_H%d_n%d.jpg',noh,h);
    full_path = fullfile(output_cortical, file_name);
    saveas(f, full_path)
    % pause()
    close(f)

end

%% Define the matrix of correlations

% Define path
file_path = '/path/to/FCH/FCH_NeuroMap'; % USER.adapt!
file_name = 'neuroMap_FCH_correlation_summary_2.xlsx'; % USER.adapt!

% Read table with correlation values
T = readtable(fullfile(file_path, file_name));

% Gradients to keep
grads_keep = {'Harmonic','Metric','FunctionalGradient_1', 'FunctionalGradient_2', 'FunctionalGradient_3', 'FunctionalGradient_4', 'AllometricScaling_NIH_', 'PC1Gene', 'PC1NeuroSynth'};
T = T(1:12, grads_keep);

corr_cols = find(strcmp(T.Metric,'correlation'));
pvalue_cols = find(strcmp(T.Metric,'p-value'));


% Define titles
titles = {
    ['Functional', newline, 'Gradient #1'],
    ['Functional', newline, 'Gradient #2'],
    ['Functional', newline, 'Gradient #3'],
    ['Functional', newline, 'Gradient #4'],
    ['Allometric Scaling', newline, '(NIH)'],
    ['PC1', newline, 'Gene expression'],
    ['PC1', newline, 'NeuroSynth']
};

%% Define reference networks' images
% Define path
ref_img_path = '/path/to/FCH/FCH_NeuroMap/neuroMap_reference_networks'; % USER.adapt!
ref_img_names = {'gradient1.png','gradient2.png','gradient3.png','gradient4.png','allometricScaling.png','geneExpression.png','neuroSynth.png'}; % USER.adapt!

%% Visualization of cortical images in the same figures
addpath(genpath(output_cortical))

% Dimension of the table inside the figure
n_psi = height(T)/2;  % corr and pvalue for each harmonic
n_grad = width(T)-2;  % discard 'harmonic' and 'metric' cols

scale = 1;
figure('Position', [100, 100, scale*(200*n_grad + 200), scale*(200*(n_psi+1))]); % +1 row for reference image

%%%% First row %%%%
% First row's parameters (with only titles)
top_row_height = 0.7 / (n_psi+1);
left_img = 0.05;
height_img = top_row_height;
width_img = 0.18;

% Parameters for reference images (from 2nd row on)
width_small = 0.7 / n_grad; % the ref image has the same width as tabular data
left_start = left_img + width_img; 

% Reference images + gradient name (1st row)
for j = 1:n_grad
    left_pos = left_start + (j-1)*width_small;
    axes('Position', [left_pos, 1 - height_img - 0.1, width_small*0.9, height_img]);
    img_ref = imread(ref_img_names{j});
    imshow(img_ref);
    axis off;
    title(titles{j}, 'Interpreter', 'tex', 'FontSize', 16);
end

%%%% Following rows %%%%
for i = 1:n_psi
    bottom_img = 1 - (i+1)/ (n_psi + 1); % i+1 because first row is with ref images
    
    % 1. Harmonic image
    axes('Position', [left_img, bottom_img, width_img, height_img]);
    img = imread(sprintf('ext_H%d_n%d.jpg', noh, i));
    imshow(img(350:700,100:1050,:)); 

    text(0, 0.5, sprintf('\\psi_{%d}', i), ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 20,'FontWeight', 'bold');
    
    % 2. Correlation and p-value
    for j = 1:n_grad
        left_text = left_start + (j-1)*width_small;
        axes('Position', [left_text, bottom_img, width_small*0.9, height_img]); 
        
        corr_idx = corr_cols(i); 
        pvalue_idx = pvalue_cols(i);

        if T{pvalue_idx, j+2} < 0.05
            text(0.5, 0.5, sprintf('r = %.2f\np = %.3f', T{corr_idx,j+2}, T{pvalue_idx, j+2}), ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 18, 'FontWeight', 'bold');
        else
            text(0.5, 0.5, sprintf('r = %.2f\np = %.3f', T{corr_idx,j+2}, T{pvalue_idx, j+2}), ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 18);
        end
        axis off;
    end
end
