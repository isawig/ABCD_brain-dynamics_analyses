%% ENIGMA script to plot FCH 
% already ordered

% Script by Isabella M. Wigley and Aurora Berto
% 05 / 2025; aurber@utu.fi, ilmawi@utu.fi

clear; close all; clc

%% Number nn
nn = 20; % USER.adapt!

%% Load data and set path
addpath(genpath('/path/to/ENIGMA/matlab')) % USER.adapt!

input_path = sprintf('/path/to/FCH/FCH_results/FCH_nn%d_H113/ALL',nn); % USER.adapt!
load(fullfile(input_path,"Harmonics.mat"))

%% Output path
% Output folders: where figures will be saved
output_dir = sprintf('/path/to/FCH/FCH_results/Connectomic_maps'); % USER.adapt!
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

noh = 6; % number of harmonics to visualize
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
    % delete(findall(gcf,'Type','ColorBar'))

    % Save the figure
    file_name = sprintf('ext_H%d_n%d.jpg',noh,h);
    full_path = fullfile(output_cortical, file_name);
    saveas(f, full_path)
    % pause()
    close(f)

end

%% Save subcortical figures
% for each harmonics eigenvector
for h = 1:noh

    % Set the harmonic name
    H_tmp = sprintf("H%d",h);

    % Select the column 
    V = Y(101:114, h+1);
    % V_ordered = V(:,order_cort);
    % V(V<0) = 0;

    % Project the results on the surface brain
    f = figure('Visible', 'off'); % not displayed
    [a, cb] = plot_subcortical(V, 'color_range',[-.15 .15],'ventricles', 'False');
    % delete(findall(gcf,'Type','ColorBar'))

    % Save the figure
    file_name = sprintf('ext_H%d_n%d.jpg',noh, h);
    full_path = fullfile(output_subcortical, file_name);
    saveas(f, full_path)
    % pause()
    close(f)

end

%% Visualization of images in the same figures - first 6 harmonics
addpath(genpath(output_cortical))
addpath(genpath(output_subcortical))

% Dimensions of the figure
nRows = noh + 1;
nCols = 2;

f = figure;
t = tiledlayout(nRows, nCols, 'TileSpacing', 'none', 'Padding', 'compact');
% t.TileSpacing = 'none';  

% Visualize all images together
for h = 1:noh % rows

        % Cortical image
        cortical_name = fullfile(output_cortical,sprintf('ext_H%d_n%d.jpg', noh, h));
        cortical_tmp = imread(cortical_name);

        % Subcortical image
        subcortical_name = fullfile(output_subcortical,sprintf('ext_H%d_n%d.jpg', noh, h));
        subcortical_tmp = imread(subcortical_name);
    
        % --- Cortical image ---
        ax1 = nexttile((h - 1) * nCols + 1);
        imshow(cortical_tmp(400:600,100:1050,:));
    
        % Add text to the left of the cortical image
        text(ax1, 0, 0.5, sprintf('\\psi_{%d}', h), ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'right', ...
            'FontSize', 18);
    
        % --- Subcortical image ---
        ax2 = nexttile((h - 1) * nCols + 2);
        imshow(subcortical_tmp(400:600,100:1050,:));
end

% --- Colorbar ---
ax_cb = nexttile([1 2]); 
pos = get(ax_cb, 'Position');
pos(3) = pos(3) * 1.5;
set(ax_cb, 'Position', pos);

imshow(cortical_tmp(630:780, 320:850 ,:));
axis off;

% Set full screen
set(f, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Save
file_name = sprintf('%dnn_cortical_H%d.jpg', nn,noh);
full_path = fullfile(output_dir, file_name);
saveas(f, full_path);

%% Visualize first 20 harmonics
addpath(genpath(output_cortical))
addpath(genpath(output_subcortical))

% Dimensions of the figure
noh = 20;
nRows = noh/2 +1;  
nCols = 4;   

f = figure;
t = tiledlayout(nRows, nCols, 'TileSpacing', 'none', 'Padding', 'compact');

% Loop through harmonics (2 for each line)
for h = 1:noh
    % Calculate row and column position
    if h <= 10
        % Columns 1–2
        row = h;
        col_offset = 0;
    else
        % Columns 3–4
        row = h - 10;
        col_offset = 2;
    end

    % Load images
    cortical_name = fullfile(output_cortical, sprintf('ext_H%d_n%d.jpg', noh, h));
    subcortical_name = fullfile(output_subcortical, sprintf('ext_H%d_n%d.jpg', noh, h));
    cortical_tmp = imread(cortical_name);
    subcortical_tmp = imread(subcortical_name);

    % Index initial tile
    base_tile = (row - 1) * nCols + 1 + col_offset;

    % --- Cortical ---
    ax1 = nexttile(base_tile);
    imshow(cortical_tmp(400:600, 100:1050, :));
    text(ax1, 0, 0.5, sprintf('\\psi_{%d}', h), ...
        'Units', 'normalized', ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'right', ...
        'FontSize', 18);

    % --- Subcortical ---
    ax2 = nexttile(base_tile + 1);
    imshow(subcortical_tmp(400:600, 100:1050, :));
end

% --- Colorbar in all the last line (columns 1-4) ---
ax_cb = nexttile([1 4]);
pos = get(ax_cb, 'Position');
pos(3) = pos(3) * 1.02;  % allarga leggermente se vuoi
set(ax_cb, 'Position', pos);
imshow(cortical_tmp(630:780, 320:850 ,:));  % Usa uno dei file già caricati
axis off;

% Full screen
set(f, 'Units', 'normalized', 'Position', [0 0 1 1]);
