%% Plot reference networks (Yeo et al.)
% Create a figure with seven reference RSN and overlap to the pyramid plot
% with LEiDA centroids visualization, color-coded according to Yeo matching
% network.

% Isabella L.C. Mariani Wigley, 04 / 2025, ilmawi@utu.fi
% Aurora Berto, 04 / 2025; aurber@utu.fi

clear; close all; clc

%% Set path
addpath(genpath('/path/to/ENIGMA/matlab')) % USER.adapt!

output_path = '/path/to/comparison_networks_Yeo_et_al/ENIGMA/Reference_networks'; % USER.adapt!

if ~exist(output_path,"dir")
    mkdir(output_path)
end

%% Assign each ROI to the correspondent network

% Load Schaefer100 file
schaefer100_path = '/path/to/comparison_networks_Yeo_et_al/Atlas/Schaefer100'; % USER.adapt!
schaefer100_name = 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv'; % USER.adapt!

schaefer100_file = readtable(fullfile(schaefer100_path, schaefer100_name));

% Extract ROIs names
schaefer100_ROIs = table2array(schaefer100_file(1:100,"ROIName"));

% Define the 7 network vectors
network_labels = zeros(size(schaefer100_ROIs)); % inizializza il vettore

network_labels(contains(schaefer100_ROIs, 'Vis')) = 1;
network_labels(contains(schaefer100_ROIs, 'SomMot')) = 2;
network_labels(contains(schaefer100_ROIs, 'DorsAttn')) = 3;
network_labels(contains(schaefer100_ROIs, 'SalVentAttn')) = 4;
network_labels(contains(schaefer100_ROIs, 'Limbic')) = 5;
network_labels(contains(schaefer100_ROIs, 'Cont')) = 6;
network_labels(contains(schaefer100_ROIs, 'Default')) = 7;

network_ids = unique(network_labels);


%% Set different color for each network 
% Using Yeo et al. as reference

colors = {
    'visual',          [120,  18, 134]/255;  % violet
    'somatomotor',     [ 70, 130, 180]/255;  % blue
    'dorsal_attention',[  0, 118,  14]/255;  % dark green
    'ventral_attention',[196, 58, 250]/255; % light violet
    'limbic',          [220, 248, 164]/255;  % lime
    'frontoparietal',  [230, 148,  34]/255;  % yellow
    'default',         [205,  62,  78]/255;  % red
};
rgb_values = colors(:,2);

% Create a map to access to colors through indices
network_names = {'visual', 'somatomotor', 'dorsal_attention', 'ventral_attention', 'limbic', 'frontoparietal', 'default'};
color_map = containers.Map(network_names, rgb_values, 'UniformValues', false);

%% Save cortical figures
% for each cluster
for n = 1:length(network_ids)

    % Set the current network
    net_id = network_ids(n);
    net_name = network_names{n};
    
    % Network indices 
    V = double(network_labels == net_id);
   
    % Map parcelled data to the cortical surface
    Schaefer100 = parcel_to_surface(V,'schaefer_100_fsa5');

    % Select current network's color
    net_color = color_map(net_name);

    % Project the results on the surface brain
    f = figure('Visible', 'on'); % not displayed
    [a, cb] = plot_cortical(Schaefer100, 'surface_name', 'fsa5','color_range',[-1 1]);
    colorbar off;
    
    % Delete left and inner cortical visualizations and colorbar
    delete(a([1, 2, 3]));
    delete(cb);

    % Set the right hemisphere in the center and zoom it
    set(a(4), 'Position', [0.2, 0.1, 0.6, 0.8]);

    % Access to the surfice patch
    p = a(4).Children(end);
    
    % Creat a new colormap
    cmap = zeros(length(p.FaceVertexCData),3);
    for v = 1:length(p.FaceVertexCData)
        if p.FaceVertexCData(v) > 0 % if p = 1 -> network
            cmap(v,:) = net_color;        
        else                        % if p = 0 -> background
            cmap(v,:) = [0.8 0.8 0.8];  
        end
    end
    
    p.FaceVertexCData = cmap;

    % Save the figure
    file_name = sprintf('network_%d_%s.jpg',net_id,net_name);
    full_path = fullfile(output_path, file_name);
    saveas(f, full_path)
    % pause()
    close(f)
end


%% Visualization of cortical images in the same figures
addpath(genpath(output_path))

% Path to store image 
output_dir = '/path/to/comparison_networks_Yeo_et_al/ENIGMA'; % USER.adapt!

% Define titles
network_titles = {'Visual', 'Somatomotor', 'Dorsal Attention', 'Ventral Attention', 'Limbic', 'Frontoparietal', 'Default'};
 
% Dimensions of the figure
nRows = 1;
nCols = 7;

f = figure;
t = tiledlayout(nRows, nCols, 'TileSpacing', 'none', 'Padding', 'compact');

% Visualize all images together
for n = 1:length(network_ids)

    % Set the current network
    net_id = network_ids(n);
    net_name = network_names{n};

    % File name
    fig_name = sprintf('network_%d_%s.jpg',net_id,net_name);
    fig_tmp = imread(fig_name);

    % Tile's position
    p = n;
    ax = nexttile(p);
    imshow(fig_tmp)
    title(sprintf('%s',network_titles{n}),'FontSize',20)
end
% Set the figure as the whole screen
set(f, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Title
sgtitle('Reference functional brain networks (Yeo et al., 2011)','FontSize',30);

% Save
file_name = 'reference_networks.jpg'; % USER.adapt!
full_path = fullfile(output_dir, file_name);
saveas(f, full_path);


%% Add legend to megaLEiDA networks' image
img_path = '/path/to/comparison_networks_Yeo_et_al/ENIGMA'; % USER.adapt!

% Load images
main_img = imread(fullfile(img_path,'cortical_k2_to_20.jpg')); % USER.adapt!
legend_img = imread(fullfile(img_path,'reference_networks.jpg')); % USER.adapt!

% Create figure
figure('Units','normalized','OuterPosition',[0 0 1 1]);

% Add the main image as background
ax1 = axes('Position',[0 0 1 1]);
imshow(main_img, 'Parent', ax1, 'InitialMagnification', 'fit');
axis off;

% Add the legend image on the top-right of the screen
legend_width = 0.5;
legend_height = 0.5;
legend_left = 1 - legend_width - 0.01; 
legend_bottom = 1 - legend_height - 0.01;  

ax2 = axes('Position',[legend_left legend_bottom legend_width legend_height]);
imshow(legend_img(300:1000,:,:), 'Parent', ax2);
axis off;

% Bring main image on the back
uistack(ax1, 'bottom');

% Save
file_name = 'megaLEiDA_networks_complete.jpg'; % USER.adapt!
full_path = fullfile(output_dir, file_name);
saveas(gcf, full_path);

