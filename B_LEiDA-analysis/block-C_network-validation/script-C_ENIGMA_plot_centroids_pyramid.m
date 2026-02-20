%% ENIGMA script to plot megaLEiDA networks using Yeo et al. classification
% already ordered

% Isabella L.C. Mariani Wigley, 04 / 2025, ilmawi@utu.fi
% Aurora Berto, 04 / 2025; aurber@utu.fi

clear; close all; clc

%% Load data and set path
addpath(genpath('/path/to/ENIGMA/matlab')) % USER.adapt!

input_path = '/path/to/pooled_LEiDA/results'; % USER.adapt!
load(fullfile(input_path,"Kmeans_results_ordered_k20.mat")) % USER.adapt!


%% Output path
% Output folders: where figures will be saved
output_dir = sprintf('/path/to/comparison_networks_Yeo_et_al/ENIGMA'); % USER.adapt!
output_cortical = fullfile(output_dir, '/Cortical/K2toK20'); % USER.adapt!

% Create output folders to store figures
if ~exist(output_cortical,'dir')
    mkdir(output_cortical)
end


%% Set parameters

K_values = 2:20;
K_names = strcat('K',string(K_values));

N_cortical_areas = 100;

%% Load network states classification

classification_path = '/path/to/comparison_networks_Yeo_et_al/Results-Yeo'; % USER.adapt!
load(fullfile(classification_path,"networks_correlation_with_Yeo.mat")) % USER.adapt!

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


%% 1. Plot brain states

%% Save cortical figures
% for each cluster
for k = K_values(1):K_values(end)

    % Set the cluster name
    K_tmp = K_names(k-1);

    % Select current network correlation and p-values with Yeo et al.
    correlation_tmp = networks_correlations_k20.dominantCorr.(K_tmp);
    network_tmp = networks_correlations_k20.dominantNetwork.(K_tmp);
    pvalues_tmp = networks_correlations_k20.p_values.(K_tmp);

    % for each state c
    for c = 1:k
    
        % Cluster indices 
        V = C_ordered.(K_tmp)(c,1:100);

        % Map parcelled data to the cortical surface
        Schaefer100 = parcel_to_surface(V,'schaefer_100_fsa5');

        % Select current network's color
        net_idx = network_tmp(c);
        net_name = network_names{net_idx};
        net_color = color_map(net_name);

        % Decide whether to use the network color or black
        if correlation_tmp(c) < 0.4
            net_color = [0 0 0]; % black for low correlation
        end
        
        % Project the results on the surface brain
        f = figure('Visible', 'off'); % not displayed
        [a, cb] = plot_cortical(Schaefer100, 'surface_name', 'fsa5','color_range',[-1 1]);
        colorbar off;
        
        % Delete left and inner cortical visualizations and colorbar
        delete(a([1, 2, 3]));
        delete(cb);

        % Set the right hemisphere in the center and zoom it
        set(a(4), 'Position', [0.2, 0.1, 0.6, 0.8]);

        % Access to the surfice patch
        p = a(4).Children(end);
        
        % Background color: gray
        cmap = repmat([0.8 0.8 0.8], length(p.FaceVertexCData), 1);
    
        % Get the original values
        vals = p.FaceVertexCData;
    
        % Find positive indices (regions to color)
        pos_idx = vals > 0;
        
        % Normalize the positive values between 0 and 1
        norm_vals = vals;
        norm_vals(~pos_idx) = 0;
        max_val = max(norm_vals);
        if max_val > 0
            norm_vals = norm_vals / max_val;
        end

        % Apply gradient from white to net_color
        for v = find(pos_idx)'
            t = norm_vals(v); % 0 (white) to 1 (net_color)
            cmap(v,:) = (1 - t) * [1.2 1.2 1.2] + t * net_color;
        end
       
        % Update surface colors
        p.FaceVertexCData = cmap;

        % Save the figure
        file_name = sprintf('shaded_k%d_c%d.jpg',k,c);
        full_path = fullfile(output_cortical, file_name);
        saveas(f, full_path)
        % pause()
        close(f)
    end
end

%% Visualization of cortical images in the same figures
addpath(genpath(output_cortical))

% number of clusters
K_clust = 2:20;

% Dimensions of the figure
nRows = K_clust(end) - K_clust(1) + 1;
nCols = K_clust(end);

f = figure;
t = tiledlayout(nRows, nCols, ...
    'TileSpacing', 'none', ... % puoi anche provare 'loose'
    'Padding', 'compact');

% Aumenta il rapporto verticale (più alta e stretta)
set(f, 'Units', 'normalized', 'Position', [0.1 0 0.7 1]); 

% Visualize all images together
for k = K_clust(1):K_clust(end) % rows

    % Set the cluster name
    K_tmp = K_names(k-1);

    % Select current network correlation and p-values with Yeo et al.
    correlation_tmp = networks_correlations_k20.dominantCorr.(K_tmp);
    network_tmp = networks_correlations_k20.dominantNetwork.(K_tmp);
    pvalues_tmp = networks_correlations_k20.p_values.(K_tmp);


    for c = 1:k % columns visible for this row

        % File name
        fig_name = sprintf('shaded_k%d_c%d.jpg', k, c);
        fig_tmp = imread(fig_name);

        % Find current correlation and p-value
        % Correlation
        network_corr = correlation_tmp(c);
        
        % P-value
        network = network_tmp(c);
        network_pval = pvalues_tmp(network,c);

        % Find exponent to visualize in the title
        exponent = floor(log10(abs(network_pval)));
        display_exponent = exponent + 1;

        % Tile's position
        p = (k - K_clust(1)) * K_clust(end) + c;
        ax = nexttile(p);
        imshow(fig_tmp(180:end, 200:950, :))
        if network_corr >= 0.4
            % title(sprintf('r=%.2f p<1e%d', network_corr, exponent), ...
            %     'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
            title(sprintf('r = %.2f', network_corr), ...
                'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end

        % If the first column, add the string "k = ..."
        if c == 1
            text(ax, -0.3, 0.5, sprintf('(%d)', k), ...
                'Units', 'normalized', ...
                'VerticalAlignment', 'middle', ...
                'HorizontalAlignment', 'right', ...
                'FontSize', 10);
        end
    end
end
% Title
% sgtitle('Overlap of cluster centroids with functional Brain Networks');

% Save
file_name = sprintf('cortical_k%d_to_%d.jpg', K_clust(1), K_clust(end));
full_path = fullfile(output_dir, file_name);
saveas(f, full_path);

