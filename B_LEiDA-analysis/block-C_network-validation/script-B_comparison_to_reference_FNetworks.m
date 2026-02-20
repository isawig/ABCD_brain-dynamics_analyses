%% Comparison to reference intrinsic functional networks
% Articles: Cabral et al. 2020, Yeo et al. 2011

% Isabella L.C. Mariani Wigley, 04 / 2025, ilmawi@utu.fi
% Aurora Berto, 04 / 2025; aurber@utu.fi

%% Rationale behind the code
% Intrinsic Functional Networks, typically assessed using correlation 
% analysis, have been consistently detected in large cohorts of 
% resting-state fMRI experiments (Yeo et al., 2011).

% Here, we verify if the centroids obtained from clustering
% BOLD phase leading eigenvectors obtained at TR resolution share
% spatial similarities with the seven cerebral intrinsic functional
% networks estimated by Yeo et al. (2011) clustering correlation-
% based functional connectivity. 

% 1. Take the mask of the Yeo parcellation into seven non-overlapping 
% functional networks defined in MNI152 space

% 2. Take the mask of the Schaefer100 parcellation in the same MNI152 space

% 3. Calculate, for each of the 100 cortical brain areas, the
% proportion of voxels assigned to each of the seven functional
% networks (output: 7 1x100 vectors representing the intrinsic functional 
% networks in Schaefer space)

% 4. Compute the Pearson’s correlation between these seven networks and 
% the centroids obtained from our clustering analysis across
% the whole range of k explored (setting all negative values of the
% centroids’ vectors to zero, to consider only the areas whose BOLD
% phase is shifted from the main orientation)

%% Calculate ROIs assigned to each functional network

clear; close all; clc

%% Schaefer 100

% Load Schaefer100 file
schaefer100_path = '/path/to/comparison_networks_Yeo_et_al/Atlas/Schaefer100'; % USER.adapt!
schaefer100_name = 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv'; % USER.adapt!

schaefer100_file = readtable(fullfile(schaefer100_path, schaefer100_name));

% Extract ROIs names
schaefer100_ROIs = table2array(schaefer100_file(1:100,"ROIName"));

% Define the 7 network vectors
vis_net_cort = contains(schaefer100_ROIs, 'Vis');   % visual
SomMot_net_cort = contains(schaefer100_ROIs, 'SomMot'); % somatomotor
DorsAttn_net_cort = contains(schaefer100_ROIs, 'DorsAttn'); % dorsal attention
SalVentAttn_net_cort = contains(schaefer100_ROIs, 'SalVentAttn'); % salience ventral attention
Limbic_net_cort = contains(schaefer100_ROIs, 'Limbic'); % limbic
Cont_net_cort = contains(schaefer100_ROIs, 'Cont'); % frontoparietal control
Default_net_cort = contains(schaefer100_ROIs, 'Default'); % default mode

% Load Schaefer100 voxel count for all networks
voxel_count = readtable(fullfile(schaefer100_path,"NetworkVoxels_Yeo7_assignment.csv"));

%% Define networks

vis_net = vis_net_cort;
SomMot_net = SomMot_net_cort;
DorsAttn_net = DorsAttn_net_cort;
SalVentAttn_net = SalVentAttn_net_cort;
Limbic_net = Limbic_net_cort;
Cont_net = Cont_net_cort;
Default_net = Default_net_cort;

%% Concatenate in the same matrix

n_net = 7;
networks = [vis_net, SomMot_net, DorsAttn_net, SalVentAttn_net, Limbic_net, Cont_net, Default_net];

%% Load megaLEiDA results

megaLEiDA_path = '/path/to/pooled_LEiDA/results'; % USER.adapt!
load(fullfile(megaLEiDA_path, "Kmeans_results_ordered_k20.mat")); % USER.adapt!

% LEiDA_path = '/Users/auroraberto/Desktop/24-25/MSc_thesis/5_LEiDA_results/Results';
% load(fullfile(LEiDA_path,"reordering_results.mat"))
% % Choose a site
% brain_states = sites_centroids_ordered.S014;

%% Set parameters

rangeK = 2:20;
K_names = strcat('K',string(rangeK));

%% Set to zero negative centroids

for k = 1:length(rangeK)
    
    % Current K
    Kname = K_names(k);
    centroids_tmp = C_ordered.(Kname);

    % Set to zero negative values
    centroids_tmp(centroids_tmp < 0) = 0;

    % Save
    brain_states_pos.(Kname) = centroids_tmp;
end

%% Calculate correlations between networks and centroids

for k = 1:length(rangeK)

    % Current K
    Ktmp = k+1;
    Kname = K_names(k);

    % Initialize matrices
    correlations = zeros(n_net, Ktmp);
    p_values = zeros(n_net, Ktmp);

    % Select current K cortical brain states
    brain_states_tmp = brain_states_pos.(Kname)(:,1:100);

    for c = 1:Ktmp % for each state
        for j = 1:n_net

            % Calculate correlation and relative p-value
            [correlations(j,c), p_values(j,c)] = corr(brain_states_tmp(c,:)', networks(:,j));
        end
    end

    % Calculate maximum correlation and relative network
    [dominantCorr, dominantNetwork] = max(correlations, [], 1);

    % Save in the structure
    networks_correlations_k20.correlations.(Kname) = correlations;
    networks_correlations_k20.p_values.(Kname) = p_values;
    networks_correlations_k20.dominantCorr.(Kname) = round(dominantCorr,3);
    networks_correlations_k20.dominantNetwork.(Kname) = dominantNetwork;
end


%% Plot heatmaps of network-centroid correlations

figure;
for k = 1:length(rangeK)
    
    % Current K
    Kname = K_names(k);
    correlations = networks_correlations_k20.correlations.(Kname);
    
    subplot(5,4,k)
    imagesc(correlations)
    colormap(parula)
    colorbar
    caxis([-1 1]) % range of Pearson's r
    
    xlabel('Centroid (State)')
    ylabel('Yeo Network')
    
    title(['K = ', num2str(rangeK(k))])
    set(gca, 'YTick', 1:n_net, ...
             'YTickLabel', {'Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default'})
end

set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);
sgtitle('Correlations between centroids and Yeo et al. 7 networks')
saveas(gcf,'Heatmap_correlation_centroids_and_Yeo','jpg')

%% Save
save('networks_correlation_with_Yeo',"networks_correlations_k20") % USER.adapt!

