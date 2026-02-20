%% MegaLEiDA reordering
% Reorder megaLEiDA states with decreasing Fractional Occupancy for each K

% Isabella L.C. Mariani Wigley, 09 / 2025, ilmawi@utu.fi
% Aurora Berto, 09 / 2025; aurber@utu.fi

clear; close all; clc

%% Input path and load data
input_path = '/path/to/pooled_LEiDA/Results'; % USER.adapt!
input_file = fullfile(input_path,'megaLEiDA_results_k20.mat'); % USER.adapt!

load(input_file, "Kmeans_results", "N_areas", "n_Condition_1", "rangeK", "P", "LT", "PT")

%% Output path
output_path = '/path/to/output/pooled_LEiDA'; % USER.adapt!
output_folder = fullfile(output_path, "Results"); % USER.adapt!

if ~exist(output_folder,"dir")
    mkdir(output_folder)
end

%% Set structures to store data

% Structure to store brain states for all clusters
K_names = strcat('K',string(rangeK));
C_ordered = cell2struct(cell(1, numel(rangeK)), K_names, 2);

% Matrices to store probabilities and lifetimes
P_ordered = zeros(size(P));
LT_ordered = zeros(size(LT));
PT_ordered = cell(size(PT));

%% Reorder brain states based on correlations
% For each cluster from K = 2
for k = rangeK(1):rangeK(end)

    % Set current K
    K_tmp = sprintf("K%d",k);

    % Select probability matrix for current K
    Pk = P(:,k-1,1:k);
    LTk = LT(:,k-1,1:k);

    % Calculate mean probability for each state
    Pk_mean = mean(Pk,1);

    % Reorder probabilities in decreasing order
    [Pk_val, idx] = sort(Pk_mean, 'descend');

    % Apply reordering to centroids, FO and DT
    C_tmp = Kmeans_results{1,k-1}.C;
    C_ordered.(K_tmp) = C_tmp(idx,:);

    P_ordered(:,k-1,1:k) = Pk(:,:,idx);
    LT_ordered(:,k-1,1:k) = LTk(:,:,idx);

    % Apply reordering to transition probabilities
    % For each subject
    for sub = 1:size(PT,1)

        % Select the current transition matrix
        PT_tmp = PT{sub,k-1};

        % Apply reordering
        PT_order_tmp = PT_tmp(idx, idx);

        % Save it in the matrix
        PT_ordered{sub,k-1} = PT_order_tmp;
    end
end

%% Save results
save(fullfile(output_folder,"Kmeans_results_ordered_k20.mat"),"C_ordered","P_ordered","LT_ordered","PT_ordered")
