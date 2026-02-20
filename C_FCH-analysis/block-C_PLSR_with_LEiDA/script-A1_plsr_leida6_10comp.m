%% PLSR model to explore relation between LEiDA and FCH
% Part A1, the script continues with code_A2_cross_validation.m

% Calculates percentage of LEiDA variance explained varying the number of PLS 
% components estimated from all harmonics (ordered and randomized) and plots the 
% predicted vs fitted response.

% Jetro J. Tuulari, jetro.tuulari@utu.fi

% Adapted by:
% Isabella L.C. Mariani Wigley, 09 / 2025; ilmawi@utu.fi
% Aurora Berto, 09 / 2025; aurber@utu.fi

clear; close all; clc

%% Set parameters for reproducibility

s = rng;
rng(s);

%% Load data
% Predictors and responses

% Load FCH results
fch_path = '/path/to/FCH/FCH_results/FCH_nn10_H113/ALL'; % USER.adapt!
load(fullfile(fch_path, "Harmonics.mat"))

% Load megaLEiDA results
leida_path = '/path/to/pooledLEiDA/Results'; % USER.adapt!
load(fullfile(leida_path, "Kmeans_results_ordered_k20.mat"))

%% Redefine variables' names

H = Y(:,2:end); % harmonics eigevectors, discard the null one
H_values = V(2:end); % harmonics eigenvalues
clear Y V

% Select LEiDA centroids with K=6
V = C_ordered.K6'; % LEiDA centroids

%% Prepare H_rand

%threshold = 0
%%V = V>threshold

H_rand = reshape(H(randperm(numel(H))), size(H));
[n,p] = size(H_rand);

%% Define parameters

% Define the state names
stateNames = {'State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6'};

%% Calculate percentage of explained variation and fitted models 

% Loop over each state
for i_state = 1:6

    % Extract the predictors and responses for this state
    %predictors H and H_rand
    Y = V(:, i_state); % responses
    
    % Perform PLSR with the optimal number of components randomized
    [XL_r, YL_r, XS_r, YS_r, BETA_r, PCTVAR_r] = plsregress(H_rand, Y, 10);
    [XL, YL, XS, YS, BETA, PCTVAR] = plsregress(H, Y, 10);

    % Calculate the predicted responses using only the first optimalNumComponents columns of XS
    yfitPLS_r = [ones(n,1) H_rand]*BETA_r;
    yfitPLS = [ones(n,1) H]*BETA;
    
    % Plot the predicted responses vs fitted response
    fig = figure;
    plot(Y, yfitPLS, 'bo', Y, yfitPLS_r, 'r^','LineWidth',1.5);
    xlabel('Observed Response','FontSize',16);
    ylabel('Fitted Response','FontSize',16);
    title(['Observed vs Fitted Response for ' stateNames{i_state}],'FontSize',18);
    legend({'PLSR with harmonics' 'PLSR with randomized harmonics'},  ...
	'location','NW','FontSize',14);
    saveas(fig, sprintf(['10comp_obv_pred_' stateNames{i_state} '.png'], i_state));
    % Close the figure (optional, to save memory)
    close(fig);

    
    % Plot the percent variance explained
    fig = figure;
    plot(1:10,cumsum(100*PCTVAR(2,:)),'-bo', 1:10, cumsum(100*PCTVAR_r(2,:)),'-r^','LineWidth',1.5);
    xlabel('Number of PLS components','FontSize',16);
    ylabel('Percent Variance Explained in Y','FontSize',16);
    ylim([0 110])
    title(['Percent Variance Explained by Each PLS Component for ' stateNames{i_state}],'FontSize',17);
    legend({'PLSR with harmonics' 'PLSR with randomized harmonics'},'location','SE','FontSize',14)
    saveas(fig, sprintf(['10comp_per_var_' stateNames{i_state} '.png'], i_state));    % Close the figure (optional, to save memory)
    close(fig);

end
