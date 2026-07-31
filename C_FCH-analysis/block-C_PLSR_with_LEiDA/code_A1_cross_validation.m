%% PLSR model to explore relation between LEiDA and FCH
% Part A1, the script continues with code_A2_perm_test.m

% Calculates cross-validated MSE and find the number of PLS components that 
% minimizes it.

% Jetro J. Tuulari, jetro.tuulari@utu.fi

% Adapted by:
% Isabella L.C. Mariani Wigley, 09 / 2025; ilmawi@utu.fi
% Aurora Berto, 09 / 2025; aurber@utu.fi

clear; close all; clc

%% Set parameters for reproducibility

s = rng;
rng(1);

%% Load data
% Predictors and responses

% Load FCH results
fch_path = '/Users/auroraberto/Desktop/25-26/UTU/ABCD_main-article/data'; % USER.adapt!
load(fullfile(fch_path, "Harmonics.mat"))

% Load megaLEiDA results
leida_path = '/Users/auroraberto/Desktop/25-26/UTU/ABCD_main-article/data'; % USER.adapt!
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

% Define the number of folds for cross-validation
numFolds = 10;
cvp = cvpartition(size(H,1), 'KFold', numFolds); % fixed partition

% Define the maximum number of PLS components to consider
maxNumComponents = 10;

% Preallocate a vector to hold the cross-validated MSE for each number of components
cvMSE = zeros(maxNumComponents, 1);
cvMSE_rand = zeros(maxNumComponents, 1);

% Preallocate a vector to store the optimal number of components selected
nComp_opt = zeros(6, 1);
figure('Color','w','Units','normalized','Position',[0.05 0.05 0.6 0.9])

for i_state = 1:6
    % Extract the predictors and responses for this state
    Y = V(:, i_state); % responses

    for i_comp = 1:maxNumComponents
        [~, ~, ~, ~, ~, ~, MSE_rand] = plsregress(H_rand, Y, i_comp, 'CV', cvp);
        cvMSE_rand(i_comp) = MSE_rand(2, i_comp+1); % take the mean of the second row of MSE
    end

    for i_comp = 1:maxNumComponents
        [~, ~, ~, ~, ~, ~, MSE] = plsregress(H, Y, i_comp, 'CV', cvp);
        cvMSE(i_comp) = MSE(2, i_comp+1); % take the mean of the second row of MSE
    end
    
% Plot the cross-validated MSE
    subplot(3,2,i_state)
    ax = gca;
    plot(1:3, cvMSE(1:3), '-o', 1:3, cvMSE_rand(1:3), '-r^','LineWidth',1.5);
% plot(1:maxNumComponents, cvMSE_rand, 'r^')
    xlabel('Number of PLS components','FontSize',13);
    ylabel('Cross-validated MSE','FontSize',13);
    ax = gca;
    ax.YAxis.Exponent = 0;
    ytickformat('%.3f')
    title(['PLS components vs. CV MSE for ' stateNames{i_state}],'FontSize',15);

% Find the number of components that gives the minimum MSE
    [~, optimalNumComponents] = min(cvMSE);
    nComp_opt(i_state,1) = optimalNumComponents;
    disp(['The optimal number of PLS components for state ' num2str(i_state) ' is ' num2str(optimalNumComponents) ' (MSE =' num2str(min(cvMSE)) ')']);
end

% Save the figure
saveas(gcf, 'cross_validation.png');
% Close the figure (optional, to save memory)
% close(gcf);
