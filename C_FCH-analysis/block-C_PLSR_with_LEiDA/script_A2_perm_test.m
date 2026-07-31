%% Permutation test
% Part A2, to run after code code_A1_cross_validation.m

% Runs a permutation test using the optimal number of PLS components 
% (selected in Script A1 based on the minimum MSE in the observed data) 
% to assess whether the MSE of the observed FCHs is lower than expected 
% under a null distribution generated from randomized harmonics.

% Isabella L.C. Mariani Wigley, 09 / 2025; ilmawi@utu.fi
% Aurora Berto, 09 / 2025; aurber@utu.fi

close all;
s = rng;
rng(1);

addpath(genpath('/Users/auroraberto/Desktop/25-26/UTU/ABCD_main-article/GitHub_codes/fdr_bh'))

numPerm = 1000;
numFolds = 10;
cvp = cvpartition(size(H,1), 'KFold', numFolds); % fixed partition

% Number of components selected from CV
optimalComponents = nComp_opt; % [1,3,1,2,1,3]

observedMSE = zeros(6,1);
permMSE = zeros(numPerm,6);

for i_state = 1:6
    fprintf('\nState %d\n',i_state);
    Y = V(:,i_state);                   % LEiDA centroid
    nComp = optimalComponents(i_state); % optimal n. components

    % Observed model (LEiDA from ordered FCHs)
    [~,~,~,~,~,~,MSE] = plsregress(H,Y,nComp,'CV',cvp);
    observedMSE(i_state) = MSE(2,nComp+1); 

    % (The first row of MSE contains mean
    % squared errors for the predictor variables in X and the second row
    % contains mean squared errors for the response variable(s) in Y.)

    % Permutations (LEiDA from randomized FCHs)
    for p = 1:numPerm
        H_perm = reshape(H(randperm(numel(H))),size(H)); % Randomize harmonics 
        [~,~,~,~,~,~,MSEperm] = plsregress(H_perm,Y,nComp,'CV',cvp);
        permMSE(p,i_state) = MSEperm(2,nComp+1); 
    end
end

%% Empirical p-values
% Smaller MSE = better model

% Sum the number of permutations (from 1 to 1000) with p-value lower than
% the one from the observed model (sum(permMSE(:,i_state) <= observedMSE(i_state))
% Then estimate the empirical p-value with the formula (count + 1)/(numPerm + 1)

p_empirical = zeros(6,1);

for i_state = 1:6
    count = sum(permMSE(:,i_state) <= observedMSE(i_state));
    p_empirical(i_state) = (count + 1)/(numPerm + 1);
end

%% FDR correction
[p_fdr,~,~,adj_p] = fdr_bh(p_empirical);

%% Results table
Results = table( ...
    string(stateNames)', ...
    optimalComponents, ...
    observedMSE, ...
    p_empirical, ...
    adj_p, ...
    p_fdr, ...
    'VariableNames',{'State','Ncomponents','ObservedMSE',...
                     'pPermutation','qFDR','Significant'});
disp(Results)
save('Permutation_results.mat',...
    'Results','permMSE','observedMSE','p_empirical','adj_p');

%% Plot null distributions

NETWORK = {'DMN','VIS','VIS','VIS','DMN','DAN'};
figure('Color','w','Units','normalized','Position',[0.05 0.05 0.6 0.9])

for i_state = 1:6
    subplot(3,2,i_state)
    histogram(permMSE(:,i_state),30,"FaceColor",[0.7 0.7 0.7])
    hold on
    xline(observedMSE(i_state),'--r','LineWidth',2);
    xlabel('Cross-validated MSE','FontSize',13)
    ylabel('Count','FontSize',13)
    ylim([0 110]); xlim([0.001 0.0075])
    if adj_p(i_state) < 0.001
        p_plot = 0.001;
    elseif adj_p(i_state) < 0.01
        p_plot = 0.01;
    elseif adj_p(i_state) < 0.05
        p_plot = 0.05;
    else 
        p_plot = adj_p(i_state);
    end
    title(sprintf('%s (%s) - p < %.3f',stateNames{i_state}, NETWORK{i_state}, p_plot),"FontSize",15)
    legend('Permutations','Observed')
end
exportgraphics(gcf,'nullDistributions.png')