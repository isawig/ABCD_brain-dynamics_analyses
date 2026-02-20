%% Script to compute LEiDA Markov Chain transition probabilities

% Jetro J. Tuulari, 01 / 2020, jetro.tuulari@utu.fi
% Isabella L.C. Mariani Wigley, 04 / 2025, ilmawi@utu.fi
% Aurora Berto, 04 / 2025, aurber@utu.fi

function script_A_LEiDA_Markov_Chain

%% Add path and load data
close all; clc

% Set path
addpath /path/to/LEiDA/results

% Load data
load megaLEiDA_results_k20.mat Kmeans_results rangeK P Time_all n_Condition_1 

%% For every fMRI scan and for each K calculate probability of transitions PT.

PT = cell(n_Subjects,size(rangeK,2));

for K=1:length(rangeK)
    for s=1:n_Subjects
            % Select the time points representing this subject
            T = (Time_all==s);
            Ctime=Kmeans_results{K}.IDX(T);

            % Probability Transition Matrix
            transferMatrix = zeros(rangeK(K));
            for tp = 2:size(Ctime,2)%-5
                transferMatrix(Ctime(tp-1),Ctime(tp)) = transferMatrix(Ctime(tp-1),Ctime(tp)) + 1;
            end
            PT{s,K} = transferMatrix./(size(Ctime,2)-1); % normalised by T-1
            PT{s,K} = PT{s,K}./squeeze(P(s,K,1:rangeK(K)));
            PT{s,K}(isnan(PT{s,K}))=0;
    end    
end

%% Save results
save megaLEiDA_results_k20.mat PT -append % removed from the saved list as there are no stats; PT_pval
