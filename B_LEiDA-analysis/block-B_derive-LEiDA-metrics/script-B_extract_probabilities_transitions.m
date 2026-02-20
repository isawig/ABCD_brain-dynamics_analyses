%% Script for extracting the probability values from LEiDA analyses

% Jetro J. Tuulari, 02/2023, jjtuul@utu.fi
% Isabella L.C. Mariani Wigley, 04 / 2025, ilmawi@utu.fi
% Aurora Berto, 04 / 2025, aurber@utu.fi

%% Load LEiDA results for further analyses
% you need this in case you have not used osl before
clear; close all; clc

% Input path
addpath(genpath('/path/to/LEiDA/results'))
addpath(genpath('/path/to/osl-core-master'))

% then you start osl
% osl_startup

% site
% site = "S090";

% you should have run the script LEiDA_Markov_Chain.m before this
% load(sprintf('%s_LEiDA_results.mat',site),'n_Condition_1','rangeK')
% load('sitesLEiDA_outcomes.mat')

load('megaLEiDA_results_k20.mat','n_Condition_1','rangeK')
load('Kmeans_results_ordered_k20.mat','PT_ordered')

% Select transition matrices for the current site
% PT_order = PT_order_all_sites.(site);

close all
Index_Condition_1=1:n_Condition_1;  % The first data are for Condition_1...
n_Subjects=n_Condition_1; % +n_Condition_2; % Total number of scans index

%% Output path

% output_path = sprintf('/path/to/LEiDA/site_results/Transitions/%s',site);
output_path = '/path/to/LEiDA/pooled_results/Transitions/K2toK20';
addpath(genpath(output_path))

% Create output folders to store results
if ~exist(output_path,'dir')
    mkdir(output_path)
end


%% Plot the results for a defined K

% disp(['Choose number of clusters between ' num2str(rangeK(1)) ' and ' num2str(rangeK(end)) ])
% K = input('Number of clusters: ');
% k=find(rangeK==K);
% print the K that was chosen
% disp(['K=' num2str(K)])

% Define K
% K = 6;

% For each K
for K = rangeK(1):rangeK(end)
    
    % Define output path
    K_output_path = fullfile(output_path,sprintf("K%d",K));
    addpath(genpath(K_output_path))

    % Create output folders to store results
    if ~exist(K_output_path,'dir')
        mkdir(K_output_path)
    end
    
    % initialise the matrix
    PT_K=zeros(n_Subjects,K,K);
    
    for s=1:n_Subjects
        PT_K(s,:,:)=PT_ordered{s,rangeK==K};
    end
    
    PT_K_Condition_1=squeeze(mean(PT_K(Index_Condition_1,:,:)));
    % PT_K_Condition_2=squeeze(mean(PT_K(Index_Condition_2,:,:)));
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    colormap(osl_colormap('rwb'))
    colorbar
    
    % Transition Matrix in Group / Condition 1
    % subplot(1,2,1)
    imagesc(PT_K_Condition_1,[0 0.04])
    for j1 = 1:K
        for j2 = 1:K
            sem_cont=std(PT_K(Index_Condition_1,j1,j2))/sqrt(numel(Index_Condition_1));
            caption = {['  ' sprintf('%.3f', PT_K_Condition_1(j1,j2))],[char(177) num2str(sem_cont,'%.3f')]};
            if j1==j2
                text(j2-0.35, j1, caption, 'FontSize', 11, 'Color', 'k','FontWeight','bold');
            else
                text(j2-0.35, j1, caption, 'FontSize', 11, 'Color', 'w','FontWeight','bold');
            end
        end
    end
    xticks(1:K);
    yticks(1:K); 
    % yticklabels({'State 1','State 2','State 3','State 4','State 5'})
    %xtickangle(45);ytickangle(45)
    set(gca,'FontSize',12,'FontWeight','bold')
    axis square
    ylabel('From PL State')
    xlabel('To PL State')
    title({'Transition Probabilities'})
    
    % save the figure
    % fig_name = sprintf('%s_k%d_0_Transition_probability_matrix.png',site,K);
    fig_name = sprintf('k%d_0_Transition_probability_matrix.png',K);
    fig_path = fullfile(output_path,fig_name);
    saveas(f, fig_path)
    close(f)

    % Save the derived measures 
    for k=1:K
        x = squeeze(PT_K(:,k,:)); % dependent variable(s) transition probabilities % b=squeeze(PT_K(:,6,:)); % gives a matrix K = 6; 6->1, 6->2, 6->3 etc.
        %x = x(:, 1:K); % use all rows, and selected amount of cluster solutions 
        e = array2table(x);
        % table_name = sprintf('%s_K%dx%d_Extracted_LEiDA_transitions.csv',site,K,k);
        table_name = sprintf('K%dx%d_Extracted_LEiDA_transitions.csv',K,k);
        table_path = fullfile(K_output_path,table_name);
        writetable(e, table_path);
    end
end
%% END
