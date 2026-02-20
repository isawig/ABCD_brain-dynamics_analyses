function nn = script_CSC_compute_leida_from_bold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
%
% This set of scripts processes, clusters and analyses BOLD data using LEiDA
% Comparing between Condition 1, Condition 2 and Condition 3
%
%
% A - Read the BOLD data from the folders and computes the BOLD phases
%   - Calculate the instantaneous BOLD synchronization matrix
%   - Compute the Leading Eigenvector at each frame from all fMRI scans
% B - Cluster the Leading Eigenvectors into recurrent Functional Networks
% C - Compute the probability and lifetimes of each FN in each session
%   - Saves the Eigenvectors, Clusters and statistics into LEiDA_results.mat
%
%% IMPORTANT NOTE: Visualization can be run independently once
%  LEiDA_results.mat has been saved:
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Joana Cabral March 2018
% joana.cabral@psych.ox.ac.uk
%
% First use in
% Cabral, et al. 2017 Scientific reports 7, no. 1 (2017): 5135.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adaptations by Jetro Tuulari, 07 / 2023 for local set up
% Adaptations by Jetro Tuulari, 08 / 2024 CSC's Puhti server
% jetro.tuulari@utu.fi

% Isabella L.C. Mariani Wigley, 04 / 2025, ilmawi@utu.fi
% Aurora Berto, 04 / 2025, aurber@utu.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User defined key settings are done here
% it is recommended to use a single source directory for all analyses
% place this and other script there :)

site = 'S090'; %QQQ

    % USER: Define here data parameters:
    n_Condition_1 = 210; %QQQ sample size
    N_areas    = 114; %QQQ THIS THE AAL ATLAS and cortical areas; adapt for the one you use
    TR         = 0.8; %QQQ 
    Tmax       = 380;  % Set here the total number of frames in each scan %QQQ
    % Note: if Tmax is different at each scan let me know to adapt.
    
    % USER: Change here to your folder where the timeseries data are
    % Directory  = '/scratch/project_2006897/bertoaur/5_megaLEiDA/All_subjs_data/'; %QQQ your working folder, input data
    Directory  = sprintf('/path/to/data/%s/temp/114_rows_380_columns/',site); %QQQ your working folder, input data

    % Filename_Condition_1=dir([Directory '*_sub-*.txt']); % fMRI timeseries data
    Filename_Condition_1=dir([Directory 'sub-*.txt']); % fMRI timeseries data

    % USER: Change here to your folder where data will be saved
    % Outdir = '/scratch/project_2006897/bertoaur/5_megaLEiDA/MegaLEiDA_results/';
    Outdir = sprintf('/path/to/data/%s/LEiDA_results/',site);

    % USER: We will concatenate data from all Sessions as it is loaded,
    % so please set here the indices that will correspond to Conditions and
    % update the indices later on as well, search / check QQQ comments..!
    Index_Condition_1=1:n_Condition_1;  % The first data are for Condition_1... 
    n_Sessions=n_Condition_1; % Total number of scan sessions index.. for easier adaptation if needed..


%% 1 - Compute the Leading Eigenvectors from the BOLD datasets

% if strcmp(mode,'run') % compute LEiDA

% % Bandpass filter settings %QQQ comment If no filtering will be done
    fnq=1/(2*TR);                 % Nyquist frequency
    flp = 0.02;                   % lowpass frequency of filter (Hz)
    fhi = 0.10;                   % highpass
    Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
    k=2;                          % 2nd order butterworth filter
    [bfilt,afilt]=butter(k,Wn);   % construct the filter
    clear fnq flp fhi Wn k

    disp('Processing the eigenvectors from BOLD data')    

    % Preallocate variables to save FC patterns and associated information
    V1_all   = zeros((Tmax-2)*n_Sessions,N_areas); % All leading eigenvectors
    t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*(Tmax-2))
    Time_all= zeros((Tmax-2)*n_Sessions,1); % Vector that links each frame to a subject

    for s=1:n_Sessions

        % USER: Adapt to load here the BOLD matrix (NxT) from each subject
        if sum(Index_Condition_1==s)
            BOLD = load([Directory Filename_Condition_1(s).name]);        
        end

% Get the BOLD phase using the Hilbert transform with filtering in place
        Phase_BOLD=zeros(N_areas,Tmax);
        for seed=1:N_areas
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
            signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end

                %QQQ apply this only when there is no filtering
                % Get the BOLD phase using the Hilbert transform 
                %         Phase_BOLD=zeros(N_areas,Tmax);
                %         for seed=1:N_areas
                %             BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
                %             Phase_BOLD(seed,:) = angle(hilbert(BOLD(seed,:)));
                %         end

        % Slide over time discarding the first and last epochs
        for t=2:Tmax-1

            %Calculate the Instantaneous FC (BOLD Phase Synchrony)
            iFC=zeros(N_areas);
            for n=1:N_areas
                for p=1:N_areas
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end

            % Get the leading eigenvector
            [V1,~]=eigs(iFC,1);
            if sum(V1)>0
                V1=-V1;
            end

            % Save V1 from all frames in all fMRI sessions
            t_all=t_all+1; % Update time
            V1_all(t_all,:)=V1;
            Time_all(t_all)=s;
        end
    end
    %%
    % save the estimated data
    % save LEiDA_extra.mat -v7.3

    save(fullfile(Outdir, 'LEiDA_vectors_full.mat'), '-v7.3');

    % then remove the variables we do not need in the immediate modeling
    clear BOLD tc_aal signal_filt iFC V1 Phase_BOLD

    % save a lean version of the estimated data for the next step
    save(fullfile(Outdir, 'LEiDA_vectors.mat'), '-v7.3');


  %% 2 - Cluster the Leading Eigenvectors

    % Load results
    load(fullfile(Outdir, 'LEiDA_vectors_full.mat'))
    
    disp('Clustering the eigenvectors into')
    % V1_all is a matrix containing all the eigenvectors:
    % Columns: N_areas are brain areas (variables)
    % Rows: (Tmax-2)*n_Sessions are all time points (independent observations)
    
    % USER: Set maximum/minimum number of clusters
    % There is no fixed number of FC states that the brain can display
    % Keep the range small for the first trials
    % Extend the range depending on the hypothesis of each work
    maxk=20; % QQQ
    mink=2; % QQQ
    rangeK=mink:maxk;
    
    % Set the parameters for Kmeans clustering
    Kmeans_results=cell(size(rangeK));
    
    for k=1:length(rangeK)
        disp(['- ' num2str(rangeK(k)) ' FC states'])
        % Distance can be cosine, cityblock, default one is sqeuclidean
        [IDX, C, SUMD, D]=kmeans(V1_all,rangeK(k),'Distance','cosine','Replicates',50,'MaxIter',5000,'Display','final','Options',statset('UseParallel',1)); % QQQ modify settings 
        [~, ind_sort]=sort(hist(IDX,1:rangeK(k)),'descend');
        [~,idx_sort]=sort(ind_sort,'ascend');
        Kmeans_results{k}.IDX=idx_sort(IDX);   % Cluster time course - numeric column vectors
        Kmeans_results{k}.C=C(ind_sort,:);       % Cluster centroids (FC patterns)
        Kmeans_results{k}.SUMD=SUMD(ind_sort); % Within-cluster sums of point-to-centroid distances
        Kmeans_results{k}.D=D(:,ind_sort);       % Distance from each point to every centroid 
    end

  %% 3 - Analyse the Clustering results
    
    % For every fMRI scan calculate probability and lifetimes of each state c.
    P=zeros(n_Sessions,maxk-mink+1,maxk);
    LT=zeros(n_Sessions,maxk-mink+1,maxk);
    
    for k=1:length(rangeK)
        for s=1:n_Sessions
            
            % Select the time points representing this subject and task
            T=(Time_all==s);
            Ctime=Kmeans_results{k}.IDX(T);
            
            for c=1:rangeK(k)
                % Probability
                P(s,k,c)=mean(Ctime==c);
                
                % Mean Lifetime
                Ctime_bin=Ctime==c;
                
                % Detect switches in and out of this state
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);
                
                % We discard the cases where state sarts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end
                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                else
                    C_Durations=0;
                end
                LT(s,k,c)=mean(C_Durations)*TR;
            end
        end
    end
    
    P_pval=zeros(3,maxk-mink+1,maxk);
   
    LT_pval=zeros(3,maxk-mink+1,maxk);
    
    % save the full output from estimation
    save(fullfile(Outdir, sprintf('LEiDA_results_k%d.mat',maxk)), '-v7.3');
  
    disp('%%%%% LEiDA SUCCESSFULLY COMPLETED %%%%%%%')
    disp('LEiDA results are now saved.')
    
    
end

% NED