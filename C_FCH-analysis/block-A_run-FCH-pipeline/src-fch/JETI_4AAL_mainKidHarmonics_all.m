%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Selen Atasoy
%
% 21/12/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%% local set up and SPM integration
%%%%%%%%%%%%% added by Jetro to make Matlab work (SPM has issues otherwise)
%%%%%%% for MacBook Pro, MAC OS Catalina
addpath(genpath('/Users/jetrojtuulari/Documents/spm12'));
addpath('/Users/jetrojtuulari/Documents/spm12_fieldtrip');

%% to get the run time
%tic

%% input folder
%folder.in = '/home/c/csatas/Research/Results/KidHarmonics/5V_2020-11-17-HARMONICS-SAMPLE/functional_harmonics';
folder.in = '/Users/jetrojtuulari/Desktop/tryout-func-harmo-5v';

%% input files
files.data.left      = 'filtered_func_data_clean_nonaggr_tri.L.func.gii';
files.data.right     = 'filtered_func_data_clean_nonaggr_tri.R.func.gii';

%% subject to process %%
subject = '7043234';
disp(['Processing subject ', subject, '...']);
    
folder.data         = fullfile(folder.in, subject);

%% if subject's folder doesnt exist, create one
if ~exist(folder.data, 'dir')
    disp('Subject directory doesnt exist, please set the path to the files.. ', files.data.left, ', and ', files.data.right);
else
    
    %% read the data
    [data_norm] = surfFMRI_readData(fullfile(folder.data, files.data.left), fullfile(folder.data, files.data.right), 1);
    % data_norm.all = data_norm.all(1:5000,:); % jjtuul: commented as it
    % was added by Selen only for testing purposes 

    %% compute the FC
    FC = pdist(data_norm.all, 'correlation');

    % to test the correlations
    %[RHO,PVAL] = corr(data_norm.all(v,:)', data_norm.all');

    %% save correlations
    % save(fullfile(folder.data, 'FC.mat'), 'FC'); 
    % jjtuul: added a version definition to save a large file '-v7.3'
    save(fullfile(folder.data, 'FC.mat'), 'FC', '-v7.3'); 
    
    %% compute the adjacency
    nn = 30 % number of nearest neighbours to create the adjacency
    [A] = computeAdjacency(squareform(FC), 'nn', nn, 0);
 
    %% save adjacency
    % jjtuul: added a version definition to save a large file '-v7.3'
    save(fullfile(folder.data, 'Adjacency.mat'), 'A', '-v7.3');
    clear FC; % clean

    %% compute the harmonics
    nr_vertices = 32492*2; 
    
    % get rid of the subcortical structures and corpus collosum for the
    % harmonic computation
    load('connRSM_CC.mat'); % load the indices of corpus callosum and other subcortical structures
    CC.inds         = setdiff([1:nr_vertices], CC.RestInds);
    A(CC.inds,:)   = [];
    A(:,CC.inds)    = [];
    
    % initialise
    Y = zeros(nr_vertices, size(A,2)-1);
    
    % harmonic estimation
    [Y(CC.RestInds,:), V] = computeLaplacianEigenmaps(A, size(A,2)-1, 'symmetric');
    % jjtuul: added a version definition to save a large file '-v7.3'
    save(fullfile(folder.data, 'Harmonics.mat'), 'Y', 'V', 'CC', '-v7.3');

%% to get the run time
%toc

end
   


    