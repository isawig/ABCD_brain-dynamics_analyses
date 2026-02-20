%% Local script A to run LEiDA on HPC Puhti (CSC)
% This script is a configuration job submitted from MATLAB.

% Jetro J. Tuulari, 01 / 2020, jetro.tuulari@utu.fi
% Isabella L.C. Mariani Wigley, 04 / 2025, ilmawi@utu.fi
% Aurora Berto, 04 / 2025, aurber@utu.fi

% The settings are reset by running the command "ConfigCluster", after which you 
% will receive a prompt to input your CSC user name.
% You may want to do this before starting new runs of analyses to empty the
% job history.

% Please always consider the parameters below carefully
% INSTRUCTIONS ONLINE: https://docs.csc.fi/#apps/matlab/

% make the tools available
addpath("~/Documents/MATLAB/mps_puhti")
% savepath

configCluster

% TIP: you can connect (ssh) to Puhti to see your home directory and also
% to monitor the queue online (this seems to update slowly to Matlab)
% after you have the ssh ready type "squeue -l -u $USER" to see the queue
% to see the available licenses (for workers) type the following
% "scontrol show lic=mdcs"

% get a handle on the cluster
c = parcluster;

% days-hh:mm:s, e.g. 2h = 02:00:0 | 2d12h = 2-12:00:0
c.AdditionalProperties.WallTime = '02:00:0';
% no. CPU's
c.AdditionalProperties.CPUsPerNode = '';
% memory usage, this will per node / worker for parallel or multiple cpu runs
% Set memory per CPU (4 GB in this case)
c.AdditionalProperties.MemPerCPU = '8g';  % Change to '8g' or desired value
c.AdditionalProperties.MemUsage = '2g';
% partition on the cluster, see: https://docs.csc.fi/#computing/running/batch-job-partitions/
c.AdditionalProperties.QueueName = 'small';
% this is the account name for HPC project
c.AdditionalProperties.ComputingProject = 'project_ID'; % your project ID
c.AdditionalProperties.AccountName = 'project_ID';      % your project ID
% email to receive notification to, default is begin and end i.e. "ALL"
c.AdditionalProperties.EmailAddress = 'your@email.com'; % your email
% Check configured values
c.AdditionalProperties
% GPU stuff
c.AdditionalProperties.GpuCard = '';
c.AdditionalProperties.GPUsPerNode = '';

% save profile for the submitted scripts
c.saveProfile
% to check currently running and finished jobs
% jobs = c.Jobs
