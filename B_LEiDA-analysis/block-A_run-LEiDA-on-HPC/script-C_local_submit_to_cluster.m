%% Local script B to run LEiDA on HPC Puhti (CSC)
% This script is to submit the work to the cluster.

% Jetro J. Tuulari, 01 / 2020, jetro.tuulari@utu.fi
% Isabella L.C. Mariani Wigley, 04 / 2025, ilmawi@utu.fi
% Aurora Berto, 04 / 2025, aurber@utu.fi

j = batch(c, @script_CSC_compute_leida_from_bold, 1, {}, 'CurrentFolder', '.', 'AutoAttachFiles',false, 'AutoAddClientPath', false, 'pool', 39);

% NED
