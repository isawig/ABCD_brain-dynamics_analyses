%% FCH validation with NeuroMap

% Compares via correlation the first 6 low-frequency FCHs with seven
% reference gradients by NeuroMap. It uses the python package nibabel.

% Isabella L.C. Mariani Wigley, 09 / 2025; ilmawi@utu.fi
% Aurora Berto, 09 / 2025; aurber@utu.fi

clear; close all; clc

% preps; source the neuromaps virtual environment
% CMD in terminal
% source /Users/jetsonite/neuromaps/venv_neuromaps/bin/activate
% then open Matlab from that terminal session
% open /Applications/MATLAB_R2022b.app 

%% Import packages and load harmonic files

% import nibabel
py.importlib.import_module('nibabel')

% load harmonics nifti file
har = py.nibabel.load('/path/to/FCH/FCH_NeuroMap/Harmonics_as_Schaefer_nifti/Harmonic_7/NEUROmap_Harmonic_7.nii.gz') % USER.adapt!

% import neuromaps
py.importlib.import_module('neuromaps')

% ans = 
% 
%   Python module with properties:
% 
%          resampling: [1×1 py.module]
%     resample_images: [1×1 py.function]
%      compare_images: [1×1 py.function]
%               utils: [1×1 py.module]
%              images: [1×1 py.module]
%               stats: [1×1 py.module]
%          transforms: [1×1 py.module]
%            datasets: [1×1 py.module]

%% margulies2016 - fcgradient01

% load the annotation; margulies2016
margulies = py.neuromaps.datasets.fetch_annotation(source='margulies2016')
py_key = py.tuple({'margulies2016', 'fcgradient01', 'fsLR', '32k'})
grad = margulies{py_key};

% import neuromaps
py.importlib.import_module('neuromaps.transforms')

% transform harmonic to fslr 32k
har_res = py.neuromaps.transforms.mni152_to_fslr(har, '32k');

% estimate the correlation
corr = py.neuromaps.stats.compare_images(har_res, grad)

% load the surface data harmonics
data = py.neuromaps.images.load_data(har_res)

rotated = py.neuromaps.nulls.alexander_bloch(data, atlas='fsLR', density='32k')

% calculate correlation p-value
corr, pval = py.neuromaps.stats.compare_images(har_res, grad, nulls=rotated)

% save the stats in matlab format
correlations = double(py.array.array('d', py.numpy.nditer(corr))); 
pvalue = double(py.array.array('d', py.numpy.nditer(pval{2})));

%% margulies2016 - fcgradient02

% load the annotation; margulies2016
margulies = py.neuromaps.datasets.fetch_annotation(source='margulies2016')
py_key = py.tuple({'margulies2016', 'fcgradient02', 'fsLR', '32k'})
grad = margulies{py_key};

% import neuromaps
py.importlib.import_module('neuromaps.transforms')

% transform harmonic to fslr 32k
har_res = py.neuromaps.transforms.mni152_to_fslr(har, '32k');

% estimate the correlation
corr = py.neuromaps.stats.compare_images(har_res, grad)

% load the surface data harmonics
data = py.neuromaps.images.load_data(har_res)

rotated = py.neuromaps.nulls.alexander_bloch(data, atlas='fsLR', density='32k')

% calculate correlation p-value
corr, pval = py.neuromaps.stats.compare_images(har_res, grad, nulls=rotated)

% save the stats in matlab format
correlations = double(py.array.array('d', py.numpy.nditer(corr))); 
pvalue = double(py.array.array('d', py.numpy.nditer(pval{2})));

%% margulies2016 - fcgradient03

% load the annotation; margulies2016
margulies = py.neuromaps.datasets.fetch_annotation(source='margulies2016')
py_key = py.tuple({'margulies2016', 'fcgradient03', 'fsLR', '32k'})
grad = margulies{py_key};

% import neuromaps
py.importlib.import_module('neuromaps.transforms')

% transform harmonic to fslr 32k
har_res = py.neuromaps.transforms.mni152_to_fslr(har, '32k');

% estimate the correlation
corr = py.neuromaps.stats.compare_images(har_res, grad)

% load the surface data harmonics
data = py.neuromaps.images.load_data(har_res)

rotated = py.neuromaps.nulls.alexander_bloch(data, atlas='fsLR', density='32k');

% calculate correlation p-value
corr, pval = py.neuromaps.stats.compare_images(har_res, grad, nulls=rotated)

% save the stats in matlab format
correlations = double(py.array.array('d', py.numpy.nditer(corr))); 
pvalue = double(py.array.array('d', py.numpy.nditer(pval{2})));


%% margulies2016 - fcgradient04

% load the annotation; margulies2016
margulies = py.neuromaps.datasets.fetch_annotation(source='margulies2016')
py_key = py.tuple({'margulies2016', 'fcgradient04', 'fsLR', '32k'})
grad = margulies{py_key};

% import neuromaps
py.importlib.import_module('neuromaps.transforms')

% transform harmonic to fslr 32k
har_res = py.neuromaps.transforms.mni152_to_fslr(har, '32k');

% estimate the correlation
corr = py.neuromaps.stats.compare_images(har_res, grad)

% load the surface data harmonics
data = py.neuromaps.images.load_data(har_res)

rotated = py.neuromaps.nulls.alexander_bloch(data, atlas='fsLR', density='32k');

% calculate correlation p-value
corr, pval = py.neuromaps.stats.compare_images(har_res, grad, nulls=rotated)

% save the stats in matlab format
correlations = double(py.array.array('d', py.numpy.nditer(corr))); 
pvalue = double(py.array.array('d', py.numpy.nditer(pval{2})));


%% reardon2018 - scalinghcp

% load the annotation; reardon2018
reardon = py.neuromaps.datasets.fetch_annotation(source='reardon2018')
py_key = py.tuple({'reardon2018', 'scalinghcp', 'civet', '41k'})
grad = reardon{py_key};

% import neuromaps
py.importlib.import_module('neuromaps.transforms')

% transform harmonic to civet 41k
har_res = py.neuromaps.transforms.mni152_to_civet(har, '41k');

% estimate the correlation
corr = py.neuromaps.stats.compare_images(har_res, grad)

% load the surface data harmonics
data = py.neuromaps.images.load_data(har_res)
rotated = py.neuromaps.nulls.alexander_bloch(data, atlas='civet', density='41k');

% calculate correlation p-value
corr, pval = py.neuromaps.stats.compare_images(har_res, grad, nulls=rotated)

% save the stats in matlab format
correlations = double(py.array.array('d', py.numpy.nditer(corr))); 
pvalue = double(py.array.array('d', py.numpy.nditer(pval{2})));


%% abagen - genepc1

% load the annotation; abagen
abagen = py.neuromaps.datasets.fetch_annotation(source='abagen')

% import neuromaps
py.importlib.import_module('neuromaps.transforms')

% transform harmonic to fsaverage 10k
har_res = py.neuromaps.transforms.mni152_to_fsaverage(har, '10k');

% estimate the correlation
corr = py.neuromaps.stats.compare_images(har_res, abagen)

% load the surface data harmonics
data = py.neuromaps.images.load_data(har_res)
rotated = py.neuromaps.nulls.alexander_bloch(data, atlas='fsaverage', density='10k');

% calculate correlation p-value
corr, pval = py.neuromaps.stats.compare_images(har_res, abagen, nulls=rotated)

% save the stats in matlab format
correlations = double(py.array.array('d', py.numpy.nditer(corr))); 
pvalue = double(py.array.array('d', py.numpy.nditer(pval{2})));


%% neurosynth

% load the annotation; neurosynth
neurosynth = py.neuromaps.datasets.fetch_annotation(source='neurosynth')

% import neuromaps
py.importlib.import_module('neuromaps.transforms')

% transform harmonic to fsavg 10k
har_res = py.neuromaps.transforms.mni152_to_fsaverage(har, '10k')

% transform neurosynth to fsavg 10k
neurosynth_res = py.neuromaps.transforms.mni152_to_fsaverage(neurosynth, '10k')

% estimate the correlation
corr = py.neuromaps.stats.compare_images(har_res, neurosynth_res)

% load the surface data harmonics
data = py.neuromaps.images.load_data(har_res)
rotated = py.neuromaps.nulls.alexander_bloch(data, atlas='fsaverage', density='10k')

% calculate correlation p-value
corr, pval = py.neuromaps.stats.compare_images(har_res, neurosynth_res, nulls=rotated)

% save the stats in matlab format
correlations = double(py.array.array('d', py.numpy.nditer(corr))); 
pvalue = double(py.array.array('d', py.numpy.nditer(pval{2})));
