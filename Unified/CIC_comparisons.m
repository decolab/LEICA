%% Comparison of Metric Distributions between States
%	


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% COMPARISON PIPELINE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%	SETUP

%% Set paths & filenames

clear; close all; clc
set(groot,'defaultLegendAutoUpdate','off');	% prevents legend from automatically adding new data

% Shuffle random seed.  Necessary in array parallelization to avoid
% repeating same random seed across arrays.
rng('shuffle');

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2},'Data');
path{5,1} = fullfile(path{2},'Results');
path{6,1} = fullfile(path{1},'MATLAB', 'spm12');
path{7,1} = fullfile(path{2},'Functions','LEICA');
path{8,1} = fullfile(path{2},'Results', 'LEICA');

% Add relevant paths
addpath(path{6});
addpath(genpath(path{7}));


%% Set file names & load data

% Define files to load
loadFile = 'LEICA90_CIC_Metrics';

% Load data
% metricFig = openfig(fullfile(path{8}, loadFile));
load(fullfile(path{8}, loadFile));

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_', S{2}, '_Comparisons');
clear loadFile S


%% Reset paths overwritten by loaded file

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2},'Data');
path{5,1} = fullfile(path{2},'Results');
path{6,1} = fullfile(path{1},'MATLAB', 'spm12');
path{7,1} = fullfile(path{3},'Functions');
path{8,1} = fullfile(path{2},'Results', 'LEICA');




%% Compare transition matrices between states

% Preallocate storage arrays


% Compute distance(s) between matrices



%% Compare activations between conditions

% Compare ROI activation time series
sig.AAL = robustTests(dFC.cond{1}, dFC.cond{2}, N.assemblies, pval.target, 'kstest2');

% Compare component activation time series
sig.IC = robustTests(activities.cond{1}, activities.cond{2}, N.ROI, pval.target, 'kstest2');



%% Compare activation metrics between conditions

% Average activation magnitude(s)
con = activities.av.subj(:,:,1); con = con(isfinite(con));
pat = activities.av.subj(:,:,2); pat = pat(isfinite(pat));
sig.av = robustTests(con, pat, N.assemblies, pval.target, 'kstest2');
clear con pat

% Activation medians(s)
con = activities.md.subj(:,:,1); con = con(isfinite(con));
pat = activities.md.subj(:,:,2); pat = pat(isfinite(pat));
sig.md = robustTests(con, pat, N.assemblies, pval.target, 'kstest2');
clear con pat

% Activation standard deviations(s)
con = activities.sd.subj(:,:,1); con = con(isfinite(con));
pat = activities.sd.subj(:,:,2); pat = pat(isfinite(pat));
sig.sd = robustTests(con, pat, N.assemblies, pval.target, 'kstest2');
clear con pat

% Subject metastabilities
con = metastable.subj{:,'Controls'}(isfinite(metastable.subj{:,'Controls'}));
pat = metastable.subj{:,'Patients'}(isfinite(metastable.subj{:,'Patients'}));
if adtest(con) && adtest(pat)
	[sig.metastable.h, sig.metastable.p] = kstest2(con, pat, 'Alpha',pval.target);
else
	[sig.metastable.h, sig.metastable.p] = ttest2(con, pat, 'Alpha',pval.target);
end
clear con pat

% Subject entropies
con = entro.subj{:,'Controls'}(isfinite(entro.subj{:,'Controls'}));
pat = entro.subj{:,'Patients'}(isfinite(entro.subj{:,'Patients'}));
if adtest(con) && adtest(pat)
	[sig.entro.h, sig.entro.p] = kstest2(con, pat, 'Alpha',pval.target);
else
	[sig.entro.h, sig.entro.p] = ttest2(con, pat, 'Alpha',pval.target);
end
clear con pat



%% Save results

% Get file list
fList = dir(fullfile(path{8}, strcat(fileName, '*')));

% Find number of previous iterations
nIter = numel(fList);
clear fList

% Edit fileName
fileName = strcat(fileName, '_Iteration', num2str(nIter));
clear nIter

% Save figures
savefig(metricFig, fullfile(path{8}, fileName), 'compact')
clear metricFig

% Save variables
save(fullfile(path{8}, fileName));

