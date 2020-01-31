%% ICA Common Assemblies: Computing Component Metrics
%	


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ASSEMBLY PIPELINE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%	SETUP

%% Set paths & filenames

clear; close all; clc

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
path{6,1} = fullfile(path{2},'Functions','LEICA');
path{7,1} = fullfile(path{2},'Results','LEICA');

% Add relevant paths
addpath(genpath(path{6}));


%% Set file names & load data

% Define files to load
loadFile = 'LEICA90_CIC_Assemblies';

% Load data
load(fullfile(path{7}, loadFile));

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_', S{2}, '_Metrics');
clear loadFile S

% Reset N.fig
N.fig = 1;


%% Reset paths (in case original paths overwritten by loaded file)

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2},'Data');
path{5,1} = fullfile(path{2},'Results');
path{6,1} = fullfile(path{2},'Functions','LEICA');
path{7,1} = fullfile(path{2},'Results','LEICA');



%% Compute activity metrics

% Preallocate storage arrays:
%	Activation magnitude means, medians, standard deviations
%	Component-wise Kuramoto order parameter & metastability
%	Subject-wise entropies
activities.av.cond = nan(N.assemblies, N.conditions);
activities.md.cond = nan(N.assemblies, N.conditions);
activities.sd.cond = nan(N.assemblies, N.conditions);
metastable.cond = nan(1, N.conditions);
kuramoto.cond = nan(N.conditions, max(T.condition));
activities.av.subj = nan(N.assemblies, max(N.subjects), N.conditions);
activities.md.subj = nan(N.assemblies, max(N.subjects), N.conditions);
activities.sd.subj = nan(N.assemblies, max(N.subjects), N.conditions);
metastable.subj = nan(max(N.subjects), N.conditions);
kuramoto.subj = nan(max(N.subjects), N.conditions, T.scan);
entro.subj = nan(N.assemblies, N.subjects, N.conditions);

for c = 1:N.conditions
	activities.av.cond(:,c) = mean(activities.cond{c}, 2, 'omitnan');
	activities.md.cond(:,c) = median(activities.cond{c}, 2, 'omitnan');
	activities.sd.cond(:,c) = std(activities.cond{c}, 0, 2, 'omitnan');
	[kuramoto.cond(c, 1:T.condition(c)), metastable.cond(c)] = findStability(activities.cond{c});
	for s = 1:N.subjects(c)
		% 
		activities.av.subj(:,s,c) = mean(activities.subj{s,c}, 2, 'omitnan');
		activities.md.subj(:,s,c) = median(activities.subj{s,c}, 2, 'omitnan');
		activities.sd.subj(:,s,c) = std(activities.subj{s,c}, 0, 2, 'omitnan');
		[kuramoto.subj(s,c,:), metastable.subj(s,c)] = findStability(activities.subj{s,c});
		for ass = 1:N.assemblies
			entro.subj(ass, s, c) = HShannon_kNN_k_estimation(activities.subj.TS{s,c}(ass,:), co);
		end
	end
end
clear c s
entro.subj = squeeze(sum(entro.subj, 1));

% Convert metrics to table format
activities.av.cond = array2table(activities.sv.cond, 'VariableNames',labels.Properties.VariableNames);
activities.md.cond = array2table(activities.md.cond, 'VariableNames',labels.Properties.VariableNames);
activities.sd.cond = array2table(activities.sd.cond, 'VariableNames',labels.Properties.VariableNames);
metastable.cond = array2table(metastable.cond, 'VariableNames', labels.Properties.VariableNames);
metastable.subj = array2table(metastable.subj, 'VariableNames', labels.Properties.VariableNames);
entro.subj = array2table(entro.subj, 'VariableNames', labels.Properties.VariableNames);



%% Save results

% Get file list
fList = dir(fullfile(path{8}, strcat(fileName, '*')));

% Find number of previous iterations
nIter = numel(fList);
clear fList

% Edit fileName
fileName = strcat(fileName, '_Iteration', num2str(nIter));
clear nIter

% Save variables
save(fullfile(path{8}, fileName));


