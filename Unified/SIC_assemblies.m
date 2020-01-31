%% ICA Extraction of Network States: Subject Assemblies
%	This script computes the neural assemblies and assembly time courses
% from the time series extracted from the BOLD data.  In addition, it
% computes the assemblywise Shannon entropy of each subject, and computes
% an upper bound of the total assemblywise entropy of each condition.
%	This version of the script calculates the assemblies of each subject
% separately.  This allows us to assess how the  assemblies change between
% subjects, and whether one group displays greater regularity than another.


%%	SETUP
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
path{6,1} = fullfile(path{1},'MATLAB','FastICA');
path{7,1} = fullfile(path{2},'Functions','LEICA');
path{8,1} = fullfile(path{2},'Results','LEICA');

% Add relevant paths
for k = 6:numel(path)-1
	addpath(genpath(path{k}));
end
clear k

% Load data
loadFile = 'LEICA90_Data_Iteration1.mat';
load(fullfile(path{8}, loadFile));
l = labels.Properties.VariableNames;

% Reset paths (in case original paths overwritten by loaded file)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');
path{4,1} = fullfile(path{2},'Data');
path{5,1} = fullfile(path{2},'Results');
path{6,1} = fullfile(path{1},'MATLAB','FastICA');
path{7,1} = fullfile(path{2},'Functions','LEICA');
path{8,1} = fullfile(path{2},'Results','LEICA');

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_SIC_Assemblies');
clear loadFile S


%% 1) Find significant components, percentage of explained variance

% Extract number of assemblies using Marcenko-Pastur distribution
N.IC = nan(max(N.subjects), N.conditions);
explained = cell(max(N.subjects), N.conditions);
tVar = nan(max(N.subjects), N.conditions);
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		N.IC(s,c) = NumberofIC(dFC.subj{s,c});				% find number of components
		[~,~,~,~,explained{s,c},~] = pca(dFC.subj{s,c}');	% find percentage of explained variance per component
		tVar(s,c) = sum(explained{s,c}(1:N.IC(s,c)));		% amount of variance captured with N.IC
	end
end
explained = tVar;
clear s c tVar

% Test for significant difference in number of ICs per group
[h,p] = kstest2(N.IC(:,1), N.IC(:,2));

% Visualization
F(N.fig) = figure;
N.fig = N.fig + 1;

% Visualize variance explained per IC
x = repmat(mean(N.IC,1, 'omitnan'), 19, 1);
y = 0:18;
col = {'b', 'r'};
subplot(1, 2, 1); hold on;
for c = 1:N.conditions
	s(c) = histogram(N.IC(:,c), 'FaceColor',col{c});
	plot(x(:,c), y, strcat(':', col{c}));
end
if h == true
	plot(mean(x(1,:)), 16, '*k')
end
title('Number of ICs');
ylabel('Counts');
legend(s, l, 'Location','northwest');

subplot(1, 2, 2); hold on;
for c = 1:N.conditions
	histogram(explained(:,c));
end
title('% Captured Variance');
ylabel('Counts');
legend(l);
clear c s



%% 2) Compute the assemblies & assembly activations

disp('Processing the ICs from BOLD data');

% Compute assembly activity timecourses and memberships
activities = cell(max(N.subjects), N.conditions);
memberships = cell(max(N.subjects), N.conditions);
W = cell(max(N.subjects), N.conditions);
ICs = cell(max(N.subjects), N.conditions);
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		[activities{s,c}, memberships{s,c}, W{s,c}] = fastica(dFC.subj{s,c}, 'numOfIC', N.IC(s,c), 'verbose','off');

		% Normalize membership weights
		for i = 1:N.IC(s,c)
			memberships{s,c}(:,i) = memberships{s,c}(:,i)./norm(memberships{s,c}(:,i));
		end

		% Compute IC matrices
		ICs{s,c} = nan(N.ROI, N.ROI, N.IC(s,c));
		for i = 1:size(memberships,2)
			ICs{s,c}(:,:,i) = memberships{s,c}(:,i) * memberships{s,c}(:,i)';
		end
		clear i
	end
end
clear s c



%% Save results

% Get file list
fList = dir(fullfile(path{8}, strcat(fileName, '*')));

% Find number of previous iterations
nIter = 1+numel(fList);
clear fList

% Edit fileName
fileName = strcat(fileName, '_Iteration', num2str(nIter));
clear nIter

% Save figures
savefig(F, fullfile(path{8}, fileName), 'compact');
clear F

% Save variables
save(fullfile(path{8}, fileName));


