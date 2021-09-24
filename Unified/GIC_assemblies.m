%% ICA Extraction of Network States: Group Assemblies
%	This script computes the neural assemblies and assembly time courses
% from the time series extracted from the BOLD data.  In addition, it
% computes the assemblywise Shannon entropy of each subject, and computes
% an upper bound of the total assemblywise entropy of each condition.
%	This version of the script calculates  assemblies using the time series
% of each condition separately.  This allows us to assess how the
% assemblies change between conditions.



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
fileName = strcat(S{1}, '_GIC_Assemblies');
clear loadFile S



%% 1) Find significant components, percentage of explained variance

% Visualize variance explained per IC
F(N.fig) = figure;
N.fig = N.fig + 1;

% Extract number of assemblies using Marcenko-Pastur distribution
N.IC = nan(N.conditions, 1);
explained = cell(N.conditions, 1);
tVar = nan(N.conditions, 1);
for c = 1:N.conditions
	N.IC(c) = NumberofIC(dFC.cond{c});				% find number of components
	[~,~,~,~,explained{c},~] = pca(dFC.cond{c}');	% find percentage of explained variance per component
	tVar(c) = sum(explained{c}(1:N.IC(c)));			% amount of variance captured with N.IC

	subplot(1, N.conditions, c);
	plot(1:numel(explained{c}), explained{c}, ':ob'); hold on;
	bar(explained{c}, 'c');
	p = plot(N.IC(c), explained{c}(N.IC(c)), '*r', 'MarkerSize',8);
	title('% Variance Captured per PC')
	xlabel('Principal Components');
	ylabel('Variance Captured');
	legend(p, [num2str(tVar(c)), '% variance captured', newline, 'with N = ', num2str(N.IC(c)), ' components.'])
end
clear explained tVar p



%% 2) Compute the assemblies & assembly activations

disp('Processing the ICs from BOLD data');

% Compute assembly activity timecourses and memberships
activities = cell(N.conditions, 1);
memberships = cell(N.conditions, 1);
W = cell(N.conditions, 1);
ICs = cell(N.conditions, 1);
for c = 1:N.conditions
	[activities{c}, memberships{c}, W{c}] = fastica(dFC.cond{c}, 'numOfIC', N.IC(c), 'verbose','off');
	
	% Normalize membership weights
	for k = 1:N.IC(c)
		memberships{c}(:,k) = memberships{c}(:,k)./norm(memberships{c}(:,k));
	end
	
	% Compute IC matrices
	ICs{c} = nan(N.ROI, N.ROI, N.IC(c));
	for i = 1:size(memberships,2)
		ICs{c}(:,:,i) = memberships{c}(:,i) * memberships{c}(:,i)';
	end
	clear i
end

% Project timeseries onto alternative components
projections = cell(nchoosek(N.conditions,1), 1);
projections{1} =  W{1}*dFC.cond{2};
projections{2} =  W{2}*dFC.cond{1};



%% 4) Visualize assemblies

l = labels.Properties.VariableNames;

% Visualize membership weights
F(N.fig) = figure;
N.fig = N.fig+1;
d = abs(N.IC - max(N.IC));
for c = 1:N.conditions
	for k = 1:N.IC(c)
		subplot(N.conditions, max(N.IC), k+(c-1)*max(N.IC)); hold on;

		I = memberships{c}(:,k) > 0;
		m = memberships{c}(:,k); m(I) = 0;
		barh(1:numel(I), m, 'r');

		I = memberships{c}(:,k) <= 0;
		m = memberships{c}(:,k); m(I) = 0;
		barh(1:numel(I), m, 'g');

		title([l{c}, ': IC ', num2str(k)]);
		ylim([0 numel(I)+1]);
		xlim([-max(abs(memberships{c}(:,k))) max(abs(memberships{c}(:,k)))]);
	end
end
clear k m n I


% Visualize component relationships
F(N.fig) = figure; colormap jet
N.fig = N.fig+1;
for c = 1:N.conditions
	i = (c-1)*N.conditions;
	v = i + [1,2];
	
	subplot(2, N.conditions*2, v);
	imagesc(activities{c}); colorbar;
	title([l{c}, ': Component Activation']);
	xlabel('Time Points');
	ylabel('Components');

	% Component spatial correlation
	subplot(2, N.conditions*2, 5+i);
	comp.spatial = corr(memberships{c});
	imagesc(comp.spatial); colorbar;
	title([l{c}, ': Component Spatial Correlation']);

	% Component temporal overlap
	subplot(2, N.conditions*2, 6+i);
	imagesc(corr(activities{c}')); colorbar;
	title([l{c}, ': IC Temporal Correlations']);
end
clear c v i


% Visualize inter-group compontent distances
F(N.fig) = figure; colormap jet
N.fig = N.fig+1;
cTypes = {'correlation', 'spearman', 'cosine'};
cTitles = {'Pearson Distance', 'Spearman Distance', 'Cosine Distance'};
title([l{1}, ' vs. ', l{2}]);
for k = 1:numel(cTypes)
	subplot(2, numel(cTypes), k);
	colormap jet
	imagesc(pdist2(memberships{1}', memberships{2}', cTypes{k})); colorbar;
	title(['Spatial ', cTitles{k}]);
	xlabel(l{2});
	ylabel(l{1});
	
	subplot(2, numel(cTypes), k+numel(cTypes));
	colormap jet
	imagesc(pdist2(activities{1}(:,1:T.scan*min(N.subjects)), activities{2}(:,1:T.scan*min(N.subjects)), cTypes{k})); colorbar;
	title(['Temporal ', cTitles{k}]);
	xlabel(l{2});
	ylabel(l{1});
end


% Compare LEICA assemblies to RSNs

% Load Yeo RSNs
load(fullfile(path{4}, 'Atlases/AAL', 'RSNsAAL.mat'));
RSNs = LR_version_symm(Yeo_AAL');
clear Yeo_AAL

% Visualize relations between ICs and RSNs: Pearson, Spearman, cosine
F(N.fig) = figure; hold on;
N.fig = N.fig+1;
for c = 1:N.conditions
	for k = 1:numel(cTypes)
		subplot(N.conditions, numel(cTypes), k+(c-1)*numel(cTypes));
		colormap jet
		imagesc(pdist2(memberships{c}', RSNs', cTypes{k})); colorbar;
		title(['ICs vs. RSNs: ', cTitles{k}]);
		xlabel('Yeo RSNs');
		ylabel([l{c}, ': ICs']);
	end
end
clear cTypes cTitles k c


% Visualize projection relationships
F(N.fig) = figure; colormap jet
N.fig = N.fig+1;
d = N.conditions:-1:1;
for c = 1:N.conditions
	i = (c-1)*N.conditions;
	v = i + [1,2];
	
	subplot(3, N.conditions*2, v);
	imagesc(projections{c}); colorbar;
	title(['Projection: ', l{d(c)}, ' onto ', l{c}, ' Components']);
	xlabel('Time Points');
	ylabel('Components');
	
	% Component temporal overlap
	subplot(3, N.conditions*2, [5 6 9 10]+i);
	imagesc(corr(activities{c}')); colorbar;
	title([l{c}, ': Projections Temporal Correlations']);
end
clear d c v i l



%% 5) Save results

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


