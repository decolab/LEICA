%% ICA Extraction of Network States: Common Assemblies
%	This script computes the neural assemblies and assembly time courses
% from the time series extracted from the BOLD data.  In addition, it
% computes the assemblywise Shannon entropy of each subject, and computes
% an upper bound of the total assemblywise entropy of each condition.
%	This version of the script calculates the assemblies using the time
% series of both conditions.  This prevents the membership of the resting
% state networks from changing between conditions, but ensures the
% comparability of the time series and Shannon entropies between
% conditions.


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
path{6,1} = fullfile(path{1},'MATLAB','FastICA');
path{7,1} = fullfile(path{2},'Functions','LEICA');
path{8,1} = fullfile(path{2},'Results','LEICA');

% Add relevant paths
for k = 6:numel(path)-1
	addpath(genpath(path{k}));
end
clear k


%% Set file names & load data

% Define files to load
loadFile = 'LEICA90_Data_Iteration1.mat';

% Load data
load(fullfile(path{8}, loadFile));

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_CIC_Assemblies');
clear loadFile S


%% Reset paths (in case original paths overwritten by loaded file)

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


%% 1) Find significant components, percentage of explained variance

% Extract number of assemblies using Marcenko-Pastur distribution
N.assemblies = NumberofIC(dFC.concat);

% Confirm with PCA
[~,~,~,~,explained,~] = pca(dFC.concat');

% Quantify amount of variance captured with N.assemblies
tVar = sum(explained(N.assemblies+1:end));

% Visualize variance explained per PCA
F(N.fig) = figure;
plot(1:numel(explained), explained, ':ob'); hold on;
bar(explained, 'c');
c = plot(N.assemblies, explained(N.assemblies), '*r', 'MarkerSize',8);
title('% Variance Captured per PC')
xlabel('Principal Components');
ylabel('Variance Captured');
legend(c, [num2str(tVar), '% variance captured', newline, 'with N = ', num2str(N.assemblies), ' components.'])

N.fig = N.fig + 1;
clear explained tVar


%% 2) Compute the assemblies & assembly activations

disp('Processing the ICs from BOLD data');

% Compute assembly activity timecourses and memberships
[activities.concat, memberships, W] = fastica(dFC.concat, 'numOfIC', N.assemblies, 'verbose','off');


%% 4) Separate assembly activations by condition & subject

% Declare storage arrays
activities.cond = cell(1, N.conditions);
activities.subj = cell(max(N.subjects), N.conditions);

% Separate assembly activations by condition & subject
for c = 1:N.conditions
	
	I = T.index(2,:) == c;
	activities.cond{c} = activities.concat(:,I);
	
	for s = 1:N.subjects(c)
		I = (T.index(2,:) == c & T.index(1,:) == s);
		activities.subj{s,c} = activities.concat(:,I);
	end
end
clear I s c


%% 5) Normalize & visualize membership weights

% Normalize memberships
for k = 1:N.assemblies
	memberships(:,k) = memberships(:,k)./norm(memberships(:,k));
end

% Visualize membership weights
F(N.fig) = figure;
for k = 1:size(memberships, 2)
	subplot(1,size(memberships,2),k)
	
	I = memberships(:,k) > 0;
	m = memberships(:,k); m(I) = 0;
	barh(1:numel(I), m, 'r'); hold on;
	
	I = memberships(:,k) <= 0;
	m = memberships(:,k); m(I) = 0;
	barh(1:numel(I), m, 'g');
	
	title(['IC ', num2str(k)]);
	ylim([0 numel(I)+1]);
	xlim([-max(abs(memberships(:,k))) max(abs(memberships(:,k)))]);
end
N.fig = N.fig+1;
clear k m n I


%% Compute LEICA assembly matrices

% Preallocate storage array
ICs = nan(N.ROI, N.ROI, N.assemblies);

% Compute IC matrices
for i = 1:size(memberships,2)
	ICs(:,:,i) = memberships(:,i) * memberships(:,i)';
end
clear i


%% Compare LEICA assemblies to RSNs

% Load Yeo RSNs
load(fullfile(path{4}, 'Atlases/AAL', 'RSNsAAL.mat'));
RSNs = LR_version_symm(Yeo_AAL');
clear Yeo_AAL

% Generate comparison types, figure titles
cTypes = {'correlation', 'spearman', 'cosine'};
cTitles = {'Pearson Distance', 'Spearman Distance', 'Cosine Distance'};

% Visualize relations between ICs and RSNs: Pearson, Spearman, cosine
F(N.fig) = figure; hold on;
for k = 1:3
	subplot(1,3,k); colormap jet
	imagesc(pdist2(memberships', RSNs', cTypes{k})); colorbar;
	title(cTitles{k});
	xlabel('Yeo RSNs');
	ylabel('LEICA ICs');
end
clear cTypes cTitles k
N.fig = N.fig+1;


% Visualize component relationships
F(N.fig) = figure; colormap jet

subplot(2,2,[1,2]);
imagesc(activities.concat); colorbar;
title('Component Activation')
xlabel('Time Points');
ylabel('Components');

% Component spatial correlation
subplot(2,2,3);
comp.spatial = corr(memberships);
F(N.fig) = figure;
imagesc(comp.spatial); colorbar;
title('Component Spatial Correlation')

% Component temporal overlap
subplot(2,2,4);
imagesc(corr(activities.concat')); colorbar;
title('Inter-Assembly Pearson Correlations');
clear k comp tempTitles

% Count
N.fig = N.fig+1;


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
savefig(F, fullfile(path{8}, fileName), 'compact');
clear F

% Save variables
save(fullfile(path{8}, fileName));


