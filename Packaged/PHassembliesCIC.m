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
loadFile = 'LEICA90_PhaseData';

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


%% 1) Find number of components

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
[activities.concat.TS, memberships, W] = fastica(dFC.concat, 'numOfIC', N.assemblies, 'verbose','off');

% Compute activation, event time courses
activities.concat.activations = activeMat(activities.concat.TS, 1);
activities.concat.events = eventMat(activities.concat.activations);


%% 4) Separate assembly activations by condition & subject

% Declare storage arrays
activities.cond.TS = cell(1, N.conditions);
activities.cond.activations = cell(1, N.conditions);
activities.cond.events = cell(1, N.conditions);
activities.subj.TS = cell(max(N.subjects), N.conditions);
activities.subj.activations = cell(max(N.subjects), N.conditions);
activities.subj.events = cell(max(N.subjects), N.conditions);

% Separate assembly activations by condition & subject
for c = 1:N.conditions
	I = T.index(2,:) == c;
	activities.cond.TS{c} = activities.concat.TS(:,I);
	activities.cond.activations{c} = activities.concat.activations(:,I);
	activities.cond.events{c} = activities.concat.events(:,I);
	
	for s = 1:N.subjects(c)
		I = (T.index(2,:) == c & T.index(1,:) == s);
		activities.subj.TS{s,c} = activities.concat.TS(:,I);
		activities.subj.activations{s,c} = activities.concat.activations(:,I);
		activities.subj.events{s,c} = activities.concat.events(:,I);
	end
end
clear I s c


%% 5) Normalize & visualize membership weights

% Normalize memberships
memberships = memberships./max(abs(memberships));

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
	xlim([-1 1]);
end
N.fig = N.fig+1;
clear k m n I


%% Visualize component activation over time

F(N.fig) = figure;
colormap jet
imagesc(activities.concat.TS); colorbar; hold on;
title('Component Activation')
xlabel('Time Points');
ylabel('Components');
N.fig = N.fig + 1;


%% Compute component spatial overlap

comp.spatial = corr(memberships);
F(N.fig) = figure;
imagesc(comp.spatial); colorbar; hold on;
title('Component Spatial Correlation')
N.fig = N.fig + 1;


%% Compute component temporal overlap

% Set titles
tempTitles = {'Pearson Correlation', 'Proportion of Co-Activity', 'Proportion of Co-Activation'};

% Compute temporal correlation & co-activity
comp.temporal{1} = corr(activities.concat.TS');
comp.temporal{2} = zeros(N.assemblies, N.assemblies, sum(T.condition));
comp.temporal{3} = zeros(N.assemblies, N.assemblies, sum(T.condition));
for t = 1:sum(T.condition)
	comp.temporal{2}(:,:,t) = activities.concat.activations(:,t) * activities.concat.activations(:,t)';
	comp.temporal{3}(:,:,t) = activities.concat.events(:,t) * activities.concat.events(:,t)';
end
clear t

% Extract activation and event probabilities
activities.concat.prob.sd = diag(std(comp.temporal{2}, 0, 3));
events.concat.prob.sd = diag(std(comp.temporal{3}, 0, 3));

% Find mean correlation and co-activation matrix over time
comp.temporal{2} = mean(comp.temporal{2}, 3);
comp.temporal{3} = mean(comp.temporal{3}, 3);

% Extract mean probability of activity & activation
activities.concat.prob.av = diag(comp.temporal{2});
events.concat.prob.av = diag(comp.temporal{3});

% Visualize correlation, coactivity, co-activation matrices
F(N.fig) = figure;
colormap jet
title('Temporal Relations Between Assemblies');
for k = 1:3
	subplot(1,3,k);
	imagesc(comp.temporal{k}); colorbar; hold on;
	title(tempTitles{k});
end
clear k comp tempTitles
N.fig = N.fig+1;


%% Save results

% Save figures
savefig(F, fullfile(path{8}, fileName), 'compact');
clear F

% Save variables
save(fullfile(path{8}, fileName));


