%% ICA Extraction of Network States: Group Assemblies
%	This script computes the neural assemblies and assembly time courses
% from the time series extracted from the BOLD data.  In addition, it
% computes the assemblywise Shannon entropy of each subject, and computes
% an upper bound of the total assemblywise entropy of each condition.
%	This version of the script calculates  assemblies using the time series
% of each condition separately.  This allows us to assess how the
% assemblies change between conditions.



%%	SETUP

% Clear workspace
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
path{6,1} = fullfile(path{2},'Results','LEICA');

% Add relevant paths
fpath{1,1} = fullfile(path{1},'MATLAB','spm12');
fpath{2,1} = fullfile(path{1},'MATLAB','FastICA');
fpath{3,1} = fullfile(path{1},'MATLAB','permutationTest');
fpath{4,1} = fullfile(path{2},'Functions');
fpath{5,1} = fullfile(path{3},'Functions');
for k = 1:numel(fpath)-1
	addpath(fpath{k});
end
addpath(genpath(fpath{numel(fpath)}));
clear fpath

% Load structural data
load(fullfile(path{4}, 'sc90.mat'));

% Set methods
phaseType = 'cosine';	% for measuring phase: cosine or exponential
compressType = 'eigenvector';	% for compressing matrix: eigenvector, average, or none

% File to save
fileName = 'ZICA90_GIC';
fList = dir(fullfile(path{6}, strcat(fileName, '_*')));	% Get file list
nIter = numel(fList); clear fList						% Find number of previous iterations
fileName = strcat(fileName, '_Iteration', num2str(nIter)); clear nIter	% Edit fileName

% Define list of patients
patientList = 'ClinicalData_OCD.xlsx';

% Define files containing BOLD data
loadDirectory = 'Subjects';
loadFiles = 'P*';

% Get file list
fList = dir(fullfile(path{4}, loadDirectory, loadFiles));

% Get number of files of interest
nFiles = numel(fList);

% Load patient data
patientData = readtable(fullfile(path{4}, patientList));

% Set data array sizes
D = load(fullfile(path{4}, loadDirectory, fList(1).name))';
N.ROI = size(D,1);
T.scan = size(D,2);

% Set data arrays
subjectBOLD = nan(N.ROI, T.scan, nFiles);
fName = cell(nFiles, 1);

% Select ROIs of interest (if necessary)
nROI = 90;
dROI = N.ROI - nROI;
I = [1:(nROI/2), (nROI/2+dROI+1):N.ROI]';

% Get functional data
for k = 1:nFiles
	
	% Rename files
	D = strsplit(fList(k).name, '_');
	fName{k} = D{1};
	
	% Visualize data being extracted
	disp(['Extracting data for subject ', fName{k}, '.']);
	
	% Load data
	subjectBOLD(:,:,k) = load(fullfile(path{4}, loadDirectory, fList(k).name))';
	
	% Convert LR-symmetric to LR-mirrored ordering
	subjectBOLD(:,:,k) = LR_version_symm(subjectBOLD(:,:,k));
	
end
clear D ans k loadDirectory loadFiles nFiles fList

% Remove cerebellar regions
if exist('I', 'var')
	subjectBOLD = subjectBOLD(I,:,:);
	N.ROI = size(subjectBOLD,1);
	clear I
end
clear nROI dROI

% Symmetrize node labels
load(fullfile(path{4},'Atlases','AAL','AAL_labels'));
labelROI = string(label90);
labelROI = LR_version_symm(labelROI);
clear label90

% Extract indices for BOLD signals of patients, controls
I(:,1) = ~ismember(fName, patientData{:, 'Code'});
I(:,2) = ismember(fName, patientData{:, 'Code'});
I = array2table(I, 'VariableNames',{'Controls','OCD'});

% Set number of conditions
N.conditions = size(I,2);

% Find number of subjects in each condition
for c = 1:N.conditions
	N.subjects(c) = nnz(I{:,c});
end

% Extract labels for patients, controls
labels = cell(max(N.subjects), N.conditions);
for c = 1:N.conditions
	labels(1:nnz(I{:,c}), c) = fName(I{:,c});
end
labels = cell2table(labels, 'VariableNames',{'Control','Patient'});
clear c fName



%% Extract dFC

% Temporal parameters
T.TR = 2.73;			% Repetition Time (seconds)

% Set bandpass filter
fnq = 1/(2*T.TR);				% Nyquist frequency
flp = 0.04;						% lowpass frequency of filter (Hz)
fhi = 0.07;						% highpass
Wn = [flp/fnq fhi/fnq];			% butterworth bandpass non-dimensional frequency
k = 2;							% 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn);	% construct the filter
clear fnq flp fhi Wn k

% Filter BOLD signal
for s = 1:sum(N.subjects)
	subjectBOLD(:,:,s) = detrend(demean(squeeze(subjectBOLD(:,:,s)))')';
	subjectBOLD(:,:,s) = filtfilt(bfilt,afilt, squeeze(subjectBOLD(:,:,s))')';
end

% z-score filtered BOLD signals
subjectZ = zscore(subjectBOLD, 0, 2);

% Indices for vectorizing lower triangle
Isubdiag = find(tril(ones(N.ROI),-1));

% Set hypothesis test parameters
pval.target = 0.05;

% Set number of neighbors to search for in KNN
co = HShannon_kNN_k_initialization(1);

% Set figure counter
N.fig = 1;

% Preallocate data arrays
TS.BOLD = cell(max(N.subjects), N.conditions);
TS.Z = cell(max(N.subjects), N.conditions);
FC = nan(N.ROI, N.ROI, max(N.subjects), N.conditions);

% Separate time series and FC matrix by subject & by condition
for c = 1:N.conditions
	ts = subjectBOLD(:,:,logical(I{:,c}));
	tz = subjectZ(:,:,logical(I{:,c}));
	for s = 1:N.subjects(c)
		TS.BOLD{s,c} = squeeze(ts(:,:,s)); TS.BOLD{s,c} = TS.BOLD{s,c}(:, 1:T.scan);
		TS.Z{s,c} = squeeze(tz(:,:,s)); TS.Z{s,c} = TS.Z{s,c}(:, 1:T.scan);
		FC(:,:,s,c) = corr(TS.BOLD{s, c}');
	end
end
clear c s ts tz subjectBOLD subjectZ



%% 2) Find significant components and activations

% Visualize variance explained per IC
F(N.fig) = figure;
N.fig = N.fig + 1;

% Extract number of assemblies using Marcenko-Pastur distribution
Z = cell(1,N.conditions);
N.IC = nan(N.conditions, 1);
explained = cell(N.conditions, 1);
tVar = nan(N.conditions, 1);
for c = 1:N.conditions
	Z{c} = cell2mat(TS.Z(:,c)');			% Get group-level z-score
	N.IC(c) = NumberofIC(Z{c});				% find number of components
	[~,~,~,~,explained{c},~] = pca(Z{c}');	% find percentage of explained variance per component
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

% Compute assembly memberships
memberships = cell(N.conditions, 1);
W = cell(N.conditions, 1);
ICs = cell(N.conditions, 1);
for c = 1:N.conditions
	disp(['Processing the ICs of group ', num2str(c), ' from BOLD data']);
	[~, memberships{c}, W{c}] = fastica(Z{c}, 'numOfIC', N.IC(c), 'verbose','off');
	
	% Normalize membership weights
	for i = 1:N.IC(c)
		memberships{c}(:,i) = memberships{c}(:,i)./norm(memberships{c}(:,i));
	end
	
	% Compute IC matrices
	ICs{c} = nan(N.ROI, N.ROI, N.IC(c));
	for i = 1:N.IC(c)
		ICs{c}(:,:,i) = memberships{c}(:,i) * memberships{c}(:,i)';
	end
end

% Compute number of pairings between ICs and time series
pairings = vertcat(repmat([1:N.conditions]', 1,2), nchoosek(1:N.conditions,2), flip(nchoosek(1:N.conditions,2),2));
N.pairings = size(pairings, 1);

% Compute IC time courses, both original and projected
activities.cond = cell(N.conditions, N.pairings-N.conditions);
activities.subj = cell(max(N.subjects), N.conditions, N.pairings-N.conditions);
for p = 1:N.pairings
	activities.cond{pairings(p,1), pairings(p,2)} = W{pairings(p,1)}*Z{pairings(p,2)};
	for s = 1:N.subjects(pairings(p,2))
		activities.subj{s,pairings(p,1),pairings(p,2)} = W{pairings(p,1)}*TS.Z{s,pairings(p,2)};
	end
end
clear I k t c s i p


%% Compute IC metrics

% Store results for IC activations:
%	Component-wise Kuramoto order parameter & metastability
%	Subject-wise entropies
activities.entro = nan(max(N.subjects), N.conditions, N.pairings-N.conditions);
activities.metastable.cond = nan(N.conditions, N.pairings-N.conditions);
activities.kuramoto.cond = nan(N.conditions, N.pairings-N.conditions, T.scan*max(N.subjects));
activities.metastable.subj = nan(max(N.subjects), N.conditions, N.pairings-N.conditions);
activities.kuramoto.subj = nan(max(N.subjects), N.conditions, N.pairings-N.conditions, T.scan);

for p = 1:N.pairings
	entro = nan(N.IC(pairings(p,1)), N.subjects(pairings(p,2)));
	[activities.kuramoto.cond(pairings(p,1),pairings(p,2), 1:T.scan*N.subjects(pairings(p,2))), activities.metastable.cond(pairings(p,1),pairings(p,2))] = findStability(activities.cond{pairings(p,1),pairings(p,2)});
	for s = 1:N.subjects(pairings(p,2))
		[activities.kuramoto.subj(s,pairings(p,1),pairings(p,2),:), activities.metastable.subj(s,pairings(p,1),pairings(p,2))] = findStability(activities.subj{s,pairings(p,1),pairings(p,2)});
		for i = 1:N.IC(pairings(p,1))
			entro(i, s) = HShannon_kNN_k_estimation(activities.subj{s,pairings(p,1),pairings(p,2)}(i,:), co);
		end
	end
	activities.entro(1:N.subjects(pairings(p,2)), pairings(p,1),pairings(p,2)) = squeeze(sum(entro, 1))';
end
clear p s i entro



%% Compare IC metrics between conditions, vs. permuted null distribution

% Define test types
ttype = {'kstest2', 'permutation'};

% Compare activations between conditions.  The second index c indicates which IC space is used
for t = 1:numel(ttype)
	disp(['Running ', ttype{t}, ' test on activations in AAL space.']);
	sig.AAL(t) = robustTests(Z{1}, Z{2}, N.ROI, 'pval',pval.target, 'testtype',ttype{t});						% Compare ROI time series
	for c = 1:N.conditions
		disp(['Running ', ttype{t}, ' test on condition ', num2str(c), ' IC space activations.']);
		sig.IC(t,c) = robustTests(activities.cond{c,1}, activities.cond{c,2}, N.IC(c), 'pval',pval.target, 'testtype',ttype{t});	% Compare IC time series
	end
end
clear c t

% Test activation metrics 
for c = 1:N.conditions
	% Subject metastabilities
	disp(['Running permutation test on condition ', num2str(c), ' IC space metastability.']);
	con = activities.metastable.subj(:,c,1); con = con(isfinite(con));
	pat = activities.metastable.subj(:,c,2); pat = pat(isfinite(pat));
	[sig.metastable.p(c), ~, sig.metastable.effsize(c)] = permutationTest(con, pat, 10000);
	
	% Subject entropies
	disp(['Running permutation test on condition ', num2str(c), ' IC space entropy.']);
	con = activities.entro(:,c,1); con = con(isfinite(con));
	pat = activities.entro(:,c,2); pat = pat(isfinite(pat));
	[sig.entro.p(c), ~, sig.entro.effsize(c)] = permutationTest(con, pat, 10000);
end
clear c con pat



%% Save results

% Save figures
if exist('F', 'var')
	savefig(F, fullfile(path{6}, fileName), 'compact');
	clear F
end

% Save all data
save(fullfile(path{6}, fileName));
