%% LEICA Extraction of Network States
%	This script extracts and transforms neuroimaging data for ICA-based
% comparisons of the network states.  The current script only computes
% phase-based synchrony measures.  It is possible to compute synchrony
% measures based on the z-score and the BOLD signal, but at the moment this
% will only add complexity.


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
clear fpath k

% Load structural data
load(fullfile(path{4}, 'sc90.mat'));

% Set methods
distType = 'cosine';	% for measuring phase: cosine or exponential
compressType = 'eigenvector';	% for compressing matrix: eigenvector, average, or none

% File to save
fileName = 'ZICA90_CIC';
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



%% Compute ICs from BOLD signal

% Concatenate z-scored signals
Z = cell(1,N.conditions);
for c = 1:N.conditions
	Z{c} = cell2mat(TS.Z(:,c)');
end
Z = cell2mat(Z);

% Extract number of assemblies using Marcenko-Pastur distribution
N.IC = NumberofIC(Z);

% Use PCA to estimate amount of variance captured with N.IC
[~,~,~,~,explained,~] = pca(Z');
explainedVar = sum(explained(N.IC+1:end));

% Compute assembly activity timecourses and memberships
disp('Processing the ICs from BOLD data');
[activities.concat, memberships, W] = fastica(Z, 'numOfIC', N.IC, 'verbose','off');

% Separate assembly activations by condition & subject
Z = cell(1, N.conditions);
activities.cond = cell(1, N.conditions);
activities.subj = cell(max(N.subjects), N.conditions);
for c = 1:N.conditions
	Z{c} = cell2mat(TS.Z(:,c)');
	activities.cond{c} = W*Z{c};
	for s = 1:N.subjects(c)
		activities.subj{s,c} = W*TS.Z{s,c};
	end
end
clear I s c

% Normalize membership weights
for k = 1:N.IC
	memberships(:,k) = memberships(:,k)./norm(memberships(:,k));
end
clear k

% Compute LEICA assembly matrices
ICs = nan(N.ROI, N.ROI, N.IC);
for i = 1:size(memberships,2)
	ICs(:,:,i) = memberships(:,i) * memberships(:,i)';
end
clear i

% Compute IC frequency spectra
L = sum(N.subjects)*T.scan;
ICspect.concat = abs(fft(activities.concat')/L)';
ICspect.concat = ICspect.concat(:, 1:round(L/2+1));
ICspect.concat(:, 2:end-1) = 2*ICspect.concat(:, 2:end-1);
for c = 1:N.conditions
	L = N.subjects(c)*T.scan;
	ICspect.cond{c} = abs(fft(activities.cond{c}')/L)';
	ICspect.cond{c} = ICspect.cond{c}(:, 1:round(L/2+1));
	ICspect.cond{c}(:, 2:end-1) = 2*ICspect.cond{c}(:, 2:end-1);
	for s = 1:N.subjects(c)
		ICspect.subj{s,c} = fft(activities.subj{s,c}')';
		ICspect.subj{s,c} = ICspect.subj{s,c}(:, 1:round(T.scan/2+1));
		ICspect.subj{s,c}(:, 2:end-1) = 2*ICspect.subj{s,c}(:, 2:end-1);
	end
end
clear c s L


%% Compute IC metrics

% Preallocate storage arrays:
%	Activation magnitude means, medians, standard deviations
%	Component-wise Kuramoto order parameter & metastability
%	Subject-wise entropies
activities.av.cond = nan(N.IC, N.conditions);
activities.md.cond = nan(N.IC, N.conditions);
activities.sd.cond = nan(N.IC, N.conditions);
metastable.cond = nan(1, N.conditions);
kuramoto.cond = nan(N.conditions, T.scan*max(N.subjects));
activities.av.subj = nan(N.IC, max(N.subjects), N.conditions);
activities.md.subj = nan(N.IC, max(N.subjects), N.conditions);
activities.sd.subj = nan(N.IC, max(N.subjects), N.conditions);
metastable.subj = nan(max(N.subjects), N.conditions);
kuramoto.subj = nan(max(N.subjects), N.conditions, T.scan);
entro.subj = nan(N.IC, max(N.subjects), N.conditions);
for c = 1:N.conditions
	activities.av.cond(:,c) = mean(activities.cond{c}, 2, 'omitnan');
	activities.md.cond(:,c) = median(activities.cond{c}, 2, 'omitnan');
	activities.sd.cond(:,c) = std(activities.cond{c}, 0, 2, 'omitnan');
	[kuramoto.cond(c, 1:T.scan*N.subjects(c)), metastable.cond(c)] = findStability(activities.cond{c});
	for s = 1:N.subjects(c)
		activities.av.subj(:,s,c) = mean(activities.subj{s,c}, 2, 'omitnan');
		activities.md.subj(:,s,c) = median(activities.subj{s,c}, 2, 'omitnan');
		activities.sd.subj(:,s,c) = std(activities.subj{s,c}, 0, 2, 'omitnan');
		[kuramoto.subj(s,c,:), metastable.subj(s,c)] = findStability(activities.subj{s,c});
		for ass = 1:N.IC
			entro.subj(ass, s, c) = HShannon_kNN_k_estimation(activities.subj{s,c}(ass,:), co);
		end
	end
end
clear c s
entro.subj = squeeze(sum(entro.subj, 1));

% Convert metrics to table format
activities.av.cond = array2table(activities.av.cond, 'VariableNames',labels.Properties.VariableNames);
activities.md.cond = array2table(activities.md.cond, 'VariableNames',labels.Properties.VariableNames);
activities.sd.cond = array2table(activities.sd.cond, 'VariableNames',labels.Properties.VariableNames);
metastable.cond = array2table(metastable.cond, 'VariableNames', labels.Properties.VariableNames);
metastable.subj = array2table(metastable.subj, 'VariableNames', labels.Properties.VariableNames);
entro.subj = array2table(entro.subj, 'VariableNames', labels.Properties.VariableNames);



%% Compare IC metrics between conditions, vs. permuted null distribution

% Define test types
ttype = {'kstest2', 'ranksum', 'permutation'};

% Test with Kolmogorov-Smirnov, permutation test
for t = 1:numel(ttype)
	disp(['Running ', ttype{t}, ' tests on activations.']);
	
	% Compare activations between conditions
	sig.AAL{t} = robustTests(Z{1}, Z{2}, N.ROI, 'p',pval.target, 'testtype',ttype{t});						% Compare ROI time series
	sig.IC{t} = robustTests(activities.cond{1}, activities.cond{2}, N.IC, 'p',pval.target, 'testtype',ttype{t});	% Compare IC time series
	
	% Average activation magnitude(s)
	con = activities.av.subj(:,:,1); con = con(isfinite(con));
	pat = activities.av.subj(:,:,2); pat = pat(isfinite(pat));
	sig.av{t} = robustTests(con, pat, N.IC, 'p',pval.target, 'testtype',ttype{t});

	% Activation medians(s)
	con = activities.md.subj(:,:,1); con = con(isfinite(con));
	pat = activities.md.subj(:,:,2); pat = pat(isfinite(pat));
	sig.md{t} = robustTests(con, pat, N.IC, 'p',pval.target, 'testtype',ttype{t});

	% Activation standard deviations(s)
	con = activities.sd.subj(:,:,1); con = con(isfinite(con));
	pat = activities.sd.subj(:,:,2); pat = pat(isfinite(pat));
	sig.sd{t} = robustTests(con, pat, N.IC, 'p',pval.target, 'testtype',ttype{t});
end
clear con pat t Z

% Subject metastabilities
disp('Running permutation tests on metastability.');
con = metastable.subj{:,'Control'}(isfinite(metastable.subj{:,'Control'}));
pat = metastable.subj{:,'Patient'}(isfinite(metastable.subj{:,'Patient'}));
sig.metastable.p = permutationTest(con, pat, 10000, 'sidedness','both');

% Subject entropies
disp('Running permutation tests on entropy.');
con = entro.subj{:,'Control'}(isfinite(entro.subj{:,'Control'}));
pat = entro.subj{:,'Patient'}(isfinite(entro.subj{:,'Patient'}));
sig.entro.p = permutationTest(con, pat, 10000, 'sidedness','both');



%% Save results

% Save figures
if exist('F', 'var')
	savefig(F, fullfile(path{6}, fileName), 'compact');
	clear F
end

% Save all data
save(fullfile(path{6}, fileName));

