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
distType = 'cosine';	% for measuring phase: cosine or exponential
compressType = 'eigenvector';	% for compressing matrix: eigenvector, average, or none

% File to save
switch compressType
	case {'LEICA', 'eigenvector'}
		fileName = 'LEICA90_CIC';
	case 'average'
		fileName = 'MICA90_CIC';
	otherwise
		fileName = 'ICA90_CIC';
end
switch distType
	case 'cosine'
		fileName = strcat(fileName, '_COS');
	case 'exponential'
		fileName = strcat(fileName, '_EXP');
end
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
FC = nan(N.ROI, N.ROI, max(N.subjects), N.conditions);

% Separate time series and FC matrix by subject & by condition
for c = 1:N.conditions
	ts = subjectBOLD(:,:,logical(I{:,c}));
	for s = 1:N.subjects(c)
		TS.BOLD{s,c} = squeeze(ts(:,:,s));
		TS.BOLD{s,c} = TS.BOLD{s,c}(:, 1:T.scan);
		FC(:,:,s,c) = corr(TS.BOLD{s, c}');
	end
end
clear c s ts subjectBOLD

% Preallocate storage arrays
TS.PH = cell(max(N.subjects), N.conditions);

% Compute BOLD phase and z-score
disp('Computing phase of BOLD signal');
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		[TS.PH{s,c}, TS.BOLD{s,c}] = regionPhase(TS.BOLD{s,c}, bfilt, afilt);
	end
end
clear s c

% Preallocate storage arrays
switch compressType
	case {'LEICA', 'eigenvector'}
		dFC.concat = zeros(N.ROI, T.scan*sum(N.subjects));
	case 'average'
		dFC.concat = zeros(N.ROI, T.scan*sum(N.subjects));
	otherwise
		dFC.concat = zeros(length(Isubdiag), T.scan*sum(N.subjects));
end

% Preallocate variables to save FC patterns and associated information
T.index = zeros(2, T.scan*sum(N.subjects));	% vector with subject nr and task at each t
t = 0;

% Compute instantaneous FC (BOLD Phase Synchrony) and leading eigenvector (V1) for each time point
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		
		% Index subject, condition of current dFC sequence
		t = t+1 : t+T.scan;
		T.index(:,t) = repmat([s c]', 1,T.scan);
		
		% Extract dFC
		dFC.concat(:,t) = LEdFC(TS.PH{s,c}, 'distType',distType, 'compressType',compressType, 'nROI',N.ROI, 'T',T.scan);
		t = t(end);
	end
end
clear afilt bfilt c s t m n iPH iZ V1 t_all Isubdiag sc90

% Preallocate storage arrays for condition-wise dFC, subject-wise dFC
dFC.cond = cell(1, N.conditions);
dFC.subj = cell(max(N.subjects), N.conditions);

% Segment dFC
for c = 1:N.conditions
	I = T.index(2,:) == c;
	dFC.cond{c} = dFC.concat(:,I);
	for s = 1:N.subjects(c)
		I = (T.index(2,:) == c & T.index(1,:) == s);
		dFC.subj{s,c} = dFC.concat(:,I);
	end
end
clear I s c

% Compute dFC frequency spectra
L = sum(N.subjects)*T.scan;
dFCspect.concat = fft(dFC.concat')';
dFCspect.concat = dFCspect.concat(:, 1:round(L/2+1));
dFCspect.concat(:, 2:end-1) = 2*dFCspect.concat(:, 2:end-1);
for c = 1:N.conditions
	L = N.subjects(c)*T.scan;
	dFCspect.cond{c} = fft(dFC.cond{c}')';
	dFCspect.cond{c} = dFCspect.cond{c}(:, 1:round(L/2+1));
	dFCspect.cond{c}(:, 2:end-1) = 2*dFCspect.cond{c}(:, 2:end-1);
	for s = 1:N.subjects(c)
		dFCspect.subj{s,c} = fft(dFC.subj{s,c}')';
		dFCspect.subj{s,c} = dFCspect.subj{s,c}(:, 1:round(T.scan/2+1));
		dFCspect.subj{s,c}(:, 2:end-1) = 2*dFCspect.subj{s,c}(:, 2:end-1);
	end
end
clear c s L



%% 2) Find significant components and activations

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

% Compute assembly memberships
memberships = cell(N.conditions, 1);
W = cell(N.conditions, 1);
ICs = cell(N.conditions, 1);
for c = 1:N.conditions
	disp(['Processing the ICs of group ', num2str(c), ' from BOLD data']);
	[~, memberships{c}, W{c}] = fastica(dFC.cond{c}, 'numOfIC', N.IC(c), 'verbose','off');
	
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
	activities.cond{pairings(p,1), pairings(p,2)} = W{pairings(p,1)}*dFC.cond{pairings(p,2)};
	for s = 1:N.subjects(pairings(p,2))
		t = T.index(1,(T.index(2,:) == pairings(p,2)));
		I = (t == s);
		activities.subj{s,pairings(p,1),pairings(p,2)} = activities.cond{pairings(p,1),pairings(p,2)}(:,I);
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
	sig.AAL(t) = robustTests(dFC.cond{1}, dFC.cond{2}, N.ROI, 'pval',pval.target, 'testtype',ttype{t});						% Compare ROI time series
	for c = 1:N.conditions
		sig.IC(t,c) = robustTests(activities.cond{c,1}, activities.cond{c,2}, N.IC(c), 'pval',pval.target, 'testtype',ttype{t});	% Compare IC time series
	end
end
clear c t

% Test activation metrics 
for c = 1:N.conditions
	% Subject metastabilities
	con = activities.metastable.subj(:,c,1); con = con(isfinite(con));
	pat = activities.metastable.subj(:,c,2); pat = pat(isfinite(pat));
	[sig.metastable(c).p, ~, sig.metastable(c).effsize] = permutationTest(con, pat, 10000);
	
	% Subject entropies
	con = activities.entro(:,c,1); con = con(isfinite(con));
	pat = activities.entro(:,c,2); pat = pat(isfinite(pat));
	[sig.entro(c).p, ~, sig.entro(c).effsize] = permutationTest(con, pat, 10000);
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
