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
if strcmpi(compressType, 'eigenvector')
	fileName = 'LEICA90';
elseif strcmpi(compressType, 'average')
	fileName = 'MICA90';
elseif strcmpi(compressType, 'none')
	fileName = 'ICA90';
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

% Save BOLD data, indices
save(fullfile(path{6},fileName));



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

% Save interim results
save(fullfile(path{6}, fileName));

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
if strcmpi(compressType, 'LEICA') ||  strcmpi(compressType, 'Eigenvector')	% leading eigenvector
	dFC.concat = zeros(T.scan*sum(N.subjects), N.ROI);
elseif strcmpi(compressType, 'average')	% average activity
	dFC.concat = zeros(T.scan*sum(N.subjects), N.ROI);
else	% vectorized synchrony data
	dFC.concat = zeros(T.scan*sum(N.subjects), length(Isubdiag));
end

% Preallocate variables to save FC patterns and associated information
T.index = zeros(2, T.scan*sum(N.subjects));	% vector with subject nr and task at each t
t_all = 0;									% Index of time (starts at 0, updated until N.subjects*T.scan)

% Compute instantaneous FC (BOLD Phase Synchrony) and leading eigenvector (V1) for each time point
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		for t = 1:T.scan
			
			% Update time
			t_all = t_all+1;
			
			% Compute instantaneous phase synchrony FC
			iPH = zeros(N.ROI, N.ROI);
			if strcmpi(phaseType, 'cosine')
				for n = 1:N.ROI
					for m = 1:N.ROI
						iPH(n,m) = cos(TS.PH{s,c}(n,t) - TS.PH{s,c}(m,t));
					end
				end
			elseif strcmpi(phaseType, 'exponential')
				for n = 1:N.ROI
					for m = 1:N.ROI
						iPH(n,m) = exp(-3*adif(TS.PH{s,c}(n,t), TS.PH{s,c}(m,t))); 
					end
				end
			end
			
			% Extract leading eigenvector
			if strcmpi(compressType, 'eigenvector')
				[dFC.concat(t_all,:), ~] = eigs(iPH,1);
			elseif strcmpi(compressType, 'average')
				dFC.concat(t_all,:) = mean(iPH);
			else
				dFC.concat(t_all,:) = iPH(Isubdiag);
			end
			T.index(:,t_all) = [s c]';	% Information that at t_all, V1 corresponds to subject s in a given task
		end
	end
end
dFC.concat = dFC.concat';

% Clear dud variables
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


% Compute assembly activity timecourses and memberships
activities.cond = cell(N.conditions, 1);
activities.subj = cell(max(N.subjects), N.conditions);
memberships = cell(N.conditions, 1);
W = cell(N.conditions, 1);
ICs = cell(N.conditions, 1);
for c = 1:N.conditions
	disp(['Processing the ICs of group ', num2str(c), ' from BOLD data']);
	[activities.cond{c}, memberships{c}, W{c}] = fastica(dFC.cond{c}, 'numOfIC', N.IC(c), 'verbose','off');
	
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
	
	% Extract subject activation sequences
	for s = 1:N.subjects(c)
		I = (T.index(2,:) == c & T.index(1,:) == s);
		activities.subj{s,c} = activities.cond{c}(:,I);
	end
end

% Project timeseries onto alternative components
N.permutations = nchoosek(N.conditions, 1);
projections.cond = cell(N.permutations, 1);
projections.subj = cell(max(N.subjects), N.permutations);
projections.cond{1} =  W{1}*dFC.cond{2};
projections.cond{2} =  W{2}*dFC.cond{1};
for c = 1:N.permutations
	for s = 1:N.subjects(c)
		I = (T.index(2,:) == c);
		I = (I(1,:) == s);
		projections.subj{s,c} = projections.cond{c}(:,I);
	end
end



%% Compute IC metrics

% Store results for IC activations:
%	Activation magnitude means, medians, standard deviations
%	Component-wise Kuramoto order parameter & metastability
%	Subject-wise entropies
activities.av.cond = nan(N.assemblies, N.conditions);
activities.md.cond = nan(N.assemblies, N.conditions);
activities.sd.cond = nan(N.assemblies, N.conditions);
activities.metastable.cond = nan(N.conditions, 1);
activities.kuramoto.cond = nan(N.conditions, T.scan*max(N.subjects));
activities.av.subj = nan(N.assemblies, max(N.subjects), N.conditions);
activities.md.subj = nan(N.assemblies, max(N.subjects), N.conditions);
activities.sd.subj = nan(N.assemblies, max(N.subjects), N.conditions);
activities.metastable.subj = nan(max(N.subjects), N.conditions);
activities.kuramoto.subj = nan(max(N.subjects), N.conditions, T.scan);
activities.entro.subj = nan(N.assemblies, max(N.subjects), N.conditions);
for c = 1:N.conditions
	activities.av.cond(:,c) = mean(activities.cond{c}, 2, 'omitnan');
	activities.md.cond(:,c) = median(activities.cond{c}, 2, 'omitnan');
	activities.sd.cond(:,c) = std(activities.cond{c}, 0, 2, 'omitnan');
	[kuramoto.cond(c, 1:T.scan*N.subjects(c)), metastable.cond(c)] = findStability(activities.cond{c});
	for s = 1:N.subjects(c)
		activities.av.subj(:,s,c) = mean(activities.subj{s,c}, 2, 'omitnan');
		activities.md.subj(:,s,c) = median(activities.subj{s,c}, 2, 'omitnan');
		activities.sd.subj(:,s,c) = std(activities.subj{s,c}, 0, 2, 'omitnan');
		[activities.kuramoto.subj(s,c,:), activities.metastable.subj(s,c)] = findStability(activities.subj{s,c});
		for ass = 1:N.assemblies
			activities.entro.subj(ass, s, c) = HShannon_kNN_k_estimation(activities.subj{s,c}(ass,:), co);
		end
	end
end
clear c s
activities.entro.subj = squeeze(sum(entro.subj, 1));

% Convert metrics to table format
activities.av.cond = array2table(activities.av.cond, 'VariableNames',labels.Properties.VariableNames);
activities.md.cond = array2table(activities.md.cond, 'VariableNames',labels.Properties.VariableNames);
activities.sd.cond = array2table(activities.sd.cond, 'VariableNames',labels.Properties.VariableNames);
activities.metastable.cond = array2table(activities.metastable.cond, 'VariableNames', labels.Properties.VariableNames);
activities.metastable.subj = array2table(activities.metastable.subj, 'VariableNames', labels.Properties.VariableNames);
activities.entro.subj = array2table(activities.entro.subj, 'VariableNames', labels.Properties.VariableNames);



%% Compute projection metrics

% Store results for switched projections:
%	Activation magnitude means, medians, standard deviations
%	Component-wise Kuramoto order parameter & metastability
%	Subject-wise entropies
projections.av.cond = nan(N.assemblies, N.conditions);
projections.md.cond = nan(N.assemblies, N.conditions);
projections.sd.cond = nan(N.assemblies, N.conditions);
projections.metastable.cond = nan(N.conditions, 1);
projections.kuramoto.cond = nan(N.conditions, T.scan*max(N.subjects));
projections.av.subj = nan(N.assemblies, max(N.subjects), N.conditions);
projections.md.subj = nan(N.assemblies, max(N.subjects), N.conditions);
projections.sd.subj = nan(N.assemblies, max(N.subjects), N.conditions);
projections.metastable.subj = nan(max(N.subjects), N.conditions);
projections.kuramoto.subj = nan(max(N.subjects), N.conditions, T.scan);
projections.entro.subj = nan(N.assemblies, max(N.subjects), N.conditions);
for c = 1:N.conditions
	projections.av.cond(:,c) = mean(projections.cond{c}, 2, 'omitnan');
	projections.md.cond(:,c) = median(projections.cond{c}, 2, 'omitnan');
	projections.sd.cond(:,c) = std(projections.cond{c}, 0, 2, 'omitnan');
	[projections.kuramoto.cond(c, 1:T.scan*N.subjects(c)), projections.metastable.cond(c)] = findStability(projections.cond{c});
	for s = 1:N.subjects(c)
		projections.av.subj(:,s,c) = mean(projections.subj{s,c}, 2, 'omitnan');
		projections.md.subj(:,s,c) = median(projections.subj{s,c}, 2, 'omitnan');
		projections.sd.subj(:,s,c) = std(projections.subj{s,c}, 0, 2, 'omitnan');
		[projections.kuramoto.subj(s,c,:), projections.metastable.subj(s,c)] = findStability(projections.subj{s,c});
		for ass = 1:N.assemblies
			projections.entro.subj(ass, s, c) = HShannon_kNN_k_estimation(projections.subj{s,c}(ass,:), co);
		end
	end
end
clear c s
projections.entro.subj = squeeze(sum(entro.subj, 1));

% Convert metrics to table format
projections.av.cond = array2table(projections.av.cond, 'VariableNames',labels.Properties.VariableNames);
projections.md.cond = array2table(projections.md.cond, 'VariableNames',labels.Properties.VariableNames);
projections.sd.cond = array2table(projections.sd.cond, 'VariableNames',labels.Properties.VariableNames);
projections.metastable.cond = array2table(projections.metastable.cond, 'VariableNames', labels.Properties.VariableNames);
projections.metastable.subj = array2table(projections.metastable.subj, 'VariableNames', labels.Properties.VariableNames);
projections.entro.subj = array2table(projections.entro.subj, 'VariableNames', labels.Properties.VariableNames);



%% Compare IC metrics between conditions, vs. permuted null distribution

% Define test types
ttype = {'kstest2', 'permutation'};

% Test with Kolmogorov-Smirnov, permutation test
for t = 1:numel(ttype)
	% Compare activations between conditions
	sig.AAL.TS{t} = robustTests(dFC.cond{1}, dFC.cond{2}, N.ROI, 'pval',pval.target, 'testtype',ttype{t});						% Compare ROI time series
	sig.act.IC.TS{t} = robustTests(activities.cond{1}, activities.cond{2}, N.assemblies, 'pval',pval.target, 'testtype',ttype{t});	% Compare IC time series
	
	% Average activation magnitude(s)
	con = activities.av.subj(:,:,1); con = con(isfinite(con));
	pat = activities.av.subj(:,:,2); pat = pat(isfinite(pat));
	sig.act.av{t} = robustTests(con, pat, N.assemblies, 'pval',pval.target, 'testtype',ttype{t});

	% Activation medians(s)
	con = activities.md.subj(:,:,1); con = con(isfinite(con));
	pat = activities.md.subj(:,:,2); pat = pat(isfinite(pat));
	sig.act.md = robustTests(con, pat, N.assemblies, 'pval',pval.target, 'testtype',ttype{t});

	% Activation standard deviations(s)
	con = activities.sd.subj(:,:,1); con = con(isfinite(con));
	pat = activities.sd.subj(:,:,2); pat = pat(isfinite(pat));
	sig.act.sd = robustTests(con, pat, N.assemblies, 'pval',pval.target, 'testtype',ttype{t});
end
clear con pat t

% Subject metastabilities
con = activities.metastable.subj{:,'Control'}(isfinite(activities.metastable.subj{:,'Control'}));
pat = activities.metastable.subj{:,'Patient'}(isfinite(activities.metastable.subj{:,'Patient'}));
sig.act.metastable.p = permutationTest(con, pat, 'exact',1, 'sidedness','both');
clear con pat

% Subject entropies
con = activities.entro.subj{:,'Control'}(isfinite(activities.entro.subj{:,'Control'}));
pat = activities.entro.subj{:,'Patient'}(isfinite(activities.entro.subj{:,'Patient'}));
sig.act.entro.p = permutationTest(con, pat, 'exact',1, 'sidedness','both');
clear con pat



%% Compare IC metrics between conditions, vs. permuted null distribution

% Define test types
ttype = {'kstest2', 'permutation'};

% Test with Kolmogorov-Smirnov, permutation test
for t = 1:numel(ttype)
	% Compare activations between conditions
	sig.proj.AAL.TS{t} = robustTests(dFC.cond{1}, dFC.cond{2}, N.ROI, 'pval',pval.target, 'testtype',ttype{t});						% Compare ROI time series
	sig.proj.IC.TS{t} = robustTests(projections.cond{1}, projections.cond{2}, N.assemblies, 'pval',pval.target, 'testtype',ttype{t});	% Compare IC time series
	
	% Average activation magnitude(s)
	con = projections.av.subj(:,:,1); con = con(isfinite(con));
	pat = projections.av.subj(:,:,2); pat = pat(isfinite(pat));
	sig.proj.av{t} = robustTests(con, pat, N.assemblies, 'pval',pval.target, 'testtype',ttype{t});

	% Activation medians(s)
	con = projections.md.subj(:,:,1); con = con(isfinite(con));
	pat = projections.md.subj(:,:,2); pat = pat(isfinite(pat));
	sig.proj.md = robustTests(con, pat, N.assemblies, 'pval',pval.target, 'testtype',ttype{t});

	% Activation standard deviations(s)
	con = projections.sd.subj(:,:,1); con = con(isfinite(con));
	pat = projections.sd.subj(:,:,2); pat = pat(isfinite(pat));
	sig.proj.sd = robustTests(con, pat, N.assemblies, 'pval',pval.target, 'testtype',ttype{t});
end
clear con pat t

% Subject metastabilities
con = projections.metastable.subj{:,'Control'}(isfinite(projections.metastable.subj{:,'Control'}));
pat = projections.metastable.subj{:,'Patient'}(isfinite(projections.metastable.subj{:,'Patient'}));
sig.proj.metastable.p = permutationTest(con, pat, 'exact',1, 'sidedness','both');
clear con pat

% Subject entropies
con = projections.entro.subj{:,'Control'}(isfinite(projections.entro.subj{:,'Control'}));
pat = projections.entro.subj{:,'Patient'}(isfinite(projections.entro.subj{:,'Patient'}));
sig.proj.entro.p = permutationTest(con, pat, 'exact',1, 'sidedness','both');
clear con pat



%% Save results

% Save figures
if exist('F', 'var')
	savefig(F, fullfile(path{6}, fileName), 'compact');
	clear F
end

% Save all data
save(fullfile(path{6}, fileName));
