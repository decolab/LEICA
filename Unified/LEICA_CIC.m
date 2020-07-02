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
distType = 'cosine';	% for measuring distance: cosine or exponential
compressType = 'none';	% for compressing matrix: eigenvector, average, or none

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
Isubdiag = find(triu(ones(N.ROI), 1));

% Set hypothesis test parameters
pval.target = 0.05;

% Set number of neighbors to search for in KNN
co = HShannon_kNN_k_initialization(1);

% Set figure counter
N.fig = 1;

% Preallocate data arrays

% Separate time series and FC matrix by subject & by condition
BOLD = cell(max(N.subjects), N.conditions);
FC = nan(N.ROI, N.ROI, max(N.subjects), N.conditions);
for c = 1:N.conditions
	ts = subjectBOLD(:,:,logical(I{:,c}));
	for s = 1:N.subjects(c)
		BOLD{s,c} = squeeze(ts(:,:,s)); BOLD{s,c} = BOLD{s,c}(:, 1:T.scan);
		FC(:,:,s,c) = corr(BOLD{s, c}');
	end
end
clear c s ts subjectBOLD

% Plot BOLD signals
F.BOLD = figure; hold on;
subplot(2,2,1); imagesc(cell2mat(BOLD(:,1)')); colorbar; title('Patient BOLD');
subplot(2,2,2); imagesc(cell2mat(BOLD(:,1)')); colorbar; title('Control BOLD');
subplot(2,2,[3 4]); hold on; histogram(cell2mat(BOLD(:,1)')); histogram(cell2mat(BOLD(:,2)')); legend('Patient', 'Control');

% Compute BOLD phase and z-score
PH = cell(max(N.subjects), N.conditions);
Z.subj = cell(max(N.subjects), N.conditions);
Z.cond = cell(1,N.conditions);
disp('Computing phase of BOLD signal');
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		[PH{s,c}, BOLD{s,c}] = regionPhase(BOLD{s,c}, bfilt, afilt);
		Z.subj{s,c} = zscore(BOLD{s,c},0,2);
	end
	Z.cond{c} = cell2mat(Z.subj(:,c)');
end
Z.concat = cell2mat(Z.cond);
clear s c

% Preallocate storage arrays
switch compressType
	case {'LEICA', 'eigenvector', 'average'}
		dFC.concat = zeros(N.ROI, T.scan*sum(N.subjects));
	otherwise
		dFC.concat = zeros(length(Isubdiag), T.scan*sum(N.subjects));
end

% Preallocate variables to save FC patterns and associated information
T.index = zeros(2, T.scan*sum(N.subjects));	% vector with subject nr and task at each t
t = 0;									% Index of time (starts at 0, updated until N.subjects*T.scan)

% Compute instantaneous FC (BOLD Phase Synchrony) and leading eigenvector (V1) for each time point
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		
		% Index subject, condition of current dFC sequence
		t = t+1 : t+T.scan;
		T.index(:,t) = repmat([s c]', 1,T.scan);
		
		% Extract dFC
		dFC.concat(:,t) = LEdFC(PH{s,c}, 'distType',distType, 'compressType',compressType, 'nROI',N.ROI, 'T',T.scan);
		t = t(end);
	end
end
clear afilt bfilt c s t m n iPH iZ V1 t_all Isubdiag sc90

% Segment dFC
dFC.cond = cell(1, N.conditions);
dFC.subj = cell(max(N.subjects), N.conditions);
for c = 1:N.conditions
	I = T.index(2,:) == c;
	dFC.cond{c} = dFC.concat(:,I);
	for s = 1:N.subjects(c)
		I = (T.index(2,:) == c & T.index(1,:) == s);
		dFC.subj{s,c} = dFC.concat(:,I);
	end
end
clear I s c

% Plot LEdFC signals
F.LEdFC = figure; hold on;
subplot(2,2,1); imagesc(dFC.cond{1}); colorbar; title('Patient LEdFC');
subplot(2,2,2); imagesc(dFC.cond{2}); colorbar; title('Control LEdFC');
subplot(2,2,[3 4]); hold on; histogram(dFC.cond{1}); histogram(dFC.cond{2}); legend('Patient', 'Control');


% Compute FCD, power spectra of dFC
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		FCD.dFC.subj{s,c} = computeFCD(dFC.subj{s,c}, 'cosine');
		pspect.dFC.subj{s,c} = pspectrum(dFC.subj{s,c}', 1/T.TR)';
	end
end

% Compute KS distances between control, patient groups
pat = reshape(cell2mat(FCD.dFC.subj(1:N.subjects(1),1)), [N.subjects(1)*175^2, 1]);
con = reshape(cell2mat(FCD.dFC.subj(1:N.subjects(2),2)), [N.subjects(2)*175^2, 1]);
[FCD.dFC.h, FCD.dFC.p, FCD.dFC.ksdist] = kstest2(con, pat);

% Visualize dFC FCD
F.FCD.dFC = figure; hold on;
subplot(2,2, 1); imagesc(cell2mat(FCD.dFC.subj(1:N.subjects(1),1))'); colorbar; title('Patient dFC FCD');
subplot(2,2, 2); imagesc(cell2mat(FCD.dFC.subj(1:N.subjects(2),2))'); colorbar; title('Control dFC FCD');
subplot(2,2, [3 4]); hold on; histogram(pat); histogram(con); legend({'Patient', 'Control'});
clear con pat



%% Compute ICs from dFC

% Extract number of assemblies using Marcenko-Pastur distribution
N.IC = NumberofIC(dFC.concat);

% Use PCA to estimate amount of variance captured with N.IC
[~,~,~,~,explained,~] = pca(dFC.concat');
explainedVar = sum(explained(N.IC+1:end));

% Compute assembly activity timecourses and memberships
disp('Processing the ICs from BOLD data');
[activities.concat, memberships, W] = fastica(dFC.concat, 'numOfIC', N.IC, 'verbose','off');

% Separate assembly activations by condition & subject
activities.cond = cell(1, N.conditions);
activities.subj = cell(max(N.subjects), N.conditions);
for c = 1:N.conditions
	I = T.index(2,:) == c;
	activities.cond{c} = activities.concat(:,I);
	for s = 1:N.subjects(c)
		I = (T.index(2,:) == c & T.index(1,:) == s);
		activities.subj{s,c} = activities.concat(:,I);
	end
end
clear I s c

% Normalize membership weights
for k = 1:N.IC
	memberships(:,k) = memberships(:,k)./norm(memberships(:,k));
end
clear k

% Compute component matrices
if strcmpi(compressType, {'LEICA', 'eigenvector', 'average'})
	ICs = nan(N.ROI, N.ROI, size(memberships,2));
	for i = 1:size(memberships,2)
		ICs(:,:,i) = memberships(:,i) * memberships(:,i)';
	end
end
clear i

% Visualize IC activations
F.LEICA = figure; hold on;
subplot(2,2, 1); imagesc(cell2mat(activities.subj(1:N.subjects(1),1)')); colorbar; title('Patient LEICA Activations');
subplot(2,2, 2); imagesc(cell2mat(activities.subj(1:N.subjects(2),2)')); colorbar; title('Control LEICA Activations');
subplot(2,2, [3 4]); hold on; histogram(cell2mat(activities.subj(1:N.subjects(1),1))); histogram(cell2mat(activities.subj(1:N.subjects(2),2))); legend({'Patient', 'Control'});



%% Compute IC metrics

% Preallocate storage arrays:
%	Activation magnitude means, medians, standard deviations
%	Component-wise Kuramoto order parameter & metastability
%	Subject-wise entropies
metastable.cond = nan(1, N.conditions);
kuramoto.cond = nan(N.conditions, T.scan*max(N.subjects));
metastable.subj = nan(max(N.subjects), N.conditions);
kuramoto.subj = nan(max(N.subjects), N.conditions, T.scan);
entro.subj = nan(N.IC, max(N.subjects), N.conditions);
for c = 1:N.conditions
	[kuramoto.cond(c, 1:T.scan*N.subjects(c)), metastable.cond(c)] = findStability(activities.cond{c});
	for s = 1:N.subjects(c)
		[kuramoto.subj(s,c,:), metastable.subj(s,c)] = findStability(activities.subj{s,c});
		for ass = 1:N.IC
			entro.subj(ass, s, c) = HShannon_kNN_k_estimation(activities.subj{s,c}(ass,:), co);
		end
	end
end
clear c s ass
entro.subj = squeeze(sum(entro.subj, 1));

% Convert metrics to table format
metastable.cond = array2table(metastable.cond, 'VariableNames', labels.Properties.VariableNames);
metastable.subj = array2table(metastable.subj, 'VariableNames', labels.Properties.VariableNames);
entro.subj = array2table(entro.subj, 'VariableNames', labels.Properties.VariableNames);

% Visualize IC metrics
F.metrics = figure; hold on;
subplot(1,2,1); hold on; histogram(metastable.subj{:,'Patient'}); histogram(metastable.subj{:,'Control'}); legend('Patients', 'Controls'); title('Metastability');
subplot(1,2,2); hold on; histogram(entro.subj{:,'Patient'}); histogram(entro.subj{:,'Control'}); legend('Patients', 'Controls'); title('Entropy');



%% Compute FCD, power spectra, and periodograms of ICs

% Compute FCD, power spectra, periodogram of subject-level ICs
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		FCD.IC.subj{s,c} = computeFCD(activities.subj{s,c}, 'cosine');
		pspect.IC.subj{s,c} = pspectrum(activities.subj{s,c}', 1/T.TR)';
		pgram.IC.subj{s,c} = periodogram(activities.subj{s,c}', [], [], 1/T.TR)';
	end
end
clear c s i

% Visualize FCD
F.FCD.IC = figure; hold on;
subplot(2,2,1); imagesc(cell2mat(FCD.IC.subj(:,1))'); colorbar; title('LEICA FCD: Patients');
subplot(2,2,2); imagesc(cell2mat(FCD.IC.subj(:,2))'); colorbar; title('LEICA FCD: Controls');
subplot(2,2,[3 4]); hold on; histogram(cell2mat(FCD.IC.subj(:,1))'); histogram(cell2mat(FCD.IC.subj(:,2))'); legend({'Patient', 'Control'});

% Compute KS distances between control, patient FCD
pat = reshape(cell2mat(FCD.IC.subj(1:N.subjects(1),1)), [N.subjects(1)*175^2, 1]);
con = reshape(cell2mat(FCD.IC.subj(1:N.subjects(2),2)), [N.subjects(2)*175^2, 1]);
[FCD.IC.h, FCD.IC.p, FCD.IC.ksdist] = kstest2(con, pat);
clear con pat

% Compute Euclidean distances between subject power spectra, periodograms
dummy.spect = nan(N.IC, size(pspect.IC.subj{1,1},2), sum(N.subjects));
dummy.gram = nan(N.IC, size(pgram.IC.subj{1,1},2), sum(N.subjects));
pspect.IC.dist = nan(sum(N.subjects), sum(N.subjects), N.IC);
pgram.IC.dist = nan(sum(N.subjects), sum(N.subjects), N.IC);
ind = 0;
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		ind = ind+1;
		dummy.spect(:,:,ind) = pspect.IC.subj{s,c};
		dummy.gram(:,:,ind) = pgram.IC.subj{s,c};
	end
end
for i = 1:N.IC
	pspect.IC.dist(:,:,i) = squareform(pdist(squeeze(dummy.spect(i,:,:))'));
	pgram.IC.dist(:,:,i) = squareform(pdist(squeeze(dummy.gram(i,:,:))'));
end
clear c s i ind dummy

% Compute mean, standard deviation of inter-subject distance
pspect.IC.ave = mean(pspect.IC.dist, 3);
pspect.IC.std = std(pspect.IC.dist, [], 3);
pgram.IC.ave = mean(pgram.IC.dist, 3);
pgram.IC.std = std(pgram.IC.dist, [], 3);

% Visualize mean, standard deviations of spectral and periodogram distances
figure; hold on;
subplot(1,2,1); imagesc(pspect.IC.ave); colorbar; title('Average Spectral Distance between Subjects')
subplot(1,2,2); imagesc(pspect.IC.std); colorbar; title('Standard Deviation in Spectral Distance between Subjects')
figure; hold on;
subplot(1,2,1); imagesc(pgram.IC.ave); colorbar; title('Average Periodogram Distance between Subjects')
subplot(1,2,2); imagesc(pgram.IC.std); colorbar; title('Standard Deviation in Periodogram Distance between Subjects')

% Separate into classes (control, patient, inter)
pspect.IC.sect{1} = pspect.IC.dist(1:N.subjects(1), 1:N.subjects(1));
pspect.IC.sect{2} = pspect.IC.dist(1:N.subjects(1), N.subjects(1)+1:sum(N.subjects));
pspect.IC.sect{3} = pspect.IC.dist(N.subjects(1)+1:sum(N.subjects), N.subjects(1)+1:sum(N.subjects));
pgram.IC.sect{1} = pgram.IC.dist(1:N.subjects(1), 1:N.subjects(1));
pgram.IC.sect{2} = pgram.IC.dist(1:N.subjects(1), N.subjects(1)+1:sum(N.subjects));
pgram.IC.sect{3} = pgram.IC.dist(N.subjects(1)+1:sum(N.subjects), N.subjects(1)+1:sum(N.subjects));

% Compute KS distances between control, patient power spectra
pat = reshape(pspect.IC.sect{1}, [1, numel(pspect.IC.sect{1})]);
con = reshape(pspect.IC.sect{3}, [1, numel(pspect.IC.sect{3})]);
inter = reshape(pspect.IC.sect{2}, [1, numel(pspect.IC.sect{2})]);
[pspect.IC.h(1), pspect.IC.p(1), pspect.IC.ksdist(1)] = kstest2(con, pat);
[pspect.IC.h(2), pspect.IC.p(2), pspect.IC.ksdist(2)] = kstest2(pat, inter);
[pspect.IC.h(3), pspect.IC.p(3), pspect.IC.ksdist(3)] = kstest2(con, inter);

% Visualize power spectral distances
edges = 0:5:50;
F.pspect = figure; hold on;
subplot(2,3,1); histogram(pat, edges); title('Patient Spectral Distances'); subplot(2,3,2); histogram(con, edges); title('Control'); subplot(2,3,3); histogram(inter, edges); title('Inter');
subplot(2,3,4); hold on; histogram(pat, edges); histogram(con, edges); title('Grouped Spectral Distances'); legend({'Patient', 'Control'});
subplot(2,3,5); hold on; histogram(pat, edges); histogram(inter, edges); title('Grouped Spectral Distances'); legend({'Patient', 'Inter'});
subplot(2,3,6); hold on; histogram(con, edges); histogram(inter, edges); title('Grouped Spectral Distances'); legend({'Control', 'Inter'});
clear con pat inter edges

% Compute KS distances between control, patient power periodograms
pat = reshape(pgram.IC.sect{1}, [1, numel(pgram.IC.sect{1})]);
con = reshape(pgram.IC.sect{3}, [1, numel(pgram.IC.sect{3})]);
inter = reshape(pgram.IC.sect{2}, [1, numel(pgram.IC.sect{2})]);
[pgram.IC.h(1), pgram.IC.p(1), pgram.IC.ksdist(1)] = kstest2(con, pat);
[pgram.IC.h(2), pgram.IC.p(2), pgram.IC.ksdist(2)] = kstest2(pat, inter);
[pgram.IC.h(3), pgram.IC.p(3), pgram.IC.ksdist(3)] = kstest2(con, inter);

% Visualize periodogram distances
edges = 0:50:1500;
F.pgram = figure; hold on;
subplot(2,3,1); histogram(pat, edges); title('Patient Periodogram Distances'); subplot(2,3,2); histogram(con, edges); title('Control Periodogram Distances'); subplot(2,3,3); histogram(inter, edges); title('Inter Periodogram Distances');
subplot(2,3,4); hold on; histogram(pat, edges); histogram(con, edges); title('Grouped Periodogram Distances'); legend({'Patient', 'Control'});
subplot(2,3,5); hold on; histogram(pat, edges); histogram(inter, edges); title('Grouped Periodogram Distances'); legend({'Patient', 'Inter'});
subplot(2,3,6); hold on; histogram(con, edges); histogram(inter, edges); title('Grouped Periodogram Distances'); legend({'Control', 'Inter'});
clear con pat inter


%% Compute coherence

% Construct temporary index of activities
dummy = nan(N.IC, T.scan, sum(N.subjects));
i = 0;
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		i = i+1;
		dummy(:,:,i) = activities.subj{s,c};
	end
end
clear c s i

% Find all possible pairings
coms = nchoosek(1:sum(N.subjects), 2);

% Test all possible pairwise coherences
coherence = cell(length(coms), 1);
for c = 1:length(coms)
	coherence{c} = mscohere(squeeze(dummy(:,:,coms(c,1)))', squeeze(dummy(:,:,coms(c,2)))')';
end
clear c dummy

% Split coherence matrices into groups
pat = nan(size(coherence{1})); ip = 0;
con = nan(size(coherence{1})); ic = 0;
inter = nan(size(coherence{1})); ii = 0;
for c = 1:length(coms)
	if coms(c,1) <= N.subjects(1) && coms(c,2) <= N.subjects(1)
		ip = ip+1;
		pat(:,:,ip) = coherence{c};
	elseif coms(c,1) > N.subjects(1) && coms(c,2) > N.subjects(1)
		ic = ic+1;
		con(:,:,ic) = coherence{c};
	else
		ii = ii+1;
		inter(:,:,ii) = coherence{c};
	end
end
clear c ip ic ii

% Compile into structure
clear coherence;
coherence.pat = pat;
coherence.con = con;
coherence.inter = inter;
clear pat inter con

% Visualize coherence averages, standard deviations
F.coherence = figure; hold on;
subplot(3,2,1); imagesc(mean(coherence.pat, 3)); colorbar; title('Average Inter-Patient Coherence'); ylabel('IC'); xlabel('Frequency');
subplot(3,2,2); imagesc(var(coherence.pat, [], 3)); colorbar; title('Variance of Inter-Patient Coherence'); ylabel('IC'); xlabel('Frequency');
subplot(3,2,3); imagesc(mean(coherence.con, 3)); colorbar; title('Average Inter-Control Coherence'); ylabel('IC'); xlabel('Frequency');
subplot(3,2,4); imagesc(var(coherence.con, [], 3)); colorbar; title('Variance of Inter-Control Coherence'); ylabel('IC'); xlabel('Frequency');
subplot(3,2,5); imagesc(mean(coherence.inter, 3)); colorbar; title('Average Inter-Group Coherence'); ylabel('IC'); xlabel('Frequency');
subplot(3,2,6); imagesc(var(coherence.inter, [], 3)); colorbar; title('Variance of Inter-Group Coherence'); ylabel('IC'); xlabel('Frequency');


%% Compare IC metrics between conditions, vs. permuted null distribution

% Define test types
ttype = {'kstest2', 'permutation'};

% Test with Kolmogorov-Smirnov, permutation test
for t = 1:numel(ttype)
	disp(['Running ', ttype{t}, ' tests on activations.']);
	
	% Compare activations between conditions
	sig.BOLD(t) = robustTests(cell2mat(BOLD(:,1)'), cell2mat(BOLD(:,2)'), N.ROI, 'p',pval.target, 'testtype',ttype{t});	% Compare ROI time series
	sig.dFC(t) = robustTests(dFC.cond{1}, dFC.cond{2}, size(dFC.concat,1), 'p',pval.target, 'testtype',ttype{t});					% Compare dFC time series
	sig.IC(t) = robustTests(activities.cond{1}, activities.cond{2}, N.IC, 'p',pval.target, 'testtype',ttype{t});		% Compare IC time series
end
clear con pat t

% Subject metastabilities
disp('Running permutation tests on metastability.');
con = metastable.subj{:,'Control'}(isfinite(metastable.subj{:,'Control'}));
pat = metastable.subj{:,'Patient'}(isfinite(metastable.subj{:,'Patient'}));
[sig.metastable(1).h, sig.metastable(1).p, sig.metastable(1).effsize] = kstest2(con, pat);
[sig.metastable(2).p, ~, sig.metastable(2).effsize] = permutationTest(con, pat, 10000, 'sidedness','both');
clear con pat

% Subject entropies
disp('Running permutation tests on entropy.');
con = entro.subj{:,'Control'}(isfinite(entro.subj{:,'Control'}));
pat = entro.subj{:,'Patient'}(isfinite(entro.subj{:,'Patient'}));
[sig.entro(1).h, sig.entro(1).p, sig.entro(1).effsize] = kstest2(con, pat);
[sig.entro(2).p, ~, sig.entro(2).effsize] = permutationTest(con, pat, 10000, 'sidedness','both');
clear con pat



%% Save results

% Save figures
if exist('F', 'var')
	savefig(F, fullfile(path{6}, fileName), 'compact');
	clear F
end

% Save all data
save(fullfile(path{6}, fileName));

