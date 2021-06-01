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
path{4,1} = fullfile(path{2}, 'OCD', 'Data');
path{5,1} = fullfile(path{2}, 'OCD', 'Results');
path{6,1} = fullfile(path{2}, 'OCD', 'Results','LEICA');

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
aType.dist = 'cosine';	% for measuring distance: cosine or exponential
aType.compress = 'eigenvector';	% for compressing matrix: eigenvector, average, or none
aType.filter = 'wideband';	% determine which type of filter to use; highpass or bandpass
aType.segment = 'kmeans';		% determine which segmentation to use: ICA, k-means, or binary k-means (only assigns states as ON or OFF)

% Set number of neighbors to search for in KNN
co = HShannon_kNN_k_initialization(1);

% Set highpass or bandpass filter
T.TR = 2.73;				% Repetition Time (seconds)
fnq = 1/(2*T.TR);			% Nyquist frequency
k = 2;						% 2nd order butterworth filter\
switch aType.filter
    case 'bandpass'
        flp = 0.04;     % lowpass frequency of filter (Hz)
        fhi = 0.07;		% highpass
        Wn = [flp/fnq fhi/fnq];	% butterworth bandpass non-dimensional frequency
        fType = 'bandpass';
    case 'wideband'
        flp = 0.01;     % lowpass frequency of filter (Hz)
        fhi = 0.1;		% highpass
        Wn = [flp/fnq fhi/fnq];	% butterworth bandpass non-dimensional frequency
        fType = 'bandpass';
    otherwise
        Wn = flp/fnq;			% butterworth highpass non-dimensional frequency
        fType = 'high';
end
clear flp fhi

% Set filename to save
switch aType.compress
	case {'LE', 'eigenvector'}
		fileName = 'LE';
	case 'average'
		fileName = 'M';
	otherwise
		fileName = 'dFC';
end
switch aType.dist
	case 'cosine'
		fileName = strcat(fileName, '_', aType.segment, '_AAL90_CIC_COS');
	case 'exponential'
		fileName = strcat(fileName, '_', aType.segment, '_AAL90_CIC_EXP');
end
fileName = strcat(fileName, '_', aType.filter, '_k', num2str(co.mult));

fList = dir(fullfile(path{6}, strcat(fileName, '*')));	% Get file list
nIter = numel(fList) + 1; clear fList;						% Find number of previous iterations
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
for n = 1:nFiles
	
	% Rename files
	D = strsplit(fList(n).name, '_');
	fName{n} = D{1};
	
	% Visualize data being extracted
	disp(['Extracting data for subject ', fName{n}, '.']);
	
	% Load data
	subjectBOLD(:,:,n) = load(fullfile(path{4}, loadDirectory, fList(n).name))';
	
	% Convert LR-symmetric to LR-mirrored ordering
	subjectBOLD(:,:,n) = LR_version_symm(subjectBOLD(:,:,n));
	
end
clear D ans n loadDirectory loadFiles nFiles fList

% Remove cerebellar regions
if exist('I', 'var')
	subjectBOLD = subjectBOLD(I,:,:);
	N.ROI = size(subjectBOLD,1);
	clear I
end
clear nROI dROI

% Symmetrize node coordinates, labels
load(fullfile(path{2},'Atlases','AAL','aal_cog.txt'));
load(fullfile(path{2},'Atlases','AAL','AAL_labels.mat'));
coords_AAL90 = LR_version_symm(aal_cog);
label_AAL90 = string(LR_version_symm(label90));
clear aal_cog label90

% Extract indices for BOLD signals of patients, controls
label_groups = ["Controls","OCD"];
I(:,1) = ~ismember(fName, patientData{:, 'Code'});
I(:,2) = ismember(fName, patientData{:, 'Code'});
I = array2table(I, 'VariableNames', label_groups);

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
labels = cell2table(labels, 'VariableNames', label_groups);
clear c fName

% Filter signal
[bfilt,afilt] = butter(k,Wn,fType);
clear fnq flp fhi Wn k fType

% Indices for vectorizing lower triangle
Isubdiag = find(triu(ones(N.ROI), 1));

% Set hypothesis test parameters
pval.target = 0.05;

% Set figure counter
N.fig = 1;



%% Extract dFC

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
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,2,1); imagesc(cell2mat(BOLD(:,1)')); colorbar; title('Patient BOLD');
subplot(2,2,2); imagesc(cell2mat(BOLD(:,2)')); colorbar; title('Control BOLD');
subplot(2,2,[3 4]); hold on; histogram(cell2mat(BOLD(:,1)')); histogram(cell2mat(BOLD(:,2)')); legend('Patient', 'Control');

% Preallocate storage arrays
PH = cell(max(N.subjects), N.conditions);
dFC.subj = cell(max(N.subjects), N.conditions);
T.index = nan(N.conditions, sum(N.subjects)*T.scan);
t = zeros(2,1);

% Compute subject-level BOLD phase and dFC
disp('Computing subject-level dFC');
dFC.cond = cell(1, N.conditions);
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		[PH{s,c}, dFC.subj{s,c}] = phasesync(BOLD{s,c}, N.ROI, T.scan, bfilt, afilt, aType);
		t(1) = t(2) + 1;
		t(2) = t(1)-1 + size(dFC.subj{s,c}, 2);
		T.index(:, t(1):t(2)) = repmat([s,c]', [1, size(dFC.subj{s,c},2)]);
	end
	dFC.cond{c} = cell2mat(dFC.subj(1:N.subjects(c),c)');
end
dFC.concat = cell2mat(dFC.cond);
clear t s c afilt bfilt Isubdiag sc90

% Plot dFC signals
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,2,1); imagesc(dFC.cond{1}); colorbar; title('Patient LEdFC');
subplot(2,2,2); imagesc(dFC.cond{2}); colorbar; title('Control LEdFC');
subplot(2,2,[3 4]); hold on; histogram(dFC.cond{1}); histogram(dFC.cond{2}); legend('Patient', 'Control');

% Compute FCD, power spectra of dFC
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		FCD.dFC.subj{s,c} = computeFCD(dFC.subj{s,c}, 'cosine');
		% pspect.dFC.subj{s,c} = pspectrum(dFC.subj{s,c}', 1/T.TR)';
	end
end


%% Compute ICs from dFC

% Extract number of assemblies using Marcenko-Pastur distribution
N.IC = NumberofIC(dFC.concat);

% Use PCA to estimate amount of variance captured with N.IC
[~,~,~,~,explained,~] = pca(dFC.concat');
explainedVar = sum(explained(N.IC+1:end));

% Check how many eigenvalues needed to explain 95% of variance
% [eVal, ~] = eig(covarianceMatrix);
% find(cumsum(eVal)/sum(eVal) > 0.95);

% Compute assembly activity timecourses and memberships
switch aType.segment
	case 'ICA'
		disp('Processing ICs from dFC data');
		[activities.concat, memberships, W] = fastica(dFC.concat, 'numOfIC', N.IC, 'verbose','off');
	case 'binary'
		disp('Processing clusters from dFC data');
		[idx, memberships] = kmeans(dFC.concat', N.IC);
        memberships = memberships';
		activities.concat = zeros(N.IC, length(idx));
		for i = 1:N.IC
			activities.concat(i, :) = (idx == i)';
		end
		clear i
	case 'kmeans'
		disp('Processing clusters from dFC data');
		[idx, memberships, ~, D] = kmeans(dFC.concat', N.IC);
        memberships = memberships'; D = D';
		activities.concat = 1./D;
		activities.concat = activities.concat./max(activities.concat, [], 'all');
end
meanActivity.concat = mean(activities.concat, 2);

% Separate assembly activations by condition & subject
activities.cond = cell(1, N.conditions);
activities.subj = cell(max(N.subjects), N.conditions);
meanActivity.cond = nan(N.IC, N.conditions);
meanActivity.subj = nan(N.IC, N.conditions, max(N.subjects));
for c = 1:N.conditions
	I = T.index(2,:) == c;
	activities.cond{c} = activities.concat(:,I);
	for s = 1:N.subjects(c)
		I = (T.index(2,:) == c & T.index(1,:) == s);
		activities.subj{s,c} = activities.concat(:,I);
		meanActivity.subj(:,c,s) = mean(activities.subj{s,c}, 2);
	end
	meanActivity.cond(:,c) = mean(activities.cond{c}, 2);
end
clear I s c

% Sort memberships by activation level and normalize weights
[meanActivity.concat, i] = sort(meanActivity.concat, 1, 'descend');
memberships = memberships(:,i);
for k = 1:N.IC
	memberships(:,k) = memberships(:,k)./norm(memberships(:,k));
end
clear k i

% Compute component matrices
if strcmpi(aType.compress, 'LEICA') || strcmpi(aType.compress, 'eigenvector') || strcmpi(aType.compress, 'average')
	ICs = nan(N.ROI, N.ROI, N.IC);
	for i = 1:size(memberships,2)
		ICs(:,:,i) = memberships(:,i) * memberships(:,i)';
	end
end
clear i

% Visualize IC activations
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,3, [1 2]); imagesc(cell2mat(activities.subj(1:N.subjects(1),1)')); colorbar; title('Patient LEICA Activations');
subplot(2,3, [4 5]); imagesc(cell2mat(activities.subj(1:N.subjects(2),2)')); colorbar; title('Control LEICA Activations');
subplot(2,3, [3 6]); hold on; histogram(cell2mat(activities.subj(1:N.subjects(1),1))); histogram(cell2mat(activities.subj(1:N.subjects(2),2))); legend({'Patient', 'Control'});

% Compute FCD of subject-level ICs
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		FCD.IC.subj{s,c} = computeFCD(activities.subj{s,c}, 'cosine');
		% pspect.IC.subj{s,c} = pspectrum(activities.subj{s,c}', 1/T.TR)';
		% pgram.IC.subj{s,c} = periodogram(activities.subj{s,c}', [], [], 1/T.TR)';
	end
end
clear c s i

% Compute correlation between static FC, weighted component average
d = nan(N.ROI, N.ROI, N.IC);
if strcmpi(aType.compress, 'none')
	m = zeros(N.ROI);
	for c = 1:N.IC
		m(find(tril(ones(N.ROI), -1))) = memberships(:,c);
		m = m + m';
		m = m./max(abs(m), [], 'all', 'omitnan') + eye(N.ROI);
		d(:,:,c) = meanActivity.concat(c).*m(:,:);
	end
else
	for c = 1:N.IC
		d(:,:,c) = meanActivity.concat(c).*(memberships(:,c)*memberships(:,c)');
	end
end
d = d./max(abs(d), [], 'all', 'omitnan');
rho = corrcoef(mean(d, 3), mean(FC, [3 4], 'omitnan'));
F(N.fig) = figure; hold on; N.fig = N.fig+1; colormap jet
subplot(1,2, 2); imagesc(mean(FC, [3 4], 'omitnan')); colorbar; title('Static FC');  xticks([]); yticks([]);
subplot(1,2, 1); imagesc(mean(d,3)); colorbar;
title('Weighted Motif Average'); yticks(1:N.ROI); yticklabels(label_AAL90); xticks([]);
clear c d i m



%% Compute global metrics (entropy, metastability)

% Preallocate storage arrays:
%	Activation magnitude means, medians, standard deviations
%	Component-wise Kuramoto order parameter & metastability
%	Subject-wise entropies
%	Subject-wise temporal complexity
metastable.IC = nan(max(N.subjects), N.conditions);
metastable.dFC = nan(max(N.subjects), N.conditions);
metastable.BOLD = nan(max(N.subjects), N.conditions);
entro.BOLD = nan(N.ROI, max(N.subjects), N.conditions);
entro.dFC = nan(N.ROI, max(N.subjects), N.conditions);
entro.IC = nan(N.IC, max(N.subjects), N.conditions);
fcomp.BOLD = nan(N.ROI, max(N.subjects), N.conditions);
fcomp.dFC = nan(N.ROI, max(N.subjects), N.conditions);
fcomp.IC = nan(N.IC, max(N.subjects), N.conditions);
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		[~, metastable.BOLD(s,c)] = findStability(BOLD{s,c});
		[~, metastable.dFC(s,c)] = findStability(dFC.subj{s,c});
		[~, metastable.IC(s,c)] = findStability(activities.subj{s,c});
		for ass = 1:N.IC
			entro.IC(ass, s, c) = HShannon_kNN_k_estimation(activities.subj{s,c}(ass,:), co);
			fcomp.IC(ass, s, c) = funccomp(activities.subj{s,c}(ass,:), []);
		end
		for roi = 1:N.ROI
			entro.BOLD(roi, s, c) = HShannon_kNN_k_estimation(BOLD{s,c}(roi,:), co);
			entro.dFC(roi, s, c) = HShannon_kNN_k_estimation(dFC.subj{s,c}(roi,:), co);
			fcomp.BOLD(roi, s, c) = funccomp(BOLD{s,c}(roi,:), []);
			fcomp.dFC(roi, s, c) = funccomp(dFC.subj{s,c}(roi,:), []);
		end
	end
end
clear c s ass roi
entro.subj = squeeze(mean(entro.IC, 1, 'omitnan'));
entro.mIC = squeeze(mean(entro.IC, 2, 'omitnan'));
fcomp.subj = squeeze(mean(fcomp.IC, 1, 'omitnan'));
fcomp.mIC = squeeze(mean(fcomp.IC, 2, 'omitnan'));

% Convert metrics to table format
metastable.BOLD = array2table(metastable.BOLD, 'VariableNames', label_groups);
metastable.dFC = array2table(metastable.dFC, 'VariableNames', label_groups);
metastable.IC = array2table(metastable.IC, 'VariableNames', label_groups);
entro.subj = array2table(entro.subj, 'VariableNames', label_groups);
entro.mIC = array2table(entro.mIC, 'VariableNames', label_groups);
fcomp.subj = array2table(fcomp.subj, 'VariableNames', label_groups);
fcomp.mIC = array2table(fcomp.mIC, 'VariableNames', label_groups);



%% Power spectral and periodogram analysis of ICs
% 
% % Compute FCD of subject-level ICs
% for c = 1:N.conditions
% 	for s = 1:N.subjects(c)
% 		pspect.IC.subj{s,c} = pspectrum(activities.subj{s,c}', 1/T.TR)';
% 		pgram.IC.subj{s,c} = periodogram(activities.subj{s,c}', [], [], 1/T.TR)';
% 	end
% end
% clear c s i
% 
% % Compute Euclidean distances between subject power spectra, periodograms
% dummy.spect = nan(N.IC, size(pspect.IC.subj{1,1},2), sum(N.subjects));
% dummy.gram = nan(N.IC, size(pgram.IC.subj{1,1},2), sum(N.subjects));
% pspect.IC.dist = nan(sum(N.subjects), sum(N.subjects), N.IC);
% pgram.IC.dist = nan(sum(N.subjects), sum(N.subjects), N.IC);
% ind = 0;
% for c = 1:N.conditions
% 	for s = 1:N.subjects(c)
% 		ind = ind+1;
% 		dummy.spect(:,:,ind) = pspect.IC.subj{s,c};
% 		dummy.gram(:,:,ind) = pgram.IC.subj{s,c};
% 	end
% end
% for i = 1:N.IC
% 	pspect.IC.dist(:,:,i) = squareform(pdist(squeeze(dummy.spect(i,:,:))'));
% 	pgram.IC.dist(:,:,i) = squareform(pdist(squeeze(dummy.gram(i,:,:))'));
% end
% clear c s i ind dummy
% 
% % Compute mean, standard deviation of inter-subject distance
% pspect.IC.ave = mean(pspect.IC.dist, 3);
% pspect.IC.std = std(pspect.IC.dist, [], 3);
% pgram.IC.ave = mean(pgram.IC.dist, 3);
% pgram.IC.std = std(pgram.IC.dist, [], 3);
% 
% % Visualize mean, standard deviations of spectral and periodogram distances
% figure; hold on;
% subplot(1,2,1); imagesc(pspect.IC.ave); colorbar; title('Average Spectral Distance between Subjects')
% subplot(1,2,2); imagesc(pspect.IC.std); colorbar; title('Standard Deviation in Spectral Distance between Subjects')
% figure; hold on;
% subplot(1,2,1); imagesc(pgram.IC.ave); colorbar; title('Average Periodogram Distance between Subjects')
% subplot(1,2,2); imagesc(pgram.IC.std); colorbar; title('Standard Deviation in Periodogram Distance between Subjects')
% 
% % Separate into classes (control, patient, inter)
% pspect.IC.sect{1} = pspect.IC.dist(1:N.subjects(1), 1:N.subjects(1));
% pspect.IC.sect{2} = pspect.IC.dist(1:N.subjects(1), N.subjects(1)+1:sum(N.subjects));
% pspect.IC.sect{3} = pspect.IC.dist(N.subjects(1)+1:sum(N.subjects), N.subjects(1)+1:sum(N.subjects));
% pgram.IC.sect{1} = pgram.IC.dist(1:N.subjects(1), 1:N.subjects(1));
% pgram.IC.sect{2} = pgram.IC.dist(1:N.subjects(1), N.subjects(1)+1:sum(N.subjects));
% pgram.IC.sect{3} = pgram.IC.dist(N.subjects(1)+1:sum(N.subjects), N.subjects(1)+1:sum(N.subjects));
% 
% % Compute KS distances between control, patient power spectra
% pat = reshape(pspect.IC.sect{1}, [1, numel(pspect.IC.sect{1})]);
% con = reshape(pspect.IC.sect{3}, [1, numel(pspect.IC.sect{3})]);
% inter = reshape(pspect.IC.sect{2}, [1, numel(pspect.IC.sect{2})]);
% [pspect.IC.h(1), pspect.IC.p(1), pspect.IC.ksdist(1)] = kstest2(con, pat);
% [pspect.IC.h(2), pspect.IC.p(2), pspect.IC.ksdist(2)] = kstest2(pat, inter);
% [pspect.IC.h(3), pspect.IC.p(3), pspect.IC.ksdist(3)] = kstest2(con, inter);
% 
% % Visualize power spectral distances
% edges = 0:5:50;
% F(N.fig) = figure; hold on; N.fig = N.fig+1;
% subplot(2,3,1); histogram(pat, edges); title('Patient Spectral Distances'); subplot(2,3,2); histogram(con, edges); title('Control'); subplot(2,3,3); histogram(inter, edges); title('Inter');
% subplot(2,3,4); hold on; histogram(pat, edges); histogram(con, edges); title('Grouped Spectral Distances'); legend({'Patient', 'Control'});
% subplot(2,3,5); hold on; histogram(pat, edges); histogram(inter, edges); title('Grouped Spectral Distances'); legend({'Patient', 'Inter'});
% subplot(2,3,6); hold on; histogram(con, edges); histogram(inter, edges); title('Grouped Spectral Distances'); legend({'Control', 'Inter'});
% clear con pat inter edges
% 
% % Compute KS distances between control, patient power periodograms
% pat = reshape(pgram.IC.sect{1}, [1, numel(pgram.IC.sect{1})]);
% con = reshape(pgram.IC.sect{3}, [1, numel(pgram.IC.sect{3})]);
% inter = reshape(pgram.IC.sect{2}, [1, numel(pgram.IC.sect{2})]);
% [pgram.IC.h(1), pgram.IC.p(1), pgram.IC.ksdist(1)] = kstest2(con, pat);
% [pgram.IC.h(2), pgram.IC.p(2), pgram.IC.ksdist(2)] = kstest2(pat, inter);
% [pgram.IC.h(3), pgram.IC.p(3), pgram.IC.ksdist(3)] = kstest2(con, inter);
% 
% % Visualize periodogram distances
% edges = 0:50:1500;
% F(N.fig) = figure; hold on; N.fig = N.fig+1;
% subplot(2,3,1); histogram(pat, edges); title('Patient Periodogram Distances'); subplot(2,3,2); histogram(con, edges); title('Control Periodogram Distances'); subplot(2,3,3); histogram(inter, edges); title('Inter Periodogram Distances');
% subplot(2,3,4); hold on; histogram(pat, edges); histogram(con, edges); title('Grouped Periodogram Distances'); legend({'Patient', 'Control'});
% subplot(2,3,5); hold on; histogram(pat, edges); histogram(inter, edges); title('Grouped Periodogram Distances'); legend({'Patient', 'Inter'});
% subplot(2,3,6); hold on; histogram(con, edges); histogram(inter, edges); title('Grouped Periodogram Distances'); legend({'Control', 'Inter'});
% clear con pat inter


%% Compute coherence
% 
% % Construct temporary index of activities
% dummy = nan(N.IC, T.scan, sum(N.subjects));
% i = 0;
% for c = 1:N.conditions
% 	for s = 1:N.subjects(c)
% 		i = i+1;
% 		dummy(:,:,i) = activities.subj{s,c};
% 	end
% end
% clear c s i
% 
% % Find all possible pairings
% coms = nchoosek(1:sum(N.subjects), 2);
% 
% % Test all possible pairwise coherences
% coherence = cell(length(coms), 1);
% for c = 1:length(coms)
% 	coherence{c} = mscohere(squeeze(dummy(:,:,coms(c,1)))', squeeze(dummy(:,:,coms(c,2)))')';
% end
% clear c dummy
% 
% % Split coherence matrices into groups
% pat = nan(size(coherence{1})); ip = 0;
% con = nan(size(coherence{1})); ic = 0;
% inter = nan(size(coherence{1})); ii = 0;
% for c = 1:length(coms)
% 	if coms(c,1) <= N.subjects(1) && coms(c,2) <= N.subjects(1)
% 		ip = ip+1;
% 		pat(:,:,ip) = coherence{c};
% 	elseif coms(c,1) > N.subjects(1) && coms(c,2) > N.subjects(1)
% 		ic = ic+1;
% 		con(:,:,ic) = coherence{c};
% 	else
% 		ii = ii+1;
% 		inter(:,:,ii) = coherence{c};
% 	end
% end
% clear c ip ic ii
% 
% % Compile into structure
% clear coherence;
% coherence.pat = pat;
% coherence.con = con;
% coherence.inter = inter;
% clear pat inter con
% 
% % Visualize coherence averages, standard deviations
% F(N.fig) = figure; hold on; N.fig = N.fig+1;
% subplot(3,2,1); imagesc(mean(coherence.pat, 3)); colorbar; title('Average Inter-Patient Coherence'); ylabel('IC'); xlabel('Frequency');
% subplot(3,2,2); imagesc(var(coherence.pat, [], 3)); colorbar; title('Variance of Inter-Patient Coherence'); ylabel('IC'); xlabel('Frequency');
% subplot(3,2,3); imagesc(mean(coherence.con, 3)); colorbar; title('Average Inter-Control Coherence'); ylabel('IC'); xlabel('Frequency');
% subplot(3,2,4); imagesc(var(coherence.con, [], 3)); colorbar; title('Variance of Inter-Control Coherence'); ylabel('IC'); xlabel('Frequency');
% subplot(3,2,5); imagesc(mean(coherence.inter, 3)); colorbar; title('Average Inter-Group Coherence'); ylabel('IC'); xlabel('Frequency');
% subplot(3,2,6); imagesc(var(coherence.inter, [], 3)); colorbar; title('Variance of Inter-Group Coherence'); ylabel('IC'); xlabel('Frequency');


%% Compare IC metrics between conditions, vs. permuted null distribution

% Determine all possible pairwise comparisons
C = nchoosek(1:N.conditions, 2);
N.comp = size(C,1);

% Define test types
ttype = {'kstest2', 'permutation'};

% Loop through all pairwise comparisons
for c = 1:N.comp
    
    % Test with Kolmogorov-Smirnov, permutation test
    for t = 1:numel(ttype)
        disp(['Running ', ttype{t}, ' tests.']);

        % Compare activations
        sig.BOLD(c,t) = robustTests(cell2mat(BOLD(:,C(c,1))'), cell2mat(BOLD(:,C(c,2))'), N.ROI, 'p',pval.target, 'testtype',ttype{t});				% Compare ROI time series
        sig.dFC(c,t) = robustTests(dFC.cond{C(c,1)}, dFC.cond{C(c,2)}, size(dFC.concat,1), 'p',pval.target, 'testtype',ttype{t});					% Compare dFC time series
        sig.IC(c,t) = robustTests(activities.cond{C(c,1)}, activities.cond{C(c,2)}, N.IC, 'p',pval.target, 'testtype',ttype{t});					% Compare IC time series
        
        % Compare entropies
        sig.entro.IC(c,t) = robustTests(squeeze(entro.IC(:,:,C(c,1))), squeeze(entro.IC(:,:,C(c,2))), N.IC, 'p',pval.target, 'testtype',ttype{t});          % Compare IC entropies
        sig.entro.BOLD(c,t) = robustTests(squeeze(entro.BOLD(:,:,C(c,1))), squeeze(entro.BOLD(:,:,C(c,2))), N.ROI, 'p',pval.target, 'testtype',ttype{t});	% Compare BOLD entropies
        
        % Compare complexities
        sig.fcomp.IC(c,t) = robustTests(squeeze(fcomp.IC(:,:,C(c,1))), squeeze(fcomp.IC(:,:,C(c,2))), N.IC, 'p',pval.target, 'testtype',ttype{t});          % Compare IC functional complexities
        sig.fcomp.BOLD(c,t) = robustTests(squeeze(fcomp.BOLD(:,:,C(c,1))), squeeze(fcomp.BOLD(:,:,C(c,2))), N.ROI, 'p',pval.target, 'testtype',ttype{t});	% Compare BOLD functional complexities
    end
    
    % Compare dFC FCD distributions
    disp('Comparing dFC FCD.');
    pat = reshape(cell2mat(FCD.dFC.subj(1:N.subjects(C(c,1)),C(c,1))), [N.subjects(C(c,1))*175^2, 1]);
    con = reshape(cell2mat(FCD.dFC.subj(1:N.subjects(C(c,2)),C(c,2))), [N.subjects(C(c,2))*175^2, 1]);
    [FCD.dFC.h(c,1), FCD.dFC.p(c,1), FCD.dFC.effsize(c,1)] = kstest2(con, pat);
    [FCD.dFC.p(c,2), ~, FCD.dFC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if FCD.dFC.p(c,2) < pval.target
        FCD.dFC.h(c,2) = 1;
    else
        FCD.dFC.h(c,2) = 0;
    end
    
    % Compare IC FCD distributions
    disp('Comparing IC FCD.');
    con = reshape(cell2mat(FCD.IC.subj(1:N.subjects(C(c,1)),C(c,1))), [N.subjects(C(c,1))*175^2, 1]);
    pat = reshape(cell2mat(FCD.IC.subj(1:N.subjects(C(c,2)),C(c,2))), [N.subjects(C(c,2))*175^2, 1]);
    [FCD.IC.h(c,1), FCD.IC.p(c,1), FCD.IC.effsize(c,1)] = kstest2(con, pat);
    [FCD.IC.p(c,2), ~, FCD.ID.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if FCD.IC.p(c,2) < pval.target
        FCD.IC.h(c,2) = 1;
    else
        FCD.IC.h(c,2) = 0;
    end
    
    % Compare subject BOLD metastabilities
    disp('Comparing BOLD metastability.');
    con = metastable.BOLD{:,label_groups(C(c,1))}(isfinite(metastable.BOLD{:,label_groups(C(c,1))}));
    pat = metastable.BOLD{:,label_groups(C(c,2))}(isfinite(metastable.BOLD{:,label_groups(C(c,2))}));
    [sig.metastable.BOLD.h(c,1), sig.metastable.BOLD.p(c,1), sig.metastable.BOLD.effsize(c,1)] = kstest2(con, pat);
    [sig.metastable.BOLD.p(c,2), ~, sig.metastable.BOLD.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.metastable.BOLD.p(c,2) < pval.target
        sig.metastable.BOLD.h(c,2) = 1;
    else
        sig.metastable.BOLD.h(c,2) = 0;
    end

    % Compare subject dFC metastabilities
    disp('Comparing dFC metastability.');
    con = metastable.dFC{:,label_groups(C(c,1))}(isfinite(metastable.dFC{:,label_groups(C(c,1))}));
    pat = metastable.dFC{:,label_groups(C(c,2))}(isfinite(metastable.dFC{:,label_groups(C(c,2))}));
    [sig.metastable.dFC.h(c,1), sig.metastable.dFC.p(c,1), sig.metastable.dFC.effsize(c,1)] = kstest2(con, pat);
    [sig.metastable.dFC.p(c,2), ~, sig.metastable.dFC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.metastable.dFC.p(c,2) < pval.target
        sig.metastable.dFC.h(c,2) = 1;
    else
        sig.metastable.dFC.h(c,2) = 0;
    end

    % Subject IC metastabilities
    disp('Comparing component metastability.');
    con = metastable.IC{:,label_groups(C(c,1))}(isfinite(metastable.IC{:,label_groups(C(c,1))}));
    pat = metastable.IC{:,label_groups(C(c,2))}(isfinite(metastable.IC{:,label_groups(C(c,2))}));
    [sig.metastable.IC.h(c,1), sig.metastable.IC.p(c,1), sig.metastable.IC.effsize(c,1)] = kstest2(con, pat);
    [sig.metastable.IC.p(c,2), ~, sig.metastable.IC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.metastable.IC.p(c,2) < pval.target
        sig.metastable.IC.h(c,2) = 1;
    else
        sig.metastable.IC.h(c,2) = 0;
    end

    % Subject entropies
    disp('Comparing subject entropies.');
    con = entro.subj{:,label_groups(C(c,1))}(isfinite(entro.subj{:,label_groups(C(c,1))}));
    pat = entro.subj{:,label_groups(C(c,2))}(isfinite(entro.subj{:,label_groups(C(c,2))}));
    [sig.entro.meansubj.h(c,1), sig.entro.meansubj.p(c,1), sig.entro.meansubj.effsize(c,1)] = kstest2(con, pat);
    [sig.entro.meansubj.p(c,2), ~, sig.entro.meansubj.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.entro.meansubj.p(c,2) < pval.target
        sig.entro.meansubj.h(c,2) = 1;
    else
        sig.entro.meansubj.h(c,2) = 0;
    end

    % IC entropies
    disp('Comparing component entropies.');
    con = entro.mIC{:,label_groups(C(c,1))}(isfinite(entro.mIC{:,label_groups(C(c,1))}));
    pat = entro.mIC{:,label_groups(C(c,2))}(isfinite(entro.mIC{:,label_groups(C(c,2))}));
    [sig.entro.meanIC.h(c,1), sig.entro.meanIC.p(c,1), sig.entro.meanIC.effsize(c,1)] = kstest2(con, pat);
    [sig.entro.meanIC.p(c,2), ~, sig.entro.meanIC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.entro.meanIC.p(c,2) < pval.target
        sig.entro.meanIC.h(c,2) = 1;
    else
        sig.entro.meanIC.h(c,2) = 0;
    end

    % Subject functional complexities
    disp('Comparing subject functional complexity.');
    con = fcomp.subj{:,label_groups(C(c,1))}(isfinite(fcomp.subj{:,label_groups(C(c,1))}));
    pat = fcomp.subj{:,label_groups(C(c,2))}(isfinite(fcomp.subj{:,label_groups(C(c,2))}));
    [sig.fcomp.meansubj.h(c,1), sig.fcomp.meansubj.p(c,1), sig.fcomp.meansubj.effsize(c,1)] = kstest2(con, pat);
    [sig.fcomp.meansubj.p(c,2), ~, sig.fcomp.meansubj.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.fcomp.meansubj.p(c,2) < pval.target
        sig.fcomp.meansubj.h(c,2) = 1;
    else
        sig.fcomp.meansubj.h(c,2) = 0;
    end

    % IC functional complexities
    disp('Comparing component functional complexity.');
    con = fcomp.mIC{:,label_groups(C(c,1))}(isfinite(fcomp.mIC{:,label_groups(C(c,1))}));
    pat = fcomp.mIC{:,label_groups(C(c,2))}(isfinite(fcomp.mIC{:,label_groups(C(c,2))}));
    [sig.fcomp.meanIC.h(c,1), sig.fcomp.meanIC.p(c,1), sig.fcomp.meanIC.effsize(c,1)] = kstest2(con, pat);
    [sig.fcomp.meanIC.p(c,2), ~, sig.fcomp.meanIC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.fcomp.meanIC.p(c,2) < pval.target
        sig.fcomp.meanIC.h(c,2) = 1;
    else
        sig.fcomp.meanIC.h(c,2) = 0;
    end
end
clear con pat t c


%% Visualize FCD Distributions

% Visualize dFC FCD
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax(2,1) = subplot(2, N.conditions, 1:N.conditions); hold on;
for c = 1:N.conditions
    histogram(ax(2,1), cell2mat(FCD.dFC.subj(1:N.subjects(c),c)), 'Normalization','probability');	% FCD histograms
    ax(1,c) = subplot(2, N.conditions, c); colormap jet             % FCD matrix
    imagesc(ax(1,c), FCD.dFC.subj{1,c}); colorbar; title(["dFC FCD of ", label_groups(c), num2str(1)]);
end
legend(ax(2,1), label_groups);  % FCD histogram legend
clear ax c

% Visualize IC dFC
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax(2,1) = subplot(2, N.conditions, 1:N.conditions); hold on;
for c = 1:N.conditions
    histogram(ax(2,1), cell2mat(FCD.IC.subj(1:N.subjects(c),c)), 'Normalization','probability');	% FCD histograms
    ax(1,c) = subplot(2, N.conditions, c); colormap jet             % FCD matrix
    imagesc(ax(1,c), FCD.IC.subj{1,c}); colorbar; title(["IC FCD of ", label_groups(c), num2str(1)]);
end
legend(ax(2,1), label_groups);  % FCD histogram legend
clear ax c


%% Visualize Metastability Distributions

% Get metastability bin sizes
f = figure; hold on;
hg{1} = histogram(metastable.dFC{:, label_groups(1)}, 'Normalization','probability');
hg{2} = histogram(metastable.dFC{:, label_groups(2)}, 'Normalization','probability');
sz(1) = min(hg{1}.BinWidth, hg{2}.BinWidth);
hg{1} = histogram(metastable.IC{:, label_groups(1)}, 'Normalization','probability');
hg{2} = histogram(metastable.IC{:, label_groups(2)}, 'Normalization','probability');
sz(2) = min(hg{1}.BinWidth, hg{2}.BinWidth);
close(f); clear hg f

% Visualize metastability
F(N.fig) = figure; hold on; sgtitle('Metastability'); N.fig = N.fig+1;
ax(1) = subplot(1,2,1); hold on; title('dFC Metastability');
ax(2) = subplot(1,2,2); hold on; title('IC Metastability');
for c = 1:N.conditions
    histogram(ax(1), metastable.dFC{:,label_groups(c)}, 'BinWidth',sz(1), 'Normalization','probability');
    histogram(ax(2), metastable.IC{:,label_groups(c)}, 'BinWidth',sz(2), 'Normalization','probability');
end
legend(ax(1), label_groups);
legend(ax(2), label_groups);


%% Visualize Entropy Distributions

% Get entropy bin sizes
f = figure; hold on;
hg{1} = histogram(entro.subj{:, label_groups(1)}, 'Normalization','probability');
hg{2} = histogram(entro.subj{:, label_groups(2)}, 'Normalization','probability');
sz(1) = min(hg{1}.BinWidth, hg{2}.BinWidth);
hg{1} = histogram(entro.mIC{:, label_groups(1)}, 'Normalization','probability');
hg{2} = histogram(entro.mIC{:, label_groups(2)}, 'Normalization','probability');
sz(2) = min(hg{1}.BinWidth, hg{2}.BinWidth);
close(f); clear hg f

% Visualize entropy
F(N.fig) = figure; hold on; sgtitle('Entropy'); N.fig = N.fig+1;
ax(1) = subplot(1,2,1); hold on; title('Mean per Subject');
ax(2) = subplot(1,2,2); hold on; title('Mean Per IC');
for c = 1:N.conditions
    histogram(ax(1), entro.subj{:, label_groups(c)}, 'BinWidth',sz(1), 'Normalization','probability');
	histogram(ax(2), entro.mIC{:, label_groups(c)}, 'BinWidth',sz(2), 'Normalization','probability');
end
legend(ax(1), label_groups);
legend(ax(2), label_groups);
clear sz


%% Visualize IC memberships
if strcmpi(aType.compress, 'none')
	for j = 1:N.IC

		% Get bin sizes
		f = figure; hold on;
		hg{1} = histogram(entro.IC(j,:,1));
		hg{2} = histogram(entro.IC(j,:,2));
		sz = min(hg{1}.BinWidth, hg{2}.BinWidth);
		close(f);

		% Scale memberships (optional)
		mships = zscore(memberships(:,j));

		% Open figure
		F(nFig) = figure; nFig = nFig + 1; hold on;

		% Connectivity
		kax = subplot(2, 4, [3 4]); hold on;
		sgtitle(['Component ', num2str(j)]);
		a = squeeze(memberships(:,j))*squeeze(memberships(:,j))';
		imagesc(a); colorbar; hold on;
		xlim([1 size(a,2)]); ylim([1 size(a,1)]);
		title('Connectivity');
		yticks(1:N.ROI); yticklabels(label_AAL90); xticks([]);
		pbaspect([1 1 1]);

		% Histogram of component entropies
		kax = subplot(2, 4, [1 2]); hold on;
		histogram(entro.IC(j,:,1), 'BinWidth',sz, 'Normalization','Probability');
		histogram(entro.IC(j,:,2), 'BinWidth',sz, 'Normalization','Probability');
		legend(labels.Properties.VariableNames);
		title('Entropy');
		ylabel('Counts'); xlabel('Mean Entropy');

		% Bar Plots
		kax = subplot(2, 4, [5 8]); hold on;
		if sum(sign(squeeze(memberships(:,j)))) >= 0
			ind(:,1) = mships < -zthresh;	% select node weights which surpass threshold
			ind(:,2) = mships > zthresh;	% select node weights which surpass threshold
			ind(:,3) = mships > -zthresh & mships < 0;
			ind(:,4) = mships < zthresh & mships > 0;
			a = mships; a(~ind(:,1)) = 0; bar(1:N.ROI, a, 'b');
			a = mships; a(~ind(:,2)) = 0; bar(1:N.ROI, a, 'r');
			a = mships; a(~ind(:,3)) = 0; bar(1:N.ROI, a, 'b', 'FaceAlpha',0.3);
			a = mships; a(~ind(:,4)) = 0; bar(1:N.ROI, a, 'r', 'FaceAlpha',0.3);
		elseif sum(sign(squeeze(memberships(:,j)))) < 0
			ind(:,1) = mships > zthresh;	% only plot node weights which surpass threshold
			ind(:,2) = mships < -zthresh;	% only plot node weights which surpass threshold
			ind(:,3) = mships < zthresh & mships > 0;
			ind(:,4) = mships > -zthresh & mships < 0;
			a = mships; a(~ind(:,1)) = 0; bar(1:N.ROI, a, 'b');
			a = mships; a(~ind(:,2)) = 0; bar(1:N.ROI, a, 'r');
			a = mships; a(~ind(:,3)) = 0; bar(1:N.ROI, a, 'b', 'FaceAlpha',0.3);
			a = mships; a(~ind(:,4)) = 0; bar(1:N.ROI, a, 'r', 'FaceAlpha',0.3);
		end
		title(['z-score threshold: ' num2str(zthresh)]);
		xticks(1:N.ROI); xticklabels(label_AAL90); xtickangle(-90);
		xlabel('z-score');
	end
end
clear mships j hg sz f ind a


%% Save results

% Save figures
if exist('F', 'var')
	savefig(F, fullfile(path{6}, fileName), 'compact');
	clear F
end

% Save all data
save(fullfile(path{6}, fileName));