%% LEICA Analysis of Network States
%	This script extracts and transforms neuroimaging data for ICA-based
% comparisons of the network states.  The current script only computes
% phase-based synchrony measures.  It is possible to compute synchrony
% measures based on the z-score and the BOLD signal, but at the moment this
% will only add complexity.


%%	SET UP FILE SYSTEM

% Clear workspace
clear; close all; clc

% Shuffle random seed.  Necessary to avoid repeating random seeds across parallel computations.
rng('shuffle');

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set data-specific subdirectories
path{4,1} = fullfile(path{2}, 'UCLA');
path{5,1} = fullfile(path{4}, 'Data');
path{6,1} = fullfile(path{4}, 'Results', 'LEICA');

% Add relevant paths
fpath{1,1} = fullfile(path{1}, 'MATLAB','spm12');
fpath{2,1} = fullfile(path{1}, 'MATLAB','FastICA');
fpath{3,1} = fullfile(path{1}, 'MATLAB','permutationTest');
fpath{4,1} = fullfile(path{2}, 'Functions');
fpath{5,1} = fullfile(path{2}, 'LEICA', 'Functions');
for k = 1:numel(fpath)-1
	addpath(fpath{k});
end
addpath(genpath(fpath{numel(fpath)}));
clear fpath k


%% Load data

% Load formatted data
load(fullfile(path{5}, 'formattedUCLA.mat'));

% Set figure counter
N.fig = 1;


%% Set analysis methods

% Methods key
aType.filter = 'wideband';  % determine which type of filter to use; highpass or bandpass
aType.dist = 'cosine';      % for measuring distance: cosine or exponential
aType.compress = 'LE';      % for compressing matrix: eigenvector, average, or none
aType.segment = 'ICA';	% determine which segmentation to use: ICA, kmeans, or binary (k-means: only assigns states as ON or OFF)

% Set number of neighbors to search for in KNN
co = HShannon_kNN_k_initialization(1);

% Set hypothesis test parameters
pval.target = 0.05;

% Determine pairwise comparisons to make
C = ["CONTROL", "BIPOLAR"];
if strcmpi(C, "all")
    C = nchoosek(1:N.conditions, 2);    % all comparisons
else
    C = find(matches(groups, C))';
end
N.comp = size(C,1);

% Define test types
ttype = {'kstest2', 'permutation'};


%% Set filename to save

% Set base filename
switch aType.compress
	case {'LE', 'LEICA', 'eigenvector'}
		fileName = 'LE';
	case {'average' , 'mean'}
		fileName = 'M';
    case 'none'
		fileName = 'dFC';
    otherwise
        error("Please select one of the supported compression methods");
end
switch aType.dist
	case 'cosine'
		fileName = strcat(fileName, '_', aType.segment, '_AAL90_CIC_COS');
	case 'exponential'
		fileName = strcat(fileName, '_', aType.segment, '_AAL90_CIC_EXP');
end
if size(C,1) == 1
    fileName = strcat(fileName, '_', groups(C(1)), 'v', groups(C(2)));
end

% Set iteration number
fList = dir(fullfile(path{6}, strcat(fileName, '*')));	% Get file list
if numel(fList) == 0
    nIter = 1;
else
    nIter = numel(fList)+1;
end
clear fList					% Find number of previous iterations

% Set full filename
fileName = strcat(fileName, '_', aType.filter, '_k', num2str(co.mult), '_Iteration', num2str(nIter));
clear nIter


%% Set filter

% Set highpass or bandpass filter
fnq = 1/(2*T.TR);			% Nyquist frequency
k = 2;						% 2nd order butterworth filter
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

% Define signal filter
[bfilt,afilt] = butter(k,Wn,fType);
clear fnq flp fhi Wn k fType


%% Compute dFC, FCD

% Preallocate storage arrays
PH = cell(max(N.subjects), N.conditions);
dFC.subj = cell(max(N.subjects), N.conditions);
T.index = nan(2, sum(N.subjects)*T.scan);
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
if size(C,1) == 1
    dFC.concat = dFC.cond{C(1)};
else
    dFC.concat = cell2mat(dFC.cond);
end
clear t s c afilt bfilt sc90

% Plot BOLD signals
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,2,1); imagesc(cell2mat(BOLD(:,1)')); colorbar; title('Patient BOLD');
subplot(2,2,2); imagesc(cell2mat(BOLD(:,2)')); colorbar; title('Control BOLD');
subplot(2,2,[3 4]); hold on; histogram(cell2mat(BOLD(:,1)')); histogram(cell2mat(BOLD(:,2)')); legend('Patient', 'Control');

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
		[~, memberships, W] = fastica(dFC.concat, 'numOfIC', N.IC, 'verbose','off');
        activities.concat = W*dFC.concat;
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
dFC.concat = cell2mat(dFC.cond);
meanActivity.concat = mean(activities.concat, 2);

% Separate assembly activations by condition & subject
activities.cond = cell(1, N.conditions);
activities.subj = cell(max(N.subjects), N.conditions);
meanActivity.cond = nan(N.IC, N.conditions);
meanActivity.subj = nan(N.IC, N.conditions, max(N.subjects));
for c = 1:N.conditions
	i = T.index(2,:) == c;
	activities.cond{c} = activities.concat(:,i);
	for s = 1:N.subjects(c)
		i = (T.index(2,:) == c & T.index(1,:) == s);
		activities.subj{s,c} = activities.concat(:,i);
		meanActivity.subj(:,c,s) = mean(activities.subj{s,c}, 2);
	end
	meanActivity.cond(:,c) = mean(activities.cond{c}, 2);
end
clear i s c

% Sort memberships by activation level and normalize weights
[meanActivity.concat, i] = sort(meanActivity.concat, 1, 'descend');
memberships = memberships(:,i);
% for k = 1:N.IC
% 	memberships(:,k) = memberships(:,k)./norm(memberships(:,k));
% end
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
title('Weighted Motif Average'); yticks(1:N.ROI); yticklabels(labels_ROI); xticks([]);
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
metastable.BOLD = array2table(metastable.BOLD, 'VariableNames', groups);
metastable.dFC = array2table(metastable.dFC, 'VariableNames', groups);
metastable.IC = array2table(metastable.IC, 'VariableNames', groups);
entro.subj = array2table(entro.subj, 'VariableNames', groups);
entro.mIC = array2table(entro.mIC, 'VariableNames', groups);
fcomp.subj = array2table(fcomp.subj, 'VariableNames', groups);
fcomp.mIC = array2table(fcomp.mIC, 'VariableNames', groups);



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

% Loop through all pairwise comparisons
for c = 1:N.comp
    
    % Test with Kolmogorov-Smirnov, permutation test
    for t = 1:numel(ttype)
        disp(['Running ', ttype{t}, ' tests.']);

        % Compare activations
        [sig.BOLD.h(:,c,t), sig.BOLD.p(:,c,t), sig.BOLD.tstat(:,c,t)] = robustTests(cell2mat(BOLD(:,C(c,1))'), cell2mat(BOLD(:,C(c,2))'), N.ROI, 'p',pval.target, 'testtype',ttype{t});				% Compare ROI time series
        [sig.dFC.h(:,c,t), sig.dFC.p(:,c,t), sig.dFC.tstat(:,c,t)] = robustTests(dFC.cond{C(c,1)}, dFC.cond{C(c,2)}, size(dFC.concat,1), 'p',pval.target, 'testtype',ttype{t});					% Compare dFC time series
        [sig.IC.h(:,c,t), sig.IC.p(:,c,t), sig.IC.tstat(:,c,t)] = robustTests(activities.cond{C(c,1)}, activities.cond{C(c,2)}, N.IC, 'p',pval.target, 'testtype',ttype{t});					% Compare IC time series
        
        % Compare complexities
        [sig.fcomp.BOLD.h(:,c,t), sig.fcomp.BOLD.p(:,c,t), sig.fcomp.BOLD.tstat(:,c,t)] = robustTests(squeeze(fcomp.BOLD(:,:,C(c,1))), squeeze(fcomp.BOLD(:,:,C(c,2))), N.ROI, 'p',pval.target, 'testtype',ttype{t});	% Compare BOLD functional complexities
        [sig.fcomp.IC.h(:,c,t), sig.fcomp.IC.p(:,c,t), sig.fcomp.IC.tstat(:,c,t)] = robustTests(squeeze(fcomp.IC(:,:,C(c,1))), squeeze(fcomp.IC(:,:,C(c,2))), N.IC, 'p',pval.target, 'testtype',ttype{t});          % Compare IC functional complexities
    end
    
    % Compare dFC FCD distributions
    disp('Comparing dFC FCD.');
    pat = reshape(cell2mat(FCD.dFC.subj(1:N.subjects(C(c,1)),C(c,1))), [N.subjects(C(c,1))*T.scan^2, 1]);
    con = reshape(cell2mat(FCD.dFC.subj(1:N.subjects(C(c,2)),C(c,2))), [N.subjects(C(c,2))*T.scan^2, 1]);
    [FCD.dFC.h(c,1), FCD.dFC.p(c,1), FCD.dFC.effsize(c,1)] = kstest2(con, pat);
    [FCD.dFC.p(c,2), ~, FCD.dFC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if FCD.dFC.p(c,2) < pval.target
        FCD.dFC.h(c,2) = 1;
    else
        FCD.dFC.h(c,2) = 0;
    end
    
    % Compare IC FCD distributions
    disp('Comparing IC FCD.');
    con = reshape(cell2mat(FCD.IC.subj(1:N.subjects(C(c,1)),C(c,1))), [N.subjects(C(c,1))*T.scan^2, 1]);
    pat = reshape(cell2mat(FCD.IC.subj(1:N.subjects(C(c,2)),C(c,2))), [N.subjects(C(c,2))*T.scan^2, 1]);
    [FCD.IC.h(c,1), FCD.IC.p(c,1), FCD.IC.effsize(c,1)] = kstest2(con, pat);
    [FCD.IC.p(c,2), ~, FCD.IC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if FCD.IC.p(c,2) < pval.target
        FCD.IC.h(c,2) = 1;
    else
        FCD.IC.h(c,2) = 0;
    end
    
    % Compare subject BOLD metastabilities
    disp('Comparing BOLD metastability.');
    con = metastable.BOLD{:,groups(C(c,1))}(isfinite(metastable.BOLD{:,groups(C(c,1))}));
    pat = metastable.BOLD{:,groups(C(c,2))}(isfinite(metastable.BOLD{:,groups(C(c,2))}));
    [sig.metastable.BOLD.h(c,1), sig.metastable.BOLD.p(c,1), sig.metastable.BOLD.effsize(c,1)] = kstest2(con, pat);
    [sig.metastable.BOLD.p(c,2), ~, sig.metastable.BOLD.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.metastable.BOLD.p(c,2) < pval.target
        sig.metastable.BOLD.h(c,2) = 1;
    else
        sig.metastable.BOLD.h(c,2) = 0;
    end

    % Compare subject dFC metastabilities
    disp('Comparing dFC metastability.');
    con = metastable.dFC{:,groups(C(c,1))}(isfinite(metastable.dFC{:,groups(C(c,1))}));
    pat = metastable.dFC{:,groups(C(c,2))}(isfinite(metastable.dFC{:,groups(C(c,2))}));
    [sig.metastable.dFC.h(c,1), sig.metastable.dFC.p(c,1), sig.metastable.dFC.effsize(c,1)] = kstest2(con, pat);
    [sig.metastable.dFC.p(c,2), ~, sig.metastable.dFC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.metastable.dFC.p(c,2) < pval.target
        sig.metastable.dFC.h(c,2) = 1;
    else
        sig.metastable.dFC.h(c,2) = 0;
    end

    % Subject IC metastabilities
    disp('Comparing component metastability.');
    con = metastable.IC{:,groups(C(c,1))}(isfinite(metastable.IC{:,groups(C(c,1))}));
    pat = metastable.IC{:,groups(C(c,2))}(isfinite(metastable.IC{:,groups(C(c,2))}));
    [sig.metastable.IC.h(c,1), sig.metastable.IC.p(c,1), sig.metastable.IC.effsize(c,1)] = kstest2(con, pat);
    [sig.metastable.IC.p(c,2), ~, sig.metastable.IC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.metastable.IC.p(c,2) < pval.target
        sig.metastable.IC.h(c,2) = 1;
    else
        sig.metastable.IC.h(c,2) = 0;
    end

    % Subject functional complexities
    disp('Comparing subject functional complexity.');
    con = fcomp.subj{:,groups(C(c,1))}(isfinite(fcomp.subj{:,groups(C(c,1))}));
    pat = fcomp.subj{:,groups(C(c,2))}(isfinite(fcomp.subj{:,groups(C(c,2))}));
    [sig.fcomp.meansubj.h(c,1), sig.fcomp.meansubj.p(c,1), sig.fcomp.meansubj.effsize(c,1)] = kstest2(con, pat);
    [sig.fcomp.meansubj.p(c,2), ~, sig.fcomp.meansubj.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.fcomp.meansubj.p(c,2) < pval.target
        sig.fcomp.meansubj.h(c,2) = 1;
    else
        sig.fcomp.meansubj.h(c,2) = 0;
    end

    % IC functional complexities
    disp('Comparing component functional complexity.');
    con = fcomp.mIC{:,groups(C(c,1))}(isfinite(fcomp.mIC{:,groups(C(c,1))}));
    pat = fcomp.mIC{:,groups(C(c,2))}(isfinite(fcomp.mIC{:,groups(C(c,2))}));
    [sig.fcomp.meanIC.h(c,1), sig.fcomp.meanIC.p(c,1), sig.fcomp.meanIC.effsize(c,1)] = kstest2(con, pat);
    [sig.fcomp.meanIC.p(c,2), ~, sig.fcomp.meanIC.effsize(c,2)] = permutationTest(con, pat, 10000, 'sidedness','both');
    if sig.fcomp.meanIC.p(c,2) < pval.target
        sig.fcomp.meanIC.h(c,2) = 1;
    else
        sig.fcomp.meanIC.h(c,2) = 0;
    end
end

% Run mutiple-comparison correction
if c > 1
    for t = 1:numel(ttype)
        
        % BOLD activations
        p_dum = reshape(squeeze(sig.BOLD.p(:,:,t)), [N.comp*N.ROI, 1]);     % reshape p-value array)
        [f, ~, S] = mCompCorr([], p_dum, pval.target);                      % run mutliple comparison tests
        sig.BOLD.FDR(:,:,t) = reshape(f, [N.ROI, N.comp]);
        sig.BOLD.Sidak(:,:,t) = reshape(S, [N.ROI, N.comp]);
        
        % dFC activations
        p_dum = reshape(squeeze(sig.dFC.p(:,:,t)), [N.comp*N.ROI, 1]);      % reshape p-value array)
        [f, ~, S] = mCompCorr([], p_dum, pval.target);                      % run mutliple comparison tests
        sig.dFC.FDR(:,:,t) = reshape(f, [N.ROI, N.comp]);
        sig.dFC.Sidak(:,:,t) = reshape(S, [N.ROI, N.comp]);
        
        % dFC activations
        p_dum = reshape(squeeze(sig.IC.p(:,:,t)), [N.comp*N.IC, 1]);	% reshape p-value array)
        [f, ~, S] = mCompCorr([], p_dum, pval.target);                  % run mutliple comparison tests
        sig.IC.FDR(:,:,t) = reshape(f, [N.IC, N.comp]);
        sig.IC.Sidak(:,:,t) = reshape(S, [N.IC, N.comp]);
        
        % BOLD functional complexity
        p_dum = reshape(squeeze(sig.fcomp.BOLD.p(:,:,t)), [N.comp*N.ROI, 1]);	% reshape p-value array)
        [f, ~, S] = mCompCorr([], p_dum, pval.target);                      % run mutliple comparison tests
        sig.fcomp.BOLD.FDR(:,:,t) = reshape(f, [N.ROI, N.comp]);
        sig.fcomp.BOLD.Sidak(:,:,t) = reshape(S, [N.ROI, N.comp]);
        
        % Component functional complexity
        p_dum = reshape(squeeze(sig.fcomp.IC.p(:,:,t)), [N.comp*N.IC, 1]);	% reshape p-value array)
        [f, ~, S] = mCompCorr([], p_dum, pval.target);                      % run mutliple comparison tests
        sig.fcomp.IC.FDR(:,:,t) = reshape(f, [N.IC, N.comp]);
        sig.fcomp.IC.Sidak(:,:,t) = reshape(S, [N.IC, N.comp]);
        
        % FCD Activations
        [FCD.IC.FDR, FCD.IC.Bonferroni, FCD.IC.Sidak] = mCompCorr([], FCD.IC.p(:,t), pval.target);
        [FCD.dFC.FDR, FCD.dFC.Bonferroni, FCD.dFC.Sidak] = mCompCorr([], FCD.dFC.p(:,t), pval.target);

        % Metastability
        [sig.metastable.BOLD.FDR, sig.metastable.BOLD.Bonferroni, sig.metastable.BOLD.Sidak] = mCompCorr([], sig.metastable.BOLD.p(:,t), pval.target);
        [sig.metastable.dFC.FDR, sig.metastable.dFC.Bonferroni, sig.metastable.dFC.Sidak] = mCompCorr([], sig.metastable.dFC.p(:,t), pval.target);
        [sig.metastable.IC.FDR, sig.metastable.IC.Bonferroni, sig.metastable.IC.Sidak] = mCompCorr([], sig.metastable.IC.p(:,t), pval.target);

        % Mean Complexities
        [sig.fcomp.meansubj.FDR, sig.fcomp.meansubj.Bonferroni, sig.fcomp.meansubj.Sidak] = mCompCorr([], sig.fcomp.meansubj.p(:,t), pval.target);
        [sig.fcomp.meanIC.FDR, sig.fcomp.meanIC.Bonferroni, sig.fcomp.meanIC.Sidak] = mCompCorr([], sig.fcomp.meanIC.p(:,t), pval.target);
    end
end
clear con pat t c pdum f S


%% Visualize FCD Distributions

% Visualize dFC FCD
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax(2,1) = subplot(2, N.conditions, N.conditions+1:2*N.conditions); hold on;
for c = 1:N.conditions
    histogram(ax(2,1), cell2mat(FCD.dFC.subj(1:N.subjects(c),c)), 'Normalization','probability');	% FCD histograms
    ax(1,c) = subplot(2, N.conditions, c); colormap jet             % FCD matrix
    imagesc(ax(1,c), FCD.dFC.subj{1,c}); colorbar; title(["dFC FCD of ", groups(c), num2str(1)]);
end
legend(ax(2,1), groups);  % FCD histogram legend
clear ax c

% Visualize IC FCD
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax(2,1) = subplot(2, N.conditions, N.conditions+1:2*N.conditions); hold on;
for c = 1:N.conditions
    histogram(ax(2,1), cell2mat(FCD.IC.subj(1:N.subjects(c),c)), 'Normalization','probability');	% FCD histograms
    ax(1,c) = subplot(2, N.conditions, c); colormap jet             % FCD matrix
    imagesc(ax(1,c), FCD.IC.subj{1,c}); colorbar; title(["IC FCD of ", groups(c), num2str(1)]);
end
legend(ax(2,1), groups);  % FCD histogram legend
clear ax c


%% Visualize Metastability Distributions

% Get metastability bin sizes
f = figure; hold on;
hg{1} = histogram(metastable.dFC{:, groups(1)}, 'Normalization','probability');
hg{2} = histogram(metastable.dFC{:, groups(2)}, 'Normalization','probability');
sz(1) = min(hg{1}.BinWidth, hg{2}.BinWidth);
hg{1} = histogram(metastable.IC{:, groups(1)}, 'Normalization','probability');
hg{2} = histogram(metastable.IC{:, groups(2)}, 'Normalization','probability');
sz(2) = min(hg{1}.BinWidth, hg{2}.BinWidth);
close(f); clear hg f

% Visualize metastability
F(N.fig) = figure; hold on; sgtitle('Metastability'); N.fig = N.fig+1;
ax(1) = subplot(1,2,1); hold on; title('dFC Metastability');
ax(2) = subplot(1,2,2); hold on; title('IC Metastability');
for c = 1:N.conditions
    histogram(ax(1), metastable.dFC{:,groups(c)}, 'BinWidth',sz(1), 'Normalization','probability');
    histogram(ax(2), metastable.IC{:,groups(c)}, 'BinWidth',sz(2), 'Normalization','probability');
end
legend(ax(1), groups);
legend(ax(2), groups);


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
		yticks(1:N.ROI); yticklabels(labels_ROI); xticks([]);
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
		xticks(1:N.ROI); xticklabels(labels_ROI); xtickangle(-90);
		xlabel('z-score');
	end
end
clear mships j hg sz f ind a


%% Save results

% % Save figures
% if exist('F', 'var')
% 	savefig(F, fullfile(path{6}, fileName), 'compact');
% 	clear F
% end
% 
% % Save all data
% save(fullfile(path{6}, fileName));