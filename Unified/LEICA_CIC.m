%% LEICA Analysis of Network States
%	This script extracts and transforms neuroimaging data for ICA-based
% comparisons of the network states.  The current script only computes
% phase-based synchrony measures.  It is possible to compute synchrony
% measures based on the z-score and the BOLD signal, but at the moment this
% will only add complexity.
%   We implement three methods for assigning time points to components.  The
% first is ICA, which assigns an activation to each time point.  The second
% is binary k-means, in which each time point is assigned a single
% time point.  The final option is probabilistic k-means in which each
% component has a certain probability of activation for each time point.


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
path{4,1} = fullfile(path{2}, 'OCD');
path{5,1} = fullfile(path{4}, 'Data');
path{6,1} = fullfile(path{4}, 'Results', 'LEICA');
path{7,1} = fullfile(path{1}, 'Project','Atlases','AAL');

% Add relevant paths
fpath{1,1} = fullfile(path{1}, 'MATLAB','spm12');
fpath{2,1} = fullfile(path{1}, 'MATLAB','FastICA');
fpath{3,1} = fullfile(path{1}, 'MATLAB','permutationTest');
fpath{4,1} = fullfile(path{1}, 'MATLAB','BrainNetViewer');
fpath{5,1} = fullfile(path{2}, 'Functions');
fpath{6,1} = fullfile(path{2}, 'LEICA', 'Functions');
for k = 1:numel(fpath)-1
	addpath(fpath{k});
end
addpath(genpath(fpath{numel(fpath)}));
addpath(fullfile(path{1},'MATLAB','spm12'));
clear fpath k


%% Load data

% Load formatted data
load(fullfile(path{5}, 'formattedOCD.mat'));

% Set figure counter
N.fig = 1;


%% Set analysis methods

% Methods key
aType.filter = 'wideband';		% determine which type of filter to use; highpass or bandpass
aType.dist = 'cosine';			% for measuring distance: cosine or exponential
aType.compress = 'eigenvector';	% for compressing matrix: eigenvector, average, or none
aType.segment = 'ICA';			% determine which segmentation to use: ICA, kmeans, or binary (k-means: only assigns states as ON or OFF)

% Set number of neighbors to search for in KNN
co = HShannon_kNN_k_initialization(1);

% Set hypothesis test parameters
pval.target = 0.05;

% Set group on which to define componnets
compGroup = "All";
% importIC = fullfile(path{6}, 'GroupIC/LE_COS_ICA_ControlIC_ALL_wideband_k1_Iteration1');

% Determine how to normalize activity
normal = 'none';    % options: 'none' (no normalization) or 'probability' (normalize time series as probabilities of activation)

% Determine pairwise comparisons to make
comps = "all";
if strcmpi(comps, "all")
    comps = nchoosek(1:N.conditions, 2);    % all comparisons
else
    comps = find(matches(groups, comps))';
end
N.comp = size(comps,1);

% Define test type(s)
ttypes = ["permutation"];


%% Set comparison parameters

% Determine whether to z-scale membership weights
z.thresh = 1:0.2:1.6;
z.scale = false;

% set color index for histograms, cortex network plots
cind.hist = [0 0 0; 0 0 1; 1 0 0; 0 1 0];
cind.node = [1 0 0; 0 0 1];
cind.conn = [1 0 0; 0 0 1];

% Set decompositions, spaces, comparisions
spaces = ["IC"];                        % space in which to compare
dim = ["Subject" "Component"];			% dimensions to compare
pTarget = 0.05;							% target p-value
prange = 0.025:0.025:0.075;				% Set range of p-values to test


%% Define filename based on parameters

% Name compression type
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

% Name distance measurement
switch aType.dist
	case 'cosine'
		fileName = strcat(fileName, '_COS');
	case 'exponential'
		fileName = strcat(fileName, '_EXP');
end

% Name group which defines ICs
switch compGroup
    case {'All', 'ALL', 'all'}
        fileName = fullfile('Common_ICs', strcat(fileName, '_', aType.segment));
    case {'Import', 'IMPORT', 'import'}
        a = strsplit(importIC, '/');
        b = strsplit(path{6}, '/');
        a = a{length(b)+1};
        fileName = fullfile(a, strcat(fileName, '_', aType.segment, '_', compGroup, 'IC'));
        clear a b
    otherwise
        fileName = fullfile(strcat(compGroup,"_ICs"), strcat(fileName, '_', aType.segment));
end

% Name preprocessing pipeline
if exist('ncorr', 'var')
    ncorr = strsplit(ncorr,'+');
    ncorr = ncorr(end);
    fileName = fullfile(ncorr, fileName);
    clear ncorr
end

% Name: pairwise or all-way comparisons?
if N.comp == 1
    fileName = strsplit(fileName, '/');
    fileName = fullfile(fileName{1}, 'Pairwise', strcat(strjoin(fileName(2:end),'/'), '_', groups(comps(1)), 'v', groups(comps(2))));
else
    fileName = strjoin([fileName; "All"], '_');
end

% Set iteration number
fList = dir(fullfile(path{6}, strcat(fileName, '*.mat')));	% Get file list
a = false(numel(fList),1);
for n = 1:numel(fList)
    a(n) = matches('entroCompare.mat', strsplit(fList(n).name, '_'));
end
nIter = numel(fList)-sum(a)+1;

% Set full filename
fileName = strcat(fileName, '_', aType.filter, '_k', num2str(co.mult), '_Iteration', num2str(nIter));
clear fList nIter a


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
if strcmpi(compGroup, "ALL") || strcmpi(compGroup, "IMPORT")
    dFC.concat = cell2mat(dFC.cond);
else
    dFC.concat = cell2mat(dFC.cond(find(matches(groups, compGroup, 'IgnoreCase',true))'));
end
clear t s c afilt bfilt sc90

% Plot BOLD signals
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,2,1); imagesc(cell2mat(BOLD(:,1)')); colorbar; title('Patient BOLD');
subplot(2,2,2); imagesc(cell2mat(BOLD(:,2)')); colorbar; title('Control BOLD');
subplot(2,2,[3 4]); hold on; legend('Patient', 'Control');
histogram(cell2mat(BOLD(:,1)'), 'Normalization','pdf');
histogram(cell2mat(BOLD(:,2)'), 'Normalization','pdf');

% Plot dFC signals
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,2,1); imagesc(dFC.cond{1}); colorbar; title('Patient LEdFC');
subplot(2,2,2); imagesc(dFC.cond{2}); colorbar; title('Control LEdFC');
subplot(2,2,[3 4]); hold on;
histogram(dFC.cond{1}, 'Normalization','pdf');
histogram(dFC.cond{2}, 'Normalization','pdf');
legend('Patient', 'Control');

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

% Compute assembly activity timecourses and memberships
switch compGroup
    case {'import', 'Import', 'IMPORT'}
        disp('Loading ICs from import file');
        e = load(importIC, 'aType', 'N');
        aType.segment = e.aType.segment;
        N.IC = e.N.IC; clear e
        switch aType.segment
            case 'ICA'
                load(importIC, 'memberships', 'W');
                activities.concat = W*dFC.concat;
            case 'binary'
                load(importIC, 'memberships', 'idx');
                activities.concat = zeros(N.IC, length(idx));
                for i = 1:N.IC
                    activities.concat(i, :) = (idx == i)';
                end
                clear i
            case 'kmeans'
                load(importIC, 'memberships', 'D');
                activities.concat = 1./D;
        end
    otherwise
        switch aType.segment
            case 'ICA'
                disp('Processing ICs from dFC data');
                [~, memberships, W] = fastica(dFC.concat, 'numOfIC', N.IC, 'verbose','off');
                dFC.concat = cell2mat(dFC.cond);
                activities.concat = W*dFC.concat;
            case 'binary'
                disp('Processing clusters from dFC data');
                [idx, memberships] = kmeans(dFC.concat', N.IC);
                if ~strcmpi(compGroup,'ALL') && ~strcmpi(compGroup,'IMPORT')
                    dFC.concat = cell2mat(dFC.cond);
                    X = vertcat(memberships, dFC.concat');
                    D = squareform(pdist(X));
                    D = D(1:N.IC, N.IC+(1:size(T.index,2)));
                    [~,idx] = min(D);
                    activities.concat = zeros(N.IC, size(dFC.concat,2));
                    for i = 1:N.IC
                        activities.concat(i, :) = (idx == i);
                    end
                else
                    activities.concat = zeros(N.IC, length(idx));
                    for i = 1:N.IC
                        activities.concat(i, :) = (idx == i)';
                    end
                end
                memberships = memberships';
            case 'kmeans'
                disp('Processing clusters from dFC data');
                [~, memberships, ~, D] = kmeans(dFC.concat', N.IC);
                D = D';
                if ~strcmpi(compGroup,'ALL') && ~strcmpi(compGroup,'IMPORT')
                    dFC.concat = cell2mat(dFC.cond);
                    X = vertcat(memberships, dFC.concat');
                    D = squareform(pdist(X));
                    D = D(1:N.IC, 1+N.IC:end);
                    activities.concat = 1./D;
                else
                    activities.concat = 1./D;
                end
                memberships = memberships';
        end
end
clear importIC D X idx i

% Normalize activity (if desired)
switch normal
    case 'probability'
        activities.concat = activities.concat./sum(activities.concat, 1);   % scale as probabilities
end
meanActivity.concat = mean(activities.concat, 2);

% Sort memberships by activation level
[meanActivity.concat, i] = sort(meanActivity.concat, 1, 'descend');
memberships = memberships(:,i);
clear i

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
subplot(2,3, [3 6]); hold on;
histogram(cell2mat(activities.subj(1:N.subjects(1),1)), 'Normalization','pdf');
histogram(cell2mat(activities.subj(1:N.subjects(2),2)), 'Normalization','pdf');
legend({'Patient', 'Control'});

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

% Run comparison with entropy
[F, h, p, tstat, FDR, Sidak] = compareComponents(path, z, cind, spaces, dim, ttypes, pTarget, prange, origin, cortex, sphereScale, rdux, ROI, entro, N, memberships, I.Properties.VariableNames, comps);



%% Visualize FCD Distributions

% Visualize dFC FCD
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax(2,1) = subplot(2, N.conditions, N.conditions+1:2*N.conditions); hold on;
for c = 1:N.conditions
    histogram(ax(2,1), cell2mat(FCD.dFC.subj(1:N.subjects(c),c)), 'Normalization','pdf');	% FCD histograms
    ax(1,c) = subplot(2, N.conditions, c); colormap jet             % FCD matrix
    imagesc(ax(1,c), FCD.dFC.subj{1,c}); colorbar; title(["dFC FCD of ", groups(c), num2str(1)]);
end
legend(ax(2,1), groups);  % FCD histogram legend
clear ax c

% Visualize IC FCD
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax(2,1) = subplot(2, N.conditions, N.conditions+1:2*N.conditions); hold on;
for c = 1:N.conditions
    histogram(ax(2,1), cell2mat(FCD.IC.subj(1:N.subjects(c),c)), 'Normalization','pdf');	% FCD histograms
    ax(1,c) = subplot(2, N.conditions, c); colormap jet             % FCD matrix
    imagesc(ax(1,c), FCD.IC.subj{1,c}); colorbar; title(["IC FCD of ", groups(c), num2str(1)]);
end
legend(ax(2,1), groups);  % FCD histogram legend
clear ax c


%% Visualize Metastability Distributions

% Get metastability bin sizes
f = figure; hold on;
hg{1} = histogram(metastable.dFC{:, groups(1)}, 'Normalization','pdf');
hg{2} = histogram(metastable.dFC{:, groups(2)}, 'Normalization','pdf');
sz(1) = min(hg{1}.BinWidth, hg{2}.BinWidth);
hg{1} = histogram(metastable.IC{:, groups(1)}, 'Normalization','pdf');
hg{2} = histogram(metastable.IC{:, groups(2)}, 'Normalization','pdf');
sz(2) = min(hg{1}.BinWidth, hg{2}.BinWidth);
close(f); clear hg f

% Visualize metastability
F(N.fig) = figure; hold on; sgtitle('Metastability'); N.fig = N.fig+1;
ax(1) = subplot(1,2,1); hold on; title('dFC Metastability');
ax(2) = subplot(1,2,2); hold on; title('IC Metastability');
for c = 1:N.conditions
    histogram(ax(1), metastable.dFC{:,groups(c)}, 'BinWidth',sz(1), 'Normalization','pdf');
    histogram(ax(2), metastable.IC{:,groups(c)}, 'BinWidth',sz(2), 'Normalization','pdf');
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
		histogram(entro.IC(j,:,1), 'BinWidth',sz, 'Normalization','pdf');
		histogram(entro.IC(j,:,2), 'BinWidth',sz, 'Normalization','pdf');
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
clear mships j hg sz f ind a c


%% Compare component metrics

% Entropy
[F, h, p, tstat, FDR, Sidak] = compareComponents(path,z,cind, spaces,dim,ttypes,pTarget,prange, origin,cortex,sphereScale,rdux, ROI,entro,N,memberships,labels,comps);


%% Save results

% Save figures
if exist('F', 'var')
	savefig(F, fullfile(path{6}, fileName), 'compact');
	clear F
end

% Save all data
save(fullfile(path{6}, fileName));