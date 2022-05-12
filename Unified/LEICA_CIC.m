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
path{4,1} = fullfile(path{2}, 'UCLA');
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
load(fullfile(path{5}, 'formattedUCLA_GMR.mat'));
cortex.file = fullfile(path{7}, cortex.file);

% Enforce: groups must be row string (necessary for boxplot grouping)
if iscolumn(groups)
    groups = groups';
end

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
compGroup = "Control";					% note: letter cases must match save directory
% importIC = fullfile(path{6}, "Common_ICs/LE_COS_ICA_ControlvOCD_wideband_k1_Iteration2");

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


%% Set comparison parameters

% Determine whether to z-scale membership weights
z.thresh = 1:0.2:1.6;
z.scale = false;

% set color index for histograms, cortex network plots
cind.hist = [0 0 0; 0 0 1; 1 0 0; 0 1 0];
cind.node = [1 0 0; 0 0 1];
cind.conn = [1 0 0; 0 0 1];

% Set tests, decompositions, spaces, comparisions
ttype = "permutation";
spaces = "IC";                          % space(s) in which to compare
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
    fileName = fullfile(fileName{1}, strcat(strjoin(fileName(2:end),'/'), '_', groups(comps(1)), 'v', groups(comps(2))));
%     fileName = fullfile(fileName{1}, 'Pairwise', strcat(strjoin(fileName(2:end),'/'), '_', groups(comps(1)), 'v', groups(comps(2))));
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
dFC.subj = cell(max(N.subjects), N.conditions);
T.index = nan(2, sum(N.subjects)*T.scan);
t = zeros(2,1);

% Compute subject-level BOLD phase and dFC
disp('Computing subject-level dFC');
dFC.cond = cell(1, N.conditions);
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		[~, dFC.subj{s,c}] = phasesync(BOLD{s,c}, N.ROI, T.scan, bfilt, afilt, aType);
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

% Compute FCD, power spectra of dFC
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		FCD.dFC.subj{s,c} = computeFCD(dFC.subj{s,c}, 'cosine');
		% pspect.dFC.subj{s,c} = pspectrum(dFC.subj{s,c}', 1/T.TR)';
	end
end


%% Visualize BOLD, dFC signals

% Plot BOLD signals
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax = subplot(N.conditions+1 ,1, N.conditions+1); hold on;
for g = 1:N.conditions
    subplot(N.conditions+1 ,1, g); imagesc(cell2mat(BOLD(:,1)')); colorbar; title(groups{g}); xlabel('Time (s)'); ylabel('ROI');
    histogram(ax, cell2mat(BOLD(:,g)'), 'Normalization','pdf');
end
legend(ax, groups);
sgtitle('BOLD');

% Plot dFC signals
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax = subplot(N.conditions+1 ,1, N.conditions+1); hold on;
for g = 1:N.conditions
    subplot(N.conditions+1 ,1, g); imagesc(dFC.cond{g}); colorbar; title(groups{g}); xlabel('Time (s)'); ylabel('ROI');
    histogram(ax, dFC.cond{g}, 'Normalization','pdf');
end
legend(ax, groups);
sgtitle('LEdFC');


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
activities.concat = activities.concat(i,:);
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


%% Visualize IC activations & relation to FC

% Visualize IC activations
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax = subplot(N.conditions+1 ,1, N.conditions+1); hold on;
for g = 1:N.conditions
    subplot(N.conditions+1 ,1, g); imagesc(cell2mat(activities.subj(1:N.subjects(g),g)')); colorbar; title(groups{g}); xlabel('Time (s)'); ylabel('ROI');
    histogram(ax, cell2mat(activities.subj(1:N.subjects(g),g)), 'Normalization','pdf');
end
legend(ax, groups);
sgtitle('LEICA Activations');
clear g

% Visualize weighted average vs. static FC
F(N.fig) = figure; hold on; N.fig = N.fig+1; colormap jet
subplot(1,2, 2); imagesc(mean(FC, [3 4], 'omitnan')); colorbar; title('Static FC');  xticks([]); yticks([]);
subplot(1,2, 1); imagesc(mean(d,3)); colorbar;
title('Weighted Motif Average'); yticks(1:N.ROI); yticklabels(ROI{:,"Label"}); xticks([]);
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

% Find mean values (stored in table format)
metastable.BOLD = array2table(metastable.BOLD, 'VariableNames', groups);
metastable.dFC = array2table(metastable.dFC, 'VariableNames', groups);
metastable.IC = array2table(metastable.IC, 'VariableNames', groups);
entro.subj = array2table(squeeze(mean(entro.IC, 1, 'omitnan')), 'VariableNames', groups);
entro.mIC = array2table(squeeze(mean(entro.IC, 2, 'omitnan')), 'VariableNames', groups);
fcomp.subj = array2table(squeeze(mean(fcomp.IC, 1, 'omitnan')), 'VariableNames', groups);
fcomp.mIC = array2table(squeeze(mean(fcomp.IC, 2, 'omitnan')), 'VariableNames', groups);

% Find joint values (necessary because sum(NaN) == 0)
entro.joint = squeeze(sum(entro.IC, 1, 'omitnan'));
entro.joint(entro.joint == 0) = NaN;
entro.joint = array2table(entro.joint, 'VariableNames', groups);


%% Visualize FCD Distributions

% Visualize dFC FCD
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax(2,1) = subplot(2, N.conditions, N.conditions+1:2*N.conditions); hold on;
for c = 1:N.conditions
    histogram(ax(2,1), cell2mat(FCD.dFC.subj(1:N.subjects(c),c)), 'Normalization','pdf');	% FCD histograms
    ax(1,c) = subplot(2, N.conditions, c); colormap jet             % FCD matrix
    imagesc(ax(1,c), FCD.dFC.subj{1,c}); colorbar; title(strjoin(["dFC FCD of", groups(c), num2str(1)])); ylabel("Time (s)"); xlabel("Time (s)");  pbaspect([1 1 1]);
end
legend(ax(2,1), groups);  % FCD histogram legend
ylabel(ax(2,1), "Count");
xlabel(ax(2,1), "Correlation");
title(ax(2,1), "FCD Scores");
clear ax c

% Visualize IC FCD
F(N.fig) = figure; hold on; N.fig = N.fig+1;
ax(2,1) = subplot(2, N.conditions, N.conditions+1:2*N.conditions); hold on;
for c = 1:N.conditions
    histogram(ax(2,1), cell2mat(FCD.IC.subj(1:N.subjects(c),c)), 'Normalization','pdf');	% FCD histograms
    ax(1,c) = subplot(2, N.conditions, c); colormap jet             % FCD matrix
    imagesc(ax(1,c), FCD.IC.subj{1,c}); colorbar; title(strjoin(["IC FCD of", groups(c), num2str(1)])); ylabel("Time (s)"); xlabel("Time (s)"); pbaspect([1 1 1]);
end
legend(ax(2,1), groups);  % FCD histogram legend
ylabel(ax(2,1), "Count");
xlabel(ax(2,1), "Correlation");
title(ax(2,1), "FCD Scores");
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

if ~strcmpi(aType.compress, 'none')
    
    % Scale memberships (optional)
    if z.scale == true || sum(z.thresh ~= 0, 'all') > 0
        mships = squeeze(zscore(memberships));
    else
        mships = squeeze(memberships);
    end
    
    % Plot components
    fDims(1,:) = [0 0 1280 1024];
    fDims(2,:) = [0 0 480 1024];
    F(N.fig:N.fig+N.IC+1) = plotComponents(mships, z.thresh(1), ROI, cortex, origin, cind, N, fDims);
    N.fig = N.fig+N.IC+2;
end
clear mships fDims


%% Compare joint metrics

% Generate dummy variables (necessary for compatibility)
ej(1,:,:) = table2array(entro.joint);
nic = N.IC; N.IC = 1;

% Run permutation test on joint entropy
[hjoint, pjoint, ~, FDRjoint, Sidakjoint] = compareComponentMetrics(ttype, pTarget, ej, N, groups, comps);
% [hjoint, pjoint, ~, ~, ~] = robustTests(entro.joint{:,1}', entro.joint{:,2}', [], 'p',pTarget, 'testtype',ttype);
% hjoint = table(hjoint, 'VariableNames', strjoin([groups(comps(1)), "v.", groups(comps(2))], ' '));
% pjoint = table(pjoint, 'VariableNames', strjoin([groups(comps(1)), "v.", groups(comps(2))], ' '));

% Plot results for joint entropy
fDims = [0 1024-256 1280 256];
F(N.fig) = plotComponentMetrics(1, ej, "Joint Entropy", comps, groups, FDRjoint, cind, [1.07 1.1 1.13], fDims, []);

% Reconvert dummy variables
N.IC = nic;
clear ej nic


%% Compare component metrics

% Component Entropies
fDims = [0 1024-256 1280 256];
for s = 1:numel(spaces)
    N.fig = N.fig + 1;
    [h, p, tstat, FDR, Sidak] = compareComponentMetrics(ttype, pTarget, entro.(spaces(s)), N, groups, comps);
    if ~strcmpi(aType.compress, 'none')
        F(N.fig) = plotComponentMetrics(N.IC, entro.(spaces(s)), "Entropy", comps, groups, FDR, cind, [1.03 1.07 1.1], fDims, []);
    end
end


%% Save results

% Save figures
if exist('F', 'var') && strcmpi(aType.compress, 'none')
	savefig(F, fullfile(path{6}, fileName), 'compact');
    saveas(F(8), fullfile(path{6}, strjoin([fileName, "JointEntropy"], '_')), 'jpeg');
	clear F ax
elseif exist('F', 'var')
	savefig(F, fullfile(path{6}, fileName), 'compact');
    saveas(F(8), fullfile(path{6}, strjoin([fileName, "SpatialMaps"], '_')), 'jpeg');
    saveas(F(9), fullfile(path{6}, strjoin([fileName, "ComponentMaps"], '_')), 'jpeg');
    for k = 10:N.IC+9
        saveas(F(k), fullfile(path{6}, strjoin([fileName, strcat("Comp", num2str(k-9))], '_')), 'jpeg');
    end
    saveas(F(N.IC+10), fullfile(path{6}, strjoin([fileName, "JointEntropy"], '_')), 'jpeg');
    saveas(F(N.IC+11), fullfile(path{6}, strjoin([fileName, "ComponentEntropy"], '_')), 'jpeg');
	clear F ax
end

% Save all data
save(fullfile(path{6}, fileName));