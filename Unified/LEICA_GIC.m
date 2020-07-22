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
clear fpath k

% Load structural data
load(fullfile(path{4}, 'sc90.mat'));

% Set methods
distType = 'cosine';	% for measuring distance: cosine or exponential
compressType = 'none';	% for compressing matrix: eigenvector, average, or none

% File to save
switch compressType
	case {'LEICA', 'eigenvector'}
		fileName = 'LEICA90_GIC';
	case 'average'
		fileName = 'MICA90_GIC';
	otherwise
		fileName = 'ICA90_GIC';
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
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,2,1); imagesc(cell2mat(BOLD(:,1)')); colorbar; title('Patient BOLD');
subplot(2,2,2); imagesc(cell2mat(BOLD(:,1)')); colorbar; title('Control BOLD');
subplot(2,2,[3 4]); hold on; histogram(cell2mat(BOLD(:,1)')); histogram(cell2mat(BOLD(:,2)')); legend('Patient', 'Control');

% Compute BOLD phase and z-score
PH = cell(max(N.subjects), N.conditions);
disp('Computing phase of BOLD signal');
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		PH{s,c} = phaseTS(BOLD{s,c}, [], [], bfilt, afilt);
	end
end
clear s c

% Preallocate storage arrays
dFC.cond = cell(N.conditions, 1);
switch compressType
	case {'LEICA', 'eigenvector', 'average'}
		for c = 1:N.conditions
			dFC.cond{c} = zeros(N.ROI, T.scan*N.subjects(c));
		end
	otherwise
		for c = 1:N.conditions
			dFC.cond{c} = zeros(length(Isubdiag), T.scan*N.subjects(c));
		end
end

% Compute instantaneous FC (BOLD Phase Synchrony) and leading eigenvector (V1) for each time point
T.index = zeros(2, T.scan*max(N.subjects), c);	% vector with subject nr and task at each t
for c = 1:N.conditions
	t = 0;			% Index of time (starts at 0, updated until N.subjects*T.scan)
	for s = 1:N.subjects(c)
		
		% Index subject, condition of current dFC sequence
		t = t+1 : t+T.scan;
		T.index(:,t,c) = repmat([s c]', 1,T.scan);
		
		% Extract dFC
		dFC.cond{c}(:,t) = LEdFC(PH{s,c}, 'distType',distType, 'compressType',compressType, 'nROI',N.ROI, 'T',T.scan);
		t = t(end);
	end
end
clear afilt bfilt c s t m n iPH iZ V1 t_all Isubdiag sc90

% Segment dFC
dFC.subj = cell(max(N.subjects), N.conditions);
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		I = (T.index(2,:,c) == c & T.index(1,:,c) == s);
		dFC.subj{s,c} = dFC.cond{c}(:,I);
	end
end
clear I s c

% Plot LEdFC signals
F(N.fig) = figure; hold on; N.fig = N.fig+1;
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
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,2, 1); imagesc(cell2mat(FCD.dFC.subj(1:N.subjects(1),1))'); colorbar; title('Patient dFC FCD');
subplot(2,2, 2); imagesc(cell2mat(FCD.dFC.subj(1:N.subjects(2),2))'); colorbar; title('Control dFC FCD');
subplot(2,2, [3 4]); hold on; histogram(pat); histogram(con); legend({'Patient', 'Control'});
clear con pat



%% Compute ICs from dFC

% Preallocate arrays
activities.cond = cell(1, N.conditions);
activities.subj = cell(max(N.subjects), N.conditions);
memberships = cell(1, N.conditions);
W = cell(1, N.conditions);
if strcmpi(compressType, {'LEICA', 'eigenvector', 'average'})
	ICs = cell(1, N.conditions);
end

% Loop through conditions
explained = nan(N.ROI, N.conditions);
explainedVar = nan(1, N.conditions);
for c = 1:N.conditions
	
	% Extract number of assemblies using Marcenko-Pastur distribution
	N.IC(c) = NumberofIC(dFC.cond{c});
	[~,~,~,~,explained(:,c),~] = pca(dFC.cond{c}');
	explainedVar(c) = sum(explained(N.IC(c)+1:end, c));
	
	% Compute assembly activity timecourses and memberships
	disp('Processing the ICs from BOLD data');
	[activities.cond{c}, memberships{c}, W{c}] = fastica(dFC.cond{c}, 'numOfIC', N.IC(c), 'verbose','off');
	
	% Normalize membership weights
	for k = 1:N.IC(c)
		memberships{c}(:,k) = memberships{c}(:,k)./norm(memberships{c}(:,k));
	end
	
	% Separate assembly activations by subject
	for s = 1:N.subjects(c)
		I = (T.index(2,:,c) == c & T.index(1,:,c) == s);
		activities.subj{s,c} = activities.cond{c}(:,I);
	end
	
	% Compute component matrices
	if strcmpi(compressType, {'LEICA', 'eigenvector', 'average'})
		ICs{c} = nan(N.ROI, N.ROI, size(memberships{c},2));
		for i = 1:size(memberships,2)
			ICs{c}(:,:,i) = memberships(:,i) * memberships(:,i)';
		end
	end
end
clear I s c k i


% Visualize IC activations
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,2, 1); imagesc(cell2mat(activities.subj(1:N.subjects(1),1)')); colorbar; title('Patient LEICA Activations');
subplot(2,2, 2); imagesc(cell2mat(activities.subj(1:N.subjects(2),2)')); colorbar; title('Control LEICA Activations');
subplot(2,2, [3 4]); hold on; histogram(cell2mat(activities.subj(1:N.subjects(1),1))); histogram(cell2mat(activities.subj(1:N.subjects(2),2))); legend({'Patient', 'Control'});



%% Compute IC metrics

% Preallocate storage arrays:
%	Component-wise metastability
%	Component-wise and subject-wise entropies
metastable.BOLD = nan(max(N.subjects), N.conditions);
metastable.PH = nan(max(N.subjects), N.conditions);
metastable.IC = nan(max(N.subjects), N.conditions);
metastable.dFC = nan(max(N.subjects), N.conditions);
entro.all = nan(max(N.IC), max(N.subjects), N.conditions);
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		[~, metastable.BOLD(s,c)] = findStability(BOLD{s,c});
		[~, metastable.PH(s,c)] = findStability(PH{s,c});
		[~, metastable.dFC(s,c)] = findStability(dFC.subj{s,c});
		for g = 1:N.conditions
			[~, metastable.IC(s,c,g)] = findStability(activities.subj{s,c,g});
			for ass = 1:N.IC(c)
				entro.all(ass,s,c,g) = HShannon_kNN_k_estimation(activities.subj{s,c,g}(ass,:), co);
			end
		end
	end
end
clear c s ass g
entro.subj = squeeze(mean(entro.all, 1, 'omitnan'));
entro.IC = squeeze(mean(entro.all, 2, 'omitnan'));

% Convert metrics to table format
metastable.BOLD = array2table(metastable.BOLD, 'VariableNames', labels.Properties.VariableNames);
metastable.PH = array2table(metastable.PH, 'VariableNames', labels.Properties.VariableNames);
metastable.dFC = array2table(metastable.dFC, 'VariableNames', labels.Properties.VariableNames);
metastable.IC = array2table(metastable.IC, 'VariableNames', labels.Properties.VariableNames);
entro.subj = array2table(entro.subj, 'VariableNames', labels.Properties.VariableNames);
entro.IC = array2table(entro.IC, 'VariableNames', labels.Properties.VariableNames);

% Visualize IC metrics
F(N.fig) = figure; hold on; N.fig = N.fig+1; title('Metastability');
subplot(1,2,1); hold on; histogram(metastable.dFC{:,'Patient'}); histogram(metastable.dFC{:,'Control'}); legend('Patients', 'Controls'); title('dFC Metastability');
subplot(1,2,2); hold on; histogram(metastable.IC{:,'Patient'}); histogram(metastable.IC{:,'Control'}); legend('Patients', 'Controls'); title('IC Metastability');
F(N.fig) = figure; hold on; N.fig = N.fig+1; title('Entropy');
subplot(1,2,1); hold on; histogram(entro.subj{:,'Patient'}); histogram(entro.subj{:,'Control'}); legend('Patients', 'Controls'); title('Mean Subject Entropy');
subplot(1,2,2); hold on; histogram(entro.IC{:,'Patient'}); histogram(entro.IC{:,'Control'}); legend('Patients', 'Controls'); title('Mean IC Entropy');


%% Compute FCD, power spectra, and periodograms of ICs

% Compute FCD, power spectra, periodogram of subject-level ICs
for c = 1:N.conditions
	for g = 1:N.conditions
		for s = 1:N.subjects(c)
			FCD.IC.subj{s,c,g} = computeFCD(activities.subj{s,c,g}, 'cosine');
			pspect.IC.subj{s,c,g} = pspectrum(activities.subj{s,c,g}', 1/T.TR)';
			pgram.IC.subj{s,c,g} = periodogram(activities.subj{s,c,g}', [], [], 1/T.TR)';
		end
	end
end
clear c s i g

% Visualize FCD
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,2,1); imagesc(cell2mat(FCD.IC.subj(:,1))'); colorbar; title('LEICA FCD: Patients');
subplot(2,2,2); imagesc(cell2mat(FCD.IC.subj(:,2))'); colorbar; title('LEICA FCD: Controls');
subplot(2,2,[3 4]); hold on; histogram(cell2mat(FCD.IC.subj(:,1))'); histogram(cell2mat(FCD.IC.subj(:,2))'); legend({'Patient', 'Control'});

% Compute KS distances between control, patient FCD
pat = reshape(cell2mat(FCD.IC.subj(1:N.subjects(1),1)), [N.subjects(1)*175^2, 1]);
con = reshape(cell2mat(FCD.IC.subj(1:N.subjects(2),2)), [N.subjects(2)*175^2, 1]);
[FCD.IC.h, FCD.IC.p, FCD.IC.ksdist] = kstest2(con, pat);
clear con pat

% Compute Euclidean distances between subject power spectra, periodograms
dummy.spect = nan(max(N.IC), size(pspect.IC.subj{1,1},2), sum(N.subjects));
dummy.gram = nan(max(N.IC), size(pgram.IC.subj{1,1},2), sum(N.subjects));
pspect.IC.dist = nan(sum(N.subjects), sum(N.subjects), max(N.IC));
pgram.IC.dist = nan(sum(N.subjects), sum(N.subjects), max(N.IC));
ind = 0;
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		ind = ind+1;
		dummy.spect(1:N.IC(c),:,ind) = pspect.IC.subj{s,c};
		dummy.gram(1:N.IC(c),:,ind) = pgram.IC.subj{s,c};
	end
end
for c = 1:N.conditions
	for i = 1:N.IC(c)
		pspect.IC.dist(:,:,i) = squareform(pdist(squeeze(dummy.spect(i,:,:))'));
		pgram.IC.dist(:,:,i) = squareform(pdist(squeeze(dummy.gram(i,:,:))'));
	end
end
clear c s i ind dummy

% Compute mean, standard deviation of inter-subject distance
pspect.IC.ave = mean(pspect.IC.dist, 3, 'omitnan');
pspect.IC.std = std(pspect.IC.dist, [], 3, 'omitnan');
pgram.IC.ave = mean(pgram.IC.dist, 3, 'omitnan');
pgram.IC.std = std(pgram.IC.dist, [], 3, 'omitnan');

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
F(N.fig) = figure; hold on; N.fig = N.fig+1;
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
F(N.fig) = figure; hold on; N.fig = N.fig+1;
subplot(2,3,1); histogram(pat, edges); title('Patient Periodogram Distances'); subplot(2,3,2); histogram(con, edges); title('Control Periodogram Distances'); subplot(2,3,3); histogram(inter, edges); title('Inter Periodogram Distances');
subplot(2,3,4); hold on; histogram(pat, edges); histogram(con, edges); title('Grouped Periodogram Distances'); legend({'Patient', 'Control'});
subplot(2,3,5); hold on; histogram(pat, edges); histogram(inter, edges); title('Grouped Periodogram Distances'); legend({'Patient', 'Inter'});
subplot(2,3,6); hold on; histogram(con, edges); histogram(inter, edges); title('Grouped Periodogram Distances'); legend({'Control', 'Inter'});
clear con pat inter


% %% Compute coherence
% 
% % Construct temporary index of activities
% dummy = cell(1, N.conditions);
% for c = 1:N.conditions
% 	dummy{c} = nan(N.IC(c), T.scan, N.subjects(c));
% 	for s = 1:N.subjects(c)
% 		dummy{c}(:,:,s) = activities.subj{s,c};
% 	end
% end
% clear c s
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

% Define test types
ttype = {'kstest2', 'permutation'};

% Test with Kolmogorov-Smirnov, permutation test
for t = 1:numel(ttype)
	disp(['Running ', ttype{t}, ' tests.']);
	
	% Compare activations between conditions
	% sig.BOLD(t) = robustTests(cell2mat(BOLD(:,1)'), cell2mat(BOLD(:,2)'), N.ROI, 'p',pval.target, 'testtype',ttype{t});				% Compare ROI time series
	% sig.dFC(t) = robustTests(dFC.cond{1}, dFC.cond{2}, size(dFC.cond,1), 'p',pval.target, 'testtype',ttype{t});					% Compare dFC time series
	% sig.IC(t) = robustTests(activities.cond{1}, activities.cond{2}, N.IC, 'p',pval.target, 'testtype',ttype{t});					% Compare IC time series
	% sig.entro.IC(t) = robustTests(squeeze(entro.all(:,:,1)), squeeze(entro.all(:,:,2)), N.IC, 'p',pval.target, 'testtype',ttype{t});	% Compare IC time series
end
clear con pat t

% Subject dFC metastabilities
disp('Running permutation tests on metastability.');
con = metastable.dFC{:,'Control'}(isfinite(metastable.dFC{:,'Control'}));
pat = metastable.dFC{:,'Patient'}(isfinite(metastable.dFC{:,'Patient'}));
[sig.metastable.dFC(1).h, sig.metastable.dFC(1).p, sig.metastable.dFC(1).effsize] = kstest2(con, pat);
[sig.metastable.dFC(2).p, ~, sig.metastable.dFC(2).effsize] = permutationTest(con, pat, 10000, 'sidedness','both');
clear con pat

% Subject IC metastabilities
disp('Running permutation tests on metastability.');
con = metastable.IC{:,'Control'}(isfinite(metastable.IC{:,'Control'}));
pat = metastable.IC{:,'Patient'}(isfinite(metastable.IC{:,'Patient'}));
[sig.metastable.IC(1).h, sig.metastable.IC(1).p, sig.metastable.IC(1).effsize] = kstest2(con, pat);
[sig.metastable.IC(2).p, ~, sig.metastable.IC(2).effsize] = permutationTest(con, pat, 10000, 'sidedness','both');
clear con pat

% Subject entropies
disp('Running permutation tests on entropy.');
con = entro.subj{:,'Control'}(isfinite(entro.subj{:,'Control'}));
pat = entro.subj{:,'Patient'}(isfinite(entro.subj{:,'Patient'}));
[sig.entro.meansubj(1).h, sig.entro.meansubj(1).p, sig.entro(1).meansubj.effsize] = kstest2(con, pat);
[sig.entro.meansubj(2).p, ~, sig.entro.meansubj(2).effsize] = permutationTest(con, pat, 10000, 'sidedness','both');
clear con pat

% IC entropies
disp('Running permutation tests on entropy.');
con = entro.IC{:,'Control'}(isfinite(entro.IC{:,'Control'}));
pat = entro.IC{:,'Patient'}(isfinite(entro.IC{:,'Patient'}));
[sig.entro.meanIC(1).h, sig.entro.meanIC(1).p, sig.entro.meanIC(1).effsize] = kstest2(con, pat);
[sig.entro.meanIC(2).p, ~, sig.entro.meanIC(2).effsize] = permutationTest(con, pat, 10000, 'sidedness','both');
clear con pat



%% Save results

% Save figures
if exist('F', 'var')
	savefig(F, fullfile(path{6}, fileName), 'compact');
	clear F
end

% Save all data
save(fullfile(path{6}, fileName));