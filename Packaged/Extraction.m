%% ICA Extraction of Network States
%	This script extracts and transforms neuroimaging data for ICA-based
% comparisons of the network states.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DATA PIPELINE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%	SETUP

%% Set paths & filenames

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
path{6,1} = fullfile(path{1},'MATLAB','spm12');
path{7,1} = fullfile(path{2},'Functions');
path{8,1} = fullfile(path{2},'Results','LEICA');

% Add relevant paths
addpath(path{6});
addpath(genpath(path{7}));


%% Set file names, analysis type

% Define list of patients
patientList = 'ClinicalData_OCD.xlsx';

% Define files containing BOLD data
loadDirectory = 'Subjects';
loadFiles = 'P*';

% Set methods
phaseType = 'cosine';	% for measuring phase: cosine or exponential
compressType = 'none';	% for compressing matrix: eigenvector, average, or none

% Files to save
if strcmpi(compressType, 'eigenvector')
	fileName = 'LEICA90_Data';				% all results
	phaseFileName = 'LEICA90_PhaseData';	% phase only
elseif strcmpi(compressType, 'average')
	fileName = 'MICA90_Data';				% all results
	phaseFileName = 'MICA90_PhaseData';		% phase only
elseif strcmpi(compressType, 'none')
	fileName = 'ICA90_Data';				% all results
	phaseFileName = 'ICA90_PhaseData';		% phase only
end


%% Get list of files

% Get file list
fList = dir(fullfile(path{4}, loadDirectory, loadFiles));

% Get number of files of interest
nFiles = numel(fList);


%% Load data

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


%% Set indices

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
labels = cell2table(labels, 'VariableNames',{'Controls','Patients'});
clear c fName

% Save BOLD data, indices
save(fullfile(path{8},fileName));


%% Set parameters

% Temporal parameters
T.TR = 2;			% Repetition Time (seconds)
T.condition(1) = T.scan*N.subjects(1);
T.condition(2) = T.scan*N.subjects(2);

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


%% Set pairwise comparisions (unnecessary for OCD because only OCD vs. controls)



%% SORTING

%% Extract BOLD time series data to cell array (TS.BOLD)

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

% Visualize interim results
load(fullfile(path{4}, 'sc90.mat'));
figure; imagesc(sc90); title('Structural');
figure; imagesc(squeeze(mean(FC(:,:,:,1),3,'omitnan'))); title('Patients');
figure; imagesc(squeeze(mean(FC(:,:,:,2),3,'omitnan'))); title('Controls');

% Save interim results
save(fullfile(path{8}, fileName));



%% CONVERSION

%% 7) Convert data to z-score, phase data

% Preallocate storage arrays
TS.PH = cell(max(N.subjects), N.conditions);
TS.Z = cell(max(N.subjects), N.conditions);

% Compute BOLD phase and z-score
disp('Computing z-score and phase of BOLD signal');
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		
		% Extract filtered time series, Hilbert phase from BOLD data
		[TS.PH{s,c}, TS.BOLD{s,c}] = regionPhase(TS.BOLD{s,c}, bfilt, afilt);
		
		% Extract metastability, entropy from BOLD data
		[kuramoto.AAL.subj(s,c,:), metastable.AAL.subj(s,c)] = findStability(TS.PH{s,c});
		entro.AAL.subj(s,c) = HShannon_kNN_k_estimation(TS.PH{s,c}, co);
		
		% Convert filtered BOLD data to z-score
		TS.Z{s,c} = zscore(TS.BOLD{s,c}')';
	end
end
clear s c


%% Extract metastability, synchrony, entropy from phase time series

% Preallocate arrays
kuramoto.AAL.concat = nan(1, sum(T.condition));
kuramoto.AAL.cond = nan(N.conditions, max(T.condition));
kuramoto.AAL.subj = nan(max(N.subjects), N.conditions, T.scan);
metastable.AAL.cond = nan(1, N.conditions);
metastable.AAL.subj = nan(max(N.subjects), N.conditions);
entro.AAL.cond = nan(1, N.conditions);
entro.AAL.subj = nan(max(N.subjects), N.conditions);

% Extract metastability, synchrony, entropy
z = TS.PH';
t = cell(2,1);

% Extract metastability, entropy from BOLD data
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		[kuramoto.AAL.subj(s,c,:), metastable.AAL.subj(s,c)] = findStability(TS.PH{s,c});
		entro.AAL.subj(s,c) = HShannon_kNN_k_estimation(TS.PH{s,c}, co);
	end
	t{c} = cell2mat(z(c,:));
	[kuramoto.AAL.cond(c, 1:T.condition(c)), metastable.AAL.cond(c)] = findStability(t{c});
	entro.AAL.cond(c) = HShannon_kNN_k_estimation(t{c}, co);
end
t = horzcat(t{1}, t{2});
[kuramoto.AAL.concat, metastable.AAL.concat] = findStability(t);
entro.AAL.concat = HShannon_kNN_k_estimation(t, co);
clear s c z t

% Convert metastability, entropy data to tables
metastable.AAL.subj = array2table(metastable.AAL.subj, 'VariableNames', {'Controls','OCD'});
metastable.AAL.cond = array2table(metastable.AAL.cond, 'VariableNames', {'Controls','OCD'});
entro.AAL.subj = array2table(entro.AAL.subj, 'VariableNames', {'Controls','OCD'});
entro.AAL.cond = array2table(entro.AAL.cond, 'VariableNames', {'Controls','OCD'});


%% 8) Convert z-score to point process

% % Preallocate arrays
% TS.PP = cell(max(N.subjects), N.conditions);
% TS.A = cell(max(N.subjects), N.conditions);
% 
% % Compute events, activations
% for c = 1:N.conditions
% 	for s = 1:N.subjects(c)
% 		TS.PP{s,c} = eventMat(TS.Z{s,c});
% 		TS.A{s,c} = activeMat(TS.Z{s,c});
% 	end
% end


%% 9) Compute dFC tensors

% Preallocate storage arrays
if strcmpi(compressType, 'LEICA') ||  strcmpi(compressType, 'Eigenvector')	% leading eigenvector
	dFC.PH.concat = zeros(sum(T.condition), N.ROI);
	dFC.Z.concat = zeros(sum(T.condition), N.ROI);
elseif strcmpi(compressType, 'average')	% average activity
	dFC.PH.concat = zeros(sum(T.condition), N.ROI);
	dFC.Z.concat = zeros(sum(T.condition), N.ROI);
else	% vectorized synchrony data
	dFC.PH.concat = zeros(sum(T.condition), length(Isubdiag));
	dFC.Z.concat = zeros(sum(T.condition), length(Isubdiag));
end
dFC.BOLD.concat = [];

% Preallocate variables to save FC patterns and associated information
T.index = zeros(2, sum(T.condition));	% vector with subject nr and task at each t
t_all = 0;								% Index of time (starts at 0, updated until N.subjects*T.scan)

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
						iPH(n,m)=exp(-3*adif(TS.PH{s,c}(n,t), TS.PH{s,c}(m,t))); 
					end
				end
			end
			
			% Compute z-score matrix
			iZ = corr(TS.Z{s,c}(:,t)')';
			
			% Extract leading eigenvector
			if strcmpi(compressType, 'eigenvector')
				[dFC.PH.concat(t_all,:), ~] = eigs(iPH,1);
				[dFC.Z.concat(t_all,:), ~] = eigs(iZ,1);
			elseif strcmpi(compressType, 'average')
				dFC.PH.concat(t_all,:) = mean(iPH);
				dFC.PH.concat(t_all,:) = mean(iZ);
			else
				dFC.PH.concat(t_all,:) = iPH(Isubdiag);
				dFC.Z.concat(t_all,:) = iZ(Isubdiag);
			end
			T.index(:,t_all) = [s c]';	% Information that at t_all, V1 corresponds to subject s in a given task
		end
		
		dFC.BOLD.concat = horzcat(dFC.BOLD.concat, TS.BOLD{s,c});
		% dFC.Z.concat(N.conditions(1)*(c-1)+(1+T.scan*(s-1):T.scan*s), :) = TS.Z{s,c}';
	end
end
dFC.PH.concat = dFC.PH.concat';
dFC.Z.concat = dFC.Z.concat';

% Clear dud variables
clear afilt bfilt c s t m n iPH iZ V1 t_all Isubdiag sc90


%% Split dFC into condition, subject signals

% Preallocate storage arrays
dFC.PH.cond = cell(1, N.conditions);
dFC.Z.cond = cell(1, N.conditions);
dFC.BOLD.cond = cell(1, N.conditions);
dFC.PH.subj = cell(max(N.subjects), N.conditions);
dFC.Z.subj = cell(max(N.subjects), N.conditions);
dFC.BOLD.subj = cell(max(N.subjects), N.conditions);

% Segment dFC
for c = 1:N.conditions
	I = T.index(2,:) == c;
	dFC.PH.cond{c} = dFC.PH.concat(:,I);
	dFC.Z.cond{c} = dFC.Z.concat(:,I);
	dFC.BOLD.cond{c} = dFC.BOLD.concat(:,I);
	
	for s = 1:N.subjects(c)
		I = (T.index(2,:) == c & T.index(1,:) == s);
		dFC.PH.subj{s,c} = dFC.PH.concat(:,I);
		dFC.Z.subj{s,c} = dFC.Z.concat(:,I);
		dFC.BOLD.subj{s,c} = dFC.BOLD.concat(:,I);
	end
end
clear I s c


%% Save extracted, packaged data

% Save all data
save(fullfile(path{8},fileName));

% Save phase-only data
dFC = dFC.PH;
save(fullfile(path{8},phaseFileName));

