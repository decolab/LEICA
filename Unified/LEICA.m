%% LEICA Extraction of Network States
%	This script extracts and transforms neuroimaging data for ICA-based
% comparisons of the network states.  The current script only computes
% phase-based synchrony measures.  It is possible to compute synchrony
% measures based on the z-score and the BOLD signal, but at the moment this
% will only add complexity.


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
path{7,1} = fullfile(path{1},'MATLAB','FastICA');
path{8,1} = fullfile(path{2},'Functions');
path{9,1} = fullfile(path{2},'Results','LEICA');

% Add relevant paths
addpath(path{6});
addpath(path{7});
addpath(genpath(path{8}));


%% Set file names, analysis type

% Define list of patients
patientList = 'ClinicalData_OCD.xlsx';

% Define files containing BOLD data
loadDirectory = 'Subjects';
loadFiles = 'P*';

% Set methods
phaseType = 'cosine';	% for measuring phase: cosine or exponential
compressType = 'eigenvector';	% for compressing matrix: eigenvector, average, or none

% Files to save
if strcmpi(compressType, 'eigenvector')
	fileName = 'LEICA90_Data';
elseif strcmpi(compressType, 'average')
	fileName = 'MICA90_Data';
elseif strcmpi(compressType, 'none')
	fileName = 'ICA90_Data';
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
figure; imagesc(squeeze(mean(FC(:,:,:,1),3,'omitnan'))); title('Mean FC: Patients');
figure; imagesc(squeeze(mean(FC(:,:,:,2),3,'omitnan'))); title('Mean FC: Controls');

% Save interim results
save(fullfile(path{8}, fileName));



%% CONVERSION

%% 7) Demean BOLD signal and convert data to phase angle time series

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


%% 8) Compute dFC tensors

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


%% Split dFC into condition, subject signals

% Preallocate storage arrays
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


%% Save extracted data

% Get file list
fList = dir(fullfile(path{8}, strcat(fileName, '*')));

% Find number of previous iterations
nIter = numel(fList);
clear fList

% Save all data
save(fullfile(path{8}, strcat(fileName, '_Iteration', num2str(nIter))));

