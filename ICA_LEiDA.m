%% ICA Extraction of Network States
%	This script extracts the "activity states" of the network using
% independent component analysis.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DATA PIPELINE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%	SETUP

%% Set paths & filenames

% Shuffle random seed.  Necessary in array parallelization to avoid
% repeating same random seed across arrays.
rng('shuffle');

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = pwd;

% Set required subdirectories
path{4,1} = fullfile(path{2},'Data');
path{5,1} = fullfile(path{2},'Results');
path{6,1} = fullfile(path{1},'MATLAB');
path{7,1} = fullfile(path{3},'Functions');
path{8,1} = fullfile(path{3},'Results');

% % Set required subdirectories
% path{4,1} = fullfile(path{1},'Serotonin/shared');
% path{5,1} = fullfile(path{2},'Serotonin/knn');
% path{6,1} = fullfile(path{2},'Serotonin/sqdistance');
% path{7,1} = fullfile(path{2},'Serotonin/Non_parametric_Entropy_Estimator');

% Add relevant paths
for k = 6:numel(path)-1
	addpath(genpath(path{k}));
end
clear k


%% Set file names, analysis type

% Define list of patients
patientList = 'ClinicalData_OCD.xlsx';

% Define files containing BOLD data
loadDirectory = 'Subjects';
loadFiles = 'P*';

% File to save
fileName = 'LEICA';

% Set type of phase information to extract: eigenvector or synchrony
phaseType = 'LEICA';


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

% Get functional data
for k = 1:nFiles
	
	% Rename files
	D = strsplit(fList(k).name, '_');
	fName{k} = D{1};
	
	% Visualize data being extracted
	disp(['Extracting data for subject ', fName{k}, '.']);
	
	% Load data
	subjectBOLD(:,:,k) = load(fullfile(path{4}, loadDirectory, fList(k).name))';
	
end
clear D ans k


%% Set indices

% Extract indices for BOLD signals of patients, controls
I(:,1) = ~ismember(fName, patientData{:, 'Code'});
I(:,2) = ismember(fName, patientData{:, 'Code'});
I = array2table(I, 'VariableNames',{'Controls','OCD'});

% Set number of conditions
N.condition = size(I,2);

% Find number of subjects in each condition
for c = 1:N.condition
	N.subjects(c) = nnz(I{:,c});
end

% Extract labels for patients, controls
labels = cell(max(N.subjects), N.condition);
for c = 1:N.condition
	labels(1:nnz(I{:,c}), c) = fName(I{:,c});
end
clear c

% Save BOLD data, indices
save(fullfile(path{8},fileName));


%% Set parameters

% Temporal parameters
TR = 2;			% Repetition Time (seconds)
T.total(1) = T.scan*N.subjects(1);
T.total(2) = T.scan*N.subjects(2);

% Set bandpass filter
[bfilt, afilt] = narrowband(TR);

% Indices for vectorizing lower triangle
Isubdiag = find(tril(ones(N.ROI),-1));

% Set hypothesis test parameters
pval.target = 0.01;

% Set number of neighbors to search for in KNN
co = HShannon_kNN_k_initialization(1);

% Set figure counter
N.fig = 1;


%% Set pairwise comparisions (invalid for OCD because only OCD vs. controls)

% nfig=1;
% for NA=1:4
%     for NB=NA+1:4


%% SORTING

%% Extract data to cell array (TC_AAL)

% Preallocate data array
TS = cell(1,N.condition);
TS_AAL = cell(max(N.subjects), N.condition);

for c = 1:N.condition
	TS{c} = subjectBOLD(:,:,logical(I{:,c}));
	for s = 1:N.subjects(c)
		TS_AAL{s, c} = squeeze(TS{c}(:,:,s));
	end
end

% Save interim results
save(fullfile(path{8},fileName));



%% CONVERSION


%% 7) Convert data to z-score, phase data

% Preallocate variables to save FC patterns and associated information
T.all = zeros(2, sum(T.total));		% vector with subject nr and task at each t
t_all = 0;							% Index of time (starts at 0, updated until N.subjects*T.scan)

% Preallocate storage arrays
timeseriedata = nan(N.ROI, T.scan);
Z.timeseries = [];
if strcmpi(phaseType, 'LEiDA') ||  strcmpi(phaseType, 'Eigenvector')
	PH.timeseries = zeros(sum(T.total), N.ROI);				% All leading eigenvectors
else
	PH.timeseries = zeros(sum(T.total), length(Isubdiag));	% vectorized synchrony data
end

disp('Converting BOLD signal to z-score, phase synchrony data');

for c = 1:N.condition
	for s = 1:N.subjects(c)
		
		% Get the size of the BOLD signals from this subject in this task
		[N.ROI, T.scan] = size(TS_AAL{s,c});
		Phase_BOLD = zeros(N.ROI, T.scan);
		
		% Get the BOLD signals from this subject in this task
		BOLD = TS_AAL{s,c};
		BOLD = BOLD(:,1:T.scan);
		for seed=1:N.ROI
			BOLD(seed,:) = BOLD(seed,:) - mean(BOLD(seed,:));		% demean BOLD signal
			signal_filt = filtfilt(bfilt, afilt, BOLD(seed,:));		% bandpass filter BOLD signal
			timeseriedata(seed,:) = zscore(signal_filt);			% z-score BOLD signal
			Phase_BOLD(seed,:) = angle(hilbert(signal_filt));		% compute phase synchrony
		end
		
		% Instantaneous FC (BOLD Phase Synchrony) for each time point
		for t = 1:T.scan
			iPH = zeros(N.ROI, N.ROI);
			for n = 1:N.ROI
				for m = 1:N.ROI
					iPH(n,m) = cos(Phase_BOLD(n,t) - Phase_BOLD(m,t));
				end
			end
			t_all=t_all+1; % Update time

			if strcmpi(phaseType, 'LEiDA') ||  strcmpi(phaseType, 'Eigenvector')
				[V1,~] = eigs(iPH,1);
				PH.timeseries(t_all,:) = V1;
			else
				PH.timeseries(t_all,:) = iPH(Isubdiag);
			end
			T.all(:,t_all) = [s c]';	% Information that at t_all, V1 corresponds to subject s in a given task
		end

		Z.timeseries = horzcat(Z.timeseries, timeseriedata);

	end
end
clear Phase_BOLD BOLD iPH signal_filt timeseriedata m n nFiles loadFiles loadDirectory patientList s seed TR

PH.timeseries = PH.timeseries';

% Save interim results
save(fullfile(path{8},fileName));


%% 8) Convert z-score to point process

% = eventMat();


% Save interim results
save(fullfile(path{8},fileName));




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ASSEMBLY PIPELINE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ASSEMBLIES

%% 1) Compute the assemblies

disp('Processing the ICs from BOLD data')

N.assemblies(1,1) = NumberofIC(Z.timeseries);
N.assemblies(2,1) = NumberofIC(PH.timeseries);

% Compute assembly activations, timeseries
[assemblies.Z.activations{1}, assemblies.Z.memberships{1}, W.Z] = fastica(Z.timeseries, 'numOfIC',N.assemblies(1), 'verbose','off');
[assemblies.PH.activations{1}, assemblies.PH.memberships{1}, W.PH] = fastica(PH.timeseries, 'numOfIC',N.assemblies(2), 'verbose','off');

% Set assemblies, activations for BOLD signal
assemblies.BOLD.memberships{1} = ones(N.ROI);
assemblies.BOLD.activations{1} = Z.timeseries;

% Save interim results
save(fullfile(path{8},fileName));


%% COMPARISON

%% 3) Compute Shannon entropy for each each subject, condition

% Compute segmented entropy, activation distribution for each extraction type
[entro.BOLD, meanvalue.BOLD, distribution.BOLD] = HShannon_segmented(assemblies.BOLD.activations{1}, T.all, N.condition, N.subjects, N.ROI, co);
[entro.Z, meanvalue.Z, distribution.Z] = HShannon_segmented(assemblies.Z.activations{1}, T.all, N.condition, N.subjects, N.assemblies(1), co);
[entro.PH, meanvalue.PH, distribution.PH] = HShannon_segmented(assemblies.PH.activations{1}, T.all, N.condition, N.subjects, N.assemblies(2), co);

% Find number of data types
datatype = fieldnames(entro);
N.datatype = numel(datatype);

% Save interim results
save(fullfile(path{8},fileName));


%% 4) Compare Shannon entropies between conditions

% Preallocate arrays
E = cell(ndims(entro.BOLD)-1, N.datatype);
p = nan(ndims(entro.BOLD)-1, N.datatype);
h = nan(ndims(entro.BOLD)-1, N.datatype);

% Run comparisons between conditions
for c = 1:ndims(entro.BOLD)-1
	for d = 1:N.datatype
		% Get entropy distribution per condition
		%	First row: sum over subjects (total entropy per assembly)
		%	Second row: sum over assemblies (total entropy per subject)
		E{c,d} = squeeze(sum(entro.(datatype{d}), c+1, 'omitnan'));
        % E{c,d} = squeeze(cat(c+1, entro.(datatype{d})));
		
		% Run ranksum comparison
		[p(c,d), h(c,d), stats(c,d)] = ranksum(E{c,d}(1,:), E{c,d}(2,:));
	end
end
clear d c

% Convert results to tables
rLabel = {'Assembly', 'Subject'};
p = array2table(p, 'RowNames',rLabel, 'VariableNames',datatype);
h = array2table(h, 'RowNames',rLabel, 'VariableNames',datatype);


%% Compare other statistics between conditions

% Get metastability distribution for each condition & data type
% Get hierarchy distribution for each condition & data type
% Get functional complexity distribution for each condition & data type


%% 5) Find significantly different assemblies & assembly activations

% Find assemblies with significant activations
[sig.BOLD, pval.BOLD] = entroSigTest(distribution.BOLD, pval.target, N.ROI);
[sig.Z, pval.Z] = entroSigTest(distribution.Z, pval.target, N.assemblies(1));
[sig.PH, pval.PH] = entroSigTest(distribution.PH, pval.target, N.assemblies(2));

% Extract significant components, activations, 
for d = 1:N.datatype
	% Extract significant components and activations
	assemblies.(datatype{d}).memberships{1,2} = assemblies.(datatype{d}).memberships{1}(:,sig.(datatype{d}).bool);
	assemblies.(datatype{d}).activations{1,2} = assemblies.(datatype{d}).activations{1}(sig.(datatype{d}).bool,:);
	
	% Extract activations per condition
	for c = 1:N.condition
		% All activations
		assemblies.(datatype{d}).activations{2,1}{c,1} = distribution.(datatype{d}){c};
		
		% Significant activations
		assemblies.(datatype{d}).activations{2,2}{c,1} = distribution.(datatype{d}){c}(sig.(datatype{d}).bool,:);
	end
	
	% Extract significant entropies & mean values
	entro.(datatype{d}) = entro.(datatype{d})(:,:,sig.(datatype{d}).bool);
	meanvalue.(datatype{d}) = meanvalue.(datatype{d})(:,:,sig.(datatype{d}).bool);
end

% Save interim results
save(fullfile(path{8},fileName));


%% 6) Compute assembly correlation matrices

% Preallocate arrays
% distribution = cell(N.datatype, N.condition);
for d = 1:N.datatype
	for n = 1:numel(assemblies.(datatype{d}).memberships)
		assemblies.(datatype{d}).corr{1,n} = nan(size(assemblies.(datatype{d}).memberships{n}, 2));
		% assemblies.(datatype{d}).corr{1,2} = nan(size(assemblies.(datatype{d}).memberships.sig, 2));
		assemblies.(datatype{d}).corr{2,n} = nan(size(assemblies.(datatype{d}).activations{1,n}, 1));
		% assemblies.(datatype{d}).corr{2,2} = nan(size(assemblies.(datatype{d}).activations.sig, 1));
	end
	assemblies.(datatype{d}).corr{3,1} = nan(length(distribution.(datatype{d})), 1);
	assemblies.(datatype{d}).corr{3,2} = nan(length(sig.(datatype{d}).ind), 1);
end
L = {'Membership', 'Activation', 'Condition'};
A = nan(max(T.total), N.condition);
tvec = zeros(N.datatype,1);
clear d n

% Find data types with multiple significant assemblies, activation profiles
for d = 1:N.datatype
	for n = 1:numel(assemblies.(datatype{d}).memberships)
		tvec(d,n) = ~isempty(assemblies.(datatype{d}).memberships{n});
	end
	tvec(d,2) = (size(assemblies.(datatype{d}).memberships{2}, 2) > 1);
end
tvec = tvec(:,1) & tvec(:,2);
[tvec, ~] = find(tvec);


% Compute membership and activation correlation matrices
for d = 1:length(tvec)
    tv = tvec(d,:);
    
    % Generate assembly index vectors
	N.assemblies(tv-1, 2) = numel(sig.(datatype{tv}).ind);		
	
	% Assembly membership & activation correlation over all subjects
	for s = 1:2	% All asssemblies: s = 1.  Only significant assemblies: s = 2.
		
		% Compute assembly correlation over all subjects and conditions
		assemblies.(datatype{tv}).corr{1,s} = corr(assemblies.(datatype{tv}).memberships{s});		% membership
		assemblies.(datatype{tv}).corr{2,s} = corr(assemblies.(datatype{tv}).activations{1,s}');	% activations
		
		% Compute assembly activation correlation for each condition
		for ass = 1:N.assemblies(tv-1, s)
			% extract relevant assembly index
			% a = ind.(datatype{tv}){tv,s}(ass);
			% extract assembly activation matrix
			for c = 1:N.condition
				A(1:T.total(c), c) = assemblies.(datatype{tv}).activations{2,s}{c}(ass,:)';
			end
			% compute correlation matrix & extract pairwise correlation coefficient
			C = corr(A, 'rows','pairwise');
			assemblies.(datatype{tv}).corr{3,s}(ass) = C(1,2);
		end
	end
end
clear C A ass c d s n a


%% 7) Compute LEiDA connectivity matrix

% Find number of required components
N.comp.PH = nnz(sig.PH.bool);
N.comp.Z = nnz(sig.Z.bool);

% Preallocate LEiDA connectivity array
ICmatrix.PH = nan(N.ROI, N.ROI, N.comp.PH);

if strcmpi(phaseType, 'LEiDA') ||  strcmpi(phaseType, 'Eigenvector')
	for ass=1:N.comp.PH
		ICmatrix.PH(:,:,ass) = assemblies.PH.memberships{2}(:,ass)*(assemblies.PH.memberships{2}(:,ass))';
	end
	for ass=1:N.comp.Z
		ICmatrix.Z(:,:,ass) = assemblies.Z.memberships{2}(:,ass)*(assemblies.Z.memberships{2}(:,ass))';
	end
	clear ass
else
	numind = size(assemblies.PH.memberships{2}, 1);
	for ass=1:N.comp.PH
		for i=1:numind
			[ii, jj] = ind2sub([N.ROI N.ROI], Isubdiag(i));
			ICmatrix.PH(ii,jj,ass) = assemblies.PH.memberships{2}(i,ass);
			ICmatrix.PH(jj,ii,ass) = ICmatrix.PH(ii,jj,ass);
		end
	end
	
	numind = size(assemblies.Z.memberships{2}, 1);
	for ass=1:N.comp.Z
		for i=1:numind
			[ii, jj] = ind2sub([N.ROI N.ROI], Isubdiag(i));
			ICmatrix.Z(ii,jj,ass) = assemblies.Z.memberships{2}(i,ass);
			ICmatrix.Z(jj,ii,ass) = ICmatrix.Z(ii,jj,ass);
		end
	end
	clear numind ass
end

% Save interim results
save(fullfile(path{8},fileName));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% VISUALIZATION PIPELINE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 8) Plot results of Shannon entropy comparison

% Open new figure
F(N.fig) = figure(N.fig);
N.fig = N.fig+1;

% Plot Wilcoxon rank-sum results of comparing entropies between conditions
for c = 1:ndims(entro.BOLD)-1
	subplot(ndims(entro.BOLD)-1, 1, c);
	bar(table2array(p(c,:))); hold on;
	plot(xlim, [0.05 0.05], '--r');	% plot 0.05 significance level
	ylim([0 1]);
	ylabel('p-value');
	set(gca, 'xticklabel', datatype);
	title(strcat('comparison of total entropy per ', rLabel{c}));
end
clear c

% Open new figure
F(N.fig) = figure(N.fig);
N.fig = N.fig+1;

% Plot p-values of assembly activations between conditions & mark
% significantly different assemblies
subplot(3,1,1);
xlim([1, N.ROI]);
bar(1:N.ROI, pval.BOLD);
hold on;
plot(sig.BOLD.ind, sig.PH.bool(sig.BOLD.bool > 0), '*k');
plot(xlim, [0.05 0.05], '--r');	% plot 0.05 significance level
ylabel('p-value');
xlabel('Assembly Index');
title('BOLD Assembly Activation Significance');

subplot(3,1,2);
xlim([1, N.assemblies(1)]);
bar(1:N.assemblies(1), pval.Z);
hold on;
plot(sig.Z.ind, sig.Z.bool(sig.Z.bool > 0), '*k');
plot(xlim, [0.05 0.05], '--r');	% plot 0.05 significance level
ylabel('p-value');
xlabel('Assembly Index');
title('Z-Scored Assembly Activation Significance');

subplot(3,1,3);
xlim([1, N.assemblies(2)]);
bar(1:N.assemblies(2), pval.PH);
hold on;
plot(sig.PH.ind, sig.PH.bool(sig.PH.bool > 0), '*k');
plot(xlim, [0.05 0.05], '--r');	% plot 0.05 significance level
ylabel('p-value');
xlabel('Assembly Index');
title('Phase Assembly Activation Significance');


%% Plot assembly membership, activation correlations

% Plot overall membership, activation correlation matrices
for d = 1:length(tvec)
    tv = tvec(d);
    
    % Plot assembly activation, membership correlation matrices
    F(N.fig) = figure;
    N.fig = N.fig+1;
    for k = 1:size(assemblies.(datatype{tv}).corr,1)-1
        subplot(length(tvec), size(assemblies.(datatype{tv}).corr,1), k+2*(d-1));
        imagesc(assemblies.(datatype{tv}).corr{k}); colorbar;
        title(strcat(L{k}, ' Correlation for ', datatype{tv}));
    end
    
    % Plot assemblywise activation correlation between conditions
    F(N.fig) = figure;
    N.fig = N.fig+1;
    bar(assemblies.(datatype{tv}).corr{numel(assemblies.(datatype{tv}).corr)});
    xlabel('Assembly Index');
    title('Activation Correlation Between Conditions');
end
clear d tv ass c A C k


%% Save results

% Save interim results
savefig(F, fullfile(path{8},fileName), 'compact');

clear c d tv ass A C k phaseType t t_all F afilt bfilt tvec

% Save interim results
save(fullfile(path{8},fileName));


