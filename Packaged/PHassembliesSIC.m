%% ICA Extraction of Network States: Separate Assemblies
%	This script computes the neural assemblies and assembly time courses
% from the time series extracted from the BOLD data.  In addition, it
% computes the assemblywise Shannon entropy of each subject, and computes
% an upper bound of the total assemblywise entropy of each condition.
%	This version of the script calculates the assemblies of each condition
% separately.  This allows us to compare the membership of the resting
% state networks between conditions, but may make comparision of the
% activation time series between conditions more difficult to interpret.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ASSEMBLY PIPELINE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%	SETUP

%% Set paths & filenames

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
path{6,1} = fullfile(path{1},'MATLAB','FastICA');
path{7,1} = fullfile(path{2},'Functions','LEICA');
path{8,1} = fullfile(path{2},'Results','LEICA');

% Add relevant paths
for k = 6:numel(path)-1
	addpath(genpath(path{k}));
end
clear k


%% Set file names & load data

% Define files to load
loadFile = 'LEICA90_PhaseData';

% Load data
load(fullfile(path{8}, loadFile));

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_SIC_Assemblies');
clear loadFile S


%% Reset paths (in case original paths overwritten by loaded file)

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2},'Data');
path{5,1} = fullfile(path{2},'Results');
path{6,1} = fullfile(path{1},'MATLAB','FastICA');
path{7,1} = fullfile(path{2},'Functions','LEICA');
path{8,1} = fullfile(path{2},'Results','LEICA');


%% 1) Compute the assemblies

disp('Processing the ICs from BOLD data')

% Set metrics to calculate
clist = {'All','TimeSeries','Probability','Metastability','Entropy','Cohesion'};

% Set storage arrays
N.assemblies = nan(1, N.condition);
activations.cond = cell(1, N.condition);
activations.events = cell(1, N.condition);
activations.prob = cell(1, N.condition);
memberships.cond = cell(1, N.condition);
W.concat = cell(1, N.condition);
entro.cond = cell(1, N.condition);
meanvalue.cond = cell(1, N.condition);
cohesiveness = cell(1,N.condition);

for c = 1:N.condition
	disp(['Analyzing condition ', num2str(c)]);
	TI = T.index(2,:)==c;
	data = timeseries(:,TI);

	% Extract number of assemblies
	N.assemblies(c) = NumberofIC(data);

	% Compute assembly activation timecourses and memberships
	[activations.cond{c}, memberships.cond{c}, W.cond{c}] = fastica(data, 'numOfIC', N.assemblies(c), 'verbose','off');

	% Normalize memberships
	memberships.cond{c} = memberships.cond{c}./max(abs(memberships.cond{c}));

	% Compute activation entropy and timecourse for each assembly
	TI = T.index(1,TI);
	[entro.cond{c}, meanvalue.cond{c}, ~] = HShannon_segmented(data, TI, N.subjects(c), N.assemblies(c), co);

	% Preallocate arrays for subjectwise statistics
    activations.subj{c} = nan(N.assemblies(c), T.scan, N.subjects(c));
    activations.events{c} = activations.subj{c};
    cohesiveness{c} = nan(N.assemblies(c), N.subjects(c));
    
    % Compute subjectwise activations, events
    for nsub = 1:N.subjects(c)
        activations.subj{c}(:,:,nsub) = activations.cond{c}(:, (1+(nsub-1)*T.scan):(nsub*T.scan));
    	activations.events{c}(:,:,nsub) = activeMat(activations.subj{c}(:,:,nsub), 1);
    end
    
    % Compute assembly activation probability per subject
    activations.prob{c} = squeeze(sum(activations.events{c},2))./nnz(TI);
    
    % Compute regional cohesiveness per subject
    cohesiveness{c} = cohesive(memberships.cond{c}, activations.prob{c});
end
clear data TI c

% Save results
save(fullfile(path{8}, fileName));


