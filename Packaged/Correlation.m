%% ICA Extraction of Network States
%	This script computes the neural assemblies and assembly time courses
% from the data processed by the extraction script and computes the Shannon
% entropies of each assembly's time course.  In addition, it computes the
% total Shannon entropy of each condition.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CORRELATION PIPELINE


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
path{6,1} = fullfile(path{1},'MATLAB', 'spm12');
path{7,1} = fullfile(path{3},'Functions');
path{8,1} = fullfile(path{2},'Results', 'LEICA');

% Add relevant paths
addpath(path{6});
addpath(genpath(path{7}));


%% Set file names & load data

% Define files to load
loadFile = 'LEiDA90_Comparisons_OCDvHealth';

% Load data
load(fullfile(path{8}, loadFile));
clear loadFile

% File to save
fileName = 'LEiDA90_Correlations_OCDvHealth';


%%

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2},'Data');
path{5,1} = fullfile(path{2},'Results');
path{6,1} = fullfile(path{1},'MATLAB', 'spm12');
path{7,1} = fullfile(path{3},'Functions');
path{8,1} = fullfile(path{2},'Results', 'LEICA');


%% 6) Preallocate assembly correlation matrices

% for d = 1:N.datatype
% 	for n = 1:size(memberships.(datatype{d}),2)
% 		memberships.(datatype{d}){2,n} = nan(size(memberships.(datatype{d}){1,n}, 2));
% 		activations.(datatype{d}){3,n} = nan(size(activations.(datatype{d}){1,n}, 1));
% 	end
% 	assemblies.(datatype{d}).corr{3,1} = nan(length(distribution.(datatype{d})), 1);
% 	assemblies.(datatype{d}).corr{3,2} = nan(length(sig.(datatype{d}).ind), 1);
% end
% L = {'Membership', 'Activation', 'Condition'};
% A = nan(max(T.condition), N.condition);
% tvec = zeros(N.datatype,1);
% clear d n


%% Find data types with multiple significant assemblies

% Find data with significant assemblies
for d = 1:N.datatype
	for n = 1:numel(memberships.(datatype{d}))
		tvec(d,n) = ~isempty(memberships.(datatype{d}){n});
	end
	tvec(d,2) = (size(memberships.(datatype{d}){1,2}, 2) > 1);
end
clear d n

% Confirm that data possess multiple significant assemblies
tvec = tvec(:,1) & tvec(:,2);
[tvec, ~] = find(tvec);


%% Compute membership and activation correlation matrices

for d = 1:length(tvec)
    
	% Extract index of interest
	tv = tvec(d,:);
	
	% Assembly membership & activation correlation over all subjects
	% All asssemblies: s = 1.  Only significant assemblies: s = 2.
	for s = 1:size(memberships.(datatype{tv}),2)
		
		% Compute assembly correlation over all subjects and conditions
		memberships.(datatype{tv}){2,s} = corr(memberships.(datatype{tv}){1,s});	% membership
		activations.(datatype{tv}){3,s} = corr(activations.(datatype{tv}){1,s}');	% activations
		
		% Compute assembly activation correlation between conditions
		for ass = 1:N.assemblies(tv, s)
			for c = 1:N.condition
				A(1:T.condition(c), c) = activations.(datatype{tv}){2,s}{c}(ass,:)';
			end
			% compute correlation matrix & extract pairwise correlation coefficient
			C = corr(A, 'rows','pairwise');
			activations.(datatype{tv}){4,s}(ass) = C(1,2);
		end
	end
end
clear C A ass c d s n a


%% 7) Compute assembly connectivity matrices

% Determine how to compute connectivity matrix
if strcmpi(phaseType, 'LEiDA') ||  strcmpi(phaseType, 'Eigenvector')
    for d = 1:numel(tvec)
		
		% find correlation matrices to compute
		tv = tvec(d);
		
        % Preallocate LEiDA connectivity array
        ICmatrix.(datatype{tv}) = nan(N.ROI, N.ROI, N.assemblies(tv,2));
        
        % Compute connectivity matrix
        for ass = 1:N.assemblies(tv,2)
            ICmatrix.(datatype{tv})(:,:,ass) = memberships.(datatype{tv}){1,2}(:,ass)*(memberships.(datatype{tv}){1,2}(:,ass))';
        end
    end
	clear ass
else
    for d = 1:numel(tvec)
		
		% find correlation matrices to compute
		tv = tvec(d);
		
        % Preallocate LEiDA connectivity array
        ICmatrix.(datatype{tv}) = nan(N.ROI, N.ROI, N.assemblies(tv,2));
        
        % Compute connectivity matrix
        numind = size(memberships.(datatype{tv}){1,2}, 1);
        for ass=1:N.assemblies(tv,2)
            for i=1:numind
                [ii, jj] = ind2sub([N.ROI N.ROI], Isubdiag(i));
                ICmatrix.(datatype{d})(ii,jj,ass) = memberships.(datatype{d}){1,2}(i,ass);
                ICmatrix.(datatype{d})(jj,ii,ass) = ICmatrix.(datatype{d})(ii,jj,ass);
            end
        end
    end
	clear numind
end
clear tv ass

% Save interim results
save(fullfile(path{8},fileName));


%% Plot assembly membership, activation correlations

% Declare title
L = {'All', 'Signficant'};

% Plot overall membership, activation correlation matrices
for d = 1:length(tvec)
	
	% Extract data type of interest
	tv = tvec(d);
	
	% Declare figure
	F(N.fig) = figure;
    N.fig = N.fig+1;
    
	% Repeat for each condition
	for s = 1:N.condition
		% Plot assembly membership correlation
		subplot(3,2,s);
		imagesc(memberships.(datatype{tv}){2,s}); colorbar;
		title({[datatype{tv}], ' Membership Correlation', [L{s}, ' Assemblies']});
	
		% Plot assembly activation correlations over all time
		subplot(3,2, 2+s);
		imagesc(activations.(datatype{tv}){3,s}); colorbar;
		title({'Activation Correlation', [L{s}, ' Assemblies']});
		xlabel('Assembly Index');
	
		% Plot assembly activation correlation across conditions
		subplot(3,2, 4+s);
		b = bar(activations.(datatype{tv}){4,s});
		if s == 2
			xticks(sig.(datatype{tv}).ind);
		end
		xlabel('Assembly Index');
		title({'Correlation over Conditions', [L{s}, ' Assemblies']});
	end
end
clear d tv k


%% Plot connectivity matrices for each data type

% Find connectivity matrices to plot
% dt = fieldnames(ICmatrix);

for f = 1:numel(tvec)
	
	% Extract data type of interest
	tv = tvec(f);
    
    % Plot figure
    F(N.fig) = figure;
    N.fig = N.fig+1;
	
    % Plot each assembly
	for ass = 1:size(ICmatrix.(datatype{tv}), 3)
		subplot(1, size(ICmatrix.(datatype{tv}), 3), ass)
		imagesc(ICmatrix.(datatype{tv})(:,:,ass)); colorbar;
		title({strcat(datatype{tv}, ' Connectivity'), strcat('Assembly ', num2str(ass))});
	end
end
clear f ass tv tvec s

% Save figure(s) results
savefig(F, fullfile(path{8},fileName), 'compact');
