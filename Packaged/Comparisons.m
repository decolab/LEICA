%% ICA Extraction of Network States
%	This script computes the neural assemblies and assembly time courses
% from the data processed by the extraction script and computes the Shannon
% entropies of each assembly's time course.  In addition, it computes the
% total Shannon entropy of each condition.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% COMPARISON PIPELINE


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
loadFile = 'LEICA90_PhaseAssemblies_CIC';

% Load data
load(fullfile(path{8}, loadFile));
clear loadFile

% File to save
fileName = 'LEICA90_PhaseComparisons_CIC';


%% Reset paths overwritten by loaded file

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


%% 1) Find assemblies with significantly different activations between conditions

% Find assemblies with significant activations
for d = 1:N.datatype
	
	% Locate assemblies with significantly different activation distributions
	a = activations.(datatype{d})(2:N.condition+1,1);
	[sig.(datatype{d}).active, pval.(datatype{d}).active] = sigTest(a, pval.target, N.assemblies(d));
	N.assemblies(d, 2) = numel(sig.(datatype{d}).active{:,'FDR'});		% List number of significantly different assemblies
	
	% Locate assemblies with significantly different entropy distributions
	for k = 1:2
		a{k} = entro.(datatype{d}){k,1}';
	end
	[sig.(datatype{d}).entro, pval.(datatype{d}).entro] = sigTest(a, pval.target, N.assemblies(d));
	N.assemblies(d, 2) = numel(sig.(datatype{d}).entro.FDR);		% List number of significantly different assemblies
	
	% Compare assembly activation probabilities between conditions
	
	
	% Compare assembly metastabilities between conditions
	
    
	% Extract significant components and activations
	m = memberships.(datatype{d});
	memberships.(datatype{d}) = cell(1,2);
	memberships.(datatype{d}){1,1} = m;
	memberships.(datatype{d}){1,2} = m(:,sig.(datatype{d}).active{:,'FDR'});
	activations.(datatype{d}){1,2} = activations.(datatype{d}){1,1}(sig.(datatype{d}).active{:,'FDR'},:);
	
	% Extract significant activations, entropies, mean values per condition
	for c = 1:N.condition
		activations.(datatype{d}){c+1,2} = activations.(datatype{d}){c+1,1}(sig.(datatype{d}).active{:,'FDR'},:);
		entro.(datatype{d}){c,2} = entro.(datatype{d}){c,1}(:,sig.(datatype{d}).active{:,'FDR'});
		meanvalue.(datatype{d}){c,2} = meanvalue.(datatype{d}){c,1}(:,sig.(datatype{d}).active{:,'FDR'});
	end
end
clear c d a m


%% 2) Compare assemblywise entropy distributions between conditions

% Preallocate arrays
E.p = nan(N.datatype,1);
E.h = nan(N.datatype,1);

% Compare assemblywise entropy distributions between conditions
for d = 1:N.datatype
	% First row: total entropy per assembly
	% Second row: total entropy per subject
	% for s = 1:2
		for c = 1:N.condition
			% Get entropy distribution per condition
			E.sum{d,1}(c,:) = squeeze(sum(entro.(datatype{d}){c,2},1, 'omitnan'));
		end
		% Run ranksum comparison
		if ~isempty(E.sum{d,s})
			[E.p(d,1), E.h(d,1)] = ranksum(E.sum{d,1}(1,:), E.sum{d,1}(2,:));
		end
	% end
end
clear d c

% Convert results to tables
rLabel = {'Assembly', 'Subject'};
E.p = array2table(E.p, 'RowNames',datatype, 'VariableNames',rLabel);
E.h = array2table(E.h, 'RowNames',datatype, 'VariableNames',rLabel);


% Save interim results
% save(fullfile(path{8},fileName));


%% 3) Compare sum of Shannon entropies between conditions


%% 3) Plot results of Shannon entropy comparison

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
	title(strcat('Comparison: Total Entropy per ', rLabel{c}));
end
clear c

% Open new figure
F(N.fig) = figure(N.fig);
N.fig = N.fig+1;

% Plot p-values of assembly activations between conditions
for d = 1:N.datatype
    subplot(N.datatype, 1, d);
    bar(1:N.assemblies(d), pval.(datatype{d}));
    xlim([1, N.assemblies(d)]);
    hold on;
    % Mark significantly different assemblies
    plot(sig.(datatype{d}).ind, sig.(datatype{d}).bool(sig.(datatype{d}).bool > 0), '*k');
    plot(xlim, [0.05 0.05], '--r');	% plot 0.05 significance level
    ylabel('p-value');
    xlabel('Assembly Index');
    title({'Assembly Activation Significance', datatype{d}});
end

% Save figure(s)
savefig(F, fullfile(path{8},fileName), 'compact');


