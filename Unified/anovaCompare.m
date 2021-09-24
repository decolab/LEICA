%% Import data

clear; close all; clc;

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2}, 'UCLA', 'Data');
path{5,1} = fullfile(path{2}, 'UCLA', 'Results','LEICA');
path{6,1} = fullfile(path{1}, 'Project','Atlases','AAL');

% Add relevant paths
fpath{1,1} = fullfile(path{1},'MATLAB','BrainNetViewer');
fpath{2,1} = fullfile(path{1},'MATLAB','permutationTest');
fpath{3,1} = fullfile(path{1},'Project','Functions');
fpath{4,1} = fullfile(path{1},'Project','LEICA','Functions');
for k = 1:numel(fpath)
	addpath(genpath(fpath{k}));
end
clear fpath k
addpath(fullfile(path{1},'MATLAB','spm12'));

% Determine whether to z-scale membership weights
zthresh = 1:0.2:1.6;
zscale = false;

% Determine which files to compare
band = 'WideBand';
k = 1;
iteration = 2;
fname = strcat('_',band, '_k',num2str(k), '_Iteration',num2str(iteration));
pfix = {'LE_ICA_AAL90_CIC_COS', 'M_ICA_AAL90_CIC_EXP'};
clear assemblies distance band k iteration

% Set decompositions, spaces, comparisions
nFig = 1;
titles = {'LEICA', 'MICA'};	% compression type
spaces = {'dFC' 'IC'};				% space in which to compare
dim = {'Subject', 'Component'};		% dimensions to compare
ttype = {'kstest2', 'permutation'};	% set test types to run
pTarget = 0.05;						% target p-value
prange = 0.025:0.025:0.1;			% Set range of p-values to test

% Load metadata
load(fullfile(path{4}, 'formattedUCLA'), 'metadata', 'I', 'ROI', 'badSubj');
labels = I.Properties.VariableNames;

% Load entropies
N = cell(numel(pfix), numel(spaces));
entro = cell(numel(pfix), numel(spaces));
FCD = cell(numel(pfix), numel(spaces));
memberships = cell(numel(pfix), 1);	% have 2 dFC / FCD spaces: ROI and IC.
for f = 1:numel(pfix)
	e = load(fullfile(path{5}, 'CommonIC', strcat(pfix{f}, fname)), 'entro','N','memberships','FCD','I','T','ROI');
	N{f,1} = e.N; N{f,1}.IC = N{f,1}.ROI;   % N{f,1} = rmfield(N{f,1}, 'ROI');
	% N{f,2} = N{f,1};
	N{f,2} = e.N;   % N{f,2} = rmfield(N{f,2}, 'ROI');
	% entro{f,1} = e.entro.BOLD;
	entro{f,1} = e.entro.dFC;
	entro{f,2} = e.entro.IC;
	memberships{f,1} = e.memberships;
	FCD{f,1} = e.FCD.dFC;
	FCD{f,2} = e.FCD.IC;
end
T = e.T;
clear e f

% Load network labels
labels_ROI = string(ROI.Properties.RowNames);
% labels_ROI = load(fullfile(path{5}, 'AAL_labels'));
% labels_ROI = string(labels_ROI.label90);
% labels_ROI = strip(LR_version_symm(labels_ROI));

% Generate MNI coordinates
coords_ROI = horzcat(ROI{:,'x'}, ROI{:,'y'}, ROI{:,'z'});
% coords_ROI = load(fullfile(path{5}, 'aal_cog.txt'), 'aal_cog');
% coords_ROI = LR_version_symm(coords_ROI);
origin = [65 45.5 35];							% center origin
MNIscale = 5.5/10;								% scale MNI coordinates
sphereScale = 2.5;								% scale sphere size
coords_ROI = MNIscale*coords_ROI;			% scale MNI coordinates
clear label90 MNIscale

% Set brain rendering parameters
cortex.file = fullfile(path{6},'MNI152_T1_2mm_brain_mask.nii');	% file containing cortical atlas
cortex.color = [0.9 0.9 0.9];					% color for cortical rendering
cortex.transparency = 0.1;						% set to 1 for opaque cortex
cortex.val = 0.3;								% set isonormal line spacing
cortex.view = [-90 90];							% set camera angle
rdux = 0.7;						% proportion of surface faces to keep (<1)


%% Set up variables, group matrices for repeated measures ANOVA

e = nan(sum(N{1,2}.subjects), N{1,2}.IC);
group = cell(sum(N{1,2}.subjects),1);
i = [0 N{1,2}.subjects];
for c = 1:N{1,2}.conditions
    e(1+sum(i(1:c)):sum(i(1:c+1)),:) = entro{1,2}(:,1:N{1,2}.subjects(c),c)';
    group(1+sum(i(1:c)):sum(i(1:c+1))) = labels(c);
end
clear i c
group = string(group);

% Convert to table
t = table(group,e(:,1),e(:,2),e(:,3),e(:,4),e(:,5),e(:,6),e(:,7),e(:,8),e(:,9),e(:,10), 'VariableNames',{'Group','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10'});
ic = table([1:10]', 'VariableNames',{'Components'});

% fit repeated measures model
rm = fitrm(t, 'C1-C10~Group', 'WithinDesign',ic);

% run repeated measures ANOVA
atbl = anova(rm);
ratbl = ranova(rm);
[matbl, A,C,D] = manova(rm, 'By','Group');
clear group


%% Set up n-way ANOVA

% Set N-way ANOVA groups: component & conditions
key{1} = repmat([1:N{1,2}.IC]', [1 size(entro{1,2},2) N{1,2}.conditions]);
k = nan(1,1,N{1,2}.conditions);
k(:,:,:) = 1:N{1,2}.conditions;
key{2} = repmat(k, [N{1,2}.IC size(entro{1,2},2) 1]);
clear k

% Reorder entropy values & keys
e = reshape(entro{1,2}, [numel(entro{1,2}) 1]);
group{2} = reshape(key{1}, [numel(key{1}) 1]);
group{1} = reshape(key{2}, [numel(key{2}) 1]);
group{2}(isnan(e)) = [];
group{1}(isnan(e)) = [];
e(isnan(e)) = [];
group{2} = string(group{2});
group{1} = string(labels(group{1}));
clear key

% Set N-way ANOVA groups: subjects
% group{3} = reshape(repmat(string(I.Properties.RowNames)', [N{1,2}.IC, 1]), [numel(group{1}) 1]);

% Set model terms & nesting matrices
terms = [1 0; 0 1; 1 1];
nest = [0 0 0; 0 0 0; 0 1 0];

% Run N-way ANOVA
[p, tbl, stats] = anovan(e, group, 'varnames',{'Diagnosis' 'Component'}, 'display','on', 'model','interaction');% terms);

% Multiple comparison correction
C = cell(size(stats.varnames));
M = cell(size(stats.varnames));
F = cell(size(stats.varnames));
groups = cell(size(stats.varnames));
for d = numel(stats.varnames):-1:1
    [C{d}, M{d}, F{d}, groups{d}] = multcompare(stats, 'Dimension',d);
end
clear e


%% Set up Kruskal-Wallis, Friedman tests

% [p, tbl, stats] = kruskalwallis();
% [p, tbl, stats] = freidman();