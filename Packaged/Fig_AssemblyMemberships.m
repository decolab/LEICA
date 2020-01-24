%% This script plots the spatial assemblies in three formats:
%	1)	As a normalized bar plot, to show the influence of each region
%	2)	As a cortical network in MNI space
%	3)	As a connectivity matrix, one per assembly
% 
% NB: current template only supports 90-region data.  Need cerebellar
% template in order to plot 116-region data


%% Set paths

% Find general path (enclosing folder of current directory)
path{2,1} = strsplit(pwd, '/');
path{3,1} = strjoin(path{2}(1:end-2),'/');
path{4,1} = strjoin(path{2}(1:end-1),'/');
path{2,1} = strjoin(path{2}(1:end-3),'/');
path{1,1} = strjoin({path{2},'MATLAB','spm12'}, '/');

% Set required subdirectories
path{5,1} = fullfile(path{4},'Functions');
path{6,1} = fullfile(path{4},'Results');
path{7,1} = fullfile(path{3},'Data');

% Add required paths
addpath(genpath(path{1}));
addpath(path{5});
addpath(path{7});


%% Set files to use

% Set template file
cortex.path = fullfile(path{7},'MNI152_T1_2mm_brain_mask.nii');

% Set file names
loadFile = 'LEiDA90_Comparisons_OCDvHealth.mat';	% file to load
fileName = 'LEiDA90_Memberships_OCDvHealth';		% file storing figure(s)


%% Load data

% Load data
load(fullfile(path{6},loadFile), 'memberships', 'N');

% Extract field names from assembly vector
fnames = fieldnames(memberships);


%% Reset paths

% Find general path (enclosing folder of current directory)
path{2,1} = strsplit(pwd, '/');
path{3,1} = strjoin(path{2}(1:end-2),'/');
path{4,1} = strjoin(path{2}(1:end-1),'/');
path{2,1} = strjoin(path{2}(1:end-3),'/');
path{1,1} = strjoin({path{2},'MATLAB','spm12'}, '/');

% Set required subdirectories
path{5,1} = fullfile(path{4},'Functions');
path{6,1} = fullfile(path{4},'Results');
path{7,1} = fullfile(path{3},'Data');


%% Extract assembly membership vectors

V = cell(numel(fnames),1);
z = ones(numel(fnames),1);
for k = 1:numel(fnames)
	V{k} = memberships.(fnames{k}){1,2};		% extract assembly membership vectors
	if isempty(memberships.(fnames{k}){1,2})	% detect empty membership arrays (no significance)
		z(k) = 0;
	end
end
z = logical(z);
clear k

% Eliminate empty membership vectors
V = V(z,:);
fnames = fnames(z,:);

% Find number of significant assemblies per data type
N.assemblies = N.assemblies(z,2);
clear z

% Set order for bar plot
Order=[1:2:N.ROI N.ROI:-2:2];


%% Plot figure

% Run through data types
for f = 1:numel(fnames)
	
	% Declare figure
	F(f) = figure;
	colormap(jet);
	
	for c = 1:N.assemblies(f)

		% Sort and normalize assembly membership vectors
		Vc = V{f}(Order,c);
		Vc = Vc/max(abs(Vc));

		% Plot horizontal bar plots of assemblies
		subplot(5, N.assemblies(f), [c N.assemblies(f)+c 2*N.assemblies(f)+c])
		hold on
		barh(Vc.*(Vc<=0), 'FaceColor', [0.2  .2  1], 'EdgeColor', 'none', 'Barwidth', .5)
		barh(Vc.*(Vc>0), 'FaceColor', [1 .2 .2], 'EdgeColor', 'none', 'Barwidth', .5)
		ylim([0 91])
		xlim([-1 1])
		grid on
		set(gca, 'YTick', 1:N.ROI, 'Fontsize', 8)    
		set(gca, 'YTickLabel', [])
		title({[fnames{f}], ['Assembly ' num2str(c)]});

		% Plot networks in cortical space
		subplot(5, N.assemblies(f), 3*N.assemblies(f)+c)       
		plot_nodes_in_cortex(V{f}(:,c)', cortex);	% requires aal_cog.mat

		% Plot connectivity matrices
		subplot(5,N.assemblies(f),4*N.assemblies(f)+c);		
		FC_V = Vc*Vc';			% compute functional connectivity
		li=max(abs(FC_V(:)));	% 
		imagesc(FC_V,[-li li]);	colorbar;	% 
		axis square
		title(['V_C_' num2str(c) '.V_C_' num2str(c) '^T ']) 
		ylabel('Brain area #')
		xlabel('Brain area #')   
	end
end


%% Save figure
savefig(F, fullfile(path{6},fileName));
