%% Visualization of LEICA Analysis
%	This script generates visual representations of the findings from the
% LEICA script.  This version is specifically intended to represent
% findings from the common independent component (CIC) version of LEICA, as
% this is the version in which component activations can be directly
% compared between groups.
%	As LEICA has several phases, each of which generate separate metrics,
% this script has multiple sections.  Each section corresponds to one
% section of the LEICA script.  It is hoped that this segementation will
% allow for greater readability and flexibilty in the use of this script.



%% SETUP

clear; close all; clc

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2}, 'OCD', 'Data');
path{5,1} = fullfile(path{2}, 'OCD', 'Results');
path{6,1} = fullfile(path{2}, 'OCD', 'Results','LEICA');

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

% Load data
load(fullfile(path{4}, 'sc90.mat'));
load(fullfile(path{6}, 'LEICA90_CIC_COS_WideBand_k1_Iteration2.mat'));

% Set N.fig
N.fig = 1;



%% dFC EXTRACTION

% Visualize interim results
F(N.fig) = figure;
N.fig = N.fig+1;
subplot(1,3,1); imagesc(sc90); title('Structural');
subplot(1,3,2); imagesc(squeeze(mean(FC(:,:,:,1),3,'omitnan'))); title('Mean FC: Patients');
subplot(1,3,3); imagesc(squeeze(mean(FC(:,:,:,2),3,'omitnan'))); title('Mean FC: Controls');



%% COMPUTING ICs FROM dFC

% Visualize assembly activation magnitude histograms
F(N.fig) = figure;
for ass = 1:N.IC
	subplot(3, ceil(N.IC/3), ass);
	for c = 1:N.conditions
		histogram(activities.cond{c}(ass,:)); hold on;
	end
	legend('Controls','Patients');
	title(['Assembly ', num2str(ass), ' Activation Magnitudes']);
	xlabel('Magnitude');
	ylabel('Counts');
end
N.fig = N.fig + 1;
clear ass c

% Visualize means, standard deviations of IC activation magnitudes per condition
F(N.fig) = figure; N.fig = N.fig + 1;
bar(table2array(activities.av.cond)); hold on;
errorbar((1:N.IC)-0.15, activities.av.cond{:,'Controls'}, activities.sd.cond{:,'Controls'}, '.b');
errorbar((1:N.IC)+0.15, activities.av.cond{:,'Patients'}, activities.sd.cond{:,'Patients'}, '.r');
legend('Controls','Patients');
title('Mean Assembly Activation Magnitude');
xlabel('Assembly');
ylabel('Mean Activation Magnitudes');
if ~isempty(sig.IC.FDR)
	set(F(nFig),'defaultLegendAutoUpdate','off');
	cax = F(nFig).CurrentAxes;
	mark = 1.3*ones(numel(sig.IC.FDR),1);
	plot(cax, sig.IC.FDR, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax

% Boxplot of IC activation magnitudes per condition
F(N.fig) = figure; hold on; N.fig = N.fig+1;
for i = 1:length(sig.TS.FDR)
	a = nan(1, max(T.condition));
	a(1:min(T.condition)) = activities.cond{1}(sig.IC.TS.FDR(i),:);
	act = cat(1, a, activities.cond.TS{2}(sig.IC.TS.FDR(i),:));
	subplot(3, ceil(length(sig.TS.FDR)/3), i);
	boxplot(act', 'Notch','on', 'Labels',{'Controls','Patients'});
	title(['IC ', num2str(sig.TS.FDR(i)), ' Activation']);
end
clear a act

% Visualize means, standard deviations of IC activation magnitudes per subject
F(N.fig) = figure; N.fig = N.fig + 1;
for c = 1:N.conditions
	subplot(N.conditions, 1, c);
	bar(squeeze(activities.av.subj(:,:,c))); hold on;
	title([labels.Properties.VariableNames{c}(1:end-1), ' Mean Assembly Activation']);
	xlabel('Assembly');
	ylabel('Mean Activation Magnitudes');
end
clear c

% Visualize time series p-values; AAL vs. ICs
figure; hold on; N.fig = N.fig+1;
subplot(1,2,1);
bar(sig.IC.TS.p);
scatter(find(sig.IC.TS.Sidak), 0.5.*ones(length(find(sig.IC.TS.Sidak)),1), 'r*');
xlabel('ICs');
ylabel('p-value');
legend('p-value','Sidak significance');
title('IC Time Series p-values: OCD vs. Controls');
subplot(1,2,2);
bar(sig.AAL.TS.p);
scatter(find(sig.AAL.TS.Sidak), 0.5.*ones(length(find(sig.AAL.TS.Sidak)),1), 'r*');
xlabel('ICs');
ylabel('p-value');
legend('p-value','Sidak significance');
title('AAL Time Series p-values: OCD vs. Controls');

% 
F(N.fig) = figure; hold on; N.fig = N.fig+1;
d = 1:ceil(length(sig.TS.FDR)/3);
d = [d, d+2*ceil(length(sig.TS.FDR)/3), d+4*ceil(length(sig.TS.FDR)/3)];
for i = 1:length(sig.TS.FDR)
	subplot(3*2, ceil(length(sig.TS.FDR)/3), d(i));
	scatter(T.TR*(1:T.condition(1)), activities.cond{1}(sig.TS.FDR(i),:), 'filled');
	title(['IC ', num2str(sig.TS.FDR(i)), ' Activations']);
	ylabel('Control');
	
	subplot(3*2, ceil(length(sig.TS.FDR)/3), d(i)+ceil(length(sig.TS.FDR)/3));
	scatter(T.TR*(1:T.condition(2)), activities.cond.TS{2}(sig.TS.FDR(i),:), 'filled');
	ylabel('OCD');
end
clear d i

% Visualize ICs with significantly different time series
F(N.fig) = figure; hold on; N.fig = N.fig+1;
colormap jet
for i = 1:length(sig.IC.FDR)
	subplot(3, ceil(length(sig.IC.FDR)/3), i);
	imagesc(ICs(:,:,sig.IC.FDR(i))); colorbar;
	title(['Significant Component ', num2str(i)]);
	xlabel('AAL Region');
	ylabel('AAL Region');
end
clear i



%% IC TIME SERIES METRICS

% Kuramoto Order Parameter
F(N.fig) = figure;
N.fig = N.fig + 1;
% Conditions
subplot(2,4,[1 3]);
plot(1:max(T.condition), kuramoto.cond); hold on;
xlim([0 max(T.condition)]);
legend(labels.Properties.VariableNames);
ylabel('Kuramoto Order Parameter');
xlabel('Time');
title('Kuramoto Order Parameter per Condition');
% Histogram
subplot(2,4,4);
for c = 1:N.conditions
	histogram(kuramoto.cond(c,:), 'orientation', 'horizontal'); hold on;
end
legend(labels.Properties.VariableNames, 'Location', 'southeast');
xlabel('Counts');
% Subjects
ind = [5 6; 7 8];
for c = 1:N.conditions
	subplot(2,4, ind(c,:));
	plot(1:T.scan, squeeze(kuramoto.subj(:,c,:))); hold on;
	xlim([0 T.scan]);
	ylabel('Kuramoto Order Parameter');
	xlabel('Time');
	title(labels.Properties.VariableNames{c});
end
clear s c ind


% Metastability
F(N.fig) = figure;
N.fig = N.fig + 1;
subplot(1,2,1);
boxplot(table2array(metastable.subj), 'Notch','on'); hold on;
xticklabels(metastable.LEICA.cond.Properties.VariableNames);
ylabel('Metastability');
title('LEICA Metastability per Condition');
ind = find(sig.metastable.h);
mark = ones(numel(ind),1);
if ~isempty(ind)
	cax = F(nFig).CurrentAxes;
	scatter(cax, ind, mark, '*', 'MarkerEdgeColor','k');
end
% Subjects
subplot(1,2,2);
for c = 1:N.conditions
	histogram(metastable.subj{:,c}); hold on;
end
xlabel('Subject Metastabilities');
ylabel('Counts');
title('LEICA Metastabilities');
legend(metastable.LEICA.cond.Properties.VariableNames);
clear c ind mark cax


% Entropy
F(N.fig) = figure;
N.fig = N.fig + 1;
% Conditions
subplot(2,1,1);
boxplot(table2array(entro.subj), 'Notch','on'); hold on;
xticklabels(entro.cond.Properties.VariableNames);
ylabel('Entropy');
title('IC Entropy per Condition');
ind = find(sig.entro.h);
mark = ones(numel(ind),1);
if ~isempty(ind)
	cax = F(nFig).CurrentAxes;
	scatter(cax, ind, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax
% Subjects
subplot(2,1,2);
for c = 1:N.conditions
	histogram(entro.subj{:,c}); hold on;
end
xlabel('Subject Entropy');
ylabel('Counts');
title('Subject Entropy per Condition');
legend(entro.cond.Properties.VariableNames);
clear c



%% Save results

% % Get file list
% fList = dir(fullfile(path{8}, strcat(fileName, '*')));
% 
% % Find number of previous iterations
% nIter = numel(fList);
% clear fList
% 
% % Edit fileName
% fileName = strcat(fileName, '_Iteration', num2str(nIter));
% clear nIter
% 
% % Save figures
% savefig(F, fullfile(path{8}, fileName), 'compact')
% clear metricFig


