%% Visualization

% Visualize assembly activation magnitude histograms
metricFig(N.fig) = figure;
for ass = 1:N.assemblies
	subplot(3, ceil(N.assemblies/3), ass);
	for c = 1:N.conditions
		histogram(activities.cond.TS{c}(ass,:)); hold on;
	end
	legend('Controls','Patients');
	title(['Assembly ', num2str(ass), ' Activation Magnitudes']);
	xlabel('Magnitude');
	ylabel('Counts');
end
N.fig = N.fig + 1;
clear ass c

% Visualize means, standard deviations of assembly activation magnitudes
metricFig(N.fig) = figure;
bar(table2array(activities.cond.TSav)); hold on;
errorbar((1:N.assemblies)-0.15, activities.cond.TSav{:,'Controls'}, activities.cond.TSsd{:,'Controls'}, '.b');
errorbar((1:N.assemblies)+0.15, activities.cond.TSav{:,'Patients'}, activities.cond.TSsd{:,'Patients'}, '.r');
legend('Controls','Patients');
title('Mean Assembly Activation Magnitude');
xlabel('Assembly');
ylabel('Mean Activation Magnitudes');
N.fig = N.fig + 1;

% Visualize assembly activation time series
metricFig(N.fig) = figure;
for c = 1:N.conditions
	subplot(N.conditions, 1, c);
	bar(squeeze(activities.subj.TSav(:,:,c))); hold on;
	title([labels.Properties.VariableNames{c}(1:end-1), ' Mean Assembly Activation']);
	xlabel('Assembly');
	ylabel('Mean Activation Magnitudes');
end
N.fig = N.fig + 1;
clear c



%% Visualize Kuramoto order parameter

metricFig(N.fig) = figure;
% Conditions
subplot(2,4,[1 3]);
plot(1:max(T.condition), kuramoto.LEICA.cond); hold on;
xlim([0 max(T.condition)]);
legend(labels.Properties.VariableNames);
ylabel('Kuramoto Order Parameter');
xlabel('Time');
title('Kuramoto Order Parameter per Condition');
% Histogram
subplot(2,4,4);
for c = 1:N.conditions
	histogram(kuramoto.LEICA.cond(c,:), 'orientation', 'horizontal'); hold on;
end
legend(labels.Properties.VariableNames, 'Location', 'southeast');
xlabel('Counts');
% Subjects
ind = [5 6; 7 8];
for c = 1:N.conditions
	subplot(2,4, ind(c,:));
	plot(1:T.scan, squeeze(kuramoto.LEICA.subj(:,c,:))); hold on;
	xlim([0 T.scan]);
	ylabel('Kuramoto Order Parameter');
	xlabel('Time');
	title(labels.Properties.VariableNames{c});
end
clear s c ind
N.fig = N.fig + 1;


%% Visualize metastability

metricFig(N.fig) = figure;

% LEICA Conditions
subplot(2,2,1);
boxplot(table2array(metastable.LEICA.subj), 'Notch','on'); hold on;
% errorbar(1:N.conditions, table2array(metastable.LEICA.cond), std(table2array(metastable.LEICA.subj), 0, 1, 'omitnan'), '.r');
xticklabels(metastable.LEICA.cond.Properties.VariableNames);
ylabel('Metastability');
title('LEICA Metastability per Condition');
% ROI Conditions
subplot(2,2,2);
boxplot(table2array(metastable.AAL.subj), 'Notch','on'); hold on;
% errorbar(1:N.conditions, table2array(metastable.AAL.cond), std(table2array(metastable.AAL.subj), 0, 1, 'omitnan'), '.r');
xticklabels(metastable.LEICA.cond.Properties.VariableNames);
ylabel('Metastability');
title('AAL Metastability per Condition');
% Subjects
subplot(2,2,3);
for c = 1:N.conditions
	histogram(metastable.LEICA.subj{:,c}); hold on;
end
xlabel('Subject Metastabilities');
ylabel('Counts');
title('LEICA Metastabilities');
legend(metastable.LEICA.cond.Properties.VariableNames);
% ROI
subplot(2,2,4);
for c = 1:N.conditions
	histogram(metastable.AAL.subj{:,c}); hold on;
end
xlabel('Subject Metastabilities');
ylabel('Counts');
title('AAL Metastability');
legend(metastable.LEICA.cond.Properties.VariableNames);
clear c
N.fig = N.fig + 1;


%% Compute component-wise and global entropy

% Visualize entropies
metricFig(N.fig) = figure;
% Conditions
subplot(2,1,1);
boxplot(table2array(entro.subj), 'Notch','on'); hold on;
% errorbar(1:N.conditions, table2array(entro.LEICA.cond), std(table2array(entro.LEICA.subj), 0, 1, 'omitnan'), '.r');
xticklabels(entro.cond.Properties.VariableNames);
ylabel('Entropy');
title('LEICA Entropy per Condition');
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
N.fig = N.fig + 1;


%% Compute component-wise and global cohesion

% Preallocate storage
cohesion.cond = nan(N.ROI, N.conditions);
cohesion.subj = nan(N.ROI, max(N.subjects), N.conditions);

% Compute entropy for each condition and each subject
cohesion.concat = cohesive(memberships, events.concat.prob.av);
for c = 1:N.conditions
	cohesion.cond(:, c) = cohesive(memberships, events.cond.prob.av(:,c));
	for s = 1:N.subjects(c)
		cohesion.subj(:, s, c) = cohesive(memberships, events.subj.prob.av(:,s,c));
	end
end
clear c s

% Visualize entropies
metricFig(N.fig) = figure;
% Conditions
subplot(2,1,1);
bar(cohesion.cond); hold on;
legend(labels.Properties.VariableNames);
xticklabels(labelROI);
ylabel('Cohesion');
title('Global Cohesion per Condition');
% Subjects
subplot(2,1,2);
for c = 1:N.conditions
	histogram(cohesion.subj(:,:,c)); hold on;
end
xlabel('Subject Cohesion');
ylabel('Counts');
title('Subject Cohesion per Condition');
legend(labels.Properties.VariableNames);
clear c
N.fig = N.fig + 1;


%% Save results

% Get file list
fList = dir(fullfile(path{8}, strcat(fileName, '*')));

% Find number of previous iterations
nIter = numel(fList);
clear fList

% Edit fileName
fileName = strcat(fileName, '_Iteration', num2str(nIter));
clear nIter

% Save figures
savefig(metricFig, fullfile(path{8}, fileName), 'compact')
clear metricFig





%%


%%


%% Visualize Comparison Results


%% Visualize comparison between IC activation magnitudes

nFig = 2;
if ~isempty(sig.IC.FDR)
	set(metricFig(nFig),'defaultLegendAutoUpdate','off');
	cax = metricFig(nFig).CurrentAxes;
	mark = 1.3*ones(numel(sig.IC.FDR),1);
	plot(cax, sig.IC.FDR, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax

% Visualize ICs with significantly different time series
N.fig = N.fig+1;
metricFig(N.fig) = figure; hold on;
colormap jet
for i = 1:length(sig.IC.FDR)
	subplot(3, ceil(length(sig.IC.FDR)/3), i);
	imagesc(ICs(:,:,sig.IC.FDR(i))); colorbar;
	title(['Significant Component ', num2str(i)]);
	xlabel('AAL Region');
	ylabel('AAL Region');
end
clear i

% Plot IC activations over time
N.fig = N.fig+1;
metricFig(N.fig) = figure; hold on;
d = 1:ceil(length(sig.TS.FDR)/3);
d = [d, d+2*ceil(length(sig.TS.FDR)/3), d+4*ceil(length(sig.TS.FDR)/3)];
for i = 1:length(sig.TS.FDR)
	subplot(3*2, ceil(length(sig.TS.FDR)/3), d(i));
	scatter(T.TR*(1:T.condition(1)), activities.cond.TS{1}(sig.TS.FDR(i),:), 'filled');
	title(['IC ', num2str(sig.TS.FDR(i)), ' Activations']);
	ylabel('Control');
	
	subplot(3*2, ceil(length(sig.TS.FDR)/3), d(i)+ceil(length(sig.TS.FDR)/3));
	scatter(T.TR*(1:T.condition(2)), activities.cond.TS{2}(sig.TS.FDR(i),:), 'filled');
	ylabel('OCD');
end
clear d i

% Plot values of IC activations
N.fig = N.fig+1;
metricFig(N.fig) = figure; hold on;
for i = 1:length(sig.TS.FDR)
	a = nan(1, max(T.condition));
	a(1:min(T.condition)) = activities.cond.TS{1}(sig.TS.FDR(i),:);
	act = cat(1, a, activities.cond.TS{2}(sig.TS.FDR(i),:));
	
	subplot(3, ceil(length(sig.TS.FDR)/3), i);
	boxplot(act', 'Notch','on', 'Labels',{'Controls','Patients'});
	title(['IC ', num2str(sig.TS.FDR(i)), ' Activation']);
end
clear a act


%% Compare mean component activation magnitudes

% Visualize results
figure; hold on;
bar(sig.TS.p);
scatter(find(sig.TS.Sidak), 0.5.*ones(length(find(sig.TS.Sidak)),1), 'r*');
xlabel('ICs');
ylabel('p-value');
legend('p-value','Sidak significance');
title('IC Time Series p-values: OCD vs. Controls');

nFig = 1;
ind = find(sig.active.FDR);
mark = ones(numel(ind),1);
if ~isempty(ind)
	cax = metricFig(nFig).CurrentAxes;
	scatter(cax, ind, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax


%% Compare metastability between states

% Run tests for each assembly
con = metastable.LEICA.subj{:,'Controls'}(isfinite(metastable.LEICA.subj{:,'Controls'}));
pat = metastable.LEICA.subj{:,'OCD'}(isfinite(metastable.LEICA.subj{:,'OCD'}));
if adtest(con) && adtest(pat)
	[sig.meta.h, sig.meta.p] = kstest2(con, pat, 'Alpha',pval.target);
else
	[sig.meta.h, sig.meta.p] = ttest2(con, pat, 'Alpha',pval.target);
end
clear con pat

% Visualize results
nFig = nFig+1;
ind = find(sig.meta.h);
mark = ones(numel(ind),1);
if ~isempty(ind)
	cax = metricFig(nFig).CurrentAxes;
	scatter(cax, ind, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax


%% Compare entropy between states

% Visualize results
nFig = nFig+1;
ind = find(sig.entro.h);
mark = ones(numel(ind),1);
if ~isempty(ind)
	cax = metricFig(nFig).CurrentAxes;
	scatter(cax, ind, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax


%% Compare cohesion between states

% Preallocate storage arrays
sig.cohesion.h = nan(N.ROI, 1);
sig.cohesion.p = nan(N.ROI, 1);

% Run tests for each assembly
for roi = 1:N.ROI
	con = cohesion.LEICA.cond(roi,1); con = con(isfinite(con));
	pat = cohesion.LEICA.cond(roi,2); pat = pat(isfinite(pat));
	[sig.cohesion.h(roi), sig.cohesion.p(roi)] = kstest2(con, pat, 'Alpha',pval.target);
end
clear con pat roi

% Run FDR correction
sig.cohesion.FDR = FDR_benjHoch(sig.cohesion.p, pval.target);

% Run Bonferroni multiple comparison correction
sig.cohesion.Bonferroni = (sig.cohesion.p < (pval.target/N.assemblies));

% Run Dunn-Sidak multiple comparison correction
alpha = 1-(1-pval.target)^(1/N.assemblies);
sig.cohesion.Sidak = (sig.cohesion.p < alpha);
clear alpha

% Visualize results
nFig = nFig+1;
ind = find(sig.cohesion.FDR);
mark = ones(numel(ind),1);
if ~isempty(ind)
	cax = metricFig(nFig).CurrentAxes;
	scatter(cax, ind, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax





