%% ICA Common Assemblies: Computing Component Metrics
%	


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ASSEMBLY PIPELINE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%	SETUP

%% Set paths & filenames

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
path{6,1} = fullfile(path{2},'Functions','LEICA');
path{7,1} = fullfile(path{2},'Results','LEICA');

% Add relevant paths
addpath(genpath(path{6}));


%% Set file names & load data

% Define files to load
loadFile = 'LEICA90_CIC_Assemblies';

% Load data
load(fullfile(path{7}, loadFile));

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_', S{2}, '_Metrics');
clear loadFile S

% Reset N.fig
N.fig = 1;


%% Reset paths (in case original paths overwritten by loaded file)

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2},'Data');
path{5,1} = fullfile(path{2},'Results');
path{6,1} = fullfile(path{2},'Functions','LEICA');
path{7,1} = fullfile(path{2},'Results','LEICA');


%% Compute component activation magnitudes

% Activation magnitude means, standard deviations
activities.cond.TSav = nan(N.assemblies, N.conditions);
activities.cond.TSsd = nan(N.assemblies, N.conditions);
activities.subj.TSav = nan(N.assemblies, max(N.subjects), N.conditions);
activities.subj.TSsd = nan(N.assemblies, max(N.subjects), N.conditions);
for c = 1:N.conditions
	activities.cond.TSav(:,c) = mean(activities.cond.TS{c}, 2, 'omitnan');
	activities.cond.TSsd(:,c) = std(activities.cond.TS{c}, 0, 2, 'omitnan');
	for s = 1:N.subjects(c)
		activities.subj.TSav(:,s,c) = mean(activities.subj.TS{s,c}, 2, 'omitnan');
		activities.subj.TSsd(:,s,c) = std(activities.subj.TS{s,c}, 0, 2, 'omitnan');
	end
end
clear c s
activities.cond.TSav = array2table(activities.cond.TSav, 'VariableNames',labels.Properties.VariableNames);
activities.cond.TSsd = array2table(activities.cond.TSsd, 'VariableNames',labels.Properties.VariableNames);

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


%% Compute component event probabilities

% Event probabilities
con = nan(N.assemblies, N.conditions, sum(T.condition));
sub = nan(N.assemblies, max(N.subjects), N.conditions, sum(T.condition));
for c = 1:N.conditions
	for t = 1:T.condition(c)
		con(:,c,t) = diag(activities.cond.events{c}(:,t) * activities.cond.events{c}(:,t)');
	end
	for s = 1:N.subjects(c)
		for t = 1:T.scan
			sub(:,s,c,t) = diag(activities.subj.events{s,c}(:,t) * activities.subj.events{s,c}(:,t)');
		end
	end
end
clear c s t
events.cond.prob.av = mean(con, 3, 'omitnan');
events.cond.prob.sd = std(con, 0, 3, 'omitnan');
events.subj.prob.av = mean(sub, 4, 'omitnan');
events.subj.prob.sd = std(sub, 0, 4, 'omitnan');
clear con sub

% Visualize means, standard deviations of event probabilities
metricFig(N.fig) = figure;
bar(events.cond.prob.av); hold on;
errorbar((1:N.assemblies)-0.15, events.cond.prob.av(:,1), events.cond.prob.sd(:,1), '.b');
errorbar((1:N.assemblies)+0.15, events.cond.prob.av(:,2), events.cond.prob.sd(:,2), '.r');
legend('Controls','Patients');
title('Assembly Event Probability');
xlabel('Assembly');
ylabel('Event Probability');
N.fig = N.fig + 1;

% Visualize event probability distributions
metricFig(N.fig) = figure;
for c = 1:N.conditions
	subplot(N.conditions, 1, c);
	bar(squeeze(events.subj.prob.av(:,:,c))); hold on;
	title([labels.Properties.VariableNames{c}(1:end-1), ' Assembly Event Probabilities']);
	xlabel('Assembly');
	ylabel('Event Probability');
end
N.fig = N.fig + 1;


%% Compute component activity probabilities

con = nan(N.assemblies, N.conditions, sum(T.condition));
sub = nan(N.assemblies, max(N.subjects), N.conditions, sum(T.condition));
for c = 1:N.conditions
	for t = 1:T.condition(c)
		con(:,c,t) = diag(activities.cond.activations{c}(:,t) * activities.cond.activations{c}(:,t)');
	end
	for s = 1:N.subjects(c)
		for t = 1:T.scan
			sub(:,s,c,t) = diag(activities.subj.activations{s,c}(:,t) * activities.subj.activations{s,c}(:,t)');
		end
	end
end
clear c s t
activities.cond.prob.av = mean(con, 3, 'omitnan');
activities.cond.prob.sd = std(con, 0, 3, 'omitnan');
activities.subj.prob.av = mean(sub, 4, 'omitnan');
activities.subj.prob.sd = std(sub, 0, 4, 'omitnan');

% Visualize means, standard deviations of event probabilities
metricFig(N.fig) = figure;
bar(activities.cond.prob.av); hold on;
errorbar((1:N.assemblies)-0.15, activities.cond.prob.av(:,1), activities.cond.prob.sd(:,1), '.b');
errorbar((1:N.assemblies)+0.15, activities.cond.prob.av(:,2), activities.cond.prob.sd(:,2), '.r');
legend('Controls','Patients');
title('Assembly Activation Probability');
xlabel('Assembly');
ylabel('Activation Probability');
N.fig = N.fig + 1;

% Visualize event probability distributions
metricFig(N.fig) = figure;
for c = 1:N.conditions
	subplot(N.conditions, 1, c);
	bar(squeeze(activities.subj.prob.av(:,:,c))); hold on;
	title([labels.Properties.VariableNames{c}(1:end-1), ' Assembly Activation Probabilities']);
	xlabel('Assembly');
	ylabel('Activation Probability');
end
N.fig = N.fig + 1;


%% Compute component lifetime histograms & means

% Find lifetimes for each assembly
for c = 1:N.conditions
	activities.cond.lifetime.hist{c} = lifetime(activities.cond.activations{c}, activities.cond.events{c});
	activities.cond.lifetime.av(:,c) = sum(squeeze(con(:,c,:)), 2, 'omitnan')./sum(activities.cond.events{c}, 2, 'omitnan');
	activities.cond.lifetime.sd(:,c) = std(activities.cond.lifetime.hist{c}, 0, 2, 'omitnan');
	for s = 1:N.subjects(c)
		activities.subj.lifetime.hist{s,c} = lifetime(activities.subj.activations{s,c}, activities.subj.events{s,c});
		activities.subj.lifetime.av(:,s,c) = sum(squeeze(sub(:,s,c,:)), 2, 'omitnan')./sum(activities.subj.events{s,c}, 2, 'omitnan');
		activities.subj.lifetime.sd(:,s,c) = std(activities.subj.lifetime.hist{s,c}, 0, 2, 'omitnan');
	end
end
clear c s con sub

% Visualize means, standard deviations of event lifetimes
metricFig(N.fig) = figure;
bar(activities.cond.lifetime.av); hold on;
errorbar((1:N.assemblies)-0.15, activities.cond.lifetime.av(:,1), activities.cond.lifetime.sd(:,1), '.b');
errorbar((1:N.assemblies)+0.15, activities.cond.lifetime.av(:,2), activities.cond.lifetime.sd(:,2), '.r');
legend(labels.Properties.VariableNames);
title('Mean Assembly Lifetimes');
xlabel('Assembly');
ylabel('Lifetime in TR');
N.fig = N.fig + 1;

% Visualize event lifetime distributions
metricFig(N.fig) = figure;
for c = 1:N.conditions
	subplot(N.conditions, 1, c);
	bar(squeeze(activities.subj.lifetime.av(:,:,c))); hold on;
	title([labels.Properties.VariableNames{c}(1:end-1), ' Assembly Lifetimes']);
	xlabel('Assembly');
	ylabel('Lifetime in TR');
end
N.fig = N.fig + 1;
clear c


%% Compute component switching probabilities




%% Compute component transition matrices




%% Compute component-wise Kuramoto order parameter & metastability

% Preallocate storage arrays for metastability
metastable.LEICA.cond = nan(1, N.conditions);
metastable.LEICA.subj = nan(max(N.subjects), N.conditions);

% Preallocate storage arrays for Kuramoto order parameter
kuramoto.LEICA.cond = nan(N.conditions, max(T.condition));
kuramoto.LEICA.subj = nan(max(N.subjects), N.conditions, T.scan);

% Extract metastability from BOLD data
[kuramoto.LEICA.concat, metastable.LEICA.concat] = findStability(activities.concat.TS);
for c = 1:N.conditions
	[kuramoto.LEICA.cond(c, 1:T.condition(c)), metastable.LEICA.cond(c)] = findStability(activities.cond.TS{c});
	for s = 1:N.subjects
		[kuramoto.LEICA.subj(s,c,:), metastable.LEICA.subj(s,c)] = findStability(activities.subj.TS{s,c});
	end
end
clear c s

% Convert metastabilities to tables
metastable.LEICA.cond = array2table(metastable.LEICA.cond, 'VariableNames', {'Controls','OCD'});
metastable.LEICA.subj = array2table(metastable.LEICA.subj, 'VariableNames', {'Controls','OCD'});

% Visualize Kuramoto order parameter
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

% Visualize metastability
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

% Preallocate storage
entro.cond = nan(1, N.conditions);
entro.subj = nan(max(N.subjects), N.conditions);

% Compute entropy for each condition and each subject.
% QUESTION: can compute entropy for each individual assembly?
entro.concat = HShannon_kNN_k_estimation(activities.concat.TS, co);
for c = 1:N.conditions
	entro.cond(c) = HShannon_kNN_k_estimation(activities.cond.TS{c}, co);
	for s = 1:N.subjects(c)
		entro.subj(s, c) = HShannon_kNN_k_estimation(activities.subj.TS{s,c}, co);
	end
end
clear n c s

% Convert entropies to tables
entro.cond = array2table(entro.cond, 'VariableNames', labels.Properties.VariableNames);
entro.subj = array2table(entro.subj, 'VariableNames', labels.Properties.VariableNames);

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

% Save figures
savefig(metricFig, fullfile(path{8}, fileName), 'compact')
clear metricFig

% Save variables
save(fullfile(path{8}, fileName));


