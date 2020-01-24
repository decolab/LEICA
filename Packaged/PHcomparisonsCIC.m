%% Comparison of Metric Distributions between States
%	


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% COMPARISON PIPELINE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%	SETUP

%% Set paths & filenames

clear; close all; clc
set(groot,'defaultLegendAutoUpdate','off');	% prevents legend from automatically adding new data

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
path{7,1} = fullfile(path{2},'Functions','LEICA');
path{8,1} = fullfile(path{2},'Results', 'LEICA');

% Add relevant paths
addpath(path{6});
addpath(genpath(path{7}));


%% Set file names & load data

% Define files to load
loadFile = 'LEICA90_CIC_Metrics';

% Load data
metricFig = openfig(fullfile(path{8}, loadFile));
load(fullfile(path{8}, loadFile));

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_', S{2}, '_Comparisons');
clear loadFile S


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


%% Compute LEICA assembly matrices

% Preallocate sotrage array
ICs = nan(N.ROI, N.ROI, size(memberships,2));

% Compute IC matrices
for i = 1:size(memberships,2)
	ICs(:,:,i) = memberships(:,i) * memberships(:,i)';
end
clear i


%% Compare LEICA assemblies to RSNs

% Load Yeo RSNs
load(fullfile(path{4}, 'Atlases/AAL', 'RSNsAAL.mat'));
RSNs = Yeo_AAL';
clear Yeo_AAL

% Visualize relations between ICs and RSNs
N.fig = N.fig+1;
metricFig(N.fig) = figure; hold on;

% Visualize Pearson distance: ICs to RSNs
subplot(1,3,1);
colormap jet
imagesc(pdist2(memberships', RSNs', 'correlation')); colorbar;
title('Pearson Distance: ICs to RSNs');
xlabel('Yeo RSNs');
ylabel('LEICA ICs');

% Visualize Spearman distance: ICs to RSNs
subplot(1,3,2);
colormap jet
imagesc(pdist2(memberships', RSNs', 'spearman')); colorbar;
title('Spearman Distance: ICs to RSNs');
xlabel('Yeo RSNs');
ylabel('LEICA ICs');

% Visualize cosine distance: ICs to RSNs
subplot(1,3,3);
colormap jet
imagesc(pdist2(memberships', RSNs', 'cosine')); colorbar;
title('Cosine Distance: ICs to RSNs');
xlabel('Yeo RSNs');
ylabel('LEICA ICs');


%% Compare transition matrices between states

% Preallocate storage arrays


% Compute distance(s) between matrices



%% Compare component activation time series

% Preallocate storage arrays
sig.TS.h = nan(N.assemblies, 1);
sig.TS.p = nan(N.assemblies, 1);

% Run tests for each assembly
for ass = 1:N.assemblies
	[sig.TS.h(ass), sig.TS.p(ass)] = kstest2(activities.cond.TS{1}(ass,:), activities.cond.TS{2}(ass,:), 'Alpha',pval.target);
end
clear con pat ass

% Run FDR correction
[sig.TS.FDR] = sort(FDR_benjHoch(sig.TS.p, pval.target));

% Run Bonferroni multiple comparison correction
sig.TS.Bonferroni = (sig.TS.p < (pval.target/N.assemblies));

% Run Dunn-Sidak multiple comparison correction
alpha = 1-(1-pval.target)^(1/N.assemblies);
sig.TS.Sidak = (sig.TS.p < alpha);
clear alpha

% Visualize results
nFig = 2;
if ~isempty(sig.TS.FDR)
	set(metricFig(nFig),'defaultLegendAutoUpdate','off');
	cax = metricFig(nFig).CurrentAxes;
	mark = 1.3*ones(numel(sig.TS.FDR),1);
	plot(cax, sig.TS.FDR, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax

% Visualize ICs with significantly different time series
N.fig = N.fig+1;
metricFig(N.fig) = figure; hold on;
colormap jet
for i = 1:length(sig.TS.FDR)
	subplot(3, ceil(length(sig.TS.FDR)/3), i);
	imagesc(ICs(:,:,sig.TS.FDR(i))); colorbar;
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

% Preallocate storage arrays
sig.active.h = nan(N.assemblies, 1);
sig.active.p = nan(N.assemblies, 1);

% Run tests for each assembly
for ass = 1:N.assemblies
	con = activities.subj.TSav(ass,:,1); con = con(isfinite(con));
	pat = activities.subj.TSav(ass,:,2); pat = pat(isfinite(pat));
	[sig.active.h(ass), sig.active.p(ass)] = kstest2(con, pat, 'Alpha',pval.target);
end
clear con pat ass

% Run FDR correction
[sig.active.FDR] = FDR_benjHoch(sig.active.p, pval.target);

% Run Bonferroni multiple comparison correction
sig.active.Bonferroni = (sig.active.p < (pval.target/N.assemblies));

% Run Dunn-Sidak multiple comparison correction
alpha = 1-(1-pval.target)^(1/N.assemblies);
sig.active.Sidak = (sig.active.p < alpha);
clear alpha

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


%% Compare component event probabilities

% Preallocate storage arrays
sig.events.h = nan(N.assemblies, 1);
sig.events.p = nan(N.assemblies, 1);

% Run tests for each assembly
for ass = 1:N.assemblies
	con = events.subj.prob.av(ass,:,1); con = con(isfinite(con));
	pat = events.subj.prob.av(ass,:,2); pat = pat(isfinite(pat));
	[sig.events.h(ass), sig.events.p(ass)] = kstest2(con, pat, 'Alpha',pval.target);
end
clear con pat ass

% Run FDR correction
[sig.events.FDR] = FDR_benjHoch(sig.events.p, pval.target);

% Run Bonferroni multiple comparison correction
sig.events.Bonferroni = (sig.events.p < (pval.target/N.assemblies));

% Run Dunn-Sidak multiple comparison correction
alpha = 1-(1-pval.target)^(1/N.assemblies);
sig.events.Sidak = (sig.events.p < alpha);
clear alpha

% Visualize results
nFig = nFig+2;
ind = find(sig.events.FDR);
mark = ones(numel(ind),1);
if ~isempty(ind)
	cax = metricFig(nFig).CurrentAxes;
	scatter(cax, ind, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax


%% Compare component activity probabilities

% Preallocate storage arrays
sig.prob.h = nan(N.assemblies, 1);
sig.prob.p = nan(N.assemblies, 1);

% Run tests for each assembly
for ass = 1:N.assemblies
	con = activities.subj.prob.av(ass,:,1); con = con(isfinite(con));
	pat = activities.subj.prob.av(ass,:,2); pat = pat(isfinite(pat));
	[sig.prob.h(ass), sig.prob.p(ass)] = kstest2(con, pat, 'Alpha',pval.target);
end
clear con pat ass

% Run FDR correction
sig.prob.FDR = FDR_benjHoch(sig.prob.p, pval.target);

% Run Bonferroni multiple comparison correction
sig.prob.Bonferroni = (sig.prob.p < (pval.target/N.assemblies));

% Run Dunn-Sidak multiple comparison correction
alpha = 1-(1-pval.target)^(1/N.assemblies);
sig.prob.Sidak = (sig.prob.p < alpha);
clear alpha

% Visualize results
nFig = nFig+2;
ind = find(sig.prob.FDR);
mark = ones(numel(ind),1);
if ~isempty(ind)
	cax = metricFig(nFig).CurrentAxes;
	scatter(cax, ind, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax


%% Compare mean dwell times per assembly

% Preallocate storage arrays
sig.lifetime.h = nan(N.assemblies, 1);
sig.lifetime.p = nan(N.assemblies, 1);

% Run tests for each assembly
for ass = 1:N.assemblies
	con = cell(N.conditions, 1);
	for c = 1:N.conditions
		for s = 1:N.subjects(c)
			con{c} = horzcat(con{c}, activities.subj.lifetime.hist{s,c}(ass,:));
		end
		con{c} = con{c}(isfinite(con{c}));
	end
	% con = activities.subj.lifetime.av(ass,:,1); con = con(isfinite(con));
	% pat = activities.subj.lifetime.av(ass,:,2); pat = pat(isfinite(pat));
	[sig.lifetime.h(ass), sig.lifetime.p(ass)] = kstest2(con{1}, con{2}, 'Alpha',pval.target);
end
clear con pat ass

% Run FDR correction
sig.lifetime.FDR = FDR_benjHoch(sig.lifetime.p, pval.target);

% Run Bonferroni multiple comparison correction
sig.lifetime.Bonferroni = (sig.lifetime.p < (pval.target/N.assemblies));

% Run Dunn-Sidak multiple comparison correction
alpha = 1-(1-pval.target)^(1/N.assemblies);
sig.lifetime.Sidak = (sig.lifetime.p < alpha);
clear alpha

% Visualize results
nFig = nFig+2;
ind = find(sig.lifetime.FDR);
mark = ones(numel(ind),1);
if ~isempty(ind)
	cax = metricFig(nFig).CurrentAxes;
	scatter(cax, ind, mark, '*', 'MarkerEdgeColor','k');
end
clear ind mark cax


%% Compare synchrony (Kuramoto order parameter) between states

% Run tests for each assembly
con = mean(squeeze(kuramoto.LEICA.subj(:,1,:)), 2, 'omitnan');
pat = mean(squeeze(kuramoto.LEICA.subj(:,2,:)), 2, 'omitnan');
[sig.sync.h, sig.sync.p] = kstest2(con, pat, 'Alpha',pval.target);
clear con pat

% Run FDR correction
% sig.sync.FDR = FDR_benjHoch(sig.sync.p, pval.target);

% Run Bonferroni multiple comparison correction
% sig.sync.Bonferroni = (sig.sync.p < (pval.target/N.assemblies));

% Run Dunn-Sidak multiple comparison correction
% alpha = 1-(1-pval.target)^(1/N.assemblies);
% sig.sync.Sidak = (sig.sync.p < alpha);
% clear alpha

% Visualize results
nFig = nFig+2;
ind = find(sig.sync.h);
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

% Run FDR correction
% sig.meta.FDR = FDR_benjHoch(sig.meta.p, pval.target);

% Run Bonferroni multiple comparison correction
% sig.meta.Bonferroni = (sig.meta.p < (pval.target/N.assemblies));

% Run Dunn-Sidak multiple comparison correction
% alpha = 1-(1-pval.target)^(1/N.assemblies);
% sig.meta.Sidak = (sig.meta.p < alpha);
% clear alpha

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

% Run tests for each assembly
con = entro.LEICA.subj{:,'Controls'}(isfinite(entro.LEICA.subj{:,'Controls'}));
pat = entro.LEICA.subj{:,'Patients'}(isfinite(entro.LEICA.subj{:,'Patients'}));
if adtest(con) && adtest(pat)
	[sig.entro.h, sig.entro.p] = kstest2(con, pat, 'Alpha',pval.target);
else
	[sig.entro.h, sig.entro.p] = ttest2(con, pat, 'Alpha',pval.target);
end
clear con pat

% Run FDR correction
% sig.entro.FDR = FDR_benjHoch(sig.entro.p, pval.target);

% Run Bonferroni multiple comparison correction
% sig.entro.Bonferroni = (sig.entro.p < (pval.target/N.assemblies));

% Run Dunn-Sidak multiple comparison correction
% alpha = 1-(1-pval.target)^(1/N.assemblies);
% sig.entro.Sidak = (sig.entro.p < alpha);
% clear alpha

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


%% Save results

% Save figures
savefig(metricFig, fullfile(path{8}, fileName), 'compact')
clear metricFig

% Save variables
save(fullfile(path{8}, fileName));







%% 


%% 1) Find assemblies which significantly differ between conditions

% % Locate assemblies with significantly different activation distributions
% [sig.activation, pval.activation] = sigTest(activations.cond, pval.target, N.assemblies(1));
% N.assemblies(1,2) = nnz(sig.activation{:,'Bonferroni'});		% List number of significantly different assemblies
% 
% % Locate assemblies with significantly different entropy distributions
% [sig.entro.ass, pval.entro.ass] = sigTest(entro.cond, pval.target, N.assemblies(1));
% N.assemblies(1,5) = nnz(sig.entro.ass{:,'Bonferroni'});			% List number of significantly different assemblies
% 
% % Locate assemblies with significantly different hierarchy distributions
% [sig.cohesive, pval.cohesive] = sigTest(cohesiveness, pval.target, N.ROI);
% N.assemblies(1,6) = nnz(sig.cohesive{:,'Bonferroni'});			% List number of significantly different assemblies
% 
% % Compare assembly activation probabilities between conditions
% [sig.prob, pval.prob] = sigTest(activations.prob, pval.target, N.assemblies(1));
% N.assemblies(1,3) = nnz(sig.prob{:,'Bonferroni'});
% 
% % Compare assembly metastabilities between conditions
% [sig.mstable, pval.mstable] = kstest2(metastable(:,1), metastable(:,2));
% N.assemblies(1,4) = nnz(sig.mstable);
% 
% % Convert number of assemblies to table
% N.assemblies = array2table(N.assemblies, 'VariableNames',clist);
% 
% % Save interim results
% save(fullfile(path{8},fileName));
% 
% 
% %% 2) Find total subjectwise entropy
% 
% % Preallocate array
% entro.subj = nan(max(N.subjects), N.condition);
% 
% % Find total entropy (upper bound) per subject
% for c = 1:N.condition
% 	entro.subj(1:N.subjects(c),c) = sum(entro.cond{c},2);
% end
% entro.subj = array2table(entro.subj, 'VariableNames',I.Properties.VariableNames);
% clear c
% 
% % Compare using KS test
% [sig.entro.subj, pval.entro.subj] = kstest2(entro.subj{:,'Controls'}, entro.subj{:,'OCD'});
% 
% % Save interim results
% save(fullfile(path{8},fileName));
% 
% 
% %% Plot assemblywise comparison results
% 
% % Plot results of activation time series
% F(N.fig) = PlotComparison(activations.cond, pval.activation, pval.target, sig.activation, N, I.Properties.VariableNames, 'TimeSeries', 'Bonferroni');
% 
% % Plot results of activation probability
% N.fig = N.fig+1;
% F(N.fig) = PlotComparison(activations.prob, pval.prob, pval.target, sig.prob, N, I.Properties.VariableNames, 'Probability', 'Bonferroni');
% 
% % Plot results of entropy comparison
% N.fig = N.fig+1;
% F(N.fig) = PlotComparison(entro.cond, pval.entro.ass, pval.target, sig.entro.ass, N, I.Properties.VariableNames, 'Entropy', 'Bonferroni');
% 
% % Plot results of entropy comparison
% N.fig = N.fig+1;
% F(N.fig) = PlotComparison(cohesiveness, pval.cohesive, pval.target, sig.cohesive, N, I.Properties.VariableNames, 'Cohesion', 'Bonferroni');
% 
% % Save figure
% savefig(F, fullfile(path{8}, fileName), 'compact');
% 
% 
% %% Plot results of metastability comparison
% 
% % Open new figure
% N.fig = N.fig+1;
% F(N.fig) = figure(N.fig);
% 
% % Plot subjectwise entropies as box plot
% boxplot(metastable, entro.subj.Properties.VariableNames, 'Notch','on');
% title('Metastability by Condition');
% ylabel('Metastabiltiy by Subject');
% if sig.mstable == true
% 	hold on;
% 	yt = get(gca, 'YTick');
% 	axis([xlim    0  ceil(max(yt)*1.2)])
% 	xt = get(gca, 'XTick');
% 	hold on
% 	plot(xt, [1 1]*max(yt)*1.05, '-k',  mean(xt), max(yt)*1.1, '*k');
% 	ylim([min(yt)*0.8 max(yt)*1.15])
% 	hold off
% end
% 
% % Save figure(s)
% savefig(F, fullfile(path{8},fileName), 'compact');
% 
% 
% %% 4) Plot subjectwise Shannon entropy
% 
% % Open new figure
% N.fig = N.fig+1;
% F(N.fig) = figure(N.fig);
% 
% % Plot subjectwise entropies as box plot
% boxplot(table2array(entro.subj), entro.subj.Properties.VariableNames, 'Notch','on');
% title('Total Subjectwise Entropy');
% ylabel('Total Entropy');
% if sig.entro.subj == true
% 	hold on;
% 	yt = get(gca, 'YTick');
% 	axis([xlim    0  ceil(max(yt)*1.2)])
% 	xt = get(gca, 'XTick');
% 	hold on
% 	plot(xt, [1 1]*max(yt)*1.05, '-k',  mean(xt), max(yt)*1.1, '*k');
% 	ylim([min(yt)*0.8 max(yt)*1.15])
% 	hold off
% end
% 
% % Save figure(s)
% savefig(F, fullfile(path{8},fileName), 'compact');
% 
% 
% %% Compute & plot transition matrices
% 
% 
% 
% 
% clear F

