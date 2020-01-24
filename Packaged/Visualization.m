%% ICA Extraction of Network States
%	This script computes the neural assemblies and assembly time courses
% from the data processed by the extraction script and computes the Shannon
% entropies of each assembly's time course.  In addition, it computes the
% total Shannon entropy of each condition.


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

% Plot p-values of assembly activations between conditions
subplot(3,1,1);
xlim([1, N.ROI]);
bar(1:N.ROI, pval.BOLD);
hold on;
% Mark significantly different assemblies
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