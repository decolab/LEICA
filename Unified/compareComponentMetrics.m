function [h, p, tstat, FDR, Sidak] = compareComponentMetrics(cind, dim, ttypes, pTarget, met, N, groups, comps)
%  Summary of this function goes here
%	

%% Compare metric across components

% % Preallocate storage arrays
% k = cell(numel(ttypes), N.comp);
% r = cell(numel(ttypes), N.comp);

% Allocate threshold storage
FDR = nan(N.IC, numel(ttypes), N.comp);
Sidak = nan(N.IC, numel(ttypes), N.comp);
p = nan(N.IC, numel(ttypes), N.comp);
tstat = nan(N.IC, numel(ttypes), N.comp);
h = nan(N.IC, numel(ttypes), N.comp);

% Test for differences
for t = 1:numel(ttypes)
    disp(strcat("Running ", ttypes(t), " tests on entropy."));
    for c = 1:N.comp
        % Run comparisons
        disp(strjoin(['Comparing', groups(comps(c,1)), 'and', groups(comps(c,2))], ' '));
        [h{t,c}, p(:,t,c), tstat(:,t,c)] = robustTests(squeeze(met(:,1:N.subjects(comps(c,1)),comps(c,1))), squeeze(met(:,1:N.subjects(comps(c,2)),comps(c,2))), [], 'p',pTarget, 'testtype',ttypes(t));
    end

    % Multiple comparison correction
    if N.comp > 1
        p_dum = reshape(squeeze(p(:,t,:)), [N.comp*N.IC, 1]);	% reshape p-value array
        [f, ~, S] = mCompCorr([], p_dum, pTarget);				% run mutliple comparison tests
        FDR(:,t,:) = reshape(f, [N.IC, N.comp]);
        Sidak(:,t,:) = reshape(S, [N.IC, N.comp]);
    end
end
clear c f S p_dum







%%






% Display test p-values
figure(F(nFig-4));
sgtitle([groups{comps(c,1)}, ' vs. ', groups{comps(c,2)}]);	% sgtitle('p-values');
ax(2) = subplot(numel(ttypes), 1, n(2)); 
n(2) = n(2)+1;
barh(ax(2), p(:,t,c)); hold on;
plot((max(p(:,t,c),[],'all')*1.2).*ones(nnz(h{t,c}{:}),1), find(h{t,c}{:}), '*r');
title(strjoin([ttypes(t), "test"]));    % , labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);
ylabel('Component'); xlabel('p-value');

% Display test effect sizes
figure(F(nFig-3));
sgtitle([groups{comps(c,1)}, ' vs. ', groups{comps(c,2)}]);	% sgtitle('Effect Sizes');
ax(1) = subplot(numel(ttypes), 1, n(1));
n(1) = n(1)+1;
barh(ax(1), tstat(:,t,c)); hold on;
plot(sign(tstat(logical(h{t,c}{:}),t,c)).*(max(abs(tstat(:,t,c)),[],'all')*1.2).*ones(nnz(h{t,c}{:}),1), find(h{t,c}{:}), '*r');
title(strjoin([ttypes(t), "test"]));    %, labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);
ylabel('Component'); xlabel('Effect Size');




end