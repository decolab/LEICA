function [F] = plotComponentComparison(h, p, tstat, groups, comps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Display test p-values
F(1) = figure;
sgtitle([groups{comps(c,1)}, ' vs. ', groups{comps(c,2)}]);	% sgtitle('p-values');
ax(2) = subplot(numel(ttypes), 1, n(2)); 
barh(ax(2), p(:,t,c)); hold on;
plot((max(p(:,t,c),[],'all')*1.2).*ones(nnz(h{t,c}{:}),1), find(h{t,c}{:}), '*r');
title(strjoin([ttypes(t), "test"]));    % , labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);
ylabel('Component'); xlabel('p-value');

% Display test effect sizes
F(2) = figure;
sgtitle([groups{comps(c,1)}, ' vs. ', groups{comps(c,2)}]);	% sgtitle('Effect Sizes');
ax(1) = subplot(numel(ttypes), 1, n(1));
barh(ax(1), tstat(:,t,c)); hold on;
plot(sign(tstat(logical(h{t,c}{:}),t,c)).*(max(abs(tstat(:,t,c)),[],'all')*1.2).*ones(nnz(h{t,c}{:}),1), find(h{t,c}{:}), '*r');
title(strjoin([ttypes(t), "test"]));    %, labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);
ylabel('Component'); xlabel('Effect Size');

end

