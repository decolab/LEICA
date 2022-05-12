function [h, p, tstat, FDR, Sidak] = compareGroupMetrics(ttype, pTarget, met, N, groups, comps)
%  Summary of this function goes here
%	

%% Compare metric across components

% Allocate threshold storage
FDR = nan(1, N.comp);
Sidak = nan(1, N.comp);
p = nan(1, N.comp);
tstat = nan(1, N.comp);
h = nan(1, N.comp);

% Test for differences
disp(strjoin(["Running", ttype, "test."], " "));
for c = 1:N.comp
    % Run comparisons
    cont(c) = strjoin([groups(comps(c,1)), "v.", groups(comps(c,2))], ' ');
    disp(strjoin(['Comparing', cont(c)], ' '));
    [h(c), p(c), tstat(c), FDR(c), Sidak(c)] = robustTests(squeeze(met(:,comps(c,1))), squeeze(met(:,comps(c,2))), [], 'p',pTarget, 'testtype',ttype);
end

% Convert to tables
FDR = table(FDR, 'VariableNames',cont);
Sidak = table(Sidak, 'VariableNames',cont);
p = table(p, 'VariableNames',cont);
tstat = table(tstat, 'VariableNames',cont);
h = table(h, 'VariableNames',cont);

% Multiple comparison correction
if N.comp > 1
    p_dum = reshape(squeeze(p(:,:)), [N.comp*N.IC, 1]);	% reshape p-value array
    [f, ~, S] = mCompCorr([], p_dum, pTarget);				% run mutliple comparison tests
    FDR(:,:) = reshape(f, [N.IC, N.comp]);
    Sidak(:,:) = reshape(S, [N.IC, N.comp]);
end
clear c f S p_dum

end