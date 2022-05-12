function [h, p, tstat, FDR, Sidak] = compareComponentMetrics(ttype, pTarget, met, N, groups, comps)
%  Summary of this function goes here
%	

%% Compare metric across components

% Allocate threshold storage
FDR = nan(N.IC, N.comp);
Sidak = nan(N.IC, N.comp);
p = nan(N.IC, N.comp);
tstat = nan(N.IC, N.comp);
h = nan(N.IC, N.comp);

% Test for differences
disp(strjoin(["Running", ttype, "test."], " "));
for c = 1:N.comp
    % Run comparisons
    cont(c) = strjoin([groups(comps(c,1)), "v.", groups(comps(c,2))], ' ');
    disp(strjoin(['Comparing', cont(c)], ' '));
    [h(:,c), p(:,c), tstat(:,c), FDR(:,c), Sidak(:,c)] = robustTests(squeeze(met(:,1:N.subjects(comps(c,1)),comps(c,1))), squeeze(met(:,1:N.subjects(comps(c,2)),comps(c,2))), [], 'p',pTarget, 'testtype',ttype);
end

% Multiple comparison correction
if N.comp > 1
    p_dum = reshape(squeeze(p(:,:)), [N.comp*N.IC, 1]);	% reshape p-value array
    [f, ~, S] = mCompCorr([], p_dum, pTarget);				% run mutliple comparison tests
    FDR(:,:) = reshape(f, [N.IC, N.comp]);
    Sidak(:,:) = reshape(S, [N.IC, N.comp]);
end

% Convert to tables
FDR = array2table(FDR, 'VariableNames',cont);
Sidak = array2table(Sidak, 'VariableNames',cont);
p = array2table(p, 'VariableNames',cont);
tstat = array2table(tstat, 'VariableNames',cont);
h = array2table(h, 'VariableNames',cont);

clear c f S p_dum

end