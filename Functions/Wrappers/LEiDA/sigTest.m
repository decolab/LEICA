function [sig, pval] = sigTest(distribution, sigvalue, nAssembly)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Confirm that distribution is in proper format
if size(distribution{1},1) ~= nAssembly
    distribution{1} = distribution{1}';
    distribution{2} = distribution{2}';
end

% Preallocate assembly significance vectors
pval = nan(nAssembly,1);
sig = logical((zeros(nAssembly, 3)));
sig = array2table(sig, 'VariableNames',{'FDR' 'Bonferroni' 'Sidak'});

% Run multiple pairwise Kolmogorov-Smirnov two-tailed tests
for ass = 1:nAssembly
	a = distribution{1}(ass,:);
	b = distribution{2}(ass,:);
	[~, pval(ass)] = kstest2(a,b);
end

% Run false discovery rate test
[FDR] = FDR_benjHoch(pval, sigvalue);
FDR = sort(FDR);
sig{:,'FDR'}(FDR) = true;

% Run Bonferroni multiple comparison correction
alpha = sigvalue/nAssembly;
sig{:,'Bonferroni'} = (pval < alpha);

% Run Dunn-Sidak multiple comparison correction
alpha = 1-(1-sigvalue)^(1/nAssembly);
sig{:,'Sidak'} = (pval < alpha);

end


