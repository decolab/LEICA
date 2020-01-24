% Compute entropy for each condition and each subject.
% QUESTION: can compute entropy for each individual assembly?
for c = 1:N.conditions
	entro.cond.all(c) = HShannon_kNN_k_estimation(activities.cond.TS{c}, co);
	for ass = 1:13
		entro.cond.ass(ass, c) = HShannon_kNN_k_estimation(activities.cond.TS{c}(ass,:), co);
		for s = 1:N.subjects(c)
			entro.subj.ass(ass, s, c) = HShannon_kNN_k_estimation(activities.subj.TS{s,c}(ass,:), co);
		end
	end
	
	
	% 
	disp([num2str(entro.cond.all(c)), ' ', num2str(sum(entro.cond.ass(:,c)))]);
	subplot(1,2,c);
	imagesc(entro.subj.ass(:,1:N.subjects(c),c));
	colorbar;
end
clear c s ass

% test assembly-wise entropy
p = nan(13,1);
h = nan(13,1);
for ass = 1:13
	[h(ass),p(ass)] = kstest2(entro.subj.ass(ass,1:N.subjects(1),1), entro.subj.ass(ass,1:N.subjects(2),2));
end
s = sort(FDR_benjHoch(p, pval.target));