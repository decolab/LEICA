%% 4) Visualize assemblies

l = labels.Properties.VariableNames;

% Visualize membership weights
F(N.fig) = figure;
N.fig = N.fig+1;
d = abs(N.IC - max(N.IC));
for c = 1:N.conditions
	for k = 1:N.IC(c)
		subplot(N.conditions, max(N.IC), k+(c-1)*max(N.IC)); hold on;

		I = memberships{c}(:,k) > 0;
		m = memberships{c}(:,k); m(I) = 0;
		barh(1:numel(I), m, 'r');

		I = memberships{c}(:,k) <= 0;
		m = memberships{c}(:,k); m(I) = 0;
		barh(1:numel(I), m, 'g');

		title([l{c}, ': IC ', num2str(k)]);
		ylim([0 numel(I)+1]);
		xlim([-max(abs(memberships{c}(:,k))) max(abs(memberships{c}(:,k)))]);
	end
end
clear k m n I


% Visualize component relationships
F(N.fig) = figure; colormap jet
N.fig = N.fig+1;
for c = 1:N.conditions
	i = (c-1)*N.conditions;
	v = i + [1,2];
	
	subplot(2, N.conditions*2, v);
	imagesc(activities{c}); colorbar;
	title([l{c}, ': Component Activation']);
	xlabel('Time Points');
	ylabel('Components');

	% Component spatial correlation
	subplot(2, N.conditions*2, 5+i);
	comp.spatial = corr(memberships{c});
	imagesc(comp.spatial); colorbar;
	title([l{c}, ': Component Spatial Correlation']);

	% Component temporal overlap
	subplot(2, N.conditions*2, 6+i);
	imagesc(corr(activities{c}')); colorbar;
	title([l{c}, ': IC Temporal Correlations']);
end
clear c v i


% Visualize inter-group compontent distances
F(N.fig) = figure; colormap jet
N.fig = N.fig+1;
cTypes = {'correlation', 'spearman', 'cosine'};
cTitles = {'Pearson Distance', 'Spearman Distance', 'Cosine Distance'};
title([l{1}, ' vs. ', l{2}]);
for k = 1:numel(cTypes)
	subplot(2, numel(cTypes), k);
	colormap jet
	imagesc(pdist2(memberships{1}', memberships{2}', cTypes{k})); colorbar;
	title(['Spatial ', cTitles{k}]);
	xlabel(l{2});
	ylabel(l{1});
	
	subplot(2, numel(cTypes), k+numel(cTypes));
	colormap jet
	imagesc(pdist2(activities{1}(:,1:T.scan*min(N.subjects)), activities{2}(:,1:T.scan*min(N.subjects)), cTypes{k})); colorbar;
	title(['Temporal ', cTitles{k}]);
	xlabel(l{2});
	ylabel(l{1});
end


% Compare LEICA assemblies to RSNs

% Load Yeo RSNs
load(fullfile(path{4}, 'Atlases/AAL', 'RSNsAAL.mat'));
RSNs = LR_version_symm(Yeo_AAL');
clear Yeo_AAL

% Visualize relations between ICs and RSNs: Pearson, Spearman, cosine
F(N.fig) = figure; hold on;
N.fig = N.fig+1;
for c = 1:N.conditions
	for k = 1:numel(cTypes)
		subplot(N.conditions, numel(cTypes), k+(c-1)*numel(cTypes));
		colormap jet
		imagesc(pdist2(memberships{c}', RSNs', cTypes{k})); colorbar;
		title(['ICs vs. RSNs: ', cTitles{k}]);
		xlabel('Yeo RSNs');
		ylabel([l{c}, ': ICs']);
	end
end
clear cTypes cTitles k c


% Visualize projection relationships
F(N.fig) = figure; colormap jet
N.fig = N.fig+1;
d = N.conditions:-1:1;
for c = 1:N.conditions
	i = (c-1)*N.conditions;
	v = i + [1,2];
	
	subplot(3, N.conditions*2, v);
	imagesc(projections{c}); colorbar;
	title(['Projection: ', l{d(c)}, ' onto ', l{c}, ' Components']);
	xlabel('Time Points');
	ylabel('Components');
	
	% Component temporal overlap
	subplot(3, N.conditions*2, [5 6 9 10]+i);
	imagesc(corr(activities{c}')); colorbar;
	title([l{c}, ': Projections Temporal Correlations']);
end
clear d c v i l