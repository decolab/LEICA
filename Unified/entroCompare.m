%% Import data

clear; close all; clc;

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2}, 'OCD', 'Results','LEICA');
path{5,1} = fullfile(path{1},'Project','Atlases','AAL');

% Add relevant paths
fpath{1,1} = fullfile(path{1},'MATLAB','BrainNetViewer');
fpath{2,1} = fullfile(path{1},'MATLAB','permutationTest');
fpath{3,1} = fullfile(path{1},'Project','Functions');
fpath{4,1} = fullfile(path{1},'Project','LEICA','Functions');
for k = 1:numel(fpath)
	addpath(genpath(fpath{k}));
end
clear fpath k
addpath(fullfile(path{1},'MATLAB','spm12'));

% Determine whether to z-scale membership weights
zthresh = 0.5:0.2:1.5;
zscale = false;

% Load network labels
label_AAL90 = load(fullfile(path{5}, 'AAL_labels'));
label_AAL90 = string(label_AAL90.label90);
label_AAL90 = strip(LR_version_symm(label_AAL90));

% Generate MNI coordinates
coords_AAL90 = load(fullfile(path{5}, 'aal_cog.txt'), 'aal_cog');
coords_AAL90 = LR_version_symm(coords_AAL90);
origin = [65 45.5 35];							% center origin
MNIscale = 5.5/10;								% scale MNI coordinates
sphereScale = 2.5;								% scale sphere size
coords_AAL90 = MNIscale*coords_AAL90;			% scale MNI coordinates
clear label90 MNIscale

% Set brain rendering parameters
cortex.file = fullfile(path{2},'Data','MNI152_T1_2mm_brain_mask.nii');	% file containing cortical atlas
cortex.color = [0.9 0.9 0.9];					% color for cortical rendering
cortex.transparency = 0.1;						% set to 1 for opaque cortex
cortex.val = 0.3;								% set isonormal line spacing
cortex.view = [-90 90];							% set camera angle
rdux = 0.7;						% proportion of surface faces to keep (<1)

% Determine which files to compare
band = 'WideBand';
k = 1;
iteration = 2;
fname = strcat('_',band, '_k',num2str(k), '_Iteration',num2str(iteration));
pfix = {'FCD90_CIC_COS', 'MICA90_CIC_EXP', 'LEICA90_CIC_COS'};
clear assemblies distance band k iteration

% Set decompositions, spaces, comparisions
nFig = 1;
titles = {'ICA', 'MICA', 'LEICA'};	% compression type
spaces = {'dFC' 'IC'};				% space in which to compare
dim = {'Subject', 'Component'};		% dimensions to compare
ttype = {'kstest2', 'permutation'};	% set test types to run
pTarget = 0.05;						% target p-value
prange = 0.025:0.025:0.15;			% Set range of p-values to test

% Load entropies
N = cell(numel(pfix), numel(spaces));
entro = cell(numel(pfix), numel(spaces));
FCD = cell(numel(pfix), numel(spaces));
memberships = cell(numel(pfix), 1);	% have 2 dFC / FCD spaces: ROI and IC.
for f = 1:numel(pfix)
	e = load(fullfile(path{4}, strcat(pfix{f}, fname)), 'entro','N','memberships','FCD','labels','T');
	labels = e.labels.Properties.VariableNames;
	T = e.T;
	N{f,1} = e.N; N{f,1}.comp = N{f,1}.ROI; N{f,1} = rmfield(N{f,1}, 'ROI');
	% N{f,2} = N{f,1};
	N{f,2} = e.N; N{f,2}.comp = N{f,2}.IC; N{f,2} = rmfield(N{f,2}, 'IC');
	% entro{f,1} = e.entro.BOLD;
	entro{f,1} = e.entro.dFC;
	entro{f,2} = e.entro.IC;
	memberships{f,1} = e.memberships;
	FCD{f,1} = e.FCD.dFC;
	FCD{f,2} = e.FCD.IC;
end
clear e f


%% Display FCDs

% Preallocate array to store KS distances
ksdist.FCD = nan(size(FCD));

% Set pre-ICA titles
ttls = {'Standard', 'Mean', 'Leading Eigenvector'};	% compression type

for s = 1:size(FCD,2)
	F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1;
	for t = 1:size(FCD,1)

		% Plot examples of patient, control FCDs
		for l = 1:numel(labels)
			subplot(numel(labels)+1, size(FCD,1), t+size(FCD,1)*(l-1));
			colormap jet;
			ind = find(~mod(round(T.TR*1:T.scan), 25));
			imagesc(FCD{t,s}.subj{1,l}); colorbar; caxis([-1 1]);
			xticks([]);
			yticks(ind); yticklabels(round(T.TR.*ind));
			title([ttls{t}, ', ', labels{l}, ' 1']);
			xlabel('Time (s)'); ylabel('Time (s)');
            pbaspect([1 1 1]);
		end
		
		% Measure KS distance between patient, control FCD histograms
		D = cell(size(FCD{t,s}.subj, 2));
		grp = D;
		for f = 1:size(FCD{t,s}.subj, 2)
			D{f} = reshape(cell2mat(FCD{t,s}.subj(:,f)), [numel(cell2mat(FCD{t,s}.subj(:,f))), 1]);
			grp{f} = string(repmat(labels{f}, [numel(cell2mat(FCD{t,s}.subj(:,f))), 1]));
		end
		[ksdist.h(t,s), ksdist.p(t,s), ksdist.FCD(t,s)] = kstest2(D{1}(:), D{2}(:));
		% d(t,s) = mean(D{1}(:)) - mean(D{2}(:));
		
		% Plot patient, control FCD histograms
		subplot(numel(labels)+1, size(FCD,1), t+size(FCD,1)*2); hold on;
		grp = vertcat(grp{:});
		D = cell2mat(D);
		boxplot(D,grp, 'Notch','on');
		
		% histogram(cell2mat(FCD{t,s}.subj(:,1))', 'Normalization','probability', 'FaceAlpha',0.2, 'FaceColor',[1 0 0]);
		% histogram(cell2mat(FCD{t,s}.subj(:,2))', 'Normalization','probability', 'FaceAlpha',0.2, 'FaceColor',[0 0 1]);
		% legend(labels, 'location','northwest');
		title([ttls{t}, ' FCD Values']);
		ylabel('Correlation');	% ylabel('Probability');
	end
	sgtitle(['FCD in ', spaces{s}, ' Space']);
end
clear s t l ind con pat ttls grp D d


%% Compare mean entropy across methods

% Preallocate arrays
mEntro = cell(size(entro,1), numel(spaces), ndims(entro{1,1})-1);
cond = cell(1, N{1}.conditions);

% Set comparison order
dType = {'Subject Mean' 'IC Mean'};

for t = 1:numel(titles)
	for s = 1:numel(spaces)
        
        % Determine all possible pairwise comparisons
        C = nchoosek(1:N{t,s}.conditions, 2);
        N{t,s}.comp = size(C,1);
        
		for d = 1:numel(dType)
			
			% Compute means
			mEntro{t,s,d} = squeeze(mean(entro{t,s}, d, 'omitnan'));
			mEntro{t,s,d} = array2table(mEntro{t,s, d}, 'VariableNames', labels);
			
			% Pairwise comparisons between conditions
            for c = 1:N{t,s}.comp
                cond{1} = mEntro{t,s,d}{:,labels{C(c,1)}}(isfinite(mEntro{t,s,d}{:,labels{C(c,1)}}));
                cond{2} = mEntro{t,s,d}{:,labels{C(c,2)}}(isfinite(mEntro{t,s,d}{:,labels{C(c,2)}}));
                
                [sig.av.h(t,s,d,c,1), sig.av.p(t,s,d,c,1), sig.av.effsize(t,s,d,c,1)] = kstest2(cond{1}, cond{2});
                [sig.av.p(t,s,d,c,2), ~, sig.av.effsize(t,s,d,c,2)] = permutationTest(cond{1}, cond{2}, 10000, 'sidedness','both');
                if sig.av.p(t,s,d,c,2) < 0.05
                    sig.av.h(t,s,d,c,2) = 1;
                else
                    sig.av.h(t,s,d,c,2) = 0;
                end
            end
		end
	end 
end
clear s t c d

% Split into comparisons
sig.mSubj.h = squeeze(sig.av.h(:,:,1,:,:));
sig.mSubj.p = squeeze(sig.av.p(:,:,1,:,:));
sig.mSubj.effsize = squeeze(sig.av.effsize(:,:,1,:,:));
sig.mIC.h = squeeze(sig.av.h(:,:,2,:,:));
sig.mIC.p = squeeze(sig.av.p(:,:,2,:,:));
sig.mIC.effsize = squeeze(sig.av.effsize(:,:,2,:,:));

% Display mean entropy histograms
F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; hold on;
col = 0;
for d = 1:numel(dType)			% dimension average to display
	for s = 1:numel(spaces)		% space to display (dFC or IC)
		col = col+1;
		for t = 1:numel(titles)
			
			% Get bin sizes
            binwidth = nan(N{t,s}.conditions,1);
			f = figure; hold on;
            for c = 1:N{t,s}.conditions
                hg = histogram(mEntro{t,s,d}{:,labels{c}}, 'Normalization','probability');
                binwidth(c) = hg.BinWidth;
            end
			sz = min(binwidth);
			close(f); clear f hg c binwidth
			
			% Plot mean entropy histograms
			nc = numel(spaces)*(ndims(entro{t,s})-1);
			row = nc*(t-1);
			subplot(size(entro,1), nc, col+row); hold on;
            for c = 1:N{t,s}.conditions
                h{c} = histogram(mEntro{t,s,d}{:,labels{c}}, 'Normalization','probability', 'FaceAlpha',0.2, 'FaceColor',[1 0 0], 'BinWidth',sz);
            end
			
			% Plot means, distances for significant differences
            for c = 1:N{t,s}.comp
                if sig.av.h(t,s,d,c,1) || sig.av.h(t,s,d,c,2)
                    [mc(1),i(1)] = max(h{C(c,1)}.BinCounts); mc(1) = mc(1)/sum(h{C(c,1)}.BinCounts);
                    [mc(2),i(2)] = max(h{C(c,2)}.BinCounts); mc(2) = mc(2)/sum(h{C(c,2)}.BinCounts);
                    mp{1} = mean([h{C(c,1)}.BinEdges(i(1):i(1)+1), h{C(c,2)}.BinEdges(i(2):i(2)+1)]);
                    mp{2} = [mean(h{C(c,1)}.BinEdges(i(1):i(1)+1)), mean(h{C(c,2)}.BinEdges(i(2):i(2)+1))];
                    if sig.av.h(t,s,d,c,1) && sig.av.h(t,s,d,c,2)
                        plot(mp{1}, 1.05*max(mc), '*g');
                        plot(mp{2}, 1.02*[max(mc),max(mc)], '-g');
                        plot(mp{2}(1), 1.02*[max(mc),max(mc)], '+g');
                        plot(mp{2}(2), 1.02*[max(mc),max(mc)], '+g');
                    elseif sig.av.h(t,s,d,c,2)
                        plot(mp{1}, 1.05*max(mc), '*b');
                        plot(mp{2}, 1.02*[max(mc),max(mc)], '-b');
                        plot(mp{2}(1), 1.02*[max(mc),max(mc)], '+b');
                        plot(mp{2}(2), 1.02*[max(mc),max(mc)], '+b');
                    else
                        plot(mp{1}, 1.05*max(mc), '*r');
                        plot(mp{2}, 1.02*[max(mc),max(mc)], '-r');
                        plot(mp{2}(1), 1.02*[max(mc),max(mc)], '+r');
                        plot(mp{2}(2), 1.02*[max(mc),max(mc)], '+r');
                    end
                end
            end
            
			title([titles{t}, ': ', dType{d}, ', ', spaces{s}, ' Space']);
			legend(labels, 'location','northwest');
			xlabel('Mean Entropy'); ylabel('Probability');
		end
	end
end
clear d s t nc col row mc mp i h hg sz

% Visualize overall entropy distribution between groups
F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; hold on;
title('Entropy per Condition');
n = 0;
for s = 1:numel(spaces)
	for t = 1:numel(titles)
		n = n+1;
        
        % Reshape entropy
        d = nan(numel(entro{t,s}(:,:,1)), N{t,s}.conditions);
        for c = 1:N{t,s}.conditions
            d(:,c) = reshape(entro{t,s}(:,:,c), [numel(entro{t,s}(:,:,c)), 1]);
        end
        
        % Run comparisons
        a = nan(N{t,s}.comp, 1); p = nan(N{t,s}.comp, 1);
        for c = 1:N{t,s}.comp
    		[a(c), ~, ~] = kstest2(d(:,C(c,1)), d(:,C(c,2)));
        	[p(c), ~, ~] = permutationTest(d(:,C(c,1)), d(:,C(c,2)), 10000, 'sidedness','both');
        end
        
		% Get bin sizes
		f = figure; hold on;
        binWidth = nan(N{t,s}.conditions, 1);
        for c = 1:N{t,s}.conditions
            hg = histogram(d(:,c));
            binWidth(c) = hg.BinWidth;
        end
		sz = min(binWidth);
		close(f); clear hg binWidth f
		
		% Plot histograms
		subplot(numel(spaces),numel(titles),n); hold on;
        for c = 1:N{t,s}.conditions
            h{c} = histogram(entro{t,s}(:,:,c), 'Normalization','probability', 'BinWidth',sz);
        end
        legend(labels, 'location','northwest');
        
        % Plot signfiicance
        for c = 1:N{t,s}.comp
            if a(c) || p(c) < pTarget
                [mc(1),i(1)] = max(h{C(c,1)}.BinCounts); mc(1) = mc(1)/sum(h{C(c,1)}.BinCounts);
                [mc(2),i(2)] = max(h{C(c,2)}.BinCounts); mc(2) = mc(2)/sum(h{C(c,2)}.BinCounts);
                mp{1} = mean([h{C(c,1)}.BinEdges(i(1):i(1)+1), h{C(c,2)}.BinEdges(i(2):i(2)+1)]);
                mp{2} = [mean(h{C(c,1)}.BinEdges(i(1):i(1)+1)), mean(h{C(c,2)}.BinEdges(i(2):i(2)+1))];
                if a(c) && p(c) < pTarget
                    plot(mp{1}, 1.05*max(mc), '*g');
                    plot(mp{2}, 1.02*[max(mc),max(mc)], '-g');
                    plot(mp{2}(1), 1.02*[max(mc),max(mc)], '+g');
                    plot(mp{2}(2), 1.02*[max(mc),max(mc)], '+g');
                elseif p(c) < pTarget
                    plot(mp{1}, 1.05*max(mc), '*b');
                    plot(mp{2}, 1.02*[max(mc),max(mc)], '-b');
                    plot(mp{2}(1), 1.02*[max(mc),max(mc)], '+b');
                    plot(mp{2}(2), 1.02*[max(mc),max(mc)], '+b');
                else
                    plot(mp{1}, 1.05*max(mc), '*r');
                    plot(mp{2}, 1.02*[max(mc),max(mc)], '-r');
                    plot(mp{2}(1), 1.02*[max(mc),max(mc)], '+r');
                    plot(mp{2}(2), 1.02*[max(mc),max(mc)], '+r');
                end
            end
        end
		title([titles{t} ' Entropy in ' spaces{s} ' Space']);
		ylabel('Probability'); xlabel('Entropy');
		
		clear d
	end
end
clear s t h p a mc mp i n sz hg C


%% Compare component-level entropy between conditions

% Preallocate storage arrays
k = cell(numel(titles), numel(spaces), numel(ttype), N{1,1}.comp);
r = cell(numel(titles), numel(spaces), numel(ttype), N{1,1}.comp);
FDR = cell(numel(titles), numel(spaces));
Sidak = cell(numel(titles), numel(spaces));

% Set index locations
ticlocs = (N{1,1}.comp/5):(N{1,1}.comp/5):N{1,1}.comp;

% Extract thresholds at which components become FDR-significant
for e = 1:numel(titles)
	for s = 1:numel(spaces)
        
        % Determine all possible pairwise comparisons
        C = nchoosek(1:N{e,s}.conditions, 2);
        N{e,s}.comp = size(C,1);

		% Allocate threshold storage
		FDR{e,s} = nan(size(entro{e,s},1), numel(ttype), N{e,s}.comp, numel(prange));
		Sidak{e,s} = nan(size(entro{e,s},1), numel(ttype), N{e,s}.comp, numel(prange));
		
		% Test for differences
        for t = 1:numel(ttype)
			disp(['Running ', ttype{t}, ' tests on ', titles{e}, ' entropy in ', spaces{s}, ' space.']);
            for c = 1:N{e,s}.comp
                % Run comparisons
                disp(["Comparing ", labels(C(c,1)), " and ", labels(C(c,2))]);
                sig.comp(e,s,t,c) = robustTests(squeeze(entro{e,s}(:,:,C(c,1))), squeeze(entro{e,s}(:,:,C(c,2))), [], 'p',prange, 'testtype',ttype{t}); % robustTests(squeeze(entro{d}.IC(:,:,1)), squeeze(entro{d}.IC(:,:,2)), N{d}.IC, 'p',prange, 'testtype',ttype{t});
                
                % Extract significant components
                FDR{e,s}(:, t, c, :) = sig.comp(e,s,t,c).FDR(:,:);
                Sidak{e,s}(:, t, c, :) = sig.comp(e,s,t,c).Sidak(:,:);
                
               % Find components with 5% significance
                [r{e,s,t,c},~] = unique(find(squeeze(FDR{e,s}(:,t,c,prange==0.05))));
                [k{e,s,t,c},~] = unique(find(squeeze(Sidak{e,s}(:,t,c,prange==0.05))));
                k{e,s,t,c} = union(r{e,s,t,c}, k{e,s,t,c});
                disp([num2str(numel(k{e,s,t,c})), ' component(s) displays significant differences.']); 
            end
        end
	end
end

% Save significant component indices to files
if size(k,4) == 1
    for e = 1:numel(titles)
        h = cell2table(squeeze(k(e,:,:,:)), 'RowNames',ttype, 'VariableNames',spaces);
        % save(fullfile(path{4}, strcat(pfix{e}, fname)), 'h', '-append');
    end
else
    for e = 1:numel(titles)
        for s = 1:numel(spaces)
            h = cell2table(squeeze(k(e,s,:,:)), 'RowNames',ttype, 'VariableNames',labels(C'));
            % save(fullfile(path{4}, strcat(pfix{e}, spaces{s}, fname)), 'h', '-append');
        end
    end
end
h = k;
clear d e s t k


%% Plot significance results

% set color index for cortex network plots
cind = [1 0 0; 0 0 1];

F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(4) = 1;
F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(3) = 1;
F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(2) = 1;
F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(1) = 1;
kFig = 1;
for e = 1:numel(titles)
	for s = 1:numel(spaces)
		for t = 1:numel(ttype)
            for c = 1:N{e,s}.comp
                % Plot FDR significance level
                figure(F(nFig-4)); sgtitle('FDR Significance Threshold');
                ax(4) = subplot(numel(titles), numel(ttype)+numel(spaces), n(4)); 
                n(4) = n(4)+1;
                imagesc(ax(4), squeeze(FDR{e,s}(:,t,c,:)), [0 1]); colormap gray; colorbar;
                xticks(ax(4), 1:numel(prange));
                xticklabels(ax(4), strsplit(num2str(prange))); hold on;
                xlabel('FDR Significant'); ylabel('Component');
                title([titles{e}, ': ', spaces{s}, ' Space, ', ttype{t}, ' test']);

                % Plot Sidak significance level
                figure(F(nFig-3)); sgtitle('Sidak Significance Threshold');
                ax(3) = subplot(numel(titles), numel(ttype)+numel(spaces), n(3));
                n(3) = n(3)+1;
                imagesc(ax(3), squeeze(Sidak{e,s}(:,t,c,:)), [0 1]);
                colormap gray; colorbar;
                xticks(ax(3), 1:numel(prange));
                xticklabels(ax(3), strsplit(num2str(prange))); hold on;
                xlabel('Sidak Significant'); ylabel('Component');
                title([titles{e}, ': ', spaces{s}, ' Space, ', ttype{t}, ' test']);

                % Display test p-values
                figure(F(nFig-2)); sgtitle('p-values');
                ax(2) = subplot(numel(titles), numel(ttype)+numel(spaces), n(2)); 
                n(2) = n(2)+1;
                barh(ax(2), sig.comp(e,s,t,c).p); hold on;
                plot(ones(numel(h{e,s,t,c}),1), h{e,s,t,c}, '*r');
                title([titles{e}, ': ', ttype{t}, ' test, ', spaces{s}, ' space']);
                ylabel('Component'); xlabel('p-value');

                % Display test effect sizes
                figure(F(nFig-1)); sgtitle('Effect Sizes');
                ax(1) = subplot(numel(titles), numel(ttype)+numel(spaces), n(1));
                n(1) = n(1)+1;
                barh(ax(1), sig.comp(e,s,t,c).tstat); hold on;
                plot(sign(sig.comp(e,s,t,c).tstat(h{e,s,t,c})).*max(abs(sig.comp(e,s,t,c).tstat)).*ones(numel(h{e,s,t,c}),1), h{e,s,t,c}, '*r');
                title([titles{e}, ': ', ttype{t}, ', ', spaces{s}, ' space, ', labels(C(c,1)), ' vs. ', labels(C(c,2))]);
                ylabel('Component'); xlabel('Effect Size');

                if ~isempty(h{e,s,t,c})
                    for j = 1:numel(h{e,s,t,c})

                        % Get bin sizes
                        f = figure; hold on;
                        hg{1} = histogram(entro{e,s}(h{e,s,t,c}(j), :, C(c,1)));
                        hg{2} = histogram(entro{e,s}(h{e,s,t,c}(j), :, C(c,2)));
                        sz = min(hg{1}.BinWidth, hg{2}.BinWidth);
                        close(f);

                        % Plot results
                        if strncmpi(spaces{s}, 'IC', 2)

                            % Scale memberships (optional)
                            if zscale == true
                                mships = zscore(memberships{e});
                            elseif zthresh
                                mships = zscore(memberships{e});
                            else
                                mships = memberships{e};
                            end

                            if zthresh
                                for z = 1:numel(zthresh)
                                    K(kFig) = figure('Position', [0 0 1280 1024]); kFig = kFig + 1; hold on;

                                    % Connectivity
                                    kax = subplot(numel(h{e,s,t,c})*2, 5, [9 10]); hold on;
                                    sgtitle([titles{e}, ' Component ', num2str(h{e,s,t,c}(j)), ' ', labels(C(c,1)), ' vs. ', labels(C(c,2))], 'FontSize',18);
                                    a = squeeze(memberships{e}(:,h{e,s,t,c}(j)))*squeeze(memberships{e}(:,h{e,s,t,c}(j)))';
                                    imagesc(a); colorbar; hold on;
                                    xlim([1 size(a,2)]); ylim([1 size(a,1)]);
                                    yticks(1:N{e,s}.ROI); set(kax, 'YTickLabel',label_AAL90, 'FontSize', 3); xticks([]);
                                    title("Connectivity", 'FontSize',16); pbaspect([1 1 1]);

                                    % Histogram of component entropies
                                    kax = subplot(numel(h{e,s,t,c})*2, 5, [4 5]); hold on;
                                    histogram(entro{e,s}(h{e,s,t,c}(j), :, C(c,1)), 'BinWidth',sz, 'Normalization','Probability');
                                    histogram(entro{e,s}(h{e,s,t,c}(j), :, C(c,2)), 'BinWidth',sz, 'Normalization','Probability');
                                    legend(labels(C(c,:)));
                                    title("Entropy", 'FontSize',16);
                                    ylabel('Counts'); xlabel('Mean Entropy');

                                    % Brain Renderings
                                    kax = subplot(numel(h{e,s,t,c})*2, 5, [2 3 7 8]); hold on;
                                    plot_nodes_in_cortex(cortex, mships(:,h{e,s,t,c}(j)), coords_AAL90, origin, sphereScale, zthresh(z), [], cind, [], [], rdux);

                                    % Bar Plots
                                    kax = subplot(numel(h{e,s,t,c})*2, 5, [1 6]); hold on;
                                    if sum(sign(squeeze(memberships{e}(:,h{e,s,t,c}(j))))) >= 0
                                        ind(:,1) = squeeze(mships(:,h{e,s,t,c}(j))) < -zthresh(z);	% select node weights which surpass threshold
                                        ind(:,2) = squeeze(mships(:,h{e,s,t,c}(j))) > zthresh(z);	% select node weights which surpass threshold
                                        ind(:,3) = squeeze(mships(:,h{e,s,t,c}(j))) > -zthresh(z) & squeeze(mships(:,h{e,s,t,c}(j))) < 0;
                                        ind(:,4) = squeeze(mships(:,h{e,s,t,c}(j))) < zthresh(z) & squeeze(mships(:,h{e,s,t,c}(j))) > 0;
                                        a = squeeze(mships(:,h{e,s,t,c}(j))); a(~ind(:,1)) = 0; barh(1:N{1,1}.comp, a, 'b');
                                        a = squeeze(mships(:,h{e,s,t,c}(j))); a(~ind(:,2)) = 0; barh(1:N{1,1}.comp, a, 'r');
                                        a = squeeze(mships(:,h{e,s,t,c}(j))); a(~ind(:,3)) = 0; barh(1:N{1,1}.comp, a, 'b', 'FaceAlpha',0.3);
                                        a = squeeze(mships(:,h{e,s,t,c}(j))); a(~ind(:,4)) = 0; barh(1:N{1,1}.comp, a, 'r', 'FaceAlpha',0.3);
                                    elseif sum(sign(squeeze(memberships{e}(:,h{e,s,t,c}(j))))) < 0
                                        ind(:,1) = squeeze(mships(:,h{e,s,t,c}(j))) > zthresh(z);	% only plot node weights which surpass threshold
                                        ind(:,2) = squeeze(mships(:,h{e,s,t,c}(j))) < -zthresh(z);	% only plot node weights which surpass threshold
                                        ind(:,3) = squeeze(mships(:,h{e,s,t,c}(j))) < zthresh(z) & squeeze(mships(:,h{e,s,t,c}(j))) > 0;
                                        ind(:,4) = squeeze(mships(:,h{e,s,t,c}(j))) > -zthresh(z) & squeeze(mships(:,h{e,s,t,c}(j))) < 0;
                                        a = squeeze(mships(:,h{e,s,t,c}(j))); a(~ind(:,1)) = 0; barh(1:N{1,1}.comp, a, 'b');
                                        a = squeeze(mships(:,h{e,s,t,c}(j))); a(~ind(:,2)) = 0; barh(1:N{1,1}.comp, a, 'r');
                                        a = squeeze(mships(:,h{e,s,t,c}(j))); a(~ind(:,3)) = 0; barh(1:N{1,1}.comp, a, 'b', 'FaceAlpha',0.3);
                                        a = squeeze(mships(:,h{e,s,t,c}(j))); a(~ind(:,4)) = 0; barh(1:N{1,1}.comp, a, 'r', 'FaceAlpha',0.3);
                                    end
                                    yticks(1:N{e,s}.ROI); set(kax, 'YTickLabel',label_AAL90, 'FontSize', 6);
                                    title("Component Membership", 'FontSize',16, 'Position',[5 93]);
                                    subtitle(['z-score threshold: ' num2str(zthresh(z))], 'FontSize',12, 'Position',[5 91]);
                                    xlabel('z-score', 'FontSize',12);

                                    % Save as png file
                                    % saveas(K(kFig-1), strcat('entroCompare', fname, '_Z', strjoin(split(num2str(zthresh(z)), '.')', '')), 'png');
                                end
                            else
                                K(kFig) = figure('Position', [0 0 1280 1024]); kFig = kFig + 1; hold on;

                                % Connectivity
                                kax = subplot(numel(h{e,s,t,c})*2, 5, [9 10]); hold on;
                                sgtitle([titles{e}, ' Component ', num2str(h{e,s,t,c}(j))], 'FontSize',18);
                                a = squeeze(memberships{e}(:,h{e,s,t,c}(j)))*squeeze(memberships{e}(:,h{e,s,t,c}(j)))';
                                imagesc(a); colorbar; hold on;
                                xlim([1 size(a,2)]); ylim([1 size(a,1)]);
                                yticks(1:N{e,s}.ROI); set(kax, 'YTickLabel',label_AAL90, 'FontSize', 3); xticks([]);
                                title("Connectivity", 'FontSize',16); pbaspect([1 1 1]);

                                % Histogram of component entropies
                                kax = subplot(numel(h{e,s,t,c})*2, 5, [4 5]); hold on;
                                histogram(entro{e,s}(h{e,s,t,c}(j), :, C(c,1)), 'BinWidth',sz, 'Normalization','Probability');
                                histogram(entro{e,s}(h{e,s,t,c}(j), :, C(c,2)), 'BinWidth',sz, 'Normalization','Probability');
                                legend(labels(C(c,:)));
                                title("Entropy", 'FontSize',16);
                                ylabel('Counts'); xlabel('Mean Entropy');

                                % Bar Plots
                                kax = subplot(numel(h{e,s,t,c})*2, 5, [1 6]); hold on;
                                if sum(sign(squeeze(memberships{e}(:,h{e,s,t,c}(j))))) >= 0
                                    ind(:,1) = (squeeze(mships(:,h{e,s,t,c}(j))) < 0);
                                    ind(:,2) = (squeeze(mships(:,h{e,s,t,c}(j))) > 0);
                                    a = squeeze(mships(:,h{e,s,t,c}(j))); a(ind(:,1)) = 0; barh(1:N{1,1}.comp, a);
                                    a = squeeze(mships(:,h{e,s,t,c}(j))); a(ind(:,2)) = 0; barh(1:N{1,1}.comp, a, 'r');
                                elseif sum(squeeze(memberships{e}(:,h{e,s,t,c}(j)))) < 0
                                    ind(:,1) = (squeeze(mships(:,h{e,s,t,c}(j))) > 0);
                                    ind(:,2) = (squeeze(mships(:,h{e,s,t,c}(j))) < 0);
                                    a = squeeze(mships(:,h{e,s,t,c}(j))); a(ind(:,1)) = 0; barh(1:N{1,1}.comp, a);
                                    a = squeeze(mships(:,h{e,s,t,c}(j))); a(ind(:,2)) = 0; barh(1:N{1,1}.comp, a, 'r');
                                end
                                yticks(1:N{e,s}.ROI); set(kax, 'YTickLabel',label_AAL90, 'FontSize', 6);
                                title("Component Membership", 'FontSize',16, 'Position',[5 93]);
                                subtitle(['z-score threshold: ' num2str(zthresh(z))], 'FontSize',12, 'Position',[5 91]);
                                xlabel('z-score', 'FontSize',12);

                                % Brain Renderings
                                kax = subplot(numel(h{e,s,t,c})*2, 5, [2 3 7 8]); hold on;
                                plot_nodes_in_cortex(cortex, mships(:,h{e,s,t,c}(j)), coords_AAL90, origin, sphereScale, [], rdux);

                                % Save as png file
                                % saveas(K(kFig-1), strcat('entroCompare', fname, '_Z', join(split(num2str(zthresh(z)), '.'))), 'png');
                            end
                        else
                            % Histogram of component entropies
                            K(kFig) = figure('Position', [0 0 1280 1024]); kFig = kFig + 1; hold on;
                            histogram(entro{e,s}(h{e,s,t,c}(j), :, C(c,1)), 'BinWidth',sz, 'Normalization','Probability');
                            histogram(entro{e,s}(h{e,s,t,c}(j), :, C(c,2)), 'BinWidth',sz, 'Normalization','Probability');
                            legend(labels(C(c,:)));
                            title(strjoin({'Entropy of', char(label_AAL90(h{e,s,t,c}(j))), 'in', spaces{s}, 'Space'}));
                            ylabel('Counts'); xlabel('Entropy');

                            % Save as png file
                            % saveas(K(kFig-1), strcat('entroCompare', fname, '_Component', num2str(h{e,s,t,c}(j))), 'png');
                        end
                    end
                end
            end
        end
	end
end
clear e s t p B n ax n k j a sz hg f mships m ind

% Add kFigs to F
for k = 1:kFig-1
	F(nFig) = K(k); nFig = nFig + 1;
end
clear K k kFig



%% Compare distances between entropy distributions
% 
% esubdist = cell(numel(titles), numel(spaces), numel(dim), 4);
% for s = 1:numel(spaces)
% 	F(nFig) = figure; nFig = nFig + 1; hold on;
% 	F(nFig) = figure; nFig = nFig + 1; hold on;
% 	title('Euclidean Distances Between Subject Entropies');
% 	for t = 1:numel(titles)
% 		e = nan(N{t,s}.comp, sum(N{t,s}.subjects));
% 		e(:, 1:N{t,s}.subjects(1)) = squeeze(entro{t,s}(:, 1:N{t,s}.subjects(1), 1));
% 		for c = 2:N{t,s}.conditions
% 			e(:, N{t,s}.subjects(c-1)+1:sum(N{t,s}.subjects(1:c))) = squeeze(entro{t,s}(:, 1:N{t,s}.subjects(c), c));
% 		end
% 		
% 		for d = 1:numel(dim)
% 			e = e';
% 			if strcmpi(dim{d}, 'Subject')
% 				C = nchoosek(1:sum(N{t,s}.subjects), 2);
% 				esubdist{t,s,d,2} = zeros(sum(N{t,s}.subjects), sum(N{t,s}.subjects));
% 				esubdist{t,s,d,3} = zeros(sum(N{t,s}.subjects), sum(N{t,s}.subjects));
% 				esubdist{t,s,d,4} = zeros(sum(N{t,s}.subjects), sum(N{t,s}.subjects));
% 			else
% 				C = nchoosek(1:N{t,s}.comp, 2);
% 				esubdist{t,s,d,2} = zeros(N{t,s}.comp, N{t,s}.comp);
% 				esubdist{t,s,d,3} = zeros(N{t,s}.comp, N{t,s}.comp);
% 				esubdist{t,s,d,4} = zeros(N{t,s}.comp, N{t,s}.comp);
% 			end
% 			
% 			% Euclidean distances
% 			esubdist{t,s,d,1} = squareform(pdist(e));
% 			figure(F(nFig-d)); subplot(3,numel(titles),t); colormap jet; hold on;
% 			imagesc(esubdist{t,s,d,1}); colorbar;
% 			xlim([1 size(esubdist{t,s,d,1}, 2)]); ylim([1 size(esubdist{t,s,d,1}, 1)]);
% 			title([titles{t},' Euclidean Distance, ', spaces{s}, ' Space']);
% 			xlabel(dim{d}); ylabel(dim{d});
% 			
% 			% KS distances & significance
% 			for i = 1:size(C,1)
% 				[esubdist{t,s,d,3}(C(i,2), C(i,1)), ~, esubdist{t,s,d,2}(C(i,1), C(i,2))] = kstest2(e(C(i,1),:), e(C(i,2),:));
% 			end
% 			esubdist{t,s,d,2} = esubdist{t,s,d,2}+esubdist{t,s,d,2}';
% 			esubdist{t,s,d,3} = esubdist{t,s,d,3}+esubdist{t,s,d,3}';
% 			figure(F(nFig-d)); subplot(3,numel(titles),t+numel(titles)); colormap jet; hold on;
% 			imagesc(esubdist{t,s,d,2}); colorbar; caxis([0 1]);
% 			xlim([1 size(esubdist{t,s,d,2}, 2)]); ylim([1 size(esubdist{t,s,d,2}, 1)]);
% 			title([titles{t},' KS Distance, ', spaces{s}, ' Space']);
% 			xlabel(dim{d}); ylabel(dim{d});
% 			
% 			% Absolute mean differences
% 			em = mean(e, 2);
% 			for i = 1:size(C,1)
% 				esubdist{t,s,d,4}(C(i,2), C(i,1)) = abs(em(C(i,1)) - em(C(i,2)));
% 			end
% 			esubdist{t,s,d,4} = esubdist{t,s,d,4}+esubdist{t,s,d,4}';
% 			figure(F(nFig-d)); subplot(3,numel(titles),t+2*numel(titles)); colormap jet; hold on;
% 			imagesc(esubdist{t,s,d,4}); colorbar; caxis([0 1]);
% 			xlim([1 size(esubdist{t,s,d,4}, 2)]); ylim([1 size(esubdist{t,s,d,4}, 1)]);
% 			title([titles{t},' Mean Distance, ', spaces{s}, ' Space']);
% 			xlabel(dim{d}); ylabel(dim{d});
% 		end
% 	end
% end
% clear t s d i C e em j c


%% Save results
% 
% savefig(F, strcat('entroCompare', fname), 'compact');
% clear F nFig
% save(strcat('entroCompare', fname));