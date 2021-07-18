%% Import data

clear; close all; clc;

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1}(1:end-2),'/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{1,1} = strjoin(path{1}(1:end-3),'/');

% Set required subdirectories
path{4,1} = fullfile(path{2}, 'UCLA', 'Results','LEICA');
path{5,1} = fullfile(path{1}, 'Project','Atlases','AAL');

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
zthresh = 1:0.2:1.6;
zscale = false;

% Determine which files to compare
band = 'WideBand';
k = 1;
iteration = 1;   % 2;
fname = strcat('_',band, '_k',num2str(k), '_Iteration',num2str(iteration));
pfix = "GroupIC/Pairwise/LE_ICA_ControlIC_COS_CONTROLvSCHZ";  % {'LE_ICA_AAL90_CIC_COS', 'M_ICA_AAL90_CIC_EXP'};
clear assemblies distance band k iteration

% Set decompositions, spaces, comparisions
nFig = 1;
titles = "LEICA"; %  {'LEICA', 'MICA'};	% compression type
spaces = {'dFC' 'IC'};				% space in which to compare
dim = {'Subject', 'Component'};		% dimensions to compare
ttype = {'kstest2', 'permutation'};	% set test types to run
pTarget = 0.05;						% target p-value
prange = 0.025:0.025:0.1;			% Set range of p-values to test

% Load entropies
N = cell(numel(pfix), numel(spaces));
entro = cell(numel(pfix), numel(spaces));
FCD = cell(numel(pfix), numel(spaces));
memberships = cell(numel(pfix), 1);	% have 2 dFC / FCD spaces: ROI and IC.
C = cell(numel(pfix),1);
comps = cell(numel(pfix),1);
for f = 1:numel(pfix)
	e = load(fullfile(path{4}, strcat(pfix{f}, fname)), 'entro','N','memberships','FCD','I','T','ROI','C');
	labels = e.I.Properties.VariableNames;
	T = e.T;
	N{f,1} = e.N; N{f,1}.IC = N{f,1}.ROI;   % N{f,1} = rmfield(N{f,1}, 'ROI');
	% N{f,2} = N{f,1};
	N{f,2} = e.N;   % N{f,2} = rmfield(N{f,2}, 'ROI');
	% entro{f,1} = e.entro.BOLD;
	entro{f,1} = e.entro.dFC;
	entro{f,2} = e.entro.IC;
	memberships{f,1} = e.memberships;
	FCD{f,1} = e.FCD.dFC;
	FCD{f,2} = e.FCD.IC;
    ROI = e.ROI;
    C{f} = e.C;
    comps{f} = unique(e.C);
end
clear e f

% Load network labels
labels_ROI = string(ROI.Properties.RowNames);
% labels_ROI = load(fullfile(path{5}, 'AAL_labels'));
% labels_ROI = string(labels_ROI.label90);
% labels_ROI = strip(LR_version_symm(labels_ROI));

% Generate MNI coordinates
coords_ROI = horzcat(ROI{:,'x'}, ROI{:,'y'}, ROI{:,'z'});
% coords_ROI = load(fullfile(path{5}, 'aal_cog.txt'), 'aal_cog');
% coords_ROI = LR_version_symm(coords_ROI);
origin = [65 45.5 35];							% center origin
MNIscale = 5.5/10;								% scale MNI coordinates
sphereScale = 2.5;								% scale sphere size
coords_ROI = MNIscale*coords_ROI;			% scale MNI coordinates
clear label90 MNIscale

% Set brain rendering parameters
cortex.file = fullfile(path{5},'MNI152_T1_2mm_brain_mask.nii');	% file containing cortical atlas
cortex.color = [0.9 0.9 0.9];					% color for cortical rendering
cortex.transparency = 0.1;						% set to 1 for opaque cortex
cortex.val = 0.3;								% set isonormal line spacing
cortex.view = [-90 90];							% set camera angle
rdux = 0.7;						% proportion of surface faces to keep (<1)


%% Display FCDs

% Preallocate array to store KS distances
ksdist.FCD = nan(horzcat(size(FCD), size(C{1},1)));

% Set pre-ICA titles
ttls = {'Leading Eigenvector', 'Mean', 'Standard'};	% compression type

for s = 1:size(FCD,2)
	F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1;
	for t = 1:size(FCD,1)

		% Plot examples of patient, control FCDs
		for l = comps{t}    % 1:numel(comps{t})
			subplot(size(FCD,1), numel(comps{t})+1, find(comps{t}==l)+(numel(C{t})+1)*(t-1));
			colormap jet;
			ind = find(~mod(round(T.TR*1:T.scan), 25));
			imagesc(FCD{t,s}.subj{1,l}); colorbar; caxis([-1 1]);
			xticks([]);
			yticks(ind); yticklabels(round(T.TR.*ind));
			title([ttls{t}, ', ', labels{l}, ' 1']);
			xlabel('Time (s)'); ylabel('Time (s)');
            pbaspect([1 1 1]);
		end
		
		% Measure KS distance between FCD histograms
		D = cell(size(comps{t}));
		grp = D;
        for f = comps{t}
			D{comps{t}==f} = reshape(cell2mat(FCD{t,s}.subj(:,f)), [numel(cell2mat(FCD{t,s}.subj(:,f))), 1]);
			grp{comps{t}==f} = string(repmat(labels{f}, [numel(cell2mat(FCD{t,s}.subj(:,f))), 1]));
        end
        for c = 1:size(C{t},1)
            a = D{comps{t}==C{t}(c,1)}(:);
            b = D{comps{t}==C{t}(c,2)}(:);
            [ksdist.h(t,s,c), ksdist.p(t,s,c), ksdist.FCD(t,s,c)] = kstest2(a(isfinite(a)), b(isfinite(b)));
        end
		% d(t,s) = mean(D{1}(:)) - mean(D{2}(:));
		
		% Plot patient, control FCD histograms
		subplot(size(FCD,1), numel(C{t})+1, (numel(C{t})+1)*t); hold on;
		grp = vertcat(grp{:});
		D = vertcat(D{:});	% cell2mat(D);
		boxplot(D,grp, 'Notch','on');
		
		% histogram(cell2mat(FCD{t,s}.subj(:,1))', 'Normalization','probability', 'FaceAlpha',0.2, 'FaceColor',[1 0 0]);
		% histogram(cell2mat(FCD{t,s}.subj(:,2))', 'Normalization','probability', 'FaceAlpha',0.2, 'FaceColor',[0 0 1]);
		% legend(labels, 'location','northwest');
		title([ttls{t}, ' FCD Values']);
		ylabel('Correlation');	% ylabel('Probability');
	end
	sgtitle(['FCD in ', spaces{s}, ' Space']);
end
clear s t l ind con pat ttls grp D d a b


%% Compare mean entropy across methods

% Preallocate arrays
mEntro = cell(size(entro,1), numel(spaces), ndims(entro{1,1})-1);
cond = cell(1, N{1}.conditions);

% Set comparison order
dType = {'Subject Mean' 'IC Mean'};

for t = 1:numel(titles)
	for s = 1:numel(spaces)
		for d = 1:numel(dType)
			
			% Compute means
			mEntro{t,s,d} = squeeze(mean(entro{t,s}, d, 'omitnan'));
			mEntro{t,s,d} = array2table(mEntro{t,s, d}, 'VariableNames', labels);
			
			% Pairwise comparisons between conditions
            for c = 1:N{t,s}.comp
                for p = 1:size(C{t},2)
                    cond{p} = mEntro{t,s,d}{:,labels{C{t}(c,p)}}(isfinite(mEntro{t,s,d}{:,labels{C{t}(c,p)}}));
                end
                [sig.av.h(t,s,d,c,1), sig.av.p(t,s,d,c,1), sig.av.effsize(t,s,d,c,1)] = kstest2(cond{1}, cond{2});
                [sig.av.p(t,s,d,c,2), ~, sig.av.effsize(t,s,d,c,2)] = permutationTest(cond{1}, cond{2}, 10000, 'sidedness','both');
                if sig.av.p(t,s,d,c,2) < 0.05
                    sig.av.h(t,s,d,c,2) = 1;
                else
                    sig.av.h(t,s,d,c,2) = 0;
                end
            end
            
            % Run mutliple-comparison corrections for conditions
            if c > 1
                [sig.av.FDR(t,s,d,:,1), sig.av.Bonferroni(t,s,d,:,1), sig.av.Sidak(t,s,d,:,1)] = mCompCorr([], squeeze(sig.av.p(t,s,d,:,1)), pTarget);
                [sig.av.FDR(t,s,d,:,2), sig.av.Bonferroni(t,s,d,:,2), sig.av.Sidak(t,s,d,:,2)] = mCompCorr([], squeeze(sig.av.p(t,s,d,:,2)), pTarget);
            end
		end
	end 
end
clear s t d

% Split into comparisons
sig.mSubj.h = squeeze(sig.av.h(:,:,1,:,:));
sig.mSubj.p = squeeze(sig.av.p(:,:,1,:,:));
sig.mSubj.effsize = squeeze(sig.av.effsize(:,:,1,:,:));
sig.mIC.h = squeeze(sig.av.h(:,:,2,:,:));
sig.mIC.p = squeeze(sig.av.p(:,:,2,:,:));
sig.mIC.effsize = squeeze(sig.av.effsize(:,:,2,:,:));
if c > 1
    sig.mSubj.FDR = squeeze(sig.av.FDR(:,:,1,:,:));
    sig.mSubj.Bonferroni = squeeze(sig.av.Bonferroni(:,:,1,:,:));
    sig.mSubj.Sidak = squeeze(sig.av.Sidak(:,:,1,:,:));
    sig.mIC.FDR = squeeze(sig.av.FDR(:,:,2,:,:));
    sig.mIC.Bonferroni = squeeze(sig.av.Bonferroni(:,:,2,:,:));
    sig.mIC.Sidak = squeeze(sig.av.Sidak(:,:,2,:,:));
end
clear c

% Display mean entropy histograms
F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; hold on;
col = 0;
for d = 1:numel(dType)			% dimension average to display
	for s = 1:numel(spaces)		% space to display (dFC or IC)
		col = col+1;
		for t = 1:numel(titles)
			
			% Plot mean entropy histograms
			nc = numel(spaces)*(ndims(entro{t,s})-1);
			row = nc*(t-1);
			subplot(size(entro,1), nc, col+row); hold on;
            for c = 1:N{t,s}.conditions
                h{c} = histogram(mEntro{t,s,d}{:,labels{c}}, 'Normalization','pdf', 'FaceAlpha',0.4);
            end
			
			% Plot means, distances for significant differences
            for c = 1:N{t,s}.comp
                if sig.av.h(t,s,d,c,1) || sig.av.h(t,s,d,c,2)
                    [mc(1),i(1)] = max(h{C{t}(c,1)}.BinCounts); mc(1) = mc(1)/sum(h{C{t}(c,1)}.BinCounts);
                    [mc(2),i(2)] = max(h{C{t}(c,2)}.BinCounts); mc(2) = mc(2)/sum(h{C{t}(c,2)}.BinCounts);
                    mp{1} = mean([h{C{t}(c,1)}.BinEdges(i(1):i(1)+1), h{C{t}(c,2)}.BinEdges(i(2):i(2)+1)]);
                    mp{2} = [mean(h{C{t}(c,1)}.BinEdges(i(1):i(1)+1)), mean(h{C{t}(c,2)}.BinEdges(i(2):i(2)+1))];
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
clear d s t nc col row mc mp i h hg


%% Overall Entropy

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
        a = nan(N{t,s}.comp, 2);
        for c = 1:N{t,s}.comp
    		[~, a(c,1), ~] = kstest2(d(:,C{t}(c,1)), d(:,C{t}(c,2)));
        	[a(c,2), ~, ~] = permutationTest(d(:,C{t}(c,1)), d(:,C{t}(c,2)), 10000, 'sidedness','both');
        end
            
        % Run mutliple-comparison corrections for conditions
        if c > 1
            [FDR(:,1), Bonferroni(:,1), Sidak(:,1)] = mCompCorr([], a(:,1), pTarget);
            [FDR(:,2), Bonferroni(:,2), Sidak(:,2)] = mCompCorr([], a(:,2), pTarget);
        end
		
		% Plot histograms
		subplot(numel(spaces),numel(titles),n); hold on;
        for c = 1:N{t,s}.conditions
            h{c} = histogram(entro{t,s}(:,:,c), 'Normalization','pdf');
        end
        legend(labels, 'location','northwest');
        
        % Plot significance
        for c = 1:N{t,s}.comp
            if exist('FDR','var') & FDR(c,:)
                [mc(1),i(1)] = max(h{C{t}(c,1)}.BinCounts); mc(1) = mc(1)/sum(h{C{t}(c,1)}.BinCounts);
                [mc(2),i(2)] = max(h{C{t}(c,2)}.BinCounts); mc(2) = mc(2)/sum(h{C{t}(c,2)}.BinCounts);
                mp{1} = mean([h{C{t}(c,1)}.BinEdges(i(1):i(1)+1), h{C{t}(c,2)}.BinEdges(i(2):i(2)+1)]);
                mp{2} = [mean(h{C{t}(c,1)}.BinEdges(i(1):i(1)+1)), mean(h{C{t}(c,2)}.BinEdges(i(2):i(2)+1))];
                if a(c,1) && a(c,2) < pTarget
                    plot(mp{1}, 1.05*max(mc), '*g');
                    plot(mp{2}, 1.02*[max(mc),max(mc)], '-g');
                    plot(mp{2}(1), 1.02*[max(mc),max(mc)], '+g');
                    plot(mp{2}(2), 1.02*[max(mc),max(mc)], '+g');
                elseif a(c,1) < pTarget
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
clear s t h p a mc mp i n hg


%% Compare component-level entropy between conditions

% Preallocate storage arrays
k = cell(numel(titles), numel(spaces), numel(ttype), N{1,1}.comp);
r = cell(numel(titles), numel(spaces), numel(ttype), N{1,1}.comp);
FDR = cell(numel(titles), numel(spaces));
Sidak = cell(numel(titles), numel(spaces));
p = cell(numel(titles), numel(spaces));
h = cell(numel(titles), numel(spaces));
tstat = cell(numel(titles), numel(spaces));

% Set index locations
ticlocs = (N{1,1}.comp/5):(N{1,1}.comp/5):N{1,1}.comp;

% Extract thresholds at which components become FDR-significant
for e = 1:numel(titles)
	for s = 1:numel(spaces)
		% Allocate threshold storage
		FDR{e,s} = nan(N{e,s}.IC, numel(ttype), N{e,s}.comp, numel(prange));
		Sidak{e,s} = nan(N{e,s}.IC, numel(ttype), N{e,s}.comp, numel(prange));
        p{e,s} = nan(N{e,s}.IC, numel(ttype), N{e,s}.comp);
        tstat{e,s} = nan(N{e,s}.IC, numel(ttype), N{e,s}.comp);
        h{e,s} = nan(N{e,s}.IC, numel(ttype), N{e,s}.comp);
		
		% Test for differences
        for t = 1:numel(ttype)
			disp(['Running ', ttype{t}, ' tests on ', titles{e}, ' entropy in ', spaces{s}, ' space.']);
            for c = 1:N{e,s}.comp
                % Run comparisons
                disp(["Comparing ", labels(C{e}(c,1)), " and ", labels(C{e}(c,2))]);
                [h{e,s}(:,t,c), p{e,s}(:,t,c), tstat{e,s}(:,t,c)] = robustTests(squeeze(entro{e,s}(:,:,C{e}(c,1))), squeeze(entro{e,s}(:,:,C{e}(c,2))), [], 'p',prange, 'testtype',ttype{t});
            end
            
            % Multiple comparison correction
            if N{e,s}.comp > 1
                p_dum = reshape(squeeze(p{e,s}(:,t,:)), [N{e,s}.comp*N{e,s}.IC, 1]);	% reshape p-value array
                [f, ~, S] = mCompCorr([], p_dum, prange);                        % run mutliple comparison tests
                FDR{e,s}(:,t,:,:) = reshape(f, [N{e,s}.IC, N{e,s}.comp, numel(prange)]);
                Sidak{e,s}(:,t,:,:) = reshape(S, [N{e,s}.IC, N{e,s}.comp, numel(prange)]);
                
                % Find components with 5% significance
                for c = 1:N{e,s}.comp
                    [r{e,s,t,c},~] = unique(find(squeeze(FDR{e,s}(:,t,c,prange == pTarget))));
                    [k{e,s,t,c},~] = unique(find(squeeze(Sidak{e,s}(:,t,c,prange == pTarget))));
                    k{e,s,t,c} = union(r{e,s,t,c}, k{e,s,t,c});
                    disp([num2str(numel(k{e,s,t,c})), ' component(s) displays significant differences.']); 
                end
            end
        end
	end
end
clear c f S p_dum

% Save significant component indices to files
vn = cell(size(C{1},1), 1);
for c = 1:size(C{1},1)
    vn{c} = [labels{C{1}(c,1)}, ' v. ', labels{C{1}(c,2)}];
end
if N{e,s}.comp > 1
    for e = 1:numel(titles)
        for s = 1:numel(spaces)
            h{e,s} = cell2table(squeeze(k(e,s,:,:)), 'RowNames',ttype, 'VariableNames',vn);
        end
    end
    save(fullfile(path{4}, strcat(pfix{e}, fname)), 'h', '-append');
    h = k;
else
    save(fullfile(path{4}, strcat(pfix{e}, fname)), 'h', '-append');
    k = h; h = cell(numel(titles), numel(spaces), numel(ttype), 1);
    for e = 1:numel(titles)
        for s = 1:numel(spaces)
            for t = 1:numel(ttype)
                h{e,s,t,1} = find(k{e,s}(:,t,1));
            end
        end
    end
end
clear d e s t k r c vn


%% Plot significance results

% set color index for cortex network plots
cind.node = [1 0 0; 0 0 1];
cind.conn = [1 0 0; 0 0 1];

kFig = 1;
for c = 1:N{1,1}.comp
    F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(4) = 1;
    F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(3) = 1;
    F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(2) = 1;
    F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(1) = 1;
	for e = 1:numel(titles)
        for s = 1:numel(spaces)
            
            % Set index for all components found across all tests
            k = [];
            
            % Plot significance test results
            for t = 1:numel(ttype)
                
                % Set index for all components found across all tests
                k = union(k, h{e, s, t, c});
                
                if N{e,s}.comp > 1
                    % Plot FDR significance level
                    figure(F(nFig-2));
                    sgtitle([labels{C{e}(c,1)}, ' vs. ', labels{C{e}(c,2)}]);	% sgtitle('FDR Significance Threshold');
                    ax(4) = subplot(numel(titles), numel(ttype)+numel(spaces), n(4)); 
                    n(4) = n(4)+1;
                    imagesc(ax(4), squeeze(FDR{e,s}(:,t,c,:)), [0 1]); colormap gray; colorbar;
                    xticks(ax(4), 1:numel(prange));
                    xticklabels(ax(4), strsplit(num2str(prange))); hold on;
                    xlabel('FDR Significant'); ylabel('Component');
                    title([titles{e}, ': ', ttype{t}, ', ', spaces{s}, ' space']);

                    % Plot Sidak significance level
                    figure(F(nFig-1));
                    sgtitle([labels{C{e}(c,1)}, ' vs. ', labels{C{e}(c,2)}]);	% sgtitle('Sidak Significance Threshold');
                    ax(3) = subplot(numel(titles), numel(ttype)+numel(spaces), n(3));
                    n(3) = n(3)+1;
                    imagesc(ax(3), squeeze(Sidak{e,s}(:,t,c,:)), [0 1]);
                    colormap gray; colorbar;
                    xticks(ax(3), 1:numel(prange));
                    xticklabels(ax(3), strsplit(num2str(prange))); hold on;
                    xlabel('Sidak Significant'); ylabel('Component');
                    title([titles{e}, ': ', ttype{t}, ', ', spaces{s}, ' space']);    % , labels{C{e}(c,1)}, ' vs. ', labels{C{e}(c,2)}]);
                end

                % Display test p-values
                figure(F(nFig-4));
                sgtitle([labels{C{e}(c,1)}, ' vs. ', labels{C{e}(c,2)}]);	% sgtitle('p-values');
                ax(2) = subplot(numel(titles), numel(ttype)+numel(spaces), n(2)); 
                n(2) = n(2)+1;
                barh(ax(2), p{e,s}(:,t,c)); hold on;
                plot(ones(numel(h{e,s,t,c}),1), h{e,s,t,c}, '*r');
                title([titles{e}, ': ', ttype{t}, ', ', spaces{s}, ' space']);    % , labels{C{e}(c,1)}, ' vs. ', labels{C{e}(c,2)}]);
                ylabel('Component'); xlabel('p-value');

                % Display test effect sizes
                figure(F(nFig-3));
                sgtitle([labels{C{e}(c,1)}, ' vs. ', labels{C{e}(c,2)}]);	% sgtitle('Effect Sizes');
                ax(1) = subplot(numel(titles), numel(ttype)+numel(spaces), n(1));
                n(1) = n(1)+1;
                barh(ax(1), tstat{e,s}(:,t,c)); hold on;
                plot(sign(tstat{e,s}(h{e,s,t,c},t,c)).*max(abs(tstat{e,s}(:,t,c))).*ones(numel(h{e,s,t,c}),1), h{e,s,t,c}, '*r');
                title([titles{e}, ': ', ttype{t}, ', ', spaces{s}, ' space']);    %, labels{C{e}(c,1)}, ' vs. ', labels{C{e}(c,2)}]);
                ylabel('Component'); xlabel('Effect Size');
            end
            
            % Plot entropy distributions for signifcant separations
            if ~isempty(k)
                for j = 1:numel(k)

                    % Plot results
                    if strncmpi(spaces{s}, 'IC', 2)

                        % Scale memberships (optional)
                        if zscale == true || sum(zthresh ~= 0, 'all') > 0
                            mships = squeeze(zscore(memberships{e}(:,k(j))));
                        else
                            mships = squeeze(memberships{e}(:,k(j)));
                        end
                        
                        % Enforce consensus: smaller community should be "positive"
                        if sum(sign(mships)) > 0
                            mships = -mships;
                        end
                        
                        % Plot component
                        if sum(zthresh ~= 0, 'all') > 0
                            for z = 1:numel(zthresh)
                                K(kFig) = figure('Position', [0 0 1280 1024]); kFig = kFig + 1; hold on;

                                % Connectivity
                                kax = subplot(2, 5, [9 10]); hold on;   % subplot(numel(h{e,s,t,c})*2, 5, [9 10]); hold on;
                                sgtitle([titles{e}, ' Component ', num2str(k(j)), ' ', labels{C{e}(c,1)}, ' vs. ', labels{C{e}(c,2)}], 'FontSize',18);
                                a = squeeze(memberships{e}(:,k(j)))*squeeze(memberships{e}(:,k(j)))';
                                imagesc(a); colorbar; hold on;
                                xlim([1 size(a,2)]); ylim([1 size(a,1)]);
                                yticks(1:N{e,s}.ROI); set(kax, 'YTickLabel',labels_ROI, 'FontSize', 3); xticks([]);
                                title("Connectivity", 'FontSize',16); pbaspect([1 1 1]);

                                % Histogram of component entropies
                                kax = subplot(2, 5, [4 5]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [4 5]); hold on;
                                for f = 1:length(C{e}(c,:))
                                    histogram(entro{e,s}(k(j), :, C{e}(c,f)), 'Normalization','pdf');
                                end
                                legend(labels(C{e}(c,:)));
                                title("Entropy", 'FontSize',16);
                                ylabel('Probability'); xlabel('Mean Entropy');

                                % Bar Plots
                                kax = subplot(2, 5, [1 6]); hold on;    % subplot(numel(h{e,s,t,c})*2, 5, [1 6]); hold on;
                                ind(:,1) = mships < -zthresh(z);	% select node weights which surpass threshold
                                ind(:,2) = mships > zthresh(z); 	% select node weights which surpass threshold
                                ind(:,3) = mships > -zthresh(z) & mships < 0;
                                ind(:,4) = mships < zthresh(z) & mships > 0;
                                a = mships; a(~ind(:,1)) = 0; barh(1:N{e,s}.ROI, a, 'b');
                                a = mships; a(~ind(:,2)) = 0; barh(1:N{e,s}.ROI, a, 'r');
                                a = mships; a(~ind(:,3)) = 0; barh(1:N{e,s}.ROI, a, 'b', 'FaceAlpha',0.3);
                                a = mships; a(~ind(:,4)) = 0; barh(1:N{e,s}.ROI, a, 'r', 'FaceAlpha',0.3);
                                yticks(1:N{e,s}.ROI); set(kax, 'YTickLabel',labels_ROI, 'FontSize', 6);
                                title("Component Membership", 'FontSize',16, 'Position',[0 93]);
                                subtitle(['z-score threshold: ' num2str(zthresh(z))], 'FontSize',12, 'Position',[0 91]);
                                xlabel('z-score', 'FontSize',12);

                                % Brain Renderings
                                kax = subplot(2, 5, [2 3 7 8]); hold on;  % subplot(numel(h{e,s,t,c})*2, 5, [2 3 7 8]); hold on;
                                plot_nodes_in_cortex(cortex, mships, coords_ROI, origin, sphereScale, zthresh(z), [], cind, [], [], rdux);
                                % Note: if want to weight node color by
                                % "strength of association, must encode weighting in cind.node

                                % Save as png file
                                saveas(K(kFig-1), fullfile(path{4}, strcat('entroCompare', fname, '_Z', strjoin(split(num2str(zthresh(z)), '.')', ''))), 'png');
                            end
                        else
                            K(kFig) = figure('Position', [0 0 1280 1024]); kFig = kFig + 1; hold on;

                            % Connectivity
                            kax = subplot(2, 5, [9 10]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [9 10]); hold on;
                            sgtitle([titles{e}, ' Component ', num2str(k(j))], 'FontSize',18);
                            a = squeeze(memberships{e}(:,k(j)))*squeeze(memberships{e}(:,k(j)))';
                            imagesc(a); colorbar; hold on;
                            xlim([1 size(a,2)]); ylim([1 size(a,1)]);
                            yticks(1:N{e,s}.ROI); set(kax, 'YTickLabel',labels_ROI, 'FontSize', 3); xticks([]);
                            title("Connectivity", 'FontSize',16); pbaspect([1 1 1]);

                            % Histogram of component entropies
                            kax = subplot(2, 5, [4 5]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [4 5]); hold on;
                            histogram(entro{e,s}(k(j), :, C{e}(c,1)), 'Normalization','pdf');
                            histogram(entro{e,s}(k(j), :, C{e}(c,2)), 'Normalization','pdf');
                            legend(labels(C{e}(c,:)));
                            title("Entropy", 'FontSize',16);
                            ylabel('Probability'); xlabel('Mean Entropy');

                            % Membership Bar Plots
                            kax = subplot(2, 5, [1 6]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [1 6]); hold on;
                            if sum(sign(squeeze(memberships{e}(:,k(j))))) >= 0
                                ind(:,1) = (mships < 0);
                                ind(:,2) = (mships > 0);
                                a = mships; a(ind(:,1)) = 0; barh(1:N{e,s}.ROI, a);
                                a = mships; a(ind(:,2)) = 0; barh(1:N{e,s}.ROI, a, 'r');
                            elseif sum(squeeze(memberships{e}(:,k(j)))) < 0
                                ind(:,1) = (mships > 0);
                                ind(:,2) = (mships < 0);
                                a = mships; a(ind(:,1)) = 0; barh(1:N{e,s}.ROI, a);
                                a = mships; a(ind(:,2)) = 0; barh(1:N{e,s}.ROI, a, 'r');
                            end
                            yticks(1:N{e,s}.ROI); set(kax, 'YTickLabel',labels_ROI, 'FontSize', 6);
                            title("Component Membership", 'FontSize',16, 'Position',[0 93]);
                            subtitle(['z-score threshold: ' num2str(zthresh(z))], 'FontSize',12, 'Position',[0 91]);
                            xlabel('z-score', 'FontSize',12);

                            % Brain Renderings
                            kax = subplot(2, 5, [2 3 7 8]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [2 3 7 8]); hold on;
                            plot_nodes_in_cortex(cortex, mships, coords_ROI, origin, sphereScale, [], rdux);

                            % Save as png file
                            saveas(K(kFig-1), strcat('entroCompare', fname, '_Z', join(split(num2str(zthresh(z)), '.'))), 'png');
                        end
                    else
                        % Histogram of component entropies
                        K(kFig) = figure('Position', [0 0 1280 1024]); kFig = kFig + 1; hold on;
                        for q = 1:size(C{e},2)
                            histogram(entro{e,s}(k(j), :, C{e}(c,q)), 'Normalization','pdf');
                        end
                        legend(labels(C{e}(c,:)));
                        title(strjoin({'Entropy of', char(labels_ROI(k(j))), 'in', spaces{s}, 'Space'}));
                        ylabel('Counts'); xlabel('Entropy');

                        % Save as png file
                        saveas(K(kFig-1), fullfile(path{4}, strcat('entroCompare', fname, '_Component', num2str(k(j)))), 'png');
                    end
                end
            end
        end
	end
end
clear e s t q B n ax n k j a hg f mships m ind

% Add kFigs to F
if N{1,1}.comp > 1
    close(F(nFig-2:nFig-1));
    nFig = nFig-2;
end
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
% 				C{t} = nchoosek(1:sum(N{t,s}.subjects), 2);
% 				esubdist{t,s,d,2} = zeros(sum(N{t,s}.subjects), sum(N{t,s}.subjects));
% 				esubdist{t,s,d,3} = zeros(sum(N{t,s}.subjects), sum(N{t,s}.subjects));
% 				esubdist{t,s,d,4} = zeros(sum(N{t,s}.subjects), sum(N{t,s}.subjects));
% 			else
% 				C{t} = nchoosek(1:N{t,s}.comp, 2);
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
% 			for i = 1:size(C{t},1)
% 				[esubdist{t,s,d,3}(C{t}(i,2), C{t}(i,1)), ~, esubdist{t,s,d,2}(C{t}(i,1), C{t}(i,2))] = kstest2(e(C{t}(i,1),:), e(C{t}(i,2),:));
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
% 			for i = 1:size(C{t},1)
% 				esubdist{t,s,d,4}(C{t}(i,2), C{t}(i,1)) = abs(em(C{t}(i,1)) - em(C{t}(i,2)));
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
% clear t s d i e em j c


%% Save results

% savefig(F, fullfile(path{4}, strcat('entroCompare', fname)), 'compact');
% clear F nFig
% save(strcat('entroCompare', fname));