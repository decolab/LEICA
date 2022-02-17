function [F, h, p, tstat, FDR, Sidak] = compareComponents(path, z, cind, spaces, dim, ttypes, pTarget, prange, origin,cortex,sphereScale,rdux, ROI,entro,N,memberships,labels,comps)
%  Summary of this function goes here
%	

% Reindex count array
nc = N.IC;
N.IC = cell(1, numel(spaces));
if strcmpi(spaces, 'dFC'); N.IC{strcmpi(spaces, 'dFC')} = N.ROI; end
if strcmpi(spaces, 'IC'); N.IC{strcmpi(spaces, 'IC')} = nc; end; clear nc
nFig = 1;

% file containing cortical atlas
cortex.file = fullfile(path{7}, cortex.file);


%% Compare mean entropy across methods

% Preallocate arrays
mEntro = cell(numel(spaces), numel(dim));
cond = cell(1, N.conditions);

% Set comparison order
for s = 1:numel(spaces)
    for d = 1:numel(dim)

        % Compute means
        mEntro{s,d} = squeeze(mean(entro.(spaces(s)), d, 'omitnan'));
        mEntro{s,d} = array2table(mEntro{s,d}, 'VariableNames', labels);

        % Pairwise comparisons between conditions
        for c = 1:N.comp
            for p = 1:size(comps,2)
                cond{p} = mEntro{s,d}{:,labels{comps(c,p)}}(isfinite(mEntro{s,d}{:,labels{comps(c,p)}}));
            end
            [sig.av.h(s,d,c,1), sig.av.p(s,d,c,1), sig.av.effsize(s,d,c,1)] = kstest2(cond{1}, cond{2});
            [sig.av.p(s,d,c,2), ~, sig.av.effsize(s,d,c,2)] = permutationTest(cond{1}, cond{2}, 10000, 'sidedness','both');
            if sig.av.p(s,d,c,2) < 0.05
                sig.av.h(s,d,c,2) = 1;
            else
                sig.av.h(s,d,c,2) = 0;
            end
        end

        % Run mutliple-comparison corrections for conditions
        [sig.av.FDR(s,d,:,1), sig.av.Bonferroni(s,d,:,1), sig.av.Sidak(s,d,:,1)] = mCompCorr([], squeeze(sig.av.p(s,d,:,1)), pTarget);
        [sig.av.FDR(s,d,:,2), sig.av.Bonferroni(s,d,:,2), sig.av.Sidak(s,d,:,2)] = mCompCorr([], squeeze(sig.av.p(s,d,:,2)), pTarget);
    end
end
clear s d cond

% Split into comparisons
sig.mSubj.h = squeeze(sig.av.h(:,1,:,:));
sig.mSubj.p = squeeze(sig.av.p(:,1,:,:));
sig.mSubj.effsize = squeeze(sig.av.effsize(:,1,:,:));
sig.mIC.h = squeeze(sig.av.h(:,2,:,:));
sig.mIC.p = squeeze(sig.av.p(:,2,:,:));
sig.mIC.effsize = squeeze(sig.av.effsize(:,2,:,:));
if c > 1
    sig.mSubj.FDR = squeeze(sig.av.FDR(:,1,:,:));
    sig.mSubj.Bonferroni = squeeze(sig.av.Bonferroni(:,1,:,:));
    sig.mSubj.Sidak = squeeze(sig.av.Sidak(:,1,:,:));
    sig.mIC.FDR = squeeze(sig.av.FDR(:,2,:,:));
    sig.mIC.Bonferroni = squeeze(sig.av.Bonferroni(:,2,:,:));
    sig.mIC.Sidak = squeeze(sig.av.Sidak(:,2,:,:));
end
clear c

% Display mean entropy histograms
for d = 1:numel(dim)			% dimension average to display
    F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; hold on;
    col = 0;
	for s = 1:numel(spaces)		% space to display (dFC or IC)
		col = col+1;

        % Plot mean entropy histograms
        sz = binWidth(mEntro{s,d}{:,:}, 2);                 % Get bin sizes
        for c = 1:N.comp
            row = numel(spaces)*(c-1);
            subplot(N.comp, numel(spaces), col+row); hold on;
            for f = comps(c,:)
                l = mEntro{s,d}{:,labels{f}};
                l = l(isfinite(l));
                h{c,f} = histogram(l, 'BinWidth',sz, 'Normalization','probability', 'FaceAlpha',0.4, 'FaceColor',cind.hist(f,:));
            end
            means(comps(c,:), h(c,:), sig.av.h(s,d,c,:));     % Plot means, distances for significant differences
            title(strjoin(["Mean Entropy per", dim(d), "in" spaces(s), "Space"]));
            subtitle(strjoin([labels(comps(c,1)) "v." labels(comps(c,2))]));
            xlabel('Mean Entropy'); ylabel('Probability');
        end
	end
end

% Display mean entropy boxplots
F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; hold on;
col = 0;
for d = 1:numel(dim)			% dimension average to display
	for s = 1:numel(spaces)		% space to display (dFC or IC)
        col = col+1;
        b(d,s) = subplot(numel(dim), numel(spaces), col);
        a = table2array(mEntro{s,d});
        l = reshape(repmat(string(labels), [size(a,1) 1]), [numel(a) 1]);
        a = reshape(a, [numel(a) 1]);
        i = isfinite(a);
        a = a(i); l = l(i);
        boxplot(a, l, 'Colors',cind.hist, 'Notch','on'); hold on
        title(strjoin(["Mean Entropy per", dim(d), "in" spaces(s), "Space"]));
        xt = get(b(d,s), 'XTick');
        yt = get(b(d,s), 'YTick');
        yl = ylim(b(d,s));
        amp = [1.02 1.03 1.07];
        
        for c = 1:size(comps,1)
            if sig.av.FDR(s,d,c,:)
                plot(xt(comps(c,:)), [1 1]*max(yt)*amp(1), '-k',  mean(xt(comps(c,:))), max(yt)*amp(2), '*k');
                l = [yl(1) yl(2)*amp(3)];
                ylim(b(d,s), l);
                amp = amp+0.03;
            end
        end
	end
end
clear d s a nc col row mc mp i h hg sz l c f p b yt xt yl amp


%% Compare component-level entropy between conditions

% Preallocate storage arrays
k = cell(numel(spaces), numel(ttypes), N.comp);
r = cell(numel(spaces), numel(ttypes), N.comp);
FDR = cell(1, numel(spaces));
Sidak = cell(1, numel(spaces));
p = cell(1, numel(spaces));
h = cell(1, numel(spaces));
tstat = cell(1, numel(spaces));

% Set index locations
ticlocs = (N.comp/5):(N.comp/5):N.comp;

% Extract thresholds at which components become FDR-significant
for s = 1:numel(spaces)
    % Allocate threshold storage
    FDR{s} = nan(N.IC{s}, numel(ttypes), N.comp, numel(prange));
    Sidak{s} = nan(N.IC{s}, numel(ttypes), N.comp, numel(prange));
    p{s} = nan(N.IC{s}, numel(ttypes), N.comp);
    tstat{s} = nan(N.IC{s}, numel(ttypes), N.comp);
    h{s} = nan(N.IC{s}, numel(ttypes), N.comp);

    % Test for differences
    for t = 1:numel(ttypes)
        disp(strcat("Running ", ttypes(t), " tests on entropy in ", spaces(s), " space."));
        for c = 1:N.comp
            % Run comparisons
            disp(strjoin(['Comparing', labels(comps(c,1)), 'and', labels(comps(c,2))], ' '));
            [h{s}(:,t,c), p{s}(:,t,c), tstat{s}(:,t,c)] = robustTests(squeeze(entro.(spaces(s))(:,1:N.subjects(comps(c,1)),comps(c,1))), squeeze(entro.(spaces(s))(:,1:N.subjects(comps(c,2)),comps(c,2))), [], 'p',prange, 'testtype',ttypes(t));
        end

        % Multiple comparison correction
        if N.comp > 1
            p_dum = reshape(squeeze(p{s}(:,t,:)), [N.comp*N.IC{s}, 1]);	% reshape p-value array
            [f, ~, S] = mCompCorr([], p_dum, prange);                   % run mutliple comparison tests
            FDR{s}(:,t,:,:) = reshape(f, [N.IC{s}, N.comp, numel(prange)]);
            Sidak{s}(:,t,:,:) = reshape(S, [N.IC{s}, N.comp, numel(prange)]);

            % Find components with 5% significance
            for c = 1:N.comp
                [r{s,t,c},~] = unique(find(squeeze(FDR{s}(:,t,c,prange == pTarget))));
                [k{s,t,c},~] = unique(find(squeeze(Sidak{s}(:,t,c,prange == pTarget))));
                k{s,t,c} = union(r{s,t,c}, k{s,t,c});
                disp([num2str(numel(k{s,t,c})), ' component(s) displays significant differences.']); 
            end
        end
    end
end
clear c f S p_dum

% Save significant component indices to files
vn = cell(size(comps,1), 1);
for c = 1:size(comps,1)
    vn{c} = [labels{comps(c,1)}, ' v. ', labels{comps(c,2)}];
end
for s = 1:numel(spaces)
    h{s} = cell2table(squeeze(k(s,:,:)), 'RowNames',ttypes, 'VariableNames',vn);
end
% save(fullfile(path{4}, fname), 'labels', 'h','p','tstat', '-append');
% h = k;
clear d e s t k r c vn


%% Plot significance results

kFig = 1;
for c = 1:N.comp
    F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(4) = 1;
    F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(3) = 1;
    F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(2) = 1;
    F(nFig) = figure('Position', [0 0 1280 1024]); nFig = nFig + 1; n(1) = 1;
    for s = 1:numel(spaces)

        % Set index for all components found across all tests
        k = table;

        % Plot significance test results
        for t = 1:numel(ttypes)

            % Set index for all components found across all tests
            k = union(k, h{s,t,c});

            if N.comp > 1
                % Plot FDR significance level
                figure(F(nFig-2));
                sgtitle([labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);	% sgtitle('FDR Significance Threshold');
                ax(4) = subplot(numel(ttypes), numel(spaces), n(4)); 
                n(4) = n(4)+1;
                imagesc(ax(4), squeeze(FDR{s}(:,t,c,:)), [0 1]); colormap gray; colorbar;
                xticks(ax(4), 1:numel(prange));
                xticklabels(ax(4), strsplit(num2str(prange))); hold on;
                xlabel('FDR Significant'); ylabel('Component');
                title(strjoin([ttypes(t), "test in", spaces(s), "space"]));

                % Plot Sidak significance level
                figure(F(nFig-1));
                sgtitle([labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);	% sgtitle('Sidak Significance Threshold');
                ax(3) = subplot(numel(ttypes), numel(spaces), n(3));
                n(3) = n(3)+1;
                imagesc(ax(3), squeeze(Sidak{s}(:,t,c,:)), [0 1]);
                colormap gray; colorbar;
                xticks(ax(3), 1:numel(prange));
                xticklabels(ax(3), strsplit(num2str(prange))); hold on;
                xlabel('Sidak Significant'); ylabel('Component');
                title(strjoin([ttypes(t), "test in", spaces(s), "space"]));    % , labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);
            end

            % Display test p-values
            figure(F(nFig-4));
            sgtitle([labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);	% sgtitle('p-values');
            ax(2) = subplot(numel(ttypes), numel(spaces), n(2)); 
            n(2) = n(2)+1;
            barh(ax(2), p{s}(:,t,c)); hold on;
            plot((max(p{s}(:,t,c),[],'all')*1.2).*ones(numel(h{s,t,c}),1), h{s,t,c}, '*r');
            title(strjoin([ttypes(t), "test in", spaces(s), "space"]));    % , labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);
            ylabel('Component'); xlabel('p-value');

            % Display test effect sizes
            figure(F(nFig-3));
            sgtitle([labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);	% sgtitle('Effect Sizes');
            ax(1) = subplot(numel(ttypes), numel(spaces), n(1));
            n(1) = n(1)+1;
            barh(ax(1), tstat{s}(:,t,c)); hold on;
            plot(sign(tstat{s}(h{s,t,c},t,c)).*(max(abs(tstat{s}(:,t,c)),[],'all')*1.2).*ones(numel(h{s,t,c}),1), h{s,t,c}, '*r');
            title(strjoin([ttypes(t), "test in", spaces(s), "space"]));    %, labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}]);
            ylabel('Component'); xlabel('Effect Size');
        end

        % Plot entropy distributions for signifcant separations
        if ~isempty(k)
            for j = 1:numel(k)

                % Get bin sizes
                d = squeeze(entro.(spaces(s))(k(j),:,comps));
                sz = binWidth(d, 2);

                % Plot results
                if strncmpi(spaces(s), 'IC', 2)

                    % Scale memberships (optional)
                    if z.scale == true || sum(z.thresh ~= 0, 'all') > 0
                        mships = squeeze(zscore(memberships(:,k(j))));
                    else
                        mships = squeeze(memberships(:,k(j)));
                    end

                    % Enforce consensus: smaller community should be "positive"
                    if sum(sign(mships)) > 0
                        mships = -mships;
                    end

                    % Plot component
                    if sum(z.thresh ~= 0, 'all') > 0
                        for i = 1:numel(z.thresh)
                            K(kFig) = figure('Position', [0 0 1280 1024]); kFig = kFig + 1; hold on;

                            % Connectivity
                            kax = subplot(2, 5, [9 10]); hold on;   % subplot(numel(h{e,s,t,c})*2, 5, [9 10]); hold on;
                            sgtitle(['Component ', num2str(k(j)), ' ', labels{comps(c,1)}, ' vs. ', labels{comps(c,2)}], 'FontSize',18);
                            a = squeeze(memberships(:,k(j)))*squeeze(memberships(:,k(j)))';
                            imagesc(a); colorbar; hold on;
                            xlim([1 size(a,2)]); ylim([1 size(a,1)]);
                            yticks(5:5:N.ROI); xticks(5:5:N.ROI);
                            title("Connectivity", 'FontSize',16); pbaspect([1 1 1]);

                            % Histogram of component entropies
                            kax = subplot(2, 5, [4 5]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [4 5]); hold on;
                            for f = 1:length(comps(c,:))
                                d = entro.(spaces(s))(k(j), :, comps(c,f));
                                d = d(isfinite(d));
                                histogram(d, 'BinWidth',sz, 'Normalization','probability');
                            end
                            legend(labels(comps(c,:)));
                            title("Entropy", 'FontSize',16);
                            ylabel('Probability'); xlabel('Mean Entropy');

                            % Bar Plots
                            kax = subplot(2, 5, [1 6]); hold on;    % subplot(numel(h{e,s,t,c})*2, 5, [1 6]); hold on;
                            ind(:,1) = mships < -z.thresh(i);	% select node weights which surpass threshold
                            ind(:,2) = mships > z.thresh(i); 	% select node weights which surpass threshold
                            ind(:,3) = mships > -z.thresh(i) & mships < 0;
                            ind(:,4) = mships < z.thresh(i) & mships > 0;
                            a = mships; a(~ind(:,1)) = 0; barh(1:N.ROI, a, 'b');
                            a = mships; a(~ind(:,2)) = 0; barh(1:N.ROI, a, 'r');
                            a = mships; a(~ind(:,3)) = 0; barh(1:N.ROI, a, 'b', 'FaceAlpha',0.3);
                            a = mships; a(~ind(:,4)) = 0; barh(1:N.ROI, a, 'r', 'FaceAlpha',0.3);
                            yticks(1:N.ROI); set(kax, 'YTickLabel',ROI{:,"Labels"}, 'FontSize', 6);
                            title("Component Membership", 'FontSize',16);                           % , 'Position',[0 93]);
                            subtitle(['z-score threshold: ' num2str(z.thresh(i))], 'FontSize',12);   % , 'Position',[0 91]);
                            xlabel('z-score', 'FontSize',12);

                            % Brain Renderings
                            kax = subplot(2, 5, [2 3 7 8]); hold on;  % subplot(numel(h{e,s,t,c})*2, 5, [2 3 7 8]); hold on;
                            plot_nodes_in_cortex(cortex, mships, ROI{:,{'x','y','z'}}, origin, sphereScale, z.thresh(i), [], cind, [], [], rdux);
                            % Note: if want to weight node color by strength of association, must encode weighting in cind.node

                            % Save as png file
                            saveas(K(kFig-1), fullfile(path{4}, strjoin([fname, strcat(labels(comps(c,1)),'v',labels(comps(c,2))), strjoin(['Component',num2str(k(j))]), strjoin(["Z", strjoin(split(num2str(z.thresh(i)), '.'),'')], '')], '_')), 'png');
                        end
                    else
                        K(kFig) = figure('Position', [0 0 1280 1024]); kFig = kFig + 1; hold on;

                        % Connectivity
                        kax = subplot(2, 5, [9 10]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [9 10]); hold on;
                        sgtitle(['Component ', num2str(k(j))], 'FontSize',18);
                        a = squeeze(memberships(:,k(j)))*squeeze(memberships(:,k(j)))';
                        imagesc(a); colorbar; hold on;
                        xlim([1 size(a,2)]); ylim([1 size(a,1)]);
                        yticks(5:5:N.ROI); xticks(5:5:N.ROI);
                        title("Connectivity", 'FontSize',16); pbaspect([1 1 1]);

                        % Histogram of component entropies
                        kax = subplot(2, 5, [4 5]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [4 5]); hold on;
                        histogram(entro.(spaces(s))(k(j), :, comps(c,1)), 'BinWidth',sz, 'Normalization','probability');
                        histogram(entro.(spaces(s))(k(j), :, comps(c,2)), 'BinWidth',sz, 'Normalization','probability');
                        legend(labels(comps(c,:)));
                        title("Entropy", 'FontSize',16);
                        ylabel('Probability'); xlabel('Mean Entropy');

                        % Membership Bar Plots
                        kax = subplot(2, 5, [1 6]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [1 6]); hold on;
                        if sum(sign(squeeze(memberships(:,k(j))))) >= 0
                            ind(:,1) = (mships < 0);
                            ind(:,2) = (mships > 0);
                            a = mships; a(ind(:,1)) = 0; barh(1:N.ROI, a);
                            a = mships; a(ind(:,2)) = 0; barh(1:N.ROI, a, 'r');
                        elseif sum(squeeze(memberships(:,k(j)))) < 0
                            ind(:,1) = (mships > 0);
                            ind(:,2) = (mships < 0);
                            a = mships; a(ind(:,1)) = 0; barh(1:N.ROI, a);
                            a = mships; a(ind(:,2)) = 0; barh(1:N.ROI, a, 'r');
                        end
                        yticks(1:N.ROI); set(kax, 'YTickLabel',ROI{:,"Labels"}, 'FontSize', 6);
                        title("Component Membership", 'FontSize',16, 'Position',[0 93]);
                        subtitle(['z-score threshold: ' num2str(z.thresh)], 'FontSize',12, 'Position',[0 91]);
                        xlabel('z-score', 'FontSize',12);

                        % Brain Renderings
                        kax = subplot(2, 5, [2 3 7 8]); hold on;	% subplot(numel(h{e,s,t,c})*2, 5, [2 3 7 8]); hold on;
                        plot_nodes_in_cortex(cortex, mships, ROI{:,{'x','y','z'}}, origin, sphereScale, z.thresh, [], cind, [], [], rdux);
                        
                        % Save as png file
                        saveas(K(kFig-1), fullfile(path{4}, strjoin([fname, strcat(labels(comps(c,1)),'v',labels(comps(c,2))), strjoin(['Component',num2str(k(j))]), strjoin(['Z', split(num2str(z.thresh(i)), '.')], '')], '_')), 'png');
                    end
                else
                    % Histogram of component entropies
                    K(kFig) = figure('Position', [0 0 1280 1024]); kFig = kFig + 1; hold on;
                    for q = 1:size(comps,2)
                        en = entro.(spaces(s))(k(j), :, comps(c,q));
                        en = en(isfinite(en));
                        histogram(en, 'BinWidth',sz, 'Normalization','probability');
                    end
                    legend(labels(comps(c,:)));
                    title(strjoin(["Entropy of", ROI{:,"Labels"}(k(j)), "in", spaces(s), "Space"]));
                    ylabel('Probability'); xlabel('Entropy');

                    % Save as png file
                    saveas(K(kFig-1), fullfile(path{4}, strcat(fname, '_', strcat(labels(comps(c,1)),'v',labels(comps(c,2))), '_Component', num2str(k(j)))), 'png');
                end
            end
        end
    end
end
clear e s t q B n ax n k j a hg f mships m ind w d sz c p_dum en i

% Add kFigs to F
if N.comp > 1
    close(F(nFig-2:nFig-1));
    nFig = nFig-2;
    F = F(1:nFig-2);
end
for k = 1:kFig-1
	F(nFig) = K(k); nFig = nFig + 1;
end
clear K k kFig kax z



%% Compare distances between entropy distributions




end