function F = PlotComparison(data, p, ptarget, sigvec, N, condition, comptype, sigtest)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Open new figure
F = figure;

% Determine how many columns are needed
if N.assemblies{:,comptype} > 0
	np = 3;
else
	np = 2;
end

% Find significantly different assemblies
x = find(sigvec{:,sigtest} == 1);
l = mat2cell(x, ones(length(x),1));
for k = 1:length(x)
	l{k} = num2str(l{k});
end
clear k

% Plot entropy per assembly
for c = 1:N.condition
	
	% Confirm that distribution is in proper format
	if size(data{c},2) ~= N.assemblies{:,'All'}
		data{c} = data{c}';
	end
	
	subplot(N.condition, np, (np*c-(np-1))); hold on;
	boxplot(data{c}, 'Notch','on', 'PlotStyle','compact');
	xlabel('Assembly');
	ylabel(comptype)
	title(strcat(comptype, ' per Assembly: ', condition{c}));
	
	if N.assemblies{:,comptype} > 0
		subplot(N.condition, np, np*c);
		boxplot(data{c}(:,x), 'Notch','on', 'PlotStyle','compact');
		xticklabels(l);
		xlabel('Assembly');
		ylabel(comptype);
		title(strcat(comptype, ' per Significant Assembly: ', condition{c}));
	end
end

if N.assemblies{:,comptype} > 0
	subplot(N.condition, np, [2,5]); hold on;
else
	subplot(N.condition, np, [2,4]); hold on;
end
bar(p);
plot(0:length(p)+1, ptarget./length(p).*ones(1,2+length(p)), '--r');
if N.assemblies{:,comptype} > 0
	scatter(x, ones(1,N.assemblies{:,comptype}), 50, 'red','*');
end
xlabel('Assembly');
ylabel('p-value'); ylim([0 1]);
title(strcat('Assembly Significance: ', condition{c}));
legend('p-values',strcat(sigtest, ' threshold'), 'significant assemblies');

end

