function [entro, meanvalue, distribution] = HShannon_segmented(data, TI, nSubjects, nAssembly, co)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Preallocate arrays
entro = nan(nSubjects, nAssembly);
meanvalue = nan(nSubjects, nAssembly);
distribution = nan(nAssembly, size(TI,2));

% Compute entropy per condition, subject, assembly
for ass = 1:nAssembly
	disp(['Analyzing assembly ', num2str(ass)]);
	lista = [];

	for s = 1:nSubjects
		disp(['Analyzing subject ', num2str(s)]);

		% Select the time points representing this subject and task
		T = (TI(1,:)==s);	% & (TI(2,:)==c));

		lista = horzcat(lista, squeeze(data(ass,T)));

		entro(s,ass) = HShannon_kNN_k_estimation(squeeze(data(ass,T)),co);
		
		meanvalue(s,ass) = mean(squeeze(data(ass,T)));
	end

	% Store assembly activations
	distribution(ass,:) = lista;
end

end

