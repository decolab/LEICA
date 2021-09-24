function [cohesiveness] = cohesive(memberships, activeprob)
% COHESIVE calculates the heirarchy in an activity pattern
%	INPUTS:
%		memberships: brain assemblies (spatial patterns) detected by
%			extraction algorithms.
%		activeprob: probability of activation of each assembly.  May be
%           vector or matrix.
%	OUTPUTS:
%		cohesiveness: cohesiveness found in the assembly patterns.  Same
%           number of samples (columns) as activeprob
%
% Reference:
%	Gustavo Deco, Josephine Cruzat, & Morten L. Kringelbach.  Brain Songs
% Framework Used for Discovering the Relevant Timescale of the Human Brain.
% Nature Communications 1:10 (2019)

% get number of patterns
[nROI, nAssemblies] = size(memberships);

% Declare data arrays
cohesiveness = nan(nROI,1);

if nAssemblies == 0
	cohesiveness = NaN;
else
	for roi = 1:nROI		% index over ROIs
		suma = 0;
		for ass = 1:nAssemblies		% index over number of assemblies
			suma = suma + abs(memberships(roi,ass)).*sum(abs(memberships(:,ass)))*activeprob(ass);
		end
    	cohesiveness(roi) = suma;
	end
end

end

