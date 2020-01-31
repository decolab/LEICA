function [terminus] = terminusMat(activations)
% terminusMat calculates the end of  activation epochs using the method
% of Tagliazucchi et. al. 2016
%	INPUTS:
%		activations: the activation time series matrix (output of
%			activeMat).  Strictly a binary matrix (active/inactive).
%	OUTPUTS:
%		terminus: terminus matrix.  Rows correspond to assemblies, columns
%			correspond to time samples.

% Fill event matrix
events2 = [activations(:,2:end) zeros(size(activations,1),1)];
terminus = (activations-events2) > 0;

% Convert event matrix to logical
terminus = logical(terminus);

end

