function [events] = eventMat(activations)
% activeMat calculates the activations in the matrix using the method of
% Tagliazucchi et. al. 2016
%	INPUTS:
%		activations: the activation time series matrix (output of
%			activeMat).  Strictly a binary matrix (active/inactive).
%	OUTPUTS:
%		events: event matrix.  Rows correspond to assemblies, columns
%			correspond to time samples.

% Fill event matrix
events2 = [zeros(size(activations,1),1) activations(:,1:end-1)];
events = (activations - events2) > 0;

% Convert event matrix to logical
events = logical(events);

end

