function [active] = activeMat(activations, threshold)
% activeMat calculates the activations in the matrix using the method of
% Tagliazucchi et. al. 2016
%	INPUTS:
%		activations: the matrix of continuous activity data
%		threshold: the z-score (number of SD's) threshold of "activity"
%	OUTPUTS:
%		active: activation matrix.  Rows correspond to assemblies, columns
%			correspond to time samples.

% Preallocate event matrices
active = zeros(size(activations));

% Compute z-score
score = zscore(activations')';

% Fill event matrix
active(abs(score) > threshold) = 1;
active = logical(active);

end

