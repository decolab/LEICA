function [lifetimes] = lifetime(activations, events)
% LIFETIME extracts the activity lifetimes from an activation matrix.  The
% activation matrix is assumed to be logical.
%	INPUTS:
%		activations: an N x T logical matrix marking the temporal epochs of
%			activity.  N is the number of nodes and T the number of samples
%			(time points).
%		events: an N x T matrix marking the start of activation epochs.  N
%			is the number of nodes and T the number of samples.
%	OUTPUTS:
%		lifetimes: an N x M array of activity lifetimes.  In this instance,
%			N is the number of nodes, whereas M is the (maximum) number of
%			events per node.  Each 

% Find sizes
N = size(activations,1);
T = size(activations,2);
M = sum(events,2);

% Preallocate storage arrays
lifetimes = nan(max(M), N);

% Count length of each activation
n = 1;
[row, col] = find(events');
activations = activations';
for k = 1:length(col)
	if k == length(col)
		lifetimes(n, col(k)) = sum(activations(row(k):T, col(k)));
	elseif col(k+1) == col(k)
		lifetimes(n, col(k)) = sum(activations(row(k):row(k+1)-1, col(k)));
		n = n+1;
	else
		lifetimes(n, col(k)) = sum(activations(row(k):T, col(k)));
		n = 1;
	end
end

% Correct orientation
lifetimes = lifetimes';


end

