function [syncdata, metastabilitydata] = findStability(phasedata)
% FINDSTABILITY finds the Kuramoto order and metastability of the time
% series phase data.
%	INPUTS
%		phasedata: 
%	OUTPUTS
%		syncdata: 
%		metastabilitydata: 

% Preallocate array
syncdata = nan(1, size(phasedata,2));

% Compute Kuramoto order parameter
for t = 1:size(phasedata,2)		% compute synchrony at each time step
	kudata = sum(complex(cos(phasedata(:,t)), sin(phasedata(:,t))))/size(phasedata,1);	% Kuramoto Order Parameter
	syncdata(t) = abs(kudata);		% the time-wise magnitude of kudata
end

% Compute metastability
metastabilitydata = std(syncdata);	% time-wise variance of the magnitude of kudata

end

