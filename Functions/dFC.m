function [iPH] = dFC(TS, varargin)
%% DFC computes the dynamic functional connectivity from the time series.
% Inputs:
%	TS: time series from which to compute dFC.  Must already be processed
%		into appropriate format.
%	distType (optional): similarity measure to use in computing dFC.
%		Supports cosine similarity and a variant of exponential distance.
%		Default: cosine similarity.
%	nROI (optional): number of regions of interest (ROIs) in time series.
%		Must be less than or equal to number of rows in TS.
%		Default: number of rows in TS.
%	T (optional): number of time points in scan.  Must be less than or
%		equal to number of columns in TS.
%		Default: number of columns in TS.


%% Set type of dFC to compute

% Set default values
distType = 'cosine';
nROI = size(TS,1);
T = size(TS,2);

% unpack varagin
for k = 1:2:length(varargin)
	switch varargin{k}
		case 'distType'
			distType = varargin{k+1};
		case 'nROI'
			nROI = varargin{k+1};
		case 'T'
			T = varargin{k+1};
	end
end


%% Compute dFC
%   Types available: cosine distance or exponential

% Preallocate storage
iPH = zeros(nROI, nROI, T);

% Compute dFC
for t = 1:T
	switch distType
		case 'cosine'
			for n = 1:nROI
				for m = 1:nROI
					iPH(n,m,t) = cos(TS(n,t) - TS(m,t));
				end
			end
		case 'exponential'
			for n = 1:nROI
				for m = 1:nROI
					iPH(n,m,t) = exp(-3*min(abs(TS(n,t)-TS(m,t)), 2*pi-abs(TS(n,t)-TS(m,t)))); 
				end
			end
	end
end

end

