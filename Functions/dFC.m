function [iPH] = dFC(TS, varargin)
% DFC 
%	


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
					iPH(n,m,t) = exp(-3*adif(TS(n,t), TS(m,t))); 
				end
			end
	end
end

end

