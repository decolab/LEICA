function [dFC] = LEdFC(TS, nROI, T, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Preallocate storage arrays
dFC = nan(T, nROI);

% Indices for vectorizing lower triangle of FC matrix
Isubdiag = find(tril(ones(N.ROI),-1));

% Set default values
distType = 'cosine';
compressType = 'eigenvector';

% unpack varagin
for k = 1:2:length(varargin)
	switch varargin{k}
		case 'phaseType'
			distType = varargin{k+1};
		case 'compressType'
			compressType = varargin{k+1};
	end
end

% Compute dFC
for t = 1:T
	iPH = zeros(nROI, nROI);
	switch distType
		case 'cosine'
			for n = 1:nROI
				for m = 1:nROI
					iPH(n,m) = cos(TS(n,t) - TS(m,t));
				end
			end
		case 'exponential'
			for n = 1:nROI
				for m = 1:nROI
					iPH(n,m) = exp(-3*adif(TS(n,t), TS(m,t))); 
				end
			end
	end

	% Compress dFC
	switch compressType
		case 'eigenvector'
			[dFC(t,:), ~] = eigs(iPH,1);
		case 'average'
			dFC(t,:) = mean(iPH);
		case 'none'
			dFC(t,:) = iPH(Isubdiag);
	end
end

dFC = dFC';

end

