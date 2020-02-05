function [iFC] = LEdFC(TS, varargin)
% LEDFC 
%	


% Set default values
distType = 'cosine';
compressType = 'eigenvector';
nROI = size(TS,1);
T = size(TS,2);

% unpack varagin
for k = 1:2:length(varargin)
	switch varargin{k}
		case 'distType'
			distType = varargin{k+1};
		case 'compressType'
			compressType = varargin{k+1};
		case 'nROI'
			nROI = varargin{k+1};
		case 'T'
			T = varargin{k+1};
	end
end

% Preallocate storage arrays
iFC = nan(T, nROI);

% Indices for vectorizing lower triangle of FC matrix
Isubdiag = find(tril(ones(nROI),-1));

% Compute dFC
iPH = dFC(TS, 'nROI',nROI, 'T',T, 'distType',distType);

% Compress dFC
for t = 1:T
	switch compressType
		case 'eigenvector'
			[iFC(t,:), ~] = eigs(squeeze(iPH(:,:,t)),1);
		case 'average'
			iFC(t,:) = mean(squeeze(iPH(:,:,t)));
		case 'none'
			d = squeeze(iPH(:,:,t));
			iFC(t,:) = d(Isubdiag);
	end
end

% Convert dFC matrix to standard format
iFC = iFC';

end

