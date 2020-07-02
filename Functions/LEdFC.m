function [iFC] = LEdFC(TS, varargin)
% LEDFC 
%	


% Set default values
distType = 'cosine';
compressType = 'eigenvector';
nROI = size(TS, 1);
T = size(TS, 2);

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

% Compute dFC
iPH = dFC(TS, 'nROI',nROI, 'T',T, 'distType',distType);

% Compress dFC
switch compressType
	case {'LEICA', 'eigenvector'}
		iFC = nan(T, nROI);
		for t = 1:T
			[iFC(t,:), ~] = eigs(squeeze(iPH(:,:,t)),1);
		end
	case 'average'
		iFC = nan(T, nROI);
		for t = 1:T
			iFC(t,:) = mean(squeeze(iPH(:,:,t)));
		end
	otherwise
		iFC = nan(T, length(find(tril(ones(nROI),-1))));
		for t = 1:T
			d = squeeze(iPH(:,:,t));
			iFC(t,:) = d(find(tril(ones(nROI),-1))');
		end
end

% Convert dFC matrix to standard format
iFC = iFC';

end

