function [iFC] = LEdFC(TS, varargin)
%% LEDFC computes the dFC from the time series TS.
%	The precise

%% Set type of dFC to compute

% Set default values
distType = 'cosine';
compressType = 'eigenvector';
N = size(TS, 1);
T = size(TS, 2);

% unpack varagin
for k = 1:2:length(varargin)
	switch varargin{k}
		case 'distType'
			distType = varargin{k+1};
		case 'compressType'
			compressType = varargin{k+1};
		case 'N'
			N = varargin{k+1};
		case 'T'
			T = varargin{k+1};
	end
end

%% Compute and compress dFC

% Compute dFC
iPH = dFC(TS, 'nROI',N, 'T',T, 'distType',distType);

% Compress dFC
switch compressType
	case {'LE', 'LEICA', 'eigenvector'}
		iFC = nan(T, N);
		for t = 1:T
			[iFC(t,:), ~] = eigs(squeeze(iPH(:,:,t)),1);
		end
	case {'average', 'mean'}
		iFC = nan(T, N);
		for t = 1:T
			iFC(t,:) = mean(squeeze(iPH(:,:,t)));
		end
    case 'none'
		iFC = nan(T, length(find(tril(ones(N),-1))));
        for t = 1:T
			d = squeeze(iPH(:,:,t));
			iFC(t,:) = d(find(tril(ones(N),-1))');
        end
    otherwise
        error("Please select one of the supported compression methods.");
end

% Convert dFC matrix to standard format
iFC = iFC';

end

