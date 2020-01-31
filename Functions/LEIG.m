function [timeserietotal, iPH, iPHmean] = LEIG(phBOLD, T, N)
%%	LEIG computes instantaneous phase connectivity & leading eigenvector
% LEIG computes the instantaneous phase connectivity, its leading
% eigenvector, and its time average for one subject.

% Initialize storage array for leading eigenvector
iPH = nan(N, N, T);
timeserietotal = nan(T, N);

% Compute instantaneous phase connectivity & leading eigenvector
for t = 1:T
	
	% Compute instantaneous phase connectivity
	for n = 1:N
		for m = 1:N
			iPH(n,m,t) = cos(phBOLD(n,t) - phBOLD(m,t));
		end
	end
	
	try
		[V1,~] = eigs(squeeze(iPH(:,:,t)), 1);
	catch
		try
			[V1,~] = eigs(squeeze(iPH(:,:,t)), eye(size(squeeze(iPH(:,:,t)))), 1);
		catch
			nA = sum(squeeze(iPH(:,:,t)).^2, 'all') / numel(squeeze(iPH(:,:,t)));
			sA = squeeze(iPH(:,:,t));
			sA(sA.^2 < 1e-10*nA) = 0;
			[V1,~] = eigs(sA, 1);
		end
	end
	
	% Compute & record leading eigenvector
	timeserietotal(t,:) = V1;
end

% find time average of instantaneous phase connectivity
iPHmean = squeeze(mean(iPH,3));

% Convert timeserietotal to standard time series format
timeserietotal = timeserietotal';

end

