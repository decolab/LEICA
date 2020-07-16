function [phBOLD, iPH] = phasesync(BOLD, N, T, bfilt, afilt, aType)
%% PHASESYNC computes the phase angle and dFC of the BOLD time series.
%	

% Set defaults
if isempty(aType)
	aType = struct;
	aType.dist = 'cosine';
	aType.compress = 'eigenvector';
end

% Compute phase of BOLD time series
phBOLD = phaseTS(BOLD, N, T, bfilt, afilt);

% Compute leading eigenvector of time series
iPH = LEdFC(phBOLD, 'distType', aType.dist, 'compressType', aType.compress, 'N', N, 'T', T);

end

