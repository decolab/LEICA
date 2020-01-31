function [phBOLD, timeserietotal, iPH, iPHmean] = phasesync(BOLD, N, T, bfilt, afilt)
%% PHASESYNC computes the phase synchrony of BOLD time series over time
%	

% Compute phase of BOLD time series
phBOLD = phaseTS(BOLD, N, T, bfilt, afilt);

% Compute leading eigenvector of time series
[timeserietotal, iPH, iPHmean] = LEIG(phBOLD, T, N);

end

