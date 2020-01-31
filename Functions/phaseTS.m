function [Phase_BOLD] = phaseTS(BOLD, nROI, T, bfilt, afilt)
%% PHASETS computes the phase time series of BOLD connectivity over time
%	

% Initialize storage arrays
Phase_BOLD = zeros(nROI, T);	% phase synchrony array

% Compute filtered BOLD signal & Hilbert phase
for seed = 1:nROI
	BOLD(seed,:) = BOLD(seed,:) - mean(BOLD(seed,:));
	signal_filt = filtfilt(bfilt, afilt, BOLD(seed,:));
	Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
end

end

