function [Phase_BOLD] = phaseTS(BOLD, N, T, bfilt, afilt)
%% PHASETS computes the Hilbert phase coherence between regions.
%   INPUTS:
%		BOLD: time series of an individual subject.  Rows correspond to
%			regions (ROIs, seeds), columns to time points.
%		nROI: number of regions of interest (ROIs) in BOLD.  Defaults to
%			number of rows in BOLD.
%		T: number of time points in time series.  Defaults to number of
%			columns in BOLD.
%		bfilt: transfer coefficient b of a Butterworth bandpass filter
%		afilt: transfer coefficient a of a Butterworth bandpass filter
%	OUTPUTS:
%		phasedata: phase angle of Hilbert-transformed data

% Set default values
if isempty(N)
	N = size(BOLD,1);
end
if isempty(T)
	T = size(BOLD,2);
end

% Initialize storage arrays
Phase_BOLD = zeros(N, T);	% phase synchrony array

% Compute filtered BOLD signal & Hilbert phase
for seed = 1:N
	BOLD(seed,:) = detrend(BOLD(seed,:) - mean(BOLD(seed,:)));
	signal_filt = filtfilt(bfilt, afilt, BOLD(seed,:));
	Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
end

end

