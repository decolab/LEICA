function [BOLD_filt, BOLD_complex, BOLD_phase] = phaseTS(BOLD, N, T, bfilt, afilt)
%% PHASETS computes the Hilbert phase coherence between regions.
%   INPUTS:
%		BOLD: time series of an individual subject.  Rows correspond to
%			regions (ROIs, seeds), columns to time points.
%		nROI: number of regions of interest (ROIs) in BOLD.  Defaults to
%			number of rows in BOLD.
%		T: number of time points in time series.  Defaults to number of
%			columns in BOLD.
%		bfilt: transfer coefficient b of a Butterworth bandpass filter (optional)
%		afilt: transfer coefficient a of a Butterworth bandpass filter (optional)
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
BOLD_filt = zeros(N, T);	% filtered BOLD signal
BOLD_complex = zeros(N, T);	% complex BOLD signal

% Compute filtered BOLD signal & Hilbert phase
if isempty(afilt) && isempty(bfilt)
	for seed = 1:N
		BOLD(seed,:) = detrend(BOLD(seed,:) - mean(BOLD(seed,:)));
		BOLD_complex(seed,:) = hilbert(BOLD(seed,:));
	end
else
	for seed = 1:N
		BOLD(seed,:) = detrend(BOLD(seed,:) - mean(BOLD(seed,:)));
		BOLD_filt(seed,:) = filtfilt(bfilt, afilt, BOLD(seed,:));
		BOLD_complex(seed,:) = hilbert(BOLD_filt(seed,:));
	end
end
BOLD_phase = angle(BOLD_complex);

end

