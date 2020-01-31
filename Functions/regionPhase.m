function [phasedata, timedata] = regionPhase(data, bfilt, afilt)
% REGIONPHASE computes the Hilbert phase coherence between regions for a
% single subject
%   INPUTS:
%		data: time series of an individual subject.  Rows correspond to
%			regions (ROIs, seeds), columns to time points.
%		bfilt: transfer coefficient b of a Butterworth bandpass filter
%		afilt: transfer coefficient a of a Butterworth bandpass filter
%	OUTPUTS:
%		phasedata: phase angle of Hilbert-transformed data
%		timedata: bandpass-filtered time series


% Declare data arrays
timedata = zeros(size(data));
phasedata = zeros(size(data));

% Run filtration
for seed = 1:size(data,1)							% loop over regions (seeds)
	x = data(seed,:) - mean(data(seed,:));			% demean data
	x = detrend(x);									% detrend data
	timedata(seed,:) = filtfilt(bfilt, afilt, x);	% zero phase filter data
end

% Extract hilbert phase
for seed = 1:size(data,1)
	Xanalytic = hilbert(timedata(seed,:));			% convert to analytic (real + Hilbert transformed) signal
	phasedata(seed,:) = angle(Xanalytic);			% compute phase angle of signal
end

end

