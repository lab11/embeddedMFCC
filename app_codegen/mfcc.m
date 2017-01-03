% Samuel Rohrer < rohrer@umich.edu >
% December 7, 2016
% mfcc: compute the MFCC coefficients for a 30 ms window and store all of 
%       these results in matrix where each row is a vector of MFCC coefficients

function [mfcc_coeff] = mfcc( currFrame ) %#codegen

	assert( isa(currFrame, 'double'));

	% framing the signal to prevent spectral leakage
	Fs = 44100;
	samples10ms = Fs / 100;
	% calculate the hamming window
	hamWindow = hamming(3 * samples10ms );
	% apply the hamming window
	currFrame_hamming = currFrame .* hamWindow' ;
	% take the DFT of the current frame
	dft_currFrame = fft(currFrame_hamming);
	
	% put into mel freq banks
	% calculate magnitude of DFT coefficients
	mag_dft_currFrame = abs(dft_currFrame).^2;
	% num freq bins in DFT
	numBins = floor(numel(mag_dft_currFrame)); 
	% remove redundant information from DFT
	mag_dft_currFrame = mag_dft_currFrame(1 : numBins);
	
	% Compute the Mel filters
	maxFreq = 3000; % in hertz
	minFreq = 100; % in hertz
	numMelBins = 18; % somewhat arbitrarily chosen
	melMax = 1125*log(1 + (maxFreq/700)); % compute max mel freq
	melMin = 1125*log(1 + (minFreq/700)); % compute min mel freq
	melBins = linspace(melMin, melMax, numMelBins); % linearly space the mel bins
	% find the mel filters
	melFilters = 700.*(exp(melBins./1125) - 1);
	% these now correspond to the index in the DFT of each bin
	melFilters = floor((numBins).*melFilters / Fs);

	% now the triangular binning
	% iterate over every bin
	mfc_spectrum_bins = [0];
	for i=2:(numel(melFilters)-1)
		left_edge = melFilters(i-1);
		center = melFilters(i);
		right_edge = melFilters(i+1);
		sum = 0;
		% iterate over left side individual bin to the center
		for j=1:(center - left_edge)
			sum = sum + ((j-1)/(center-left_edge))*mag_dft_currFrame(left_edge+(j-1)); 
		end
		% add the middle term
		sum = sum + mag_dft_currFrame(center);
		% iterate over left side individual bin to the center
		diff_r = right_edge - center;
		for j=1:diff_r
			sum = sum + ((diff_r - j)/diff_r)*mag_dft_currFrame(center+(j)); 
		end			
		mfc_spectrum_bins = [mfc_spectrum_bins sum];
	end

	% compute the log-spectral vector
	log_spectral = log(mfc_spectrum_bins);
	% get the mel coefficients by take the IFFT
	%   in this case reduces to DCT as signal real
	mfcc_coeff = fft(real(log_spectral));

end % mfcc
