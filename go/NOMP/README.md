# README #

### What is this repository for? ###

* MATLAB code to estimate gains and frequencies in a mixture of sinusoids
* Frequency estimation on the continuum. The frequencies can take any value on [0, 2π)
* Both compressive and regular measurements are supported. So are spectral weighting functions like hamming, hanning, etc.,

### CONTENTS ###

* extractSpectrum.m - frequency estimation algorithm. Takes as input the measurements, the measurement matrix and parameter tau (larger values of tau return SPARSER solutions)
* generateMeasMat.m - generates commonly used measurement matrices including compressive measurement matrices and those with weighing functions
* example.m - illustrates how to use extractSpectrum.m for line spectral estimation
* example_CS.m - same as above for compressive line spectral estimation

### ALGORITHM OVERVIEW ###

* SCENARIO: sinusoid( ω )  = exp(1j * ω * (0:(N-1))')/sqrt(N) with frequency ω ∈ [0, 2π) and length N, 
	S is a M times N measurement matrix and noise is complex Gaussian with covariance matrix sigma2 * eye(M) (white noise)
	
* MEASUREMENTS: Mixture of sinusoids: y = S * ( gain(1) * sinusoid(omega(1)) + ... + gain(K) * sinusoid(omega(K)) ) + noise

* USAGE: [gainList, omegaList, residueList] = extractSpectrum(y, S, tau)
	
* COMMENTS: When tau is large the explanation will be sparse (number of sinusoids estimated will be small). Conversely, if tau is small, the explanation will include many sinusoids (number of sinusoids estimated will be large)

* http://dx.doi.org/10.1109/TSP.2016.2580523, http://dx.doi.org/10.1109/TSP.2014.2306180

### Dowload this repository ###
https://bitbucket.org/wcslspectralestimation/continuous-frequency-estimation/downloads
