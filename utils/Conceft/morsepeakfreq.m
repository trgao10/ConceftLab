function [peakAF, peakCF] = morsepeakfreq(ga,be)
% peakfreq = morsepeakfreq(ga,be) returns the peak frequency for the 
%   zero-th order member of the Morse wavelet family parameterized by ga
%   (gamma) and be (beta): % $(\frac{\beta}{\gamma})^{1/\gamma}$
%   This peak frequency has nothing to do with the order of Morse wavelets.
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Jan 29, 2017
% 

peakAF = exp(1/ga*(log(be)-log(ga)));

%%% obtain the peak frequency in cyclical frequency (mostly useless)
peakCF = peakAF/(2*pi);

end
