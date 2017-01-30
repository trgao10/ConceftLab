function [psift,frequencies] = morsewavft(omega,scales,ga,be,k)


psift = zeros(length(scales),length(omega));

for jj = 1:length(scales)
    psift(jj,:) = genMorseWavelet(scales(jj)*omega,ga,be,k);
end

peakradfreq = morsepeakfreq(ga,be);
frequencies = (1/(2*pi))*peakradfreq./scales;

end
