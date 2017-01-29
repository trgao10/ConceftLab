function [cwtResult, cwtTic, sqcwtResult, sqcwtTic] = sqcwt(t, x, freqlow, freqhigh, DF, opts)
%
% Example: [~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);
%	time: 	time of the signal
%	xm: 	the signal to be analyzed
%	[lowfreq, highfreq]: the frequency range in the output time-frequency
%                        representation. For the sake of computational efficiency.
%	alpha:	the frequency resolution in the output time-frequency representation
%	opts:	parameters for the CWT analysis. See below
%	tfr/tfrtic:	the CWT and its scale tic
%	tfrsq/tfrsqtic: the SST-CWT and its frequency tic
%

Ex = mean(abs(x).^2);
Gamma = 1.0e-8*Ex;  % originally it was 1e-6*Ex
dt = t(2) - t(1);

[cwtResult, cwtTic] = CWT(t, x, opts);
Dtfr = (-1i/2/pi/dt)*[cwtResult(2:end,:) - cwtResult(1:end-1,:); cwtResult(end,:)-cwtResult(end-1,:)];

Dtfr((abs(cwtResult) < Gamma)) = NaN;
omega = Dtfr./cwtResult;

[sqcwtResult, sqcwtTic] = SQ(cwtResult, omega, freqlow, freqhigh, DF);
cwtResult = cwtResult';
sqcwtResult = sqcwtResult';

end

%=====================================================================
	%% function for CWT-based SST
function [tfrsq, tfrsqtic] = SQ(cwtResult, omega, lowFreq, highFreq, DF)

nvoice = 32;
scale = 2;

omega = abs(omega);
[n, nscale] = size(cwtResult);

nalpha = floor((highFreq - lowFreq)./DF);
tfrsq = zeros(n, nalpha);
tfrsqtic = (1:1:nalpha)*DF + lowFreq;
ntfrsqtic = length(tfrsqtic);
	
for b = 1:n
    for kscale = 1:nscale
        qscale = scale .* (2^(kscale/nvoice));
        if (isfinite(omega(b, kscale)) && (omega(b, kscale)>0))
            k = floor( ( omega(b,kscale) - lowFreq )./ DF )+1;
            if (isfinite(k) && (k > 0) && (k < ntfrsqtic-1))
                ha = tfrsqtic(k+1)-tfrsqtic(k);
                tfrsq(b,k) = tfrsq(b,k) + log(2)*cwtResult(b,kscale)*sqrt(qscale)./ha/nvoice;
            end
        end
    end
end

end
