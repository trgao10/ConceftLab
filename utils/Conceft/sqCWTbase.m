function [tfr, tfrtic, tfrsq, tfrsqtic, sqPsift, omega, sqDPsift, dtfr] = sqCWTbase(t, x, freqlow, freqhigh, alpha, opts, Smooth, Hemi)
%
% Synchrosqueezing transform ver 0.5 (2015-03-09)
% You can find more information in 
%	http://sites.google.com/site/hautiengwu/
%
% Example: [~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);
%	time: 	time of the signal
%	xm: 	the signal to be analyzed
%	[lowfreq, highfreq]: the frequency range in the output time-frequency
%                        representation. For the sake of computational efficiency.
%	alpha:	the frequency resolution in the output time-frequency representation
%	opts:	parameters for the CWT analysis. See below
%	tfr/tfrtic:	the CWT and its SCALE tic
%	tfrsq/tfrsqtic: the SST-CWT and its frequency tic
%
% by Hau-tieng Wu v0.1 2011-06-20 (hauwu@math.princeton.edu)
%		  v0.2 2011-09-10
%		  v0.3 2012-03-03
%		  v0.4 2012-12-12
%		  v0.5 2015-03-09

	%% you can play with these 4 parameters, but the results might not be
	%% that different if the change is not crazy
% Ex = mean(abs(x).^2);
Gamma = 1.0e-8;
% Gamma = 1.0e-8*Ex;  % originally it was 1e-6*Ex
%Gamma = 1e-8 ;
% nvoice = 128;
% scale = 1;
% oct = 1;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% noctave = floor(log2(length(x))) - oct;
dt = t(2) - t(1);

    
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %% Continuous wavelet transform 
% x = x - mean(x);
% [tfr, tfrtic, sqPsift] = CWT(t, x, opts);
[tfr, tfrtic, sqPsift, sqDPsift, dtfr] = CWTDW(x, opts);
% Dtfr = (-1i/2/pi/dt)*[tfr(2:end,:) - tfr(1:end-1,:); tfr(end,:)-tfr(end-1,:)];
% Dtfr((abs(tfr) < Gamma)) = NaN;
% omega = abs(Dtfr./tfr);
dtfr = dtfr/dt;
omega = imag(dtfr./tfr);
omega(abs(omega)<Gamma | isinf(abs(omega))) = NaN;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    	%% Synchro-squeezing transform

[tfrsq, tfrsqtic] = SQ(tfr, omega, freqlow, freqhigh, alpha, dt, Smooth, Hemi);
tfr = transpose(tfr);
tfrsq = transpose(tfrsq);
dtfr = transpose(dtfr);

end

%=====================================================================
	%% function for CWT-based SST
function [tfrsq, tfrsqtic] = SQ(tfd, omega, tfrsqticlow, tfrsqtichigh, alpha, dt, Smooth, Hemi)

nvoice = 32;
scale = 1;
oct = 1;

omega = abs(omega);
[n, nscale] = size(tfd);

% nalpha = floor((tfrsqtichigh - tfrsqticlow)./alpha); %%% very stupid and ad-hoc, making no sense whatsoever
tfrsqticlow = max(tfrsqticlow,(1/n)/dt);
tfrsqtichigh = min(tfrsqtichigh, (1/2)/dt);
% delta_omega = log2((tfrsqtichigh/tfrsqticlow))/nalpha;
% tfrsqtic = 2.^(delta_omega*(1:1:nalpha))*tfrsqticlow;

nalpha = floor((tfrsqtichigh - tfrsqticlow)./alpha);
tfrsqtic = (1:1:nalpha)*alpha + tfrsqticlow;

tfrsq = zeros(n, nalpha);
ntfrsqtic = length(tfrsqtic);
	
Mid = round(length(tfrsqtic)/2) ;
Delta = 20*(tfrsqtic(2)-tfrsqtic(1)).^2 ;
weight = exp(-(tfrsqtic(Mid-10:Mid+10)-tfrsqtic(Mid)).^2/Delta) ;
weight = weight ./ sum(weight) ;
weightIDX = [Mid-10:Mid+10] - Mid ;

for b = 1:n             %% Synchro-
    for kscale = 1: nscale       %% -Squeezing
        
        qscale = scale .* (2^(kscale/nvoice));
        
        if (isfinite(omega(b, kscale)) && (omega(b, kscale)>0))
%             k = floor((log2(omega(b,kscale))-log2(tfrsqticlow))/delta_omega) + 1;
            k = floor( ( omega(b,kscale) - tfrsqticlow )./ alpha )+1;
            
            if (isfinite(k) && (k > 0) && (k < ntfrsqtic-1))
                ha = tfrsqtic(k+1)-tfrsqtic(k);
                
                
                if Smooth
                    IDXb = find((k+weightIDX < ntfrsqtic-1) & (k+weightIDX > 0)) ;
                    IDXa = k+weightIDX(IDXb) ;
                    
                    if Hemi
                        if real(tfd(b,kscale))>0
                            tfrsq(b,IDXa) = tfrsq(b,IDXa) + weight(IDXb)*log(2)*tfd(b,kscale)*sqrt(qscale)./ha/nvoice;
                        else
                            tfrsq(b,IDXa) = tfrsq(b,IDXa) - weight(IDXb)*log(2)*tfd(b,kscale)*sqrt(qscale)./ha/nvoice;
                        end
                    else
                        tfrsq(b,IDXa) = tfrsq(b,IDXa) + log(2)*tfd(b,kscale)*sqrt(qscale)./ha/nvoice;
                    end
                else
                    
                    if Hemi
                        if real(tfd(b,kscale))>0
                            tfrsq(b,k) = tfrsq(b,k) + log(2)*tfd(b,kscale)*sqrt(qscale)./ha/nvoice;
                        else
                            tfrsq(b,k) = tfrsq(b,k) - log(2)*tfd(b,kscale)*sqrt(qscale)./ha/nvoice;
                        end
                    else
                        tfrsq(b,k) = tfrsq(b,k) + log(2)*tfd(b,kscale)*sqrt(qscale)./ha/nvoice;
                    end
                    
                end
                
            end
        end
    end
end
end
