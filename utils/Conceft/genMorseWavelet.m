function [psihat,psi] = genMorseWavelet(omega,ga,be,k)
% GENMORSEWAVELET Find the time domain, psi, and frequency domain, psihat
% of the k-th order Morse wavelet
%   [psi,psihat] = morsewavelet(omega,be,ga,k);
%   
%   omega --- the "angular frequency", i.e. with 2*pi multiplied
%   k  --- the order of the Morse wavelet, k=0,1,2,...
%   be --- parameter "beta" in Generalized Morse Wavelet, e.g. be = 8
%   gamma --- parameter "gamma" in Generalized Morse Wavelet, e.g. ga = 3
%   
%   Tingran Gao (trgao10@math.duke.edu)
%   last modified: Jan 30, 2017
%

if nargin < 4
    k = 0;
end

fo = morsepeakfreq(ga,be);
basicmorse =2*exp(-be.*log(fo)+fo.^ga+be.*log(abs(omega))-abs(omega).^ga).*(omega>0);
%%% In principle one should normalize each wavelet appropriately by
%%% correctly computing coefficients Akbg; for SST this really only causes
%%% serious numerical issues. We fix Akbg = 1 just in this routine.
Akbg = 1;
% Akbg = morsenormconstant(be,ga,k);

if k == 0
    psihat = Akbg*basicmorse;
else
    %%% when k>0, need to compute Generalized Laguerre Polynomials
    arr = 2*((abs(omega)).^ga);
    psihat = Akbg*basicmorse.*laggen(arr,k,(2*be+1)/ga-1);
end

%%% Only invert if the second output is requested.
%%% Note that J. Lilly's code actually performs an additional rotation
%%% before takign ifft --- should check if ifftshift does the same job.
if nargout == 2
    psi = ifftshift(ifft(psihat));
end

end

% %-------------------------------------------------------------
% function Akbg = morsenormconstant(be,ga,k)
% % Returns the Morse wavelet normalization constant for the k-th order
% % Morse wavelet based on the parameters, \beta and \gamma
% r = (2*be+1)/ga;
% Akbg = double((ga.*(2.^r).*exp(gammaln(k+1)-gammaln(k+r))).^(1/2));
% % Akbg = double((2*pi*ga.*(2.^r).*exp(gammaln(k+1)-gammaln(k+r))).^(1/2));
% 
% end
% %-------------------------------------------------------------------

function Lkc = laggen(x,k,c)
% LAGGEN compute generalized Laguerre poly L_k^c(x)
% 
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Jan 29, 2017
%

Lkc = zeros(size(x));
sn = -1;
for m=0:k
    sn = -sn;
    ga = exp(gammaln(k+c+1)-gammaln(k-m+1)-gammaln(c+m+1)-gammaln(m+1));
    Lkc = Lkc+sn.*ga.*x.^m;
end

end

