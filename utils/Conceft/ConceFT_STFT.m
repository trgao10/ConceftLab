function [tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_STFT(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT, Smooth, Hemi)
%
% Usage: 
% 	[tfrsq, ConceFT, tfrsqtic] = sqSTFT(t, x, lowFreq, highFreq, alpha, WinLen, dim, supp, MT)
%
% MT = 1: ordinary SST; MT > 1: ConceFT
% alpha: resolution in the frequency axis
% WinLen, dim, supp: for hermf.m
%
% Example:
% 	[tfrsq, ConceFT, tfrsqtic] = sqSTFT([1:length(y)]', y, 0,0.5, 0.0002, 121, 4, 6, 10);
%

% N = length(x);

%%%% Multitapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate the window for short time Fourier transform (STFT)
[h, Dh, ~] = hermf(WinLen, dim, supp) ;

%=======================================

fprintf(['Run ordinary STFT-SST (Smooth = ',num2str(Smooth),', Hemi = ',num2str(0),')\n']) ;
[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', Smooth, 0);

%=======================================
ConceFTall = zeros(size(tfrsq,1), size(tfrsq,2), MT);

if MT == 1
    ConceFT = ConceFTall;
else
	%%%% Conceft
    parfor ii = 1:MT
		rv = randn(1, dim); rv = rv ./ norm(rv);
		rh = rv * h; 
		rDh = rv * Dh;

		[~, ~, ConceFTall(:,:,ii)] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, rh', rDh', Smooth, Hemi);
    end

    ConceFT = mean(ConceFTall,3);
end

end
