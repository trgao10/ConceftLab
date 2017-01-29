function [stftResult, stftTic, sqstftResult, sqstftTic] = sqstft(x, lowFreq, highFreq, DF, tDS, Win, dWin)
%
% Synchrosqueezing --- Short-Time Fourier Transform (STFT) version
%
%	x            : signal to be analyzed
%	[lowFreq, highFreq] : frequency range within [0, 0.5]
%	DF           : the resolution in the frequency axis
%	tDS          : the time axis in the final TF representation is downsampled by tDS (set it to 1 unless you know what you are doing)
%	H            : frequency smoothing window, H(0) being forced to 1
%   DH           : differentiation of H	
%	stftResult   : STFT result
%	sqstftResult : synchrosqueezed STFT result
%
%

%%% number of bins along the frequency axis
%%% alwasys assume signal is sampled at period 1 second, so the frequecy is
%%% bounded by [0,1]; for symmetry (and centrality), shift (0.5,1] to (-0.5,0]
NF = length(-0.5+DF:DF:0.5);

[xrow,xcol] = size(x);
t = 1:length(x);
tLen = length(t(1:tDS:length(x)));

%%% since the signal is real, only concerned about the positive frequency
%%% half interval
LIdx = round((NF/2)*(lowFreq/0.5)) + 1;
HIdx = round((NF/2)*(highFreq/0.5));
fLen = HIdx - LIdx + 1;

%====================================================================
%%%%% check input signals
if (xcol~=1)
    error('X must have only one column');
elseif highFreq > 0.5
    error('TopFreq must be a value in [0, 0.5]');
elseif (tDS < 1) || rem(tDS,1)
    error('tDS must be an integer value >= 1');
end

%%% Win is the window function (Win and dWin have the same size)
[hrow,hcol] = size(Win);
hWinLen = (hrow-1)/2; %% half window length
if (hcol~=1) || (rem(hrow,2)==0)
    error('H must be a smoothing window with odd length');
end

%====================================================================
%%%%% run STFT and reassignment rule
stftResult = zeros(NF/2, tLen);
stftTic = linspace(0, 0.5, NF/2)';
sqstftResult = zeros(fLen, tLen);
sqstftTic = linspace(lowFreq, highFreq, fLen)';

Ex = mean(abs(x).^2);
threshold = 1.0e-8*Ex;  % originally it was 1e-6*Ex

for tIdx = 1:tLen
    ti = t((tIdx-1)*tDS+1); 
    tau = -min([round(NF/2)-1,hWinLen,ti-1]):min([round(NF/2)-1,hWinLen,xrow-ti]);
    indices = rem(NF+tau,NF)+1;
    norm_h = norm(Win(hWinLen+1+tau));

    %%%% truncate signal using the window function
	tf0 = zeros(NF,1); tf1 = zeros(NF,1);
    tf0(indices) = x(ti+tau).*conj(Win(hWinLen+1+tau)) / norm_h;
    tf1(indices) = x(ti+tau).*conj(dWin(hWinLen+1+tau)) / norm_h;
    %%% after FFT the frequencies in tf0, tf1 ranges in [0,1], and we only
    %%% pick the first half interval
    %%% tf0, tf1 are ordered from low-frequency to high-frequency
    tf0 = fft(tf0);
    tf0 = tf0(1:NF/2);
    tf1 = fft(tf1);
    tf1 = tf1(1:NF/2);

	%%%% get the first order omega
	omega = zeros(size(tf1)) ;
    avoid_warn = find(abs(tf0) > threshold);
    %%% Instantaneous frequencies are not bounded in [0,1]; instead, it can
    %%% be as large as 1/dt. The omega below is already transformed into
    %%% indices --- all instantaneous frequencies are multiplied by NF and
    %%% rounded-off to index into corresponding frequency bins
    %%% More concretely, this should really be rounding-off the following:
    %%% (NF/2)*[tf1(avoid_warn)./(tf0(avoid_warn)+eps)/(2.0*pi)]/0.5
	omega(avoid_warn) = round(imag(NF*tf1(avoid_warn)./(tf0(avoid_warn)+eps)/(2.0*pi)));

	sst = zeros(fLen,1);
    
    for jcol = 1:NF/2
        if abs(tf0(jcol)) > threshold
            %%% additional shift of phase needed here because the algorithm
            %%% uses modified STFT instead of the standard STFT
            jcolhat = jcol-omega(jcol);
            if (jcolhat <= HIdx) && (jcolhat >= LIdx)
            	sst(jcolhat-LIdx+1) = sst(jcolhat-LIdx+1) + tf0(jcol);
            end
        end
    end

	stftResult(:, tIdx) = tf0(1:NF/2);
	sqstftResult(:, tIdx) = sst * 2 * (stftTic(2)-stftTic(1));
end

