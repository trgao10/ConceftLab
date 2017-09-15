function [sst,f,stftcfs,phasetf,stftfreqs,reassntRule,dstftcfs] = fsstgao_bak(x,varargin)
%STFT Synchrosqueezed Transform (adapted from the wsst function shipped
%with the MATLAB Wavelet Toolbox)
%   Tingran Gao (trgao10@math.duke.edu)
%   Feb 15, 2017
%   

narginchk(1,15);
nbSamp = numel(x);
% x = x(:)';
x = x(:);
validateattributes(x,{'double'},{'column','finite','real'},'fsst','X');
if numel(x)<4
    error(message('Wavelet:synchrosqueezed:NumInputSamples'));
end
params = parseinputs(varargin{:});

% %%% If sampling frequency is specified, dt = 1/fs
% if (isempty(params.fs) && isempty(params.Ts))
%     % The default is 1 for normalized frequency
%     dt = params.dt;
%     Units = '';
% elseif (~isempty(params.fs) && isempty(params.Ts))
%     % Accept the sampling frequency in hertz
%     fs = params.fs;
%     dt = 1/fs;
%     Units = '';
% elseif (isempty(params.fs) && ~isempty(params.Ts))
%     % Get the dt and Units from the duration object
%     [dt,Units] = getDurationandUnits(params.Ts);
% else
% end

%%% Construct time series to analyze, pad if necessary
meanSIG = mean(x);
x = x - meanSIG;
NumExten = 0;

if params.pad
    %%% Pad the time series symmetrically
    np2 = nextpow2(nbSamp);
    NumExten = 2^np2-nbSamp;
    x = wextend('1d','symw',x,NumExten,'b');
end

%%% Record data length plus any extension
N = numel(x);

%%% Generate window function
[winfunc,dwinfunc] = sstwinfunc(params.WIN,params.winparam);

%%% Create frequency vector for output
% freq = params.FreqBounds(1):params.FreqRes:params.FreqBounds(2);
% NbFreq = length(freq);
NbFreq = numel(-0.5*params.fs+params.FreqRes:params.FreqRes:0.5*params.fs);

%%% winfunc and dwinfunc should have exactly the same size
[hrow,hcol] = size(winfunc);
hWinLen = (hrow-1)/2; %% half window length
if (hcol~=1) || (rem(hrow,2)==0)
    error('Window must be a smoothing window with odd length');
end

%%% Since the signal is real, only concerned about the positive frequency
%%% half interval
LIdx = round((NbFreq/2)*(params.FreqBounds(1)/0.5/params.fs)) + 1;
HIdx = round((NbFreq/2)*(params.FreqBounds(2)/0.5/params.fs));
fLen = HIdx - LIdx + 1;
tLen = N;

%%% run STFT and reassignment rule
stftcfs = zeros(NbFreq/2, tLen);
dstftcfs = zeros(NbFreq/2, tLen);
phasetf = zeros(NbFreq/2, tLen);
stftfreqs = linspace(0, 0.5*params.fs, NbFreq/2)';
sst = zeros(fLen, tLen);
f = linspace(params.FreqBounds(1), params.FreqBounds(2), fLen)';
reassntRule = nan(fLen, tLen);

% keyboard

for tIdx = 1:tLen
    tau = -min([round(NbFreq/2)-1,hWinLen,tIdx-1]):min([round(NbFreq/2)-1,hWinLen,N-tIdx]);
    indices = rem(NbFreq+tau,NbFreq)+1;
    norm_h = norm(winfunc(hWinLen+1+tau));

    %%%% truncate signal using the window function
	tf0 = zeros(NbFreq,1); tf1 = zeros(NbFreq,1);
    tf0(indices) = x(tIdx+tau).*conj(winfunc(hWinLen+1+tau)) / norm_h;
    tf1(indices) = x(tIdx+tau).*conj(dwinfunc(hWinLen+1+tau)) / norm_h;
    %%% After FFT the frequencies in tf0, tf1 ranges in [0,1], and we only
    %%% pick the first half interval.
    %%% tf0, tf1 are ordered from low-frequency to high-frequency
    tf0 = fft(tf0);
    tf0 = tf0(1:NbFreq/2);
    tf1 = fft(tf1);
    tf1 = tf1(1:NbFreq/2);

	%%%% get the first order omega
	omega = zeros(size(tf1)) ;
    avoid_warn = find(abs(tf0) > params.thr);
    %%% Instantaneous frequencies are not bounded in [0,1]; instead, it can
    %%% be as large as 1/dt. The omega below is already transformed into
    %%% indices --- all instantaneous frequencies are multiplied by NF and
    %%% rounded-off to index into corresponding frequency bins.
    %%% More concretely, this should really be rounding-off the following:
    %%% (NF/2)*[tf1(avoid_warn)./(tf0(avoid_warn)+eps)/(2.0*pi)]/0.5
	omega(avoid_warn) = round(imag(NbFreq*tf1(avoid_warn)./(tf0(avoid_warn)+eps)/(2.0*pi)));

	sstLocal = zeros(fLen,1);
    
    for jcol = 1:NbFreq/2
        if abs(tf0(jcol)) > params.thr
            %%% additional shift of phase needed here because the algorithm
            %%% uses modified STFT instead of the standard STFT
            jcolhat = jcol-omega(jcol);
            if (jcolhat <= HIdx) && (jcolhat >= LIdx)
                %%% elements in "jcol" go to elements in "jcolhat-LIdx+1"
            	sstLocal(jcolhat-LIdx+1) = sstLocal(jcolhat-LIdx+1) + tf0(jcol);
                reassntRule(jcol,tIdx) = jcolhat-LIdx+1;
%                 reassntRule(jcol,tIdx) = jcol-(jcolhat-LIdx+1);
            end
        end
    end

	stftcfs(:, tIdx) = tf0(1:NbFreq/2);
    dstftcfs(:, tIdx) = tf1(1:NbFreq/2);
    phasetf(:, tIdx) = omega(1:NbFreq/2);
	sst(:, tIdx) = sstLocal * 2 * params.FreqRes;
end


% %%% Create frequency vector for STFT computation
% omega = (1:fix(N/2));
% omega = omega*((2.*pi)/N); %%% omage so obtained are angular frequencies (with 2*pi multiplied)
% omega = [0., omega, -omega(fix((N-1)/2):-1:1)];

% %%% Compute FFT of the (padded) time series
% xdft = fft(x);


%%% Though we only need to obtain the plots within params.FreqBounds, the
%%% STFT has to be carried out for the full range of spectrum
% NbFreq = numel(-0.5+params.FreqRes:params.FreqRes:0.5);
% keyboard
% 
% [psift,dpsift,stftfreqs] = sstwaveft(params.WIN,omega,scales,params.winparam);
% keyboard

% %%% Obtain STFT coefficients and derivatives
% stftcfs = ifft(repmat(xdft,NbFreq,1).*psift,[],2);
% dstftcfs = ifft(repmat(xdft,NbFreq,1).*dpsift,[],2);

%%% Remove padding if any
stftcfs = stftcfs(:,NumExten+1:end-NumExten);
dstftcfs = dstftcfs(:,NumExten+1:end-NumExten);
phasetf = phasetf(:,NumExten+1:end-NumExten);

% %%% Compute the phase transform
% %%% The factor 2*pi here comes from the specific form of Discrete Fourier
% %%% transform, not the 2*pi in the angular frequency
% phasetf = imag(dstftcfs./stftcfs)./(2*pi);
% 
% %%% Threshold for synchrosqueezing
% phasetf(abs(phasetf)<params.thr | isinf(abs(phasetf))) = NaN;
% 
% 
% Tx = 1/numVoices*sstalgo(stftcfs,phasetf,params.thr,params.sqType,freq,dt);
% 
% if (nargout == 0)
%     plotsst(Tx,freq,dt,params.engunitflag,params.normalizedfreq,Units);
% else
%     sst = Tx;
%     f = freq;
% end

end

%-------------------------------------------------------------------------
function [winfunc,dwinfunc] = sstwinfunc(WIN,winparam)

if strcmpi(WIN,'hermite')
    NumPts = winparam.NumPts;
    Order = winparam.Order;
    HalfWinSpt = winparam.HalfWinSpt; %%% half-time support
    disp('Hermite Window Parameters:');
    disp(['NumPts = ' num2str(NumPts)]);
    disp(['Order = ' num2str(Order)]);
    disp(['HalfWinSpt = ' num2str(HalfWinSpt)]);
    [h, Dh, ~] = hermf(NumPts, Order, HalfWinSpt);
    winfunc = h(Order,:)';
    dwinfunc = Dh(Order,:)';
end

end

%-------------------------------------------------------------------------
function plotsst(Tx,F,dt,engunitflag,isfreqnormalized,Units)

if ~isempty(Units)
    freqUnits = Units(1:end-1);
end

t = 0:dt:(size(Tx,2)*dt)-dt;
if engunitflag && isfreqnormalized
    frequnitstrs = wgetfrequnitstrs;
    freqlbl = frequnitstrs{1};
    xlbl = 'Samples';
elseif engunitflag && ~isfreqnormalized
    [F,~,uf] = wengunits(F,'unicode');
    freqlbl = wgetfreqlbl([uf 'Hz']);
    [t,~,ut] = wengunits(t,'unicode','time');
    xlbl = [getString(message('Wavelet:getfrequnitstrs:Time')) ' (' ut ')']; 
else
    freqlbl = getString(message('Wavelet:synchrosqueezed:FreqLabel'));
    freqlbl = ...
        [freqlbl '/' freqUnits ')'];
    xlbl = getString(message('Wavelet:synchrosqueezed:Time'));
    xlbl = [xlbl ' (' Units ')'];
end

h = pcolor(t,F,qclamp(abs(Tx), 0.002, 0.998));
h.EdgeColor = 'none';
shading interp;
colormap(1-gray);
ylabel(freqlbl); xlabel(xlbl);
title(getString(message('Wavelet:synchrosqueezed:SynchrosqueezedTitle')));

end

%-------------------------------------------------------------------------
function params = parseinputs(varargin)
% Set defaults.
params.fs = [];
params.dt = 1;
params.Ts = [];
params.sampinterval = false;
params.engunitflag = true;
params.WIN = 'hermite';
params.thr = 1e-8;
params.pad = false;
params.normalizedfreq = true;

% Error out if there are any calendar duration objects
tfcalendarDuration = cellfun(@iscalendarduration,varargin);
if any(tfcalendarDuration)
    error(message('Wavelet:FunctionInput:CalendarDurationSupport'));
end

tfsampinterval = cellfun(@isduration,varargin);

if (any(tfsampinterval) && nnz(tfsampinterval) == 1)
    params.sampinterval = true;
    params.Ts = varargin{tfsampinterval>0};
    if (numel(params.Ts) ~= 1 ) || params.Ts <= 0 || isempty(params.Ts)
        error(message('Wavelet:FunctionInput:PositiveScalarDuration'));
    end
    
    params.engunitflag = false;
    params.normalizedfreq = false;
    varargin(tfsampinterval) = [];
end

%%% Look for Name-Value pairs
freqbounds = find(strncmpi('FreqBounds',varargin,1));
if any(freqbounds)
    params.FreqBounds= varargin{freqbounds+1};
    varargin(freqbounds:freqbounds+1) = [];
    if isempty(varargin)
        return;
    end
end


freqres = find(strncmpi('freqres',varargin,1));

if any(freqres)
    params.FreqRes = varargin{freqres+1};
    %validate the value is numeric
    validateattributes(params.FreqRes,{'numeric'},{'positive','scalar'},...
        'fsstgao','FreqRes');
    varargin(freqres:freqres+1) = [];
    if isempty(varargin)
        return;
    end
end


extendsignal = find(strncmpi('extendsignal',varargin,1));

if any(extendsignal)
    params.pad = varargin{extendsignal+1};
    %validate the value is logical
    if ~isequal(params.pad,logical(params.pad))
        error(message('Wavelet:FunctionInput:Logical'));
    end
    varargin(extendsignal:extendsignal+1) = [];
    if isempty(varargin)
        return;
    end
end


windowparameters = find(strncmpi('windowparameters',varargin,1));
if any(windowparameters)
    params.winparam = varargin{windowparameters+1};
    
    varargin(windowparameters:windowparameters+1) = [];
    if isempty(varargin)
        return;
    end
end


%%% Only scalar left must be sampling frequency
tfsampfreq = cellfun(@(x) (isscalar(x) && isnumeric(x)),varargin);

if (any(tfsampfreq) && (nnz(tfsampfreq) == 1) && ~params.sampinterval)
    params.fs = varargin{tfsampfreq};
    validateattributes(params.fs,{'numeric'},{'positive'},'fsstgao','Fs');
    params.normalizedfreq = false;
    params.engunits = true;
elseif any(tfsampfreq) && params.sampinterval
    error(message('Wavelet:FunctionInput:SamplingIntervalOrDuration'));
elseif nnz(tfsampfreq)>1
    error(message('Wavelet:FunctionInput:Invalid_ScalNum'));
end

%%% Only char variable left must be window
tfwin = cellfun(@ischar,varargin);
if (nnz(tfwin) == 1)
    params.WIN = varargin{tfwin>0};
    params.WIN = validatestring(params.WIN,{'hermite'},'fsstgao','WIN');
elseif nnz(tfwin)>1
    error(message('Wavelet:FunctionInput:InvalidChar'));
end

if strncmpi(params.WIN,'hermite',1)
    if ~isfield(params.winparam, 'NumPts')
        params.winparam.NumPts = 377;
    end
    if ~isfield(params.winparam, 'Order')
        params.winparam.Order = 1;
    end
    if ~isfield(params.winparam, 'HalfWinSpt')
        params.winparam.HalfWinSpt = 10;
    end
end

end

%------------------------------------------------------------------------
function Tx = sstalgo(stftcfs,phasetf,gamma,sqType,freqs,dt)
%%% The validity of freq should have been checked earlier in the main
%%% rountine "fsstgao"
%%% phasetf is really the "angular frequency", with 2*pi already multiplied
%%% The right way to handle freqs is first to normalize them appropriately
%%% so they belong to the unit interval, then multiply them by 2*pi

% M = size(cwtcfs,1);
N = size(stftcfs,2);
numFreqs = length(freqs);
phasetf = phasetf / (2*pi);
freqs = freqs*dt; %%% this is to divide each frequency by Fs

if strcmpi(sqType, 'log')
    log2Fund = log2(max(1/N,min(freqs)));
    log2Nyquist = log2(min(1/2,max(freqs)));
    iRow = real(1 + floor(numFreqs/(log2Nyquist-log2Fund)*(log2(phasetf)-log2Fund)));
elseif strcmpi(sqType, 'linear')
    FundFreq = max(1/N,min(freqs));
    NyquistFreq = min(1/2,max(freqs));
    iRow = real(1 + floor(numFreqs/(NyquistFreq-FundFreq)*(phasetf-FundFreq)));
end

idxphasetf = find(iRow>0 & iRow<=numFreqs & ~isnan(iRow));
idxcwtcfs = find(abs(stftcfs)>gamma);
idx = intersect(idxphasetf,idxcwtcfs);
iCol = repmat(1:size(stftcfs,2),size(stftcfs,1),1);
Tx = accumarray([iRow(idx) iCol(idx)],stftcfs(idx),[numFreqs, size(stftcfs,2)]);

end
