function [sst,f] = wsstgao(x,varargin)
%Wavelet Synchrosqueezed Transform (adapted from the wsst function shipped
%with the MATLAB Wavelet Toolbox)
%   Tingran Gao (trgao10@math.duke.edu)
%   Jan 27, 2017 (Lunar Chinese New Year's Eve)
%   
%   SST = wsstgao(X) returns the wavelet synchrosqueezed transform for the
%   1-D real-valued signal, X. X must have at least 4 samples. The
%   synchrosqueezed transform uses nv (32 by default) voices per octave and
%   the number of octaves is no greater than floor(log2(numel(X)))-1.
%   The transform uses the analytic Morse wavelet by default.
%   SST is a Na-by-N matrix where Na is the total number of scales, i.e.
%   no greater than nv*(floor(log2(numel(X)))-1), and N is the number of
%   samples in signal X.
%
%   [SST,F] = wsstgao(X) returns a vector of frequencies, F, in
%   cycles/sample corresponding to the rows of SST.
%
%   [...] = wsstgao(X,Fs) specifies the sampling frequency, Fs, in hertz as
%   a positive scalar. If you specify the sampling frequency, WSST returns
%   the frequencies in hertz.
%
%   [...] = wsstgao(...,WAV) uses the analytic wavelet specified by WAV to
%   compute the synchrosqueezed transform. Valid choices for WAV are
%   'amor' and 'bump' for the analytic Morlet and bump wavelet. If
%   unspecified, WAV defaults to 'amor'.
%
%   [...] = wsstgao(...,'VoicesPerOctave',NV) specifies the number of
%   voices per octave as a positive even integer between 10 and 48.
%   The number of scales is the product of the number of voices per octave
%   and the number of octaves. If unspecified, NV defaults to 32 voices
%   per octave.
%
%   [...] = wsstgao(...,'ExtendSignal',EXTENDFLAG) specifies whether to
%   symmetrically extend the signal by reflection to mitigate boundary
%   effects. EXTENDFLAG can be one of the following options [ true |
%   {false}]. If unspecified, EXTENDSIGNAL defaults to false.  You can
%   specify the 'ExtendSignal' name-value pair anywhere in the input
%   argument list after the signal X.
%
%   wsstgao(...) with no output arguments plots the wavelet synchrosqueezed
%   transform as a function of time and frequency. If you do not specify a
%   sampling frequency or interval, the synchrosqueezed transform is
%   plotted in cycles/sample. If you supply a sampling frequency, Fs, the
%   synchrosqueezed transform is plotted in hertz. If you supply a
%    <a href="matlab:help duration">duration</a> as a sampling interval,
%   the synchrosqueezed transform is plotted in cycles/unit time where the
%   time unit is the same as the duration.
%
%
%

narginchk(1,15);
nbSamp = numel(x);
x = x(:)';
validateattributes(x,{'double'},{'row','finite','real'},'wsst','X');
if numel(x)<4
    error(message('Wavelet:synchrosqueezed:NumInputSamples'));
end
params = parseinputs(nbSamp,varargin{:});
numVoices = params.nv;
numOctaves = params.noct;
%%% Create scale vector
numScales = numOctaves*params.nv;


%%% If sampling frequency is specified, dt = 1/fs
if (isempty(params.fs) && isempty(params.Ts))
    % The default is 1 for normalized frequency
    dt = params.dt;
    Units = '';
elseif (~isempty(params.fs) && isempty(params.Ts))
    % Accept the sampling frequency in hertz
    fs = params.fs;
    dt = 1/fs;
    Units = '';
elseif (isempty(params.fs) && ~isempty(params.Ts))
    % Get the dt and Units from the duration object
    [dt,Units] = getDurationandUnits(params.Ts);
else
end

% T = nbSamp*dt;

a0 = 2^(1/numVoices);
scales = a0.^(1:numScales);
NbSc = numel(scales);

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

%%% Create frequency vector for CWT computation
omega = (1:fix(N/2));
omega = omega*((2.*pi)/N); %%% omage so obtained are angular frequencies (with 2*pi multiplied)
omega = [0., omega, -omega(fix((N-1)/2):-1:1)];

%%% Compute FFT of the (padded) time series
xdft = fft(x);
[psift,dpsift,cwtfreq] = sstwaveft(params.WAV,omega,scales,params.wavparam);

cwtscales = scales;

%%% Obtain CWT coefficients and derivatives
cwtcfs = ifft(repmat(xdft,NbSc,1).*psift,[],2);
dcwtcfs = ifft(repmat(xdft,NbSc,1).*dpsift,[],2);

%%% Remove padding if any
cwtcfs = cwtcfs(:,NumExten+1:end-NumExten);
dcwtcfs = dcwtcfs(:,NumExten+1:end-NumExten);

%%% Compute the phase transform
%%% The factor 2*pi here comes from the specific form of Discrete Fourier
%%% transform, not the 2*pi in the angular frequency
phasetf = imag(dcwtcfs./cwtcfs)./(2*pi);

%%% Threshold for synchrosqueezing
phasetf(abs(phasetf)<params.thr | isinf(abs(phasetf))) = NaN;

%%% Create frequency vector for output
if isfield(params, 'FreqBounds') && isfield(params, 'FreqRes')
    if strcmpi(params.sqType, 'log')
        log2LowFreq = log2(max(min(params.FreqBounds), 1/(nbSamp*dt)));
        log2HighFreq = log2(min(max(params.FreqBounds), fs/2));
        numFreqs = fix(abs(log2HighFreq-log2LowFreq) / (params.FreqRes*dt))+1;
        freq = 2.^linspace(log2LowFreq,log2HighFreq,numFreqs);
    else
        LowFreq = max(min(params.FreqBounds), 1/(nbSamp*dt));
        HighFreq = min(max(params.FreqBounds), fs/2);
        numFreqs = fix(abs(HighFreq-LowFreq) / params.FreqRes)+1;
        freq = linspace(LowFreq, HighFreq, numFreqs);
    end
else
    params.sqType = 'log';
    log2Nyquist = log2(fs/2);
    log2Fund = log2(1/(nbSamp*dt));
    freq = 2.^linspace(log2Fund,log2Nyquist,numScales);
end

Tx = 1/numVoices*sstalgo(cwtcfs,phasetf,params.thr,params.sqType,freq,dt);

if (nargout == 0)
    plotsst(Tx,freq,dt,params.engunitflag,params.normalizedfreq,Units);
else
    sst = Tx;
    f = freq;
end

end

%-------------------------------------------------------------------------
function [wft,dwft,freq] = sstwaveft(WAV,omega,scales,wavparam)
%   Admissible wavelets are:
%   - MORSE wavelet - 'morse':
%       PSI_HAT(s) = (involving Laguerre polynomials)
%       Parameters: be (beta), ga (gamma), default: be = 30, ga = 3.
%   - MORLET wavelet (A) - 'amor':
%       PSI_HAT(s) = exp(-(s-s0).^2/2) * (s>0)
%       Parameter: s0, default: s0 = 6.
%   - Bump wavelet:  'bump':
%       PSI_HAT(s) = exp(1-(1/((s-mu)^2./sigma^2))).*(abs((s-mu)/sigma)<1)
%       Parameters: mu,sigma, default: mu=5, sigma = 1.
%   Normalized to have unit magnitude at the peak frequency of the wavelet

NbSc = numel(scales);

if strcmpi(WAV,'morse')
    ga = wavparam.ga;
    be = wavparam.be;
    k = wavparam.k;
    disp('Morse Wavelet Parameters:');
    disp(['ga = ' num2str(ga)]);
    disp(['be = ' num2str(be)]);
    disp(['k = ' num2str(k)]);
    [wft,freq] = morsewavft(omega,scales,ga,be,k);
else
    [wft,freq] = waveft(WAV,omega,scales);
end

%%% Compute derivatives
omegaMatrix = repmat(omega,NbSc,1);
dwft = 2*pi*1j*omegaMatrix.*wft;

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
function params = parseinputs(nbSamp,varargin)
% Set defaults.
params.fs = [];
params.dt = 1;
params.Ts = [];
params.sampinterval = false;
params.engunitflag = true;
params.WAV = 'amor';
params.thr = 1e-8;
params.nv = 32;
params.noct = floor(log2(nbSamp))-1;
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
        'wsstgao','VoicesPerOctave');
    varargin(freqres:freqres+1) = [];
    if isempty(varargin)
        return;
    end
end


numvoices = find(strncmpi('voicesperoctave',varargin,1));

if any(numvoices)
    params.nv = varargin{numvoices+1};
    %validate the value is numeric
    validateattributes(params.nv,{'numeric'},{'positive','scalar',...
        'even','>=',10,'<=',1024},'wsstgao','VoicesPerOctave');
    varargin(numvoices:numvoices+1) = [];
    if isempty(varargin)
        return;
    end
end


sqtype = find(strncmpi('SqType',varargin,1));

if any(sqtype)
    params.sqType = varargin{sqtype+1};
    params.sqType = validatestring(params.sqType,{'linear','log'},'wsstgao','sqType');
    varargin(sqtype:sqtype+1) = [];
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


waveletparameters = find(strncmpi('waveletparameters',varargin,1));
if any(waveletparameters)
    params.wavparam = varargin{waveletparameters+1};
    
    varargin(waveletparameters:waveletparameters+1) = [];
    if isempty(varargin)
        return;
    end
end


%%% Only scalar left must be sampling frequency
tfsampfreq = cellfun(@(x) (isscalar(x) && isnumeric(x)),varargin);

if (any(tfsampfreq) && (nnz(tfsampfreq) == 1) && ~params.sampinterval)
    params.fs = varargin{tfsampfreq};
    validateattributes(params.fs,{'numeric'},{'positive'},'wsstgao','Fs');
    params.normalizedfreq = false;
    params.engunits = true;
elseif any(tfsampfreq) && params.sampinterval
    error(message('Wavelet:FunctionInput:SamplingIntervalOrDuration'));
elseif nnz(tfsampfreq)>1
    error(message('Wavelet:FunctionInput:Invalid_ScalNum'));
end

%%% Only char variable left must be wavelet
tfwav = cellfun(@ischar,varargin);
if (nnz(tfwav) == 1)
    params.WAV = varargin{tfwav>0};
    params.WAV = validatestring(params.WAV,{'bump','amor','morse'},'wsstgao','WAV');
elseif nnz(tfwav)>1
    error(message('Wavelet:FunctionInput:InvalidChar'));
end

if strncmpi(params.WAV,'morse',1)
    if ~isfield(params.wavparam, 'be')
        params.wavparam.be = 30;
    end
    if ~isfield(params.wavparam, 'ga')
        params.wavparam.ga = 3;
    end
    if ~isfield(params.wavparam, 'k')
        params.wavparam.k = 0;
    end
elseif strncmpi(params.WAV,'amor',1)
    if ~isfield(params.wavparam, 'cf')
        params.wavparam.cf = 6;
    end
elseif strncmpi(params.WAV,'bump',1)
    if ~isfield(params.wavparam, 'mu')
        params.wavparam.mu = 5;
    end
    if ~isfield(params.wavparam, 'sigma')
        params.wavparam.sigma = 1;
    end
end

end

%------------------------------------------------------------------------
function Tx = sstalgo(cwtcfs,phasetf,gamma,sqType,freqs,dt)
%%% The validity of freq should have been checked earlier in the main
%%% rountine "wsstgao"
%%% phasetf is really the "angular frequency", with 2*pi already multiplied
%%% The right way to handle freqs is first to normalize them appropriately
%%% so they belong to the unit interval, then multiply them by 2*pi

% M = size(cwtcfs,1);
N = size(cwtcfs,2);
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
idxcwtcfs = find(abs(cwtcfs)>gamma);
idx = intersect(idxphasetf,idxcwtcfs);
iCol = repmat(1:size(cwtcfs,2),size(cwtcfs,1),1);
Tx = accumarray([iRow(idx) iCol(idx)],cwtcfs(idx),[numFreqs, size(cwtcfs,2)]);

end

