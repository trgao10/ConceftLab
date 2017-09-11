%% preparation
close all;clc;
clearvars
path(pathdef);
addpath(genpath(strcat(pwd, '/utils')), path);
scrsz = get(groot,'ScreenSize');

%%% also use the original code accompanying the critical paper
addpath(genpath(strcat(pwd, '/external')), path);

%% Simulation 1: Resolution of Two Tones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% set up parameters
Fs = 20;
T = 500;

nu1 = 2;
nu2 = 2.25;
r = 0.5;
phi1 = 0;
phi2 = 0;

%%%% cosntruct signal
tSamples = (1/Fs:1/Fs:T)';
s = cos(2*pi*nu1*tSamples+phi1)+r*cos(2*pi*nu2*tSamples+phi2);

%%%% window parameter
f0 = 0.4; %% the critical paper used f0 in [7, 2, 1, 0.4]

%%%% stft frequency range and resolution
FreqBounds = [0,0.5];
FreqRes = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use matlab built-in stft
win_truncate_threshold = 3;
hWinLen = fix(Fs*f0*sqrt(2*win_truncate_threshold));
gwin = @(x) exp(-x.^2/(2*f0^2));
wfunc = gwin(-hWinLen/Fs:1/Fs:hWinLen/Fs);
wfunc = wfunc / max(abs(wfunc));
nfft = fix(0.5/FreqRes); %% determine frequency resolution

[stft,freqLabel,timeLabel] = spectrogram(s,wfunc,length(wfunc)-1,2*nfft+1,Fs,'yaxis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use the code from the critical paper
[WFT,freq]=wft(s,Fs,'fstep',FreqRes*Fs,...
    'fmin',FreqBounds(1)*Fs,'fmax',FreqBounds(2)*Fs,'f0',f0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make comparison plots
figure('Name', 'Two Tones');

subplot(1,2,1);
pc1 = pcolor(tSamples,freq,abs(WFT));
set(pc1,'EdgeColor','none');
title('critical paper');
xlabel(sprintf('num of cols: %d',size(WFT,2)));
ylabel(sprintf('num of rows: %d',size(WFT,1)));
hold on
plot(tSamples,nu1*ones(size(tSamples)),'r');
plot(tSamples,nu2*ones(size(tSamples)),'m');
set(gca,'YLim',[0,4]);
set(gca,'XLim',[100,120]);

subplot(1,2,2);
pc2 = pcolor(timeLabel,freqLabel,abs(stft));
set(pc2,'EdgeColor','none');
title('matlab built-in');
xlabel(sprintf('num of cols: %d',size(stft,2)));
ylabel(sprintf('num of rows: %d',size(stft,1)));
hold on
plot(tSamples,nu1*ones(size(tSamples)),'r');
plot(tSamples,nu2*ones(size(tSamples)),'m');
set(gca,'YLim',[0,4]);
set(gca,'XLim',[100,120]);

%% Simulation 2: Amplitude Modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% set up parameters
Fs = 20;
T = 500;

nu1 = 0.25;
nu2 = 2;
r = 0.5;

%%%% cosntruct signal
tSamples = (1/Fs:1/Fs:T)';
s = (1+0.5*cos(2*pi*nu1*tSamples)).*cos(2*pi*nu2*tSamples);

%%%% window parameter
f0 = 0.4; %% the critical paper used f0 in [7, 2, 1, 0.4]

%%%% stft frequency range and resolution
FreqBounds = [0,0.5];
FreqRes = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use matlab built-in stft
win_truncate_threshold = 3;
hWinLen = fix(Fs*f0*sqrt(2*win_truncate_threshold));
gwin = @(x) exp(-x.^2/(2*f0^2));
wfunc = gwin(-hWinLen/Fs:1/Fs:hWinLen/Fs);
wfunc = wfunc / max(abs(wfunc));
nfft = fix(0.5/FreqRes); %% determine frequency resolution

[stft,freqLabel,timeLabel] = spectrogram(s,wfunc,length(wfunc)-1,2*nfft+1,Fs,'yaxis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use the code from the critical paper
[WFT,freq]=wft(s,Fs,'fstep',FreqRes*Fs,...
    'fmin',FreqBounds(1)*Fs,'fmax',FreqBounds(2)*Fs,'f0',f0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make comparison plots
figure('Name', 'Amplitude Modulation');

subplot(1,2,1);
pc1 = pcolor(tSamples,freq,abs(WFT));
set(pc1,'EdgeColor','none');
title('critical paper');
xlabel(sprintf('num of cols: %d',size(WFT,2)));
ylabel(sprintf('num of rows: %d',size(WFT,1)));
set(gca,'YLim',[0.1,4]);
set(gca,'XLim',[40,60]);

subplot(1,2,2);
pc2 = pcolor(timeLabel,freqLabel,abs(stft));
set(pc2,'EdgeColor','none');
title('matlab built-in');
xlabel(sprintf('num of cols: %d',size(stft,2)));
ylabel(sprintf('num of rows: %d',size(stft,1)));
set(gca,'YLim',[0.1,4]);
set(gca,'XLim',[40,60]);

%% Simulation 3: Frequency Modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% set up parameters
Fs = 20;
T = 500;

nu1 = 2;
nu2 = 0.25;

%%%% cosntruct signal
tSamples = (1/Fs:1/Fs:T)';
s = cos(2*pi*nu1*tSamples + sin(2*pi*nu2*tSamples));

%%%% window parameter
f0 = 7; %% the critical paper used f0 in [7, 2, 1, 0.4]

%%%% stft frequency range and resolution
FreqBounds = [0,0.5];
FreqRes = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use matlab built-in stft
win_truncate_threshold = 3;
hWinLen = fix(Fs*f0*sqrt(2*win_truncate_threshold));
gwin = @(x) exp(-x.^2/(2*f0^2));
wfunc = gwin(-hWinLen/Fs:1/Fs:hWinLen/Fs);
wfunc = wfunc / max(abs(wfunc));
nfft = fix(0.5/FreqRes); %% determine frequency resolution

[stft,freqLabel,timeLabel] = spectrogram(s,wfunc,length(wfunc)-1,2*nfft+1,Fs,'yaxis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use the code from the critical paper
[WFT,freq]=wft(s,Fs,'fstep',FreqRes*Fs,...
    'fmin',FreqBounds(1)*Fs,'fmax',FreqBounds(2)*Fs,'f0',f0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make comparison plots
figure('Name', 'Frequency Modulation');

subplot(1,2,1);
pc1 = pcolor(tSamples,freq,abs(WFT));
set(pc1,'EdgeColor','none');
title('critical paper');
xlabel(sprintf('num of cols: %d',size(WFT,2)));
ylabel(sprintf('num of rows: %d',size(WFT,1)));
set(gca,'YLim',[0.1,4]);
set(gca,'XLim',[40,60]);

subplot(1,2,2);
pc2 = pcolor(timeLabel,freqLabel,abs(stft));
set(pc2,'EdgeColor','none');
title('matlab built-in');
xlabel(sprintf('num of cols: %d',size(stft,2)));
ylabel(sprintf('num of rows: %d',size(stft,1)));
set(gca,'YLim',[0.1,4]);
set(gca,'XLim',[40,60]);

%% Simulation 4: Noise and Its Effects
%%%% set up parameters
Fs = 50;
T = 100;

nu = 2;
noiseLevel = 0.5;

%%%% cosntruct signal
tSamples = (1/Fs:1/Fs:T)';
s = cos(2*pi*nu*tSamples) + (noiseLevel/sqrt(2))*wgn(length(tSamples),1,0);

%%%% window parameter
f0 = 7; %% the critical paper used f0 in [7, 2, 1, 0.4]

%%%% stft frequency range and resolution
FreqBounds = [0,0.5];
FreqRes = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use matlab built-in stft
win_truncate_threshold = 3;
hWinLen = fix(Fs*f0*sqrt(2*win_truncate_threshold));
gwin = @(x) exp(-x.^2/(2*f0^2));
wfunc = gwin(-hWinLen/Fs:1/Fs:hWinLen/Fs);
wfunc = wfunc / max(abs(wfunc));
nfft = fix(0.5/FreqRes); %% determine frequency resolution

[stft,freqLabel,timeLabel] = spectrogram(s,wfunc,length(wfunc)-1,2*nfft+1,Fs,'yaxis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use the code from the critical paper
[WFT,freq]=wft(s,Fs,'fstep',FreqRes*Fs,...
    'fmin',FreqBounds(1)*Fs,'fmax',FreqBounds(2)*Fs,'f0',f0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make comparison plots
figure('Name', 'Frequency Modulation');

subplot(1,2,1);
pc1 = pcolor(tSamples,freq,abs(WFT));
set(pc1,'EdgeColor','none');
title('critical paper');
xlabel(sprintf('num of cols: %d',size(WFT,2)));
ylabel(sprintf('num of rows: %d',size(WFT,1)));
set(gca,'YLim',[0.1,4]);
set(gca,'XLim',[40,60]);

subplot(1,2,2);
pc2 = pcolor(timeLabel,freqLabel,abs(stft));
set(pc2,'EdgeColor','none');
title('matlab built-in');
xlabel(sprintf('num of cols: %d',size(stft,2)));
ylabel(sprintf('num of rows: %d',size(stft,1)));
set(gca,'YLim',[0.1,4]);
set(gca,'XLim',[40,60]);
