%% preparation
close all;clc;
clearvars
path(pathdef);
addpath(genpath(strcat(pwd, '/utils')), path);
scrsz = get(groot,'ScreenSize');

%% set up parameters
Fs = 100;
T = 16;
tSamples = (1/Fs:1/Fs:T)';
N = length(tSamples);
FreqBounds = [0,0.2];
FreqRes = 2e-4;

%% generate random signal using the model in Conceft paper
initstate(1);

am1 = smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess');
am1 = 2 + am1 ./ max(abs(am1));
am2 = smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess');
am2 = 2 + am2 ./ max(abs(am2));
am1(1:500) = 0;
am2(end-600:end) = 0;

if1 = smooth(cumsum(randn(N,1)) ./ Fs, 400, 'loess');
if1 = 10 + 6 * if1 ./ max(abs(if1));
if2 = smooth(cumsum(randn(N,1)) ./ Fs, 300, 'loess');
if2 = pi + 3 * if2 ./ max(abs(if2));
phi1 = cumsum(if1) / Fs;
phi2 = cumsum(if2) / Fs;

s1 = am1 .* cos(2*pi*phi1); 
s2 = am2 .* cos(2*pi*phi2); 
clean = s1 + s2;

if1(1:500) = nan;
if2(end-600:end) = nan;
am1(1:500) = nan;
am2(end-600:end) = nan;

sigma = 1;
noise = random('T',4,N,1);
noise = sigma * noise; 
fprintf('var(noise) = %f\n', var(noise));
snrdb = 20 * log10(std(clean)./std(noise));
fprintf('snrdb = %f\n',snrdb);

xm = clean + noise;

targetSig = xm;

% plot(targetSig);

%% SST-CWT
tic;
[sstResult, instFreqTic] = wsstgao(targetSig, Fs, 'morse',...
    'VoicesPerOctave', 32, 'ExtendSignal', false,...
    'WaveletParameters', struct('be',80,'ga',1,'k',0),...
    'SqType', 'linear', 'FreqBounds', FreqBounds*Fs, 'FreqRes', FreqRes*Fs);
toc;

itvPS = abs(Fs*sstResult/2).^2;
itvPS_logscale = qclamp(log(1+itvPS), 0.001);

%% visualize and compare SST results

hq = figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);

pcolor(tSamples, instFreqTic, itvPS_logscale);
shading interp;
colormap(1-gray);
axis([0 inf FreqBounds(1)*Fs FreqBounds(2)*Fs]);
xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
title('SST-CWT', 'Interpreter', 'latex', 'fontsize', 20);
hold on
plot(tSamples, if1, 'r', 'linewidth', 1);
plot(tSamples, if2, 'b', 'linewidth', 1);
