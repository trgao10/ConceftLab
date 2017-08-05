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
FreqBounds = [0,0.5];
FreqRes = 2e-4;

%% generate random signal using the model in Conceft paper
initstate(1);

am1 = ones(N,1);
% am1 = smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess');
% am1 = 2 + am1 ./ max(abs(am1));
% if1 = 30*ones(N,1);
if1 = linspace(0,50,N)';
% if1 = linspace(10,40,N)';
phi1 = cumsum(if1) / Fs;

s1 = am1 .* cos(2*pi*phi1); 
clean = s1;

sigma = 1;
noise = random('T',4,N,1);
noise = sigma * noise; 
fprintf('var(noise) = %f\n', var(noise));
snrdb = 20 * log10(std(clean)./std(noise));
fprintf('snrdb = %f\n',snrdb);

xm = clean + noise;

targetSig = clean;

% plot(targetSig);

%% SST-CWT
tic;
[sstResult, instFreqTic, stftcfs, phasetf, stftfreqs, reassntRule] =...
    fsstgao(targetSig, Fs, 'hermite',...
    'ExtendSignal', false,...
    'WindowParameters', struct('NumPts',377,'Order',5,'HalfWinSpt',10),...
    'FreqBounds', FreqBounds*Fs, 'FreqRes', FreqRes*Fs);
toc;

itvPS = abs(Fs*sstResult/2).^2;
itvPS_logscale = qclamp(log(1+itvPS), 0.001);

itvPS_stft = abs(Fs*stftcfs/2).^2;
itvPS_stft_logscale = qclamp(log(1+itvPS_stft), 0.0001);

%% visualize and compare SST results

hq = figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);

subplot(1,2,1);
pcolor(tSamples, instFreqTic, itvPS_stft_logscale);
shading interp;
colormap(1-gray);
axis([2 14 FreqBounds(1)*Fs FreqBounds(2)*Fs]);
xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
title('STFT', 'Interpreter', 'latex', 'fontsize', 20);
hold on
plot(tSamples, if1, 'r', 'linewidth', 1);
grid on

subplot(1,2,2);
pcolor(tSamples, instFreqTic, itvPS_logscale);
shading interp;
colormap(1-gray);
axis([2 14 FreqBounds(1)*Fs FreqBounds(2)*Fs]);
xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
title('SST-CWT', 'Interpreter', 'latex', 'fontsize', 20);
hold on
plot(tSamples, if1, 'r', 'linewidth', 1);
grid on

drawnow

%% visualize reassignment rule
reassntRuleMask = nan(size(reassntRule));
visStripWidth = 500;
for j=1:length(tSamples)
    cIdx = round((if1(j)/(Fs/2))*length(instFreqTic));    
    reassntRuleMask(max(1,(cIdx-visStripWidth)):min(length(instFreqTic),(cIdx+visStripWidth)),j) = 1;
end

reassntRuleVis = reassntRule.*reassntRuleMask;
% reassntRuleVis(isnan(reassntRuleVis)) = 0;
% reassntRuleVis(isnan(reassntRuleVis)) = 0;
reassntRuleVis_rescaled = (reassntRuleVis-mean(reassntRuleVis(~isnan(reassntRuleVis))))/std(reassntRuleVis(~isnan(reassntRuleVis)));

figure;
pcolor(tSamples, instFreqTic, reassntRuleVis_rescaled);
shading interp;
colormap(jet);
axis([2 14 FreqBounds(1)*Fs FreqBounds(2)*Fs]);
set(gca,'YDir','normal');
hold on
plot(tSamples, if1, 'r', 'linewidth', 1);
colorbar

