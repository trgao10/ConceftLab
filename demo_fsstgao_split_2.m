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

am1 = smooth(cumsum(randn(N,1)) ./ Fs, 200, 'loess');
am1 = 2 + am1 ./ max(abs(am1));
% if1 = 25*ones(N,1);
if1 = linspace(0,50,N)';
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
[sstResult, instFreqTic] = fsstgao(targetSig, Fs, 'hermite',...
    'ExtendSignal', false,...
    'WindowParameters', struct('NumPts',377,'Order',2,'HalfWinSpt',10),...
    'FreqBounds', FreqBounds*Fs, 'FreqRes', FreqRes*Fs);
toc;

itvPS = abs(Fs*sstResult/2).^2;
itvPS_logscale = qclamp(log(1+itvPS), 0.001);

%% visualize and compare SST results

hq = figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);

segIdx = 201:1400;

subplot(1,2,1);
pcolor(tSamples(segIdx), instFreqTic, itvPS_logscale(:,segIdx));
shading interp;
colormap(1-gray);
axis([tSamples(segIdx(1)) tSamples(segIdx(end)) FreqBounds(1)*Fs FreqBounds(2)*Fs]);
xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
title('SST-CWT', 'Interpreter', 'latex', 'fontsize', 20);
rectangle('Position', [6,16,3,12.5], 'EdgeColor', 'b', 'LineStyle', ':');

subplot(1,2,2);
pcolor(tSamples(segIdx), instFreqTic, itvPS_logscale(:,segIdx));
shading interp;
colormap(1-gray);
% axis([tSamples(segIdx(1)) tSamples(segIdx(end)) FreqBounds(1)*Fs FreqBounds(2)*Fs]);
xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
title('SST-CWT', 'Interpreter', 'latex', 'fontsize', 20);
axis([6,9,16,28.5]);

hold on
omega0 = if1(1);
omega1 = (if1(end)-if1(1))/T;
instFreq = omega0+omega1*linspace(0,T,N);
plot(tSamples, instFreq, 'g', 'linewidth', 1.0);
plot(tSamples, if1+omega1^2/sqrt(1+omega1^2)/(2*pi), 'r', 'linewidth', 1.0);
plot(tSamples, if1-omega1^2/sqrt(1+omega1^2)/(2*pi), 'r', 'linewidth', 1.0);
legend({'','$\omega_0+\omega_1t$','$\displaystyle\omega_0+\omega_1t\pm \frac{\omega_1^2}{1+\omega_1^2}\sqrt{1+\omega_1^2}$'}, 'Interpreter', 'latex', 'Fontsize', 15, 'Location', 'Southeast');
