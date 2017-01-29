%% preparation
close all;clc;
clearvars
path(pathdef);
addpath(genpath(strcat(pwd, '/utils')), path);
scrsz = get(groot,'ScreenSize');

%% set up parameters
Fs = 100;
tSamples = (1/Fs:1/Fs:16)';
N = length(tSamples);
FreqBounds = [0,0.5];
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
% if1 = 35 + 6 * if1 ./ max(abs(if1));
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

targetSig = clean;

% plot(targetSig);

%% SST-CWT
tic;
[sstResult, instFreqTic, cwtcfs, cwtscales, ~, psift, nativeOmega, dpsift, dcwtcfs] = wsstgao(targetSig, Fs, 'morse',...
    'VoicesPerOctave', 32, 'ExtendSignal', false,...
    'WaveletParameters', struct('be',80,'ga',1));
toc;

itvPS = abs(Fs*sstResult/2).^2;
itvPS_logscale = qclamp(log(1+itvPS), 0.005);

%% SST-CWT with Gaurav-Wu implementation
opts = struct('motherwavelet', 'morse', 'k', 0, 'beta', 80, 'gam', 1);
% FreqRes = (Fs/2)/length(instFreqTic);

tic;
% [CWTResult, cwt_tic, sqCWTResult, sqcwt_tic] = sqCWTbase(tSamples, targetSig, 0, Fs/2, FreqRes, opts, 0, 1);
[CWTResult, cwt_tic, sqCWTResult, sqcwt_tic, sqPsift, sqOmega, sqDPsift, DCWTResult] = sqCWTbase(tSamples, targetSig, FreqBounds(1)*Fs, FreqBounds(2)*Fs, FreqRes*Fs, opts, 0, 1);
toc;

Conceft_itvPS = abs(Fs*sqCWTResult/2).^2;
Conceft_itvPS_logscale = qclamp(log(1+Conceft_itvPS), 0.001);

%% compare DFT of wavelets
figure('Position',[1 scrsz(4)/2 scrsz(3)*2 scrsz(4)], 'Name', 'DFT of Wavelets');
subplot(1,2,1);
imagesc(abs(psift));
title('internal', 'Interpreter', 'latex', 'fontsize', 20);
subplot(1,2,2);
imagesc(abs(sqPsift));
% axis xy
title('Gaurav-Wu', 'Interpreter', 'latex', 'fontsize', 20);

%% compare DFT of derivatives of wavelets
figure('Position',[1 scrsz(4)/2 scrsz(3)*2 scrsz(4)], 'Name', 'DFT of Derivatives of Wavelets');
subplot(1,2,1);
imagesc(abs(dpsift));
title('internal', 'Interpreter', 'latex', 'fontsize', 20);
subplot(1,2,2);
imagesc(abs(sqDPsift));
title('Gaurav-Wu', 'Interpreter', 'latex', 'fontsize', 20);

%% compare phase transforms
figure('Position',[1 scrsz(4)/2 scrsz(3)*2 scrsz(4)]);
subplot(1,2,1);
imagesc(nativeOmega);
title('internal', 'Interpreter', 'latex', 'fontsize', 20);
subplot(1,2,2);
imagesc(sqOmega');
title('Gaurav-Wu', 'Interpreter', 'latex', 'fontsize', 20);

%% compare wavelet coefficients
figure('Position',[1 scrsz(4)/2 scrsz(3)*2 scrsz(4)],'Name','CWT coefficients');
subplot(1,2,1);
imagesc(abs(cwtcfs));
% pcolor(tSamples, cwtscales, abs(cwtcfs));
% imagesc(abs(cwtcfs));
% axis xy
% shading interp
title('internal', 'Interpreter', 'latex', 'fontsize', 20);
subplot(1,2,2);
imagesc(abs(CWTResult));
% imagesc(abs(CWTResult(end:-1:1,:)));
% axis xy
% shading interp
title('Gaurav-Wu', 'Interpreter', 'latex', 'fontsize', 20);

%% compare derivaties of wavelet coefficients
figure('Position', [1 scrsz(4)/2 scrsz(3)*2 scrsz(4)], 'Name', 'Derivatives of CWT coefficients');
subplot(1,2,1);
imagesc(abs(dcwtcfs));
title('internal', 'Interpreter', 'latex', 'fontsize', 20);
subplot(1,2,2);
imagesc(abs(DCWTResult));
title('Gaurav-Wu', 'Interpreter', 'latex', 'fontsize', 20);

%% visualize and compare SST results
% close(gcf);

hq = figure('Position',[10 10 2*scrsz(3)/3 scrsz(4)/3]);

subplot(2,1,1);
% pcolor(tSamples, instFreqTic, log(1+abs(itvPS)));
% [~,freqUBIdx] = min(abs(instFreqTic-0.2*Fs));
% pcolor(tSamples, instFreqTic(1:freqUBIdx), itvPS_logscale(1:freqUBIdx,:));
pcolor(tSamples, instFreqTic, itvPS_logscale);
shading interp;
colormap(1-gray);
% axis([0 inf FreqBounds(1)*Fs FreqBounds(2)*Fs]);
xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
title('SST-CWT', 'Interpreter', 'latex', 'fontsize', 20);
hold on
plot(tSamples, if1, 'r', 'linewidth', 1);
plot(tSamples, if2, 'b', 'linewidth', 1);

subplot(2,1,2);
% pcolor(tSamples, sqcwt_tic, log(1+abs(Conceft_itvPS)));
% [~,freqUBIdx] = min(abs(sqcwt_tic-0.2*Fs));
% pcolor(tSamples, sqcwt_tic(1:freqUBIdx), Conceft_itvPS_logscale(1:freqUBIdx,:));
pcolor(tSamples, sqcwt_tic, Conceft_itvPS_logscale);
% pcolor(tSamples, sqcwt_tic, abs(sqCWTResult));
% pcolor(tSamples,(1:length(sqcwt_tic))/length(sqcwt_tic), abs(sqCWTResult));
shading interp
colormap(1-gray);
xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
hold on
plot(tSamples, if1, 'r', 'linewidth', 1);
plot(tSamples, if2, 'b', 'linewidth', 1);
