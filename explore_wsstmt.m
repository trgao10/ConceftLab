%% preparation
close all;clc;
clearvars
path(pathdef);
addpath(genpath(strcat(pwd, '/utils')), path);
scrsz = get(groot,'ScreenSize');

%% set up parameters
Fs = 100;
T = 16;
MT = 4;
tSamples = (1/Fs:1/Fs:T)';
N = length(tSamples);
FreqBounds = [0,0.2];
FreqRes = 1e-4;

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

targetSig = clean;

% plot(targetSig);

%% SST-CWT
tic;
[sstmtResult, instFreqTic, sstmtCell, cwtcfsCell, phasetfCell, cwtfreqCell] = wsstmt(targetSig, MT, Fs, 'morse',...
    'VoicesPerOctave', 128, 'ExtendSignal', true,...
    'WaveletParameters', struct('be',30,'ga',3),...
    'SqType', 'linear', 'FreqBounds', FreqBounds*Fs, 'FreqRes', FreqRes*Fs);
toc;

%% visualize phase plots
hq = figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
segIdx = 201:1400;
[~,yAxisHigh] = min(abs(cwtfreqCell{1}*Fs-20));
yAxisIdx = yAxisHigh:length(cwtfreqCell{1});
for k = 1:MT
    subplot(ceil(MT/2),2,k);
    drawObj = cwtcfsCell{k}(yAxisIdx,segIdx);
    mask = abs(cwtcfsCell{k}(yAxisIdx,segIdx));
    mask = mask > mean(mask(:)+0.5*std(mask(:)));
    pcolor(tSamples(segIdx), (cwtfreqCell{k}(yAxisIdx))*Fs, angle(drawObj).*mask);
%     pcolor(tSamples(segIdx), (cwtfreqCell{k}(yAxisIdx))*Fs, abs(mask));

%     phaseTransPlot = phasetfCell{k}(yAxisIdx,segIdx)/(2*pi);
%     phaseTransPlot(phaseTransPlot < 0) = 0;
%     phaseTransPlot(abs(phaseTransPlot) > 1) = 0;
% %     cwtcfsPhasePlot = angle(cwtcfsCell{k}(yAxisIdx,segIdx));
% %     cwtcfsPhasePlot(phaseTransPlot < 0 | abs(phaseTransPlot) > 0.2) = 0;
%     phasePlotLogVals = log(1+log(1+phaseTransPlot));
% %     phaseTransPlot(phasePlotLogVals > mean(phasePlotLogVals(:)+std(phasePlotLogVals(:)))) = 0;
% %     pcolor(tSamples(segIdx), cwtfreqCell{k}*Fs, log(1+log(1+phasePlot)));
%     pcolor(tSamples(segIdx), (cwtfreqCell{k}(yAxisIdx))*Fs, phaseTransPlot);
% %     pcolor(tSamples(segIdx), (cwtfreqCell{k}(yAxisIdx))*Fs, cwtcfsPhasePlot);
% %     axis([6.5,11.5,6,14]);
    shading interp;
    colormap(1-gray);
    hold on
    plot(tSamples(segIdx), if1(segIdx), 'r', 'linewidth', 1);
    plot(tSamples(segIdx), if2(segIdx), 'b', 'linewidth', 1);
    title(sprintf('Morse Wavelet $k$=%d', k-1), 'Interpreter', 'latex', 'fontsize', 20);
%     print(sprintf('./result/zebra%03d', k-1), '-djpeg');
end

%% save all images
% mtarray = cat(3,sstmtCell{:});
% segIdx = 201:1400;
% figure('Position',[scrsz(1) scrsz(2) scrsz(3)/2 scrsz(4)/2]);
% for j=1:size(mtarray,3)
%     itvPSmt = abs(Fs*(sum(mtarray(:,segIdx,1:j),3) / j).^2);
%     itvPS_logscale = qclamp(log(1+itvPSmt), 0.001);
%     pcolor(tSamples(segIdx), instFreqTic, itvPS_logscale);
%     shading interp;
%     colormap(1-gray);
%     xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
%     ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
%     title(sprintf('SST-CWT-MT, MT=%d', j), 'Interpreter', 'latex', 'fontsize', 20);
%     print(sprintf('./result/MT%03d', j), '-djpeg');
% end

%% save results for each orthogonal wavelet
% segIdx = 201:1400;
% mtarray = cat(3,sstmtCell{:});
% figure('Position',[scrsz(1) scrsz(2) scrsz(3)/2 scrsz(4)/2]);
% for j=1:size(mtarray,3)
%     itvPSmt = abs(Fs*mtarray(:,segIdx,j).^2);
%     itvPS_logscale = qclamp(log(1+itvPSmt), 0.001);
%     pcolor(tSamples(segIdx), instFreqTic, itvPS_logscale);
%     shading interp;
%     colormap(1-gray);
%     xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
%     ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
%     title(sprintf('SST-CWT, Morse Wavelet k=%d', j-1), 'Interpreter', 'latex', 'fontsize', 20);
%     print(sprintf('./result/Morse%03d', j-1), '-djpeg');
% end

%% visualize and compare SST results
% segIdx = 201:1400;
% itvPS = abs(Fs*sstmtResult(:,segIdx)/2).^2;
% itvPS_logscale = qclamp(log(1+itvPS), 0.001);
% 
% hq = figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
% 
% pcolor(tSamples(segIdx), instFreqTic, itvPS_logscale);
% shading interp;
% colormap(1-gray);
% % axis([0 inf FreqBounds(1)*Fs FreqBounds(2)*Fs]);
% xlabel('Time (sec)', 'Interpreter', 'latex', 'fontsize', 20);
% ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'fontsize', 20);
% title('SST-CWT', 'Interpreter', 'latex', 'fontsize', 20);
% hold on
% plot(tSamples(segIdx), if1(segIdx), 'r', 'linewidth', 1);
% plot(tSamples(segIdx), if2(segIdx), 'b', 'linewidth', 1);

