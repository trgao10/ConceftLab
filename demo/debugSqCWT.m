%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% profiling implementations of SST: MATLAB built-in vs. original %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clc;
clearvars
path(pathdef);
addpath(genpath(strcat(pwd, '/utils')), path);
scrsz = get(groot,'ScreenSize');

Hz = 100;
tSamples = (1/Hz:1/Hz:16)';
N = length(tSamples);
MT = 100; %%% number of basis functions in multitaperring
J = 2; %%% number of basis functions in conceft
Smooth = 0;
Hemi = 0;
FreqBounds = [0,0.5];
FreqRes = 2e-4;
% FreqRes = diff(FreqBounds) / ((floor(log2(numel(tSamples)))-1)*48);

initstate(1);

%% build random signal using the model in Conceft paper
am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess');
am1 = 2 + am1 ./ max(abs(am1));
am2 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess');
am2 = 2 + am2 ./ max(abs(am2));
am1(1:500) = 0;
am2(end-600:end) = 0;

if1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess');
if1 = 25 + 6 * if1 ./ max(abs(if1));
if2 = smooth(cumsum(randn(N,1)) ./ Hz, 300, 'loess');
if2 = pi + 3 * if2 ./ max(abs(if2));
phi1 = cumsum(if1) / Hz;
phi2 = cumsum(if2) / Hz;

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

% initstate(1);

segIdx = 201:1400; %%% segment of the signal to be plotted

%% compare two different implementations of CWT
%%%% H.-T. Wu's original implementation
opts = struct('motherwavelet', 'morse', 'k', 0, 'beta', 30, 'gam', 1, 'dim', 8, 'rrnd', 0);
tic;[GauravWuCWTResult, GauravWuCWTTic] = CWT(tSamples, targetSig, opts);toc;

%%%% native cwt routine shippted with MATLAB
tic;[NativeCWTResult,NativeFreq] = cwt(targetSig, 'morse', Hz, 'TimeBandwidth', 100, 'NumOctaves', floor(log2(length(clean)))-1, 'VoicesPerOctave', 32);toc;

figure('Position',[1 scrsz(4)/2 scrsz(3)*2 scrsz(4)]);
subplot(1,2,1);
pcolor(tSamples, GauravWuCWTTic/Hz, abs(GauravWuCWTResult)');
shading interp
% imagesc(TSample, GauravWuCWTTic, abs(GauravWuCWTResult)');
% axis xy
title('scalogram - Gaurav-Wu', 'Interpreter', 'latex');
subplot(1,2,2);
pcolor(tSamples, NativeFreq, abs(NativeCWTResult));
shading interp
% imagesc(TSample, NativeFreq, abs(NativeCWTResult));
% axis xy
title('spectrogram - MATLAB native', 'Interpreter', 'latex');

% keyboard

[CWTResult, cwt_tic, sqCWTResult, sqcwt_tic] = sqCWTbase(tSamples, targetSig, FreqBounds(1)*Hz, FreqBounds(2)*Hz, FreqRes*Hz, opts, 0, 0);
DeltaT = cwt_tic(2)-cwt_tic(1);
% sqCWT_itvPS = abs(sqCWTResult);
sqCWT_itvPS = abs((sqCWTResult./(2 * DeltaT))).^2;
% sqCWT_itvPS = sqCWT_itvPS / sum(sqCWT_itvPS(:))*1e10;

% %% version under debugging
% [dCWTResult, dcwt_tic, dsqCWTResult, dsqcwt_tic] = sqCWTbase(tSamples, targetSig, FreqBounds(1)*Hz, FreqBounds(2)*Hz, FreqRes*Hz, opts, 0, 0);
% % DeltaT = dcwt_tic(2)-dcwt_tic(1);
% % dsqCWT_itvPS = abs(dsqCWTResult);
% dsqCWT_itvPS = abs((dsqCWTResult./(2 * DeltaT))).^2;
% % dsqCWT_itvPS = dsqCWT_itvPS / sum(dsqCWT_itvPS(:));

%% native wsst routine shippted with MATLAB
[NativeSSTResult,NativeSSTFreq] = wsst(targetSig, Hz, 'amor', 'VoicesPerOctave', 48);
% [NativeSSTResult,NativeSSTFreq] = wsst(clean, Hz, 'amor', 'NumOctaves', floor(log2(numel(clean)))-1, 'VoicesPerOctave', 48);
% [NativeSSTResult,NativeSSTFreq] = wsst(clean, Hz, 'amor');
% NativeSST_itvPS = abs(NativeSSTResult);
NativeSST_itvPS = abs((NativeSSTResult./(2 * DeltaT))).^2;
%%%% ad-hoc adjustment bringing equal the total sums of entries in 
%%%%% NativeSST_itvPS and in sqCWT_itvPS, just for visualization
NativeSST_itvPS = NativeSST_itvPS * (sum(sqCWT_itvPS(:))/sum(NativeSST_itvPS(:)));
% NativeSST_itvPS = NativeSST_itvPS / sum(NativeSST_itvPS(:));

%% plot and compare results
LLcwt(1) = quantile(log(1+sqCWT_itvPS(:)),0.002); 
HHcwt(1) = quantile(log(1+sqCWT_itvPS(:)),0.998);
% LLcwt(2) = quantile(log(1+dsqCWT_itvPS(:)),0.002); 
% HHcwt(2) = quantile(log(1+dsqCWT_itvPS(:)),0.998);
LLcwt(3) = quantile(log(1+NativeSST_itvPS(:)),0.002); 
HHcwt(3) = quantile(log(1+NativeSST_itvPS(:)),0.998);
% LLcwt(3) = quantile(log(1+NativeSST_itvPS(:)),0.2); 
% HHcwt(3) = quantile(log(1+NativeSST_itvPS(:)),0.8);

close(gcf);

hq = figure('Position',[1 scrsz(4)/2 scrsz(3)*2 scrsz(4)]);

% subplot(3,1,1);
% imagesc(TSample(segIdx)-TSample(segIdx(1)-1), sqcwt_tic, log(1+abs(sqCWT_itvPS)), [LLcwt(1) HHcwt(1)]); colormap(1-gray);
% axis xy; set(gca,'fontsize', 20); axis([0 inf FreqBounds(1)*Hz FreqBounds(2)*Hz]);
% subplot(3,1,1);
subplot(2,1,1);
pcolor(tSamples, sqcwt_tic, log(1+abs(sqCWT_itvPS)));
shading interp;
colormap(1-gray);
axis([0 inf FreqBounds(1)*Hz FreqBounds(2)*Hz]);
set(gca, 'xtick', []); ylabel('Freq (Hz)');
title('synchrosqueezed clean signal', 'Interpreter', 'latex');
hold on
plot(tSamples, if1, 'r', 'linewidth', 1) ;
plot(tSamples, if2, 'b', 'linewidth', 1) ;

% % subplot(3,1,2);
% % imagesc(TSample(segIdx)-TSample(segIdx(1)-1), dsqcwt_tic, log(1+abs(dsqCWT_itvPS)), [LLcwt(2) HHcwt(2)]); colormap(1-gray);
% % axis xy; set(gca,'fontsize', 20); axis([0 inf FreqBounds(1)*Hz FreqBounds(2)*Hz]);
% subplot(3,1,2);
% pcolor(tSamples, dsqcwt_tic, log(1+abs(dsqCWT_itvPS)));
% shading interp;
% colormap(1-gray);
% axis([0 inf FreqBounds(1)*Hz FreqBounds(2)*Hz]);
% set(gca, 'xtick', []); ylabel('Freq (Hz)');
% title('debugging: synchrosqueezed clean signal', 'Interpreter', 'latex');
% hold on
% plot(tSamples, if1, 'r', 'linewidth', 1) ;
% plot(tSamples, if2, 'b', 'linewidth', 1) ;

% subplot(3,1,3);
% imagesc(TSample(segIdx)-TSample(segIdx(1)-1), NativeSSTFreq, log(1+abs(NativeSST_itvPS)), [LLcwt(3) HHcwt(3)]); colormap(1-gray);
% axis xy; set(gca,'fontsize', 20); axis([0 inf FreqBounds(1)*Hz FreqBounds(2)*Hz]);
% subplot(3,1,3);
subplot(2,1,2);
pcolor(tSamples, NativeSSTFreq, log(1+abs(NativeSST_itvPS)));
% pcolor(tSamples, NativeSSTFreq, log(1+abs(NativeSST_itvPS)));
shading interp;
colormap(1-gray);
axis([0 inf FreqBounds(1)*Hz FreqBounds(2)*Hz]);
set(gca, 'xtick', []); ylabel('Freq (Hz)');
title('MATLAB wsst: synchrosqueezed clean signal', 'Interpreter', 'latex');
hold on
plot(tSamples, if1, 'r', 'linewidth', 1) ;
plot(tSamples, if2, 'b', 'linewidth', 1) ;

