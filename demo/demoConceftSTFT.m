%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% profiling implementations of SST: MATLAB built-in vs. original %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clc;
clearvars
path(pathdef);
addpath(genpath(strcat(pwd, '/utils')), path);
scrsz = get(groot,'ScreenSize');

Hz = 100;
time = (1/Hz:1/Hz:16)';
N = length(time);
MT = 100; %%% number of basis functions in multitaperring
J = 2; %%% number of basis functions in conceft

initstate(1);

%% build random signal using the model in Conceft paper
am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess');
am1 = 2 + am1 ./ max(abs(am1));
am2 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess');
am2 = 2 + am2 ./ max(abs(am2));
am1(1:500) = 0;
am2(end-600:end) = 0;

if1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess');
if1 = 10 + 6 * if1 ./ max(abs(if1));
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

%% Conceft
NPts = 377;
MOrd = 2;
HTS = 10; %%% half-time support
[h, Dh, ~] = hermf(NPts, MOrd, HTS);

%%%% tfr --- Short-Time Fourier Transform
%%%% tfrsqClean --- Synchrosqueezed Short-Time Fourier Transform (no noise)
%%%% tfrtic = linspace(0, 0.5, N/2)';
%%%% *_tic --- = linspace(lowFreq, highFreq, fLen)';
[STFTRep, tfrtic, sqSTFTRepClean, sqSTFT_tic] = sqstft(clean, 0, 1/5, 2e-4, 1, h(1,:)', Dh(1,:)');
[~, ~, ~, ConceftSTFTRepClean, ConceftSTFT_tic] = conceft_stft(clean, 0, 1/5, 2e-4, 1, NPts, J, HTS, MT);

deltaT = tfrtic(2)-tfrtic(1);
sqSTFT_itvPS = abs((sqSTFTRepClean./( 2 * deltaT))).^2;
Conceft_itvPS = abs((ConceftSTFTRepClean./( 2 * deltaT))).^2;

%% plot and compare results
LLcwt(1) = quantile(log(1+sqSTFT_itvPS(:)),0.002); 
HHcwt(1) = quantile(log(1+sqSTFT_itvPS(:)),0.998);
LLcwt(2) = quantile(log(1+Conceft_itvPS(:)),0.002); 
HHcwt(2) = quantile(log(Conceft_itvPS(:)),0.998);

segIdx = 201:1400;
hq = figure('Position',[1 scrsz(4)/2 scrsz(3)*2 scrsz(4)]);

subplot(2,1,1);
imagesc(time(segIdx)-time(segIdx(1)-1), sqSTFT_tic*Hz, log(1+abs(sqSTFT_itvPS)), [LLcwt(1) HHcwt(1)]); colormap(1-gray);
axis xy; set(gca,'fontsize', 20); axis([0 inf 0 20]);
set(gca, 'xtick', []); ylabel('Freq (Hz)');
title('synchrosqueezed clean signal', 'Interpreter', 'latex');

subplot(2,1,2);
imagesc(time(segIdx)-time(segIdx(1)-1), ConceftSTFT_tic*Hz, log(1+abs(Conceft_itvPS)), [LLcwt(2) HHcwt(2)]); colormap(1-gray);
axis xy; set(gca,'fontsize', 20); axis([0 inf 0 20]);
set(gca, 'xtick', []); ylabel('Freq (Hz)');
title('conceft clean signal', 'Interpreter', 'latex');

