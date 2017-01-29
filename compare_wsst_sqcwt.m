%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% profiling implementations of SST: MATLAB built-in vs. original %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clc;
clearvars
path(pathdef);
addpath(genpath(strcat(pwd, '/utils')), path);
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% load chirp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load quadchirp;
% Fs = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% construct a linear chirp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 1000;
T = 4; %%% signal duration time
tSamples = 1/Fs:1/Fs:T;
% tSamples = 0:1/Fs:T;
omega_0 = 400;
% instFreq = omega_0+tSamples*(Fs/2-omega_0)/T;
% instFreq = omega_0+tSamples*(Fs/2-omega_0)/T;
instFreq = omega_0*ones(size(tSamples));
AmT = ones(size(tSamples));
% AmT = (TSamples-mean(TSamples)).^2;
sigClean = AmT.*cos(2*pi*cumsum(instFreq)/Fs);

% plot(tSamples,sigClean);

% %%%% compare this with a chirp generated using MATLAB routines
% NativeLinearChirp = chirp(tSamples, omega_0, tSamples(end), instFreq(end));
% figure;
% plot(tSamples,NativeLinearChirp);
% keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MATLAB wsst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[NativeSSTResult,NativeSSTFreq] = wsstgao(sigClean, Fs, 'morse',...
    'VoicesPerOctave', 64, 'ExtendSignal', false,...
    'WaveletParameters', struct('be',80,'ga',1));
% [NativeSSTResult,NativeSSTFreq] = wsst(sigClean, Fs, 'amor', 'VoicesPerOctave', 48);
toc;

% sqCWT_itvPS = abs(sqCWTResult);
Native_itvPS = abs(Fs*NativeSSTResult/2).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Conceft wsst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = struct('motherwavelet', 'morse', 'k', 0, 'beta', 80, 'gam', 1);
FreqRes = (Fs/4)/(length(NativeSSTFreq)-1);

tic;
[CWTResult, cwt_tic, sqCWTResult, sqcwt_tic] = sqCWTbase(tSamples', sigClean', 0, Fs/2, FreqRes, opts, 0, 1);
% [CWTResult, cwt_tic, sqCWTResult, sqcwt_tic] = sqCWTbase(tSamples(2:end)', sigClean(2:end)', 0, Fs/2, FreqRes, opts, 0, 0);
toc;

Conceft_itvPS = abs(Fs*sqCWTResult/2).^2;
Conceft_itvPS = Conceft_itvPS * (sum(Native_itvPS)/sum(Conceft_itvPS));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% plot and compare sst results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [NativeQLow, NativeQHigh] = qclamp(abs(NativeSSTResult).^2, 0.002, 0.998);

% figure;
% subplot(1,2,1);
% imagesc(tSamples, NativeSSTFreq, log(1+abs(NativeSSTResult).^2), [NativeQLow, NativeQHigh]);
% % colormap(1-gray);
% axis xy;
% axis([0 inf NativeSSTFreq(1) NativeSSTFreq(end)]);
% % set(gca, 'xtick', []); ylabel('Freq (Hz)');
% % title('', 'Interpreter', 'latex');

% keyboard

figure;
subplot(1,2,1);
% pcolor(tSamples, NativeSSTFreq, abs(NativeSSTResult));
pcolor(tSamples, NativeSSTFreq, log(1+abs(NativeSSTResult)));
% pcolor(tSamples, NativeSSTFreq, log(1+abs(NativeSSTResult).^2));
% imagesc(tSamples, NativeSSTFreq, log(1+abs(NativeSSTResult).^2));
% axis xy
shading interp
colormap(1-gray);
% hold on
% plot(tSamples, instFreq, 'r', 'linewidth', 1);
% xlabel('Seconds'); ylabel('Frequency');

subplot(1,2,2);
pcolor(tSamples, sqcwt_tic, log(1+abs(sqCWTResult)));
% pcolor(tSamples, sqcwt_tic, abs(sqCWTResult));
% pcolor(tSamples(2:end), sqcwt_tic, abs(sqCWTResult));
% imagesc(tSamples, sqcwt_tic, log(1+abs(sqCWTResult).^2));
% axis xy
shading interp
colormap(1-gray);
% hold on
% plot(tSamples, instFreq, 'r', 'linewidth', 1);
% xlabel('Seconds'); ylabel('Frequency');

