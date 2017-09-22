%% Example 0: square bump function
L = 10;
N = 50;

tWinFunc = zeros(1,N);
tWinFunc(1:L) = 1;

fWinFunc = fft(tWinFunc,N);
omega = 2*pi*(0:N-1)/N;
trueFWinFunc = sin(omega*L/2)./sin(omega/2).*exp(-1i*omega/2*(L-1));

figure;

subplot(1,3,1);
plot(0:(N-1),tWinFunc);

subplot(1,3,2);
plot(0:(N-1),abs(fWinFunc),'+');
hold on
plot(0:(N-1),abs(fWinFunc),'b-');
plot(0:(N-1),abs(trueFWinFunc),'o');

subplot(1,3,3);
plot(0:(N-1),angle(fWinFunc),'+');
hold on
plot(0:(N-1),angle(fWinFunc),'b-');

angleGroundTruth = angle(trueFWinFunc);
angleGroundTruth(1) = 0;
angleGroundTruth(6:5:length(trueFWinFunc)) = angleGroundTruth(6:5:length(trueFWinFunc))+pi;
plot(0:(N-1),angleGroundTruth,'o-');

%%% the confusion about phases
figure;

subplot(1,2,1);
plot(0:(N-1),angleGroundTruth,'+');
hold on
plot(0:(N-1),angle(exp(-1i*omega/2*(L-1))),'o');

subplot(1,2,2);
plot(0:(N-1),angle(exp(-1i*omega/2*(L-1)))-angleGroundTruth,'*');
hold on
angleFromSincPart = angle(sin(omega*L/2)./sin(omega/2));
angleFromSincPart(1) = 0;
plot(0:(N-1),angleFromSincPart,'o');
legend('angle discrepancy','angle contributed from SINC');
% angleGroundTruth2 = angle(exp(-1i*omega/2*(L-1)));
% plot(0:(N-1),angle(exp(-1i*omega*(L-1)/2)),'kx');
% plot(0:(N-1),rem(-omega*(L-1)/2+20*pi,2*pi),'kx');

%% Example 1: standard Gaussian
%%%% assume our 

gau = @(x) exp(-x.^2/2);
% gau = @(x) ones(size(x));
% gau = @(x) cos(x);
Fs = 20;
T = 1;
ts = -T/2:1/Fs:(T/2-1/Fs);  %%% this way we make sure the finite sample of the
                        %%% signal is of even length

tWinFunc = gau(ts);
mask = zeros(size(tWinFunc));
mask((length(tWinFunc)/2+1):end) = 1;
tWinFunc = tWinFunc.*mask;
%%% ifft is precisely the inverse of fft!
% figure;
% plot(ifft(fft(tWinFunc)),'bo')
% hold on
% plot(tWinFunc,'r+')

%%% view the result of fft directly
figure;
subplot(1,3,1);
% plot(ts,tWinFunc);
plot(fftshift(tWinFunc));
subplot(1,3,2);
% plot(abs(fftshift(fft(tWinFunc))));
plot(abs(fft(tWinFunc)));
subplot(1,3,3);
plot(angle(fftshift(fft(tWinFunc))));
% plot(abs(fftshift(fft(tWinFunc))/Fs));
% hold on
% plot(imag(fftshift(fft(tWinFunc))/Fs));

%%
DFTFreqTicks = 2*pi*((1:length(ts))-1)/length(ts)*Fs;
% fWinFunc = abs(fftshift(fft(tWinFunc)*0.1));

figure;
plot(fftshift(real(fft(tWinFunc)/Fs)));
hold on
plot(fftshift(gau(DFTFreqTicks)));