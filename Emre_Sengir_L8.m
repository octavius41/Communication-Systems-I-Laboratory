clear all;close all;clc;

%% creation
img = imread("testpat1.png");
y = im2double(img);
Fs = numel(y);
Fsr = size(y,1);
dt = 1/Fs;
t = 0:dt:((Fs-1)/Fs);
fc = 10000;

mt = reshape(y,1,Fs);
w = linspace(-Fs/2,Fs/2,Fs);
N = Fs;

Mf = fft(mt,N);
Mf = fftshift(Mf);
Mf = abs(Mf)/N;

figure(1)
plot(w,Mf); title("magnitude response of message"); xlabel("frequency"); ylabel("magnitude"); legend("|Mf|");xlim([-20000 20000])

%% modulation

kf1 = 100;
kf2 =400;

mkf1 = 2*pi*kf1*cumsum(mt)*dt;
mkf2 = 2*pi*kf2*cumsum(mt)*dt;

s1t = cos(2*pi*fc.*t + mkf1);
s2t = cos(2*pi*fc.*t + mkf2);

S1f = fft(s1t,N);
S1f = fftshift(S1f);
S1f = abs(S1f)/N;

S2f = fft(s2t,N);
S2f = fftshift(S2f);
S2f = abs(S2f)/N;

figure(2)
subplot(211)
plot(w,S1f);title("magnitude response of modulated signal w.r.t. kf = 100"); xlabel("frequency"); ylabel("magnitude"); legend("|S1f|");xlim([-20000 20000])

subplot(212)
plot(w,S2f);title("magnitude response of modulated signal w.r.t. kf = 400"); xlabel("frequency"); ylabel("magnitude"); legend("|S2f|");xlim([-20000 20000])

%% noise_imp

s1t0 = awgn(s1t,0,'measured');
s2t0 = awgn(s2t,0,'measured');

s1t25 = awgn(s1t,25,'measured');
s2t25 = awgn(s2t,25,'measured');

s1t50 = awgn(s1t,50,'measured');
s2t50 = awgn(s2t,50,'measured');

S1f0 = fft(s1t0,N);
S1f0 = fftshift(S1f0);
S1f0 = abs(S1f0)/N;

S2f0 = fft(s2t0,N);
S2f0 = fftshift(S2f0);
S2f0 = abs(S2f0)/N;

figure(3)
subplot(211)
plot(w,S1f0);title("magnitude response of noise implemented signal w.r.t. kf = 100 and snr=0"); xlabel("frequency"); ylabel("magnitude"); legend("|S1f|");xlim([-20000 20000])

subplot(212)
plot(w,S2f0);title("magnitude response of noise implemented signal w.r.t. kf = 400 snr=0"); xlabel("frequency"); ylabel("magnitude"); legend("|S2f|");xlim([-20000 20000])

%% demod

ds1t0 = fmdemod(s1t0,fc,Fs,kf1);
ds2t0 = fmdemod(s2t0,fc,Fs,kf2);

ds1t25 = fmdemod(s1t25,fc,Fs,kf1);
ds2t25 = fmdemod(s2t25,fc,Fs,kf2);

ds1t50 = fmdemod(s1t50,fc,Fs,kf1);
ds2t50 = fmdemod(s2t50,fc,Fs,kf2);

DS1f25 = fft(ds1t25,N);
DS1f25 = fftshift(DS1f25);
DS1f25 = abs(DS1f25)/N;

DS2f25 = fft(ds2t25,N);
DS2f25 = fftshift(DS2f25);
DS2f25 = abs(DS2f25)/N;

figure(4)
subplot(211)
plot(w,DS1f25);title("magnitude response of demodulated signal w.r.t. kf = 100 and snr=25"); xlabel("frequency"); ylabel("magnitude"); legend("|S1f|");xlim([-20000 20000])

subplot(212)
plot(w,DS2f25);title("magnitude response of demodulated signal w.r.t. kf = 400 snr=25"); xlabel("frequency"); ylabel("magnitude"); legend("|S2f|");xlim([-20000 20000])

ps1nr0 = psnr(ds1t0,mt);
ps2nr0 = psnr(ds2t0,mt);

ps1nr25 = psnr(ds1t25,mt);
ps2nr25 = psnr(ds2t25,mt);

[b,a] = butter(5,(fc+5000)/(Fs/2),'low');
f1d50 = filter(b,a,ds1t50);
f2d50 = filter(b,a,ds2t50);

ps1nr50 = psnr(f1d50,mt);
ps2nr50 = psnr(f2d50,mt);

r1p0 = reshape(ds1t0,Fsr,Fsr);
r2p0 = reshape(ds2t0,Fsr,Fsr);

r1p25 = reshape(ds1t25,Fsr,Fsr);
r2p25 = reshape(ds2t25,Fsr,Fsr);

r1p50 = reshape(f1d50,Fsr,Fsr);
r2p50 = reshape(f2d50,Fsr,Fsr);

figure(5)
subplot(231)
imshow(r1p0); title("demod img of kf=100 snr=0")

subplot(232)
imshow(r2p0); title("demod img of kf=100 snr=25")

subplot(233)
imshow(r1p25); title("demod img of kf=100 snr=50")

subplot(234)
imshow(r2p25); title("demod img of kf=400 snr=0")

subplot(235)
imshow(r1p50); title("demod img of kf=400 snr=25")

subplot(236)
imshow(r2p50); title("demod img of kf=400 snr=50")

snrv = [0 25 50];
psnrv1 = [ps1nr0 ps1nr25 ps1nr50];
psnrv2 = [ps2nr0 ps2nr25 ps2nr50];

figure(6)
plot(snrv,psnrv1); 
hold on;
plot(snrv,psnrv2); title("snr-psnr relation");xlabel("snr values");ylabel("psnr values");legend("kf1 = 100","kf2 = 400");
hold off;

ms1e0 = immse(r1p0,y);
ms2e0 = immse(r2p0,y);

ms1e25 = immse(r1p25,y);
ms2e25 = immse(r2p25,y);

ms1e50 = immse(r1p50,y);
ms2e50 = immse(r2p50,y);

mse1v = [ms1e0 ms1e25 ms1e50];
mse2v = [ms2e0 ms2e25 ms2e50];

figure(7)
plot(snrv,mse1v);
hold on;
plot(snrv,mse2v);title("mse-snr relation");xlabel("snr values");ylabel("mse values");legend("kf1 = 100","kf2 = 400");
hold off;

MAX = max(mt);











