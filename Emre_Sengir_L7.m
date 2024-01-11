clc; clear all; close all;

%% creation
img = imread('cameraman.tif');
y = im2double(img);
Fs = numel(y);
Fsr = size(y,1);
dt = 1/Fs;
t = 0:dt:(Fs-1)/Fs;
fc = 20000;

mt = reshape(y,1,Fs);
ct = cos(2*pi*fc.*t);
mm = max(mt);

%% mod

ka = 0.9;
st = (1+ka.*mt).*ct;

%% noise_imp

snrv = [0 10 30];

Pm = sum(abs(mt).^2)/Fs;

snrlin1 = power(10,0.1*snrv(1));
snrlin2 = power(10,0.1*snrv(2));
snrlin3 = power(10,0.1*snrv(3));

var1 = Pm/snrlin1;
var2 = Pm/snrlin2;
var3 = Pm/snrlin3;

var = [var1 var2 var3];

n1 = sqrt(var(1)).*randn(1,Fs);
n2 = sqrt(var(2)).*randn(1,Fs);
n3 = sqrt(var(3)).*randn(1,Fs);

r1 = st + n1;
r2 = st + n2;
r3 = st + n3;

%% mag_resp

Rf = fft(r2,Fs);
Rf = fftshift(Rf);
Rf = abs(Rf)/Fs;

w = linspace(-Fs/2,Fs/2,Fs);


%% demod

dmd1 = r1.^2;
dmd2 = r2.^2;
dmd3 = r3.^2;

[b,a] = butter(3,(4*fc/5)/(Fs/2),"low");

fdmd1 = filter(b,a,dmd1);
fdmd2 = filter(b,a,dmd2);
fdmd3 = filter(b,a,dmd3);

sdmd1 = sqrt(fdmd1);
sdmd2 = sqrt(fdmd2);
sdmd3 = sqrt(fdmd3);

rp1 = (sdmd1 - 1)/ka;
rp2 = (sdmd2 - 1)/ka;
rp3 = (sdmd3 - 1)/ka;

rsp1 = reshape(rp1,Fsr,Fsr);
rsp2 = reshape(rp2,Fsr,Fsr);
rsp3 = reshape(rp3,Fsr,Fsr);

figure(1)
subplot(211)
plot(t,rp3);
subplot(212)
plot(t,mt)

figure(2)

subplot(221)
imshow(y); title('original');

subplot(222)
imshow(rsp1); title('SNR = 0');

subplot(223)
imshow(rsp2); title('SNR = 10');

subplot(224)
imshow(rsp3); title('SNR = 30');


MSE1 = 1/Fs*(sum(sum((y-rsp1).^2),2));
MSE2 = 1/Fs*(sum(sum((y-rsp2).^2),2));
MSE3 = 1/Fs*(sum(sum((y-rsp3).^2),2));

MSE = [MSE1 MSE2 MSE3];

figure(3)
plot(snrv, MSE); xlabel('SNR');ylabel('MSE');title('mean-square error according to SNR ratio');legend('trace')











