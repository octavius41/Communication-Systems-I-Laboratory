clear all;
close all;
clc;

%% creation

A1 = 1;
A2 = 2;
Ac = 2;
f1 = 10;
f2 = 15;
fc = 500;
Fs = 4000;
T = 1/Fs;
t = 0:T:0.4;

mt =A1*sin(2*pi*f1.*t) + A2*sin(2*pi*f2.*t);
ct = Ac*sin(2*pi*fc.*t);

figure(1)
plot(t,mt); title("message signal"); ylabel('amplitude'); xlabel('time'); legend('m(t)')

%% modulation

ka1 = 0.2;
ka2 = 0.6;

ka3 = 0.35;

s1t = (1 + ka1.*mt).*ct;
s2t = (1 + ka2.*mt).*ct;

s3t = (1 + ka3.*mt).*ct;

figure(2)
subplot(311)
plot(t,s1t); title("ka = 0.2 modulated signal"); ylabel('amplitude'); xlabel('time'); legend('s1(t)')

subplot(312)
plot(t,s2t); title("ka = 0.6 modulated signal"); ylabel('amplitude'); xlabel('time'); legend('s2(t)')

subplot(313)
plot(t,s3t); title("ka = 0.35 modulated signal"); ylabel('amplitude'); xlabel('time'); legend('s3(t)')

N = length(t);
w = linspace(-Fs/2,Fs/2,N);

S1f = fft(s1t,N);
S1f = fftshift(S1f);
S1f = abs(S1f/N);

S2f = fft(s2t,N);
S2f = fftshift(S2f);
S2f = abs(S2f/N);

S3f = fft(s3t,N);
S3f = fftshift(S3f);
S3f = abs(S3f/N);

figure(3)
subplot(311)
plot(w,S1f); title("magnitude response of s1(t)"); ylabel('amplitude'); xlabel('frequency'); legend('S1(f)')
xlim([-600 600])

subplot(312)
plot(w,S2f); title("magnitude response of s2(t)"); ylabel('amplitude'); xlabel('frequency'); legend('S2(f)')
xlim([-600 600])

subplot(313)
plot(w,S3f); title("magnitude response of s3(t)"); ylabel('amplitude'); xlabel('frequency'); legend('S3(f)')
xlim([-600 600])

%% demod

s1td = s1t.^2;
s2td = s2t.^2;

s3td = s3t.^2;

nl = 5;
wcl = fc/(Fs/2);
[bl , al] = butter(nl,wcl,'low');
hl = freqz(bl,al,Fs/2);

s1d = filter(bl,al,s1td);
s2d = filter(bl,al,s2td);

s3d = filter(bl,al,s3td);

S1df = fft(s1d,N);
S1df = fftshift(S1df);
S1df = abs(S1df/N);

S2df = fft(s2d,N);
S2df = fftshift(S2df);
S2df = abs(S2df/N);

S3df = fft(s3d,N);
S3df = fftshift(S3df);
S3df = abs(S3df/N);

figure(4)
subplot(311)
plot(t,s1d)


subplot(312)
plot(t,s2d)


subplot(313)
plot(t,s3d)



m1d = s1d.^(1/2);
m2d = s2d.^(1/2);

m3d = s3d.^(1/2);

m1ed = 4*(m1d - 1.5);
m2ed = 2*(m2d - 2);

m3ed = 2*(m3d - 1.4);

figure(5)
subplot(311)
plot(t,m1d)


subplot(312)
plot(t,m2d)


subplot(313)
plot(t,m3d)

figure(6)
subplot(311)
plot(t,m1ed,'r');
hold on;
plot(t,mt,'g'); title("demodulated signal for ka=0.2 and message signal"); ylabel('amplitude'); xlabel('time'); legend('m1ed(t)','m(t)')
hold off;

subplot(312)
plot(t,m2ed,'r');
hold on;
plot(t,mt,'g'); title("demodulated signal for ka=0.6 and message signal"); ylabel('amplitude'); xlabel('time'); legend('m2ed(t)','m(t)')
hold off;

subplot(313)
plot(t,m3ed,'r'); 
hold on;
plot(t,mt,'g'); title("demodulated signal for ka=0.35 and message signal"); ylabel('amplitude'); xlabel('time'); legend('m3ed(t)','m(t)')
hold off;


%%assuming that peak amplitude of the message signal was 2.85, best value for ka to obtain mü = 1 is ka = 0.35 but that causing-
%%distortion in the output wave so by trying it out I found that ka = 0.18














