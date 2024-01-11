clear all; close all; clc;

%% create
Fs = 2000;
d = 0.2;
fc = 200;
f1 = 50;
f2 = 10;
Ac = 1;
T = 1/Fs;
t = 0:T:d;
ct = Ac*cos(2*pi*fc.*t);
m1t = cos(2*pi*f1.*t);
m2t = cos(2*pi*f2.*t);
m3t = m1t + m2t;

kp1 = pi/4;
kp2 = pi/2;
kp3 = (3*pi)/4;
%% mod
st = Ac*cos(2*pi*fc.*t + kp1.*m1t);
s2t = Ac*cos(2*pi*fc.*t + kp1.*m2t);
s3t = Ac*cos(2*pi*fc.*t + kp1.*m3t);

figure(1)
subplot(211)
plot(t,s3t); title("modulated signal for m3(t) = m1(t) + m2(t)"); ylabel('amplitude'); xlabel('time(s)'); legend('s3(t)')

subplot(212);
plot(t,(st + s2t)); title("addition of seperate modulated signals m1(t) and m2(t)"); ylabel('amplitude'); xlabel('time(s)'); legend('s1(t)+s2(t)')

%%they are not the same because angle modulation is a nonlinear modulation

%% kf_diff

spm1t = Ac*cos(2*pi*fc.*t + kp2.*m1t);
spm2t = Ac*cos(2*pi*fc.*t + kp3.*m1t);

figure(2)
subplot(411)
plot(t,m1t); title("message signal m1(t)"); ylabel('amplitude'); xlabel('time(s)'); legend('m1(t)')

subplot(412)
plot(t,st); title("pm modulated m1(t) signal for kp=pi/4"); ylabel('amplitude'); xlabel('time(s)'); legend('s1(t)')

subplot(413)
plot(t,spm1t);title("pm modulated m1(t) signal for kp=pi/2"); ylabel('amplitude'); xlabel('time(s)'); legend('s2(t)')

subplot(414)
plot(t,spm2t);title("pm modulated m1(t) signal for kp=3pi/4"); ylabel('amplitude'); xlabel('time(s)'); legend('s3(t)')

%% demod

spmh = hilbert(st);
Oi1 = atan(spmh/st);
sd1t = unwrap(angle(spmh),Oi1);
ds1t = (sd1t - 2*pi*fc*t)/kp1;

spm1h = hilbert(spm1t);
Oi2 = atan(spm1h/spm1t);
sd2t = unwrap(angle(spm1h),Oi2);
ds2t = (sd2t - 2*pi*fc*t)/kp2;

spm2h = hilbert(spm2t);
Oi3 = atan(spm2h/st);
sd3t = unwrap(angle(spm2h),Oi3);
ds3t = (sd3t - 2*pi*fc*t)/kp3;

N = length(t);
S1TF = fft(ds1t,N);
S1TF = fftshift(S1TF);
S1TF = abs(S1TF/N);
w = linspace(-Fs/2,Fs/2,N);

figure(3)
subplot(311)
plot(t,m1t);
hold on;
plot(t,ds1t,'.'); title("message signal and output demodulated signal for kp1 = pi/4"); ylabel('amplitude'); xlabel('time(s)'); legend('m1(t)','spmd1(t)')
ylim([-1 1])
hold off;

subplot(312)
plot(t,m1t);
hold on;
plot(t,ds2t,'.'); title("message signal and output demodulated signal for kp2 = pi/2"); ylabel('amplitude'); xlabel('time(s)'); legend('m1(t)','spmd2(t)')
ylim([-1 1])
hold off;

subplot(313)
plot(t,m1t);
hold on;
plot(t,ds3t,'.'); title("message signal and output demodulated signal for kp3 = 3pi/4"); ylabel('amplitude'); xlabel('time(s)'); legend('m1(t)','spmd3(t)')
ylim([-1 1])
hold off;








