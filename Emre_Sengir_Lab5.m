clear all;close all;clc;

%% creation
fm = 50;
Am = 1;
Fs = 2000;
T = 1/Fs;
d = (8*fm*T);
fc = 200;
Ac = 1;

t = 0:T:d;

mt = Am*sawtooth(2*pi*fm.*t,0.5);

ct = Ac*cos(2*pi*fc.*t);

figure(1)
subplot(211)
plot(t,mt); title("message signal m(t)"); ylabel('amplitude'); xlabel('time(s)'); legend('m(t)')

subplot(212)
plot(t,ct); title("carrier signal c(t)"); ylabel('amplitude'); xlabel('time(s)'); legend('c(t)')

%% mod

kf1 = 25;
kf2 = 75;

int1 = 2*pi*kf1*cumsum(mt)*T;
int2 = 2*pi*kf2*cumsum(mt)*T;

st1 = cos(2*pi*fc.*t + int1);
st2 = cos(2*pi*fc.*t + int2);

figure(2)
subplot(211)
plot(t,st1); title("kf1 B = 0.5B modulated signal"); ylabel('amplitude'); xlabel('time(s)'); legend('s1(t)')
xlim([0.06 0.14])
subplot(212)
plot(t,st2); title("kf1 B = 1.5B modulated signal"); ylabel('amplitude'); xlabel('time(s)'); legend('s2(t)')
xlim([0.06 0.14])
%% mag_resp

N = length(t);
w = linspace(-Fs/2,Fs/2,N);

Mf = fft(mt,N);
Mf = fftshift(Mf);
Mf = abs(Mf/N);

Cf = fft(ct,N);
Cf = fftshift(Cf);
Cf = abs(Cf/N);

Sf1 = fft(st1,N);
Sf1 = fftshift(Sf1);
Sf1 = abs(Sf1/N);

Sf2 = fft(st2,N);
Sf2 = fftshift(Sf2);
Sf2 = abs(Sf2/N);

figure(3)
subplot(211)
plot(w,Mf); title("magnitude response of m(t)"); ylabel('amplitude'); xlabel('frequency'); legend('|M(f)|')

subplot(212)
plot(w,Cf); title("magnitude response of c(t)"); ylabel('amplitude'); xlabel('frequency'); legend('|C(f)|')


figure(4)
subplot(211)
plot(w,Sf1); title("magnitude response of s1(t)"); ylabel('amplitude'); xlabel('frequency'); legend('|S1(f)|')
xlim([-300 300])

subplot(212)
plot(w,Sf2); title("magnitude response of s2(t)"); ylabel('amplitude'); xlabel('frequency'); legend('|S2(f)|')
xlim([-800 800])

%% demod

ds1t = fmdemod(st1,fc,Fs,kf1);
ds2t = fmdemod(st2,fc,Fs,kf2);

figure(5)
subplot(211)
plot(t,ds1t);
hold on;
plot(t,mt); title("kf1 B = 0.5B demodulated signal and message signal"); ylabel('amplitude'); xlabel('time(s)'); legend('s1(t)','m(t)')
ylim([-1 1])
hold off;

subplot(212)
plot(t,ds2t);
hold on;
plot(t,mt);  title("kf2 B = 1.5B demodulated signal and message signal"); ylabel('amplitude'); xlabel('time(s)'); legend('s2(t)','m(t)')
ylim([-1 1])
hold off;





