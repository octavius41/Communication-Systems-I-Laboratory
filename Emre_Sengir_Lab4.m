clear all;
clc;
close all;
%% creation
Am = 2;
fm = 50;
Ac = 1;
fc = 1000;
phasec = pi/6;
Fs = 50000;
T = 1/Fs;
d = 0.1;
t=0:T:d;

mt = Am*cos(2*pi*fm.*t);
ct = Ac*cos(2*pi*fc.*t + phasec);
%% mod
st = mt.*ct;

figure(1)
plot(t,mt,'r'); 

hold on;

plot(t,ct,'g'); title("message signal and carrier signal"); ylabel('amplitude'); xlabel('time(s)'); legend('m(t)','c(t)');

hold off;

figure(2)
plot(t,st); title("dsb-sc modulated signal"); ylabel('amplitude'); xlabel('time(s)'); legend('s(t)')
%% demod
Acp = 2;
phasecp = 0;
lol = Acp*cos(2*pi*fc.*t + phasecp);

u = st.*lol;

nl = 7;
wcl = fc/(Fs/2);
[bl , al] = butter(nl,wcl,'low');

y = filter(bl,al,u);

N = length(t);
w = linspace(-Fs/2,Fs/2,N);

Mf = fft(mt,N);
Mf = fftshift(Mf);
Mf = abs(Mf/N);

Sf = fft(st,N);
Sf = fftshift(Sf);
Sf = abs(Sf/N);

Uf = fft(u,N);
Uf = fftshift(Uf);
Uf = abs(Uf/N);

figure(3);
subplot(311)
plot(w,Mf); title("magnitude response of m(t)"); ylabel('amplitude'); xlabel('frequency'); legend('M(f)')
xlim([-200 200])
subplot(312)
plot(w,Sf); title("magnitude response of s(t)"); ylabel('amplitude'); xlabel('frequency'); legend('S(f)')
xlim([-2000 2000])
subplot(313)
plot(w,Uf); title("magnitude response of u(t)"); ylabel('amplitude'); xlabel('frequency'); legend('U(f)')
xlim([-5000 5000])

figure(4)
plot(t,mt);
hold on;
plot(t,y); title("original message signal and demodulated equivalent signal"); ylabel('amplitude'); xlabel('time(s)'); legend('m(t)',"m'(t)");
hold off;

%% ssb

nbp = 5;
wcbp = [1000 2000]/(Fs/2);
[bbp ,abp] = butter(nbp,wcbp,'bandpass');
hbp = freqz(bbp,abp,Fs/2);

sut = filter(bbp,abp,st);



Suf = fft(sut,N);
Suf = fftshift(Suf);
Suf = abs(Suf/N);

figure(5)
subplot(211)
plot(abs(hbp)); title("magnitude response of ssb filter"); ylabel('amplitude'); xlabel('frequency'); legend('H(f)')
xlim([0 5000])
subplot(212)
plot(w,Suf); title("magnitude response of ssb filtered modulated signal"); ylabel('amplitude'); xlabel('frequency'); legend('Su(f)')
xlim([-5000 5000])




%% m(t).cos(2π.fc.t)    u(t) = m(t).cos(2π.fc.t).(Am/Ac).cos(2π.fc.t+θ)     m(t).[(1/2.(Am/Ac)).(cos(2a+θ) + cos(θ)) for θc = 0;

