clear all;
close all;
clc;

f1 = 20;
f2 = 50;
f3 = 1000;
f4 = 5000;

Fs = 20000;
Ts = 1/Fs;

t = 0:Ts:0.5;

pv = 2*pi.*t;

xt = 5*sin(f1.*pv) + cos(f2.*pv) + 15*cos(f3.*pv) + 10*cos(f4.*pv);

N = length(t);

XF = fft(xt,N);
XF = fftshift(XF);
XF = abs(XF)/N;

W = linspace(-Fs/2,Fs/2,N);

figure(1)

subplot(211)
plot(t,xt,'b');
title('Given Signal x(t)')
xlabel('time(s)')
ylabel('Amplitude')
legend('x(t)')

subplot(212)
plot(W,XF,'Color','black');
title('Frequency response of the given signal X(f)')
xlabel('frequency(Hz)')
ylabel('Amplitude (Normalized)')
legend('X(f)')


no1 = 5;
wc1 = f1/(Fs/2);
[b1 , a1] = butter(no1,wc1,'low');
h1 = freqz(b1,a1,Fs/2);

no2 = 80;
wc2 = f4/(Fs/2);
[b2 ,a2] = butter(no2,wc2,'high');
h2 = freqz(b2,a2,Fs/2);

no3 = 3;
wc3 = [f3-250 f3+250]/(Fs/2);
[b3 ,a3] = butter(no3,wc3,'bandpass');
h3 = freqz(b3,a3,Fs/2);

figure(2)
subplot(311)
plot(abs(h1),'color','r')
title('Low Pass Filter with cutoff f1')
ylabel('Amplitude')
xlabel('units')
legend('hlp')
xlim([0 100])

subplot(312)
plot(abs(h2),'r')
title('High Pass Filter with cutoff f4')
ylabel('Amplitude')
xlabel('units')
legend('hhp')


subplot(313)
plot(abs(h3),'r')
title('Bandpass Filter with cutoff f3-250 - f3+250')
ylabel('Amplitude')
xlabel('units')
legend('hbp')
xlim([0 5000])


y1t = filter(b1,a1,xt);
y2t = filter(b2,a2,xt);
y3t = filter(b3,a3,xt);

y1f = fft(y1t,N);
y1f = fftshift(y1f);
y1f = abs(y1f/N);

y2f = fft(y2t,N);
y2f = fftshift(y2f);
y2f = abs(y2f/N);

y3f = fft(y3t,N);
y3f = fftshift(y3f);
y3f = abs(y3f/N);

figure(3)
subplot(221)
plot(t,xt)
title('Original Signal x(t)')
xlabel('time(s)')
ylabel('Amplitude')
legend('x(t)')
xlim([0 0.5])

subplot(222)
plot(W,y1f,'g')
title("LP filtered frequency response of x(t) with f1")
xlabel('frequency(Hz)')
ylabel('Amplitude')
legend('X(f)')
xlim([-500 500])

subplot(223)
plot(W,y2f,'g')
title("HP filtered frequency response of x(t) with f4")
xlabel('frequency(Hz)')
ylabel('Amplitude')
legend('X(f)')



subplot(224)
plot(W,y3f,'g')
title("BP filtered frequency response of x(t) with f3")
xlabel('frequency(Hz)')
ylabel('Amplitude')
legend('X(f)')















