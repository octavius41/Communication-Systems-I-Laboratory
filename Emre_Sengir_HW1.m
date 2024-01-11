clc
close all
clear all

fs = 1;
T = 1/fs;


t1 = 0:T:5;
t2 = (5+T):T:(15-T);
t3 = 15:T:25;


t = [t1 t2 t3];

e = exp(1);
p = power(e,-t3/5);

xt1 = 2.*t1 - 5;
xt2 = 5*cos(2*pi.*t2);
xt3 = 100*p;

xt = [xt1 xt2 xt3];

x2t = cos(4*pi.*t);

y1t = xt.*x2t;

figure(1)
plot(t,xt)
title("x(t)")
legend('x(t) signal')
xlabel('time(s)')
ylabel('Amplitude')

figure(2)

plot(t,y1t)
title('output y(t)')
legend('y(t) signal')
xlabel('time(s)')
ylabel('Amplitude')


N = length(t);


XTF = fft(xt,N);
XTF = fftshift(XTF);
XTF = abs(XTF/N);

YTF = fft(y1t,N);
YTF = fftshift(YTF);
YTF = abs(YTF/N);

W = linspace(-fs/2,fs/2,N);

figure(3)

subplot(211)
plot(W,XTF)
title("Frequency Response of x(t), |X(w)|")
legend('X(w)')
xlabel("frequency (w)")
ylabel("Amplitude")

subplot(212)
plot(W,YTF)
title("Frequency Response of y(t), |Y(w)|")
legend('Y(w)')
xlabel("frequency (w)")
ylabel("Amplitude")


N = length(xt);
W = linspace(-fs/2,fs/2,N);

XTF = fft(xt,N);

X2TF = fft(x2t,N);

Y2TF = XTF.*X2TF;
Y2TFt = fftshift(Y2TF);
Y2TFt = abs(Y2TFt/N);



figure(4)

plot(W,Y2TFt)
title('Frequency spectrum of output |Y2(w)|')
legend('Y2(w)')
xlabel("Frequency (w)")
ylabel("Amplitude")

y2t = ifft(Y2TF,N);
tp = 0:T:51+T;

figure(5)
plot(t,y2t)
title('Output y2(t) obtained from multiplication on w domain')
legend('y2(t)')
xlabel("time(s)")
ylabel("Amplitude")
xlim([0 50])

y2tt = conv(xt,x2t);



figure(6)
plot(t,y2tt)
title('Output y2(t) obtained from convulation in time domain')
legend('y2(t)')
xlabel("time(s)")
ylabel("Amplitude")
xlim([0 50])








