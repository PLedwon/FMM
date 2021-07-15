clc

Fs = 1000; T=1/Fs; L=64; t=(0:L-1)*T; f=Fs*(0:(L/2))/L; f0=25; sig=.05;
T = sin(2*pi*f0*t); 
F = fft(T); FS=fftshift(F);
P2 = abs(FS/L);
P1=P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1);

figure(1)
plot(Fs*[1:L]/L-500,log10(P2))

plot(t,T,'.-')




