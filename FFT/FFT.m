clc
T=400; N=256; df=1/T; dt=T/N; Fm=1/dt; t=[0:dt:T-dt]; funshifted=[0:df:Fm-df]-Fm/2; %geradzahlig
RecSamVec=[0:dt/5:2*T]; f0=1/25; f=[-Fm/2:df:Fm/2-df] 

fodd=1/T*[-(N-1)/2:(N-1)/2];
feven=1/T*[-N/2:N/2-1];
if mod(N,2)==1
    f=fodd;
else
    f=feven;
end

fun = @(x) exp(-((x-50)/20).^2).*sin(2*pi*f0*x);
%fun = @(x) 0.5*(sign(x-50)-sign(x-70));

Ff=fftshift(ifft(fun(t)));

figure(1)
plot(f,real(Ff),f,imag(Ff))
figure(4)
semilogy(f,abs(Ff))

W=exp(1i*RecSamVec'*f*2*pi);
g=W*Ff';


figure(2)
plot(RecSamVec,fun(RecSamVec),RecSamVec,g)
figure(3)
plot(RecSamVec,fun(RecSamVec))



