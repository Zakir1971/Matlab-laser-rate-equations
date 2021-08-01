unction MML2
clc;
clear all
m=5;% Number of modes 10
n=m+1;
tspan=linspace(0,5e-9, 10000);%[0  1e-8];
y0=1e-6*ones(1,n);
[T,Y]= ode45(@(t,y) rate_equation2(t,y,n),tspan,y0);
figure;
for s=2:n
    plot(T,Y(:,s))
hold on
end
hold off
figure;
n1=length(T);
Ts=T(n1)/(n1-1); 
fs=1/Ts;
f=-fs/2:fs/(n1-1):fs/2;
f=f/fs;
c=3e8; lambda=980e-9;Vp=3e-12;
F=1; vg=1e10; h=6.63e-34; f1=c/lambda;
Power=F*vg*h*f1*Y(:,1)*Vp*1e3;
y=fftshift(fft(Power));
z=log(abs(y));
%m=max(z);
plot(f,z,'b');
%axis([-0.2,0.2,0,0.5]);
xlabel('Notmalized frequency');ylabel('Power (arb. units)')
%xlim([-.08 .08])
T=T*1e9;
figure;
plot(T,(Power),'k');
xlabel('Time (ns)');ylabel('Output power (mW)'); 
 figure;
%Stability_window = T/2:T; % Time window with stable dynamics
y = Y(:,1);
L = length(T);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
Y = fftshift(Y);
Y=Y/(max(Y));
lambda=linspace(1546,1554,NFFT);
window = 1000; % Part of spectrum to see
plot(lambda(NFFT/2-window:NFFT/2+window),log(abs(Y(NFFT/2-window:NFFT/2+window))),'b','LineWidth',1)
xlabel('\lambda(nm)')
ylabel('Power (dB)')
set(gca,'FontSize',12,'FontWeight','Bold')
%ylim([-10 1])
return
function dy=rate_equation2(t,y,n)
dy = zeros(n,1);
Jth=0.5e-9;q=1.602176487e-19*Jth; d=0.3e-6; D1=5e-6;L=250e-6;ts=3e-9;c=3e8; ng=4; nr=3.4;eta=.5; N0=8.25e5;
alpha0=2000; lambda0=1310e-9; Dlc=0.845e-9 ; Dld=60e-9;R=0.3;alpha=alpha0+(1.0/L)*log(1.0/R);
J=15e-3;
    sum1=0.0;
    lambda=zeros(n,1);
    gamma=zeros(n,1);
    D=zeros(n,1);
    g=zeros(n,1);
  
for i=2:n
 lambda(i)=lambda0+(i-1)*Dlc;
% other equations
end

return