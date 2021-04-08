clear;
clc;
T=1;                                   %final time       
T0=0.01;
alpha=1.5;
h=1/160;
dt = 0.5*h^alpha; 
dtdx = dt/h;
t=T0:dt:T;
M=round((T-T0)/dt);

load Tem1evy15T1.mat;

Xt=zeros(M,1);
Yt=zeros(M,1);
Xt(1)=-1;
Yt(1)=-1;
for k=1:M-1
Xt(k+1)=Xt(k)+(Xt(k)-Xt(k)^3)*dt+Z1(k+1)-Z1(k);
Yt(k+1) = Yt(k) + cos(Xt(k))*dt + randn(1,1)*sqrt(dt);
end

save TemBM1o515T1 Yt
subplot(2,1,1)
 plot(t(2:end),Xt,'b')
% hold on
subplot(2,1,2)
plot(t(2:end),Yt,'b')