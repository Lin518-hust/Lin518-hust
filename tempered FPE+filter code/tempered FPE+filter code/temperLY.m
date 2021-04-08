% Author Ernest Jum
% Department of Mathematics
% University of Tennessee
% Knoxville
% Simulating alpha   stable and
% tempered alpha   stable processes
% August 2013
clc;
clear;
h=1/160;
to=0.01; tf=1;
alpha=1.5; %index of stability
 dt = 0.5*h^alpha; 
m=(tf-to)/dt+1;
% m=1000;
% dt=(tf-to)/(m-1);
t=to:dt:tf;
epsilon=0.001; %truncation level (precision)
lambda=.01; %tempering parameter
kappa=1/(2*abs(gamma(-alpha))); %constant C_alpha
tau=1/epsilon^(alpha); %truncation level
k=0;
T(1)= log(rand);
while sum(T(1:k)) < tau
k=k+1;
T(k)=-log(rand);
U(k)=rand;
s=rand;
if s<0.5
    V(k)=1;
else
    V(k)=-1;
end
    eta(k)= -log(rand)/lambda;
    xi(k)=rand;
    eta_xi(k)=eta(k)*((xi(k))^(1/alpha));
end
%% Jump sizes process %%%
const=(alpha/(2*kappa*tf))^(-1/alpha);
%%
Z(1)=0; Z1(1)=0;
for n=2:length(t)
    for i=1:length(T)
    gamma(i)=sum(T(1:i));
        if U(i) <=t(n)
%         J(i)=V(i)*(gamma(i)^(-1/alpha));
        J1(i)=V(i)*(min(const*((gamma(i))^(-1/alpha)),eta_xi(i)));
        else
%     J(i)=0;
    J1(i)=0;
    end
end
% Z(n)=sum(J(1:length(T)));
Z1(n)=sum(J1(1:length(T)));
end
save Tem1evy15T1 Z1
figure(1)
% plot(t, Z,'k', t, Z1)
% hleg1 = legend('\alpha-stable','tempered-\alpha-stable');
plot(t, Z1)
hleg1 = legend('tempered-\alpha-stable');
xlabel('time')
ylabel('smaple path')