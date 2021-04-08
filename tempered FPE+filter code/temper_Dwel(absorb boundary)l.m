% compute pdf of fpe corr to Levy noise with absorbing BC 
% in a bounded, symmetric domain (-r,r);
% central differencing for 2nd-order deriv 
% and using one-sided differencing near the boundaries;
% 3rd-order WENO for 1st-order deriv
% 3rd-order TVD RK in time
% 


% function [X,U]=temper_Dwell(Alp,T,x0)
% clear all;
% clf

alpha=1.5;
T=1;
x0=0;
eps=1;
d=0;
rr=1;             %right
rl=1;             %left
r=4;
f=@(x) x-x.^3;
fpmax=f(1/sqrt(3));
lam=.1*r;

h=2/160;
J=rr/h;             %right
L=rl/h;             %left
dt = 0.5*h^alpha; 
dtdx = dt/h;
%T=2;

Jt=L+J-1;     %Jt=total number of unknowns  
x=-(rl+rr):h:(rl+rr);
C=-zeta(alpha-1)*h^(2-alpha);          % correction term u''(x)
cons=1/(r^alpha*2*abs(gamma(-alpha)));
% coeff for the diffusion and the correction term
Chh=( d/2 + cons*C*eps )/h^2;
% c1=eps*cons/alpha;
% c2=eps*cons*h;
c1=eps*cons;
c2=eps*cons*h;

b=zeros(Jt,1);
a=zeros(Jt,1);
c=zeros(Jt,1); 

W1=lam^(0.5*(alpha-1))*(1+x(J+2:3*J)).^(-0.5*(alpha+1)).*exp(-0.5*lam*(1+x(J+2:3*J))).*whittakerW(-0.5*(alpha+1),-0.5*alpha,lam*(1+x(J+2:3*J)));
W2=lam^(0.5*(alpha-1))*(1-x(J+2:3*J)).^(-0.5*(alpha+1)).*exp(-0.5*lam*(1-x(J+2:3*J))).*whittakerW(-0.5*(alpha+1),-0.5*alpha,lam*(1-x(J+2:3*J)));
%nonintegral part

% coefficient of U_j
b(2:Jt-1) = (-2*Chh - c1*(W1(2:Jt-1)+W2(2:Jt-1)))';
% one-sided diff near boundaries
b(1)  = 2*Chh -  c1*(W1(1)+W2(1)); 
b(Jt) = 2*Chh - c1*(W1(Jt)+W2(Jt)); % one-sided diff
a= Chh*ones(Jt,1);  % coefficient of U_(j-1)
c= Chh*ones(Jt,1);  % coefficient of U_(j+1) 
c(1)  = -5*Chh; % one-sided diff
a(Jt) = -5*Chh; % one-sided diff
vp2 = zeros(Jt,1); vp2(3) = 4*Chh;  % one-sided diff
vp3 = zeros(Jt,1); vp3(4) =  -Chh; % one-sided diff 
vm2 = zeros(Jt,1); vm2(Jt-2) = 4*Chh;  % one-sided diff
vm3 = zeros(Jt,1); vm3(Jt-3) =  -Chh; % one-sided diff 




for j=-L+1:J-1
   b(j+L)= b(j+L) - c2*( sum(1./(exp(lam*abs(x(J+2-j:L+J))).*abs(x(J+2-j:L+J)).^(1+alpha))) ...
                       + sum(1./(exp(lam*abs(x(L+J+2:2*J+L-j))).*abs(x(L+J+2:2*J+L-j)).^(1+alpha))) ...
         + .5/(exp(lam*abs(x(J+1-j)))*(abs(x(J+1-j))^(1+alpha))) + .5/(exp(lam*abs(x(2*J+L+1-j)))*abs(x(2*J+L+1-j))^(1+alpha) ));  
end 

A=spdiags([vm3 vm2 [a(2:end); 0] b ...
           [0; c(1:end-1)] vp2 vp3],-3:3,Jt,Jt);

% coefficient of u_(j+k) 
B=zeros(size(A));
for j=-L+1:J-1
  B(L+j,:)=[1./(exp(lam*abs(x(J+2-j:L+J))).*abs(x(J+2-j:L+J)).^(1+alpha))  0  1./(exp(lam*abs(x(L+J+2:2*J+L-j))).*abs(x(L+J+2:2*J+L-j)).^(1+alpha))];
end


%A = diag(ones(Jt,1))+ dt*(A+c2*B); % the iterative matrix for each time-step

X=-r:r*h:r;

UU=sqrt(40/pi)*exp(-40*(X+x0).^2); %gaussian

U=UU(2:end-1)';
Un=U; U1=U; U2=U;
%U=UU(2:end-1);


nft=round(T/dt);
%figure
%hold on
nu = length(U);
data = zeros(nu+4,1);
load TemBM1o515T1.mat;
for nt=1:nft-1

    U1 = U + dt*(A+c2*B)*U;
    % global Lax-Friedrichs(LF) flux splitting
    data(3:nu+2) = (f(X(2:end-1)) + fpmax)'.*U/2;

    fx1 = derWENOr2_minus(data,h);

   data(3:nu+2) = (f(X(2:end-1)) - fpmax)'.*U/2;
    fx2 = derWENOr2_plus(data,h);
    U1 = U1 - dtdx*(fx1+fx2)/r;
    
    U2 = 0.75*U + U1/4 + (dt/4)*(A+c2*B)*U1;
    data(3:nu+2) = (f(X(2:end-1)) + fpmax)'.*U1/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(X(2:end-1)) - fpmax)'.*U1/2;
    fx2 = derWENOr2_plus(data,h);
    U2 = U2 - (dtdx/4)*(fx1+fx2)/r;
    
    Un = U/3 + 2*U2/3 + (2*dt/3)*(A+c2*B)*U2;
    data(3:nu+2) = (f(X(2:end-1)) + fpmax)'.*U2/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(X(2:end-1)) - fpmax)'.*U2/2;
    fx2 = derWENOr2_plus(data,h);
    Un = Un - (2*dtdx/3)*(fx1+fx2)/r;
  

            
         U=Un;
     %t=nt*dt;
    
     % F(nt) = getframe;
         
end
%save('plotdiffT_uniform','P1','P2','P3','P4')
%title(['final time t= ', num2str(t)])
% Ux = U; 
% 
% a=max(U);
% m=find(U==a)
 hold on;
 plot(X(2:end-1)',U,'b-.')

 
