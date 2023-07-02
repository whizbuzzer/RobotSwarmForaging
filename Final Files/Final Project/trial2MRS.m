clear all 
clc
%Here we assume population after certain time period
N=1;
x0= [N/4 N/4 N/N N/5 N/5 N/N 0]'
[t,x]= ode45(@model_ode,[0 400], x0);

figure(1)
hold on
grid on
plot(t,x(:,2),'k', 'Linewidth',2)
plot(t,x(:,4),'b', 'Linewidth',2)
plot(t,x(:,3),'r', 'Linewidth',2)
plot(t,x(:,5),'m', 'Linewidth',2)
xlabel('Time')
ylabel('Population in each state')
%%
%Matrix for reaction vectors (mji=mj-mi)
R= [1 -1 0 0 0 0 0;
    -1 1 0 0 0 0 0;
    1 0 1 -1 0 0 0;
    0 1 0 1 -1 0 0;
    -1 0 0 1 0 1 0;
    1 0 0 0 1 0 -1];

a1=2; a2=4;
a3=1.2; b1=0.7;
ki=1.5; kj=0.4;

Ntot=30; %Total number of robots

c1=ki;
c2=kj;
c3=a1/Ntot;
c4=b1/Ntot;
c5=a2/Ntot;
c6=a3/Ntot;
t=0;
tfinal=1000;
N=Ntot*x0;
N(3)=0;
N(6)=0;
tvec=t;
Nvec=N/Ntot;

while t< tfinal
    % Reaction propensties
a(1)=c1*N(1);
a(2)=c2*N(2);
a(3)=c3*N(1)*N(3);
a(4)=c4*N(4)*N(2);
a(5)=c5*N(4)*N(6);
a(6)=c6*N(1)*N(5);

asum=sum(a);
j=min(find(rand<cumsum(a/asum)));
tau=log(1/rand)/asum;
N=N+R(j,:);
t=t+tau;
tvec(end+1)=t;
Nvec(end+1,:)=N/Ntot;

end
figure(1)
hold on
plot (tvec, Nvec(:,2),'k', 'Linewidth',2)
plot (tvec, Nvec(:,4),'b', 'Linewidth',2)
plot (tvec, Nvec(:,3),'r', 'Linewidth',2)
plot (tvec, Nvec(:,5),'m', 'Linewidth',2)

%%
function dx = model_ode(t,x)
M=[ 1 0 1 0 0 0 0 1 0; %Wm
    0 1 0 0 1 0 0 0 0; %Ba
    0 0 1 0 0 0 0 0 0; %F
    0 0 0 1 1 0 1 0 0; %Wn
    0 0 0 0 0 1 0 1 0; %Bb
    0 0 0 0 0 0 1 0 0; %N
    0 0 0 0 0 0 0 0 1];%Wo

%Reaction coefficients corresponding the equations
%Taking arbitrary values
a1=2; a2=4;
a3=1.2; b1=0.7;
ki=1.5; kj=0.4;
% This should be a 9x9 Matrix
K=[ki -kj 0 0 0 0 -a2 0 0 ;
   -ki kj 0 0 0 0 0 0 0;
    0 0 a1 0 0 0 0 0 0;
    0 0 -a1 0 0 0 0 0 0;
    0 0 0 0 b1 0 0 0 0;
    0 0 0 0 -b1 0 0 0 0;
    0 0 0 0 0 0 a2 0 0;
    0 0 0 0 0 0 0 -a3 0
    0 0 0 0 0 0 0 -a3 0];
   
% % y(x) is a 11x1 matrix
% Species Population
Wm= x(1);
Ba= x(2);
F = x(3);
Wn= x(4);
Bb= x(5);
N = x(6);
Wo= x(7);

Y=[Wm Ba Wm*F Wn Wn*Ba Bb Wn*N Wm*Bb Wo].';
dx= -M*K*Y;
end