function [] = code3()
% Function for LAB3 
% Authors: David Weicker and Florentin Goyens
close all;


% TEST kkkkkkkkkkkk

%no speed and varible N
N = [10 20 40 80];
v = 0;
T = {0,0,0,0};
z = {0,0,0,0};
for i = 1:4
    [T{i},z{i}] = temp(N(i),v);
end
figure;
plot(z{1},T{1},'b',z{2},T{2},'g',z{3},T{3},'r',z{4},T{4},'k');title('Temperature in the pipe for v=0 and various stepsizes','FontSize',14)
xlabel('z-axis [m]','FontSize',14);ylabel('Temperature [K]','FontSize',14);legend('N=10','N=20','N=40','N=80');

%fixed N=40 and variable speed
N = 40;
v = [0.1 0.5 1 10];
for i = 1:4
    [T{i},z{i}] = temp(N,v(i));
end
figure;
plot(z{1},T{1},'b',z{2},T{2},'g',z{3},T{3},'r',z{4},T{4},'k');title('Temperature in the pipe for N=40 and various speeds','FontSize',14);
xlabel('z-axis [m]','FontSize',14);ylabel('Temperature [K]','FontSize',14);legend('v=0.1','v=0.5','v=1','v=10');

%large v=10 and variable N
N = [10 20 40];
v = 10;
for i = 1:3
    [T{i},z{i}] = temp(N(i),v);
end
figure;
plot(z{1},T{1},'b',z{2},T{2},'g',z{3},T{3},'r');title('Temperature in the pipe for large speed (v=10) and various stepsizes','FontSize',14);
xlabel('z-axis [m]','FontSize',14);ylabel('Temperature [K]','FontSize',14);legend('N=10','N=20','N=40');

end

function q = Q(z)
% Right-hand side of the system
% z is the discretized axis
[a,b,~,~,~,~,~,Q0,~,~,~]=param();

q = (a<=z) & (z<=b);
q = Q0*q.*sin(pi*(z-a)/(b-a));
end

function [T,z] = temp(N,v)
%Solve the problem for N+1 points between 0 and 10 and for speed v
[~,~,L,k,kv,rho,C,~,T0,Tout,~]=param();

h = L/N;
z = linspace(0,L,N+1);
q = Q(z(2:end))*h*h;

alpha = -k-v*rho*C*h/2; %tridiagonal terms
beta = 2*k;
gamma = -k+v*rho*C*h/2;

delta = 2*h*kv/k; %boundary conditions
q(1) = q(1) - alpha*T0;
q(end) = q(end) - gamma*delta*Tout;

e = ones(N,1);
A = spdiags([alpha*e beta*e gamma*e],-1:1,N,N); %sparse matrix
A(N,N-1) = A(N,N-1) + gamma;
A(N,N) = A(N,N) - gamma*delta;

T = T0*ones(N+1,1);
T(2:end) = A\q';
end

function [a,b,L,k,kv,rho,C,Q0,T0,Tout,v] = param()
%Initialize all the parameters for the model
L=10;
k=0.5;
kv=10;
v=0;
rho=1;
C=1;
a=1;
b=3;
Q0=50;
T0=400;
Tout=300;
end









