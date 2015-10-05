function [] = code3()
L=10;
a=1;
b=3;
Q0=50;
k=0.5;
kv=10;
rho=1;
C=1;
T_out=300;
T_0=400;
N=40;
h=L/N;
x=linspace(0,L,N+1);

A=sparse(N+1);
b=zeros(N+1,1);

%T(0)=T_0
A(1,1)=1;
b(1,1)=T_0;

%T(L) condition
A(end,end)=-(k/h + kv);
A(end,end-1)=k/h;
b(end,1)=-kv*T_out;

for i=2:N
	A(i,i)=2*k;
	A(i,i-1)=-k - v*rho*C*h/2;
	A(i,i+1)=-k + v*rho*C*h/2;
end

index = find(a<x<b,1);% should be ok

b(index)=Q(x(index));% TODO Q must accept vector input

solution = A\b;

plot(x,solution);
end

function y = Q(x)
% x belongs to [0,10]
a=1;
b=3;
Q0=50;
y=zeros(size(x));
	y=Q0*sin(pi*(x-a)./(b-a));
end











