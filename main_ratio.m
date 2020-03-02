clear;
clc;

xmin=0;
xmax=10;
ymin=0;
ymax=10;

global p q r1 r2 D
p=1.1;
q=0.2;
r1=1;
r2=1;
D=0.01;

fxy=@(x)[1-x(1)*x(2);p*x(2)*(x(1)-(1+q)/(q+x(2)))];
% xnode1=fsolve(fxy,[1;1]);
% xnode2=fsolve(fxy,[1;0]);
% xsad=fsolve(fxy,[0;0]);

f1=@(x,y)1-x*y;
f2=@(x,y)p*y*(x-(1+q)/(q+y));

N=10002;
nstep=100000;
x0=[linspace(0.8,0.86,N);linspace(1,1,N)];
x1=x0;
h=0.001;
t0=0;
xm=zeros(1,N);
pos=1:N;
for i=1:nstep-1
    x2=rk4(t0,h,x1);
    I=x2(1,:)<x1(1,:);
    if isempty(I)==0
        xm(pos(I))=x1(1,I);
        x2(:,I)=[];
        pos(I)=[];
    end
    x1=x2;
    if isempty(x1)
        break;
    end
end

step=0.06/(N-1);
ratio=zeros(1,N-2);
for i=2:N-1
    ratio(i-1)=(xm(i-1)-xm(i+1))/step/2;
end
x=x0(1,2:N-1);

figure;
plot(x,ratio);

