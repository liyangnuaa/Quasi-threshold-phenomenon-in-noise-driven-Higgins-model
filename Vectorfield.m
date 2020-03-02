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

[x,y]=meshgrid(xmin:0.3:xmax,ymin:0.3:ymax);
u=1-x.*y;
v=p*y.*(x-(1+q)./(q+y));

N=4;
nstep=100000;
% x1=[linspace(0.8,0.95,N);linspace(1,1,N)];
% x1=[linspace(0.82,0.85,N);linspace(1,1,N)];
x1=[linspace(1,1,N);linspace(0.82,0.85,N)];
X=zeros(N,nstep);
Y=zeros(N,nstep);
X(:,1)=x1(1,:)';
Y(:,1)=x1(2,:)';
h=0.001;
t0=0;
for i=1:nstep-1
    x2=rk4(t0,h,x1);
    x1=x2;
    X(:,i+1)=x2(1,:)';
    Y(:,i+1)=x2(2,:)';
end

figure;
for i=1:N
    plot(X(i,:),Y(i,:));
    hold on
end
hold off

figure;
quiver(x,y,u,v);
hold on
ezplot(f1);
axis([xmin xmax ymin ymax])
hold on
ezplot(f2);
% hold on
% plot(x1(1,:),x1(2,:),'m');
hold off

figure;
streamslice(x,y,u,v);
hold on
ezplot(f1);
axis([xmin xmax ymin ymax])
hold on
ezplot(f2);
% hold on
% plot(x1(1,:),x1(2,:),'m');
hold off

% figure;
% plot(xnode(1),xnode(2),'r+');



function xout=rk4(t0,h,x0)
k1=h*fun(t0,x0);
k2=h*fun(t0+h/2,x0+0.5*k1);
k3=h*fun(t0+h/2,x0+0.5*k2);
k4=h*fun(t0+h,x0+k3);
xout=x0+(k1+2*k2+2*k3+k4)/6;


function y=fun(~,x)

global p q

y=zeros(size(x));
y(1,:)=1-x(1,:).*x(2,:);
y(2,:)=p*x(2,:).*(x(1,:)-(1+q)./(q+x(2,:)));
