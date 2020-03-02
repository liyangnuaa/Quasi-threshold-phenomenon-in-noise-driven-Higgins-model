clear;
clc;

xmin=-0.005;
xmax=10;
ymin=-0.005;
ymax=10;

% global p q r1 r2 D Nnew noise
p=1.1;
q=0.2;
r1=1;
r2=1;
invD=1/3e-4;
D=1/invD;

fxy=@(x)[1-x(1)*x(2);p*x(2)*(x(1)-(1+q)/(q+x(2)))];
xnode1=fsolve(fxy,[1;1]);
% xnode2=fsolve(fxy,[1;0]);
% xsad=fsolve(fxy,[0;0]);
xnode=xnode1;

separatrix= dlmread('separatrix3.txt');

nT=1e9;
Nphi=100;
h=0.0003;
exittime=zeros(1,Nphi);
exitpoint=zeros(2,Nphi);
pos=1:1:Nphi;
t0=0;

parfor j=1:Nphi
    fxy=@(x)[1-x(1)*x(2);p*x(2)*(x(1)-(1+q)/(q+x(2)))];
    xnode1=fsolve(fxy,[1;1]);
    % xnode2=fsolve(fxy,[1;0]);
    % xsad=fsolve(fxy,[0;0]);
    xnode=xnode1;
    
    separatrix= dlmread('separatrix3.txt');
    
    x2=xnode;
    
    for i=1:nT
        x1=x2;
        noise=sqrt(D*h)*randn(2,1);
        x2=rk4_3(t0,h,x1,noise);
        
        [m,I2]=min(abs(separatrix(1,:)-x2(1)));
        if x2(2)<separatrix(2,I2)
            break;
        end
    end
    
    exittime(j)=h*i;
    exitpoint(:,j)=x2;
end

MFPT=sum(exittime)/Nphi;

% 离出位置分布
L=length(separatrix(1,:));
arclength=zeros(1,L);
for i=2:L
    arclength(i)=arclength(i-1)+norm(separatrix(:,i)-separatrix(:,i-1));
end
x1=exitpoint(1,:);
E=[];
L=length(x1);
for i=1:L
    I=separatrix(1,:)>x1(i);
    n=length(separatrix(1,I));
    E=[E arclength(n)];
end
[EDPSTO,arc]=ksdensity(E);
figure;
plot(arc,EDPSTO);
