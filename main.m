clear;
clc;

xmin=-0.005;
xmax=10;
ymin=-0.005;
ymax=10;

global p q r1 r2 D
p=1.1;
q=0.2;
r1=1;
r2=1;
D=3e-4;

fxy=@(x)[1-x(1)*x(2);p*x(2)*(x(1)-(1+q)/(q+x(2)))];
xnode1=fsolve(fxy,[1;1]);
% xnode2=fsolve(fxy,[1;0]);
% xsad=fsolve(fxy,[0;0]);
xnode=xnode1;

separatrix= dlmread('separatrix3.txt');

A=[-xnode(2) -xnode(1);p*xnode(2) p*(xnode(1)-(1+q)/(q+xnode(2)))+p*xnode(2)*(1+q)/(q+xnode(2))^2];
B=zeros(4,4);
B(1:2,1:2)=A;
B(1:2,3:4)=[r1^2 0;0 r2^2];
B(3:4,3:4)=-A';
[Bv,Be]=eig(B);
Bv1=Bv(1:2,3:4);
Bv2=Bv(3:4,3:4);
M=real(Bv2/Bv1);

Nphi=20000;            % 环上划分精度
Sphi=zeros(1,Nphi);    % 所有的phi中S的最小值
% phi=linspace(0,2*pi,Nphi);
phi=linspace(1.84,2.15,Nphi);
% phi=linspace(1.570735367683842,1.571035517758880,Nphi);

Dp=[r1^2 0;0 r2^2];
C=[-xnode(2) -xnode(1);p*xnode(2) p*(xnode(1)-(1+q)/(q+xnode(2)))+p*xnode(2)*(1+q)/(q+xnode(2))^2];
A=[2*C(1,1) 2*C(1,2) 0;C(2,1) C(1,1)+C(2,2) C(1,2);0 2*C(2,1) 2*C(2,2)];
B=-[Dp(1,1);Dp(1,2);Dp(2,2)];
s=A\B;
Z=inv([s(1) s(2);s(2) s(3)]);

R=1e-2;
Nmap=1;
tf=200;
h=0.001;
nT=tf/h;
Np=zeros(1,Nphi);
xlamS=zeros(10,Nphi);
xlamS(1:2,:)=[xnode(1)+R*cos(phi);xnode(2)+R*sin(phi)];
% xlamS(3:4,:)=0;
xlamS(3:4,:)=M*[R*cos(phi);R*sin(phi)];
xlamS(5,:)=Z(1,1);
xlamS(6,:)=Z(1,2);
xlamS(7,:)=Z(2,1);
xlamS(8,:)=Z(2,2);
xlamS(9,:)=1/2*(Z(1,1)*(xlamS(1,:)-xnode(1)).^2+Z(2,2)*(xlamS(2,:)-xnode(2)).^2+2*Z(1,2)*(xlamS(1,:)-xnode(1)).*(xlamS(2,:)-xnode(2)));
xlamS(10,:)=1;

xyp1p2SK=[];
pos=1:1:Nphi;
delta=0;

for j=1:nT
    t0=(j-1)*h;
    xlamS2=rk4(t0,h,xlamS);

%%% 到separatrix终止
    I=[];
    for k=1:length(pos)
        [m,I2]=min(abs(separatrix(1,:)-xlamS2(1,k)));
        if xlamS2(2,k)<separatrix(2,I2)
            I=[I k];
        end
    end
    xyp1p2SK=[xyp1p2SK [xlamS2(1:4,I);xlamS2(9:10,I)]];
    
    I1=find((xlamS2(1,:)>4.5-delta)|(xlamS2(1,:)<0.25)|(xlamS2(2,:)<ymin+delta)|(xlamS2(2,:)>ymax-delta));
    I1=[I1 I];
    if isempty(I1)==0
%         xyp1p2SK(1:4,pos(I1))=xlamS2(1:4,I1);
%         xyp1p2SK(5:6,pos(I1))=xlamS2(9:10,I1);
        pos(I1)=[];
        xlamS2(:,I1)=[];
    end
    
%     Np(pos)=Np(pos)+1;
%     
%     x1(pos,j)=xlamS2(1,:)';
%     x2(pos,j)=xlamS2(2,:)';
%     x3(pos,j)=xlamS2(3,:)';
%     x4(pos,j)=xlamS2(4,:)';
%     x5(pos,j)=xlamS2(5,:)';
    
    if isempty(pos)
        break;
    end
    xlamS=xlamS2;
end

D=3e-4;
L=length(separatrix(1,:));
arclength=zeros(1,L);
for i=2:L
    arclength(i)=arclength(i-1)+norm(separatrix(:,i)-separatrix(:,i-1));
end

EDP=[];
L=length(xyp1p2SK(1,:));
Wmodify=[];
for i=1:L
    I=separatrix(1,:)>xyp1p2SK(1,i);
    n=length(separatrix(1,I));
    m=length(separatrix(1,:))+1-n;
    k=(separatrix(2,m+1)-separatrix(2,m-1))/(separatrix(1,m+1)-separatrix(1,m-1));
    normal=[k;-1]/norm([k;-1]);
%     P=exp(-xyp1p2SK(5,i)/D);
    P=(normal(1)*xyp1p2SK(3,i)+normal(2)*xyp1p2SK(4,i))*xyp1p2SK(6,i)*exp(-xyp1p2SK(5,i)/D);
    EDP=[EDP [arclength(n);P]];
    
    Wmodify=[Wmodify xyp1p2SK(5,i)-D*xyp1p2SK(6,i)];
end

[EDP1,i0,i1]=unique(EDP(1,:));
EDP2=EDP(2,i0);
Wmod2=Wmodify(i0);

S=trapz(EDP1,EDP2);
prefactor=1/S;
EDP2=prefactor*EDP2;

figure;
plot(EDP1,EDP2,'r-');

figure;
plot(EDP1,Wmod2,'g');
% figure;
% plot(EDP(1,:),EDP(2,:),'r.');

% figure;
% plot(xyp1p2SK(1,:),EDP(2,:),'r.');

% figure;
% plot(phi,xyp1p2SK(5,:));

figure;
plot(xyp1p2SK(1,:),xyp1p2SK(5,:),'g.');



I=separatrix(1,:)>3.7861;
separatrix0=separatrix(:,I);
L=length(separatrix0(1,:));
arclength=zeros(1,L);
for i=2:L
    arclength(i)=arclength(i-1)+norm(separatrix0(:,i)-separatrix0(:,i-1));
end
