clear;
clc;

xmin=-0.005;
xmax=10;
ymin=-0.005;
ymax=10;

global p q r1 r2
p=1.1;
q=0.2;
r1=1;
r2=1;
% D=0.01;

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

Nphi=80000;            % 环上划分精度
Sphi=zeros(1,Nphi);    % 所有的phi中S的最小值
phi=linspace(0,2*pi,Nphi);
% phi=linspace(1.893,1.894,Nphi);

D=[r1^2 0;0 r2^2];
C=[-xnode(2) -xnode(1);p*xnode(2) p*(xnode(1)-(1+q)/(q+xnode(2)))+p*xnode(2)*(1+q)/(q+xnode(2))^2];
A=[2*C(1,1) 2*C(1,2) 0;C(2,1) C(1,1)+C(2,2) C(1,2);0 2*C(2,1) 2*C(2,2)];
B=-[D(1,1);D(1,2);D(2,2)];
s=A\B;
Z=inv([s(1) s(2);s(2) s(3)]);

R=1e-2;
Nmap=1;
tf=200;
h=0.002;
nT=tf/h;
Np=zeros(1,Nphi);
xlamS=zeros(5,Nphi);
xlamS(1:2,:)=[xnode(1)+R*cos(phi);xnode(2)+R*sin(phi)];
% xlamS(3:4,:)=0;
xlamS(3:4,:)=M*[R*cos(phi);R*sin(phi)];
xlamS(5,:)=1/2*(Z(1,1)*(xlamS(1,:)-xnode(1)).^2+Z(2,2)*(xlamS(2,:)-xnode(2)).^2+2*Z(1,2)*(xlamS(1,:)-xnode(1)).*(xlamS(2,:)-xnode(2)));
% x1=zeros(Nphi,nT);
% x2=zeros(Nphi,nT);
% x3=zeros(Nphi,nT);
% x4=zeros(Nphi,nT);
% x5=zeros(Nphi,nT);
xyS=zeros(3,Nphi);
xyS2=zeros(3,Nphi);
pos=1:1:Nphi;
pos0=1:1:length(separatrix(1,:));
delta=0;
m1=min(separatrix(1,:));
m2=max(separatrix(1,:));

for j=1:nT
    t0=(j-1)*h;
    xlamS2=rk4(t0,h,xlamS);

%%% 到separatrix终止
    I=[];
    for k=1:length(pos)
%         [m,I2]=min(abs(separatrix(1,:)-xlamS2(1,k)));
%         if xlamS2(2,k)<separatrix(2,I2)
%             I=[I k];
%         end
        
        if (xlamS2(1,k)>m1)&&(xlamS2(1,k)<m2)
            J=separatrix(1,:)>xlamS2(1,k);
            J0=pos0(J);
            J1=J0(1);
            zx=xlamS2(1,k);
            z1=separatrix(:,J1-1);
            z2=separatrix(:,J1);
            zy=z1(2)+(z2(2)-z1(2))/(z2(1)-z1(1))*(zx-z1(1));
            if xlamS2(2,k)<zy
                I=[I k];
            end
        end
    end
    
    I1=find((xlamS2(1,:)>xmax-delta)|(xlamS2(1,:)<delta)|(xlamS2(2,:)<ymin+delta)|(xlamS2(2,:)>ymax-delta));
    I1=[I1 I];
    if isempty(I1)==0
        xyS2(1,pos(I1))=xlamS2(1,I1);
        xyS2(2,pos(I1))=xlamS2(2,I1);
        xyS2(3,pos(I1))=xlamS2(5,I1);
        
        xyS(1,pos(I1))=xlamS(1,I1);
        xyS(2,pos(I1))=xlamS(2,I1);
        xyS(3,pos(I1))=xlamS(5,I1);
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

[m,Iopt]=min(xyS(3,:));

figure;
plot(phi,xyS(3,:));

figure;
plot(xyS(1,:),xyS(3,:),'g.');
