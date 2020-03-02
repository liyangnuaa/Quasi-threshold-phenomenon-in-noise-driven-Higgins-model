clear;
clc;

xmin=-0.01;
xmax=6;
ymin=-0.01;
ymax=6;

p=1.1;
q=0.2;
r1=1;
r2=1;
D=0.01;                                                                       % 参数

fxy=@(x)[1-x(1)*x(2);p*x(2)*(x(1)-(1+q)/(q+x(2)))];
xnode1=fsolve(fxy,[1;1]);
% xnode2=fsolve(fxy,[1;0]);
% xsad=fsolve(fxy,[0;0]);
xnode=xnode1;

N=5000;                                                                      % x、y方向离散化精度
Nc=N^2;
Ntheta=10;
xlin=linspace(xmin,xmax,N);
ylin=linspace(ymin,ymax,N);
h1=(xmax-xmin)/(N-1);
h2=(ymax-ymin)/(N-1);
h=norm([h1;h2]);                                                           % 相邻点最大距离
K=20;                                                                       % update radius
[X,Y]=meshgrid(xlin,ylin);

U=zeros(1,Nc);                                                             % quasipotential
U(:)=inf;
Far=1:Nc;                                                                  % 未计算节点
lx=floor((xnode(1)-xmin)/h1)+1;
ly=floor((xnode(2)-ymin)/h2)+1;
labelxnode=(ly-1)*N+lx;

n=3;
level=n;
labelx=lx-level+1:lx+level;
labely=ly-level+1:ly+level;
Accepted=[];
for i=1:2*level
    Accepted=[Accepted (labely(i)-1)*N+labelx];
end
Acceptedlogic=false(N,N);
Acceptedlogic(labelx,labely)=true;

level=n+1;
labelx=lx-level+1:lx+level;
labely=ly-level+1:ly+level;
AcceptedFront=[];
for i=1:2*level
    AcceptedFront=[AcceptedFront (labely(i)-1)*N+labelx];
end
AcceptedFrontlogic=false(N,N);
AcceptedFrontlogic(labelx,labely)=true;

level=n+2;
labelx=lx-level+1:lx+level;
labely=ly-level+1:ly+level;
Considered=[];
for i=1:2*level
    Considered=[Considered (labely(i)-1)*N+labelx];
end
Consideredlogic=false(N,N);
Consideredlogic(labelx,labely)=true;

Farlogic=true(N,N);
Farlogic(labelx,labely)=false;
Consideredlogic=Consideredlogic&(~AcceptedFrontlogic);
AcceptedFrontlogic=AcceptedFrontlogic&(~Acceptedlogic);

L=ismember(Considered,AcceptedFront);
Considered(L)=[];
L=ismember(AcceptedFront,Accepted);
AcceptedFront(L)=[];

unFar=[Accepted AcceptedFront Considered];
L=ismember(Far,unFar);
Far(L)=[];                                                                 % Far的初始化

D=[1 0;0 1];
C=[-xnode(2) -xnode(1);p*xnode(2) p*(xnode(1)-(1+q)/(q+xnode(2)))+p*xnode(2)*(1+q)/(q+xnode(2))^2];
A=[2*C(1,1) 2*C(1,2) 0;C(2,1) C(1,1)+C(2,2) C(1,2);0 2*C(2,1) 2*C(2,2)];
B=-[D(1,1);D(1,2);D(2,2)];
s=A\B;
S=[s(1) s(2);s(2) s(3)];
L=[Accepted AcceptedFront];
for k=1:length(L)                                                          % AcceptedFront的quasipotential的初始化
    xlabel=mod(L(k)-1,N)+1;
    ylabel=floor((L(k)-1)/N)+1;
    xn=[xlin(xlabel);ylin(ylabel)];
    U(L(k))=1/2*(xn-xnode)'/S*(xn-xnode);
end

LC=length(Considered);
WC=zeros(1,LC);
for k=1:LC                                                                 % Considered的quasipotential的初始化
    z=Considered(k);
    WC(k)=CalculateW(z,N,xlin,ylin,AcceptedFront,U,Ntheta,h,K,p,q);
end

for k=1:Nc
    [M,I]=min(WC);
    z=Considered(I);
    U(z)=M;                                                                % 将最小W赋给U
    
%     if z==labelsaddle+N-1
%         zpre=zp;
%     end
%     zp=z;
    
    zxlabel=mod(z-1,N)+1;
    zylabel=floor((z-1)/N)+1;
    if (zxlabel==1)||(zxlabel==N)||(zylabel==1)||(zylabel==N)||(M>10)      % 终止条件
        break;
    end
    
    WC(I)=[];
    Considered(I)=[];
    Consideredlogic(zxlabel,zylabel)=false;
    Con=Considered;
    AcceptedFront=[AcceptedFront z];                                       % 将具有最小W的Considered赋给AcceptedFront
    AcceptedFrontlogic(zxlabel,zylabel)=true;
    
%     zFarad=neighbor(z,Far,N);                                              % 计算Far里面与z邻接的点的W
%     In=ismember(Far,zFarad);
%     Far(In)=[];
    zFarad=[];
    L=[zxlabel-1 zxlabel zxlabel+1 zxlabel-1 zxlabel+1 zxlabel-1 zxlabel zxlabel+1;
        zylabel-1 zylabel-1 zylabel-1 zylabel zylabel zylabel+1 zylabel+1 zylabel+1];
    for i=1:8
        if Farlogic(L(1,i),L(2,i))
            zFarad=[zFarad N*(L(2,i)-1)+L(1,i)];
            Farlogic(L(1,i),L(2,i))=false;
            Consideredlogic(L(1,i),L(2,i))=true;
        end
    end
    Considered=[Considered zFarad];
    
%     zAFad=neighbor(z,AcceptedFront,N);                                     % 更新AcceptedFront
    zAFad=[];
    L=[zxlabel-1 zxlabel zxlabel+1 zxlabel-1 zxlabel+1 zxlabel-1 zxlabel zxlabel+1;
        zylabel-1 zylabel-1 zylabel-1 zylabel zylabel zylabel+1 zylabel+1 zylabel+1];
    for i=1:8
        if AcceptedFrontlogic(L(1,i),L(2,i))
            zAFad=[zAFad N*(L(2,i)-1)+L(1,i)];
        end
    end
    LzA=length(zAFad);
    for i=1:LzA
        zA=zAFad(i);
        zAxlabel=mod(zA-1,N)+1;
        zAylabel=floor((zA-1)/N)+1;
        L=[zAxlabel-1 zAxlabel+1 zAxlabel zAxlabel;zAylabel zAylabel zAylabel-1 zAylabel+1];
%         In=ismember(L,Considered);
%         if isempty(L(In))
%             I=(AcceptedFront==zA);
%             AcceptedFront(I)=[];
%             Accepted=[Accepted zA];
%         end
        if ~((Consideredlogic(L(1,1),L(2,1)))||(Consideredlogic(L(1,2),L(2,2)))||(Consideredlogic(L(1,3),L(2,3)))||(Consideredlogic(L(1,4),L(2,4))))
            I=(AcceptedFront==zA);
            AcceptedFront(I)=[];
            AcceptedFrontlogic(zAxlabel,zAylabel)=false;
            Accepted=[Accepted zA];
            Acceptedlogic(zAxlabel,zAylabel)=true;
        end
    end
    
%     zFarad=neighbor(z,Far,N);                                              % 计算Far里面与z邻接的点的W
    LzF=length(zFarad);
    Wextra=zeros(1,LzF);
    for i=1:LzF
        Wextra(i)=CalculateW(zFarad(i),N,xlin,ylin,AcceptedFront,U,Ntheta,h,K,p,q);
    end
    
    zCKh=withinKh(z,Con,K,h,N,xlin,ylin);                           % 更新WC
    LzC=length(zCKh);
    W=zeros(1,LzC);
    for i=1:LzC
        W(i)=CalculateW(zCKh(i),N,xlin,ylin,AcceptedFront,U,Ntheta,h,K,p,q);
    end
    
    for i=1:LzC
        In= Con==zCKh(i);
        if W(i)<WC(In)
            WC(In)=W(i);
        end
    end
    
%     In=ismember(Far,zFarad);
%     Far(In)=[];
%     Considered=[Considered zFarad];
    WC=[WC Wextra];
    
%     M
end

Uend=zeros(N,N);
for i=1:N
    Uend(i,:)=U((i-1)*N+1:i*N);
end

% global p q X Y Ux Uy
% x0=[3.7860806776;0.0841647096];
% % x0=[3.2223406;0.12958947];
% [Ux,Uy]=gradient(Uend);
% Ux=Ux/h1;
% Uy=Uy/h2;
% MPEP=x0;
% ht=0.001;
% T=10;
% Nstep=floor(T/ht);
% for i=1:Nstep
%     x1=rk4(0,ht,x0);
%     if norm(x1-xnode)<=2*h
%         break;
%     end
%     x0=x1;
%     MPEP=[x0 MPEP];
% end
% figure;
% plot(MPEP(1,:),MPEP(2,:),'m-');

figure;
mesh(X,Y,Uend);

figure;
contour(X,Y,Uend);

h1 = openfig('Fig. 7 quasipotential by OUM.fig','reuse'); % open figure
D1=get(gca,'Children'); %get the handle of the line object
XData1=get(D1,'XData'); %get the x data
YData1=get(D1,'YData'); %get the y data
ZData1=get(D1,'ZData'); %get the z data

figure;
contour(XData1,YData1,ZData1,45);

Nxy=10000;
xyl=linspace(1,3.157,Nxy);
Uxyl=zeros(1,Nxy);
xd=XData1(1,:);
yd=YData1(:,1)';
for i=1:Nxy
    Ix=find(xd>xyl(i),1,'first');
    Iy=find(yd>xyl(i),1,'first');
    z0=ZData1(Ix-1,Iy-1);
    z1=ZData1(Ix,Iy-1);
    z2=ZData1(Ix-1,Iy);
    Uxyl(i)=z0+(z1-z0)*(xyl(i)-xd(Ix-1))/(xd(Ix)-xd(Ix-1))+(z2-z0)*(xyl(i)-yd(Iy-1))/(yd(Iy)-yd(Iy-1));
end
figure;
plot(xyl,Uxyl);

figure;
plot(xyl,xyl,'k-');


%%% 插值S点拟势
XData1=X;
YData1=Y;
ZData1=Uend;

x0=[3.65690576;0.09241667];
xd=XData1(1,:);
yd=YData1(:,1)';
Ix=find(xd>x0(1),1,'first');
Iy=find(yd>x0(2),1,'first');
z0=ZData1(Iy-1,Ix-1);
z1=ZData1(Iy-1,Ix);
z2=ZData1(Iy,Ix-1);
U0=z0+(z1-z0)*(x0(1)-xd(Ix-1))/(xd(Ix)-xd(Ix-1))+(z2-z0)*(x0(2)-yd(Iy-1))/(yd(Iy)-yd(Iy-1));
