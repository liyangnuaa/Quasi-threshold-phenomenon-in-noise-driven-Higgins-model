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
D=3e-4;

fxy=@(x)[1-x(1)*x(2);p*x(2)*(x(1)-(1+q)/(q+x(2)))];
xnode1=fsolve(fxy,[1;1]);
% xnode2=fsolve(fxy,[1;0]);
% xsad=fsolve(fxy,[0;0]);
xnode=xnode1;

separatrix= dlmread('separatrix3.txt');

nT=1e9;
Nphi=100;
ninitial=1e5;
h=0.0005;
x=zeros(Nphi,ninitial);
y=zeros(Nphi,ninitial);
totlenoisex=zeros(Nphi,ninitial);
totlenoisey=zeros(Nphi,ninitial);
pos=1:1:Nphi;
t0=0;

parfor j=1:Nphi
    fxy=@(x)[1-x(1)*x(2);p*x(2)*(x(1)-(1+q)/(q+x(2)))];
    xnode1=fsolve(fxy,[1;1]);
    % xnode2=fsolve(fxy,[1;0]);
    % xsad=fsolve(fxy,[0;0]);
    xnode=xnode1;
    
    separatrix= dlmread('separatrix3.txt');
    
    X=zeros(1,ninitial);
    Y=zeros(1,ninitial);
    noiseX=zeros(1,ninitial);
    noiseY=zeros(1,ninitial);
    X(1)=xnode(1);
    Y(1)=xnode(2);
    N=1;
    
    for i=1:nT
        x1=[X(N);Y(N)];
        noise=sqrt(D*h)*randn(2,1);
        x2=rk4_3(t0,h,x1,noise);
        if N==ninitial
            N=N-1;
            X(1:N)=X(2:N+1);
            Y(1:N)=Y(2:N+1);
            X(N+1)=x2(1);
            Y(N+1)=x2(2);
            noiseX(1:N)=noiseX(2:N+1);
            noiseY(1:N)=noiseY(2:N+1);
            noiseX(N+1)=noise(1);
            noiseY(N+1)=noise(2);
        else
            X(N+1)=x2(1);
            Y(N+1)=x2(2);
            noiseX(N+1)=noise(1);
            noiseY(N+1)=noise(2);
        end
        
        [m,I2]=min(abs(separatrix(1,:)-x2(1)));
        if x2(2)<separatrix(2,I2)
            break;
        end
        
        N=N+1;
    end
    
    x(j,:)=X;
    y(j,:)=Y;
    totlenoisex(j,:)=noiseX;
    totlenoisey(j,:)=noiseY;
end
Iw=find(x(:,end)==0);
x(Iw,:)=[];
y(Iw,:)=[];

Xfinal=x;
Yfinal=y;
Tfinal=((1:ninitial)-1)*h;
Tfinal=Tfinal-max(Tfinal);

n=length(Tfinal)-1;
% T0=zeros(1,n+1);
q1=zeros(100,n+1);
q2=zeros(100,n+1);
pdf1=zeros(100,n+1);
pdf2=zeros(100,n+1);
xopt=zeros(2,n+1);
for i=1:n+1
    [k1,k2]=ksdensity(Xfinal(:,i));
    [~,I]=max(k1);
    xopt(1,i)=k2(I);
    q1(:,i)=k2';
    pdf1(:,i)=k1';
    [k1,k2]=ksdensity(Yfinal(:,i));
    [~,I]=max(k1);
    xopt(2,i)=k2(I);
    q2(:,i)=k2';
    pdf2(:,i)=k1';
    
%     [k1,k2]=ksdensity(T(:,i));
%     [~,I]=max(k1);
%     T0(i)=k2(I);
end

% % Àë³öÎ»ÖÃ·Ö²¼
% L=length(separatrix(1,:));
% arclength=zeros(1,L);
% for i=2:L
%     arclength(i)=arclength(i-1)+norm(separatrix(:,i)-separatrix(:,i-1));
% end
% x1=Xfinal(:,end);
% E=[];
% L=length(x1);
% for i=1:L
%     I=separatrix(1,:)>x1(i);
%     n=length(separatrix(1,I));
%     E=[E arclength(n)];
% end
% [EDPSTO,arc]=ksdensity(E);
% figure;
% plot(arc,EDPSTO);

% ÔëÉùÊµÏÖ
NOISEX=zeros(1,length(x(1,:)));
NOISEY=zeros(1,length(x(1,:)));
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.003,'DesignMethod','butter');
for k=1:Nphi
    noisex= filtfilt(d1,totlenoisex(k,:)/h);
    noisey= filtfilt(d1,totlenoisey(k,:)/h);
    NOISEX=NOISEX+noisex;
    NOISEY=NOISEY+noisey;
end
NOISEX=NOISEX/Nphi;
NOISEY=NOISEY/Nphi;
% d1 = designfilt('lowpassiir','FilterOrder',12, ...
%     'HalfPowerFrequency',0.005,'DesignMethod','butter');
% noisex= filtfilt(d1,NOISEX);
% noisey= filtfilt(d1,NOISEY);
figure;
plot(Tfinal,NOISEX);
% hold on
% plot(Tfinal,noisex);
% hold off

figure;
plot(Tfinal,NOISEY);
% hold on
% plot(Tfinal,noisey);
% hold off

z=linspace(6,6,n+1);
figure;
plot3(Tfinal,xopt(1,:),z,'r-');

figure;
plot(xopt(1,:),xopt(2,:),'c-');
xlabel('x1');
ylabel('x2');

tp=Tfinal+max(abs(Tfinal));
figure;
mesh(tp,q1,pdf1);
hold on
plot3(tp,xopt(1,:),z,'r-');
hold off

figure;
mesh(Tfinal,q2,pdf2);

figure;
plot(q1(:,end),pdf1(:,end));


function xout=rk4_3(t0,h,x0,noise)
% k1=h*fun(t0,x0);
% k2=h*fun(t0+h/2,x0+0.5*k1);
% k3=h*fun(t0+h/2,x0+0.5*k2);
% k4=h*fun(t0+h,x0+k3);
% xout=x0+(k1+2*k2+2*k3+k4)/6;
% xout(2,:)=xout(2,:)+x0(3,:);

% D=2e-4;

% noise=sqrt(D*h)*randn(2,1);
k1=h*fun3(t0,x0);
x1=x0+k1+noise;
k2=h*fun3(t0+h,x1);
xout=x0+(k1+k2)/2+noise;


function fout=fun3(t0,x0)

% global p q
p=1.1;
q=0.2;

fout=[1-x0(1)*x0(2);p*x0(2)*(x0(1)-(1+q)/(q+x0(2)))];
