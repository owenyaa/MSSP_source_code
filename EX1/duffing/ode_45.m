% load 'xx1.mat';
clear;%close all;
global M A r1 r4 f w
r1=2;r4=5;f=5;w=0.5;
N_dof=5;
m1=1;m2=1;m3=1;m4=1;m5=1;
c1=0.02*pi;c2=0.05*pi;c3=0.03*pi;c4=0.07*pi;c5=0.09*pi;
c6=0.01*pi;c7=0.04*pi;c8=0.06*pi;c9=0.08*pi;c10=0.1*pi;
k1=0.2*pi^2;k2=0.8*pi^2;k3=0.5*pi^2;k4=0.7*pi^2;k5=1.2*pi^2;
k6=0.4*pi^2;k7=1*pi^2;k8=0.3*pi^2;k9=1.1*pi^2;k10=0.6*pi^2;

M=[m1,0,0,0,0;0,m2,0,0,0;0,0,m3,0,0;0,0,0,m4,0;0,0,0,0,m5];
C=[c2+c3+c4,0,-c4,-c3,0;0,c7+c8,-c8,0,0;-c4,-c8,c4+c5+c6+c8+c9+c10,-c5,-c10;...
    -c3,0,-c5,c1+c3+c5,0;0,0,-c10,0,c10];
K=[k2+k3+k4,0,-k4,-k3,0;0,k7+k8,-k8,0,0;-k4,-k8,k4+k5+k6+k8+k9+k10,-k5,-k10;...
    -k3,0,-k5,k1+k3+k5,0;0,0,-k10,0,k10];

A=[zeros(N_dof),eye(N_dof);-M\K,-M\C];
Tdata=0:0.01:100*2*pi/w;
% ini_x=[x(1,1);x(2,1);dx(1,1);dx(2,1)];
ini_x=[100;200;300;0;0;0;0;0;0;0];
[t,num]=ode45('sub_ode45',Tdata,ini_x);
figure;
plot(Tdata,num(:,1),'r-','LineWidth',1.5);
hold on;
plot(Tdata,num(:,2),'k-','LineWidth',1.5);
hold on;
plot(Tdata,num(:,3),'b-','LineWidth',1.5);
hold on;
plot(Tdata,num(:,4),'g-','LineWidth',1.5);
hold on;
plot(Tdata,num(:,5),'g--','LineWidth',1.5);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% h=0.01;
% data=num(2000:end,5);
% N=length(data);
% N_fft=2^14;
% Y=fft(data,N_fft);
% Pyy=2*abs(Y(1:N_fft/2+1))/N_fft;
% f1=1/h*(0:N_fft/2)/N_fft;
% figure;
% plot(2*pi*f1,(Pyy(1:(N_fft/2+1))),'k-')
figure;
plot(num(40000:end,1),num(40000:end,6),'r-','LineWidth',1.5);
hold on;
plot(num(40000:end,2),num(40000:end,7),'k-','LineWidth',1.5);
hold on;
plot(num(40000:end,3),num(40000:end,8),'b-','LineWidth',1.5);
hold on;
plot(num(40000:end,4),num(40000:end,9),'g-','LineWidth',1.5);
hold on;
plot(num(40000:end,5),num(40000:end,10),'g--','LineWidth',1.5);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



